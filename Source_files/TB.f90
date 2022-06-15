! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2022 Nikita Medvedev
!
! XTANT is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Although we endeavour to ensure that the code XTANT and results delivered are correct,
! no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions.
! We shall not be liable for any damage arising from the use of this code or its parts
! or any results produced with it, or from any action or decision taken
! as a result of using this code or any related material.
!
! This code is distributed as is for non-commercial peaceful purposes only,
! such as research and education. It is explicitly prohibited to use the code,
! its parts, its results or any related material for military-related and other than peaceful purposes.
!
! By using this code or its materials, you agree with these terms and conditions.
!
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with general aspects of TB

MODULE TB
use Universal_constants
use Objects
use Variables
use Algebra_tools
use Little_subroutines
use Atomic_tools
use Electron_tools
use Nonadiabatic
use TB_Fu
use TB_Pettifor
use TB_Molteni
use TB_NRL
use TB_DFTB, only : Construct_Vij_DFTB, construct_TB_H_DFTB, get_Erep_s_DFTB, get_dHij_drij_DFTB, &
                     Attract_TB_Forces_Press_DFTB, dErdr_s_DFTB, dErdr_Pressure_s_DFTB, Complex_Hamil_DFTB
use TB_3TB, only : get_Erep_s_3TB, dErdr_s_3TB, dErdr_Pressure_s_3TB, &
                     Construct_Vij_3TB, construct_TB_H_3TB, get_Mjs_factors
use TB_BOP, only : Construct_Vij_BOP, construct_TB_H_BOP, get_Erep_s_BOP
use TB_xTB, only : Construct_Vij_xTB, construct_TB_H_xTB, get_Erep_s_xTB, identify_xTB_orbitals_per_atom
use Van_der_Waals
use Coulomb
use Exponential_wall

implicit none

 contains

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Common TB tools:


! Attractive part:
subroutine get_Hamilonian_and_E(Scell, numpar, matter, which_fe, Err, t)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter ! material parameters
   integer, intent(in) :: which_fe ! which method is used to get electron distribution
   type(Error_handling), intent(inout) :: Err	! error save
   real(8), intent(in) :: t ! [fs] timestep
   !-----------------------
   integer NSC, i, Nat
   real(8), dimension(:,:,:), allocatable :: M_Vij	! matrix of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dVij	! matrix of derivatives of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_SVij	! matrix of Overlap functions for Overlap matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dSVij	! matrix of derivatives of Overlap functions for Overlap matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_E0ij     ! matrix of functions for on-site energies for all pairs of atoms, all orbitals (used for BOP)
   real(8), dimension(:,:,:), allocatable :: M_dE0ij     ! matrix of derivatives of functions for on-site energies for all pairs of atoms, all orbitals (used for BOP)
   real(8), dimension(:,:,:), allocatable :: M_Lag_exp	! matrix of Laguerres for 3-body radial funcs (for 3TB)
   real(8), dimension(:,:,:), allocatable :: M_d_Lag_exp	! matrix of derivatives of Laguerres (for 3TB)
   real(8), dimension(:,:,:), allocatable :: M_lmn	! matrix of directional cosines l, m, n
   real(8), dimension(:,:,:), allocatable :: M_x1	! matrix of x1 elements, used for forces calculations
   real(8), dimension(:,:,:), allocatable :: M_xrr	! matrix of xrr elements, used for forces calculations
   real(8), dimension(:,:,:), allocatable :: Mjs    ! matrix of K-S part of overlaps with s-orb. (for 3TB)
   real(8), dimension(:,:), allocatable :: M_Aij_x_Ei		! matrix of multipliers for non-orthogonal forces calculation
   real(8), dimension(:,:,:), allocatable :: M_drij_dsk  ! matrix of derivatives of distances between atoms i and j
   real(8), dimension(:,:,:), allocatable :: M_dlmn   ! matrix of derivatives of directional cosines

   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
      ! Create and diaganalize Hamiltonain (in Pettifor form):
      SC:do NSC = 1, size(Scell) ! for all supercells

         ! Get lists of nearest neighbors:
         call get_near_neighbours(Scell, numpar) ! see "Atomic_tools"
!          call print_time('Starting H', ind=1)
         
         ! Get and save in matrices often-used elements for forces calculations:
         ! Here, get the matrix of directional cosines and associated expressions:
         call Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn) ! see below 
         
         ! Electronic TB Hamiltonian part:
         ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types
            !print*, 'before construct_TB_H'
            
            ! Create and diagonalize TB Hamiltonian:
            select type(ARRAY)
            type is (TB_H_Pettifor) ! TB parametrization according to Pettifor
               call construct_TB_H_Pettifor(numpar, matter, ARRAY, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Pettifor"
            type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni
               call construct_TB_H_Molteni(numpar, matter, ARRAY, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Molteni"
            type is (TB_H_Fu)  ! TB parametrization accroding to Fu
               call construct_TB_H_Fu(numpar, matter, ARRAY, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Fu"
            type is (TB_H_NRL)  ! TB parametrization accroding to NRL
               call Construct_Vij_NRL(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_NRL"
               call construct_TB_H_NRL(numpar, matter, ARRAY, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_NRL"
            type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB
               call Construct_Vij_DFTB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_DFTB"
               call construct_TB_H_DFTB(numpar, matter, ARRAY, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_DFTB"

            type is (TB_H_3TB)  ! TB parametrization accroding to 3TB
               ! Get the overlaps between orbitals and ficticios s orbital (for 3-body parts):
               call get_Mjs_factors(numpar%N_basis_size, Scell(NSC), M_lmn, Mjs, M_drij_dsk, M_dlmn)   ! module "TB_3TB"
               ! Get the overlaps and reusable functions:
               call Construct_Vij_3TB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_Lag_exp, M_d_Lag_exp) ! module "TB_3TB"
               ! Construct the Hamiltonian, diagonalize it, get the energy:
               call construct_TB_H_3TB(numpar, matter, ARRAY, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Scell, NSC, Err) ! module "TB_3TB"

            type is (TB_H_BOP)  ! TB parametrization accroding to BOP (incomplete)
               call Construct_Vij_BOP(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_E0ij, M_dE0ij)    ! module "TB_BOP"
               call construct_TB_H_BOP(numpar, ARRAY, matter, M_Vij, M_SVij, M_E0ij, M_lmn, Scell, NSC, Err)    ! module "TB_BOP"
            type is (TB_H_xTB)  ! TB parametrization accroding to xTB (NOT READY)
!                call Construct_Vij_xTB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_xTB"
!                call construct_TB_H_xTB(numpar, matter, ARRAY, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_xTB"
            end select
            
            ! Get the DOS weights for each energy level, if required:
            !call get_DOS_weights(numpar%DOS_splitting, numpar%mask_DOS, numpar%DOS_weights, Hij=Scell(NSC)%Ha) ! below
            call get_DOS_weights(1, numpar%mask_DOS, numpar%DOS_weights, Hij=Scell(NSC)%Ha) ! below
            
            ! Fill corresponding energy levels + repulsive energy cotribution:
            select case (which_fe)
            case (1) ! distribution for given Te:
               call set_initial_fe(Scell, matter, Err) ! module "Electron_tools"
               call get_new_global_energy(Scell(NSC), Scell(NSC)%nrg) ! module "Electron_tools"
            case (2) ! distribution for given Ee + repulsive energy:
               call get_new_energies(Scell, matter, numpar, t, Err) ! below
               !call update_fe(Scell, matter, numpar, t, Err) ! module "Electron_tools"
            end select

            ! Get the matrix of coefficients used for calculation of forces:
            call Construct_Aij(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Aij) ! see below
            !print*, 'before dHij_s'

            if (numpar%do_atoms) then ! atoms are allowed to be moving:
               ! Construct forces:
               select type(ARRAY)
               type is (TB_H_Pettifor) ! TB parametrization according to Pettifor
                  ! Get attractive forces for atoms from the derivatives of the Hamiltonian:
                  call dHij_s(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij, M_x1, M_xrr) ! module "TB_Pettifor"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij) ! module "TB_Pettifor"
               type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni
                  ! Get attractive forces for atoms from the derivatives of the Hamiltonian:
                  call dHij_s_M(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij, M_x1, M_xrr) ! module "TB_Molteni"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_M(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij) ! module "TB_Molteni"
               type is (TB_H_Fu)
                  ! Get attractive forces for atoms from the derivatives of the Hamiltonian:
                  call dHij_s_F(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij, M_x1, M_xrr) ! module "TB_Fu"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_F(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij) ! module "TB_Fu"
               type is (TB_H_NRL)
                  ! Get the energy-weighted density matrix:
                  call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  ! Get the derivatives of the Hamiltonian:
                  call get_dHij_drij_NRL(ARRAY, Scell, NSC, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei)	! module "TB_NRL"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_NRL(ARRAY, Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_NRL"
               type is (TB_H_DFTB)
                  ! Get the energy-weighted density matrix:
                  call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  ! Get the derivatives of the Hamiltonian:
                  call get_dHij_drij_DFTB(numpar, Scell, NSC, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei)	! module "TB_DFTB"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_DFTB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_DFTB"
               type is (TB_H_3TB)
                  ! Get the energy-weighted density matrix:
                  call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  ! Get the derivatives of the Hamiltonian:
                  call get_dHij_drij_3TB(numpar, Scell, NSC, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, &
                        M_lmn, M_Aij_x_Ei, Mjs, M_Lag_exp, M_d_Lag_exp)	! module "TB_3TB"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_3TB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_3TB"
               type is (TB_H_BOP)
                  ! Get the energy-weighted density matrix:
                  call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  ! Get the derivatives of the Hamiltonian:
                  !
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  !
               type is (TB_H_xTB)
                  ! Get the energy-weighted density matrix:
                  call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  ! Get the derivatives of the Hamiltonian:
                  !call get_dHij_drij_xTB(numpar, Scell, NSC, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei)	! module "TB_xTB"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  !call Attract_TB_Forces_Press_xTB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_xTB"
               end select
            else
               Nat = size(Scell(NSC)%MDatoms)
               forall (i=1:Nat) Scell(NSC)%MDatoms(i)%forces%att(:) = 0.0d0	! no need in calculating forces if atoms are not moving anyway
            endif !(numpar%do_atoms)
         END ASSOCIATE
         !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
         !print*, 'before dErdr_s'
         
         ! Repulsive TB Hamiltonian part:
         ASSOCIATE (ARRAY2 => Scell(NSC)%TB_Repuls(:,:))
            select type(ARRAY2)
            type is (TB_Rep_Pettifor) ! TB parametrization according to Pettifor
               ! Get repulsive forces acting on all atoms:
               call dErdr_s(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s; module "TB_Pettifor"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_Pettifor"
            type is (TB_Rep_Molteni)  ! TB parametrization accroding to Molteni
               ! Get repulsive forces acting on all atoms:
               call dErdr_s_M(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s; module "TB_Molteni"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s_M(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_Molteni"
            type is (TB_Rep_Fu) ! TB parametrization according to Pettifor
               ! Get repulsive forces acting on all atoms:
               call dErdr_s_F(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s; module "TB_Fu"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s_F(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_Fu"
            type is (TB_Rep_NRL) ! TB parametrization according to Mehl (NRL)
               ! Get repulsive forces acting on all atoms:
               call dErdr_s_NRL(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s; module "TB_NRL"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s_NRL(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_NRL"
            type is (TB_Rep_DFTB) ! TB parametrization according to DFTB
               ! Get repulsive forces acting on all atoms:
               call dErdr_s_DFTB(ARRAY2, Scell, NSC) ! derivatives of the repulsive energy by s; module "TB_DFTB"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s_DFTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_DFTB"
            type is (TB_Rep_3TB) ! TB parametrization according to 3TB
               ! Get repulsive forces acting on all atoms:
               call dErdr_s_3TB(ARRAY2, Scell, NSC) ! derivatives of the repulsive energy by s; module "TB_3TB"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s_3TB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_3TB"
            type is (TB_Rep_BOP)    ! TB parametrization according to BOP
               ! Get repulsive forces acting on all atoms:
               ! Get repulsive forces acting on the supercell:
            type is (TB_Rep_xTB) ! TB parametrization according to xTB
               ! Get repulsive forces acting on all atoms:
!                call dErdr_s_xTB(ARRAY2, Scell, NSC) ! derivatives of the repulsive energy by s; module "TB_xTB"
               ! Get repulsive forces acting on the supercell:
!                call dErdr_Pressure_s_xTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_xTB"
            end select
         END ASSOCIATE
         !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         !print*, 'before vdW_forces'
         
         ! van der Waals part with TB Hamiltonian:
         call vdW_forces(Scell(NSC)%TB_Waals, Scell, NSC, numpar) ! get all the van der Waals forces, see below

         ! Coulomb potential part for modelling Coulomb explosion of a finite system:
         call Coulomb_forces(Scell(NSC)%TB_Coul, Scell, NSC, numpar) ! get Coulomb forces, see below
         
         ! Exponential wall potential part:
         call Exponential_wall_forces(Scell(NSC)%TB_Expwall, Scell, NSC, numpar) ! get Exponential wall forces, see below
         !cccccccccccccccccccccccccccccccccccccccccccccc
         !print*, 'before total_forces'


         ! Get total forces for all atoms:
         call total_forces(Scell(NSC)%MDatoms) ! sum up repulsive and attractive forces, module "Atomic_tools"
         !print*, 'before Potential_super_cell_forces'

         ! Get forces for the supercell:
         call Potential_super_cell_forces(numpar, Scell, NSC, matter)  ! module "Atomic_tools"
         !print*, 'before super_cell_forces'

         ! Adding kinetic part and pressure to super-cell forces:
         call super_cell_forces(numpar, Scell, NSC, matter, Scell(NSC)%SCforce) ! module "Atomic_tools"

         ! Get new volume of the supercell:
         call Det_3x3(Scell(NSC)%supce, Scell(NSC)%V) ! finding initial volume of the super-cell, module "Algebra_tools"
         !print*, 'cycle done'

      enddo SC
   endif DO_TB
   
   ! Get new energies of the system: potential, kinetic, global:
   call get_new_energies(Scell, matter, numpar, t, Err) ! module "TB"
   !print*, 'done get_Hamilonian_and_E'

   ! Clean up a bit:
   if (allocated(M_x1)) deallocate(M_x1)
   if (allocated(M_xrr)) deallocate(M_xrr)
   if (allocated(M_Vij)) deallocate(M_Vij)
   if (allocated(M_dVij)) deallocate(M_dVij)
   if (allocated(M_SVij)) deallocate(M_SVij)
   if (allocated(M_dSVij)) deallocate(M_dSVij)
   if (allocated(M_E0ij)) deallocate(M_E0ij)
   if (allocated(M_lmn)) deallocate(M_lmn)
   if (allocated(M_Aij_x_Ei)) deallocate(M_Aij_x_Ei)
end subroutine get_Hamilonian_and_E





subroutine get_DOS(numpar, matter, Scell, Err) ! optical coefficients, module "Optical_parameters"
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function 
   type(Solid), intent(in) :: matter     ! material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err	! error save
   !=====================================
   integer :: NSC, Nsiz, Ei_siz, i, n_types
   real(8) :: dE, Estart
   real(8), dimension(:), allocatable :: Ei_cur
   
!    print*, 'get_DOS test 0'
   
   do NSC = 1, size(Scell) ! for all super-cells
      if (numpar%save_DOS) then	! only calculate DOS if the user chose to do so:
!          call print_time_step('Before DOS:', 1.0, msec=.true.)   ! module "Little_subroutines"
         ! Set grid for DOS:
         if (.not.allocated(Scell(NSC)%DOS)) then	! it's the first time, set it:
            dE = 0.1d0	! [eV] uniform energy grid step for DOS
            Ei_siz = size(Scell(NSC)%Ei)	! the uppermost energy level at the start
            Nsiz = CEILING( (Scell(NSC)%Ei(Ei_siz) - Scell(NSC)%Ei(1) + 20.0d0*numpar%Smear_DOS)/dE )
            allocate(Scell(NSC)%DOS(2,Nsiz))
            Scell(NSC)%DOS = 0.0d0
            ! Partial DOS if needed:
            select case (numpar%DOS_splitting)
            case (1)
               n_types = size(numpar%mask_DOS,2)
               allocate(Scell(NSC)%partial_DOS(matter%N_KAO, n_types, Nsiz))
               Scell(NSC)%partial_DOS = 0.0d0
            case default
               ! No need to sort DOS per orbitals
            endselect
            
            Estart = Scell(NSC)%Ei(1) - 10.0d0*numpar%Smear_DOS	! [eV]
            !$omp PARALLEL private(i)
            !$omp do
            do i = 1, Nsiz	! set energy grid [eV]
               Scell(NSC)%DOS(1,i) = Estart + dE*dble(i-1)
            enddo
            !$omp end do
            !$omp end parallel
         endif
         
!          print*, 'get_DOS test 1'
         
         ! Now calculate the DOS:
         select case (ABS(numpar%optic_model))	! use multiple k-points, or only gamma
         case (2)	! multiple k points
            ! Partial DOS if needed:
            select case (numpar%DOS_splitting)
            case (1)
               call get_DOS_sort_complex(numpar, Scell, NSC, Scell(NSC)%DOS, numpar%Smear_DOS, Err, Scell(NSC)%partial_DOS, numpar%mask_DOS)	! see below
            case default    ! No need to sort DOS per orbitals
               call get_DOS_sort_complex(numpar, Scell, NSC, Scell(NSC)%DOS, numpar%Smear_DOS, Err)	! see below
            endselect
         case default	! gamma point
            ! Partial DOS if needed:
            
            
            select case (numpar%DOS_splitting)
            case (1)
!                print*, 'get_DOS test 2a'
               call get_DOS_sort(Scell(NSC)%Ei, Scell(NSC)%DOS, numpar%Smear_DOS, Scell(NSC)%partial_DOS, numpar%mask_DOS, Hij = Scell(NSC)%Ha)	! module "Electron_tools"
!                print*, 'get_DOS test 3a'
            case default    ! No need to sort DOS per orbitals
!                print*, 'get_DOS test 2b'
               call get_DOS_sort(Scell(NSC)%Ei, Scell(NSC)%DOS, numpar%Smear_DOS)	! module "Electron_tools"
!                print*, 'get_DOS test 3b'
            endselect
         end select
!          call print_time_step('After DOS:', 1.0, msec=.true.)   ! module "Little_subroutines"
      endif	!  (numpar%save_DOS) 
   enddo	! NSC = 1, size(Scell) ! for all super-cells
!    print*, 'get_DOS test end'
end subroutine get_DOS



subroutine get_DOS_sort_complex(numpar, Scell, NSC, DOS, smearing, Err, partial_DOS, masks_DOS)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), intent(inout) :: DOS	! [eV] grid; [a.u.] DOS
   real(8), intent(in) :: smearing	! [eV] smearing used for DOS calculations
   type(Error_handling), intent(inout) :: Err	! error save
   real(8), dimension(:,:,:), allocatable, intent(inout), optional :: partial_DOS    ! partial DOS made of each orbital type, if required to be constructed
   logical, dimension(:,:,:), intent(in), optional :: masks_DOS   ! partial DOS made of each orbital type, if required to be constructed
   !---------------------------   
   real(8), dimension(:), allocatable :: Ei	! [eV] current energy levels at the selected k point
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   complex, dimension(:,:), allocatable :: CSij	! overlap matrix of the nonorthogonal hamiltonian
   real(8), dimension(:,:), allocatable :: DOS_cur	! [eV] grid; [a.u.] DOS
   real(8), dimension(:,:,:), allocatable :: part_DOS_cur	! [a.u.] DOS
   real(8) :: kx, ky, kz
   integer :: ix, iy, iz, ixm, iym, izm, Ne, schem, Nsiz, Ngp
   real(8), pointer :: ARRAY
   logical :: do_partial
   
   Ne = size(Scell(1)%Ei)	! total number of energy levels
   do_partial = (present(partial_DOS) .and. present(masks_DOS))
   
   ! Number of k-points along each k-vector
   ixm = numpar%ixm
   iym = numpar%iym
   izm = numpar%izm
   if (allocated(numpar%k_grid)) then
      schem = 1	! user-defined grid is present
      Nsiz = size(numpar%k_grid,1)	! size of the user provided grid
   else
      schem = 0	! no user-defined grid
      Nsiz = ixm*iym*izm
   endif
   
! !$omp PARALLEL private(ix, iy, iz, kx, ky, kz, CHij, CSij, Ei, DOS_cur)
   allocate(Ei(Ne))
   allocate(CHij(Ne,Ne))
   allocate(CSij(Ne,Ne))
   allocate(DOS_cur(2,size(DOS,2)))
   DOS_cur = 0.0d0
   DOS_cur(1,:) = DOS(1,:)
   if (do_partial) then
      allocate(part_DOS_cur(size(partial_DOS,1), size(partial_DOS,2), size(partial_DOS,3)))
      part_DOS_cur = 0.0d0
   endif
   
! !$omp do reduction( + : DOS)
   do ix = 1, ixm
      do iy = 1, iym
         do iz = 1, izm
            Ngp = (ix-1)*iym*ixm + (iy-1)*ixm + iz	! number of grid point for user defined grid
            if (Ngp <= Nsiz) then
               call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid)	! below
               write(*,'(i3,i3,i3,f,f,f,a)') ix, iy, iz, kx, ky, kz, ' DOS'

               ! Construct complex Hamiltonian from the real one for the given k-point:
               call get_complex_Hamiltonian(numpar, Scell, NSC,  CHij, CSij, Ei, kx, ky, kz, Err)	! see below

               ! And get the DOS for this k point:
               if (do_partial) then
                  call get_DOS_sort(Ei, DOS_cur, smearing, part_DOS_cur, masks_DOS, CHij = CHij)	! module "Electron_tools"
                  partial_DOS(:,:,:) = partial_DOS(:,:,:) + part_DOS_cur(:,:,:)  ! save contributions of partial DOS
               else
                  call get_DOS_sort(Ei, DOS_cur, smearing)	! module "Electron_tools"
               endif
               DOS(2,:) = DOS(2,:) + DOS_cur(2,:)	! save contributions from all k points
            endif
         enddo ! iz
      enddo ! iy
   enddo ! ix
! !$omp end do
   deallocate (Ei, CHij, CSij, DOS_cur)
! !$omp END PARALLEL

   DOS(2,:) = DOS(2,:)/dble(Nsiz)
   if (do_partial) then 
      partial_DOS = partial_DOS/dble(Nsiz)
      deallocate (part_DOS_cur)
   endif
end subroutine get_DOS_sort_complex



subroutine k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, UPG)
   integer, intent(in) :: schem	! scheme for k-points samping: 0=Monkhorst-Pack; 1=user provided grid
   integer, intent(in) :: ix, iy, iz, ixm, iym, izm	! number of grid points and sizes of the grids along the reciprocal vectors
   real(8), dimension(:,:), allocatable, intent(in) :: UPG	! user provided grid
   real(8), intent(out) :: kx, ky, kz	! coordinates of the k-points in the units of the reciprocal vectors
   !------------------------------
   integer :: Ngp, Nsiz
   select case (schem)
   case (1)	! user defined grid:
      UG:if (.not.allocated(UPG)) then	! use Monkhorst-Pack grid
         kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
         ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
         kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
      else UG	! user provided the grid:
         !Ngp = ix*iy*iz	! number of grid point for user defined grid
         !Ngp = (iz-1)*iym*ixm + (iy-1)*ixm + ix	! number of grid point for user defined grid
         Ngp = (ix-1)*iym*ixm + (iy-1)*ixm + iz	! number of grid point for user defined grid
         Nsiz = size(UPG,1)	! size of the user provided grid
         if (Ngp > Nsiz) then	! If user didn't match the number of k-points to the grid provided, just use Gamma point for extra points:
            kx = 0.0d0
            ky = 0.0d0
            kz = 0.0d0
         else	! read from the user provided array:
            kx = UPG(Ngp,1)
            ky = UPG(Ngp,2)
            kz = UPG(Ngp,3)
         endif
      endif UG
   case default	! Monkhorst-Pack [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]:
      kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
      ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
      kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
   end select
end subroutine k_point_choice


subroutine get_DOS_masks(Scell, matter, numpar, only_coupling, do_cartesian)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter     ! material parameters
   type (Numerics_param), intent(inout) :: numpar    ! numerical parameters, including drude-function 
   logical, intent(in), optional :: only_coupling   ! do only coupling parameter, not the DOS?
   logical, intent(in), optional :: do_cartesian    ! cartesian basis set
   !--------------------------------
   integer :: nat, Nsiz, norb, n_types, i, NSC, j, i_orb, nkoa
   integer, pointer :: KOA
   logical :: do_DOS_masks, do_cart
   
   do_DOS_masks = .true.    ! by default, do the masks
   if (present(only_coupling)) then
      if (only_coupling) then
         do_DOS_masks = .false.    ! no need in masks
      endif
   endif
   
   NSC = 1  ! for the moment, only 1 supercell

   if (present(do_cartesian)) then
      do_cart = do_cartesian    ! user defines if use Cartesian or Spherical
   else
      do_cart = .false. ! Spherical by default
   endif
   ! Additionally, if xTB is used, it is always Cartesian:
   ASSOCIATE (ARRAY => Scell(1)%TB_Hamil(:,:))  ! only 1st supercell for now
      select type(ARRAY)
      type is (TB_H_xTB)
         do_cart = .true.
      endselect
   END ASSOCIATE

   nat = size(Scell(NSC)%MDatoms) ! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals

   BS:if (do_cart) then ! Cartesian basis set:
      ! Find number of different orbital types:
      norb = identify_xTB_orbitals_per_atom(numpar%N_basis_size) ! module "TB_xTB"
      n_types = number_of_types_of_orbitals(norb, cartesian=.true.)  ! module "Little_subroutines"

      ! Define parameters of coupling arrays:
      nkoa = size(Scell(NSC)%TB_Hamil,1)   ! kind of different atoms
      if (.not.allocated(Scell(NSC)%G_ei_partial)) allocate( Scell(NSC)%G_ei_partial(nkoa*n_types, nkoa*n_types) )

      ! Define parameters of the DOS masks:
      if (do_DOS_masks) then ! only if user requested
         if (.not.allocated(numpar%mask_DOS)) allocate(numpar%mask_DOS(matter%N_KAO, n_types, Nsiz))
         numpar%mask_DOS = .false.    ! default starting point

         do i = 1, nat
            KOA => Scell(NSC)%MDatoms(i)%KOA
            select case (numpar%N_basis_size)   ! identify basis set
            case (0) ! s
               numpar%mask_DOS(KOA, 1, i) = .true.   ! s
            case (1) ! ss*
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               numpar%mask_DOS(KOA, 2, i_orb + 2) = .true.   ! s*
            case (2) ! sp3
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
            case (3) ! sp3s*
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
               numpar%mask_DOS(KOA, 3, i_orb + 5) = .true.   ! s*
            case (4) ! sp3d6
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
               do j = 5,10
                  numpar%mask_DOS(KOA, 3, i_orb + j) = .true.   ! d
               enddo
            case (5) ! sp3d6s*
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
               do j = 5,10
                  numpar%mask_DOS(KOA, 3, i_orb + j) = .true.   ! d
               enddo
               numpar%mask_DOS(KOA, 4, i_orb + 11) = .true.   ! s*
            endselect
         enddo
      endif

   else BS ! Spherical basis set:

      ! Find number of orbitals per atom:
      norb =  Nsiz/nat ! orbitals per atom
      ! Find number of different orbital types:
      n_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"
   
      ! Define parameters of coupling arrays:
      nkoa = size(Scell(NSC)%TB_Hamil,1)   ! kind of different atoms
      if (.not.allocated(Scell(NSC)%G_ei_partial)) allocate( Scell(NSC)%G_ei_partial(nkoa*n_types, nkoa*n_types) )
   
      ! Define parameters of the DOS masks:
      if (do_DOS_masks) then ! only if user requested
         if (.not.allocated(numpar%mask_DOS)) allocate(numpar%mask_DOS(matter%N_KAO, n_types, Nsiz))
         numpar%mask_DOS = .false.    ! default starting point
   
         do i = 1, nat
            KOA => Scell(NSC)%MDatoms(i)%KOA  ! Correct order, checked by cohesive energy minimum
            select case (norb)   ! identify basis set
            case (1) ! s
               numpar%mask_DOS(KOA, 1, i) = .true.   ! s
            case (2) ! ss*
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               numpar%mask_DOS(KOA, 2, i_orb + 2) = .true.   ! s*
            case (4) ! sp3
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
            case (5) ! sp3s*
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
               numpar%mask_DOS(KOA, 3, i_orb + 5) = .true.   ! s*
            case (9) ! sp3d5
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
               do j = 5,9
                  numpar%mask_DOS(KOA, 3, i_orb + j) = .true.   ! d
               enddo
            case (10) ! sp3d5s*
               i_orb = (i-1)*norb
               numpar%mask_DOS(KOA, 1, i_orb + 1) = .true.   ! s
               do j = 2,4
                  numpar%mask_DOS(KOA, 2, i_orb + j) = .true.   ! p
               enddo
               do j = 5,9
                  numpar%mask_DOS(KOA, 3, i_orb + j) = .true.   ! d
               enddo
               numpar%mask_DOS(KOA, 4, i_orb + 10) = .true.   ! s*
            endselect
         enddo
      endif
   endif BS

   nullify(KOA)
end subroutine get_DOS_masks


 subroutine get_DOS_weights(ind, masks_DOS, DOS_weights, Hij, CHij)
   integer, intent(in) :: ind   ! which model to use
   logical, dimension(:,:,:), intent(in) :: masks_DOS   ! partial DOS made of each orbital type, if required to be constructed
   real(8), dimension(:,:,:), allocatable, intent(inout) :: DOS_weights     ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:,:), intent(in), optional :: Hij      ! real eigenvectors
   complex, dimension(:,:), intent(in), optional :: CHij ! complex eigenvectors
   !-------------------------------
   real(8), dimension(size(masks_DOS,3)) :: temp_vec
   real(8) :: temp
   integer :: j, Nsiz, N_at, N_types, i_at, i_types
   if (ind == 1) then ! get weighted DOS:
      N_at = size(masks_DOS,1)
      N_types = size(masks_DOS,2)
      Nsiz = size(masks_DOS,3)
      if (.not.allocated(DOS_weights)) allocate(DOS_weights(N_at,N_types,Nsiz))

      !$omp PARALLEL private(j, temp_vec, temp, i_at, i_types)
      !$omp do
      do j = 1, Nsiz	! for all energy levels
         if (present(Hij)) then ! real H
            temp_vec(:) = Hij(:,j) * Hij(:,j)
         elseif (present(CHij)) then    ! complex H
            temp_vec(:) = dble(conjg(CHij(:,j)) * CHij(:,j))
         else   ! undefined
            temp_vec(:) = 0.0d0
         endif
         temp = SUM(temp_vec)
         do i_at = 1, N_at
            do i_types = 1, N_types
               DOS_weights(i_at, i_types,j) = SUM(temp_vec(:), MASK = masks_DOS(i_at, i_types, :))/temp
            enddo ! i_types
        enddo ! i_at
   !       write(*,'(i4,f,f,f)') j, DOS_weights(:,:,j)
      enddo !  j = 1, Nsiz
      !$omp end do
      !$omp end parallel
   endif ! (ind == 1)
end subroutine get_DOS_weights


subroutine get_Mulliken(Mulliken_model, masks_DOS, DOS_weights, Hij, fe, MDatoms, mulliken_Ne)
   integer, intent(in) :: Mulliken_model   ! which model to use
   logical, dimension(:,:,:), intent(in) :: masks_DOS   ! partial DOS made of each orbital type, if required to be constructed
   real(8), dimension(:,:,:), intent(inout) :: DOS_weights     ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:,:), intent(in) :: Hij      ! real eigenvectors
   real(8), dimension(:), intent(in) :: fe    ! electron distribution
   type(Atom), dimension(:), intent(in) :: MDAtoms
   real(8), dimension(:), intent(out) :: mulliken_Ne   ! Mulliken charges
   !-------------------------------
   real(8) :: temp
   integer :: j, Nsiz, N_at, N_types, i_at, i_types, Nat
   if (Mulliken_model >= 1) then ! get Mulliken charges
      N_at = size(masks_DOS,1)
      N_types = size(masks_DOS,2)
      ! Check atomic charges:
      mulliken_Ne(:) = 0.0d0  ! to start from
      !$omp PARALLEL private(i_at, i_types, j)
      !$omp do
      do i_at = 1, N_at
         do i_types = 1, N_types
            mulliken_Ne(i_at) = mulliken_Ne(i_at) + SUM(fe(:) * DOS_weights(i_at, i_types,:))
         enddo   ! i_types
      enddo ! i_at
      !$omp end do
      !$omp end parallel
      ! Normalize to the number of atoms:
      do i_at = 1, N_at
         ! how many atoms of this kind are in the supercell:
         Nat = COUNT(MASK = (MDatoms(:)%KOA == i_at))
         mulliken_Ne(i_at) = mulliken_Ne(i_at) / dble(Nat)
      enddo
   endif ! (Mulliken_model >= 1) 
end subroutine get_Mulliken





subroutine get_complex_Hamiltonian(numpar, Scell, NSC,  CHij, CSij, Ei, kx, ky, kz, Err)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(inout), allocatable :: Ei	! [eV] current energy levels at the selected k point
   complex, dimension(:,:), intent(inout), allocatable :: CHij	! eigenvectors of the hamiltonian
   complex, dimension(:,:), intent(inout), allocatable :: CSij	! overlap matrix of the nonorthogonal hamiltonian
   real(8), intent(in) :: kx, ky, kz	! k point
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------
   ! Construct complex Hamiltonian from the real one for the given k-point:
   if ((abs(kx) < 1.0d-14) .AND. (abs(ky) < 1.0d-14) .AND. (abs(kz) < 1.0d-14)) then ! Gamma point:
      Ei = Scell(NSC)%Ei    !already known
   else ! any other point:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
         select type(ARRAY)
            type is (TB_H_Pettifor)
               call Complex_Hamil_tot(numpar, Scell, NSC, Scell(NSC)%MDAtoms, ARRAY, CHij=CHij, ksx=kx, ksy=ky, ksz=kz) ! "TB_Pettifor"
               call sym_diagonalize(CHij, Ei, Err%Err_descript) ! modeule "Algebra_tools"
            type is (TB_H_Molteni)
               call Complex_Hamil_tot_Molteni(numpar, Scell, NSC, Scell(NSC)%MDAtoms, ARRAY, CHij=CHij, ksx=kx, ksy=ky, ksz=kz) ! "TB_Molteni"
               call sym_diagonalize(CHij, Ei, Err%Err_descript) ! modeule "Algebra_tools"
            type is (TB_H_Fu)
               call Complex_Hamil_tot_F(numpar, Scell, NSC, Scell(NSC)%MDAtoms, ARRAY, CHij=CHij, ksx=kx, ksy=ky, ksz=kz) ! "TB_Fu"
               call sym_diagonalize(CHij, Ei, Err%Err_descript) ! modeule "Algebra_tools"
            type is (TB_H_NRL)
               call Complex_Hamil_NRL(numpar, Scell, NSC, CHij, CSij, Ei, kx, ky, kz, Err) ! "TB_NRL"
            type is (TB_H_DFTB)
               call Complex_Hamil_DFTB(numpar, Scell, NSC, CHij, CSij, Ei, kx, ky, kz, Err) ! "TB_DFTB"   
            type is (TB_H_3TB)
               call Complex_Hamil_DFTB(numpar, Scell, NSC, CHij, CSij, Ei, kx, ky, kz, Err) ! "TB_DFTB", the are the same
            type is (TB_H_xTB)
!                call Complex_Hamil_xTB(numpar, Scell, NSC, CHij, CSij, Ei, kx, ky, kz, Err) ! "TB_xTB"
         end select	!  type of Hamiltonain parameterization
      END ASSOCIATE
   endif
end subroutine get_complex_Hamiltonian


subroutine construct_complex_Hamiltonian(numpar, Scell, NSC, H_non, CHij, Ei, ksx, ksy, ksz, Err, cPRRx, cPRRy, cPRRz, Sij, CSij)   ! CHECKED
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), intent(in) :: H_non	! Non-diagonalized real hamiltonian, must be provided
   complex, dimension(:,:), intent(out), allocatable :: CHij	! complex the hamiltonian, to be constructed
   real(8), dimension(:), intent(inout), allocatable :: Ei	! [eV] current energy levels at the selected k point
   real(8), intent(in) ::  ksx, ksy, ksz	! k point [relative coordinates]
   type(Error_handling), intent(inout) :: Err	! error save   
   complex, dimension(:,:), intent(out), allocatable, optional :: cPRRx, cPRRy, cPRRz ! complex momentum operators
   real(8), dimension(:,:), intent(in), optional :: Sij	   ! real overlap matrix of the nonorthogonal hamiltonian, must be provided for nonorthogonal case
   complex, dimension(:,:), intent(out), allocatable, optional :: CSij  ! complex overlap matrix of the nonorthogonal hamiltonian, to be created
   !-------------------------------------------------------------------------------
   complex, dimension(:,:), allocatable :: CHij_temp, CHij_non, CSij_save, CHij_orth
   integer :: Nsiz, j, nat, m, atom_2, i, j1, l, i1, k, norb
   real(8), dimension(3,3) :: k_supce   ! reciprocal supercell vectors
   real(8) :: temp, kx, ky, kz
   real(8), target :: nol
   real(8), pointer :: x1, y1, z1
   complex(8) :: expfac, SH_1
   character(200) :: Error_descript
   nol = 0.0d0
   Error_descript = ''
   nat = size(Scell(NSC)%MDatoms) ! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
   norb =  Nsiz/nat ! orbitals per atom
   
   ! Allocate complex parameters:
   if (.not.allocated(CHij)) allocate(CHij(Nsiz,Nsiz))
   CHij = dcmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij_temp)) allocate(CHij_temp(Nsiz,Nsiz))
   if (.not.allocated(CHij_non)) allocate(CHij_non(Nsiz,Nsiz))
   CHij_non = dcmplx(0.0d0,0.0d0)	! to start with
   CHij_temp = dcmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij_orth)) allocate(CHij_orth(Nsiz,Nsiz))  ! orthogonalized Hamiltonian
   CHij_orth = dcmplx(0.0d0,0.0d0)	! to start with
   
   if (present(Sij)) then   ! it is nonorthogonal case:
      if (.not.allocated(CSij)) allocate(CSij(Nsiz,Nsiz))
      if (.not.allocated(CSij_save)) allocate(CSij_save(Nsiz,Nsiz))
      CSij = dcmplx(0.0d0,0.0d0)	! to start with
      CSij_save = dcmplx(0.0d0,0.0d0)	! to start with
   endif
   
   call Reciproc(Scell(NSC)%supce, k_supce) ! create reciprocal super-cell, module "Algebra_tools"
   call Convert_reciproc_rel_to_abs(ksx, ksy, ksz, k_supce, kx, ky, kz) ! get absolute k-values, molue "Atomic_tools"

   ! 1) Construct complex Hamiltonian and overlap:
   !$omp parallel
   !$omp do private(j, m, atom_2, i, x1, y1, z1, expfac, j1, l, i1, k)
   do j = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one  
         if (atom_2 == 0) then
            i = j
            x1 => nol
            y1 => nol
            z1 => nol
         else
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1)	! at this distance, X
            y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2)	! at this distance, Y
            z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3)	! at this distance, Z
         endif ! (atom_2 .EQ. 0)
         
         if ((abs(kx) < 1.0d-14) .AND. (abs(ky) < 1.0d-14) .AND. (abs(kz) < 1.0d-14)) then
            expfac = dcmplx(1.0d0,0.0d0)
         else
            expfac = exp(g_CI*dcmplx(kx*x1 + ky*y1 + kz*z1,0.0d0))
         endif

         do j1 = 1,norb ! all orbitals
            l = (j-1)*norb+j1
            do i1 = 1,norb ! all orbitals
               k = (i-1)*norb+i1

               CHij_temp(k,l) = DCMPLX(H_non(k,l),0.0d0)*expfac
               if ((isnan(real(CHij_temp(k,l)))) .OR. isnan(aimag(CHij_temp(k,l)))) then
                  Error_descript = 'CHij_temp ISNAN in construct_complex_Hamiltonian'
                  call Save_error_details(Err, 6, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  print*, i, j, k, l, CHij_temp(k,l)
               endif
               
               if (present(Sij)) then
                  CSij(k,l) = DCMPLX(Sij(k,l),0.0d0)*expfac
                  if ((isnan(real(CSij(k,l)))) .OR. isnan(aimag(CSij(k,l)))) then
                     Error_descript = 'CHij_temp ISNAN in construct_complex_Hamiltonian'
                     call Save_error_details(Err, 6, Error_descript)
                     print*, trim(adjustl(Error_descript))
                     print*, i, j, k, l, CSij(k,l)
                  endif
               endif
            enddo ! i1
         enddo ! j1
      enddo ! atom_2
   enddo ! j
   !$omp end do 
   !$omp end parallel

   ! Temporarily save nonorthogonal Hamiltonian and overlap matrix:
   CHij_non = CHij_temp
   if (present(Sij)) then 
      CSij_save = CSij
      call diagonalize_complex_Hamiltonian(CHij_temp, Ei, Err, CSij, CHij_orth)    ! below
   else
      call diagonalize_complex_Hamiltonian(CHij_temp, Ei, Err)    ! below
   endif
   CHij = CHij_temp ! save for output
   
   ! Momentum operators:
   if (numpar%optic_model .EQ. 2) then ! create matrix element:
      if (.not.allocated(cPRRx)) allocate(cPRRx(Nsiz,Nsiz))
      if (.not.allocated(cPRRy)) allocate(cPRRy(Nsiz,Nsiz))
      if (.not.allocated(cPRRz)) allocate(cPRRz(Nsiz,Nsiz))
      cPRRx = dcmplx(0.0d0,0.0d0)	! to start with
      cPRRy = dcmplx(0.0d0,0.0d0)	! to start with
      cPRRz = dcmplx(0.0d0,0.0d0)	! to start with
      ! Optical matrix elements for non-orthogonal TB are taken from:
      ! arXiv:1805.08918v1 -- https://128.84.21.199/abs/1805.08918
      !$omp parallel
      !$omp do private(j, m, atom_2, i, x1, y1, z1, j1, l, i1, k, SH_1)
      do j = 1,nat	! all atoms
         m = Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one
            if (atom_2 == 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
               x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1)	! at this distance, X
               y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2)	! at this distance, Y
               z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3)	! at this distance, Z
            endif ! (atom_2 .EQ. 0)
         
            if (i > 0) then
               do j1 = 1,norb	! all orbitals for sp3d5
                  l = (j-1)*norb+j1
                  do i1 = 1,norb	! all orbitals for sp3d5
                     k = (i-1)*norb+i1
                     if (i == j) then ! contribution of the same atom, according to Trani:

                        if (j1 == i1) then
                           ! Skip the same orbital, no overlap
                        else

                           ! [1] Nonorthogonal expression:
                           !SH_1 = DCMPLX(0.27d0,0.0d0) * (CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l))
                           SH_1 = DCMPLX(1.5d0,0.0d0) * (CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l))

                           ! Testing alternative orthogonalized Hamiltonian:
                           !SH_1 = DCMPLX(0.27d0,0.0d0) * CHij_orth(k,l)   ! wrong

                           cPRRx(k,l) = SH_1
                           cPRRy(k,l) = SH_1
                           cPRRz(k,l) = SH_1
                        endif

                     else	! different atoms at distance {x,y,z}:
                        SH_1 = CHij_non(k,l)
!                         SH_1 = CHij_non(l,k)    ! testing
                        if (present(Sij)) then ! nonorthogonal Hamiltonian:
                           ! [1] Nonorthogonal expression:
                           SH_1 = SH_1 - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l)  ! Correct
                           ! Testing alternative orthogonalized Hamiltonian:
                           !SH_1 = CHij_orth(k,l)  ! wrong
                        endif

                        cPRRx(k,l) = DCMPLX(x1,0.0d0)*SH_1
                        cPRRy(k,l) = DCMPLX(y1,0.0d0)*SH_1
                        cPRRz(k,l) = DCMPLX(z1,0.0d0)*SH_1
                     endif
                     if (dble(cPRRx(k,l)) .GT. 1d10) write(*,'(i5,i5,es,es, es,es, es,es, es, es)') i, j, cPRRx(k,l),  CHij_non(k,l), CSij_save(k,l),  Ei(k), x1
                     if (dble(cPRRy(k,l)) .GT. 1d10) print*, i, j, cPRRy(k,l)
                     if (dble(cPRRz(k,l)) .GT. 1d10) print*, i, j, cPRRz(k,l)
                  enddo ! i1
               enddo ! j1
            endif ! (i > 0)
         enddo ! atom_2
      enddo ! j
      !$omp end do 
      !$omp end parallel
      
      ! Normalization factors:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
      select type(ARRAY)
      type is (TB_H_Pettifor)
         temp = 1.0d0
      type is (TB_H_Molteni)
         temp = 1.0d0
      type is (TB_H_Fu)
         temp = 1.0d0
      type is (TB_H_NRL)
         temp = 1.0d0   ! 0.5d0 for testing
      type is (TB_H_DFTB)
         temp = 1.0d0
      type is (TB_H_3TB)
         temp = 1.0d0
      type is (TB_H_xTB)
         temp = 1.0d0
      end select
      END ASSOCIATE
      cPRRx = cPRRx * temp
      cPRRy = cPRRy * temp
      cPRRz = cPRRz * temp
   endif
   
   deallocate(CHij_temp, CHij_non)
   if (allocated(CSij_save)) deallocate(CSij_save)
   nullify(x1, y1, z1)
end subroutine construct_complex_Hamiltonian



subroutine diagonalize_complex_Hamiltonian(CHij, Ei, Err, CSij, CHij_orth)
   complex, dimension(:,:), intent(inout) :: CHij	! complex hermitian Hamiltonian
   real(8), dimension(:), intent(out), allocatable :: Ei ! eigenvalues [eV]
   type(Error_handling), intent(inout) :: Err	! error save   
   complex, dimension(:,:), intent(inout), optional ::  CSij    ! overlap matrix, for nonorthogonal Hamiltonain
   complex, dimension(:,:), intent(inout), optional ::  CHij_orth ! orthogonalized Hamiltonian, if requested
   !----------------------------
   complex, dimension(size(CHij,1), size(CHij,2)) :: CHij_temp
   integer :: Nsiz, j
   character(200) :: Error_descript
   Error_descript = ''  ! to start with, no error
   Nsiz = size(CHij,1)
   if (.not.allocated(Ei)) allocate(Ei(Nsiz))
   ORTH: if (.not.present(CSij)) then ! orthogonal:
      ! Direct diagonalization:
      call sym_diagonalize(CHij, Ei, Err%Err_descript) ! modeule "Algebra_tools"
   else ORTH ! nonorthogonal
      ! Solve linear eigenproblem:
      ! 1) Orthogonalize the Hamiltonian using Loewdin procidure:
      ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
      CHij_temp = CHij
      call Loewdin_Orthogonalization_c(Nsiz, CSij, CHij_temp, Err)	! module "TB_NRL"

      if (present(CHij_orth)) then  ! Save orthogonalized Hamiltonian (for optical coefficients below)
         CHij_orth = CHij_temp
      endif

      ! 2) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
      call sym_diagonalize(CHij_temp, Ei, Error_descript)   ! module "Algebra_tools"
      if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
         Error_descript = 'diagonalize_complex_Hamiltonian: '//trim(adjustl(Error_descript))
         call Save_error_details(Err, 6, Error_descript)    ! module "Objects"
         print*, trim(adjustl(Error_descript))
      endif
      ! 3) Convert the eigenvectors back into the non-orthogonal basis:
      call mkl_matrix_mult('N', 'N', CSij, CHij_temp, CHij)	! module "Algebra_tools"
      
!       ! 4) We need to renormalize the wave functions, as they are not normalized to 1 after this procidure:
!       do j = 1, Nsiz
!         CHij(:,j) = CHij(:,j) / SQRT(SUM( dconjg(CHij(:,j)) * CHij(:,j) ))
!       enddo

!        do j = 1, size(Ei)
!           write(*,'(i5,es,es,es,es,es)') j, Ei(j), CHij_temp(j,1), Scell(NSC)%Ei(j)
!        enddo
!        PAUSE 'Ei pause'
   endif ORTH
end subroutine diagonalize_complex_Hamiltonian



subroutine Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn)
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_x1	! matrix of coefficients x1
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_xrr	! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_lmn	! matrix of direction cosines l, m, n
   !---------------------------
   real(8), pointer :: x, y, z, r
   real(8) :: r2
   integer i, nat, m, j, atom_2
   
   nat = Scell(NSC)%Na ! number of atoms
   if (.not.allocated(M_x1)) allocate(M_x1(3,nat,nat))
   if (.not.allocated(M_xrr)) allocate(M_xrr(3,nat,nat))
   if (.not.allocated(M_lmn)) allocate(M_lmn(3,nat,nat))
   !$OMP WORKSHARE
   M_x1 = 0.0d0
   M_xrr = 0.0d0
   M_lmn = 0.0d0
   !$OMP END WORKSHARE
   
!$omp PARALLEL
!$omp do schedule(dynamic) private(i, m, atom_2, j, x, y, z, r)
   do i = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one  
         j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
!          if (j > nat) then
!             print*, 'ATOM j', i, j, atom_2
!             print*, 'HUGE!', Scell(NSC)%Near_neighbor_list(i,:)
!             pause 'CHECK-PAUSE'
!          endif
         
         if (j .GT. 0) then
            x => Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
            y => Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
            z => Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
            r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
!             r2 = r*r

            M_x1(1,i,j) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
            M_x1(2,i,j) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
            M_x1(3,i,j) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
!             M_x1(1,i,j) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
!             M_x1(2,i,j) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
!             M_x1(3,i,j) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)
            
            M_lmn(1,i,j) = x/r	! l
            M_lmn(2,i,j) = y/r	! m
            M_lmn(3,i,j) = z/r	! n
!             write(*, '(a,i2,i2,f,f,f)' ) 'cos: ', i, j, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j)
            
            M_xrr(1,i,j) = M_lmn(1,i,j)/r	!xrr
            M_xrr(2,i,j) = M_lmn(2,i,j)/r	!yrr
            M_xrr(3,i,j) = M_lmn(3,i,j)/r	!zrr
            
!             M_xrr(1,i,j) = x/r2 !xrr
!             M_xrr(2,i,j) = y/r2 !yrr
!             M_xrr(3,i,j) = z/r2 !zrr
         endif
      enddo
   enddo
!$omp end do
!$omp END PARALLEL

   nullify(x, y, z, r)
end subroutine Construct_M_x1


subroutine Exponential_wall_forces(TB_Expwall, Scell, NSC, numpar) ! get Exponential wall forces
   class(TB_Exp_wall), allocatable, dimension(:,:), intent(in) :: TB_Expwall	! exponential wall
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of super-cell
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !----------------------------------
    if (allocated(TB_Expwall)) then ! if we have vdW potential defined
      select type (TB_Expwall)
      type is (TB_Exp_wall_simple)
          ! Forces for all atoms:
         call d_Exp_wall_pot_s(Scell, NSC, TB_Expwall, numpar)	! module "Exponential_wall"
         ! Forces for the super-cell:
         call d_Exp_wall_Pressure_s(Scell, NSC, TB_Expwall, numpar)	! module "Exponential_wall"
      end select
   else
      ! No additional forces, if there is no exponential wall parameters given
   endif
end subroutine Exponential_wall_forces



subroutine Coulomb_forces(TB_Coul, Scell, NSC, numpar)
   class(TB_Coulomb), allocatable, dimension(:,:), intent(in) :: TB_Coul ! van der Waals within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of super-cell
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !----------------------------------
   real(8), dimension(:,:,:), allocatable :: Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, XijSupce, YijSupce, ZijSupce
   integer :: Nx, Ny, Nz
   
   if (allocated(TB_Coul)) then ! if we have vdW potential defined
      select type (TB_Coul)
      type is (TB_Coulomb_cut) ! so far, it is the only type we have
         ! Get multipliers used many times into temporary arrays:
         call Construct_B_C(TB_Coul, Scell, NSC, Scell(NSC)%MDatoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Coulomb"
         ! Forces for all atoms:
!          call dCoulombdr_s(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Nx, Ny, Nz) ! module "Coulomb"
         call d_Forces_s(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Nx, Ny, Nz) ! below
         ! Forces for the super-cell:
!          call dCoulombdr_Pressure_s(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz)  ! module "Coulomb"
         call d_Forces_Pressure(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! below

         if (allocated(Bij))   deallocate(Bij)
         if (allocated(A_rij)) deallocate(A_rij)
         if (allocated(Xij))   deallocate(Xij)
         if (allocated(Yij))   deallocate(Yij)
         if (allocated(Zij))   deallocate(Zij)
         if (allocated(SXij))  deallocate(SXij)
         if (allocated(SYij))  deallocate(SYij)
         if (allocated(SZij))  deallocate(SZij)
         if (allocated(XijSupce))   deallocate(XijSupce)
         if (allocated(YijSupce))   deallocate(YijSupce)
         if (allocated(ZijSupce))   deallocate(ZijSupce)
      end select
   else
      ! No additional forces, if there is no van der Waals parameters given
   endif
end subroutine Coulomb_forces


subroutine vdW_forces(TB_Waals, Scell, NSC, numpar)
   class(TB_vdW), allocatable, dimension(:,:), intent(in) :: TB_Waals ! van der Waals within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of super-cell
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !----------------------------------
   real(8), dimension(:,:,:), allocatable :: Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, XijSupce, YijSupce, ZijSupce
   integer :: Nx, Ny, Nz

   if (allocated(TB_Waals)) then ! if we have vdW potential defined
      select type (TB_Waals)
      type is (TB_vdW_Girifalco) ! so far, it is the only type we have
!          call dvdWdr_s_OLD(TB_Waals, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! module "Van_der_Waals"
!          call dvdWdr_Pressure_s_OLD(TB_Waals, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! module "Van_der_Waals"

         ! Get multipliers used many times into temporary arrays:
         call Construct_B(TB_Waals, Scell, NSC, numpar, Scell(NSC)%MDatoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Van_der_Waals"
         
         ! Forces for all atoms:
!          call dvdWdr_s(TB_Waals, Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Nx, Ny, Nz) ! module "Van_der_Waals"
         call d_Forces_s(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Nx, Ny, Nz) ! below
         
         ! Forces for the super-cell:
!          call dvdWdr_Pressure_s(TB_Waals, Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Van_der_Waals"
         call  d_Forces_Pressure(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! below

         if (allocated(Bij))   deallocate(Bij)
         if (allocated(A_rij)) deallocate(A_rij)
         if (allocated(Xij))   deallocate(Xij)
         if (allocated(Yij))   deallocate(Yij)
         if (allocated(Zij))   deallocate(Zij)
         if (allocated(SXij))  deallocate(SXij)
         if (allocated(SYij))  deallocate(SYij)
         if (allocated(SZij))  deallocate(SZij)
         if (allocated(XijSupce))   deallocate(XijSupce)
         if (allocated(YijSupce))   deallocate(YijSupce)
         if (allocated(ZijSupce))   deallocate(ZijSupce)
      end select
   else
      ! No additional forces, if there is no van der Waals parameters given
   endif
end subroutine vdW_forces




! Derivatives of the vdW energy by s:
subroutine d_Forces_s(atoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, Nx, Ny, Nz) 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: Bij, A_rij, Xij, Yij, Zij ! intermediate calculations matrices
   integer, intent(in) :: Nx, Ny, Nz ! number of super-cells to go through
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   REAL(8), DIMENSION(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, x, y, z, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, j1, ian, dik, djk, n, m, atom_2
   integer :: x_cell, y_cell, z_cell, coun_cell
   integer, DIMENSION(3) :: zb
   real(8), dimension(:,:), allocatable :: Erx_s
   logical :: origin_cell
   n = size(atoms) ! total number of atoms
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

   !$omp PARALLEL private(i1,j1,ian,dik,djk,dpsi,x1,b,a_r,ddlta,b_delta, x_cell, y_cell, z_cell, coun_cell, zb, origin_cell)
   !$omp DO
   do ian = 1, n  ! Forces for all atoms
      !Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with   
      do i1 = 1, n ! contribution from all atoms
         if (ian == i1) then
            dik = 1
         else
            dik = 0
         endif
         dpsi = 0.0d0
         do j1 = 1,n ! for each pair of atoms
            if (ian == j1) then
               djk = 1
            else
               djk = 0
            endif
            cos_if:if ((dik-djk) /= 0) then ! without it, it gives ERROR
               XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
                  YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
                     ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
                        zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
                        origin_cell = ALL(zb==0) ! if it is the origin cell
                        cell_if:if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                           ! contribution from this image of the cell:
                           coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"
                     
                           x1(1) = Xij(coun_cell,i1,j1) ! x*supce(1,1) + y*supce(1,2) + z*supce(1,3)
                           x1(2) = Yij(coun_cell,i1,j1) ! x*supce(2,1) + y*supce(2,2) + z*supce(2,3)
                           x1(3) = Zij(coun_cell,i1,j1) ! x*supce(3,1) + y*supce(3,2) + z*supce(3,3)
                           b = Bij(coun_cell,i1,j1) ! dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))
                           a_r = A_rij(coun_cell,i1,j1) ! call shortest_distance(Scell, NSC, atoms, i1, j1, A_rij(coun_cell,i1,j1), x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz) 

                           ddlta = real(dik - djk)/a_r
                           b_delta = b*ddlta
                           dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                           dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                           dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
                        endif cell_if
                     enddo ZC2
                  enddo YC2
               enddo XC2
            endif cos_if
         enddo ! j1

         Erx_s(1,ian) = Erx_s(1,ian) + dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      ! Add van der Waals force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)
!       write(*,'(a,i3,es,es,es)') 'NEW', ian, Erx_s(:,ian)
   enddo ! ian
   !$omp end do
   !$omp end parallel

   deallocate(Erx_s)
END subroutine d_Forces_s


! Derivatives of the vdW energy by h:
subroutine d_Forces_Pressure(atoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij ! intermediate calculations matrices
   integer, intent(in) :: Nx, Ny, Nz ! number of super-cells to go through
   !===============================================
   REAL(8), DIMENSION(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i,j,k,l,n
   integer x_cell, y_cell, z_cell, zb(3), coun_cell
   real(8) r, rcur(3), scur(3), PForce(3,3)
   real(8) df_psy, psi, dpsy 
   logical :: origin_cell

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = size(atoms)

      PForce = 0.0d0 ! to start with
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         do j = 1, n ! do for all pairs of atoms
            XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
               YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
                  ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
                     zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
                     origin_cell = ALL(zb==0) ! if it is the origin cell
                     cell_if:if ((j /= i) .or. (.not.origin_cell)) then ! exclude self-interaction within original super cell
                        ! contribution from this image of the cell:
                        coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"

                        rcur(1) = Xij(coun_cell,i,j) ! x
                        rcur(2) = Yij(coun_cell,i,j) ! y
                        rcur(3) = Zij(coun_cell,i,j) ! z
                        scur(1) = SXij(coun_cell,i,j) ! SX
                        scur(2) = SYij(coun_cell,i,j) ! SY
                        scur(3) = SZij(coun_cell,i,j) ! SZ
                        r = A_rij(coun_cell,i,j)     ! R
                        dpsy = Bij(coun_cell,i,j)    ! dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))

                        do k = 1,3 ! supce indices: a,b,c
                           do l = 1,3  ! supce indices: x,y,z
                              Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                           enddo ! l
                        enddo ! k
                     endif cell_if
                  enddo ZC2
               enddo YC2
            enddo XC2
         enddo ! j

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               !Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + Rep_Pr(l,k) !*0.5d0
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)
            enddo ! l
         enddo ! k
      enddo ! i
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add vdW part to existing TB part
!       print*, 'NEW P', PForce
   endif
end subroutine d_Forces_Pressure


!-------------------------------------

subroutine change_r_cut_TB_Hamiltonian(d0, TB_Repuls, TB_Hamil, TB_Waals) ! used only for testing of total energy vs super-cell size
   real(8), intent(in) :: d0
   class(TB_repulsive), dimension(:,:), intent(inout), optional :: TB_Repuls  ! parameters of the repulsive part of TB
   class(TB_Hamiltonian), dimension(:,:), intent(inout), optional ::  TB_Hamil ! parameters of the
   class(TB_vdW), allocatable, dimension(:,:), intent(inout), optional :: TB_Waals ! van der Waals within TB
   integer :: i, j

   if (present(TB_Hamil)) then
     select type (TB_Hamil)
      type is (TB_H_Pettifor)
         do i = 1,size(TB_Hamil,1)
            TB_Hamil(i,:)%rm = TB_Hamil(i,:)%rm*(1.0d0 + 1.30d0*d0/TB_Hamil(i,:)%r1)
            TB_Hamil(i,:)%r1 = 1.30d0*d0
         enddo
      type is (TB_H_Molteni)
         do i = 1,size(TB_Hamil,1)
            TB_Hamil(i,:)%rcut = 1.30d0*d0
         enddo
      type is (TB_H_Fu)
         do i = 1,size(TB_Hamil,1)
            TB_Hamil(i,:)%rm = TB_Hamil(i,:)%rm*(1.0d0 + 1.30d0*d0/TB_Hamil(i,:)%r1)
            TB_Hamil(i,:)%r1 = 1.30d0*d0
         enddo
     end select
   endif
   
   if (present(TB_Repuls)) then
     select type (TB_Repuls)
       type is (TB_Rep_Pettifor)
         do i = 1,size(TB_Hamil,1)
            TB_Repuls(i,:)%dm = TB_Repuls(i,:)%dm*(1.0d0 + 1.30d0*d0/TB_Repuls(i,:)%d1)
            TB_Repuls(i,:)%d1 = 1.30d0*d0
         enddo
       type is (TB_Rep_Molteni)
         do i = 1,size(TB_Repuls,1)
            TB_Repuls(i,:)%rcut = 1.30d0*d0
         enddo
      type is (TB_Rep_Fu)
         do i = 1,size(TB_Hamil,1)
            TB_Repuls(i,:)%dm = TB_Repuls(i,:)%dm*(1.0d0 + 1.30d0*d0/TB_Repuls(i,:)%d1)
            TB_Repuls(i,:)%d1 = 1.30d0*d0
         enddo
     endselect
   endif
   
   if (present(TB_Waals)) then
     if (allocated(TB_Waals)) then ! if we have vdW potential defined
      select type (TB_Waals)
      type is (TB_vdW_Girifalco)
         do i = 1,size(TB_Waals,1)
            TB_Waals(i,:)%r_S = TB_Waals(i,:)%r_S*(1.30d0*d0/TB_Waals(i,:)%r_L)
            TB_Waals(i,:)%r_L =  1.30d0*d0
         enddo
      end select
     endif
   endif
end subroutine change_r_cut_TB_Hamiltonian


! Repulsive part:
subroutine get_pot_nrg(Scell, matter, numpar)	! Repulsive potential energy
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !type(Error_handling), intent(inout) :: Err	! error save
   !========================================================
   real(8) Erepuls, Pot_phi
   integer NSC
   ! Calculations of the contribution of electrons into potential energy of atoms,
   ! first term in Eq.(2.44) from H.Jeschke PhD Thesis, Page 41
   !call set_total_el_energy(Ei,fe,Eelectr) ! module "Electron_tools"
   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
   
      call get_low_e_energy(Scell, matter) ! module "Electron_tools"
      !print*, 'get_pot_nrg 0'

      ! Calculations of the contribution of repulsive energies into potential energy of atoms,
      ! second term in Eq.(2.44) from H.Jeschke PhD Thesis, Page 41
      do NSC = 1, size(Scell) ! for all supercells
         
         Erepuls = Erep_s(Scell(NSC)%TB_Repuls(:,:), Scell, NSC, numpar) ! sum over all atom pairs
         !print*, 'get_pot_nrg 1'
         
         !Pot_phi = Scell(NSC)%nrg%El_tot + Erepuls + TB_Repuls%E0_TB*Scell(NSC)%Na ! [eV] potential energy, Eq.(2.44) in H.Jeschke Thesis, Page 41
         Pot_phi = Scell(NSC)%nrg%El_low + Erepuls ! [eV] potential energy, Eq.(2.44) in H.Jeschke Thesis, Page 41
         
         Scell(NSC)%nrg%At_pot = Pot_phi/real(Scell(NSC)%Na) ! [eV/atom]
         
         !Scell(NSC)%nrg%E_rep = Erepuls + TB_Repuls%E0_TB * Scell(NSC)%Na  ! [eV/atom]
         Scell(NSC)%nrg%E_rep = Erepuls ! [eV]
         
         ! van der Waals potential energy:
         Scell(NSC)%nrg%E_vdW = vdW_s(Scell(NSC)%TB_Waals, Scell, NSC, numpar)/real(Scell(NSC)%Na) ! [eV/atom], function below
         
         ! Coulomb potential energy:
         Scell(NSC)%nrg%E_coul = Coulomb_s(Scell(NSC)%TB_Coul, Scell, NSC, numpar)/real(Scell(NSC)%Na) ! [eV/atom], function below
         
         ! Exponential wall potential energy:
         Scell(NSC)%nrg%E_expwall = Exponential_wall_s(Scell(NSC)%TB_Expwall, Scell, NSC, numpar)/real(Scell(NSC)%Na) ! [eV/atom], function below
         !print*, 'get_pot_nrg 2'
      enddo
      
   endif DO_TB
end subroutine get_pot_nrg


function Exponential_wall_s(TB_Expwall, Scell, NSC, numpar) result(Pot)
   real(8) :: Pot	! Exponential wall energy [eV]
   class(TB_Exp_wall), allocatable, dimension(:,:), intent(in) :: TB_Expwall	! exponential wall
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8) a
   if (allocated(TB_Expwall)) then ! if we have Exponential wall potential defined
      select type(TB_Expwall)
      type is (TB_Exp_wall_simple)
         call get_Exp_wall_s(TB_Expwall, Scell, NSC, numpar, a)   ! module "Exponential_wall"
      end select
   else !For this material exponential wall class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   Pot = a ! [eV]
end function Exponential_wall_s



function Coulomb_s(TB_Coul, Scell, NSC, numpar)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_Coulomb), dimension(:,:), allocatable, intent(inout):: TB_Coul ! Coulomb parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8) :: Coulomb_s ! Coulomb energy [eV]
   real(8) a
   if (allocated(TB_Coul)) then ! if we have Coulomb potential defined
      select type(TB_Coul)
      type is (TB_Coulomb_cut)
         call get_Coulomb_s(TB_Coul, Scell, NSC, numpar, a) ! Coulomb energy, module "Coulomb"
      end select
   else !For this material Coulomb class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   Coulomb_s = a ! [eV]
end function Coulomb_s


function vdW_s(TB_Waals, Scell, NSC, numpar)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_vdW), dimension(:,:), allocatable, intent(inout):: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8) :: vdW_s ! van der Waals energy [eV]
   real(8) a
   if (allocated(TB_Waals)) then ! if we have vdW potential defined
      select type(TB_Waals)
      type is (TB_vdW_Girifalco)
         call get_vdW_s(TB_Waals, Scell, NSC, numpar, a)   ! van der Waals energy, module "Van_der_Waals"
      type is (TB_vdW_Dumitrica) ! UNFINISHED! DO NOT USE!
         call get_vdW_s_D(TB_Waals, Scell, NSC, numpar, a)   ! van der Waals energy, module "Van_der_Waals"
      end select
   else !For this material vdW class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   vdW_s = a ! [eV]
end function vdW_s



function vdW_interplane(TB_Waals, Scell, NSC, numpar, matter)
   real(8) :: vdW_interplane ! van der Waals energy [eV]
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_vdW), dimension(:,:), allocatable, intent(inout):: TB_Waals ! van der Waals parameters within TB
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8) a
   if (allocated(TB_Waals)) then ! if we have vdW potential defined
      select type(TB_Waals)
      type is (TB_vdW_Girifalco)
         call get_vdW_interlayer(TB_Waals, Scell, NSC, matter, numpar, a)   ! vdW interplane energy, module "Van_der_Waals"
      end select
   else !For this material vdW class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   vdW_interplane = a ! [eV]
end function vdW_interplane


FUNCTION Erep_s(TB_Repuls, Scell, NSC, numpar)   ! repulsive energy as a function of a distance
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !type(TB_Rep_Pettifor), dimension(:), intent(in) :: TB_Repuls   ! parameters of the repulsive part of TB-H
   class(TB_repulsive), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8) :: Erep_s
   !=====================================================
   real(8) a_r, a, b
   INTEGER(4) i1, j1, m, atom_2, NumTB
   real(8), DIMENSION(3) :: zb

   select type (TB_Repuls)
   type is (TB_Rep_Pettifor)
      call get_Erep_s(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Pettifor"
   type is (TB_Rep_Molteni)
      call get_Erep_s_M(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Molteni"
   type is (TB_Rep_Fu)
      call get_Erep_s_F(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Fu"
   type is (TB_Rep_NRL)
      call get_Erep_s_NRL(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_NRL"
   type is (TB_Rep_DFTB)
      call get_Erep_s_DFTB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_DFTB"
   type is (TB_Rep_3TB)
      call get_Erep_s_3TB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_3TB"
   type is (TB_Rep_BOP)
      call get_Erep_s_BOP(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_BOP"
   type is (TB_Rep_xTB)
!       call get_Erep_s_xTB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_xTB"
   end select
   Erep_s = a
END FUNCTION Erep_s



subroutine update_nrg_after_change(Scell, matter, numpar, time, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(in) :: time ! current timestep [fs]
   type(Error_handling), intent(inout) :: Err ! error save
   !=========================================
   integer NSC
   
    ! Update correspondingly the electron distribution for the new given total energy:
   call update_fe(Scell, matter, numpar, time, Err, do_E_tot=.true.) ! module "Electron_tools"
      
   do NSC = 1, size(Scell)
         ! Update correspondingly the electron distribution:
         !call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, Scell(NSC)%mu, Scell(NSC)%TeeV)
         !Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]
         !call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution
      if (ABS(Scell(NSC)%nrg%E_high_heating) > 0.0d0) then ! there was energy exchange
         ! Add energy to atoms recieved from elastic scattering of high-energy electrons:
         call Rescale_atomic_velocities(Scell(NSC)%nrg%E_high_heating, matter, Scell, NSC, Scell(NSC)%nrg) ! module "Atomic_tools"
         call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg) ! module "Atomic_tools"
         call get_new_energies(Scell, matter, numpar, time, Err) ! module "TB"
      endif
   enddo
end subroutine update_nrg_after_change



subroutine get_new_energies(Scell, matter, numpar, t, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(in) :: t ! current timestep [fs]
   type(Error_handling), intent(inout) :: Err ! error save
   !=========================================
   !print*, 'get_new_energies 0'

   ! Update electron energy according to updated distribution:
   call get_low_e_energy(Scell, matter) ! module "Electron_tools"
   !print*, 'get_new_energies 1'

   ! Get potential energy of the supercell (Repulsive part is defined inside):
   call get_pot_nrg(Scell, matter, numpar) ! see above
   !print*, 'get_new_energies 2'

   ! Get kinetic energy of atoms and of the supercell (for Parrinello-Rahman):
   call get_Ekin(Scell, matter) ! module "Atomic_tools"
   !print*, 'get_new_energies 3'

   ! And may be we want to update electron energy:
   call get_total_el_energy(Scell, matter, numpar, t, Err) ! module "Electron_tools"
   !print*, 'get_new_energies 4'
end subroutine get_new_energies



subroutine Electron_ion_coupling(t, matter, numpar, Scell, Err)
   real(8), intent(in) :: t ! [fs] current time
   type(solid), intent(inout) :: matter ! materil parameters
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err ! error save
   !=========================================
   real(8), dimension(:,:), allocatable :: Mij ! matrix element for electron-ion coupling
   real(8) :: dE_nonadiabat, E0
   integer NSC, ixm, iym, izm, ix, iy, iz, i, j

!    print*, 'Electron_ion_coupling 0'
   
   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
      SC:do NSC = 1, size(Scell) ! for all supercells
         dE_nonadiabat = 0.0d0
         if ((numpar%Nonadiabat) .AND. (t .GT. numpar%t_NA)) then ! electron-coupling included
            ! Ensure Fermi distribution:
            call update_fe(Scell, matter, numpar, t, Err) ! module "Electron_tools"
            
!             print*, 'Electron_ion_coupling 1'
            
            !--------------------------------
            ! For multiple k-points:
            ! Number of k-points to be sampled:
            ! Now it is only set for testing, to be changed later when the model works!!!
            ixm = 1
            iym = 1
            izm = 1

            do ix = 1, ixm
               do iy = 1, iym
                  do iz = 1, izm
           
                     ! Calculate nonadiabatic-coupling matrix element:
                     ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types
                        ! Different expressions for orthogonal and non-orthogonal bases:
                        select type(ARRAY)
                        type is (TB_H_Pettifor) ! TB parametrization according to Pettifor: orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                !call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij) ! module "Nonadiabatic"
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind) ! module "Nonadiabatic"
                                
!                                 do i = 1, size(Mij,1)
!                                    do j = 1, size(Mij,2)
!                                       i = 100
!                                       j = 60
!                                       print*, 'n', i, j, Mij(i,j)
!                                     enddo
!                                 enddo
!                                 call Electron_ion_coupling_Mij_OLD(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind) ! module "Nonadiabatic"
!                                 do i = 1, size(Mij,1)
!                                    do j = 1, size(Mij,2)
!                                       print*, 'o', i, j, Mij(i,j)
!                                     enddo
!                                 enddo
                                
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni: orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                !call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij) ! module "Nonadiabatic"
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind) ! module "Nonadiabatic"
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        type is (TB_H_Fu)  ! TB parametrization accroding to Fu: orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind) ! module "Nonadiabatic"
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        type is (TB_H_NRL)  ! TB parametrization accroding to NRL method: non-orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind, Scell(NSC)%Sij) ! module "Nonadiabatic"
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB: non-orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind, Sij=Scell(NSC)%Sij) ! module "Nonadiabatic"
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        type is (TB_H_3TB)  ! TB parametrization accroding to 3TB: non-orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind, Sij=Scell(NSC)%Sij) ! module "Nonadiabatic"
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        type is (TB_H_xTB)  ! TB parametrization accroding to xTB: non-orthogonal
                            if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
                                call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, numpar%NA_kind, Sij=Scell(NSC)%Sij) ! module "Nonadiabatic"
                            else ! non-gamma k-point, complex Hamiltonian needed:
                                call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, numpar%NA_kind) ! module "Nonadiabatic"
                            endif
                        end select
                     END ASSOCIATE
                     
                     
                     !sssssssssssssssssssssss
                     ! Project specific subroutine:
                     ! Test SAKUREI'S EXPRESSION:
!                     call Landau_vs_Sakurei_test(Mij, Scell(NSC)%Ha, Scell(NSC)%Ha0, Ha_non=Scell(NSC)%H_non, Ha_non0=Scell(NSC)%H_non0, wr=Scell(NSC)%Ei, wr0=Scell(NSC)%Ei0) ! module "Nonadiabatic"
                     !sssssssssssssssssssssss
                    
                     ! Calculate the electron-ion collision integral and energy exchange via it:
!                      print*, 'Electron_ion_coupling 4', numpar%DOS_splitting
!                      E0 = Scell(NSC)%nrg%El_low
!                      print*, 'before', E0
                     select case (numpar%DOS_splitting) ! include analysis of partial coupling (contribution of atomic shells) or not:
                     case (1)
                         call Electron_ion_collision_int(Scell(NSC), numpar, Scell(NSC)%nrg, Mij, Scell(NSC)%Ei, Scell(NSC)%Ei0, Scell(NSC)%fe, &
                            dE_nonadiabat, numpar%NA_kind, numpar%DOS_weights, Scell(NSC)%G_ei_partial) ! module "Nonadiabatic"
                     case default    ! No need to sort per orbitals
                        call Electron_ion_collision_int(Scell(NSC), numpar, Scell(NSC)%nrg, Mij, Scell(NSC)%Ei, Scell(NSC)%Ei0, Scell(NSC)%fe, &
                            dE_nonadiabat, numpar%NA_kind) ! module "Nonadiabatic"
                     endselect
!                      print*, 'Electron_ion_coupling 5'
                  enddo ! iz
               enddo ! iy
            enddo ! ix
            ! Averaging over multiple k-points

!             print*, 'middle', Scell(NSC)%nrg%El_low+dE_nonadiabat, E0-(Scell(NSC)%nrg%El_low+dE_nonadiabat)
            
            ! Redistribute this energy exchanged between electrons and atoms:
            ! New electron distribution:
            call update_fe(Scell, matter, numpar, t, Err) ! module "Electron_tools"

            ! New electron potential energy (Repulsive part is defined inside):
            call get_pot_nrg(Scell, matter, numpar) ! see below
            ! New atomic velosities:
            call Rescale_atomic_velocities(dE_nonadiabat, matter, Scell, NSC, Scell(NSC)%nrg) ! module "Atomic_tools"
            call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg) ! module "Atomic_tools"

            ! Update the last time-step data accordingly:
            call save_last_timestep(Scell) ! module "Atomic_tools"
            ! Calculate electron-ion coupling parameter:
            
!             print*, 'Electron_ion_coupling 6'
            call get_G_ei(Scell, NSC, numpar, dE_nonadiabat) ! module "Nonadiabatic"
            
!             print*, 'Electron_ion_coupling 7'
            
         endif
      enddo SC
   endif DO_TB
end subroutine Electron_ion_coupling



subroutine Construct_Aij(Ha, distre, Aij)
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha  ! diagonilized Hamiltonian/eigenvectors
   REAL(8), DIMENSION(:), INTENT(in) :: distre  ! electron distribution function
   REAL(8), DIMENSION(:,:), allocatable, INTENT(inout) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis
   integer i, k, N
   real(8), dimension(:,:), allocatable :: Atemp
   N = size(distre)
   if (.not. allocated(Aij)) allocate(Aij(N,N))
   allocate(Atemp(N,N))

   !$omp PARALLEL private(i)
   !$omp do !!reduction( + : Aij1)
   do i = 1, N
      Atemp(i,:) = distre(:)*Ha(i,:)
   enddo ! i
   !$omp end do
   !$omp end parallel

   call dgemm ('N','T', N, N, N, 1.0d0, Atemp, N, Ha, N, 0.0d0, Aij, N) ! mkl
   deallocate(Atemp)
   
   ! test:
!    do i = 1, N
!       do k = 1, N
!          write(*,'(a,i4,i4,f,f,f)') 'Aij:', i, k, Ha(k,i), Aij(i,k), Aij(k,i)
!       enddo
!    enddo
!    PAUSE 'Construct_Aij'
end subroutine Construct_Aij


subroutine Construct_Aij_x_En(Ha, distre, En, Aij)
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha 	! diagonilized Hamiltonian/eigenvectors
   REAL(8), DIMENSION(:), INTENT(in) :: distre 	! electron distribution function
   REAL(8), DIMENSION(:), INTENT(in) :: En	! eigenstates [eV]
   REAL(8), DIMENSION(:,:), allocatable, INTENT(inout) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis
   integer i, k, N
   real(8), dimension(:,:), allocatable :: Atemp
   N = size(distre)
   if (.not. allocated(Aij)) allocate(Aij(N,N))
   allocate(Atemp(N,N))

   !$omp PARALLEL private(i)
   !$omp do
   do i = 1, N
      Atemp(i,:) = distre(:)*En(:)*Ha(i,:)
   enddo ! i
   !$omp end do
   !$omp end parallel

   call dgemm ('N','T', N, N, N, 1.0d0, Atemp, N, Ha, N, 0.0d0, Aij, N) ! mkl
   deallocate(Atemp)
end subroutine Construct_Aij_x_En


subroutine Construct_Aij_old(Ha, distre, Aij)
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha  ! diagonilized Hamiltonian/eigenvectors
   REAL(8), DIMENSION(:), INTENT(in) :: distre  ! electron distribution function
   REAL(8), DIMENSION(:,:), allocatable, INTENT(inout) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis 
   integer i, k, nat4, j
   real(8) Vect1(size(Ha,2)), Vect2(size(Ha,2))
   real(8), dimension(:,:), allocatable :: Aij1
   nat4 = size(Ha,2)
   if (.not. allocated(Aij)) allocate(Aij(nat4,nat4))
   Aij = 0.0d0
   if (minval(distre) .LE. 1.0d-16) then
      j = 1
      do while (distre(j) .GT. 1.0d-16)
         j = j + 1
      enddo
   else
      j = nat4
   endif
   !$omp PARALLEL private(i,k,Vect1)
   !$omp do !!reduction( + : Aij1)
   do i = 1,nat4
      Vect1(:) = distre(:)*Ha(i,:)
      do k = i,nat4
         !Aij(i,k) = SUM(distre(2,:)*Ha(i,:)*Ha(k,:), MASK = distre(2,:) .GT. 0)
         
!          Aij(i,k) = DOT_PRODUCT (Vect1(1:j),Ha(k,1:j))
         Aij(i,k) = SUM(Vect1(1:j)*Ha(k,1:j))
         if (k .NE. i) then
            Aij(k,i) = Aij(i,k) ! it's simmetric
         endif
      enddo ! k

   enddo ! i
   !$omp end do
   !$omp end parallel
end subroutine Construct_Aij_old


!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Analysis subroutines:

subroutine Get_configurational_temperature(Scell, numpar, Tconf)
   type(Super_cell), dimension(:), intent(in) :: Scell	! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: Tconf	! [K] configurational temperature
   !--------------------------------------------
   real(8), dimension(:,:), allocatable :: F, dF	! forces and derivatives
   integer :: Nat
   Nat = size(Scell(1)%MDAtoms)	! number of atoms
   allocate(F(3,Nat))
   allocate(dF(3,Nat))
   F = 0.0d0
   dF = 0.0d0
   
   ! Forces and their derivatives:
   call get_derivatives_and_forces_r(Scell, numpar, F, dF)	! see below

   ! Configurational temperature:
   Tconf = SUM( (F(1,:)*F(1,:) + F(2,:)*F(2,:) + F(3,:)*F(3,:)) ) / SUM( (dF(1,:)+dF(2,:)+dF(3,:)) )	! [eV]
   Tconf = Tconf*g_kb	! [eV] -> [K]
   
!     print*, '===================='
!     print*, Tconf, SUM( (F(1,:)*F(1,:) + F(2,:)*F(2,:) + F(3,:)*F(3,:))), SUM((dF(1,:)+dF(2,:)+dF(3,:)) )
!     print*, '===================='
!     pause 'Tconfig'
   
   ! Clean up:
   deallocate(F, dF)
end subroutine Get_configurational_temperature


subroutine get_derivatives_and_forces_r(Scell, numpar, F, dF)
   type(Super_cell), dimension(:), intent(in) :: Scell	! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), dimension(:,:), allocatable, intent(inout):: F, dF	! forces and derivatives [eV/A], [eV/A^2]
   !-----------------------------------------------------------
   real(8), dimension(:,:), allocatable :: Frep, Fatr, dFrep, dFatr	! forces and derivatives [eV/A], [eV/A^2]
   real(8), dimension(:,:,:), allocatable :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), allocatable :: M_dVs ! matrix of functions dVs
   real(8), dimension(:,:,:), allocatable :: M_d2Vs ! matrix of functions d2Vs
   real(8), dimension(:,:,:), allocatable :: M_cos	! matrix of directional cosines
   integer :: i, N
   
!    open(UNIT = 7774, FILE = 'OUTPUT_Forces_for_Tconfig.dat') !<-
   
   N = size(Scell(1)%MDAtoms)	! number of atoms
   if (.not.allocated(F)) allocate(F(3,N))
   if (.not.allocated(dF)) allocate(dF(3,N))
   
   call Construct_M_cos(Scell(1), M_cos)	! see below
   
   ! Electronic TB Hamiltonian part:
   ASSOCIATE (ARRAY => Scell(1)%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types      
      ! Construct forces and their derivatives:
      select type(ARRAY)
      type is (TB_H_Pettifor) ! TB parametrization according to Pettifor
      
         ! Construct array of functions Vs and dVs for all pairs of atoms to use for forces:
         call Construct_M_Vs(Scell, 1, ARRAY, M_Vs, M_dVs, M_d2Vs) ! module "TB_Pettifor"

         ! Get attractive forces for atoms from the derivatives of the Hamiltonian:
         call dHij_r(ARRAY, Scell(1)%MDatoms, Scell, numpar, M_Vs, M_dVs, M_d2Vs, M_cos, Fatr, dFatr) ! module "TB_Pettifor"
      type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni
         print*, 'Configurational temperature calculations are not implemented for Molteni: Attractive'
      type is (TB_H_Fu)  ! TB parametrization accroding to Fu
         print*, 'Configurational temperature calculations are not implemented for Fu: Attractive'
      type is (TB_H_NRL)  ! TB parametrization accroding to NRL
         print*, 'Configurational temperature calculations are not implemented for NRL: Attractive'
      type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB
         print*, 'Configurational temperature calculations are not implemented for DFTB: Attractive'
      type is (TB_H_3TB)  ! TB parametrization accroding to 3TB
         print*, 'Configurational temperature calculations are not implemented for 3TB: Attractive'
      type is (TB_H_xTB)  ! TB parametrization accroding to xTB
         print*, 'Configurational temperature calculations are not implemented for xTB: Attractive'
      end select
   END ASSOCIATE

!     ! Test of calculations:
!     do i = 1, N
!        write(*,'(a)') 'Attractive'
!        write(*,'(a,i3,f,f,f)') 'X', i, Scell(1)%MDAtoms(i)%forces%att(1)/Scell(1)%supce(1,1), Fatr(1,i), dFatr(1,i)
!        write(*,'(a,i3,f,f,f)') 'Y', i, Scell(1)%MDAtoms(i)%forces%att(2)/Scell(1)%supce(2,2), Fatr(2,i), dFatr(2,i)
!        write(*,'(a,i3,f,f,f)') 'Z', i, Scell(1)%MDAtoms(i)%forces%att(3)/Scell(1)%supce(3,3), Fatr(3,i), dFatr(3,i)
!     enddo

   ! Repulsive TB Hamiltonian part:
   ASSOCIATE (ARRAY2 => Scell(1)%TB_Repuls(:,:))
      select type(ARRAY2)
      type is (TB_Rep_Pettifor) ! TB parametrization according to Pettifor
         ! Second derivatives of the repulsive energy by r:
         call dE2rep_dr2(ARRAY2, Scell(1)%MDAtoms, Scell, numpar, Frep, dFrep) ! module "TB_Pettifor"
      type is (TB_Rep_Molteni)  ! TB parametrization accroding to Molteni
         print*, 'Configurational temperature calculations are not implemented for Molteni: Repulsive'
      type is (TB_Rep_Fu)  ! TB parametrization accroding to Fu
         print*, 'Configurational temperature calculations are not implemented for Fu: Repulsive'
      type is (TB_Rep_NRL)  ! TB parametrization accroding to NRL
         print*, 'Configurational temperature calculations are not implemented for NRL: Repulsive'
      type is (TB_Rep_DFTB)  ! TB parametrization accroding to DFTB
         print*, 'Configurational temperature calculations are not implemented for DFTB: Repulsive'
      type is (TB_Rep_3TB)  ! TB parametrization accroding to 3TB
         print*, 'Configurational temperature calculations are not implemented for 3TB: Repulsive'
      type is (TB_Rep_BOP)  ! TB parametrization accroding to BOP
         print*, 'Configurational temperature calculations are not implemented for BOP: Repulsive'
      type is (TB_Rep_xTB)  ! TB parametrization accroding to xTB
         print*, 'Configurational temperature calculations are not implemented for xTB: Repulsive'
      end select
   END ASSOCIATE !    
   
   ! Test of calculations:
!      do i = 1, N
!         write(*,'(a)') 'Repulsive'
!         write(*,'(a,i3,f,f,f)') 'X', i, Scell(1)%MDAtoms(i)%forces%rep(1)/Scell(1)%supce(1,1), Frep(1,i), dFrep(1,i)
!         write(*,'(a,i3,f,f,f)') 'Y', i, Scell(1)%MDAtoms(i)%forces%rep(2)/Scell(1)%supce(2,2), Frep(2,i), dFrep(2,i)
!         write(*,'(a,i3,f,f,f)') 'Z', i, Scell(1)%MDAtoms(i)%forces%rep(3)/Scell(1)%supce(3,3), Frep(3,i), dFrep(3,i)
!      enddo

!      do i = 1, N
!         write(7774,'(i3,es,es,es,es,es,es,es,es,es,es,es,es)') i, Frep(:,i), Fatr(:,i),  dFrep(:,i), dFatr(:,i)
!      enddo
!      pause 'WRINTING OUT FORCES AND DERIVATIVES'
!    close(7774)
   
   ! Combine attractive and repulsive parts of forces and derivatives:
   F = Frep + Fatr	! forces [eV/A]
   dF = dFrep + dFatr	! derivatives of forces [eV/A^2]
   
   ! Clean up:
   if (allocated(Frep)) deallocate(Frep)
   if (allocated(Fatr)) deallocate(Fatr)
   if (allocated(dFrep)) deallocate(dFrep)
   if (allocated(dFatr)) deallocate(dFatr)
   if (allocated(M_Vs)) deallocate(M_Vs)
   if (allocated(M_dVs)) deallocate(M_dVs)
   if (allocated(M_d2Vs)) deallocate(M_d2Vs)
   if (allocated(M_cos)) deallocate(M_cos)
end subroutine get_derivatives_and_forces_r




subroutine Construct_M_cos(Scell,  M_cos)
   type(Super_cell), intent(in), target :: Scell	! supercell with all the atoms as one object
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_cos	! matrix of directional cosines
   !---------------------------
   real(8), pointer :: r, x, y, z
   integer, pointer :: nat,  m, j
   integer i, atom_2
   nat => Scell%Na ! number of atoms

   if (.not.allocated(M_cos)) allocate(M_cos(nat,nat,3))	! for 3 directions (x, y, z)  and pairs of atoms
   M_cos = 0.0d0
   
   !$omp PARALLEL DO private(i, m, atom_2, j, x, y, z, r)
   do i = 1,nat	! all atoms
      m => Scell%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            r => Scell%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
            x => Scell%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
            y => Scell%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
            z => Scell%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
            
            ! Directional cosines:
            M_cos(i,j,1) = x/r
            M_cos(i,j,2) = y/r
            M_cos(i,j,3) = z/r
         endif !  (j .GT. 0) 
      enddo ! atom_2 = 1,m 
   enddo ! do i = 1,nat
   !$omp END PARALLEL DO 
   nullify(nat, m, j, r, x, y, z)
end subroutine Construct_M_cos





subroutine Get_pressure(Scell, numpar, matter, P, stress_tensor_OUT)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter ! material parameters
   real(8), intent(out) :: P	! Pressure [Pa]
   real(8), dimension(3,3), intent(out), optional :: stress_tensor_OUT	! Stress tensor [Pa]
   !-------------------------------------
   real(8), dimension(3,3) :: forces_saved, forces_saved0, sigma_tensor, sigma_inversed, stress_tensor
   integer :: NSC
   logical :: P_const_save
    real(8), dimension(:,:,:), allocatable :: M_Vij	! matrix of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dVij	! matrix of derivatives of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_SVij	! matrix of Overlap functions for Overlap matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dSVij	! matrix of derivatives of Overlap functions for Overlap matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_lmn ! matrix of directional cosines l, m, n
   real(8), dimension(:,:,:), allocatable :: M_x1  ! matrix of x1 elements, used for forces calculations
   real(8), dimension(:,:,:), allocatable :: M_xrr ! matrix of xrr elements, used for forces calculations
   real(8), dimension(:,:), allocatable :: M_Aij_x_Ei ! matrix of multipliers for non-orthogonal forces calculation
   real(8), dimension(:,:,:), allocatable :: M_Lag_exp   ! matrix of Laguerres for 3-body radial funcs (for 3TB)
   real(8), dimension(:,:,:), allocatable :: M_d_Lag_exp   ! matrix of derivatives of Laguerres for 3-body radial funcs
   real(8), dimension(:,:,:), allocatable :: Mjs      ! matrix of K-S part of overlaps with s-orb. (for 3TB)
   real(8), dimension(:,:,:), allocatable :: M_drij_dsk  ! matrix of derivatives of distances between atoms i and j
   real(8), dimension(:,:,:), allocatable :: M_dlmn   ! matrix of derivatives of directional cosines

   
   ! so far we only have one supercell:
   NSC = 1
   ! Save the values that the forces had before calculations of pressure:
   forces_saved = Scell(NSC)%SCforce%total
   forces_saved0 = Scell(NSC)%SCforce%total0
   P_const_save = numpar%p_const ! P=const vs V=const as defined by the user
   ! Temporarely set it P=const so that the subroutines would calculate the supercell forces:
   numpar%p_const = .true.	! P=const
   
   ! Electronic TB Hamiltonian part:
   ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types
      ! Construct forces acting on the supercell:
      select type(ARRAY)
      type is (TB_H_Pettifor) ! TB parametrization according to Pettifor
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Attract_TB_Forces_Press(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij)	! module "TB_Pettifor"
      type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Attract_TB_Forces_Press_M(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij)	! module "TB_Molteni"
      type is (TB_H_Fu) ! TB parametrization according to Fu
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Attract_TB_Forces_Press_F(ARRAY, Scell(NSC)%MDatoms, Scell, NSC, numpar, Scell(NSC)%Aij)	! module "TB_Fu"
      type is (TB_H_NRL) ! TB parametrization according to NRL
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn) ! see below 
         call Construct_Vij_NRL(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_NRL"
         call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
         call Attract_TB_Forces_Press_NRL(ARRAY, Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_NRL"
      type is (TB_H_DFTB) ! TB parametrization according to DFTB
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn) ! see below 
         call Construct_Vij_DFTB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_DFTB"
         call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
         call Attract_TB_Forces_Press_DFTB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_DFTB"
      type is (TB_H_3TB) ! TB parametrization according to 3TB
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn) ! see below
         call get_Mjs_factors(numpar%N_basis_size, Scell(NSC), M_lmn, Mjs, M_drij_dsk, M_dlmn)   ! module "TB_3TB"
         call Construct_Vij_3TB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_Lag_exp, M_d_Lag_exp)	! module "TB_3TB"
         call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
         call Attract_TB_Forces_Press_3TB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_3TB"
      type is (TB_H_xTB) ! TB parametrization according to xTB
         ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
         call Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn) ! see below
!          call Construct_Vij_xTB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_xTB"
         call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
!          call Attract_TB_Forces_Press_xTB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_xTB"
      end select
   END ASSOCIATE
   
   ! Repulsive TB Hamiltonian part:
   ASSOCIATE (ARRAY2 => Scell(NSC)%TB_Repuls(:,:))
      select type(ARRAY2)
      type is (TB_Rep_Pettifor) ! TB parametrization according to Pettifor
         ! Get repulsive forces acting on the supercell:
         call dErdr_Pressure_s(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_Pettifor"
      type is (TB_Rep_Molteni)  ! TB parametrization accroding to Molteni
         ! Get repulsive forces acting on the supercell:
         call dErdr_Pressure_s_M(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_Molteni"
      type is (TB_Rep_Fu)  ! TB parametrization accroding to Fu
         ! Get repulsive forces acting on the supercell:
         call dErdr_Pressure_s_F(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_Fu"
      type is (TB_Rep_NRL) ! TB parametrization according to Mehl (NRL)
         call dErdr_Pressure_s_NRL(ARRAY2, Scell(NSC)%MDatoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_NRL"
      type is (TB_Rep_DFTB) ! TB parametrization according to DFTB
         call dErdr_Pressure_s_DFTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_DFTB"
      type is (TB_Rep_3TB) ! TB parametrization according to 3TB
         call dErdr_Pressure_s_3TB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_3TB"
      type is (TB_Rep_BOP) ! TB parametrization according to BOP
!          call dErdr_Pressure_s_DFTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_BOP"
      type is (TB_Rep_xTB) ! TB parametrization according to xTB
!          call dErdr_Pressure_s_xTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_xTB"
      end select
   END ASSOCIATE
   
   ! Get forces for the supercell:
   call Potential_super_cell_forces(numpar, Scell, NSC, matter)  ! module "Atomic_tools"

   ! Adding kinetic part and pressure to super-cell forces:
   call super_cell_forces(numpar, Scell, NSC, matter, Scell(NSC)%SCforce, sigma_tensor) ! module "Atomic_tools"
   
   ! Invert sigma tensor:
   call  Invers_3x3(sigma_tensor, sigma_inversed, 'Get_pressure')	! module "Algebra_tools"
   
   ! Calculate the stress tensor (factor 1.040 is to convert from [kg/A/fs^2] to [kg/m/s^2]):
   stress_tensor(:,:) = Scell(NSC)%SCforce%total(:,:) * matter%W_PR * sigma_inversed(:,:) * 1.0d40	! [kg/m/s^2]
   
   ! OUTPUT, pressure and stress tensor:
   P = 1.0d0/3.0d0 * ( stress_tensor(1,1) + stress_tensor(2,2) + stress_tensor(3,3) )
   if (present(stress_tensor_OUT)) stress_tensor_OUT = stress_tensor

   ! Restore the forces to the values they had before calculations of pressure:
   Scell(NSC)%SCforce%total = forces_saved
   Scell(NSC)%SCforce%total0 = forces_saved0
   numpar%p_const = P_const_save
end subroutine Get_pressure



!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Project-specific analysis of C60:

subroutine C60_vdW_vs_Coulomb(Scell, numpar, matter, layers)
   type(Super_cell), dimension(:), intent(in) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(solid), intent(in) :: matter	! material parameters
   integer, intent(in), optional :: layers ! number of C60 layers in the crystal 
   type(Super_cell), dimension(:), allocatable :: TEMP_Scell ! temporary supercell with all the atoms as one object
   integer :: i, j, Nat
   if (present(layers)) then
      Nat = size(Scell(1)%MDatoms) ! number of atoms in one C60
      Nat = Nat*layers ! number of atoms in all C60s
      allocate(TEMP_Scell(1)) ! only one cell always
      TEMP_Scell(1)%Na = Nat  ! number of atoms
      allocate(TEMP_Scell(1)%MDatoms(Nat)) ! all atoms are storred in here
      forall (i=1:Nat, j=1:3) TEMP_Scell(1)%MDatoms(i)%R(j) = 0.0d0   ! just to start

      ! Construct C60 crystall:
      call C60_set_coordinates(Scell, TEMP_Scell(1)%MDatoms, matter, layers) ! set coordinates for all atoms, see below
      call C60_super_cell(TEMP_Scell(1)) ! set supercell size and relative coords, see below
      
      ! Transfere parameters of potentials to this temporary super-cell:
      call Transfere_Temp_conditions(Scell, TEMP_Scell(1), Nat) ! see below
      
      ! Find the threshold charge for which Coulomb potential overcomes van der Waals for the uppermost C60 cage:
      call Coulomb_beats_vdW(TEMP_Scell, numpar, size(Scell(1)%MDatoms), Nat) ! see below
   endif
   pause 'C60_vdW_vs_Coulomb'
end subroutine C60_vdW_vs_Coulomb



subroutine Coulomb_beats_vdW(TEMP_Scell, numpar, Nat_C60, Nat) 
   type(Super_cell), dimension(:), intent(inout) :: TEMP_Scell ! temporary supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   integer, intent(in), optional :: Nat_C60, Nat
   integer i, j, coun
   real(8) :: Charge ! unballanced charge
   real(8) :: Sum_energy, vdW_nrg, Coulomb_nrg
   real(8) a, b ! parameters for bisection method
   a = 0.0d0
   b = 1.0d0
   
   ! van der Waals energy:
   vdW_nrg = vdW_s(TEMP_Scell(1)%TB_Waals, TEMP_Scell, 1, numpar) ! vdW energy
   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   coun = 0
   do while (abs(a-b)/abs(b) >= 1.0d-4)
      coun = coun + 1
      Charge = (a+b)/2.0d0
      TEMP_Scell(1)%Q = Charge !set unballanced charge per atom
      ! Coulomb potential part for modelling Coulomb explosion of a finite system:
      Coulomb_nrg = Coulomb_s(TEMP_Scell(1)%TB_Coul, TEMP_Scell, 1, numpar) ! get Coulomb energy

      Sum_energy = vdW_nrg + Coulomb_nrg ! total energy vdW + Coulomb
      write(*,'(a,f,f,f)') 'Total energy between C60:', a,b, Sum_energy
      if (Sum_energy > 0.0d0) then
         b = Charge
      else
         a = Charge
      endif
      if (coun > 1000) then
         print*, 'Too many iterations in the bisection loop:', coun
         exit
      endif
   enddo
end subroutine Coulomb_beats_vdW



subroutine Transfere_Temp_conditions(Scell, TEMP_Scell, Nat)
   type(Super_cell), dimension(:), intent(in) :: Scell ! supercell with all the atoms as one object
   type(Super_cell), intent(inout) :: TEMP_Scell ! temporary supercell with all the atoms as one object
   integer, intent(in) :: Nat ! number of atoms inside of the temporary cell
   integer i, j

   ! Set vdW parameters:
   ASSOCIATE (ARRAY => Scell(1)%TB_Waals(:,:))
      select type (ARRAY)
         type is (TB_vdW_Girifalco)
         allocate(TB_vdW_Girifalco::TEMP_Scell%TB_Waals(1,1)) ! make it for Girifalco parametrization
         ASSOCIATE (TEMP_ARRAY => TEMP_Scell%TB_Waals(:,:))
         select type (TEMP_ARRAY)
            type is (TB_vdW_Girifalco)
            TEMP_ARRAY = ARRAY ! copy from the main super-cell
          end select
         END ASSOCIATE
      end select
   END ASSOCIATE
   
   ! Set Coulomb (with soft cut-off) parameters:
   ASSOCIATE (ARRAY2 => Scell(1)%TB_Coul(:,:))
      select type (ARRAY2)
         type is (TB_Coulomb_cut)
         allocate(TB_Coulomb_cut::TEMP_Scell%TB_Coul(1,1)) ! make it for Coulomb parametrization
         ASSOCIATE (TEMP_ARRAY => TEMP_Scell%TB_Coul(:,:))
         select type (TEMP_ARRAY)
            type is (TB_Coulomb_cut)
            TEMP_ARRAY = ARRAY2 ! copy from the main super-cell
            ! Replace the Coulomb cut-off by the proper value:
            TEMP_ARRAY(1,1)%dm = maxval(TEMP_Scell%MDatoms(:)%R(3)) - minval(TEMP_Scell%MDatoms(:)%R(3)) ! allow all atoms to interact with each other
         end select
         END ASSOCIATE
      end select
   END ASSOCIATE
   
   
   ! TESTING:
   ASSOCIATE (ARRAY => Scell(1)%TB_Waals(1,1))
      select type (ARRAY)
         type is (TB_vdW_Girifalco)
         ASSOCIATE (TEMP_ARRAY => TEMP_Scell%TB_Waals(1,1))
         select type (TEMP_ARRAY)
            type is (TB_vdW_Girifalco)
            ASSOCIATE (ARRAY2 => Scell(1)%TB_Coul(1,1))
            select type (ARRAY2)
               type is (TB_Coulomb_cut)
               ASSOCIATE (TEMP_ARRAY2 => TEMP_Scell%TB_Coul(1,1))
               select type (TEMP_ARRAY2)
               type is (TB_Coulomb_cut)
                  call print_param_temp(ARRAY,TEMP_ARRAY,ARRAY2,TEMP_ARRAY2)
               end select
               END ASSOCIATE
            end select
            END ASSOCIATE   
          end select
         END ASSOCIATE
      end select
   END ASSOCIATE
end subroutine Transfere_Temp_conditions


subroutine print_param_temp(ScellTB_Waals, TB_Waals, ScellTB_Coul, TB_Coul)
   type(TB_vdW_Girifalco), intent(in) :: ScellTB_Waals, TB_Waals
   type(TB_Coulomb_cut), intent(in) :: ScellTB_Coul, TB_Coul
   print*, 'REAL vdW:', ScellTB_Waals
   print*, 'TEMP vdW:', TB_Waals
   print*, 'REAL Coulomb:', ScellTB_Coul
   print*, 'TEMP Coulomb:', TB_Coul
end subroutine print_param_temp


! This subroutine prepares coordinates for array of C60 molecules along -direction,
! separated by the equilibrium van der Waals dstance.
! This subroutine works only for the input material: C60 (only ONE C60!)
subroutine C60_set_coordinates(Scell, Atoms, matter, layers) ! there is specific arrengement for all atoms (in this project)
   type(Super_cell), dimension(:), intent(in) :: Scell ! supercell with all the atoms as one object
   type(Atom), dimension(:), intent(inout) :: Atoms ! Atoms in all the C60 molecules to be analyzed
   type(solid), intent(in) :: matter	! material parameters
   integer, intent(in) :: layers ! number of C60 layers in the crystal 
   real(8) :: Z_C60, X_min, Y_min, Z_min, Full_size
   integer :: i, j, k, Nat, NatC60
   NatC60 = size(Scell(1)%MDatoms) ! number of atoms in one C60
   Nat = size(Atoms) ! number of atoms in one C60
   Z_C60 = maxval(Scell(1)%MDatoms(:)%R(3)) - minval(Scell(1)%MDatoms(:)%R(3)) ! size of C60 in Z direction
   ! Starting points of the C60 molecule:
   X_min = minval(Scell(1)%MDatoms(:)%R(1))
   Y_min = minval(Scell(1)%MDatoms(:)%R(2))
   Z_min = minval(Scell(1)%MDatoms(:)%R(3))
   k = 0
   do i = 1, layers ! for all C60 molecules
      do j = 1, NatC60
         k = k + 1 ! curent atom in all C60s
         Atoms(k)%R(1:2) = Scell(1)%MDatoms(j)%R(1:2) ! absolute atomic coordinates X and Y
         Atoms(k)%R(3) = Scell(1)%MDatoms(j)%R(3) + dble(i-1)*(Z_C60 + 3.35d0)  ! Shift atoms in Z direction by the distance of equilibrium vdW potential
         Atoms(k)%KOA = Scell(1)%MDatoms(j)%KOA ! kind of atom
      enddo ! j
   enddo ! j
   Full_size = maxval(Atoms(:)%R(3)) - minval(Atoms(:)%R(3)) ! size of all C60s in Z direction
!    write(*,'(a)') 'Atomic coordinates:'
   do i=1,Nat ! for all atoms of C60s, shift coordinates to form crystall:
      Atoms(i)%R(1) = Atoms(i)%R(1) - X_min + 3.35d0/2.0d0
      Atoms(i)%R(2) = Atoms(i)%R(2) - Y_min + 3.35d0/2.0d0
      Atoms(i)%R(3) = Atoms(i)%R(3) - Z_min + Full_size
!       write(*,'(i4,f,f,f)') Atoms(i)%KOA, Atoms(i)%R
   end do

end subroutine C60_set_coordinates


subroutine C60_super_cell(Scell) ! set supercell size and relative coords
   type(Super_cell), intent(inout) :: Scell ! temporary supercell with all the atoms as one object
   integer :: i,Nat
   real(8) X_C60, Y_C60, Z_C60
   Nat = size (Scell%MDatoms(:)) ! total number of atoms
   X_C60 = maxval(Scell%MDatoms(:)%R(1)) - minval(Scell%MDatoms(:)%R(1)) ! size of C60 in X direction
   Y_C60 = maxval(Scell%MDatoms(:)%R(2)) - minval(Scell%MDatoms(:)%R(2)) ! size of C60 in Y direction
   Z_C60 = maxval(Scell%MDatoms(:)%R(3)) - minval(Scell%MDatoms(:)%R(3)) ! size of C60 in Z direction
   Scell%supce(:,:) = 0.0d0
   Scell%supce(1,1) = X_C60 + 3.35d0 ! image cells are separated by equilibrium vdW distance
   Scell%supce(2,2) = Y_C60 + 3.35d0 ! image cells are separated by equilibrium vdW distance
   Scell%supce(3,3) = Z_C60*3.0d0    ! free boundaries along Z-direction
!    print*, 'Size of C60:'
!    print*, X_C60, Y_C60, Z_C60
   write(*,'(a)') 'Super cell vectors coordinates:'
   write(*,'(f,f,f)') Scell%supce(1,:)
   write(*,'(f,f,f)') Scell%supce(2,:)
   write(*,'(f,f,f)') Scell%supce(3,:)
   write(*,'(a)') 'Relative atomic coordinates:'
   do i=1,Nat ! for all atoms of C60s, shift coordinates to form crystall:
      Scell%MDatoms(i)%S(1) = Scell%MDatoms(i)%R(1)/Scell%supce(1,1)
      Scell%MDatoms(i)%S(2) = Scell%MDatoms(i)%R(2)/Scell%supce(2,2)
      Scell%MDatoms(i)%S(3) = Scell%MDatoms(i)%R(3)/Scell%supce(3,3)
      write(*,'(i4,f,f,f)') Scell%MDatoms(i)%KOA, Scell%MDatoms(i)%S
      Scell%MDatoms(i)%forces%rep(:) = 0.0d0 ! we start from zero and then calculate the forces
   end do
end subroutine C60_super_cell


END MODULE TB
