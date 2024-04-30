! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT-3
! available at: https://doi.org/10.48550/arXiv.2307.03953
! or at: https://github.com/N-Medvedev/XTANT-3
!
! Developed by Nikita Medvedev
!
! XTANT-3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Although we endeavour to ensure that the code XTANT-3 and results delivered are correct,
! no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions.
! We shall not be liable for any damage arising from the use of this code or its parts
! or any results produced with it, or from any action or decision taken
! as a result of using this code or any related material.
!
! This code is distributed as is for non-commercial peaceful purposes only,
! such as research and education. The code, its parts, its results or any related material
! should never be used for military-related and other than peaceful purposes.
!
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with general aspects of TB

MODULE TB
use Universal_constants
use Objects
use Algebra_tools, only : Det_3x3, get_eigenvalues_from_eigenvectors, fit_parabola_to_3points, sym_diagonalize, Reciproc, &
                        mkl_matrix_mult, Invers_3x3, mkl_matrix_mult_c8, c8_diagonalize
use Little_subroutines, only : number_of_types_of_orbitals, count_3d, deallocate_array
use Atomic_tools, only : get_near_neighbours, total_forces, Potential_super_cell_forces, super_cell_forces, &
                        Convert_reciproc_rel_to_abs, Rescale_atomic_velocities, get_kinetic_energy_abs, &
                        get_Ekin, save_last_timestep, Potential_super_cell_forces, &
                        make_time_step_atoms, make_time_step_supercell, make_time_step_atoms_Y4, make_time_step_supercell_Y4, &
                        make_time_step_atoms_M, distance_to_given_cell, shortest_distance
use Electron_tools, only : set_initial_fe, update_fe, get_new_global_energy, find_band_gap, get_DOS_sort, &
                     get_electronic_heat_capacity, electronic_entropy, Diff_Fermi_E, get_low_e_energy, get_total_el_energy, &
                     get_orbital_resolved_data
use Nonadiabatic, only : Electron_ion_coupling_Mij, Electron_ion_coupling_Mij_complex, Electron_ion_collision_int, get_G_ei
use TB_Fu, only : dHij_s_F, Attract_TB_Forces_Press_F, dErdr_s_F, dErdr_Pressure_s_F, construct_TB_H_Fu, &
                     Complex_Hamil_tot_F, Attract_TB_Forces_Press_F, get_Erep_s_F, dErdr_Pressure_s_F
use TB_Pettifor, only : dHij_s, Attract_TB_Forces_Press, dErdr_s, dErdr_Pressure_s, construct_TB_H_Pettifor, &
                     Complex_Hamil_tot, Attract_TB_Forces_Press, get_Erep_s, dErdr_Pressure_s, Construct_M_Vs, &
                     dHij_r, dE2rep_dr2
use TB_Molteni, only : dHij_s_M, Attract_TB_Forces_Press_M, dErdr_s_M, dErdr_Pressure_s_M, construct_TB_H_Molteni, &
                     Complex_Hamil_tot_Molteni, Attract_TB_Forces_Press_M, get_Erep_s_M, dErdr_Pressure_s_M
use TB_NRL, only : get_dHij_drij_NRL, dErdr_s_NRL, Construct_Vij_NRL, construct_TB_H_NRL, Complex_Hamil_NRL, &
                     get_Erep_s_NRL, dErdr_Pressure_s_NRL, Loewdin_Orthogonalization_c, m_two_third, &
                     Attract_TB_Forces_Press_NRL, Loewdin_Orthogonalization_c8
use TB_DFTB, only : Construct_Vij_DFTB, construct_TB_H_DFTB, get_Erep_s_DFTB, get_dHij_drij_DFTB, &
                     Attract_TB_Forces_Press_DFTB, dErdr_s_DFTB, dErdr_Pressure_s_DFTB, Complex_Hamil_DFTB, &
                     identify_DFTB_orbitals_per_atom, dErdr_s_DFTB_no, get_Erep_s_DFTB_no, dErdr_Pressure_s_DFTB_no
use TB_3TB, only : get_Erep_s_3TB, dErdr_s_3TB, dErdr_Pressure_s_3TB, Attract_TB_Forces_Press_3TB, &
                     Construct_Vij_3TB, construct_TB_H_3TB, get_Mjs_factors, get_dHij_drij_3TB
use TB_BOP, only : Construct_Vij_BOP, construct_TB_H_BOP, get_Erep_s_BOP
use TB_xTB, only : Construct_Vij_xTB, construct_TB_H_xTB, get_Erep_s_xTB, identify_xTB_orbitals_per_atom
use Van_der_Waals, only : Construct_B, get_vdW_s, get_vdW_s_D, get_vdW_interlayer, d_vdW_forces
use Coulomb, only: m_k, m_sqrtPi, Coulomb_Wolf_pot, get_Coulomb_Wolf_s, cut_off_distance, Construct_B_C, get_Coulomb_s, &
                     Coulomb_Wolf_self_term, d_Coulomb_Wolf_pot, d_Coulomb_forces
use Exponential_wall, only : get_Exp_wall_s, d_Exp_wall_pot_s, d_Exp_wall_Pressure_s, &
                     get_short_range_rep_s, d_Short_range_pot_s, d_Short_range_Pressure_s, d_Exponential_wall_forces

#ifdef OMP_inside
   USE OMP_LIB, only : OMP_GET_THREAD_NUM
#endif

implicit none
PRIVATE

public :: get_new_energies, get_DOS, Get_pressure, get_electronic_thermal_parameters, &
         vdW_interplane, Electron_ion_coupling, update_nrg_after_change, get_DOS_masks, k_point_choice, &
         construct_complex_Hamiltonian, get_Hamilonian_and_E, MD_step, get_Mullikens_all, get_coupling_matrix_elements, &
         Get_configurational_temperature

 contains

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!Common TB tools:


subroutine MD_step(Scell, matter, numpar, time, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   real(8), intent(in) :: time ! [fs] timestep
   !------------------------------------------
   !real(8), dimension(size(Scell(1)%MDAtoms),3) :: V0
   integer :: i, Nat


   ! Save the velosities:
   Nat = size(Scell(1)%MDAtoms)
   do i = 1, Nat
      !V0(i,:) = Scell(1)%MDAtoms(i)%V(:)     ! before update
      Scell(1)%MDatoms(i)%accel0(:) = Scell(1)%MDatoms(i)%accel(:)   ! save from last step
      Scell(1)%MDatoms(i)%accel(:) = 0.0d0   ! reset to start over
   enddo

   ! Choose which MD propagator to use:
   select case(numpar%MD_algo)
   !00000000000000000000000000000000000000000
   case default  ! velocity Verlet (2d order):
      ! Atomic Verlet step:
      call make_time_step_atoms(Scell, matter, numpar, 2)    ! module "Atomic_tools"
      ! Supercell Verlet step:
      call make_time_step_supercell(Scell, matter, numpar, 2) ! supercall Verlet step, module "Atomic_tools"
      ! Update Hamiltonian after the atomic and supercell motion:
      call get_Hamilonian_and_E(Scell, numpar, matter, 2, Err, time) ! below

      ! Save numerical acceleration:
!       do i = 1, Nat
!          Scell(1)%MDatoms(i)%accel(:) = (Scell(1)%MDatoms(i)%V(:) - V0(i,:)) / numpar%halfdt ! [A/fs^2]
!       enddo
   !11111111111111111111111111111111111111111
   case (1)  ! Yoshida (4th order)
      ! First step of Yoshida for atomic coordinates and for supercell vectors:
      call make_time_step_atoms_Y4(Scell, matter, numpar, 1, 1)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 1, 1) ! supercall Verlet step, module "Atomic_tools"
      ! Update forces/accelerations after the first coordinates step:
      call get_Hamilonian_and_E(Scell, numpar, matter, 2, Err, time) ! module "TB"
      ! First step of Yoshida for velosities:
      call make_time_step_atoms_Y4(Scell, matter, numpar, 1, 2)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 1, 2) ! supercall Verlet step, module "Atomic_tools"

      ! Second step of Yoshida
      call make_time_step_atoms_Y4(Scell, matter, numpar, 2, 1)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 2, 1) ! supercall Verlet step, module "Atomic_tools"
      ! Update forces/accelerations after the first coordinates step:
      call get_Hamilonian_and_E(Scell, numpar, matter, 2, Err, time) ! module "TB"
      ! Second step of Yoshida for velosities:
      call make_time_step_atoms_Y4(Scell, matter, numpar, 2, 2)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 2, 2) ! supercall Verlet step, module "Atomic_tools"

      ! Third step of Yoshida
      call make_time_step_atoms_Y4(Scell, matter, numpar, 3, 1)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 3, 1) ! supercall Verlet step, module "Atomic_tools"
      ! Update forces/accelerations after the first coordinates step:
      call get_Hamilonian_and_E(Scell, numpar, matter, 2, Err, time) ! module "TB"
      ! Third step of Yoshida for velosities:
      call make_time_step_atoms_Y4(Scell, matter, numpar, 3, 2)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 3, 2) ! supercall Verlet step, module "Atomic_tools"

      ! Fourth step of Yoshida
      call make_time_step_atoms_Y4(Scell, matter, numpar, 4, 1)    ! module "Atomic_tools"
      call make_time_step_supercell_Y4(Scell, matter, numpar, 4, 1) ! supercall Verlet step, module "Atomic_tools"
      ! Update forces/accelerations after the first coordinates step:
      call get_Hamilonian_and_E(Scell, numpar, matter, 2, Err, time) ! module "TB"
      ! Fourth step of Yoshida for velosities is absent, V4=V3.

      ! Save numerical acceleration:
!       do i = 1, Nat
!          Scell(1)%MDatoms(i)%accel(:) = (Scell(1)%MDatoms(i)%V(:) - V0(i,:)) / numpar%dt ! [A/fs^2]
!       enddo
   !22222222222222222222222222222222222222222
   case (2)  ! Martyna algorithm (4th order):
      ! a) New coordinate:
      call make_time_step_atoms_M(Scell, matter, numpar, 1)    ! module "Atomic_tools"
      ! b) New potential:
      call get_Hamilonian_and_E(Scell, numpar, matter, 2, Err, time) ! module "TB"
      ! c) New velocity:
      call make_time_step_atoms_M(Scell, matter, numpar, 2)    ! module "Atomic_tools"
      ! d) New effective force:
      call make_time_step_atoms_M(Scell, matter, numpar, 3)    ! module "Atomic_tools"
      ! e) New effective force velocity:
      call make_time_step_atoms_M(Scell, matter, numpar, 4)    ! module "Atomic_tools"
      ! f) New effective force acceleration:
      call make_time_step_atoms_M(Scell, matter, numpar, 5)    ! module "Atomic_tools"

      ! For Supercell, use Verlet:
      call make_time_step_supercell(Scell, matter, numpar, 2) ! supercall Verlet step, module "Atomic_tools"

      ! Save numerical acceleration:
!       do i = 1, Nat
!          Scell(1)%MDatoms(i)%accel(:) = (Scell(1)%MDatoms(i)%V(:) - V0(i,:)) / numpar%dt ! [A/fs^2]
!       enddo
   endselect

end subroutine MD_step


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
   real(8), dimension(:,:), allocatable :: HperS   ! H_1/S for SCC calculations

   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
      ! Create and diaganalize Hamiltonain (in Pettifor form):
      SC:do NSC = 1, size(Scell) ! for all supercells

         ! Get lists of nearest neighbors:
         call get_near_neighbours(Scell, numpar) ! see "Atomic_tools"
         
         ! Get and save in matrices often-used elements for forces calculations:
         ! Here, get the matrix of directional cosines and associated expressions:
         call Construct_M_x1(Scell, NSC, M_x1, M_xrr, M_lmn) ! see below 
         
         ! Electronic TB Hamiltonian part:
         ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types

            ! Create and diagonalize TB Hamiltonian:
            call create_and_diagonalize_H(Scell, NSC, numpar, matter, ARRAY, which_fe, t, &
                                    M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Mjs, M_Lag_exp, M_d_Lag_exp, &
                                    M_E0ij, M_dE0ij, HperS, Err)   ! Below

            ! Get the matrix of coefficients used for calculation of forces:
            call Construct_Aij(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Aij) ! see below

            ! Construct forces:
            if (numpar%do_atoms) then ! atoms are allowed to be moving:
               ! Attractive part:
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
                  ! If we need SCC, get the charges:
                  if (numpar%scc) then ! include SCC term
                     call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei, HperS) ! see below
                  else ! no SCC term
                     call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  endif
                  ! Get the derivatives of the Hamiltonian:
                  call get_dHij_drij_DFTB(numpar, Scell, NSC, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_DFTB"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_DFTB(Scell, NSC, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei) ! module "TB_DFTB"

               type is (TB_H_3TB)
                  ! Get the energy-weighted density matrix:
                  if (numpar%scc) then ! include SCC term
                     call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei_scc_part, M_Aij_x_Ei, HperS) ! below [B 0]
                     !call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei, HperS) ! below [B 1]
                  else ! no SCC term
                     call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
                  endif
                  ! Get the derivatives of the Hamiltonian:
                  call get_dHij_drij_3TB(numpar, Scell, NSC, ARRAY, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, &
                        M_lmn, M_Aij_x_Ei, Mjs, M_Lag_exp, M_d_Lag_exp) ! module "TB_3TB"
                  ! Get attractive forces for supercell from the derivatives of the Hamiltonian:
                  call Attract_TB_Forces_Press_3TB(Scell, NSC, ARRAY, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, &
                        M_lmn, M_Aij_x_Ei, Mjs, M_Lag_exp, M_d_Lag_exp) ! module "TB_3TB"

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
            type is (TB_Rep_DFTB_no) ! TB parametrization according to DFTB
               ! Get repulsive forces acting on all atoms:
               call dErdr_s_DFTB_no(ARRAY2, Scell, NSC) ! derivatives of the repulsive energy by s; module "TB_DFTB"
               ! Get repulsive forces acting on the supercell:
               call dErdr_Pressure_s_DFTB_no(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_DFTB"
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

         ! SCC Coulomb contribution:
         call Coulomb_force_from_SCC(numpar, matter, Scell, NSC) ! see below

         ! van der Waals part with TB Hamiltonian:
         call vdW_forces(Scell(NSC)%TB_Waals, Scell, NSC, numpar) ! get all the van der Waals forces, see below

         ! Coulomb potential part for modelling Coulomb explosion of a finite system:
         call Coulomb_forces(Scell(NSC)%TB_Coul, Scell, NSC, numpar) ! get Coulomb forces, see below
         
         ! Exponential wall potential part:
         call Exponential_wall_forces(Scell(NSC)%TB_Expwall, Scell, NSC, matter, numpar) ! get Exponential wall forces, see below
         !cccccccccccccccccccccccccccccccccccccccccccccc


         ! Get total forces for all atoms:
         call total_forces(Scell(NSC)%MDatoms) ! sum up repulsive and attractive forces, module "Atomic_tools"

         ! Get forces for the supercell:
         call Potential_super_cell_forces(numpar, Scell, NSC, matter)  ! module "Atomic_tools"

         ! Adding kinetic part and pressure to super-cell forces:
         call super_cell_forces(numpar, Scell, NSC, matter, Scell(NSC)%SCforce) ! module "Atomic_tools"

         ! Get new volume of the supercell:
         call Det_3x3(Scell(NSC)%supce, Scell(NSC)%V) ! finding initial volume of the super-cell, module "Algebra_tools"

      enddo SC
   endif DO_TB
   
   ! Get new energies of the system: potential, kinetic, global:
   call get_new_energies(Scell, matter, numpar, t, Err) ! below

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
   if (allocated(M_dE0ij)) deallocate(M_dE0ij)
   if (allocated(M_Lag_exp)) deallocate(M_Lag_exp)
   if (allocated(M_d_Lag_exp)) deallocate(M_d_Lag_exp)
   if (allocated(Mjs)) deallocate(Mjs)
   if (allocated(M_drij_dsk)) deallocate(M_drij_dsk)
   if (allocated(M_dlmn)) deallocate(M_dlmn)
   if (allocated(HperS)) deallocate(HperS)
end subroutine get_Hamilonian_and_E




subroutine create_and_diagonalize_H(Scell, NSC, numpar, matter, TB_Hamil, which_fe, t, &
                                    M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Mjs, M_Lag_exp, M_d_Lag_exp, &
                                    M_E0ij, M_dE0ij, HperS, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC
   type(Numerics_param), intent(inout) :: numpar   ! all numerical parameters
   type(Solid), intent(inout) :: matter ! material parameters
   class(TB_Hamiltonian), dimension(:,:), intent(in) :: TB_Hamil
   integer, intent(in) :: which_fe ! which method is used to get electron distribution
   real(8), intent(in) :: t ! [fs] timestep
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Vij    ! Overlap functions for H, all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dVij   ! derivatives of Overlap functions for H orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_SVij   ! Overlap matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dSVij  ! derivatives of Overlap matrix,  all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_E0ij   ! on-site energies (used for BOP)
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dE0ij     ! derivatives on-site energies (used for BOP)
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Lag_exp   ! Laguerres for 3-body radial funcs (for 3TB)
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_d_Lag_exp ! derivatives of Laguerres (for 3TB)
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_lmn    ! directional cosines l, m, n
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Mjs      ! K-S part of overlaps with s-orb. (for 3TB)
   real(8), dimension(:,:), allocatable, intent(inout) :: HperS   ! H_1/S
   type(Error_handling), intent(inout) :: Err	! error save
   !------------------------
   real(8), dimension(:,:), allocatable :: H_scc_0, H_scc_1
   real(8), dimension(:,:), allocatable :: gam_ij
   real(8), dimension(size(matter%Atoms)) :: q, q0, q00, q1
   real(8) :: eps, iteralpha, E_coul
   integer :: i, Niter


   select type(TB_Hamil)
      type is (TB_H_Pettifor) ! TB parametrization according to Pettifor
         call construct_TB_H_Pettifor(numpar, matter, TB_Hamil, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Pettifor"
      type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni
         call construct_TB_H_Molteni(numpar, matter, TB_Hamil, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Molteni"
      type is (TB_H_Fu)  ! TB parametrization accroding to Fu
         call construct_TB_H_Fu(numpar, matter, TB_Hamil, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Fu"
      type is (TB_H_NRL)  ! TB parametrization accroding to NRL
         call Construct_Vij_NRL(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_NRL"
         call construct_TB_H_NRL(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_NRL"
      type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB
         call Construct_Vij_DFTB(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_DFTB"
         call construct_TB_H_DFTB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_DFTB"

      type is (TB_H_3TB)  ! TB parametrization accroding to 3TB
         ! Get the overlaps between orbitals and ficticios s orbital (for 3-body parts):
         call get_Mjs_factors(numpar%N_basis_size, Scell(NSC), M_lmn, Mjs)   ! module "TB_3TB"
         ! Get the overlaps and reusable functions:
         call Construct_Vij_3TB(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_Lag_exp, M_d_Lag_exp) ! module "TB_3TB"
         ! Construct the Hamiltonian, diagonalize it, get the energy:
         call construct_TB_H_3TB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Scell, NSC, Err) ! module "TB_3TB"

      type is (TB_H_BOP)  ! TB parametrization accroding to BOP (incomplete)
         call Construct_Vij_BOP(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_E0ij, M_dE0ij)    ! module "TB_BOP"
         call construct_TB_H_BOP(numpar, TB_Hamil, matter, M_Vij, M_SVij, M_E0ij, M_lmn, Scell, NSC, Err)    ! module "TB_BOP"
      type is (TB_H_xTB)  ! TB parametrization accroding to xTB (NOT READY)
!         call Construct_Vij_xTB(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_xTB"
!         call construct_TB_H_xTB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_xTB"
   end select
   if (numpar%verbose) print*, 'Hamiltonian constructed and diagonalized succesfully'

   ! Get the DOS weights for each energy level, if required:
   if (allocated(Scell(NSC)%Sij)) then ! nonorthogonal Hamiltonian
      call get_DOS_weights(1, numpar%mask_DOS, numpar%DOS_weights, Hij=Scell(NSC)%Ha, &
                           Sij=Scell(NSC)%Sij, eigen_S=Scell(NSC)%eigen_S) ! below
   else ! orthogonal Hamiltonian
      call get_DOS_weights(1, numpar%mask_DOS, numpar%DOS_weights, Hij=Scell(NSC)%Ha) ! below
   endif

   ! Fill corresponding energy levels + repulsive energy cotribution:
   select case (which_fe)
      case (1) ! distribution for given Te (or provided in the file):
         if (numpar%fe_input_exists) then ! distribution provieded by the user:
            Scell(NSC)%fe = 0.0d0 ! to start with
            Scell(NSC)%fe(1 : min( size(numpar%fe_input),size(Scell(NSC)%fe) ) ) = &
                 numpar%fe_input(1 : min( size(numpar%fe_input),size(Scell(NSC)%fe) ) )
            deallocate(numpar%fe_input)   ! no longer needed
         else ! given temperature, assuming Fermi distribution:
            call set_initial_fe(Scell, matter, Err) ! module "Electron_tools"
            call get_new_global_energy(Scell(NSC), Scell(NSC)%nrg) ! module "Electron_tools"
            if (numpar%scc) call get_Mulliken(numpar%Mulliken_model, numpar%mask_DOS, numpar%DOS_weights, Scell(NSC)%Ha, &
                     Scell(NSC)%fe, matter, Scell(NSC)%MDAtoms, matter%Atoms(:)%mulliken_Ne, matter%Atoms(:)%mulliken_q) ! below
         endif
      case (2) ! distribution for given Ee + repulsive energy:
         call get_new_energies(Scell, matter, numpar, t, Err) ! below
   end select

   ! If we need SCC, get the charges:
   if (numpar%scc) then
      ! Save Hamiltonian prior to scc cycles (H_0 term):
      if (.not.allocated(H_scc_0)) allocate(H_scc_0(size(Scell(NSC)%Ha,1),size(Scell(NSC)%Ha,2)))
      H_scc_0 = Scell(NSC)%H_non ! save the H_0 (non-SCC) part of the Hamiltonin
      if (.not.allocated(H_scc_1)) allocate(H_scc_1(size(Scell(NSC)%Ha,1),size(Scell(NSC)%Ha,2)),source = 0.0d0)


      iteralpha = numpar%scc_mix ! mixing coefficient in SCC
      Niter = 250    ! maximal number of scc iterations
      i = 0          ! to start counting iterations

      ! Define deviations of Mulliken charges from atomic ones:
      ! Define first Mulliken charges, if SCC is required:
      !q(:) = matter%Atoms(:)%NVB - matter%Atoms(:)%mulliken_Ne ! Mulliken charge of all elements [M 0]
      q(:) = matter%Atoms(:)%mulliken_q ! Mulliken charge of all elements [M 0]

      q0 = 0.0d0     ! to start with
      q00 = 1.0d10    ! to start with
      E_coul = 1.0d10 ! to start with
      Scell(NSC)%nrg%E_coul_scc = 0.0d0   ! to start with
      if (numpar%verbose) write(*,'(a,f10.5,f10.5)') ' Initial Mulliken charges calculated ', q(:)

      ! Get the gamma parameters for each pair of atoms:
      call get_all_gam_ij(Scell(NSC), matter, numpar, gam_ij) ! below
      if (numpar%verbose) print*, 'SCC gamma parameter obtained succesfully'

      ! SCC iterations:
      eps = 1.0d-6   ! precision
      SCC_ITER:do while ( abs(maxval(q(:) - q0(:))) > eps*abs(maxval(q(:))) )
      !SCC_ITER:do while ( abs(E_coul - Scell(NSC)%nrg%E_coul_scc)/max(abs(E_coul),abs(Scell(NSC)%nrg%E_coul_scc)) > eps )

         i = i + 1   ! count iterations
         if (i > Niter) then
            print*, ' Problem in SCC: number of iterations exceeds the limit: ', Niter
            print*, ' Charges 2-before-last iteration:', q00(:)
            print*, ' Charges before last iteration  :', q0(:)
            print*, ' Charges on the last iteration  :', q(:)
!             print*, ' Coulomb energy before last iteration:', E_coul
!             print*, ' Coulomb energy at the last iteration:', Scell(NSC)%nrg%E_coul_scc
            exit SCC_ITER  ! too many iterations
         endif

         ! Get the H_1_ij/S_ij:
         call get_HperS(Scell(NSC), numpar, gam_ij, q, HperS)  ! below

         ! Construct scc corrections (H_1 term):
         call create_second_order_scc_term_H(Scell(NSC), matter, numpar, Scell(NSC)%Sij, HperS, H_scc_1) ! below

         ! Get the second-order Hamiltonian, diagonalize it:
         select type(TB_Hamil)
            type is (TB_H_Pettifor) ! TB parametrization according to Pettifor
               !call construct_TB_H_Pettifor(numpar, matter, TB_Hamil, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Pettifor"
               numpar%scc = .false.
               print*, 'SCC calculations are not supportd with this TB parameterization'
               print*, 'Continue at zero-order non-SCC calculations'
            type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni
               !call construct_TB_H_Molteni(numpar, matter, TB_Hamil, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Molteni"
               numpar%scc = .false.
               print*, 'SCC calculations are not supportd with this TB parameterization'
               print*, 'Continue at zero-order non-SCC calculations'
            type is (TB_H_Fu)  ! TB parametrization accroding to Fu
               !call construct_TB_H_Fu(numpar, matter, TB_Hamil, Scell, NSC, Scell(NSC)%Ha, Err) ! module "TB_Fu"
               numpar%scc = .false.
               print*, 'SCC calculations are not supportd with this TB parameterization'
               print*, 'Continue at zero-order non-SCC calculations'
            type is (TB_H_NRL)  ! TB parametrization accroding to NRL
               !call construct_TB_H_NRL(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_NRL"
               numpar%scc = .false.
               print*, 'SCC calculations are not supportd with this TB parameterization'
               print*, 'Continue at zero-order non-SCC calculations'
            type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB
               ! Construct the Hamiltonian, diagonalize it, get the energy:
               call construct_TB_H_DFTB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err, &
                  scc=1, H_scc_0=H_scc_0, H_scc_1=H_scc_1) ! module "TB_DFTB"
            type is (TB_H_3TB)  ! TB parametrization accroding to 3TB
               ! Construct the Hamiltonian, diagonalize it, get the energy:
               call construct_TB_H_3TB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, &
                  Scell, NSC, Err, scc=1, H_scc_0=H_scc_0, H_scc_1=H_scc_1) ! module "TB_3TB"
            type is (TB_H_BOP)  ! TB parametrization accroding to BOP (incomplete)
               !call construct_TB_H_BOP(numpar, TB_Hamil, matter, M_Vij, M_SVij, M_E0ij, M_lmn, Scell, NSC, Err)    ! module "TB_BOP"
               numpar%scc = .false.
               print*, 'SCC calculations are not supportd with this TB parameterization'
               print*, 'Continue at zero-order non-SCC calculations'
            type is (TB_H_xTB)  ! TB parametrization accroding to xTB (NOT READY)
!              call Construct_Vij_xTB(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_xTB"
!              call construct_TB_H_xTB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err) ! module "TB_xTB"
               numpar%scc = .false.
               print*, 'SCC calculations are not supportd with this TB parameterization'
               print*, 'Continue at zero-order non-SCC calculations'
         end select

         ! Get the DOS weights for each energy level, if required:
         if (allocated(Scell(NSC)%Sij)) then
            call get_DOS_weights(1, numpar%mask_DOS, numpar%DOS_weights, Hij=Scell(NSC)%Ha, &
                           Sij=Scell(NSC)%Sij, eigen_S=Scell(NSC)%eigen_S) ! below
         else
            call get_DOS_weights(1, numpar%mask_DOS, numpar%DOS_weights, Hij=Scell(NSC)%Ha) ! below
         endif

         ! Update deviations of Mulliken charges from atomic ones:
         call get_Mulliken(numpar%Mulliken_model, numpar%mask_DOS, numpar%DOS_weights, Scell(NSC)%Ha, &
                            Scell(NSC)%fe, matter, Scell(NSC)%MDAtoms, matter%Atoms(:)%mulliken_Ne, q1) ! below
         q00 = q0 ! save for the next step of scc cycle
         q0 = q   ! save for the next step of scc cycle
         ! Mixing of the fraction of old and new charges:
         q(:) = q0(:)*(1.0d0 - iteralpha) + q1(:)*iteralpha


         ! Get the Coulomb contribution to energy from the charge redistribution:
         !E_coul = Scell(NSC)%nrg%E_coul_scc  ! save for the next iteration
         !call get_Coulomb_scc_energy(Scell, NSC, matter, gam_ij, Scell(NSC)%nrg%E_coul_scc)   ! below
         !Scell(NSC)%nrg%E_coul_scc = 0.0d0 ! [E 1]

         ! Check if we are moving towards convergence, or away from it:
         if ( ( abs(q(1) - q0(1)) > abs(q0(1) - q00(1)) ) .or. &
              ( abs(q(1) - q00(1)) < abs(q(1) - q0(1)) ) ) then
            iteralpha = max(0.01d0, iteralpha * 0.75d0)
            if (numpar%verbose) print*, 'SCC convergence may not be reached, reducing mixing:', iteralpha
         endif
      enddo SCC_ITER

      if (numpar%verbose) then
         print*, 'SCC cycles completed with iteration #', i
         print*, ' Charges before last iteration:', q0(:)
         print*, ' Charges on the last iteration:', q(:)
      endif

      ! Save Mulliken charges for the next time-step:
      matter%Atoms(:)%mulliken_q = q(:)

      ! Restore the unpertorbed H_0 (non-SCC) part of the Hamiltonin:
      Scell(NSC)%H_non = H_scc_0
      ! Save eigenvalues of the total Hamiltonian:
      Scell(NSC)%Ei_scc_part = Scell(NSC)%Ei

      ! Update band gap:
      call find_band_gap(Scell(NSC)%Ei_scc_part, Scell(NSC), matter, numpar) ! module "Electron_tools"

      ! Update the energy in the electronic system:
      ! 1) Update the band structure, using the wave-functions corrected after the SCC procidure:
      !call get_eigenvalues_from_eigenvectors(Scell(NSC)%H_non, Scell(NSC)%Ha, Scell(NSC)%Ei_scc_part) ! module "Algebra_tools"
      call get_eigenvalues_from_eigenvectors(Scell(NSC)%H_non, Scell(NSC)%Ha, Scell(NSC)%Ei) ! module "Algebra_tools"

      ! 2) Get the Coulomb contribution to energy from the charge redistribution:
      call get_Coulomb_scc_energy(Scell, NSC, matter, gam_ij, Scell(NSC)%nrg%E_coul_scc)   ! below

      ! Fill corresponding energy levels + repulsive energy cotribution:
      select case (which_fe)
         case (1) ! distribution for given Te:
            call set_initial_fe(Scell, matter, Err) ! module "Electron_tools"
            call get_new_global_energy(Scell(NSC), Scell(NSC)%nrg) ! module "Electron_tools"
         case (2) ! distribution for given Ee + repulsive energy:
            call get_new_energies(Scell, matter, numpar, t, Err) ! below
      end select

      ! Clear up:
      deallocate(H_scc_0, H_scc_1)
   else
      Scell(NSC)%nrg%E_coul_scc = 0.0d0 ! no energy associated with charge redistribution
      Scell(NSC)%Ei_scc_part = Scell(NSC)%Ei ! no difference between eigenstates and SCC-eigenstates
   endif
end subroutine create_and_diagonalize_H


subroutine get_Coulomb_scc_energy(Scell, NSC, matter, gam_ij, E_coulomb)   ! Coulomb energy
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(solid), intent(in), target :: matter   ! material parameters
   real(8), dimension(:,:), intent(in) :: gam_ij  ! effective energy values [eV]
   real(8), intent(out) :: E_coulomb  ! Total Coulomb energy of all atoms [eV]
   !=====================================================
   real(8) :: sum_a, Coul_pot, q(size(matter%Atoms))
   integer :: j, nat, atom_2, i
   integer, pointer :: m, KOA1, KOA2

   nat = Scell(NSC)%Na  ! number of atoms
   sum_a = 0.0d0        ! to start with
   E_coulomb = 0.0d0    ! to start with

   !q(:) = matter%Atoms(:)%NVB - matter%Atoms(:)%mulliken_Ne ! Mulliken charges of all elements [M 0]
   q(:) = matter%Atoms(:)%mulliken_q
!    q(:) = -q(:)

   !$omp PARALLEL private( j, KOA1, m, atom_2, i, KOA2, Coul_pot )
   !$omp do reduction( + : sum_a)
   do j = 1, nat  ! atom #1
      KOA1 => Scell(NSC)%MDatoms(j)%KOA   ! kind of atom #1
      m => Scell(NSC)%Near_neighbor_size(j)  ! number of nearest neighbors of atom #1
      do atom_2 = 1, m  ! only for atoms close to #1
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! atom #2
         KOA2 => Scell(NSC)%MDatoms(i)%KOA   ! kind of atom #2

         ! Energy:
         Coul_pot = gam_ij(j,i) * q(KOA1) * q(KOA2)

         ! Add to total:
         sum_a = sum_a + Coul_pot
      enddo ! atom_2
   enddo ! j
   !$omp end do
   !$omp end parallel
   ! Total Coulomb energy, excluding double-counting:
   !E_coulomb = sum_a * 0.5d0   ! [eV] ! [A 0]
   E_coulomb = -sum_a * 0.25d0   ! [eV] ! [A 1]

!    print*, 'E_coulomb', E_coulomb/dble(nat)
   nullify(KOA1, KOA2, m)
end subroutine get_Coulomb_scc_energy



subroutine Coulomb_force_from_SCC(numpar, matter, Scell, NSC) ! see below
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC
   type(Solid), intent(in) :: matter ! material parameters
   !-------------------------
   if (numpar%scc) then ! calculate SCC Coulomb contribution
      ! Forces for all atoms:
      call Coulomb_force_from_SCC_s(Scell, NSC, matter, numpar)   ! below
      ! Forces for the super-cell:
      call Coulomb_force_from_SCC_Pressure_s(Scell, NSC, matter, numpar)   ! below
   endif
end subroutine Coulomb_force_from_SCC




subroutine Coulomb_force_from_SCC_s(Scell, NSC, matter, numpar)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC
   type(Solid), intent(in) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   !---------------------------------------
   real(8), dimension(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, r1, x0, y0, z0, a, b, ddlta, b_delta, q(size(matter%Atoms)), r_cut, alpha
   integer i, j, k, ik, i1, ian, dik, djk, n, atom_2
   integer, pointer :: KOA1, KOA2, m, j1
   real(8), pointer ::  x, y, z
   n = Scell(NSC)%Na ! number of atoms

   ! Get the cot off distance
   r_cut = cut_off_distance(Scell(NSC)) ! module "Coulomb"
   alpha = 3.0d0/(4.0d0*r_cut) ! Wolf's parameter chosen according to optimal value from [4]

   ! Get the charges for all different elements in the material:
   q(:) = matter%Atoms(:)%mulliken_q

   !$omp PARALLEL private(ian, m, KOA1, dpsi, atom_2, j1, KOA2, x, y, z, x1, a_r, b)
   !$omp DO
   do ian = 1, n  ! Forces for all atoms
      m => Scell(NSC)%Near_neighbor_size(ian)
      KOA1 => Scell(NSC)%MDatoms(ian)%KOA
      dpsi = 0.0d0
      do atom_2 = 1, m   ! contribution from neighboring atoms
         j1 => Scell(NSC)%Near_neighbor_list(ian, atom_2)	! this is the list of such close atoms
         KOA2 => Scell(NSC)%MDatoms(j1)%KOA

         a_r = Scell(NSC)%Near_neighbor_dist(ian,atom_2,4) ! at this distance, R
         x  => Scell(NSC)%Near_neighbor_dist(ian,atom_2,1) ! at this distance, X
         y  => Scell(NSC)%Near_neighbor_dist(ian,atom_2,2) ! at this distance, Y
         z  => Scell(NSC)%Near_neighbor_dist(ian,atom_2,3) ! at this distance, Z

         x1(:) = x*Scell(NSC)%supce(:,1) + y*Scell(NSC)%supce(:,2) + z*Scell(NSC)%supce(:,3)

         call d_get_gamma_scc(matter%Atoms(KOA1)%Hubbard_U, matter%Atoms(KOA2)%Hubbard_U, &
                                 a_r, ian, j1, b, r_cut, alpha, numpar%scc_gam_ind)  ! below
         ! include charges to construct the projections of the forces:
         b = b*q(KOA1)*q(KOA2)
         dpsi(:) = dpsi(:) + b*x1(:)/a_r
      enddo ! atom_2

      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + dpsi(:)  ! [F 0]

   enddo ! ian
   !$omp end do
   !$omp end parallel

   nullify(j1, m, KOA1, KOA2, x, y, z)
end subroutine Coulomb_force_from_SCC_s



! Derivatives of the SCC-Coulomb by h:
subroutine Coulomb_force_from_SCC_Pressure_s(Scell, NSC, matter, numpar)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !===============================================
   real(8), dimension(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i, k, l, n, atom_2
   integer, pointer :: KOA1, KOA2, m, j
   real(8) r, rcur(3), scur(3), PForce(3,3)
   real(8) df_psy, psi, dpsy, r_cut, alpha, q(size(matter%Atoms))

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = Scell(NSC)%Na ! number of atoms
      r_cut = cut_off_distance(Scell(NSC)) ! module "Coulomb"
      alpha = 3.0d0/(4.0d0*r_cut) ! Wolf's parameter chosen according to optimal value from [4]

      ! Get the charges for all different elements in the material:
      !q(:) = matter%Atoms(:)%NVB - matter%Atoms(:)%mulliken_Ne ! Mulliken charge [M 0]
      q(:) = matter%Atoms(:)%mulliken_q
!       q(:) = -q(:)

      PForce = 0.0d0 ! to start with
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         KOA1 => Scell(NSC)%MDatoms(i)%KOA
         m => Scell(NSC)%Near_neighbor_size(i)
         do atom_2 = 1,m		! do only for atoms close to that one
            j => Scell(NSC)%Near_neighbor_list(i, atom_2)	! this is the list of such close atoms
            if (j > 0) then
               KOA2 => Scell(NSC)%MDatoms(j)%KOA
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
               call d_get_gamma_scc(matter%Atoms(KOA1)%Hubbard_U, matter%Atoms(KOA2)%Hubbard_U, &
                                 r, i, j, dpsy, r_cut, alpha, numpar%scc_gam_ind)  ! below
               ! include the charges:
               dpsy = dpsy*q(KOA1)*q(KOA2)

               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! (j > 0)
         enddo ! atom_2

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)*0.5d0   ! factor 0.5 to compensate for double-counting
            enddo ! l
         enddo ! k
      enddo ! i
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add "exponential wall" part to existing TB part
   endif
   nullify(KOA1, KOA2, m, j)
end subroutine Coulomb_force_from_SCC_Pressure_s





subroutine create_second_order_scc_term_H(Scell, matter, numpar, Sij, HperS, H_scc_1)
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), dimension(:,:), intent(in) :: Sij, HperS  ! overlap materix, gamma parameteres
   real(8), dimension(:,:), intent(inout) :: H_scc_1  ! second-order scc correction to Hamiltonian
   !------------------
   integer, pointer :: nat, m
   integer :: j, atom_2, i, j1, l, i1, k, n_orb

   nat => Scell%Na	! number of atoms in the supercell
   n_orb = identify_DFTB_orbitals_per_atom(numpar%N_basis_size)  ! module "TB_DFTB"

!$omp parallel private(j, m, atom_2, i, j1, l, i1, k)
!$omp do
   do j = 1,nat	! atom #1
      m => Scell%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one
      !do atom_2 = 1,nat ! for all atoms, no cut off

         if (atom_2 == 0) then ! the same atom
            i = j    ! atom #2 = atom #1, onsite
         else  ! different atoms
            i = Scell%Near_neighbor_list(j,atom_2) ! atom #2
            !i = atom_2 ! atom #2 - no cut off
         endif

         IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
            do j1 = 1,n_orb ! all orbitals of atom #1
               l = (j-1)*n_orb+j1   ! atom #1 (j)
               do i1 = 1,n_orb ! all orbitals of atom #2
                  k = (i-1)*n_orb+i1   ! atom #2 (i)
                  ! We fill the upper triangle here
                  H_scc_1(l,k) = Sij(l,k) * HperS(l,k)
               enddo ! i1
            enddo ! j1
         endif IJ
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel

   ! b) Construct lower triangle - use symmetry:
!$omp parallel
!$omp do  private(j, m, atom_2, i, j1, l, i1, k)
   do j = 1,nat	! all atoms
      m => Scell%Near_neighbor_size(j)
      do atom_2 = 1,m ! do only for atoms close to that one
      !do atom_2 = 1,nat ! for all atoms, no cut off
         i = Scell%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         !i = atom_2 ! this is the list of such close atoms - no cut off
         if (i < j) then ! lower triangle
            do j1 = 1,n_orb ! all orbitals of atom #1
               l = (j-1)*n_orb+j1   ! atom #1
               do i1 = 1,n_orb ! all orbitals of atom #2
                  k = (i-1)*n_orb+i1   ! atom #2
                  H_scc_1(l,k) = H_scc_1(k,l)
               enddo ! i1
            enddo ! j1
         endif
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel

   nullify(nat, m)
end subroutine create_second_order_scc_term_H



subroutine get_HperS(Scell, numpar, gam_ij, q, HperS)
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), dimension(:,:), intent(in) :: gam_ij   ! gamma parameters
   real(8), dimension(:), intent(in) :: q   ! Mulliken charges (deviation from neutral)
   real(8), dimension(:,:), allocatable, intent(inout) :: HperS
   !------------------
   integer, pointer :: nat, m, m2, KOA1, KOA2, KOA3
   integer :: j, atom_2, i, j1, l, i1, k, atom_3, n_orb, n
   real(8) :: H_ij, H_ij_1, H_ij_2

   if (.not.allocated(HperS)) allocate(HperS(size(Scell%Ha,1),size(Scell%Ha,2)))
   HperS = 0.0d0  ! to start with

   nat => Scell%Na	! number of atoms in the supercell
   n_orb = identify_DFTB_orbitals_per_atom(numpar%N_basis_size)  ! module "TB_DFTB"

!$omp parallel private(j, m, atom_2, i, KOA1, KOA2, KOA3, m2, j1, l, i1, k, n, atom_3, H_ij_1, H_ij_2)
!$omp do
   do j = 1,nat	! atom #1
      KOA1 => Scell%MDatoms(j)%KOA   ! atom #1
      m => Scell%Near_neighbor_size(j)

      H_ij_1 = 0.0d0   ! to start with
      do atom_3 = 0, m  ! sum up interactions between atom #1 and atom #3
      !do atom_3 = 1, nat  ! sum up interactions between atom #1 and atom #3 - no cut off
         if (atom_3 == 0) then ! the same atom
            n = j    ! atom #3 = atom #1, onsite
         else  ! different atoms
            n = Scell%Near_neighbor_list(j,atom_3) ! index of atom #3
            !n = atom_3 ! index of atom #3 - no cut off
         endif
         KOA3 => Scell%MDatoms(n)%KOA   ! type of atom #3
         ! Construct a part of the second-order correction to the Hamiltonian:
         H_ij_1 = H_ij_1 + gam_ij(j,n) * q(KOA3)  ! First half of the last term in Eq.(10) [1]
      enddo ! atom_3

      do atom_2 = 0,m ! do only for atoms close to that one
      !do atom_2 = 1,nat ! do only for atoms close to that one - no cut off
         if (atom_2 == 0) then ! the same atom
            i = j    ! atom #2 = atom #1, onsite
         else  ! different atoms
            i = Scell%Near_neighbor_list(j,atom_2) ! atom #2
            !i = atom_2 ! atom #2 - no cut off
         endif

         IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
            KOA2 => Scell%MDatoms(i)%KOA   ! atom #2

            m2 => Scell%Near_neighbor_size(i)
            H_ij_2 = 0.0d0 ! to start with
            do atom_3 = 0, m2 ! sum up interactions between atom #2 and atom #3
            !do atom_3 = 1, nat ! sum up interactions between atom #2 and atom #3 - no cut off
               if (atom_3 == 0) then ! the same atom
                  n = i    ! atom #3 = atom #2, onsite
               else  ! different atoms
                  n = Scell%Near_neighbor_list(i,atom_3) ! index of atom #3
                  !n = atom_3 ! index of atom #3 - no cut off
               endif
               KOA3 => Scell%MDatoms(n)%KOA   ! type of atom #3
               ! Construct a part of the second-order correction to the Hamiltonian:
               H_ij_2 = H_ij_2 + gam_ij(i,n) * q(KOA3)  ! Second half of the last term in Eq.(10) [1]
            enddo ! atom_3

            ! Save it into all shells of the two atoms:
            do j1 = 1,n_orb ! all orbitals of atom #1
               l = (j-1)*n_orb+j1   ! atom #1 (j)
               do i1 = 1,n_orb ! all orbitals of atom #2
                  k = (i-1)*n_orb+i1   ! atom #2 (i)
                  ! We fill the upper triangle here
                  HperS(l,k) = 0.5d0*(H_ij_1+H_ij_2)
               enddo ! i1
            enddo ! j1

         endif IJ
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel

   ! b) Construct lower triangle - use symmetry:
!$omp parallel
!$omp do  private(j, m, atom_2, i, j1, l, i1, k)
   do j = 1,nat	! all atoms
      m => Scell%Near_neighbor_size(j)
      do atom_2 = 1,m ! do only for atoms close to that one
      !do atom_2 = 1,nat ! do only for atoms close to that one - no cut off
         i = Scell%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         !i = atom_2 ! this is the list of such close atoms - no cut off
         if (i < j) then ! lower triangle
            do j1 = 1,n_orb ! all orbitals of atom #1
               l = (j-1)*n_orb+j1   ! atom #1
               do i1 = 1,n_orb ! all orbitals of atom #2
                  k = (i-1)*n_orb+i1   ! atom #2
                  HperS(l,k) = HperS(k,l)
               enddo ! i1
            enddo ! j1
         endif
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel

   nullify(nat, m, KOA1, KOA2, KOA3)
end subroutine get_HperS




subroutine get_all_gam_ij(Scell, matter, numpar, gam_ij)
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   real(8), dimension(:,:), allocatable, intent(inout) :: gam_ij
   !------------------
   integer, pointer :: nat, m, KOA1, KOA2
   integer :: j, atom_2, i
   real(8) :: r, r_cut, alpha

   nat => Scell%Na	! number of atoms in the supercell
   if(.not.allocated(gam_ij)) allocate(gam_ij(nat,nat))
   gam_ij = 0.0d0 ! to start with

   ! Get the cot off distance
   r_cut = cut_off_distance(Scell) ! module "Coulomb"
   alpha = 3.0d0/(4.0d0*r_cut) ! Wolf's parameter chosen according to optimal value from [4]


!$omp parallel private(j, m, atom_2, i, KOA1, KOA2, r)
!$omp do
   do j = 1,nat	! atom #1
      KOA1 => Scell%MDatoms(j)%KOA   ! atom #1
      m => Scell%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one
      !do atom_2 = 1,nat ! do only for atoms close to that one - no cut off

         if (atom_2 == 0) then ! the same atom
            i = j    ! atom #2 = atom #1, onsite
            r = 0.0d0   ! same atom, no distance
         else  ! different atoms
            i = Scell%Near_neighbor_list(j,atom_2) ! atom #2
            !i = atom_2 ! atom #2 - no cut off
            r = Scell%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]
            !call shortest_distance(Scell, j, i, r) ! at this distance, R [A] - no cut off
         endif

         IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
            KOA2 => Scell%MDatoms(i)%KOA   ! atom #2

            ! Get the gamma parameter for scc calculations:
            call get_gamma_scc(matter%Atoms(KOA1)%Hubbard_U, matter%Atoms(KOA2)%Hubbard_U, &
                                 r, j, i, gam_ij(j,i), r_cut, alpha, numpar%scc_gam_ind)  !  below

         endif IJ
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel

   ! b) Construct lower triangle - use symmetry:
!$omp parallel
!$omp do  private(j, m, atom_2, i)
      do j = 1,nat	! all atoms
         m => Scell%Near_neighbor_size(j)
         do atom_2 = 1,m ! do only for atoms close to that one
         !do atom_2 = 1,nat ! do only for atoms close to that one - no cut off
            i = Scell%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            !i = atom_2 ! this is the list of such close atoms - no cut off
            if (i < j) then ! lower triangle
               gam_ij(j,i) = gam_ij(i,j)
            endif
         enddo ! j
      enddo ! i
!$omp end do
!$omp end parallel

   nullify(nat, m, KOA1, KOA2)
end subroutine get_all_gam_ij


subroutine get_gamma_scc(Ui, Uj, Rij, i, j, gam_ij, r_cut, alpha, ind_gam)
! [1]  https://arxiv.org/pdf/2112.11585.pdf
! [2]  DOI: 10.1021/jp071338j
   integer, intent(in) :: i, j   ! indices of atoms #1 and #2
   real(8), intent(in) :: Ui, Uj ! Hubbard parameteres for atoms i and j
   real(8), intent(in) :: Rij    ! [A] distance between atoms i and j
   real(8), intent(out) :: gam_ij   ! gamma_i,j
   real(8), intent(in) :: r_cut, alpha ! [A] Wolf's method cot off distance and alpha parameter
   integer, intent(in) :: ind_gam   ! which formula to use for gamma (use only option 0, others are not ready!)
   real(8) :: Cij, R0, WCoul, SCoul
   if (i == j) then  ! the same atom
      gam_ij = Ui   ! Hubbard parameters
   else  ! two different atoms
      select case (ind_gam)
      case default   ! Use Wolf method together with Eqs. (7) and (8) from [1]
         Cij = sqrt( g_half_Pi / (1.0d0 / Ui**2 + 1.0d0/Uj**2) ) * g_ev2au ! [eV] -> [a.u.]
         ! Get the Coulomb potential according to Wolf's method:
         WCoul = Coulomb_Wolf_pot(1.0d0, 1.0d0, Rij, r_cut, alpha)  ! module "Coulomb"
         SCoul = Coulomb_Wolf_self_term(1.0d0, r_cut, alpha)  ! module "Coulomb"
         ! Collect the terms:
         R0 = Rij/g_a0  ! dimensionless units
         gam_ij = erf(Cij*R0) * (WCoul + SCoul) ! [eV]

      case (-1)   ! [1] Eqs. (7) and (8)
         Cij = sqrt( g_half_Pi/(1.0d0 / Ui**2 + 1.0d0/Uj**2) ) * g_ev2au ! [eV] -> [a.u.]
         R0 = Rij/g_a0  ! dimensionless units
         gam_ij = m_k * erf(Cij*R0) / Rij ! [eV]

      case (1) ! Klopman-Ohno [2], Eq.(13)
         R0 = Rij/g_a0
         gam_ij = g_au2ev/sqrt(R0**2 + 0.25d0 *(g_au2ev/Ui + g_au2ev/Uj)**2 ) ! -> energy units [eV]

      case (2) ! Mataga-Nishimoto [2], Eq.(20)
         R0 = Rij/g_a0
         gam_ij = g_au2ev/(R0 + 0.5d0 *(g_au2ev/Ui + g_au2ev/Uj) ) ! -> energy units [eV]

      endselect
   endif
end subroutine get_gamma_scc



subroutine d_get_gamma_scc(Ui, Uj, Rij, i, j, d_gam_ij, r_cut, alpha, ind_gam)  ! d_gamma/d_Rij
! [1]  https://arxiv.org/pdf/2112.11585.pdf
! [2]  DOI: 10.1021/jp071338j
   integer, intent(in) :: i, j   ! indices of atoms #1 and #2
   real(8), intent(in) :: Ui, Uj ! Hubbard parameteres for atoms i and j
   real(8), intent(in) :: Rij    ! [A] distance between atoms i and j
   real(8), intent(out) :: d_gam_ij   ! gamma_i,j
   real(8), intent(in) :: r_cut, alpha ! [A] Wolf's method cot off distance and alpha parameter
   integer, intent(in) :: ind_gam   ! which formula to use for gamma (use only option 0, others are not ready!)
   real(8) :: Cij, R0, WCoul, dWCoul, SCoul, d_erf
   if (i == j) then  ! the same atom
      d_gam_ij = 0.0d0
   else  ! two different atoms
      select case (ind_gam)
      case default   ! Use Wolf method together with Eqs. (7) and (8) from [1]
         Cij = sqrt( g_half_Pi / (1.0d0 / Ui**2 + 1.0d0/Uj**2) ) * g_ev2au ! [eV] -> [a.u.]
         ! Get the Coulomb potential according to Wolf's method:
         WCoul = Coulomb_Wolf_pot(1.0d0, 1.0d0, Rij, r_cut, alpha)  ! module "Coulomb"
         SCoul = Coulomb_Wolf_self_term(1.0d0, r_cut, alpha)  ! module "Coulomb"
         dWCoul = d_Coulomb_Wolf_pot(1.0d0, 1.0d0, Rij, r_cut, alpha)  ! module "Coulomb"
         R0 = Rij/g_a0  ! dimensionless units
         d_erf = 2.0d0*Cij/(m_sqrtPi * g_a0) * exp(-(Cij*R0)**2)  ! [a.u./a.u.] -> [1/A]
         ! Collect the terms:
         d_gam_ij = erf(Cij*R0)*dWCoul + d_erf*(WCoul+SCoul) ! [eV/A]
         !d_gam_ij = d_gam_ij * 100.0   ! test
      case (-1)   ! [1] Eqs. (7) and (8)
         Cij = sqrt( g_half_Pi/(1.0d0 / Ui**2 + 1.0d0/Uj**2) ) * g_ev2au ! [eV] -> [a.u.]
         R0 = Rij/g_a0  ! dimensionless units
         d_erf = 2.0d0*Cij/(m_sqrtPi*g_a0) * exp(-(Cij*R0)**2) ! [a.u./a.u.] -> [1/A]
         d_gam_ij = m_k * (-erf(Cij*R0) / Rij**2 + d_erf/Rij) ! [eV/A]

      case (1) ! Klopman-Ohno [2], Eq.(13)
         R0 = Rij/g_a0
         d_gam_ij = -g_au2ev*R0/(R0**2 + 0.25d0 *(g_au2ev/Ui + g_au2ev/Uj)**2)**(1.5d0) / g_a0 ! [eV/A]

      case (2) ! Mataga-Nishimoto [2], Eq.(20)
         R0 = Rij/g_a0
         d_gam_ij = -g_au2ev/(R0 + 0.5d0 *(g_au2ev/Ui + g_au2ev/Uj) )**2 / g_a0 ! [eV/A]
      endselect
   endif
end subroutine d_get_gamma_scc





subroutine k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, UPG)
   integer, intent(in) :: schem  ! scheme for k-points samping: 0=Monkhorst-Pack; 1=user provided grid
   integer, intent(in) :: ix, iy, iz, ixm, iym, izm   ! number of grid points and sizes of the grids along the reciprocal vectors
   real(8), dimension(:,:), allocatable, intent(in) :: UPG  ! user provided grid
   real(8), intent(out) :: kx, ky, kz  ! coordinates of the k-points in the units of the reciprocal vectors
   !------------------------------
   integer :: Ngp, Nsiz
   select case (schem)
   case (1) ! user defined grid:
      UG:if (.not.allocated(UPG)) then ! use Monkhorst-Pack grid
         kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
         ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
         kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
      else UG  ! user provided the grid:
         Ngp = (ix-1)*iym*ixm + (iy-1)*ixm + iz ! number of grid point for user defined grid
         Nsiz = size(UPG,1)   ! size of the user provided grid
         if (Ngp > Nsiz) then ! If user didn't match the number of k-points to the grid provided, just use Gamma point for extra points:
            kx = 0.0d0
            ky = 0.0d0
            kz = 0.0d0
         else  ! read from the user provided array:
            kx = UPG(Ngp,1)
            ky = UPG(Ngp,2)
            kz = UPG(Ngp,3)
         endif
      endif UG
   case default   ! Monkhorst-Pack [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]:
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

   BS:if (do_cart) then ! Cartesian basis set (UNUSED FOR NOW):
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
            KOA => Scell(NSC)%MDatoms(i)%KOA  ! atom #i
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


 subroutine get_DOS_weights(ind, masks_DOS, DOS_weights, Hij, CHij, Sij, eigen_S)
   integer, intent(in) :: ind   ! which model to use
   logical, dimension(:,:,:), intent(in) :: masks_DOS   ! partial DOS made of each orbital type, if required to be constructed
   real(8), dimension(:,:,:), allocatable, intent(inout) :: DOS_weights     ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:,:), intent(in), optional :: Hij     ! real eigenvectors
   complex, dimension(:,:), intent(in), optional :: CHij    ! complex eigenvectors
   real(8), dimension(:,:), intent(in), optional :: Sij     ! Overlap for nonorthogonal H
   real(8), dimension(:), intent(in), optional :: eigen_S   ! eigenvalues of Sij (to control for linear-dependent raws)
   !-------------------------------
   real(8), dimension(size(masks_DOS,3)) :: temp_vec
   real(8) :: temp
   integer :: j, Nsiz, N_at, N_types, i_at, i_types, k
   logical :: include_Sij


   if (ind == 1) then ! get weighted DOS:
      N_at = size(masks_DOS,1)
      N_types = size(masks_DOS,2)
      Nsiz = size(masks_DOS,3)
      if (.not.allocated(DOS_weights)) allocate(DOS_weights(N_at,N_types,Nsiz))

      ! Do we need to include overlap matrix:
      include_Sij = .false.   ! to start with
      ! Check #1: the overlap must be present to do that:
      if (present(Sij)) include_Sij = .true.
      ! Check #2: numerical consistency: Sij must have positive eigenvalues,
      ! if it doesn't, something went wrong, that means,
      ! it is worse to use a wrong Sij than not to use it at all:
      if (present(eigen_S)) then
         if (ANY(eigen_S(:) < 0.0d0)) include_Sij = .false.
      endif
!       print*, 'get_DOS_weights', include_Sij

      ! Get the contributions of all shells and elements:
      ORTH:if (include_Sij) then ! non-orthogonal, with overlap matrix Sij
         !$omp PARALLEL private(j, k, temp_vec, temp, i_at, i_types)
         !$omp do
         do j = 1, Nsiz ! for all energy levels
            if (present(Hij)) then ! real H
               do k = 1, Nsiz ! for all energy levels
                  temp_vec(k) = Hij(k,j) * SUM(Hij(:,j) * Sij(k,:))
               enddo
            elseif (present(CHij)) then    ! complex H
               do k = 1, Nsiz ! for all energy levels
                  temp_vec(k) = 0.5d0*dble( conjg(CHij(k,j)) * SUM(CHij(:,j) * Sij(k,:)) + &
                                            CHij(k,j) * SUM(conjg(CHij(:,j)) * Sij(k,:)) )
               enddo
            else   ! undefined
               temp_vec(:) = 0.0d0
            endif
            temp = SUM(temp_vec)
            if (temp == 0.0d0) temp = 1.0d0 ! if it's undefined, avoid crushing
            do i_at = 1, N_at
               do i_types = 1, N_types
                  DOS_weights(i_at, i_types, j) = SUM(temp_vec(:), MASK = masks_DOS(i_at, i_types, :))/temp
               enddo ! i_types
            enddo ! i_at
         enddo !  j = 1, Nsiz
         !$omp end do
         !$omp end parallel
      else ORTH  ! orthogonal (or inconsistent Sij)
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
            if (temp == 0.0d0) temp = 1.0d0  ! if it's undefined, avoid crushing
            do i_at = 1, N_at
               do i_types = 1, N_types
                  DOS_weights(i_at, i_types, j) = SUM(temp_vec(:), MASK = masks_DOS(i_at, i_types, :))/temp
               enddo ! i_types
            enddo ! i_at
         enddo !  j = 1, Nsiz
         !$omp end do
         !$omp end parallel
      endif ORTH
   endif ! (ind == 1)
end subroutine get_DOS_weights



subroutine get_Mullikens_all(Scell, matter, numpar)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter     ! material parameters
   type (Numerics_param), intent(in) :: numpar    ! numerical parameters, including drude-function
   !-------------------------------
   ! 1) Do the avereage Mullikens:
   call get_Mulliken(numpar%Mulliken_model, numpar%mask_DOS, numpar%DOS_weights, Scell%Ha, &
         Scell%fe, matter, Scell%MDAtoms, matter%Atoms(:)%mulliken_Ne, matter%Atoms(:)%mulliken_q) ! below

   ! 2) Do the Mullikens for each atom:
   call get_Mulliken_each_atom(numpar%Mulliken_model, Scell, matter, numpar)   ! below
end subroutine get_Mullikens_all


subroutine get_Mulliken(Mulliken_model, masks_DOS, DOS_weights, Hij, fe, matter, MDatoms, mulliken_Ne, mulliken_q)
   integer, intent(in) :: Mulliken_model   ! which model to use
   logical, dimension(:,:,:), intent(in) :: masks_DOS   ! partial DOS made of each orbital type, if required to be constructed
   real(8), dimension(:,:,:), intent(in) :: DOS_weights     ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:,:), intent(in) :: Hij      ! real eigenvectors
   real(8), dimension(:), intent(in) :: fe    ! electron distribution
   type(Solid), intent(in) :: matter     ! material parameters
   type(Atom), dimension(:), intent(in) :: MDAtoms ! all MD atoms
   real(8), dimension(:), intent(out) :: mulliken_Ne   ! Mulliken electron populations
   real(8), dimension(:), intent(out), optional :: mulliken_q   ! Mulliken charges
   !-------------------------------
   real(8) :: temp
   integer :: j, Nsiz, N_at, N_types, i_at, i_types, Nat
   if (Mulliken_model >= 1) then ! get Mulliken populations and charges
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
         if (Nat <= 0) Nat = 1   ! in case there are no atoms of this kind
         ! Mulliken electron populations:
         mulliken_Ne(i_at) = mulliken_Ne(i_at) / dble(Nat)
      enddo
      ! Mulliken charges:
      if (present(mulliken_q)) mulliken_q(:) = matter%Atoms(:)%NVB - mulliken_Ne(:)
      !if (present(mulliken_q)) mulliken_q(:) = mulliken_Ne(:) - matter%Atoms(:)%NVB
   else  ! just atomic electrons
      ! Mulliken electron populations:
      mulliken_Ne(:) = matter%Atoms(:)%NVB
      ! Mulliken charges:
      if (present(mulliken_q)) mulliken_q(:) = 0.0d0
   endif ! (Mulliken_model >= 1)
end subroutine get_Mulliken



subroutine get_Mulliken_each_atom(Mulliken_model, Scell, matter, numpar)
   integer, intent(in) :: Mulliken_model   ! which model to use
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter     ! material parameters
   type (Numerics_param), intent(in) :: numpar    ! numerical parameters, including drude-function
   !-------------------------------
   real(8), dimension(size(Scell%Ha,1), size(Scell%Ha,2)) :: D
   real(8), dimension(:), allocatable :: mulliken_Ne
   integer :: N_at, i_at, i_orb, j, Nsiz, N_orb, k

   if (numpar%save_XYZ .and. numpar%save_XYZ_extra(2)) then ! get Mulliken populations and charges
      !print*, 'get_Mulliken_each_atom:', Mulliken_model
      !print*, allocated(Scell%MDAtoms)
      N_at = size(Scell%MDAtoms) ! total number of atoms
      Nsiz = size(Scell%Ha,1) ! total number of orbitals
      N_orb = Nsiz/N_at ! orbitals per atom

      if (allocated(Scell%Sij)) then
         !$omp PARALLEL private(j, k)
         !$omp do
         do j = 1, Nsiz ! for all energy levels
            do k = 1, Nsiz ! for all energy levels
               D(k,j) = Scell%Ha(k,j) * SUM(Scell%Ha(:,j) * Scell%Sij(k,:))   ! the density matrix without occupations
            enddo
         enddo
         !$omp end do
         !$omp end parallel
      else ! orthogonal
         !$omp PARALLEL private(j)
         !$omp do
         do j = 1, Nsiz ! for all energy levels
            D(:,j) = Scell%Ha(:,j) * Scell%Ha(:,j)
         enddo
         !$omp end do
         !$omp end parallel
      endif


      allocate(mulliken_Ne(N_at), source = 0.0d0)   ! absolute charges

      !$omp PARALLEL private(i_at, i_orb, j)
      !$omp do
      do i_at = 1, N_at ! all atoms
         do i_orb = 1, N_orb  ! all orbitals of each atom
            j = (i_at-1)*N_orb + i_orb ! current orbital among all
            mulliken_Ne(i_at) = mulliken_Ne(i_at) + SUM(Scell%fe(:) * D(j,:))
         enddo   ! i_orb

         ! Mulliken charge as a deviation from the normal electron population:
         Scell%MDAtoms(i_at)%q = matter%Atoms( Scell%MDAtoms(i_at)%KOA )%NVB - mulliken_Ne(i_at)
         !print*, i_at, Scell%MDAtoms(i_at)%q, matter%Atoms( Scell%MDAtoms(i_at)%KOA )%mulliken_q
      enddo ! i_at
      !$omp end do
      !$omp end parallel

      ! Testing:
      if (numpar%verbose) then
         do i_at = 1, size(matter%Atoms)
            N_at = COUNT(MASK = (Scell%MDatoms(:)%KOA == i_at))
            if (N_at <= 0) N_at = 1   ! in case there are no atoms of this kind
            write(*,'(a,f,f)') ' Mulliken average charge: '//trim(adjustl(matter%Atoms(i_at)%Name))//':', &
                     SUM(Scell%MDAtoms(:)%q, MASK = (Scell%MDatoms(:)%KOA == i_at)) / dble(N_at), &
                     matter%Atoms(i_at)%mulliken_q
         enddo
      endif

      deallocate(mulliken_Ne)
   else  ! just atomic electrons
      ! Mulliken charges:
      Scell%MDAtoms(:)%q = 0.0d0
   endif ! (Mulliken_model >= 1)
end subroutine get_Mulliken_each_atom


subroutine band_potential_energy_atom(Scell)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   !-------------------------------
   real(8), dimension(size(Scell%Ha,1), size(Scell%Ha,2)) :: D
   real(8), dimension(:), allocatable :: mulliken_Ne
   integer :: N_at, i_at, i_orb, j, Nsiz, N_orb, k

   N_at = size(Scell%MDAtoms) ! total number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   N_orb = Nsiz/N_at ! orbitals per atom

   if (allocated(Scell%Sij)) then
      !$omp PARALLEL private(j, k)
      !$omp do
      do j = 1, Nsiz ! for all energy levels
         do k = 1, Nsiz ! for all energy levels
            D(k,j) = Scell%Ha(k,j) * SUM(Scell%Ha(:,j) * Scell%Sij(k,:))   ! the density matrix without occupations
         enddo
      enddo
      !$omp end do
      !$omp end parallel
   else ! orthogonal
      !$omp PARALLEL private(j)
      !$omp do
      do j = 1, Nsiz ! for all energy levels
         D(:,j) = Scell%Ha(:,j) * Scell%Ha(:,j)
      enddo
      !$omp end do
      !$omp end parallel
   endif

   allocate(mulliken_Ne(N_at), source = 0.0d0)   ! absolute charges

   !$omp PARALLEL private(i_at, i_orb, j)
   !$omp do
   do i_at = 1, N_at ! all atoms
      do i_orb = 1, N_orb  ! all orbitals of each atom
         j = (i_at-1)*N_orb + i_orb ! current orbital among all
         ! Sum up potential energy for all orbitals belonging to this atom:
         mulliken_Ne(i_at) = mulliken_Ne(i_at) + SUM( Scell%Ei(:) * Scell%fe(:) * D(j,:) )
      enddo   ! i_orb
      ! Add Potential energy of this atom (band contribution) to the precalculated repulsive part:
      Scell%MDAtoms(i_at)%Epot = Scell%MDAtoms(i_at)%Epot + mulliken_Ne(i_at)
   enddo ! i_at
   !$omp end do
   !$omp end parallel

   deallocate(mulliken_Ne)
end subroutine band_potential_energy_atom




subroutine get_electronic_thermal_parameters(numpar, Scell, NSC, matter, Err)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(solid), intent(in) :: matter   ! materil parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !----------------------------
   ! 1) Electron heat capacity:
   call get_electronic_heat_capacity(Scell, NSC, Scell(NSC)%Ce, numpar%do_kappa, numpar%DOS_weights, Scell(NSC)%Ce_part) ! module "Electron_tools"

   !----------------------------
   ! 2) Electronic entropy:
   !print*, 'Before Se'
   call electronic_entropy(Scell(NSC)%fe, Scell(NSC)%Se) ! module "Electron_tools"
   ! and equivalent equilibrium entropy:
   call electronic_entropy(Scell(NSC)%fe_eq, Scell(NSC)%Se_eq) ! module "Electron_tools"
   ! If needed, get partial ones (band-resolved) too:
   if (numpar%do_partial_thermal) then
      ! Valence band:
      call electronic_entropy(Scell(NSC)%fe, Scell(NSC)%Se_VB, i_end=Scell(NSC)%N_Egap) ! module "Electron_tools"
      call electronic_entropy(Scell(NSC)%fe_eq_VB, Scell(NSC)%Se_eq_VB, i_end=Scell(NSC)%N_Egap) ! module "Electron_tools"
      ! Conduction band:
      call electronic_entropy(Scell(NSC)%fe, Scell(NSC)%Se_CB, i_start=Scell(NSC)%N_Egap+1) ! module "Electron_tools"
      call electronic_entropy(Scell(NSC)%fe_eq_CB, Scell(NSC)%Se_eq_CB, i_start=Scell(NSC)%N_Egap+1) ! module "Electron_tools"
   endif

   !----------------------------
   ! 3) Orbital-resolved electronic parameters:
   call get_orbital_resolved_data(Scell(NSC), matter, numpar%DOS_weights, numpar)  ! module "Electron_tools"

   !----------------------------
   ! Electron heat conductivity, if required (does not work well...).
   ! Exclude it from here; instead, try to use Onsager coefficients in Optical module.
   !call get_electron_heat_conductivity(Scell, NSC, matter, numpar, Err) ! below
   !print*, 'Done get_electronic_thermal_parameters'
end subroutine get_electronic_thermal_parameters




subroutine get_electron_heat_conductivity(Scell, NSC, matter, numpar, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(solid), intent(in) :: matter	! materil parameters
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------
   real(8), dimension(3,3) :: kappa_e, kappa_ei, L_EE, L_ET, L_TT, kappa_ee
   real(8), dimension(:,:,:), allocatable :: kappa_e_part, L_EE_part, L_ET_part, L_TT_part
   real(8), dimension(:,:), allocatable :: v_e
   real(8), dimension(:), allocatable :: dfdE
   real(8) :: coef, temp(3,3), tau, eps, Emu, Etemp(3,3), E2temp(3,3)
   integer :: Nsiz, i, j, k, N_at, N_types, i_at, i_types, i_G1

   if (numpar%do_kappa) then ! only if requested
      if (Scell(NSC)%Te < 35.0d0) then ! effectively zero Te
         Scell(NSC)%kappa_e = 0.0d0

      else ! non-zero Te
         eps = 1.0d-8   ! precision
         Nsiz = size(Scell(NSC)%Ei)
         if (.not.allocated(Scell(NSC)%kappa_e_part)) allocate(Scell(NSC)%kappa_e_part(size(Scell(NSC)%G_ei_partial,1)))
         N_at = size(numpar%DOS_weights,1)    ! number of kinds of atoms
         N_types = size(numpar%DOS_weights,2) ! number of atomic shells (basis set size)
         ! Band resolved kappa calculations:
         allocate(kappa_e_part(3,3,size(Scell(NSC)%kappa_e_part)), source=0.0d0)
         allocate(   L_EE_part(3,3,size(Scell(NSC)%kappa_e_part)), source=0.0d0)
         allocate(   L_ET_part(3,3,size(Scell(NSC)%kappa_e_part)), source=0.0d0)
         allocate(   L_TT_part(3,3,size(Scell(NSC)%kappa_e_part)), source=0.0d0)


!          call print_time_step('Start: get_electron_heat_conductivity:', 1.0, msec=.true.)   ! module "Little_subroutines"

         !-------------------------
         ! Electron-ion contribution, according to [https://doi.org/10.1002/aenm.202200657], Eqs.(1-5):
         ! 1) Get the electron band velosities:
         allocate(v_e(Nsiz,3))   ! v_x, v_y, v_z
         call get_electron_band_velosities(numpar, Scell, NSC, v_e, Err) ! below

         ! 2) Get the derivative of Fermi:
         allocate(dfdE(Nsiz))
         do i = 1, Nsiz
            dfdE(i) = Diff_Fermi_E(Scell(NSC)%TeeV, Scell(NSC)%mu, Scell(NSC)%Ei(i)) ! module "Electron_tools"
         enddo

         ! 3) Collect the terms:
         L_EE = 0.0d0   ! to start with
         L_ET = 0.0d0   ! to start with
         L_TT = 0.0d0   ! to start with
         do i = 1, Nsiz

            Emu = Scell(NSC)%Ei(i)-Scell(NSC)%mu

            !if ( (abs(Scell(NSC)%Ei(i)) > eps) .and. (abs(Scell(NSC)%I_ij(i)) > eps) .and. (abs(dfdE(i)) > eps) ) then ! [O 1]
            if ( (abs(Emu) > eps) .and. (abs(dfdE(i)) > eps) ) then  ! [O 2]
               ! Scattering time:
               !tau = (Scell(NSC)%Te-Scell(NSC)%Ta)/Scell(NSC)%Ei(i) * Scell(NSC)%Ce_i(i)/Scell(NSC)%I_ij(i) ! [O 1]
               if (abs(Scell(NSC)%I_ij(i)) > eps) then
                  tau = (Scell(NSC)%Te-Scell(NSC)%Ta)/Emu * Scell(NSC)%Ce_i(i)/Scell(NSC)%I_ij(i) ! [O 2]
               else
                  tau = 0.0d0
               endif
!                print*, 'tau', i, Scell(NSC)%Ei(i), Scell(NSC)%I_ij(i), tau
               tau = abs(tau)

               ! Velosities:
               forall (j=1:3, k=1:3)
                  !temp(j,k) = v_e(i,j)*v_e(i,k) * dfdE(i) * tau  ! [O 1]
                  temp(j,k) = v_e(i,j)*v_e(i,k) * Scell(NSC)%Ce_i(i) * tau  ! [O 2]
               end forall
               Etemp(:,:) = Emu*temp(:,:)
               E2temp(:,:) = Emu*Etemp(:,:)

               ! L terms:
               L_TT = L_TT + temp
               L_ET = L_ET + Etemp
               L_EE = L_EE + E2temp

               do i_at = 1, N_at ! all elements
                  do i_types = 1, N_types  ! all shells of each element
                     i_G1 = (i_at-1) * N_types + i_types
                     L_TT_part(:,:,i_G1) = L_TT_part(:,:,i_G1) + temp * numpar%DOS_weights(i_at, i_types, i)
                     L_ET_part(:,:,i_G1) = L_ET_part(:,:,i_G1) + Etemp * numpar%DOS_weights(i_at, i_types, i)
                     L_EE_part(:,:,i_G1) = L_EE_part(:,:,i_G1) + E2temp * numpar%DOS_weights(i_at, i_types, i)
                  enddo ! i_types
               enddo ! i_at

            endif ! (abs(Scell(NSC)%Ei(i)) > eps)
         enddo

         ! Kappa e-i:
         ! [O 1]:
!          where (abs(L_TT(:,:)) > eps)
!             kappa_ei = 1.0d0/(Scell(NSC)%V * Scell(NSC)%Te) * (L_EE - L_ET**2/L_TT) ! [eV*m^2/(A^3*K)]
!          elsewhere
!             kappa_ei = 1.0d0/(Scell(NSC)%V * Scell(NSC)%Te) * L_EE ! [eV*m^2/(A^3*K)]
!          endwhere
!          where (abs(L_TT_part(:,:,:)) > eps)
!             kappa_e_part = 1.0d0/(Scell(NSC)%V * Scell(NSC)%Te) * (L_EE_part - L_ET_part**2/L_TT_part) ! [eV*m^2/(A^3*K)]
!          elsewhere
!             kappa_e_part = 1.0d0/(Scell(NSC)%V * Scell(NSC)%Te) * L_EE_part ! [eV*m^2/(A^3*K)]
!          endwhere
         ! [O 2]:
         kappa_ei = L_TT/Scell(NSC)%V           ! [eV*m^2/(A^3*K*s)]
         kappa_e_part = L_TT_part/Scell(NSC)%V  ! [eV*m^2/(A^3*K*s)]

         ! Convert into SI units:
         coef = 1.0d30*g_e ! [eV/A^3] -> [J/m^3]
         kappa_ei = kappa_ei * coef ! [W/(m*K)]
         kappa_e_part = kappa_e_part * coef ! [W/(m*K)]

         !-------------------------
         ! Electron-electron scattering, using CDF cross-sections:
         kappa_ee = 1.0d30
         do i = 1, 3 ! diagonal elements have non-zero contributions
            kappa_ee(i,i) = 1.0d30
         enddo

         !-------------------------
         ! Contributions from electron-ion and electron-electron scattering,
         ! accordint to [10.1103/PhysRevD.74.043004], Eq.(2):
         where (kappa_ei(:,:) > eps)
            kappa_e = 1.0d0/(1.0d0/kappa_ei + 1.0d0/kappa_ee)  ! total electron conductivity
         elsewhere
            kappa_e = kappa_ei   ! total electron conductivity
         endwhere
         ! Save average electron conductivity:
         Scell(NSC)%kappa_e = 1.0d0/3.0d0 * (kappa_e(1,1) + kappa_e(2,2) + kappa_e(3,3))  ! trace of the matrix
         Scell(NSC)%kappa_e_part(:) = 1.0d0/3.0d0 * (kappa_e_part(1,1,:) + kappa_e_part(2,2,:) + kappa_e_part(3,3,:))

         ! Clean up:
         deallocate(v_e, dfdE)

!          call print_time_step('End :  get_electron_heat_conductivity:', 2.0, msec=.true.)   ! module "Little_subroutines"
      endif ! (Scell(NSC)%Te < 30.0d0)
   endif ! (numpar%do_kappa)
end subroutine get_electron_heat_conductivity


subroutine get_electron_band_velosities(numpar, Scell, NSC, v_e, Err)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), intent(inout) :: v_e
   type(Error_handling), intent(inout) :: Err	! error save
   !--------------------------------------
   complex, dimension(:,:), allocatable :: CHij, CSij	! eigenvectors of the hamiltonian
   real(8), dimension(:), allocatable :: Ei	! [eV] current energy levels at the selected k point
   real(8), dimension(:,:), allocatable :: E_c ! [eV]
   real(8) :: k_shift, k_shift2, coef, A, B, C, k_add
   integer :: Nsiz, i

!    call print_time_step('get_electron_band_velosities:', 1.0, msec=.true.)   ! module "Little_subroutines"

   k_shift = 0.05d0
   k_shift2 = 2.0d0 * k_shift
   k_add = k_shift + k_shift2
   Nsiz = size(Scell(NSC)%Ei)

   ! Get the energy levels:
   allocate(E_c(Nsiz,6))

   ! Get energy levels at given k-points (function below):
   call get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, k_shift, 0.0d0, 0.0d0, Err)  ! below
   E_c(:,1) = Ei  ! -kx
!    call print_time_step('get_electron_band_velosities:', 2.0, msec=.true.)   ! module "Little_subroutines"
!    print*, 'Ec 1:', E_c(1,1), Scell(NSC)%Ei(1)

   call get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, k_shift2, 0.0d0, 0.0d0, Err)  ! below
   E_c(:,2) = Ei  ! +kx
!    call print_time_step('get_electron_band_velosities:', 3.0, msec=.true.)   ! module "Little_subroutines"
!    print*, 'Ec 2:', E_c(1,2), Scell(NSC)%Ei(1)

   call get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, 0.0d0, k_shift, 0.0d0, Err)  ! below
   E_c(:,3) = Ei  ! -ky
   call get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, 0.0d0, k_shift2, 0.0d0, Err)  ! below
   E_c(:,4) = Ei  ! +ky
   call get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, 0.0d0, 0.0d0, k_shift, Err)  ! below
   E_c(:,5) = Ei  ! -kz
   call get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, 0.0d0, 0.0d0, k_shift2, Err)  ! below
   E_c(:,6) = Ei  ! +kz

   ! Now get the derivetives dEi/dk:
   do i = 1, Nsiz
      call fit_parabola_to_3points(0.0d0, Scell(NSC)%Ei(i), k_shift, E_c(i,1), k_shift2, E_c(i,2), A, B, C)  ! module "Algebra_tools"
      !v_e(i,1) = B   ! v_x
      v_e(i,1) = m_two_third*(A * k_add) + B   ! v_x
      call fit_parabola_to_3points(0.0d0, Scell(NSC)%Ei(i), k_shift, E_c(i,3), k_shift2, E_c(i,4), A, B, C)  ! module "Algebra_tools"
      !v_e(i,2) = B   ! v_y
      v_e(i,2) = m_two_third*(A * k_add) + B   ! v_y
      call fit_parabola_to_3points(0.0d0, Scell(NSC)%Ei(i), k_shift, E_c(i,5), k_shift2, E_c(i,6), A, B, C)  ! module "Algebra_tools"
      !v_e(i,3) = B   ! v_z
      v_e(i,3) = m_two_third*(A * k_add) + B   ! v_z
   enddo

   ! Convert into SI units:
   coef = (g_e*1.0d-10)/g_h ! [eV*A]/h -> [J*m]/h = [m/s]
   v_e = v_e * coef
   ! Check for consistency:
   where(abs(v_e) > g_cvel)  ! cannot be bigger than the speed of light:
      v_e = sign(g_cvel,v_e)
   endwhere
!    print*, 'Ve=', v_e
!    pause 'get_electron_band_velosities'

   ! Clean up:
   deallocate(Ei, E_c)
   call deallocate_array(CHij)   ! module "Little_sobroutine"
   call deallocate_array(CSij)   ! module "Little_sobroutine"
end subroutine get_electron_band_velosities




subroutine get_complex_Hamiltonian_new(numpar, Scell, NSC,  CHij, CSij, Ei, kx, ky, kz, Err)
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
         type is (TB_H_Pettifor)	! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err) ! below
         type is (TB_H_Molteni)	! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err)  ! below
         type is (TB_H_Fu)		! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err) ! below
         type is (TB_H_NRL)	! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, &
               Sij=Scell(NSC)%Sij, CSij_out=CSij) ! below
         type is (TB_H_DFTB) 	! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, &
               Sij=Scell(NSC)%Sij, CSij_out=CSij) ! below
         type is (TB_H_3TB) 	! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, &
               Sij=Scell(NSC)%Sij, CSij_out=CSij) ! below
         type is (TB_H_xTB) 	! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, &
               Sij=Scell(NSC)%Sij, CSij_out=CSij) ! below
         end select
      END ASSOCIATE
   endif
end subroutine get_complex_Hamiltonian_new


subroutine get_complex_Hamiltonian(numpar, Scell, NSC,  CHij, CSij, Ei, kx, ky, kz, Err, flag_old)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(inout), allocatable :: Ei	! [eV] current energy levels at the selected k point
   complex, dimension(:,:), intent(inout), allocatable :: CHij	! eigenvectors of the hamiltonian
   complex, dimension(:,:), intent(inout), allocatable :: CSij	! overlap matrix of the nonorthogonal hamiltonian
   real(8), intent(in) :: kx, ky, kz	! k point
   type(Error_handling), intent(inout) :: Err	! error save
   logical, intent(in) :: flag_old  ! to eventually make a call of the new subroutine
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




subroutine construct_complex_Hamiltonian(numpar, Scell, NSC, H_non, CHij, Ei, ksx, ksy, ksz, &
            Err, cPRRx, cPRRy, cPRRz, Sij, CSij_out, cTnn)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), intent(in) :: H_non	! Non-diagonalized real hamiltonian, must be provided
   complex, dimension(:,:), intent(out), allocatable :: CHij	! complex the hamiltonian, to be constructed
   real(8), dimension(:), intent(inout), allocatable :: Ei	! [eV] current energy levels at the selected k point
   real(8), intent(in) ::  ksx, ksy, ksz	! k point [relative coordinates]
   type(Error_handling), intent(inout), optional :: Err	! error save
   complex, dimension(:,:), intent(out), allocatable, optional :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   real(8), dimension(:,:), intent(in), optional :: Sij	   ! real overlap matrix of the nonorthogonal hamiltonian, must be provided for nonorthogonal case
   complex, dimension(:,:), intent(out), allocatable, optional :: CSij_out  ! complex overlap matrix of the nonorthogonal hamiltonian, to be created
   real(8), dimension(:,:,:,:), intent(out), allocatable, optional :: cTnn ! kinetic energy-related operators
   !-------------------------------------------------------------------------------
   !complex(8), dimension(:,:), allocatable :: CHij_temp, CHij_non, CSij_save, CHij_orth, CWF_orth, cPPRx_0, cPPRy_0, cPPRz_0, CSij
   !complex(8), dimension(:,:,:,:), allocatable :: cTnn_0, cTnn_c
   !complex(8) :: expfac, SH_1
   complex, dimension(:,:), allocatable :: CHij_temp, CHij_non, CSij_save, CHij_orth, CWF_orth, cPPRx_0, cPPRy_0, cPPRz_0, CSij
   complex, dimension(:,:,:,:), allocatable :: cTnn_0, cTnn_c
   complex :: expfac, SH_1
   real(8), dimension(:), allocatable :: Norm1
   real(8), dimension(3,3,3,3), target :: distances_to_image_cells
   integer :: Nsiz, j, nat, m, atom_2, i, j1, l, i1, k, norb, n, nn, cell_x, cell_y, cell_z
   real(8), dimension(3,3) :: k_supce   ! reciprocal supercell vectors
   real(8) :: temp, kx, ky, kz, zb(3), R, x, y, z
   real(8), target :: nol, dx, dy, dz
   real(8), pointer :: x1, y1, z1
   character(200) :: Error_descript

   nol = 0.0d0
   Error_descript = ''
   nat = size(Scell(NSC)%MDatoms) ! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
   norb =  Nsiz/nat ! orbitals per atom
   
   ! Allocate complex parameters:
   if (.not.allocated(CHij)) allocate(CHij(Nsiz,Nsiz))
   CHij = cmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij_temp)) allocate(CHij_temp(Nsiz,Nsiz))
   if (.not.allocated(CHij_non)) allocate(CHij_non(Nsiz,Nsiz))
   CHij_non = cmplx(0.0d0,0.0d0)	! to start with
   CHij_temp = cmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij_orth)) allocate(CHij_orth(Nsiz,Nsiz))  ! orthogonalized Hamiltonian
   CHij_orth = cmplx(0.0d0,0.0d0)	! to start with
   !select case (abs(numpar%optic_model))
   !case (4:5) ! Graf-Vogl or Kubo-Greenwood
   if (numpar%do_kappa .or. (abs(numpar%optic_model) == 4) .or. (abs(numpar%optic_model) == 5)) then ! if requested
      if (.not.allocated(CWF_orth)) allocate(CWF_orth(Nsiz,Nsiz), source = cmplx(0.0d0,0.0d0))  ! orthogonalized Hamiltonian
      if (.not.allocated(Norm1)) allocate(Norm1(Nsiz), source = 0.0d0)  ! normalization of WF
      if (.not.allocated(cPPRx_0)) allocate(cPPRx_0(Nsiz,Nsiz), source = cmplx(0.0d0,0.0d0))  ! temporary
      if (.not.allocated(cPPRy_0)) allocate(cPPRy_0(Nsiz,Nsiz), source = cmplx(0.0d0,0.0d0))  ! temporary
      if (.not.allocated(cPPRz_0)) allocate(cPPRz_0(Nsiz,Nsiz), source = cmplx(0.0d0,0.0d0))  ! temporary
      if (present(cTnn)) then
         if (.not.allocated(cTnn_0)) allocate(cTnn_0(Nsiz,Nsiz,3,3), source = cmplx(0.0d0,0.0d0))  ! temporary
         if (.not.allocated(cTnn_c)) allocate(cTnn_c(Nsiz,Nsiz,3,3), source = cmplx(0.0d0,0.0d0))  ! temporary
      endif
   end if


   if (.not.allocated(CSij_save)) allocate(CSij_save(Nsiz,Nsiz), source=cmplx(0.0d0,0.0d0))
   if (present(Sij)) then   ! it is nonorthogonal case:
      if (.not.allocated(CSij)) allocate(CSij(Nsiz,Nsiz))
      CSij = cmplx(0.0d0,0.0d0)	! to start with
      if (present(CSij_out)) then
         if (.not.allocated(CSij_out)) allocate(CSij_out(Nsiz,Nsiz), source=cmplx(0.0d0,0.0d0))
      endif
   endif
   
   ! Get reciprocal vectors:
   call Reciproc(Scell(NSC)%supce, k_supce) ! create reciprocal super-cell, module "Algebra_tools"
   call Convert_reciproc_rel_to_abs(ksx, ksy, ksz, k_supce, kx, ky, kz) ! get absolute k-values, molue "Atomic_tools"

   ! Get distances to image cells:
   do i = 1, 3 ! x
      zb(1) = dble(i-2) ! -1:1
      do j = 1, 3 ! y
         zb(2) = dble(j-2) ! -1:1
         do k = 1, 3 ! z
            zb(3) = dble(k-2) ! -1:1
            ! Get the distance to the replica of an atom #1 in the emage cell:
            call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDAtoms, zb, 1, 1, R, x, y, z)  ! module "Atomic_tools"
            ! Save it into an array:
            distances_to_image_cells(i,j,k,1) = x
            distances_to_image_cells(i,j,k,2) = y
            distances_to_image_cells(i,j,k,3) = z
         enddo ! k
      enddo ! j
   enddo ! i


   ! 1) Construct complex Hamiltonian and overlap:
!    !$omp parallel
!    !$omp do schedule(dynamic) private(j, m, atom_2, i, x1, y1, z1, R, cell_x, cell_y, cell_z, expfac, j1, l, i1, k)
   do j = 1,nat   ! all atoms
      m = Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one  
         if (atom_2 == 0) then
            i = j
            x1 => nol
            y1 => nol
            z1 => nol
         else
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
!             x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1) ! at this distance, X
!             y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2) ! at this distance, Y
!             z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3) ! at this distance, Z

            ! Find the used image cell indices:
            !call shortest_distance(Scell(NSC), j, i, R, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z) ! module "Atomic_tools"
            call shortest_distance(Scell(NSC), i, j, R, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z) ! module "Atomic_tools"
            ! convert from cell index to array index:
            cell_x = cell_x+2
            cell_y = cell_y+2
            cell_z = cell_z+2
            ! get the distance to the given cell:
            x1 => distances_to_image_cells(cell_x,cell_y,cell_z,1)
            y1 => distances_to_image_cells(cell_x,cell_y,cell_z,2)
            z1 => distances_to_image_cells(cell_x,cell_y,cell_z,3)
         endif ! (atom_2 .EQ. 0)
         
         if ((abs(kx) < 1.0d-14) .AND. (abs(ky) < 1.0d-14) .AND. (abs(kz) < 1.0d-14)) then
            expfac = cmplx(1.0d0,0.0d0)
         else
            expfac = exp(g_CI*cmplx(kx*x1 + ky*y1 + kz*z1,0.0d0))
         endif

         do j1 = 1,norb ! all orbitals
            l = (j-1)*norb+j1
            do i1 = 1,norb ! all orbitals
               k = (i-1)*norb+i1

               CHij_temp(k,l) = CMPLX(H_non(k,l),0.0d0)*expfac
               if ((isnan(real(CHij_temp(k,l)))) .OR. isnan(aimag(CHij_temp(k,l)))) then
                  Error_descript = 'CHij_temp ISNAN in construct_complex_Hamiltonian'
                  if (present(Err)) call Save_error_details(Err, 6, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  print*, i, j, k, l, CHij_temp(k,l)
               endif
               
               if (present(Sij)) then
                  CSij(k,l) = CMPLX(Sij(k,l),0.0d0)*expfac
                  if ((isnan(real(CSij(k,l)))) .OR. isnan(aimag(CSij(k,l)))) then
                     Error_descript = 'CHij_temp ISNAN in construct_complex_Hamiltonian'
                     if (present(Err)) call Save_error_details(Err, 6, Error_descript)
                     print*, trim(adjustl(Error_descript))
                     print*, i, j, k, l, CSij(k,l)
                  endif
               endif
            enddo ! i1
         enddo ! j1
      enddo ! atom_2
   enddo ! j
!    !$omp end do
!    !$omp end parallel


   ! Temporarily save nonorthogonal Hamiltonian and overlap matrix:
   CHij_non = CHij_temp
   if (present(Sij)) then  ! nonorthogonal
      CSij_save = CSij
      !select case (abs(numpar%optic_model))
      !case (2)   ! Trani
      if (abs(numpar%optic_model) == 2) then
         call diagonalize_complex_Hamiltonian(CHij_temp, Ei, CSij)    ! below
      !case (4:5)   ! Graf-Vogl or KG
      elseif (numpar%do_kappa .or. (abs(numpar%optic_model) == 4) .or. (abs(numpar%optic_model) == 5)) then ! if requested
         !call diagonalize_complex8_Hamiltonian(CHij_temp, Ei, CSij, CHij_orth, CWF_orth)    ! below
         call diagonalize_complex_Hamiltonian(CHij_temp, Ei, CSij, CHij_orth, CWF_orth)    ! below
      else  ! gamma point, no need to diagonalize
         allocate(Ei(size(Scell(NSC)%Ei)), source=Scell(NSC)%Ei)
         !CHij_temp = Scell(NSC)%Hij
      end if
   else  ! orthogonal
      CSij_save = 0.0d0 ! exclude nonorthogonal contribution
      !select case (abs(numpar%optic_model))
      !case (2)   ! Trani
      if (abs(numpar%optic_model) == 2) then
         call diagonalize_complex_Hamiltonian(CHij_temp, Ei)    ! below
      !case (4:5)  ! Graf-Vogl or KG
      elseif (numpar%do_kappa .or. (abs(numpar%optic_model) == 4) .or. (abs(numpar%optic_model) == 5)) then ! if requested
         call diagonalize_complex_Hamiltonian(CHij_temp, Ei, CHij_orth=CHij_orth, CWF_orth=CWF_orth)    ! below
         !CWF_orth = TRANSPOSE(CWF_orth)   ! test
      else  ! gamma point, no need to diagonalize
         allocate(Ei(size(Scell(NSC)%Ei)), source=Scell(NSC)%Ei)
         !CHij_temp = Scell(NSC)%Hij
      end if
   endif
   CHij = CHij_temp ! save for output
   if (present(CSij_out)) then ! diagonalized Hamiltonian (eigenfunctions)
      CSij_out = CHij ! output
   endif

   !---------------------------------------
   ! Effective momentum operators:
   !select case (abs(numpar%optic_model))
   !case (2)  ! create matrix element for Trani method
   if (abs(numpar%optic_model) == 2) then
      if (.not.allocated(cPRRx)) allocate(cPRRx(Nsiz,Nsiz))
      if (.not.allocated(cPRRy)) allocate(cPRRy(Nsiz,Nsiz))
      if (.not.allocated(cPRRz)) allocate(cPRRz(Nsiz,Nsiz))
      cPRRx = cmplx(0.0d0,0.0d0)	! to start with
      cPRRy = cmplx(0.0d0,0.0d0)	! to start with
      cPRRz = cmplx(0.0d0,0.0d0)	! to start with
      ! Optical matrix elements for non-orthogonal TB are taken from:
      ! arXiv:1805.08918v1 -- https://128.84.21.199/abs/1805.08918
!       !$omp parallel
!       !$omp do private(j, m, atom_2, i, x1, y1, z1, j1, l, i1, k, SH_1)
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
                           SH_1 = CMPLX(0.27d0,0.0d0) * (CHij_non(k,l) - CMPLX(Ei(k),0.0d0)*CSij_save(k,l))
                           !SH_1 = DCMPLX(1.50d0,0.0d0) * (CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l))

                           cPRRx(k,l) = SH_1
                           cPRRy(k,l) = SH_1
                           cPRRz(k,l) = SH_1
                        endif

                     else  ! different atoms at distance {x,y,z}:
                        SH_1 = CHij_non(k,l)
                        if (present(Sij)) then ! nonorthogonal Hamiltonian:
                           ! [1] Nonorthogonal expression:
                           SH_1 = SH_1 - CMPLX(Ei(k),0.0d0)*CSij_save(k,l)  ! Correct
                        endif

                        cPRRx(k,l) = CMPLX(x1,0.0d0)*SH_1
                        cPRRy(k,l) = CMPLX(y1,0.0d0)*SH_1
                        cPRRz(k,l) = CMPLX(z1,0.0d0)*SH_1
                     endif
                     if (dble(cPRRx(k,l)) .GT. 1d10) write(*,'(i5,i5,es,es, es,es, es,es, es, es)') i, j, cPRRx(k,l),  CHij_non(k,l), CSij_save(k,l),  Ei(k), x1
                     if (dble(cPRRy(k,l)) .GT. 1d10) print*, i, j, cPRRy(k,l)
                     if (dble(cPRRz(k,l)) .GT. 1d10) print*, i, j, cPRRz(k,l)
                  enddo ! i1
               enddo ! j1
            endif ! (i > 0)
         enddo ! atom_2
      enddo ! j
!       !$omp end do
!       !$omp end parallel
      
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
         temp = 0.5d0   ! factor of 2 missing somewhere...
      type is (TB_H_DFTB)
         temp = 0.5d0   ! factor of 2 missing somewhere...
      type is (TB_H_3TB)
         temp = 1.0d0
      type is (TB_H_xTB)
         temp = 1.0d0
      end select
      END ASSOCIATE
      cPRRx = cPRRx * temp
      cPRRy = cPRRy * temp
      cPRRz = cPRRz * temp

   !---------------------------------------
   !case (4:5) ! Graf-Vogl or Kubo-Greenwood
   elseif (numpar%do_kappa .or. (abs(numpar%optic_model) == 4) .or. (abs(numpar%optic_model) == 5)) then ! if requested
      if (.not.allocated(cPRRx)) allocate(cPRRx(Nsiz,Nsiz))
      if (.not.allocated(cPRRy)) allocate(cPRRy(Nsiz,Nsiz))
      if (.not.allocated(cPRRz)) allocate(cPRRz(Nsiz,Nsiz))
      if (present(cTnn)) then
         if (.not.allocated(cTnn)) allocate(cTnn(Nsiz,Nsiz,3,3), source = 0.0d0)
      endif
      cPRRx = cmplx(0.0d0,0.0d0)	! to start with
      cPRRy = cmplx(0.0d0,0.0d0)	! to start with
      cPRRz = cmplx(0.0d0,0.0d0)	! to start with

!       !$omp parallel
!       !$omp do private(j, m, atom_2, i, x1, y1, z1, dx, dy, dz, j1, l, i1, k, SH_1, n, nn)
      do j = 1,nat   ! all atoms
         !m = Scell(NSC)%Near_neighbor_size(j)
         !do atom_2 = 0,m ! do only for atoms close to that one
            !if (atom_2 == 0) then
            !   i = j
            !else
            !   i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            !   x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1) ! at this distance, X
            !   y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2) ! at this distance, Y
            !   z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3) ! at this distance, Z
            !endif ! (atom_2 .EQ. 0)
         do i = 1,nat   ! all pairs of atoms
            call shortest_distance(Scell(NSC), j, i, R, x1=dx, y1=dy, z1=dz) ! module "Atomic_tools"
            x1 => dx
            y1 => dy
            z1 => dz

            if (i > 0) then
               do j1 = 1,norb	! all orbitals for sp3d5
                  l = (j-1)*norb+j1
                  do i1 = 1,norb	! all orbitals for sp3d5
                     k = (i-1)*norb+i1
                     if (i == j) then ! contribution of the same atom, according to Trani:
                        if (j1 == i1) then
                           ! Skip the same orbital, no overlap
                        else
                           ! i*(R-R') terms:
                           if (numpar%optic_model == 4) then ! orthogonal expression:
                              ! Orthogonalized Hamiltonian:
                              SH_1 = g_CI * CMPLX(0.27d0, 0.0d0) * CHij_orth(k,l)   ! testing
                           else  ! nonorthogonal
                              ! [1] Nonorthogonal expression:
                              SH_1 = g_CI * CMPLX(0.27d0, 0.0d0) * (CHij_non(k,l) - CMPLX(Ei(k),0.0d0)*CSij_save(k,l))
                           endif

                           cPRRx(k,l) = SH_1
                           cPRRy(k,l) = SH_1
                           cPRRz(k,l) = SH_1

                           if (present(cTnn)) then
                              SH_1 = -0.27d0**2
                              cTnn_c(k,l,:,:) = SH_1 * CHij_orth(k,l)
                           endif
                        endif
                     else  ! different atoms at distance {x,y,z}:
                        ! i*(R-R') terms:
                        if (numpar%optic_model == 4) then ! orthogonalized expression:
                           SH_1 = CHij_orth(k,l)
                        else  ! nonorthogonal
                           SH_1 = CHij_non(k,l)
                           if (present(Sij)) then ! nonorthogonal Hamiltonian:
                              ! [1] Nonorthogonal expression:
                              SH_1 = SH_1 - CMPLX(Ei(k),0.0d0)*CSij_save(k,l)  ! Correct
                           endif
                        endif

                        cPRRx(k,l) = g_CI * CMPLX(x1, 0.0d0)*SH_1
                        cPRRy(k,l) = g_CI * CMPLX(y1, 0.0d0)*SH_1
                        cPRRz(k,l) = g_CI * CMPLX(z1, 0.0d0)*SH_1

                        if (present(cTnn)) then
                           cTnn_c(k,l,1,1) = -(x1*x1) * CHij_orth(k,l)
                           cTnn_c(k,l,2,2) = -(y1*y1) * CHij_orth(k,l)
                           cTnn_c(k,l,3,3) = -(z1*z1) * CHij_orth(k,l)
                           SH_1 = -(x1*y1) * CHij_orth(k,l)
                           cTnn_c(k,l,1,2) = SH_1
                           cTnn_c(k,l,2,1) = SH_1 !cTnn(k,l,1,2) !-(y1*x1) * CHij_orth(k,l)
                           SH_1 = -(x1*z1) * CHij_orth(k,l)
                           cTnn_c(k,l,1,3) = SH_1 !-(x1*z1) * CHij_orth(k,l)
                           cTnn_c(k,l,3,1) = SH_1 !cTnn(k,l,1,3) !-(z1*x1) * CHij_orth(k,l)
                           SH_1 = -(y1*z1) * CHij_orth(k,l)
                           cTnn_c(k,l,2,3) = SH_1 !-(y1*z1) * CHij_orth(k,l)
                           cTnn_c(k,l,3,2) = SH_1 !cTnn(k,l,2,3) !-(z1*y1) * CHij_orth(k,l)
                        endif
                     endif

                     if (dble(cPRRx(k,l)) .GT. 1d10) write(*,'(i5,i5,es,es, es,es, es,es, es, es)') i, j, cPRRx(k,l),  CHij_non(k,l), CSij_save(k,l),  Ei(k), x1
                     if (dble(cPRRy(k,l)) .GT. 1d10) print*, i, j, cPRRy(k,l)
                     if (dble(cPRRz(k,l)) .GT. 1d10) print*, i, j, cPRRz(k,l)
                  enddo ! i1
               enddo ! j1
            endif ! (i > 0)
         enddo ! atom_2
      enddo ! j
!       !$omp end do
!       !$OMP BARRIER

      if (numpar%optic_model /= 4) then ! use nonorthogonal wave-functions:
         CWF_orth = CHij_temp ! copy them into this array used below
      endif

!       !$omp do
      do i = 1, Nsiz ! ensure WF normalization to 1
         Norm1(i) = SQRT( SUM( conjg(CWF_orth(:,i)) * CWF_orth(:,i) ) )
         !print*, i, Norm1(i)  ! not 1 for non-orthogonal Hamiltonians
      enddo
!       !$omp end do
!       !$OMP BARRIER
!       !$omp do
      do i = 1, Nsiz
         do nn = 1, Nsiz
            cPPRx_0(i,nn) = SUM(cPRRx(i,:)*CWF_orth(:,nn)) / Norm1(nn)
            cPPRy_0(i,nn) = SUM(cPRRy(i,:)*CWF_orth(:,nn)) / Norm1(nn)
            cPPRz_0(i,nn) = SUM(cPRRz(i,:)*CWF_orth(:,nn)) / Norm1(nn)
         enddo ! j
      enddo ! i
!       !$omp end do
!       !$OMP BARRIER


      if (present(cTnn)) then
!          !$omp do
         do i = 1, Nsiz ! upper triangle
            do nn = 1, Nsiz
               cTnn_0(i,nn,1,1) = SUM(cTnn_c(i,:,1,1)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,1,2) = SUM(cTnn_c(i,:,1,2)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,1,3) = SUM(cTnn_c(i,:,1,3)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,2,2) = SUM(cTnn_c(i,:,2,2)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,2,3) = SUM(cTnn_c(i,:,2,3)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,3,3) = SUM(cTnn_c(i,:,3,3)*CWF_orth(:,nn)) / Norm1(nn)
            enddo ! j
         enddo ! i
!          !$omp end do
!          !$OMP BARRIER
!          !$omp do
         do i = 1, Nsiz ! lower triangle
            do nn = 1, Nsiz
               cTnn_0(i,nn,2,1) = cTnn_0(i,nn,1,2) !SUM(cTnn(i,:,2,1)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,3,1) = cTnn_0(i,nn,1,3) !SUM(cTnn(i,:,3,1)*CWF_orth(:,nn)) / Norm1(nn)
               cTnn_0(i,nn,3,2) = cTnn_0(i,nn,2,3) !SUM(cTnn(i,:,3,2)*CWF_orth(:,nn)) / Norm1(nn)
            enddo ! j
         enddo ! i
!          !$omp end do
!          !$OMP BARRIER
      endif

!       !$omp do
      do n = 1, Nsiz
         do nn = 1, Nsiz
            cPRRx(n,nn) = SUM(conjg(CWF_orth(:,n))*cPPRx_0(:,nn)) / Norm1(n)   ! correct
            cPRRy(n,nn) = SUM(conjg(CWF_orth(:,n))*cPPRy_0(:,nn)) / Norm1(n)
            cPRRz(n,nn) = SUM(conjg(CWF_orth(:,n))*cPPRz_0(:,nn)) / Norm1(n)
!             cPRRx(n,nn) = SUM(conjg(CWF_orth(n,:))*cPPRx_0(:,nn)) / Norm1(n)  ! wrong
!             cPRRy(n,nn) = SUM(conjg(CWF_orth(n,:))*cPPRy_0(:,nn)) / Norm1(n)
!             cPRRz(n,nn) = SUM(conjg(CWF_orth(n,:))*cPPRz_0(:,nn)) / Norm1(n)
         enddo ! nn
      enddo ! n
!       !$omp end do
!       !$OMP BARRIER


      if (present(cTnn)) then
!          !$omp do
         do n = 1, Nsiz ! upper triangle
            do nn = 1, Nsiz
               cTnn_c(n,nn,1,1) = SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,1,1)) / Norm1(n)
               cTnn_c(n,nn,1,2) = SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,1,2)) / Norm1(n)
               cTnn_c(n,nn,1,3) = SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,1,3)) / Norm1(n)
               cTnn_c(n,nn,2,2) = SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,2,2)) / Norm1(n)
               cTnn_c(n,nn,2,3) = SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,2,3)) / Norm1(n)
               cTnn_c(n,nn,3,3) = SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,3,3)) / Norm1(n)
            enddo ! nn
         enddo ! n
!          !$omp end do
!          !$OMP BARRIER
!          !$omp do
         do n = 1, Nsiz ! lower triangle
            do nn = 1, Nsiz
               cTnn_c(n,nn,2,1) = cTnn_c(n,nn,1,2)  !SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,2,1)) / Norm1(n)
               cTnn_c(n,nn,3,1) = cTnn_c(n,nn,1,3)  !SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,3,1)) / Norm1(n)
               cTnn_c(n,nn,3,2) = cTnn_c(n,nn,2,3)  !SUM(conjg(CWF_orth(:,n))*cTnn_0(:,nn,3,2)) / Norm1(n)
            enddo ! j
         enddo ! i
!          !$omp end do
      endif
!       !$omp end parallel

      ! Convert in SI units of momentum:
      temp = g_me/g_h * g_e*1.0d-10 ! [kg*m/s]
      if ( present(Sij) .and. (numpar%optic_model /= 4) ) then ! for non-orthogonal calculations only
         temp = temp*0.5d0 ! correcting factor
      endif
      cPRRx = cPRRx * temp
      cPRRy = cPRRy * temp
      cPRRz = cPRRz * temp

      ! Term with k, according to
      ! [B. Holst, M. French, R. Redmer, Phys. Rev. B 83, 235120 (2011)], Eq.(18):
!       do n = 1, Nsiz
!          cPRRx(n,n) = cPRRx(n,n) + DCMPLX( sqrt(kx**2 + ky**2 + kz**2)*1.0d10*g_h, 0.0d0 )
!          cPRRy(n,n) = cPRRy(n,n) + DCMPLX( sqrt(kx**2 + ky**2 + kz**2)*1.0d10*g_h, 0.0d0 )
!          cPRRz(n,n) = cPRRz(n,n) + DCMPLX( sqrt(kx**2 + ky**2 + kz**2)*1.0d10*g_h, 0.0d0 )
!       enddo

      ! And kinetic energy:
      if (present(cTnn)) then
         cTnn = dble(cTnn_c) * temp/g_h*1.0d-10  ! [dimensionless]
      endif
   end if

   deallocate(CHij_temp, CHij_non)
   if (allocated(CSij_save)) deallocate(CSij_save)
   if (allocated(CHij_orth)) deallocate(CHij_orth)
   if (allocated(CWF_orth)) deallocate(CWF_orth)
   if (allocated(CSij)) deallocate(CSij)
   if (allocated(cPPRx_0)) deallocate(cPPRx_0)
   if (allocated(cPPRy_0)) deallocate(cPPRy_0)
   if (allocated(cPPRz_0)) deallocate(cPPRz_0)
   if (allocated(Norm1)) deallocate(Norm1)
   if (allocated(cTnn_0)) deallocate(cTnn_0)
   if (allocated(cTnn_c)) deallocate(cTnn_c)
   nullify(x1, y1, z1)
end subroutine construct_complex_Hamiltonian





subroutine diagonalize_complex_Hamiltonian(CHij, Ei, CSij, CHij_orth, CWF_orth)
   complex, dimension(:,:), intent(inout) :: CHij	! complex hermitian Hamiltonian
   real(8), dimension(:), intent(out), allocatable :: Ei ! eigenvalues [eV]
   !type(Error_handling), intent(inout) :: Err	! error save
   complex, dimension(:,:), intent(inout), optional ::  CSij    ! overlap matrix, for nonorthogonal Hamiltonain
   complex, dimension(:,:), intent(inout), optional ::  CHij_orth ! orthogonalized Hamiltonian, if requested
   complex, dimension(:,:), intent(inout), optional ::  CWF_orth  ! Wave functions of orthogonalized Hamiltonian, if requested
   !----------------------------
   !complex, dimension(size(CHij,1), size(CHij,2)) :: CHij_temp
   complex, dimension(:,:), allocatable :: CHij_temp
   integer :: Nsiz, j
   character(200) :: Error_descript


   Error_descript = ''  ! to start with, no error
   Nsiz = size(CHij,1)
   if (.not.allocated(Ei)) allocate(Ei(Nsiz))

   ORTH: if (.not.present(CSij)) then ! orthogonal:
      if (present(CHij_orth)) then  ! Save orthogonal Hamiltonian (for optical coefficients below)
         CHij_orth = CHij
      endif

      ! Direct diagonalization:
      call sym_diagonalize(CHij, Ei, Error_descript) ! modeule "Algebra_tools"

      if (present(CWF_orth)) then   ! Save WF of orthogonal Hamiltonian
         CWF_orth = CHij
      endif

   else ORTH ! nonorthogonal
      ! Solve linear eigenproblem:
      ! 1) Orthogonalize the Hamiltonian using Loewdin procidure:
      ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:

      allocate(CHij_temp(Nsiz,Nsiz))

      CHij_temp = CHij
      call Loewdin_Orthogonalization_c(Nsiz, CSij, CHij_temp)	! module "TB_NRL"

      if (present(CHij_orth)) then  ! Save orthogonalized Hamiltonian (for optical coefficients below)
         CHij_orth = CHij_temp
      endif

      ! 2) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
      call sym_diagonalize(CHij_temp, Ei, Error_descript)   ! module "Algebra_tools"
      if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
         Error_descript = 'diagonalize_complex_Hamiltonian: '//trim(adjustl(Error_descript))
         !call Save_error_details(Err, 6, Error_descript)    ! module "Objects"
         print*, trim(adjustl(Error_descript))
      endif
      if (present(CWF_orth)) then   ! Save WF of orthogonalized Hamiltonian
         CWF_orth = CHij_temp
      endif

      ! 3) Convert the eigenvectors back into the non-orthogonal basis:
      call mkl_matrix_mult('N', 'N', CSij, CHij_temp, CHij)	! module "Algebra_tools"
      
!       ! 4) If we need to renormalize the wave functions (in case they are not normalized to 1 after this procidure):
!       do j = 1, Nsiz
!         CHij(:,j) = CHij(:,j) / SQRT(SUM( dconjg(CHij(:,j)) * CHij(:,j) ))
!       enddo

!        do j = 1, size(Ei)
!           write(*,'(i5,es,es,es,es,es)') j, Ei(j), CHij_temp(j,1), Scell(NSC)%Ei(j)
!        enddo
!        PAUSE 'Ei pause'
      deallocate(CHij_temp)   ! clean up
   endif ORTH
end subroutine diagonalize_complex_Hamiltonian


subroutine diagonalize_complex8_Hamiltonian(CHij, Ei, CSij, CHij_orth, CWF_orth)
   complex, dimension(:,:), intent(inout) :: CHij	! complex hermitian Hamiltonian
   real(8), dimension(:), intent(out), allocatable :: Ei ! eigenvalues [eV]
   complex(8), dimension(:,:), intent(inout), optional ::  CSij    ! overlap matrix, for nonorthogonal Hamiltonain
   complex(8), dimension(:,:), intent(inout), optional ::  CHij_orth ! orthogonalized Hamiltonian, if requested
   complex(8), dimension(:,:), intent(inout), optional ::  CWF_orth  ! Wave functions of orthogonalized Hamiltonian, if requested
   !----------------------------
   !complex(8), dimension(size(CHij,1), size(CHij,1)) :: CHij_use ! This way of defining size breaks down!
   complex(8), dimension(:,:), allocatable :: CHij_use, CHij_temp
   integer :: Nsiz, j
   character(200) :: Error_descript

#ifdef OMP_inside
   print*, OMP_GET_THREAD_NUM(), 'diagonalize_complex8_Hamiltonian start'
#else
   print*, 'diagonalize_complex8_Hamiltonian start'
#endif

   Error_descript = ''  ! to start with, no error
   Nsiz = size(CHij,1)

   print*, 'Nsiz=', Nsiz

   if (.not.allocated(Ei)) allocate(Ei(Nsiz))
   allocate(CHij_use(Nsiz,Nsiz))

   ORTH: if (.not.present(CSij)) then ! orthogonal:

#ifdef OMP_inside
      print*, OMP_GET_THREAD_NUM(), 'diagonalize_complex8_Hamiltonian orthogonal'
#else
      print*, 'diagonalize_complex8_Hamiltonian orthogonal'
#endif

      ! Direct diagonalization:
      call sym_diagonalize(CHij, Ei, Error_descript) ! modeule "Algebra_tools"
   else ORTH ! nonorthogonal
      ! Solve linear eigenproblem:
      ! 1) Orthogonalize the Hamiltonian using Loewdin procidure:
      ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:

#ifdef OMP_inside
      print*, OMP_GET_THREAD_NUM(), 'diagonalize_complex8_Hamiltonian nonorthogonal'
#else
      print*, 'diagonalize_complex8_Hamiltonian nonorthogonal'
#endif

      CHij_use = CHij
      call Loewdin_Orthogonalization_c8(Nsiz, CSij, CHij_use) ! module "TB_NRL"

#ifdef OMP_inside
      print*, OMP_GET_THREAD_NUM(), 'Loewdin_Orthogonalization_c8 done'
#else
      print*, 'Loewdin_Orthogonalization_c8 done'
#endif

      if (present(CHij_orth)) then  ! Save orthogonalized Hamiltonian (for optical coefficients below)
         CHij_orth = CHij_use
      endif

      ! 2) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
      !call sym_diagonalize(CHij_use, Ei, Error_descript)   ! module "Algebra_tools"
      call c8_diagonalize(CHij_use, Ei, Error_descript)   ! module "Algebra_tools"
      if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
         Error_descript = 'diagonalize_complex8_Hamiltonian: '//trim(adjustl(Error_descript))
         !call Save_error_details(Err, 6, Error_descript)    ! module "Objects"
         print*, trim(adjustl(Error_descript))
      endif
      if (present(CWF_orth)) then   ! Save WF of orthogonalized Hamiltonian
         CWF_orth = CHij_use
      endif


#ifdef OMP_inside
      print*, OMP_GET_THREAD_NUM(), 'sym_diagonalize done'
#else
      print*, 'sym_diagonalize done'
#endif

      ! 3) Convert the eigenvectors back into the non-orthogonal basis:
      !call mkl_matrix_mult('N', 'N', CSij, CHij_use, CHij)	! module "Algebra_tools"
      allocate( CHij_temp(Nsiz,Nsiz), source = dcmplx( real(CHij), aimag(CHij)) )
      call mkl_matrix_mult_c8('N', 'N', CSij, CHij_use, CHij_temp)  ! module "Algebra_tools"
      CHij = cmplx( real(CHij_temp), aimag(CHij_temp) )  ! convert to cmplx (important for unix-based compileres)

#ifdef OMP_inside
      print*, OMP_GET_THREAD_NUM(), 'mkl_matrix_mult done'
#else
      print*, 'mkl_matrix_mult done'
#endif
   endif ORTH
end subroutine diagonalize_complex8_Hamiltonian



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


subroutine Exponential_wall_forces(TB_Expwall, Scell, NSC, matter, numpar) ! get Exponential wall forces
   class(TB_Exp_wall), allocatable, dimension(:,:), intent(in) :: TB_Expwall	! exponential wall
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of super-cell
   type(solid), intent(in) :: matter   ! materil parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !----------------------------------
    if (allocated(TB_Expwall)) then ! if we have short-range potential defined
      select type (TB_Expwall)
      type is (TB_Exp_wall_simple)
         ! Forces for all atoms:
         call d_Exp_wall_pot_s(Scell, NSC, TB_Expwall, numpar) ! module "Exponential_wall"
         ! Forces for the super-cell:
         call d_Exp_wall_Pressure_s(Scell, NSC, TB_Expwall, numpar)  ! module "Exponential_wall"
      type is (TB_Short_Rep)
         ! Forces for all atoms:
         call d_Short_range_pot_s(Scell, NSC, matter, TB_Expwall, numpar) ! module "Exponential_wall"
         ! Forces for the super-cell:
         call d_Short_range_Pressure_s(Scell, NSC, TB_Expwall, matter, numpar) ! module "Exponential_wall"
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
   
   if (allocated(TB_Coul)) then ! if we have Coupomb potential defined
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
      ! Get multipliers used many times into temporary arrays:
      call Construct_B(TB_Waals, Scell, NSC, numpar, Scell(NSC)%MDatoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, &
                        Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Van_der_Waals"
         
      ! Forces for all atoms:
      call d_Forces_s(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Nx, Ny, Nz) ! below
         
      ! Forces for the super-cell:
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
   allocate(Erx_s(3,n), source = 0.0d0) ! x,y,z-forces for each atoms

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
            XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
               YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
                  ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
                     zb = (/x_cell, y_cell, z_cell/) ! vector of image of the super-cell
                     origin_cell = ALL(zb==0) ! if it is the origin cell

                     cell_if:if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                        ! contribution from this image of the cell:
                        coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"

                        x1(1) = Xij(coun_cell,i1,j1)
                        x1(2) = Yij(coun_cell,i1,j1)
                        x1(3) = Zij(coun_cell,i1,j1)
                        a_r = A_rij(coun_cell,i1,j1) ! shortest distance

                        b = Bij(coun_cell,i1,j1) ! dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),a_r)

                        ddlta = dble(dik - djk)/a_r
                        b_delta = b*ddlta
                        dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                        dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                        dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
                     endif cell_if
                     !endif cos_if
                  enddo ZC2
               enddo YC2
            enddo XC2
         enddo ! j1
         !print*, allocated(Erx_s), size(Erx_s,2)
         Erx_s(1,ian) = Erx_s(1,ian) + dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      ! Add van der Waals force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)*0.5d0
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
      XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
         YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
            ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
               zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
               origin_cell = ALL(zb==0) ! if it is the origin cell
               do i = 1, n ! Forces from all atoms
                  Rep_Pr = 0.0d0 ! to start
                  dpsy = 0.0d0
                  do j = 1, n ! do for all pairs of atoms

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
                        dpsy = Bij(coun_cell,i,j)    ! dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),r)

                        do k = 1,3 ! supce indices: a,b,c
                           do l = 1,3  ! supce indices: x,y,z
                              Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                           enddo ! l
                        enddo ! k
                     endif cell_if
                  enddo ! j

                  do k = 1,3 ! supce indices
                     do l = 1,3  ! supce indices
                        !Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + Rep_Pr(l,k) !*0.5d0
                        PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)*0.5d0
                     enddo ! l
                  enddo ! k
               enddo ! i
            enddo ZC2
         enddo YC2
      enddo XC2
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add vdW part to existing TB part
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
   real(8) :: Erepuls, Pot_phi, Na_inv, E_vdW, E_coul, E_expwall
   integer :: NSC
   ! Calculations of the contribution of electrons into potential energy of atoms,
   ! first term in Eq.(2.44) from H.Jeschke PhD Thesis, Page 41
   !call set_total_el_energy(Ei,fe,Eelectr) ! module "Electron_tools"
   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
   
      ! Get the energy associated with the electrons populating band structure:
      if (numpar%scc) then ! energy for SCC is defined by for H_0:
!          call get_low_e_energy(Scell, matter, numpar) ! module "Electron_tools"
         call get_low_e_energy(Scell, matter) ! module "Electron_tools"
      else  ! no SCC, the usual expression then
         call get_low_e_energy(Scell, matter) ! module "Electron_tools"
      endif

      ! Calculations of the contribution of repulsive energies into potential energy of atoms,
      ! second term in Eq.(2.44) from H.Jeschke PhD Thesis, Page 41
      do NSC = 1, size(Scell) ! for all supercells
         ! Real number of atoms inversed:
         Na_inv = 1.0d0 / dble(Scell(NSC)%Na)

         ! Restart calculations of the potential energy of each atom:
         Scell(NSC)%MDAtoms(:)%Epot = 0.0d0
         
         !Erepuls = Erep_s(Scell(NSC)%TB_Repuls(:,:), Scell, NSC, numpar) ! below
         call Erep_s(Scell(NSC)%TB_Repuls(:,:), Scell, NSC, numpar, Erepuls) ! below
         
         Pot_phi = Scell(NSC)%nrg%El_low + Erepuls ! [eV] potential energy, Eq.(2.44) in H.Jeschke Thesis, Page 41
         
         Scell(NSC)%nrg%At_pot = Pot_phi * Na_inv ! [eV/atom]
         
         Scell(NSC)%nrg%E_rep = Erepuls ! [eV]
         
         ! van der Waals (vdW) potential energy:
         call vdW_s(Scell(NSC)%TB_Waals, Scell, NSC, numpar, E_vdW) ! below
         Scell(NSC)%nrg%E_vdW = E_vdW * Na_inv ! [eV/atom]
         
         ! Coulomb potential energy:
         call Coulomb_s(Scell(NSC)%TB_Coul, Scell, NSC, numpar, E_coul) ! below
         Scell(NSC)%nrg%E_coul = E_coul * Na_inv ! [eV/atom]
         
         ! Exponential wall potential energy:
         call Exponential_wall_s(Scell(NSC)%TB_Expwall, Scell, NSC, matter, numpar, E_expwall) ! below
         Scell(NSC)%nrg%E_expwall = E_expwall * Na_inv ! [eV/atom], below

         ! Save the additionall band energy of each atom:
         call band_potential_energy_atom(Scell(NSC)) ! above
      enddo

   endif DO_TB
end subroutine get_pot_nrg


subroutine Exponential_wall_s(TB_Expwall, Scell, NSC, matter, numpar, Pot)
   real(8), intent(out) :: Pot	! Exponential wall energy [eV]
   class(TB_Exp_wall), allocatable, dimension(:,:), intent(in) :: TB_Expwall  ! exponential wall
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in) :: matter   ! all material parameters
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   real(8) a
   if (allocated(TB_Expwall)) then ! if we have Exponential wall potential defined
      select type(TB_Expwall)
      type is (TB_Exp_wall_simple)
         call get_Exp_wall_s(TB_Expwall, Scell, NSC, numpar, a)   ! module "Exponential_wall"
      type is (TB_Short_Rep)
         call get_short_range_rep_s(TB_Expwall, Scell, NSC, matter, numpar, a)   ! module "Exponential_wall"
      end select
   else !For this material exponential wall class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   Pot = a ! [eV]
end subroutine Exponential_wall_s



subroutine Coulomb_s(TB_Coul, Scell, NSC, numpar, Coulomb_out)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_Coulomb), dimension(:,:), allocatable, intent(inout):: TB_Coul ! Coulomb parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: Coulomb_out ! Coulomb energy [eV]
   !--------------
   real(8) a
   if (allocated(TB_Coul)) then ! if we have Coulomb potential defined
      select type(TB_Coul)
      type is (TB_Coulomb_cut)
         call get_Coulomb_s(TB_Coul, Scell, NSC, numpar, a) ! Coulomb energy, module "Coulomb"
      end select
   else !For this material Coulomb class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   Coulomb_out = a ! [eV]
end subroutine Coulomb_s


subroutine vdW_s(TB_Waals, Scell, NSC, numpar, out_vdW_s)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_vdW), dimension(:,:), allocatable, intent(inout):: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8) :: out_vdW_s ! van der Waals energy [eV]
   !---------
   real(8) a

   if (allocated(TB_Waals)) then ! if we have vdW potential defined
!       select type(TB_Waals)
!       type is (TB_vdW_Girifalco)
         call get_vdW_s(TB_Waals, Scell, NSC, numpar, a)   ! van der Waals energy, module "Van_der_Waals"
!       type is (TB_vdW_Dumitrica) ! UNFINISHED! DO NOT USE!
!          call get_vdW_s_D(TB_Waals, Scell, NSC, numpar, a)   ! van der Waals energy, module "Van_der_Waals"
!       end select
   else !For this material vdW class is undefined
      a = 0.0d0 ! no energy for no potential
   endif
   out_vdW_s = a ! [eV]
end subroutine vdW_s



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


subroutine Erep_s(TB_Repuls, Scell, NSC, numpar, E_rep)   ! repulsive energy as a function of a distance
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_repulsive), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: E_rep  ! [eV] repulsive energy
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

   type is (TB_Rep_DFTB_no)
      call get_Erep_s_DFTB_no(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_DFTB"

   type is (TB_Rep_3TB)
      call get_Erep_s_3TB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_3TB"

   type is (TB_Rep_BOP)
      call get_Erep_s_BOP(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_BOP"

   type is (TB_Rep_xTB)
!       call get_Erep_s_xTB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_xTB"
   end select
   E_rep = a
end subroutine Erep_s



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
      endif

      call get_new_energies(Scell, matter, numpar, time, Err) ! module "TB"
   enddo
end subroutine update_nrg_after_change



subroutine get_new_energies(Scell, matter, numpar, t, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(in) :: t ! current timestep [fs]
   type(Error_handling), intent(inout) :: Err ! error save
   !=========================================
   ! Get potential energy of the supercell (Repulsive part is defined inside):
   call get_pot_nrg(Scell, matter, numpar) ! see above

   ! Update electron energy according to updated distribution:
   if (numpar%scc) then ! energy for SCC is defined by for H_0:
!       call get_low_e_energy(Scell, matter, numpar) ! module "Electron_tools"
      call get_low_e_energy(Scell, matter) ! module "Electron_tools"
   else  ! no SCC, the usual expression then
      call get_low_e_energy(Scell, matter) ! module "Electron_tools"
   endif
   ! Get kinetic energy of atoms and of the supercell (for Parrinello-Rahman):
   call get_Ekin(Scell, matter) ! module "Atomic_tools"

   ! And may be we want to update electron energy:
   call get_total_el_energy(Scell, matter, numpar, t, Err) ! module "Electron_tools"
end subroutine get_new_energies



subroutine Electron_ion_coupling(t, matter, numpar, Scell, Err)
! This subroutine calculates the nonadiabatic electron-ion (electron-phonon) coupling
! and dynamical electronic heat conductivity (if required)
   real(8), intent(in) :: t ! [fs] current time
   type(solid), intent(inout) :: matter ! materil parameters
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err ! error save
   !=========================================
   real(8), dimension(:,:), allocatable :: Mij ! matrix element for electron-ion coupling
   real(8) :: dE_nonadiabat, E0
   integer NSC, ixm, iym, izm, ix, iy, iz, i, j

   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
      SC:do NSC = 1, size(Scell) ! for all supercells
         dE_nonadiabat = 0.0d0
         if ( ( (numpar%Nonadiabat) .or. (numpar%do_kappa_dyn)) .AND. (t .GT. numpar%t_NA)) then ! electron-coupling included
            ! Ensure Fermi distribution (unless nonequilibrium simulation is used):
            call update_fe(Scell, matter, numpar, t, Err) ! module "Electron_tools"
            
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
                     call get_coupling_matrix_elements(numpar%NA_kind, Scell, NSC, ix, iy, iz, ixm, iym, izm, Mij)  ! below
                     
                     !sssssssssssssssssssssss
                     ! Project-specific subroutine:
                     ! Test SAKUREI'S EXPRESSION:
!                     call Landau_vs_Sakurei_test(Mij, Scell(NSC)%Ha, Scell(NSC)%Ha0, Ha_non=Scell(NSC)%H_non, Ha_non0=Scell(NSC)%H_non0, wr=Scell(NSC)%Ei, wr0=Scell(NSC)%Ei0) ! module "Nonadiabatic"
                     !sssssssssssssssssssssss
                    
                     ! Calculate the electron-ion collision integral and energy exchange via it:
                     if (numpar%Nonadiabat) then
                        call Electron_ion_collision_int(Scell(NSC), numpar, Scell(NSC)%nrg, Mij, &
                                 Scell(NSC)%Ei, Scell(NSC)%Ei0, Scell(NSC)%fe, &
                                 dE_nonadiabat, numpar%NA_kind, numpar%DOS_weights, Scell(NSC)%G_ei_partial) ! module "Nonadiabatic"
                     endif

                     ! Save Mij to calculate the dynamical electornic heat conductivity:
!                      if (numpar%do_kappa_dyn) then
!                         Scell(NSC)%Mij = Mij
!                      endif
                  enddo ! iz
               enddo ! iy
            enddo ! ix
            ! Averaging over multiple k-points

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
            call get_G_ei(Scell, NSC, numpar, dE_nonadiabat) ! module "Nonadiabatic"
            
         endif
      enddo SC
   endif DO_TB
end subroutine Electron_ion_coupling



subroutine get_coupling_matrix_elements(NA_kind, Scell, NSC, ix, iy, iz, ixm, iym, izm, Mij)
   integer, intent(in) :: NA_kind   ! kind of model for matrix elements
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   real(8), dimension(:,:), allocatable, intent(inout) :: Mij ! matrix element for electron-ion coupling
   integer, intent(in) :: NSC, ixm, iym, izm, ix, iy, iz
   !=====================================

   ! Calculate nonadiabatic-coupling matrix element:
   ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types
      ! Different expressions for orthogonal and non-orthogonal bases:
      select type(ARRAY)
      type is (TB_H_Pettifor) ! TB parametrization according to Pettifor: orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            !call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij) ! module "Nonadiabatic"
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni: orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            !call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij) ! module "Nonadiabatic"
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      type is (TB_H_Fu)  ! TB parametrization accroding to Fu: orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      type is (TB_H_NRL)  ! TB parametrization accroding to NRL method: non-orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind, Scell(NSC)%Sij) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB: non-orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind, Sij=Scell(NSC)%Sij) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      type is (TB_H_3TB)  ! TB parametrization accroding to 3TB: non-orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind, Sij=Scell(NSC)%Sij) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      type is (TB_H_xTB)  ! TB parametrization accroding to xTB: non-orthogonal
         if ( (ix == 1) .and. (iy == 1) .and. (iz == 1) ) then ! Gamma-point, real Hamiltonian is sufficient:
            call Electron_ion_coupling_Mij(Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%Ha0, Mij, NA_kind, Sij=Scell(NSC)%Sij) ! module "Nonadiabatic"
         else ! non-gamma k-point, complex Hamiltonian needed:
            call Electron_ion_coupling_Mij_complex(Scell(NSC)%H_non, Scell(NSC)%CHa, Scell(NSC)%CHa0, Mij, ix, iy, iz, ixm, iym, izm, NA_kind) ! module "Nonadiabatic"
         endif
      end select
   END ASSOCIATE
end subroutine get_coupling_matrix_elements


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


subroutine Construct_Aij_x_En(Ha, distre, En, Aij, H_ssc_1_per_S)
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha 	! diagonilized Hamiltonian/eigenvectors
   REAL(8), DIMENSION(:), INTENT(in) :: distre 	! electron distribution function
   REAL(8), DIMENSION(:), INTENT(in) :: En	! eigenstates [eV]
   REAL(8), DIMENSION(:,:), allocatable, INTENT(inout) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis
   REAL(8), DIMENSION(:,:), INTENT(in), optional :: H_ssc_1_per_S ! contribution from SCC perturbation, Eq.(23) [1]
   ! [1] Elstner et al., PRB 58, 7260 (1998)
   !----------------------------------------
   integer i, k, N
   real(8), dimension(:,:), allocatable :: Atemp
   N = size(distre)
   if (.not. allocated(Aij)) allocate(Aij(N,N))
   allocate(Atemp(N,N))

   if (present(H_ssc_1_per_S)) then ! include SCC H_1 term:
      !$omp PARALLEL private(i)
      !$omp do
      do i = 1, N
         Atemp(i,:) = distre(:)*(En(:) - H_ssc_1_per_S(i,:))*Ha(i,:)   ! original [C 0]
         !Atemp(i,:) = distre(:)*(En(:) + H_ssc_1_per_S(i,:))*Ha(i,:) ! test [C 1]
         !Atemp(i,:) = distre(:)*(En(:))*Ha(i,:) ! test [C 2]
      enddo ! i
      !$omp end do
      !$omp end parallel
   else ! no SSC term:
      !$omp PARALLEL private(i)
      !$omp do
      do i = 1, N
         Atemp(i,:) = distre(:)*En(:)*Ha(i,:)
      enddo ! i
      !$omp end do
      !$omp end parallel
   endif

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


subroutine Get_configurational_temperature(Scell, numpar, matter)
   type(Super_cell), dimension(:), intent(inout) :: Scell	! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(solid), intent(in), target :: matter	! materil parameters
   !--------------------------------------------
   real(8), dimension(:,:), allocatable :: F, dF	! forces and derivatives
   real(8), dimension(:,:), allocatable :: Frep, Fatr, dFrep, dFatr	! forces and derivatives [eV/A], [eV/A^2]
   real(8) :: Tconf, F_sum, dF_sum, acc(3), Ftest(3), dF_temp, Te, F_sum2, dF_sum2
   integer :: Nat, i
   real(8), pointer :: Mass

   Nat = size(Scell(1)%MDAtoms)	! number of atoms
   allocate(F(3,Nat))
   allocate(dF(3,Nat))
   F = 0.0d0
   dF = 0.0d0
   
   ! Forces and their derivatives:
   call get_derivatives_and_forces_r(Scell, numpar, matter, F, dF, Frep, Fatr, dFrep, dFatr)	! see below

   ! Configurational temperature:
   ! 1) Parameters to be used for evaluation of the configurational temperatures:
   !F_sum = SUM( (F(1,:)*F(1,:) + F(2,:)*F(2,:) + F(3,:)*F(3,:)) )
   dF_sum = SUM( (dF(1,:) + dF(2,:) + dF(3,:)) )

!    ! This definition only works for Ta=Te, the standard approximation; in our case, does not work, see corrected expression below!
!    if (abs(dF_sum) <= abs(F_sum) * 1.0d-10) then ! undifined, or infinite
!       Tconf = 0.0d0  ! [K]
!    else  ! defined:
!       Tconf = F_sum / dF_sum    ! [eV]
!       Tconf = Tconf*g_kb	! [eV] -> [K]
!    endif
   !write(*,'(a,f,f,f)') '1:', F_sum, dF_sum, Tconf

   ! 2) Configurational temperature (B=F) derived definition for Ta /= Te:
   F_sum = SUM( (Frep(1,:) + 0.5d0*Fatr(1,:))*F(1,:) + (Frep(2,:) + 0.5d0*Fatr(2,:))*F(2,:) + (Frep(3,:) + 0.5d0*Fatr(3,:))*F(3,:) )
   dF_temp = dF_sum
   ! Make sure the case of Te<<Egap in dielectrics makes some sense:
   Te = max( Scell(1)%TeeV, 0.1d0*Scell(1)%E_gap ) ! exclude Te->0
   !dF_sum = dF_temp + 0.5d0*SUM( F(1,:)*Fatr(1,:) + F(2,:)*Fatr(2,:) + F(3,:)*Fatr(3,:) ) / Scell(1)%TeeV
   dF_sum = dF_temp + 0.5d0*SUM( F(1,:)*Fatr(1,:) + F(2,:)*Fatr(2,:) + F(3,:)*Fatr(3,:) ) / Te
   if (abs(dF_sum) <= abs(F_sum) * 1.0d-10) then ! undifined, or infinite
      Tconf = 0.0d0  ! [K]
   else  ! defined:
      Tconf = F_sum / dF_sum    ! [eV]
      Tconf = Tconf*g_kb	! [eV] -> [K]
   endif
   Scell(1)%Tconf = Tconf
   !write(*,'(a,f,f,f,f)') '2:', F_sum, dF_temp, 0.5d0*SUM( F(1,:)*Fatr(1,:) + F(2,:)*Fatr(2,:) + F(3,:)*Fatr(3,:) ) / Scell(1)%TeeV, Tconf

   ! 3) Hyperconfigurational (second moment configurational) temperature (B=F^2*F):
   F_sum2 = SUM( (F(1,:)*F(1,:) + F(2,:)*F(2,:) + F(3,:)*F(3,:)) * &
               ( (Frep(1,:) + 0.5d0*Fatr(1,:))*F(1,:) + (Frep(2,:) + 0.5d0*Fatr(2,:))*F(2,:) + (Frep(3,:) + 0.5d0*Fatr(3,:))*F(3,:) ) )
   dF_temp = SUM( (3.0d0*F(1,:)*F(1,:) +       F(2,:)*F(2,:) +       F(3,:)*F(3,:)) * dF(1,:) + &
                  (F(1,:)*F(1,:)       + 3.0d0*F(2,:)*F(2,:) +       F(3,:)*F(3,:)) * dF(2,:) + &
                  (F(1,:)*F(1,:)       +       F(2,:)*F(2,:) + 3.0d0*F(3,:)*F(3,:)) * dF(3,:) )
   dF_sum2 = dF_temp + 0.5d0*SUM( (F(1,:)*F(1,:) + F(2,:)*F(2,:) + F(3,:)*F(3,:)) * &
                                       F(1,:)*Fatr(1,:) + F(2,:)*Fatr(2,:) + F(3,:)*Fatr(3,:) ) / Te
   if (abs(dF_sum2) <= abs(F_sum2) * 1.0d-10) then ! undifined, or infinite
      Scell(1)%Tconf2 = 0.0d0  ! [K]
   else  ! defined:
      Scell(1)%Tconf2 = F_sum2 / dF_sum2    ! [eV]
      Scell(1)%Tconf2 = Scell(1)%Tconf2*g_kb	! [eV] -> [K]
   endif


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Testing:
!    do i = 1, Nat
!       ! Convert acceleration into SI units:
!       acc(:) = Scell(1)%MDAtoms(i)%accel(:) * 1.0d20 ! [A/fs^2] -> [m/s^2]
!       Mass => matter%Atoms(Scell(1)%MDatoms(i)%KOA)%Ma ! atomic mass [kg]
!
!       ! Get the force:
!       Ftest(:) = Mass * acc(:) ! [N]
!
!       Ftest(:) = Ftest(:) * 1.0d-10 / g_e ! [eV/A]
!       write(*,'(i4, es,es,es,es,es,es)') i, Ftest(:), -F(:,i)   ! Tested, correct
!    enddo
    !print*, '===================='
    !print*, Tconf, dble(Nat)/SUM( (dF(1,:)+dF(2,:)+dF(3,:)) / (F(1,:)*F(1,:) + F(2,:)*F(2,:) + F(3,:)*F(3,:)) )*g_kb
    !print*, '===================='
    !pause 'Tconfig'
   
   ! Clean up:
   deallocate(F, dF)
end subroutine Get_configurational_temperature


subroutine get_derivatives_and_forces_r(Scell, numpar, matter, F, dF, Frep_out, Fatr_out, dFrep_out, dFatr_out)
   type(Super_cell), dimension(:), intent(inout) :: Scell	! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(solid), intent(in) :: matter	! materil parameters
   real(8), dimension(:,:), allocatable, intent(inout):: F, dF	! forces and derivatives [eV/A], [eV/A^2]
   real(8), dimension(:,:), allocatable, intent(inout), optional :: Frep_out, Fatr_out, dFrep_out, dFatr_out	! forces and derivatives [eV/A], [eV/A^2]
   !-----------------------------------------------------------
   real(8), dimension(:,:), allocatable :: Frep, Fatr, dFrep, dFatr  ! forces and derivatives [eV/A], [eV/A^2]
   real(8), dimension(:,:,:), allocatable :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), allocatable :: M_dVs ! matrix of functions dVs
   real(8), dimension(:,:,:), allocatable :: M_d2Vs ! matrix of functions d2Vs
   real(8), dimension(:,:,:), allocatable :: M_cos	! matrix of directional cosines
   real(8), dimension(:,:), allocatable :: Frep_vdW, dFrep_vdW      ! vdW forces and derivatives [eV/A], [eV/A^2]
   real(8), dimension(:,:), allocatable :: Frep_Coul, dFrep_Coul    ! Coulomb forces and derivatives [eV/A], [eV/A^2]
   real(8), dimension(:,:), allocatable :: Frep_wall, dFrep_wall          ! Short-range repulsive forces and derivatives [eV/A], [eV/A^2]
   integer :: i, N
   
   N = size(Scell(1)%MDAtoms)	! number of atoms
   if (.not.allocated(F)) allocate(F(3,N))
   if (.not.allocated(dF)) allocate(dF(3,N))
   F = 0.0d0   ! restart
   dF = 0.0d0  ! restart
   
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
         !print*, 'Configurational temperature calculations are not implemented for Molteni: Attractive'
      type is (TB_H_Fu)  ! TB parametrization accroding to Fu
         !print*, 'Configurational temperature calculations are not implemented for Fu: Attractive'
      type is (TB_H_NRL)  ! TB parametrization accroding to NRL
         !print*, 'Configurational temperature calculations are not implemented for NRL: Attractive'
      type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB
         !print*, 'Configurational temperature calculations are not implemented for DFTB: Attractive'
      type is (TB_H_3TB)  ! TB parametrization accroding to 3TB
         !print*, 'Configurational temperature calculations are not implemented for 3TB: Attractive'
      type is (TB_H_xTB)  ! TB parametrization accroding to xTB
         !print*, 'Configurational temperature calculations are not implemented for xTB: Attractive'
      end select
   END ASSOCIATE


   ! Repulsive TB Hamiltonian part:
   ASSOCIATE (ARRAY2 => Scell(1)%TB_Repuls(:,:))
      select type(ARRAY2)
      type is (TB_Rep_Pettifor) ! TB parametrization according to Pettifor
         ! Second derivatives of the repulsive energy by r:
         call dE2rep_dr2(ARRAY2, Scell(1)%MDAtoms, Scell, numpar, Frep, dFrep) ! module "TB_Pettifor"
      type is (TB_Rep_Molteni)  ! TB parametrization accroding to Molteni
         !print*, 'Configurational temperature calculations are not implemented for Molteni: Repulsive'
      type is (TB_Rep_Fu)  ! TB parametrization accroding to Fu
         !print*, 'Configurational temperature calculations are not implemented for Fu: Repulsive'
      type is (TB_Rep_NRL)  ! TB parametrization accroding to NRL
         !print*, 'Configurational temperature calculations are not implemented for NRL: Repulsive'
      type is (TB_Rep_DFTB)  ! TB parametrization accroding to DFTB
         !print*, 'Configurational temperature calculations are not implemented for DFTB: Repulsive'
      type is (TB_Rep_3TB)  ! TB parametrization accroding to 3TB
         !print*, 'Configurational temperature calculations are not implemented for 3TB: Repulsive'
      type is (TB_Rep_BOP)  ! TB parametrization accroding to BOP
         !print*, 'Configurational temperature calculations are not implemented for BOP: Repulsive'
      type is (TB_Rep_xTB)  ! TB parametrization accroding to xTB
         !print*, 'Configurational temperature calculations are not implemented for xTB: Repulsive'
      end select
   END ASSOCIATE !


   !cccccccccccccccccccccccccccccccccccccccccccccc
   ! Classical potentials contributions:
   ! SCC Coulomb contribution:
   !call Coulomb_force_from_SCC(numpar, matter, Scell, NSC) ! NOT READY
   !print*, 'before d_vdW_forces'

   ! van der Waals forces:
   call d_vdW_forces(Scell, 1, numpar, Frep_vdW, dFrep_vdW) ! module "Van_der_Waals"
   !print*, 'before d_Coulomb_forces'

   ! Coulomb potential part for modelling Coulomb explosion of a finite system:
   call d_Coulomb_forces(Scell, 1, numpar, Frep_Coul, dFrep_Coul) ! module "Coulomb"
   !print*, 'before d_Exponential_wall_forces'

   ! Exponential wall potential part:
   call d_Exponential_wall_forces(Scell, 1, matter, numpar, Frep_wall, dFrep_wall) ! module "Exponential_wall"

   !print*, 'get_derivatives_and_forces_r:'
   !print*, maxval(Frep), maxval(Frep_vdW), maxval(Frep_Coul), maxval(Frep_wall)

   ! Add all the optional forces:
   if (allocated(Frep)) then  ! if they were calculated
      Frep = Frep + Frep_vdW + Frep_Coul + Frep_wall
   endif
   if (allocated(dFrep)) then ! if they were calculated
      dFrep = dFrep + dFrep_vdW + dFrep_Coul + dFrep_wall
   endif
   !cccccccccccccccccccccccccccccccccccccccccccccc

   
   ! Combine attractive and repulsive parts of forces and derivatives:
   F = Frep + Fatr      ! forces [eV/A]
   dF = dFrep + dFatr   ! derivatives of forces [eV/A^2]


   ! Save partial contributions:
   if (present(Frep_out)) then   ! if requested
      if (.not.allocated(Frep_out)) allocate(Frep_out(3,N))
      if (allocated(Frep)) then  ! if they were calculated
         Frep_out = Frep
      else
         Frep_out = 0.0d0
      endif
   endif
   if (present(Fatr_out)) then ! if requested
      if (.not.allocated(Fatr_out)) allocate(Fatr_out(3,N))
      if (allocated(Fatr)) then  ! if they were calculated
         Fatr_out = Fatr
      else
         Fatr_out = 0.0d0
      endif
   endif
   if (present(dFrep_out)) then ! if requested
      if (.not.allocated(dFrep_out)) allocate(dFrep_out(3,N))
      if (allocated(dFrep)) then ! if they were calculated
         dFrep_out = dFrep
      else
         dFrep_out = 0.0d0
      endif
   endif
   if (present(dFatr_out)) then ! if requested
      if (.not.allocated(dFatr_out)) allocate(dFatr_out(3,N))
      if (allocated(dFatr)) then ! if they were calculated
         dFatr_out = dFatr
      else
         dFatr_out = 0.0d0
      endif
   endif

   ! Clean up:
   if (allocated(Frep)) deallocate(Frep)
   if (allocated(Fatr)) deallocate(Fatr)
   if (allocated(dFrep)) deallocate(dFrep)
   if (allocated(dFatr)) deallocate(dFatr)
   if (allocated(M_Vs)) deallocate(M_Vs)
   if (allocated(M_dVs)) deallocate(M_dVs)
   if (allocated(M_d2Vs)) deallocate(M_d2Vs)
   if (allocated(M_cos)) deallocate(M_cos)
   if (allocated(Frep_vdW)) deallocate(Frep_vdW)
   if (allocated(dFrep_vdW)) deallocate(dFrep_vdW)
   if (allocated(Frep_Coul)) deallocate(Frep_Coul)
   if (allocated(dFrep_Coul)) deallocate(dFrep_Coul)
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
            x => Scell%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
            y => Scell%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
            z => Scell%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
            r => Scell%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
            
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
   integer :: NSC, i, j
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
   ! For vdW and Coulomb contribution:
   real(8), dimension(:,:,:), allocatable :: Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, XijSupce, YijSupce, ZijSupce
   integer :: Nx, Ny, Nz

   
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
         if (numpar%verbose) print*, 'Get_pressure 3TB : Construct_M_x1 succesful'
         call get_Mjs_factors(numpar%N_basis_size, Scell(NSC), M_lmn, Mjs)   ! module "TB_3TB"
         if (numpar%verbose) print*, 'Get_pressure 3TB : get_Mjs_factors succesful'
         call Construct_Vij_3TB(numpar, ARRAY, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_Lag_exp, M_d_Lag_exp)	! module "TB_3TB"
         if (numpar%verbose) print*, 'Get_pressure 3TB : Construct_Vij_3TB succesful'
         call Construct_Aij_x_En(Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%Ei, M_Aij_x_Ei) ! see below
         if (numpar%verbose) print*, 'Get_pressure 3TB : Construct_Aij_x_En succesful'
         call Attract_TB_Forces_Press_3TB(Scell, NSC, ARRAY, numpar, Scell(NSC)%Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Aij_x_Ei, Mjs, M_Lag_exp, M_d_Lag_exp) ! module "TB_3TB"
         if (numpar%verbose) print*, 'Get_pressure 3TB : Attract_TB_Forces_Press_3TB succesful'
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
      type is (TB_Rep_DFTB_no) ! TB parametrization according to DFTB
         call dErdr_Pressure_s_DFTB_no(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_DFTB"
      type is (TB_Rep_3TB) ! TB parametrization according to 3TB
         call dErdr_Pressure_s_3TB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_3TB"
      type is (TB_Rep_BOP) ! TB parametrization according to BOP
!          call dErdr_Pressure_s_DFTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_BOP"
      type is (TB_Rep_xTB) ! TB parametrization according to xTB
!          call dErdr_Pressure_s_xTB(ARRAY2, Scell, NSC, numpar) ! derivatives of the repulsive energy by h; module "TB_xTB"
      end select
   END ASSOCIATE

   ! Other contributions to forces/pressure:
   ! van der Waals (vdW) part with TB Hamiltonian:
   if (allocated(Scell(NSC)%TB_Waals)) then ! if we have vdW potential defined
!       ASSOCIATE (ARRAY2 => Scell(NSC)%TB_Waals(:,:))
!       select type (ARRAY2)
!       type is (TB_vdW_Girifalco)
         ! Get multipliers used many times into temporary arrays:
         call Construct_B(Scell(NSC)%TB_Waals, Scell, NSC, numpar, Scell(NSC)%MDatoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Van_der_Waals"
         ! Forces for the super-cell:
         call d_Forces_Pressure(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! below
!       type is (TB_vdW_LJ_cut)
!          ! Get multipliers used many times into temporary arrays:
!          call Construct_B(Scell, NSC, numpar, Scell(NSC)%MDatoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Van_der_Waals"
!          ! Forces for the super-cell:
!          call  d_Forces_Pressure(Scell(NSC)%MDatoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! below
!       end select
!       END ASSOCIATE
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
   endif


   if (allocated(Scell(NSC)%TB_Coul)) then ! if we have Coulomb potential defined
      ASSOCIATE (ARRAY2 => Scell(NSC)%TB_Coul(:,:))
      select type (ARRAY2)
      type is (TB_Coulomb_cut) ! so far, it is the only type we have
         ! Get multipliers used many times into temporary arrays:
         call Construct_B_C(ARRAY2, Scell, NSC, Scell(NSC)%MDatoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) ! module "Coulomb"
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
      END ASSOCIATE
   endif

   ! Exponential wall potential part:
   if (allocated(Scell(NSC)%TB_Expwall)) then ! if we have short-range potential defined
      ASSOCIATE (ARRAY2 => Scell(NSC)%TB_Expwall)
      select type (ARRAY2)
      type is (TB_Exp_wall_simple)
         ! Forces for the super-cell:
         call d_Exp_wall_Pressure_s(Scell, NSC, ARRAY2, numpar)  ! module "Exponential_wall"
      type is (TB_Short_Rep)
         ! Forces for the super-cell:
         call d_Short_range_Pressure_s(Scell, NSC, ARRAY2, matter, numpar) ! module "Exponential_wall"
      end select
      END ASSOCIATE
   endif
   !cccccccccccccccccccccccccccccccccccccccccccccc

   ! Get forces for the supercell:
   call Potential_super_cell_forces(numpar, Scell, NSC, matter)  ! module "Atomic_tools"

   ! Adding kinetic part and pressure to super-cell forces:
   call super_cell_forces(numpar, Scell, NSC, matter, Scell(NSC)%SCforce, sigma_tensor) ! module "Atomic_tools"
   
   ! Invert sigma tensor:
   call Invers_3x3(sigma_tensor, sigma_inversed, 'Get_pressure')	! module "Algebra_tools"
   
   ! Calculate the stress tensor (factor 1.0e40 is to convert from [kg/A/fs^2] to [kg/m/s^2]):
   stress_tensor(:,:) = Scell(NSC)%SCforce%total(:,:) * matter%W_PR * sigma_inversed(:,:) * 1.0d40	! [kg/m/s^2]
   
   ! OUTPUT, pressure and stress tensor:
   P = 1.0d0/3.0d0 * ( stress_tensor(1,1) + stress_tensor(2,2) + stress_tensor(3,3) )
   if (present(stress_tensor_OUT)) stress_tensor_OUT = stress_tensor
   ! Save potential contributions to the pressure (to use leter for configurational temperature):
   do i = 1,3
      do j = 1,3
            Scell(NSC)%Pot_Stress(j,i) = Scell(NSC)%SCforce%att(i,j) + Scell(NSC)%SCforce%rep(j,i) ! [kg/A/fs^2]
      enddo ! j
   enddo ! i
   Scell(NSC)%Pot_Stress = Scell(NSC)%Pot_Stress * matter%W_PR * sigma_inversed(:,:) * 1.0d40   ! [kg/m/s^2]
   Scell(NSC)%Pot_Pressure = 1.0d0/3.0d0 * (Scell(NSC)%Pot_Stress(1,1) + Scell(NSC)%Pot_Stress(2,2) + Scell(NSC)%Pot_Stress(3,3) )

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
   call vdW_s(TEMP_Scell(1)%TB_Waals, TEMP_Scell, 1, numpar, vdW_nrg) ! vdW energy
   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   coun = 0
   do while (abs(a-b)/abs(b) >= 1.0d-4)
      coun = coun + 1
      Charge = (a+b)/2.0d0
      TEMP_Scell(1)%Q = Charge !set unballanced charge per atom
      ! Coulomb potential part for modelling Coulomb explosion of a finite system:
      call Coulomb_s(TEMP_Scell(1)%TB_Coul, TEMP_Scell, 1, numpar, Coulomb_nrg) ! get Coulomb energy

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
   
   ! Set Coulomb (with soft cot off) parameters:
   ASSOCIATE (ARRAY2 => Scell(1)%TB_Coul(:,:))
      select type (ARRAY2)
         type is (TB_Coulomb_cut)
         allocate(TB_Coulomb_cut::TEMP_Scell%TB_Coul(1,1)) ! make it for Coulomb parametrization
         ASSOCIATE (TEMP_ARRAY => TEMP_Scell%TB_Coul(:,:))
         select type (TEMP_ARRAY)
            type is (TB_Coulomb_cut)
            TEMP_ARRAY = ARRAY2 ! copy from the main super-cell
            ! Replace the Coulomb cot off by the proper value:
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


!================================================
! Obsolete subroutines:

! This subroutine is superceeded by the module "TB_complex"
subroutine get_DOS(numpar, matter, Scell, Err) ! Obsolete subroutine
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function
   type(Solid), intent(in) :: matter     ! material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err	! error save
   !=====================================
   integer :: NSC, Nsiz, Ei_siz, i, n_types
   real(8) :: dE, Estart, Emax
   real(8), dimension(:), allocatable :: Ei_cur

   do NSC = 1, size(Scell) ! for all super-cells
      if (numpar%save_DOS) then	! only calculate DOS if the user chose to do so:
!           call print_time_step('Before DOS:', 1.0, msec=.true.)   ! module "Little_subroutines"
         ! Set grid for DOS:
         if (.not.allocated(Scell(NSC)%DOS)) then	! it's the first time, set it:
            dE = 0.1d0	! [eV] uniform energy grid step for DOS
            Ei_siz = size(Scell(NSC)%Ei)	! the uppermost energy level at the start
            Emax = min(Scell(NSC)%Ei(Ei_siz),100.0)   ! no need to trace levels higher than 100 eV
            Nsiz = CEILING( (Emax - Scell(NSC)%Ei(1) + 20.0d0*numpar%Smear_DOS)/dE )
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

!           print*, 'get_DOS test 1'

         ! Now calculate the DOS:
         select case (ABS(numpar%optic_model))	! use multiple k-points, or only gamma!
         case (2,4,5)   ! multiple k points

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
               call get_DOS_sort(Scell(NSC)%Ei, Scell(NSC)%DOS, numpar%Smear_DOS, Scell(NSC)%partial_DOS, numpar%mask_DOS, Hij = Scell(NSC)%Ha)	! module "Electron_tools"
            case default    ! No need to sort DOS per orbitals
               call get_DOS_sort(Scell(NSC)%Ei, Scell(NSC)%DOS, numpar%Smear_DOS)	! module "Electron_tools"
            endselect
         end select
!           call print_time_step('After DOS:', 1.0, msec=.true.)   ! module "Little_subroutines"
      endif !  (numpar%save_DOS)
   enddo ! NSC = 1, size(Scell) ! for all super-cells
end subroutine get_DOS


! This subroutine is superceeded by the module "TB_complex"
subroutine get_DOS_sort_complex(numpar, Scell, NSC, DOS, smearing, Err, partial_DOS, masks_DOS) ! Obsolete subroutine
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
               call get_complex_Hamiltonian(numpar, Scell, NSC,  CHij, CSij, Ei, kx, ky, kz, Err, .true.)	! see below

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



END MODULE TB
