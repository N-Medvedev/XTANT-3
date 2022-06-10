! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2022 Nikita Medvedev
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
! This module contains subroutines to deal with 3TB hamiltonian: https://github.com/usnistgov/ThreeBodyTB.jl
! [1] https://arxiv.org/pdf/2112.11585.pdf
! [2] https://github.com/usnistgov/ThreeBodyTB.jl/blob/master/src/CalcTB_laguerre.jl

MODULE TB_3TB

use Universal_constants
use Objects
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity, Laguerre_up_to_6, d_Laguerre_up_to_6
use Atomic_tools, only : Reciproc_rel_to_abs, shortest_distance
use Little_subroutines, only : linear_interpolation, Fermi_function, d_Fermi_function, Find_in_array_monoton
use Electron_tools, only : find_band_gap
use TB_Koster_Slater
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, &
                  test_orthogonalization_c, Loewdin_Orthogonalization, Loewdin_Orthogonalization_c
use TB_DFTB, only : identify_DFTB_basis_size, identify_DFTB_orbitals_per_atom, Get_overlap_S_matrix_DFTB
use Dealing_with_3TB, only: find_3bdy_ind


implicit none


real(8), parameter :: m_a = 2.0d0   ! exponential decay parameter according to [1]



 contains



subroutine Construct_Vij_3TB(numpar, TB, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_Lag_exp)
   type(Numerics_param), intent(in), target :: numpar 	! all numerical parameters
   type(TB_H_3TB), dimension(:,:), intent(in), target :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Vij	! matrix of hoppings for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dVij	! matrix of derivatives of hoppings for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_SVij	! S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dSVij  ! derivatives of S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)

   !----------------------------------
   real(8) :: x, y, z, r, fcut, d_fcut, Fermi, dFermi, Laguer(6), d_Laguer(6), rIK, rJK, exp_ad
   integer :: i, j, atom_2, ki, N_bs, ihop, atom_3, k, at_ind, sh1, sh2
   real(8), pointer :: rm
   integer, pointer :: nat, m, KOA1, KOA2, KOA3, mm

   nat => Scell(NSC)%Na	! number of atoms in the supercell
   ! number of hopping integrals for this basis set in 3TB:
   N_bs = identify_DFTB_basis_size(numpar%N_basis_size)  ! module "TB_DFTB"

   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,N_bs))	    ! each pair of atoms, all  V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,N_bs))	! each pair of atoms, all  dV functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,N_bs))	! each pair of atoms, all S functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,N_bs))	! each pair of atoms, all dS functions
   if (.not.allocated(M_Lag_exp)) allocate(M_Lag_exp(nat,nat,6))  ! each pair of atoms, 6 polynomials

   !$OMP WORKSHARE
   M_Vij = 0.0d0
   M_dVij = 0.0d0
   M_SVij = 0.0d0
   M_dSVij = 0.0d0
   !$OMP END WORKSHARE


   ! 2-body interaction terms:

   ! Construct matrix of all the radial functions for each pair of atoms:
!$omp PARALLEL
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, Laguer, d_Laguer, exp_ad, ihop, Fermi, dFermi)
   AT1:do j = 1,nat	! all atoms

      KOA1 => Scell(NSC)%MDatoms(j)%KOA

      m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms

         KOA2 => Scell(NSC)%MDatoms(i)%KOA

         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! Get Laguerre polynomials and its derivatives to be reused below:
         call get_Laguerres(r, Laguer) ! below
         call get_d_Laguerres(r, d_Laguer) ! below

         ! Get the exponential term:
         exp_ad = Laguerre_exponent(r) ! below
         ! Include it into the solution:
         Laguer(:) = Laguer(:) * exp_ad

         ! Save [Laguerre * exp(-a*r_ij)]
         M_Lag_exp(j,i,:) = Laguer(:)

         ! All radial functions for Hamiltonian (functions below):
         M_Vij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, 1, 5)   ! (s s sigma)
         M_SVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, 1, 6)  ! (s s sigma)
         select case (numpar%N_basis_size)
         case (1)    ! sp3
            M_Vij(j,i,2) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, 2, 5)   ! (s p sigma)
            M_SVij(j,i,2) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, 2, 6) ! (s p sigma)
            M_Vij(j,i,3) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, 4, 5)   ! (p p sigma)
            M_SVij(j,i,3) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, 4, 6) ! (p p sigma)
            M_Vij(j,i,4) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, 5, 5)   ! (p p pi)
            M_SVij(j,i,4) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, 5, 6) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_Vij(j,i,ihop) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, ihop, 5)
               M_SVij(j,i,ihop) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, ihop, 6)
            enddo
         endselect

         ! All derivatives of the radial functions and the radial functions for Overlap matrix:
         M_dVij(j,i,1) = d_radial_function_3TB(r, M_Vij(j,i,1), TB(KOA1,KOA2)%Vrfx, d_Laguer, 1, 5)   ! (s s sigma)
         M_dSVij(j,i,1) = d_radial_function_3TB(r, M_SVij(j,i,1), TB(KOA1,KOA2)%Srfx, d_Laguer, 1, 6) ! (s s sigma)
         select case (numpar%N_basis_size)
         case (1)    ! sp3
            M_dVij(j,i,2) = d_radial_function_3TB(r, M_Vij(j,i,2), TB(KOA1,KOA2)%Vrfx, d_Laguer, 2, 5)   ! (s p sigma)
            M_dSVij(j,i,2) = d_radial_function_3TB(r, M_SVij(j,i,2), TB(KOA1,KOA2)%Srfx, d_Laguer, 2, 6) ! (s p sigma)
            M_dVij(j,i,3) = d_radial_function_3TB(r, M_Vij(j,i,3), TB(KOA1,KOA2)%Vrfx, d_Laguer, 4, 5)   ! (p p sigma)
            M_dSVij(j,i,3) = d_radial_function_3TB(r, M_SVij(j,i,3), TB(KOA1,KOA2)%Srfx, d_Laguer, 4, 6) ! (p p sigma)
            M_dVij(j,i,4) = d_radial_function_3TB(r, M_Vij(j,i,4), TB(KOA1,KOA2)%Vrfx, d_Laguer, 5, 5)   ! (p p pi)
            M_dSVij(j,i,4) = d_radial_function_3TB(r, M_SVij(j,i,4), TB(KOA1,KOA2)%Srfx, d_Laguer, 5, 6) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_dVij(j,i,ihop) = d_radial_function_3TB(r, M_Vij(j,i,ihop), TB(KOA1,KOA2)%Vrfx, d_Laguer, ihop, 5)
               M_dSVij(j,i,ihop) = d_radial_function_3TB(r, M_SVij(j,i,ihop), TB(KOA1,KOA2)%Srfx, d_Laguer, ihop, 6)
            enddo
         endselect

         ! Add Fermi-like smoothing approach to zero:
         Fermi = Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         dFermi = d_Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         M_dVij(j,i,:) = M_dVij(j,i,:)*Fermi + M_Vij(j,i,:)*dFermi
         M_dSVij(j,i,:) = M_dSVij(j,i,:)*Fermi + M_SVij(j,i,:)*dFermi
         M_Vij(j,i,:) = M_Vij(j,i,:)*Fermi
         M_SVij(j,i,:) = M_SVij(j,i,:)*Fermi
      enddo AT2
   enddo AT1
!$omp end do
!$omp END PARALLEL

   nullify(rm, nat, m, KOA1, KOA2, KOA3)	! clean up at the end
end subroutine Construct_Vij_3TB



! Tight Binding Hamiltonian within 3TB parametrization:
subroutine construct_TB_H_3TB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Scell, NSC, Err)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! materil parameters
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   type(Super_cell), dimension(:), intent(inout) :: Scell		! supercell with all the atoms as one object
   integer, intent(in) :: NSC		! number of supercell
   type(Error_handling), intent(inout) :: Err	! error save
   character(200) :: Error_descript
   Error_descript = ''

!$OMP WORKSHARE
   Scell(NSC)%Ha0 = Scell(NSC)%Ha	! save Hamiltonial from previous time-step
   Scell(NSC)%Ei0 = Scell(NSC)%Ei	! save energy levels for the next timestep
   Scell(NSC)%H_non0 = Scell(NSC)%H_non	! save non-diagonalized Hamiltonian from last time-step
!$OMP END WORKSHARE

    ! Construct TB Hamiltonian (with DFTB parameters),  orthogonalize it,  and then solve generalized eigenvalue problem:
   call Hamil_tot_3TB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Err)	! see below

!    ! Test (comment out for release):
!    call test_nonorthogonal_solution(Scell(NSC)) ! module "TB_NRL"

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_3TB



subroutine Hamil_tot_3TB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------------------
   real(8), dimension(:,:), allocatable :: Hij	 ! Hamiltonian
   real(8), dimension(:,:), allocatable :: Sij  ! Overlap
   real(8), dimension(size(Scell(NSC)%Ha,1)) :: Evec, EvecS
   real(8), dimension(:,:), allocatable :: Hij1, Sij1
   integer :: nat, Nsiz, n_orb
   integer, target :: j
   integer :: j1, i1, k, l, atom_2, FN, i
   real(8) :: temp, epsylon, Ev, SH_1
   real(8), pointer :: x, y, z
   integer, pointer :: KOA1, KOA2, m
   character(200) :: Error_descript

   Error_descript = ''
   epsylon = 1d-12	! acceptable tolerance : how small an overlap integral can be, before we set it to zero
   ! size of the basis set per atom, below:
   n_orb = identify_DFTB_orbitals_per_atom(numpar%N_basis_size)  ! module "TB_DFTB"
   nat = Scell(NSC)%Na  ! number of atoms in the supercell
   Nsiz = size(Scell(NSC)%Ha,1) ! size of the total basis set

   if (.not.allocated(Sij)) allocate(Sij(Nsiz,Nsiz))
   if (.not.allocated(Hij)) allocate(Hij(Nsiz,Nsiz))
   Sij = 0.0d0
   Hij = 0.0d0

!    print*, 'Hamil_tot_3TB test 0'

   !-----------------------------------
   ! 1) Construct non-orthogonal Hamiltonian H and Overlap matrix S in 2 steps:
!$omp parallel private(j, m, atom_2, i, KOA1, KOA2, j1, l, i1, k, Hij1, Sij1)
if (.not.allocated(Hij1)) allocate(Hij1(n_orb,n_orb))
if (.not.allocated(Sij1)) allocate(Sij1(n_orb,n_orb))
!$omp do
   do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)

      do atom_2 = 0,m ! do only for atoms close to that one

         if (atom_2 == 0) then ! the same atom
            i = j
         else  ! different atoms
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         endif

         IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
            KOA1 => Scell(NSC)%MDatoms(j)%KOA
            KOA2 => Scell(NSC)%MDatoms(i)%KOA
            ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
            ! Contruct a block-hamiltonian:
            call Hamilton_one_3TB(numpar%N_basis_size, Scell(NSC), j, i, TB_Hamil, Hij1, &
                                    M_Vij, M_Lag_exp, M_lmn, Mjs)   ! below

            ! Construct overlap matrix for this pair of atoms:
            call Get_overlap_S_matrix_DFTB(numpar%N_basis_size, j, i, Sij1, M_SVij, M_lmn)  ! module "TB_DFTB"

            do j1 = 1,n_orb ! all orbitals
               l = (j-1)*n_orb+j1
               do i1 = 1,n_orb ! all orbitals
                  k = (i-1)*n_orb+i1
                  Hij(k,l) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian
                  Sij(k,l) = Sij1(i1,j1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
                  if (ABS(Sij(k,l)) <= epsylon) Sij(k,l) = 0.0d0
               enddo ! i1
            enddo ! j1
         endif IJ
      enddo ! j
   enddo ! i
!$omp end do
deallocate(Hij1, Sij1)
!$omp end parallel

   ! b) Construct lower triangle - use symmetry:
!$omp parallel
!$omp do  private(j, m, atom_2, i, j1, l, i1, k)
   do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         if (i < j) then   ! lower triangle, copy from the upper one
            do j1 = 1,n_orb ! all orbitals
               l = (j-1)*n_orb+j1
               do i1 = 1,n_orb ! all orbitals
                  k = (i-1)*n_orb+i1
                  Hij(k,l) = Hij(l,k)
                  Sij(k,l) =  Sij(l,k)
               enddo ! i1
            enddo ! j1
         endif
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel


! print*, 'Hij(k,l)', Hij(1,2)   ! Printout for testing


   ! 2)    ! Save the non-orthogonalized Hamiltonian:
   !$OMP WORKSHARE
   Scell(NSC)%H_non = Hij		! nondiagonalized Hamiltonian
   Scell(NSC)%Sij = Sij		! save Overlap matrix
   !$OMP END WORKSHARE

   !-----------------------------------
   ! 3) Orthogonalize the Hamiltonian using Lowedin procidure
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization(Nsiz, Sij, Hij, Err)	! module "TB_NRL"

   !$OMP WORKSHARE
   Scell(NSC)%Hij = Hij ! save orthogonalized but non-diagonalized Hamiltonian
   !$OMP END WORKSHARE

   forall (i = 1:size(Hij,1),  j = 1:size(Hij,2), (ABS(Hij(i,j)) < 1.0d-10))
      Hij(i,j) = 0.0d0
   endforall
   !-----------------------------------
   ! 4) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
   call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript) ! module "Algebra_tools"
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Hamil_tot_DFTB: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif
   Scell(NSC)%Hij_sol = Hij		! eigenvectors of nondiagonalized Hamiltonian

   !-----------------------------------
   ! 5) Convert the eigenvectors back into the non-orthogonal basis:
   call mkl_matrix_mult('N', 'N', Sij, Hij, Scell(NSC)%Ha)	! module "Algebra_tools"


   !ccccccccccccccccccccccccc
   ! Optical coefficients:
   if (numpar%optic_model .EQ. 3) then ! create matrix element:
   ! Optical matrix elements for non-orthogonal TB are taken from:
   ! arXiv:1805.08918v1 -- https://128.84.21.199/abs/1805.08918
      if (.not.allocated(Scell(NSC)%PRRx)) allocate(Scell(NSC)%PRRx(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%PRRy)) allocate(Scell(NSC)%PRRy(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%PRRz)) allocate(Scell(NSC)%PRRz(Nsiz,Nsiz))
      Scell(NSC)%PRRx = 0.0d0
      Scell(NSC)%PRRy = 0.0d0
      Scell(NSC)%PRRz = 0.0d0
      ! Generalized Trani expression for non-orthogonal Hamiltonian:

      !$omp parallel
      !$omp do private(j, m, atom_2, i, x, y, z, j1, l, i1, k, SH_1)
      do j = 1,nat	! all atoms
         m => Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
               x => Scell(NSC)%Near_neighbor_dist(j,atom_2,1) ! at this distance, X
               y => Scell(NSC)%Near_neighbor_dist(j,atom_2,2) ! at this distance, Y
               z => Scell(NSC)%Near_neighbor_dist(j,atom_2,3) ! at this distance, Z
            endif
            if (i > 0) then ! if there really is a nearest neighbor
               do j1 = 1,n_orb	! all orbitals
                  l = (j-1)*n_orb+j1
                  do i1 = 1,n_orb	! all orbitals
                     k = (i-1)*n_orb+i1
                     if (i == j) then ! for the diagonal elements according to Trani:
                        !SH_1 =  0.0d0
                        SH_1 =  0.270d0*(Scell(NSC)%H_non(k,l) - Scell(NSC)%Ei(k) * Scell(NSC)%Sij(k,l))
                        Scell(NSC)%PRRx(k,l) = SH_1
                        Scell(NSC)%PRRy(k,l) = SH_1
                        Scell(NSC)%PRRz(k,l) = SH_1
                     else
                        SH_1 = (Scell(NSC)%H_non(k,l) - Scell(NSC)%Ei(k) * Scell(NSC)%Sij(k,l))
                        Scell(NSC)%PRRx(k,l) = x*SH_1
                        Scell(NSC)%PRRy(k,l) = y*SH_1
                        Scell(NSC)%PRRz(k,l) = z*SH_1
                     endif	! note that the PRR are not diagonal
                     if (Scell(NSC)%PRRx(k,l) .GT. 1d10) print*, i, j, Scell(NSC)%PRRx(k,l)
                     if (Scell(NSC)%PRRy(k,l) .GT. 1d10) print*, i, j, Scell(NSC)%PRRy(k,l)
                     if (Scell(NSC)%PRRz(k,l) .GT. 1d10) print*, i, j, Scell(NSC)%PRRz(k,l)
                  enddo ! i1
               enddo ! j1
            endif ! (i > 0)
         enddo ! j
      enddo ! i
      !$omp end do
      !$omp end parallel
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0

      Scell(NSC)%PRRx = Scell(NSC)%PRRx * temp
      Scell(NSC)%PRRy = Scell(NSC)%PRRy * temp
      Scell(NSC)%PRRz = Scell(NSC)%PRRz * temp
   endif ! (numpar%optic_model .EQ. 3)

   nullify(KOA1, KOA2, m, x, y, z)
   deallocate(Hij, Sij)
end subroutine Hamil_tot_3TB



subroutine Hamilton_one_3TB(basis_ind, Scell, i, j, TB, Hij, M_Vij, M_Lag_exp, M_lmn, Mjs)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: i, j   ! atoms indices
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	! all tight binding parameters
   real(8), dimension(:,:), intent(out) :: Hij  ! hamiltonian, all orbitals in sp3d5 basis set
   real(8), dimension(:,:,:), intent(in) :: M_Vij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   !---------------------------------------/
   if (i == j) then ! Onsite contributions
      call Onsite_3TB(basis_ind, i, Scell, TB, M_Lag_exp, Mjs, Hij)	! subroutine below
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:
      call Hopping_3TB(basis_ind, i, j, Scell, TB, M_Vij(i,j,:), M_Vij(j,i,:), M_lmn(:,i,j), M_Lag_exp, Mjs, Hij)	! subroutine below
   endif
end subroutine Hamilton_one_3TB



subroutine Hopping_3TB(basis_ind, i, j, Scell, TB_Hamil, M_Vij12, M_Vij21, M_lmn, M_Lag_exp, Mjs, ts)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   integer, intent(in) :: i, j ! indices of the atoms
   type(Super_cell), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:), intent(in) :: M_Vij12, M_Vij21	! Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:), intent(in) :: M_lmn	! cosine directions
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of Laguerre*exp for 3-body-interactions
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   real(8), dimension(:,:), intent(out) :: ts	! overlap integerals [eV]
   !=============================================
   real(8), dimension(:), allocatable :: vec_M_SVij12, vec_M_SVij21
   real(8), dimension(size(ts,1),size(ts,2)) :: H_3bdy_part   ! 3-body contribution
   real(8) :: G_IJK(3,3), H_temp(6)
   integer :: k, at_ind, sh1, sh2, atom_3
   integer, pointer :: m, KOA1, KOA2, KOA3

   ! Constructing total hamiltonian includes 2-body and 3-body terms:

   ! 1) Construct the 2-body overlap integrals including angular part:
   select case (basis_ind)
   case (0) ! for s basis set:
      allocate(vec_M_SVij12(1))
      vec_M_SVij12(:) = M_Vij12(:)
      call KS_s(vec_M_SVij12, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   case (1) ! for sp3 basis set:
      allocate(vec_M_SVij12(4))
      vec_M_SVij12(:) = M_Vij12(:)
      allocate(vec_M_SVij21(4))
      vec_M_SVij21(:) = M_Vij21(:)
      call KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   case default ! for sp3d5 basis set:
      allocate(vec_M_SVij12(10))
      vec_M_SVij12(:) = M_Vij12(:)
      allocate(vec_M_SVij21(10))
      vec_M_SVij21(:) = M_Vij21(:)
      call KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   end select
   if (allocated(vec_M_SVij12)) deallocate(vec_M_SVij12)
   if (allocated(vec_M_SVij21)) deallocate(vec_M_SVij21)


   ! 2) Construct the 3-body terms:
   H_3bdy_part = 0.0d0  ! to start with

   KOA1 => Scell%MDatoms(i)%KOA  ! kind of atoms #1

   if (TB_Hamil(KOA1,KOA1)%include_3body) then  ! only if user defined it to include

      ! types of atoms:
      KOA2 => Scell%MDatoms(j)%KOA  ! kind of atoms #2
      m => Scell%Near_neighbor_size(i) ! how many near neighbors

!     goto 5000   ! For testing, exclude 3-body terms

!!!$omp PARALLEL ! it is already inside parallelized region, do not parallelize again!
!!!$omp do private(atom_3, k, KOA3, at_ind, sh1, sh2)
      AT3:do atom_3 = 1,m ! do only for atoms close to that one
         k = Scell%Near_neighbor_list(i,atom_3) ! this is the list of such close atoms
         ! Make sure the third atoms is not the second atom:
         if (k /= j) then
            KOA3 => Scell%MDatoms(k)%KOA   ! kind of atom #3

            ! Find the combination-of-atoms index:
            at_ind = find_3bdy_ind(KOA1, KOA2, KOA3)  ! module "Dealing_with_3TB"

            ! Get the radial function for 3-body interaction:
            do sh1 = 1, 1+basis_ind
               do sh2 = 1, 1+basis_ind
                  G_IJK(sh1,sh2) = TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 1) * M_Lag_exp(i,k,1) * M_Lag_exp(j,k,1) + &
                              TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 2) * M_Lag_exp(i,k,1) * M_Lag_exp(j,k,2) + &
                              TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 3) * M_Lag_exp(i,k,2) * M_Lag_exp(j,k,1) + &
                              TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 4) * M_Lag_exp(i,k,1) * M_Lag_exp(j,k,1) * M_Lag_exp(i,j,1)
               enddo
            enddo
            ! Include angular parts:
            H_3bdy_part(1,1) = H_3bdy_part(1,1) + G_IJK(1,1)  ! s-s * s-s

            if (basis_ind > 0) then ! p3 orbitals:

               H_3bdy_part(1,2:4) = H_3bdy_part(1,2:4) + G_IJK(2,1) * Mjs(j,k,2:4)   ! s-s * p[x,y,z]-s
               H_3bdy_part(2:4,1) = H_3bdy_part(2:4,1) + G_IJK(1,2) * Mjs(i,k,2:4)   ! p[x,y,z]-s * s-s

               ! Calculate repeating part the K-S matrix elements:
               H_temp(1) = G_IJK(2,2) * Mjs(i,k,2) ! G * px-s
               H_temp(2) = G_IJK(2,2) * Mjs(i,k,3) ! G * px-s
               H_temp(3) = G_IJK(2,2) * Mjs(i,k,4) ! G * px-s
               ! Add it into the Hamiltonian part:
               H_3bdy_part(2,2:4) = H_3bdy_part(2,2:4) + H_temp(1) * Mjs(j,k,2:4) ! p[x,y,z]-s
               H_3bdy_part(3,2:4) = H_3bdy_part(3,2:4) + H_temp(2) * Mjs(j,k,2:4) ! p[x,y,z]-s
               H_3bdy_part(4,2:4) = H_3bdy_part(4,2:4) + H_temp(3) * Mjs(j,k,2:4) ! p[x,y,z]-s

               if (basis_ind > 1) then ! d5 orbitals:
                  H_3bdy_part(1,5) = H_3bdy_part(1,5) + G_IJK(1,3) * Mjs(j,k,5) ! s-s * dxy-s
                  H_3bdy_part(1,6) = H_3bdy_part(1,6) + G_IJK(1,3) * Mjs(j,k,6) ! s-s * dxz-s
                  H_3bdy_part(1,7) = H_3bdy_part(1,7) + G_IJK(1,3) * Mjs(j,k,7) ! s-s * dxz-s
                  H_3bdy_part(1,8) = H_3bdy_part(1,8) + G_IJK(1,3) * Mjs(j,k,8) ! s-s * (dx2-y2)-s
                  H_3bdy_part(1,9) = H_3bdy_part(1,9) + G_IJK(1,3) * Mjs(j,k,9) ! s-s * (d3z2-r2)-s

                  H_3bdy_part(5,1) = H_3bdy_part(5,1) + G_IJK(3,1) * Mjs(i,k,5) ! dxy-s * s-s
                  H_3bdy_part(6,1) = H_3bdy_part(6,1) + G_IJK(3,1) * Mjs(i,k,6) ! dxz-s * s-s
                  H_3bdy_part(7,1) = H_3bdy_part(7,1) + G_IJK(3,1) * Mjs(i,k,7) ! dxz-s * s-s
                  H_3bdy_part(8,1) = H_3bdy_part(8,1) + G_IJK(3,1) * Mjs(i,k,8) ! (dx2-y2)-s * s-s
                  H_3bdy_part(9,1) = H_3bdy_part(9,1) + G_IJK(3,1) * Mjs(i,k,9) ! (d3z2-r2)-s * s-s

                  ! Calculate repeating part the K-S matrix elements:
                  H_temp(1) = Mjs(i,k,2) * Mjs(j,k,5) ! px-s * dxy-s
                  H_temp(2) = Mjs(i,k,2) * Mjs(j,k,6) ! px-s * dxz-s
                  H_temp(3) = Mjs(i,k,2) * Mjs(j,k,7) ! px-s * dxz-s
                  H_temp(4) = Mjs(i,k,2) * Mjs(j,k,8) ! px-s * (dx2-y2)-s
                  H_temp(5) = Mjs(i,k,2) * Mjs(j,k,9) ! px-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(2,5:9) = H_3bdy_part(2,5:9) + H_temp(1:5) * G_IJK(2,3)
                  H_3bdy_part(5:9,2) = H_3bdy_part(5:9,2) + H_temp(1:5) * G_IJK(3,2)

                  ! Calculate repeating part the K-S matrix elements:
                  H_temp(1) = Mjs(i,k,3) * Mjs(j,k,5) ! py-s * dxy-s
                  H_temp(2) = Mjs(i,k,3) * Mjs(j,k,6) ! py-s * dxz-s
                  H_temp(3) = Mjs(i,k,3) * Mjs(j,k,7) ! py-s * dxz-s
                  H_temp(4) = Mjs(i,k,3) * Mjs(j,k,8) ! py-s * (dx2-y2)-s
                  H_temp(5) = Mjs(i,k,3) * Mjs(j,k,9) ! py-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(3,5:9) = H_3bdy_part(3,5:9) + H_temp(1:5) * G_IJK(2,3)
                  H_3bdy_part(5:9,3) = H_3bdy_part(5:9,3) + H_temp(1:5) * G_IJK(3,2)

                  ! Calculate repeating part the K-S matrix elements:
                  H_temp(1) = Mjs(i,k,4) * Mjs(j,k,5) ! dxy-s * pz-s
                  H_temp(2) = Mjs(i,k,4) * Mjs(j,k,6) ! dxz-s * pz-s
                  H_temp(3) = Mjs(i,k,4) * Mjs(j,k,7) ! dyz-s * pz-s
                  H_temp(4) = Mjs(i,k,4) * Mjs(j,k,8) ! (dx2-y2)-s * pz-s
                  H_temp(5) = Mjs(i,k,4) * Mjs(j,k,9) ! (d3z2-r2)-s * pz-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(4,5:9) = H_3bdy_part(4,5:9) + H_temp(1:5) * G_IJK(2,3)
                  H_3bdy_part(5:9,4) = H_3bdy_part(5:9,4) + H_temp(1:5) * G_IJK(3,2)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,5) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dxy-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dxy-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dxy-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dxy-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dxy-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(5,5:9) = H_3bdy_part(5,5:9) + H_temp(1:5)
                  H_3bdy_part(5:9,5) = H_3bdy_part(5:9,5) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,6) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dxz-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dxz-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dxz-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dxz-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dxz-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(6,5:9) = H_3bdy_part(6,5:9) + H_temp(1:5)
                  H_3bdy_part(5:9,6) = H_3bdy_part(5:9,6) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,7) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dyz-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dyz-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dyz-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dyz-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dyz-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(7,5:9) = H_3bdy_part(7,5:9) + H_temp(1:5)
                  H_3bdy_part(5:9,7) = H_3bdy_part(5:9,7) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,8) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! (dx2-y2)-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! (dx2-y2)-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! (dx2-y2)-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! (dx2-y2)-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! (dx2-y2)-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(8,5:9) = H_3bdy_part(8,5:9) + H_temp(1:5)
                  H_3bdy_part(5:9,8) = H_3bdy_part(5:9,8) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,9) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! (d3z2-r2)-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! (d3z2-r2)-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! (d3z2-r2)-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! (d3z2-r2)-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! (d3z2-r2)-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(9,5:9) = H_3bdy_part(9,5:9) + H_temp(1:5)
                  H_3bdy_part(5:9,9) = H_3bdy_part(5:9,9) + H_temp(1:5)

               endif ! (basis_ind > 1)
            endif ! (basis_ind > 0)
         endif ! (k /= j)
      enddo AT3
!!!$omp end do
!!!$omp END PARALLEL
   endif ! (TB(KOA1,KOA1)%include_3body)

   5000 continue

   ! Add up 2-body and 3-body parts:
   ts = ts + H_3bdy_part

   nullify(m, KOA1, KOA2, KOA3)
end subroutine Hopping_3TB




subroutine Onsite_3TB(basis_ind, i, Scell, TB, M_Lag_exp, Mjs, Hij)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   integer, intent(in) :: i         ! indix of the atom
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   real(8), dimension(:,:), intent(out) :: Hij	! overlap integerals [eV]
   !---------------------
   integer :: i1, Bsiz, atom_2, sh1, sh2, j, k, at_ind, atom_3
   real(8) :: r, term_s, term_p, term_d
   real(8) :: matr_spd(3,3), H_cf_temp(9)
   real(8), dimension(:,:), allocatable :: E_onsite, H_avg, H_cf, H_3bdy
   integer, pointer :: m, KOA1, KOA2, KOA3

   ! Kind of atom:
   KOA1 => Scell%MDatoms(i)%KOA
   ! Number of the nearest neighbors
   m => Scell%Near_neighbor_size(i)

   ! Size of the basis:
   Bsiz = size(Hij,1)
   ! To start with:
   allocate(E_onsite(Bsiz,Bsiz), source=0.0d0)
   allocate(H_avg(Bsiz,Bsiz), source=0.0d0)
   allocate(H_cf(Bsiz,Bsiz), source=0.0d0)
   allocate(H_3bdy(Bsiz,Bsiz), source=0.0d0)

   ! The onsite energies are constructed out of 4 terms [1]:
   ! 1) onsite eigenvalues
   ! 2) average contribution of atoms around
   ! 3) crystall field contribution
   ! 4) 3-body contributions
   ! Let's calculate them all:


   !-----------------
   ! 1) Onsite energies of spin-unpolirized orbital values:
   E_onsite(1,1) = TB(KOA1,KOA1)%Es ! s orbital
   if (basis_ind > 0) then ! p3 orbitals:
      do i1 = 2, 4
         E_onsite(i1,i1) = TB(KOA1,KOA1)%Ep
      enddo

      if (basis_ind > 1) then ! d5 orbitals:
         do i1 = 5, 9
            E_onsite(i1,i1) = TB(KOA1,KOA1)%Ed
         enddo
      endif
   endif


!    print*, 'Onsite_3TB test 0', m
! goto 5001   ! For testing, exclude environmental contributions: seems to work...

!goto 5002

   !-----------------
   ! 2) Average term:
!!!$omp parallel private(atom_2, j, i1, term_s, term_p, term_d) ! in is called from parallelized region, don't parallelize it again
!!!$omp do reduction( + : H_avg)
   do atom_2 = 1, m ! do only for atoms close to that one
      j = Scell%Near_neighbor_list(i, atom_2) ! this is the list of such close atoms
      KOA2 => Scell%MDatoms(j)%KOA

      term_s = SUM( TB(KOA1,KOA2)%Hhavg(1,1:4) * M_Lag_exp(i,j,1:4) )  ! s orbitals
      H_avg(1,1) = H_avg(1,1) + term_s

      if (basis_ind > 0) then ! p3 orbitals:
         term_p = SUM( TB(KOA1,KOA2)%Hhavg(2,1:4) * M_Lag_exp(i,j,1:4) )  ! p3 orbitals
         do i1 = 2, 4
            H_avg(i1,i1) = H_avg(i1,i1) + term_p
         enddo

         if (basis_ind > 1) then ! d5 orbitals:
            term_d = SUM( TB(KOA1,KOA2)%Hhavg(3,1:4) * M_Lag_exp(i,j,1:4) )  ! d5 orbitals
            do i1 = 5, 9
               H_avg(i1,i1) = H_avg(i1,i1) + term_d
            enddo
         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)
   enddo
!!!$omp enddo
!!!$omp end parallel

5002 continue

! print*, 'Onsite_3TB test 1', m
! goto 5001  ! For testing, exclude crystal field and 3-body terms

!  goto 5003

   !-----------------
   ! 3) Crystal field:
!!!$omp parallel private(atom_2, j, sh1, sh2, i1, matr_spd, H_cf_temp)
!!!$omp do
   do atom_2 = 1, m ! do only for atoms close to that one
      j = Scell%Near_neighbor_list(i, atom_2) ! this is the list of such close atoms
      KOA2 => Scell%MDatoms(j)%KOA

      ! Radial parts:
      matr_spd(1,1) = SUM( TB(KOA1,KOA2)%Hhcf(1,1,1:4) * M_Lag_exp(i,j,1:4) )  ! s-s orbitals
      if (basis_ind > 0) then ! p3 orbitals:
         matr_spd(1,2) = SUM( TB(KOA1,KOA2)%Hhcf(1,2,1:4) * M_Lag_exp(i,j,1:4) )  ! s-p orbitals
         matr_spd(2,1) = SUM( TB(KOA1,KOA2)%Hhcf(2,1,1:4) * M_Lag_exp(i,j,1:4) )  ! p-s orbitals
         matr_spd(2,2) = SUM( TB(KOA1,KOA2)%Hhcf(2,2,1:4) * M_Lag_exp(i,j,1:4) )  ! p-p orbitals
         if (basis_ind > 1) then ! d5 orbitals:
            matr_spd(1,3) = SUM( TB(KOA1,KOA2)%Hhcf(1,3,1:4) * M_Lag_exp(i,j,1:4) )  ! s-d orbitals
            matr_spd(3,1) = SUM( TB(KOA1,KOA2)%Hhcf(3,1,1:4) * M_Lag_exp(i,j,1:4) )  ! d-s orbitals
            matr_spd(2,3) = SUM( TB(KOA1,KOA2)%Hhcf(2,3,1:4) * M_Lag_exp(i,j,1:4) )  ! p-d orbitals
            matr_spd(3,2) = SUM( TB(KOA1,KOA2)%Hhcf(3,2,1:4) * M_Lag_exp(i,j,1:4) )  ! d-p orbitals
            matr_spd(3,3) = SUM( TB(KOA1,KOA2)%Hhcf(3,3,1:4) * M_Lag_exp(i,j,1:4) )  ! d-d orbitals
         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)

      ! Include angular parts:
      H_cf(1,1) = H_cf(1,1) + matr_spd(1,1)  ! s-s * s-s

!  goto 5004

      if (basis_ind > 0) then ! p3 orbitals:

         H_cf(1,2) = H_cf(1,2) + matr_spd(2,1) * Mjs(i,j,2)   ! s-s * px-s
         H_cf(1,3) = H_cf(1,3) + matr_spd(2,1) * Mjs(i,j,3)   ! s-s * py-s
         H_cf(1,4) = H_cf(1,4) + matr_spd(2,1) * Mjs(i,j,4)   ! s-s * pz-s

         H_cf(2,1) = H_cf(2,1) + matr_spd(1,2) * Mjs(i,j,2)   ! px-s * s-s
         H_cf(3,1) = H_cf(3,1) + matr_spd(1,2) * Mjs(i,j,3)   ! py-s * s-s
         H_cf(4,1) = H_cf(4,1) + matr_spd(1,2) * Mjs(i,j,4)   ! pz-s * s-s

         ! Calculate repeating part the K-S matrix elements:
         H_cf_temp = 0.0d0 ! reinitialize
         H_cf_temp(1) = matr_spd(2,2) * Mjs(i,j,2)
         H_cf_temp(2) = H_cf_temp(1) * Mjs(i,j,2) ! px-s * px-s
         H_cf_temp(3) = H_cf_temp(1) * Mjs(i,j,3) ! px-s * py-s
         H_cf_temp(4) = H_cf_temp(1) * Mjs(i,j,4) ! px-s * pz-s
         H_cf_temp(5) = matr_spd(2,2) * Mjs(i,j,3)
         H_cf_temp(6) = H_cf_temp(5) * Mjs(i,j,4)
         ! Add it into the Hamiltonian part:


         ! Diagonal part:
         H_cf(2,2) = H_cf(2,2) + H_cf_temp(2) / 3.0d0                               ! px-s * px-s
         H_cf(3,3) = H_cf(3,3) + H_cf_temp(5) / 3.0d0 * Mjs(i,j,3)                  ! py-s * py-s
         H_cf(4,4) = H_cf(4,4) + matr_spd(2,2) / 3.0d0 * Mjs(i,j,4) * Mjs(i,j,4)    ! pz-s * pz-s


         ! Off-diagonal part:
         H_cf(2,3) = H_cf(2,3) + H_cf_temp(3)  ! px-s * py-s
         H_cf(2,4) = H_cf(2,4) + H_cf_temp(4)  ! px-s * pz-s
         H_cf(3,2) = H_cf(3,2) + H_cf_temp(3)  ! py-s * px-s
         H_cf(3,4) = H_cf(3,4) + H_cf_temp(6)  ! py-s * pz-s
         H_cf(4,2) = H_cf(4,2) + H_cf_temp(4)  ! pz-s * px-s
         H_cf(4,3) = H_cf(4,3) + H_cf_temp(6)  ! pz-s * py-s

         if (basis_ind > 1) then ! d5 orbitals:
            H_cf(1,5) = H_cf(1,5) + Mjs(i,j,5) * matr_spd(1,3) ! s-s * dxy-s
            H_cf(1,6) = H_cf(1,6) + Mjs(i,j,6) * matr_spd(1,3) ! s-s * dxz-s
            H_cf(1,7) = H_cf(1,7) + Mjs(i,j,7) * matr_spd(1,3) ! s-s * dxz-s
            H_cf(1,8) = H_cf(1,8) + Mjs(i,j,8) * matr_spd(1,3) ! s-s * (dx2-y2)-s
            H_cf(1,9) = H_cf(1,9) + Mjs(i,j,9) * matr_spd(1,3) ! s-s * (d3z2-r2)-s

            H_cf(5,1) = H_cf(5,1) + Mjs(i,j,5) * matr_spd(3,1) ! dxy-s * s-s
            H_cf(6,1) = H_cf(6,1) + Mjs(i,j,6) * matr_spd(3,1) ! dxz-s * s-s
            H_cf(7,1) = H_cf(7,1) + Mjs(i,j,7) * matr_spd(3,1) ! dxz-s * s-s
            H_cf(8,1) = H_cf(8,1) + Mjs(i,j,8) * matr_spd(3,1) ! (dx2-y2)-s * s-s
            H_cf(9,1) = H_cf(9,1) + Mjs(i,j,9) * matr_spd(3,1) ! (d3z2-r2)-s * s-s

            ! Calculate repeating part the K-S matrix elements:
            H_cf_temp(1) = Mjs(i,j,2) * Mjs(i,j,5) ! px-s * dxy-s
            H_cf_temp(2) = Mjs(i,j,2) * Mjs(i,j,6) ! px-s * dxz-s
            H_cf_temp(3) = Mjs(i,j,2) * Mjs(i,j,7) ! px-s * dxz-s
            H_cf_temp(4) = Mjs(i,j,2) * Mjs(i,j,8) ! px-s * (dx2-y2)-s
            H_cf_temp(5) = Mjs(i,j,2) * Mjs(i,j,9) ! px-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(2,5:9) = H_cf(2,5:9) + H_cf_temp(1:5) * matr_spd(2,3)
            H_cf(5:9,2) = H_cf(5:9,2) + H_cf_temp(1:5) * matr_spd(3,2)

            ! Calculate repeating part the K-S matrix elements:
            H_cf_temp(1) = Mjs(i,j,3) * Mjs(i,j,5) ! py-s * dxy-s
            H_cf_temp(2) = Mjs(i,j,3) * Mjs(i,j,6) ! py-s * dxz-s
            H_cf_temp(3) = Mjs(i,j,3) * Mjs(i,j,7) ! py-s * dxz-s
            H_cf_temp(4) = Mjs(i,j,3) * Mjs(i,j,8) ! py-s * (dx2-y2)-s
            H_cf_temp(5) = Mjs(i,j,3) * Mjs(i,j,9) ! py-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(3,5:9) = H_cf(3,5:9) + H_cf_temp(1:5) * matr_spd(2,3)
            H_cf(5:9,3) = H_cf(5:9,3) + H_cf_temp(1:5) * matr_spd(3,2)

            ! Calculate repeating part the K-S matrix elements:
            H_cf_temp(1) = Mjs(i,j,4) * Mjs(i,j,5) ! dxy-s * pz-s
            H_cf_temp(2) = Mjs(i,j,4) * Mjs(i,j,6) ! dxz-s * pz-s
            H_cf_temp(3) = Mjs(i,j,4) * Mjs(i,j,7) ! dyz-s * pz-s
            H_cf_temp(4) = Mjs(i,j,4) * Mjs(i,j,8) ! (dx2-y2)-s * pz-s
            H_cf_temp(5) = Mjs(i,j,4) * Mjs(i,j,9) ! (d3z2-r2)-s * pz-s
            ! Add it into the Hamiltonian part:
            H_cf(4,5:9) = H_cf(4,5:9) + H_cf_temp(1:5) * matr_spd(2,3)
            H_cf(5:9,4) = H_cf(5:9,4) + H_cf_temp(1:5) * matr_spd(3,2)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,5) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! dxy-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! dxy-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! dxy-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! dxy-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! dxy-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(5,5:9) = H_cf(5,5:9) + H_cf_temp(1:5)
            H_cf(5:9,5) = H_cf(5:9,5) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,6) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! dxz-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! dxz-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! dxz-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! dxz-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! dxz-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(6,5:9) = H_cf(6,5:9) + H_cf_temp(1:5)
            H_cf(5:9,6) = H_cf(5:9,6) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,7) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! dyz-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! dyz-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! dyz-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! dyz-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! dyz-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(7,5:9) = H_cf(7,5:9) + H_cf_temp(1:5)
            H_cf(5:9,7) = H_cf(5:9,7) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,8) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! (dx2-y2)-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! (dx2-y2)-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! (dx2-y2)-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! (dx2-y2)-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! (dx2-y2)-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(8,5:9) = H_cf(8,5:9) + H_cf_temp(1:5)
            H_cf(5:9,8) = H_cf(5:9,8) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,9) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! (d3z2-r2)-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! (d3z2-r2)-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! (d3z2-r2)-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! (d3z2-r2)-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! (d3z2-r2)-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(9,5:9) = H_cf(9,5:9) + H_cf_temp(1:5)
            H_cf(5:9,9) = H_cf(5:9,9) + H_cf_temp(1:5)

         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)

5004 continue

   enddo ! atom_2 = 0, m
!!!$omp enddo
!!!$omp end parallel

5003 continue



! print*, 'Onsite_3TB test 2', m
!goto 5001   ! For testing, exclude 3-body terms

   !-----------------
   ! 4) 3-body contributions:
   if (TB(KOA1,KOA1)%include_3body) then  ! only if user defined it to include
      atom_2 = 0  ! to restart
!!!$omp PARALLEL
!!!$omp do private(atom_2, j, KOA2, atom_3, k, KOA3, at_ind, term_s, sh1)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         KOA2 => Scell%MDatoms(j)%KOA   ! index of the second atom

!          print*, 'Onsite_3TB test 2.5', j

         ! To start summing up
         AT3:do atom_3 = 1,m ! do only for atoms close to that one
            k = Scell%Near_neighbor_list(i,atom_3) ! this is the list of such close atoms
            ! Make sure the third atoms is not the second atom:
            if (k /= j) then
               KOA3 => Scell%MDatoms(k)%KOA   ! kind of atom #3

               ! Find the combination-of-atoms index:
               at_ind = find_3bdy_ind(KOA1, KOA2, KOA3)  ! module "Dealing_with_3TB"

               ! Get the 3-body elements for each shells combination:

               term_s = TB(KOA1,KOA2)%Hh3bdy(at_ind, 1) * M_Lag_exp(i,k,1) * M_Lag_exp(j,k,1) * M_Lag_exp(i,j,1) + &
                        TB(KOA1,KOA2)%Hh3bdy(at_ind, 2) * M_Lag_exp(i,k,2) * M_Lag_exp(j,k,1) * M_Lag_exp(i,j,1) + &
                        TB(KOA1,KOA2)%Hh3bdy(at_ind, 3) * M_Lag_exp(i,k,1) * M_Lag_exp(j,k,2) * M_Lag_exp(i,j,1) + &
                        TB(KOA1,KOA2)%Hh3bdy(at_ind, 4) * M_Lag_exp(i,k,1) * M_Lag_exp(j,k,1) * M_Lag_exp(i,j,2) ! Eq.(22) in [1]

               do sh1 = 1, Bsiz  ! for all shells
                  H_3bdy(sh1,sh1) = H_3bdy(sh1,sh1) + term_s   ! no orbital resolution
               enddo ! sh1
            endif ! (k /= i)
         enddo AT3
      enddo AT2
!!!$omp end do
!!!$omp END PARALLEL
   endif ! (TB(KOA1,KOA1)%include_3body)


! print*, 'Onsite_3TB test 3', m

5001 continue

   ! Collect all the terms into Hamiltonian:
   Hij = E_onsite + H_avg + H_cf + H_3bdy

   deallocate(E_onsite, H_avg, H_cf, H_3bdy)
   nullify(m, KOA1, KOA2, KOA3)
end subroutine Onsite_3TB




!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Multiplication of the two Koster-Slater angular functions (Mis*Mjs):

subroutine get_Mjs_factors(basis_ind, Scell, M_lmn, Mjs)
   integer, intent(in) :: basis_ind ! size of the basis set
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Mjs ! K-S multiplication with dummy arguments for radial function

   integer :: atom_2, j, i, Nsiz, nat
   integer, pointer :: m

   nat = Scell%Na ! number of atoms
   ! Find the number of orbitals:
   Nsiz = identify_DFTB_orbitals_per_atom(basis_ind)  ! module "TB_DFTB"

   if (.not.allocated(Mjs)) allocate( Mjs(nat,nat,Nsiz) )
   Mjs = 0.0d0 ! to start with

!$omp parallel
!$omp do private(i, j, m, atom_2)
   do i = 1, nat	! all atoms
      m => Scell%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms

         Mjs(i,j,1) = 1.0d0  ! s-s

         if (basis_ind > 0) then ! p3 orbitals:
            Mjs(i,j,2) = -M_lmn(1,i,j)   ! px-s
            Mjs(i,j,3) = -M_lmn(2,i,j)   ! py-s
            Mjs(i,j,4) = -M_lmn(3,i,j)   ! pz-s

            if (basis_ind > 1) then ! d5 orbitals (functions from module "TB_Koster_Slater"):
               Mjs(i,j,5) = t_s_dab(M_lmn(1,i,j),M_lmn(2,i,j),1.0d0)        ! dxy-s
               Mjs(i,j,6) = t_s_dab(M_lmn(1,i,j),M_lmn(3,i,j),1.0d0)        ! dxz-s
               Mjs(i,j,7) = t_s_dab(M_lmn(2,i,j),M_lmn(3,i,j),1.0d0)        ! dxz-s
               Mjs(i,j,8) = t_s_dx2_y2(M_lmn(1,i,j),M_lmn(2,i,j),1.0d0)     ! (dx2-y2)-s
               Mjs(i,j,9) = t_s_dz2_r2(M_lmn(1,i,j),M_lmn(2,i,j),M_lmn(3,i,j),1.0d0)    ! (d3z2-r2)-s
            endif
         endif
      enddo
   enddo
!$omp enddo
!$omp end parallel
nullify(m)
end subroutine get_Mjs_factors



!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Tools for radial function calculation within 3TB model:


subroutine get_Laguerres(r_given, Laguer)
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(6), intent(out) :: Laguer ! polynomials up to 6
   real(8) :: d, ad
   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au   ! [A] -> [Bohr]
   !d = r_given    ! [A] -> [Bohr] ! Testing
   ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]
   !ad = d ! for some reason, it is multiplied with the scaling parameter [2] ! Testing
   !ad = m_a*g_A2au*d ! for some reason, it is multiplied with the scaling parameter [2] ! Testing

   call Laguerre_up_to_6(ad, Laguer)   ! module "Algebra_tools"
end subroutine get_Laguerres


pure function Laguerre_exponent(r_given) result(exp_L)
   real(8) exp_L
   real(8), intent(in) :: r_given   ! distance [A]
   real(8) :: d, ad

   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au   ! [A] -> [Bohr]
   !d = r_given       ! [A] -> [Bohr] ! Testing
   ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]
   !ad = m_a*g_A2au*d ! for some reason, it is multiplied with the scaling parameter [2] ! Testing

   ! Get the exponential term:
   exp_L = exp(-0.5d0 * ad)
   !exp_L = exp(-ad)   ! Testing
end function Laguerre_exponent




subroutine get_d_Laguerres(r_given, Laguer)
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(6), intent(out) :: Laguer ! derivatives of the polynomials up to 6
   real(8) :: d, ad
   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au   ! [A] -> [Bohr]
   !ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]
   ad = m_a*g_A2au*d ! for some reason, it is multiplied with the scaling parameter [2] ! Testing

   call d_Laguerre_up_to_6(ad, Laguer)   ! module "Algebra_tools"
end subroutine get_d_Laguerres



function radial_function_3TB(r_given, fx_coefs, Laguer, sh_ind, num_lag, no_exp) result(f_out)
   real(8) f_out    ! value found in the array and interpolated
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(:,:), intent(in) :: fx_coefs    ! coefficients before Leguerres in 2-body Hamiltonian
   real(8), dimension(:), intent(in) :: Laguer   ! precalculated Laguerre polinomials
   integer, intent(in) :: sh_ind, num_lag    ! which shell, how many polynomials
   logical, intent(in), optional :: no_exp   ! indicate if Laguer = Laguerre * exp(-a*r_IJ), or not
   !---------------------
   real(8) :: d, ad, exp_ad
   integer :: i
   logical :: there_is_exp

   there_is_exp = .true.  ! by default, there is exp. provided inside of Laguer
   if (present(no_exp)) then
      if (no_exp) then
         there_is_exp = .false.   ! if it is not provided, recalculate it
      endif
   endif

   ! Construct the function:
   f_out = 0.0d0  ! to start with
   do i = 1, num_lag ! for all polynomials required
      f_out = f_out + fx_coefs(sh_ind, i) * Laguer(i)
   enddo

   ! In case the Laguerre polynomials did not contain exponent, calculate it:
   if (.not.there_is_exp) then
      ! Get the exponential term:
      exp_ad = Laguerre_exponent(r_given) ! below

      ! Include it into the solution:
      f_out = exp_ad * f_out
   endif
end function radial_function_3TB



function d_radial_function_3TB(r_given, M_Vij, fx_coefs, d_Laguer, sh_ind, num_lag) result(f_out)
   real(8) f_out    ! value found in the array and interpolated
   real(8), intent(in) :: r_given   ! [A] given radial coordinate
   real(8), intent(in) :: M_Vij     ! radial function
   real(8), dimension(:,:), intent(in) :: fx_coefs    ! coefficients before Leguerres in 2-body Hamiltonian
   real(8), dimension(:), intent(in) :: d_Laguer   ! precalculated derivatives of the Laguerre polinomials
   integer, intent(in) :: sh_ind, num_lag    ! which shell, how many polynomials
   !---------------------
   real(8) :: d, ad, exp_ad
   integer :: i

   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au   ! [A] -> [Bohr]
   ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]

   ! Get the exponential term:
   exp_ad = exp(-0.5d0 * ad)

   ! Construct the function:
   f_out = 0.0d0  ! to start with
   ! Skip i=1, since dL_0/dr = 0
   do i = 2, num_lag ! for all polynomials required
      f_out = f_out + fx_coefs(sh_ind, i) * d_Laguer(i)
   enddo
   f_out = exp_ad * f_out - (m_a*0.5d0/g_a0)*M_Vij
end function d_radial_function_3TB





!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Repulsive part of 3TB is absent, so just set zeros everywhere:

subroutine get_Erep_s_3TB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_3TB), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   a = 0.0d0	! There is no repulsive term in 3TB parameterization
end subroutine get_Erep_s_3TB


subroutine dErdr_s_3TB(TB_Repuls, Scell, NSC) ! derivatives of the repulsive energy by s
   type(TB_Rep_3TB), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer ian, n
   n = Scell(NSC)%Na ! number of atoms
   !$omp PARALLEL DO private(ian)
   do ian = 1, n  ! Forces for all atoms
     Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! no added repulsion
   enddo ! ian
   !$OMP END PARALLEL DO
END subroutine dErdr_s_3TB


subroutine dErdr_Pressure_s_3TB(TB_Repuls, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_3TB), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      Scell(NSC)%SCforce%rep = 0.0d0   ! no added repulsion
   endif
end subroutine dErdr_Pressure_s_3TB


END MODULE TB_3TB
