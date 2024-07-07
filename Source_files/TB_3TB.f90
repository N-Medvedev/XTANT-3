! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT-3
! available at: https://zenodo.org/doi/10.5281/zenodo.8392569
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
! This module contains subroutines to deal with 3TB hamiltonian: https://github.com/usnistgov/ThreeBodyTB.jl
! [1] https://arxiv.org/pdf/2112.11585.pdf
! [2] https://github.com/usnistgov/ThreeBodyTB.jl/blob/master/src/CalcTB_laguerre.jl
! So far, without the 3-body-terms; 2-body-terms are included, but problematic (not recommended to use)

MODULE TB_3TB

use Universal_constants
use Objects
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity, &
                        Laguerre_up_to_6, d_Laguerre_up_to_6, check_symmetry
use Atomic_tools, only : Reciproc_rel_to_abs, shortest_distance
use Little_subroutines, only : linear_interpolation, Fermi_function, d_Fermi_function, Find_in_array_monoton
use Electron_tools, only : find_band_gap
use TB_Koster_Slater
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, &
                  test_orthogonalization_c, Loewdin_Orthogonalization, Loewdin_Orthogonalization_c
use TB_DFTB, only : identify_DFTB_basis_size, identify_DFTB_orbitals_per_atom, Get_overlap_S_matrix_DFTB
use Dealing_with_3TB, only: find_3bdy_ind

#ifdef MPI_USED
   use MPI_subroutines, only : do_MPI_Allreduce
#endif


implicit none
PRIVATE

real(8), parameter :: m_a = 2.0d0   ! exponential decay parameter according to [1]

public :: get_Erep_s_3TB, dErdr_s_3TB, dErdr_Pressure_s_3TB, Attract_TB_Forces_Press_3TB, &
          Construct_Vij_3TB, construct_TB_H_3TB, get_Mjs_factors, get_dHij_drij_3TB

 contains



subroutine Construct_Vij_3TB(numpar, TB, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_Lag_exp, M_d_Lag_exp)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Vij	! matrix of hoppings for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dVij	! matrix of derivatives of hoppings for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_SVij	! S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dSVij  ! derivatives of S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Lag_exp, M_d_Lag_exp ! matrix of laguerre * exp(-a*r_ij)

   !----------------------------------
   real(8) :: x, y, z, r, fcut, d_fcut, Fermi, dFermi, Laguer(6), d_Laguer(6), rIK, rJK, exp_ad, d_exp_ad
   integer :: i, j, atom_2, ki, N_bs, ihop, atom_3, k, at_ind, sh1, sh2
   real(8), pointer :: rm
   integer, pointer :: nat, m, KOA1, KOA2, KOA3, mm
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   nat => Scell(NSC)%Na	! number of atoms in the supercell
   ! number of hopping integrals for this basis set in 3TB:
   N_bs = identify_DFTB_basis_size(numpar%basis_size_ind)  ! module "TB_DFTB"

   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,N_bs))	    ! each pair of atoms, all  V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,N_bs))	! each pair of atoms, all  dV functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,N_bs))	! each pair of atoms, all S functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,N_bs))	! each pair of atoms, all dS functions
   if (.not.allocated(M_Lag_exp)) allocate(M_Lag_exp(nat,nat,6))  ! each pair of atoms, 6 polynomials
   if (.not.allocated(M_d_Lag_exp)) allocate(M_d_Lag_exp(nat,nat,6))  ! each pair of atoms, 6 polynomials

   M_Vij = 0.0d0
   M_dVij = 0.0d0
   M_SVij = 0.0d0
   M_dSVij = 0.0d0


   ! 2-body interaction terms:

   ! Construct matrix of all the radial functions for each pair of atoms:
#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = nat
   ! Do the cycle (parallel) calculations:
   AT10:do j = Nstart, Nend, N_incr  ! each process does its own part
   !AT10:do j = 1,nat	! atom #1
      m => Scell(NSC)%Near_neighbor_size(j)
      AT20:do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! atom #2
         ! Kinds of atoms (elements indices):
         KOA1 => Scell(NSC)%MDatoms(j)%KOA   ! atom #1
         KOA2 => Scell(NSC)%MDatoms(i)%KOA   ! atom #2

         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! Get Laguerre polynomials and its derivatives to be reused below:
         call get_Laguerres(r, TB(KOA1,KOA2)%rc, Laguer) ! below
         call get_d_Laguerres(r, TB(KOA1,KOA2)%rc, d_Laguer) ! below

         ! Also get the cut-off functions, to include directly into radial dependence:
         Fermi = Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         dFermi = d_Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"

         ! Get the exponential term:
         exp_ad = Laguerre_exponent(r, TB(KOA1,KOA2)%rc) ! below
         d_exp_ad = d_Laguerre_exponent(r, TB(KOA1,KOA2)%rc) ! below

         ! Include the exponenet and the cut-off in the solution:
         ! 1) First, construct the derivative via chain rule:
         d_Laguer(:) = d_Laguer(:)*exp_ad*Fermi + Laguer(:)*(d_exp_ad*Fermi + exp_ad*dFermi)
         ! 2) Now, update the Laguerres themselves:
         Laguer(:) = Laguer(:)*exp_ad*Fermi

         ! Save [Laguerre * exp(-a*r_ij) * cutoff]:
         M_Lag_exp(j,i,:) = Laguer(:)
         M_d_Lag_exp(j,i,:) = d_Laguer(:)

         ! All radial functions for Hamiltonian (functions below):
         M_Vij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, 1, 5)   ! (s s sigma)
         M_SVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, 1, 6)  ! (s s sigma)
         select case (numpar%basis_size_ind)
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

         ! All derivatives of the radial functions and the overlap matrix:
         M_dVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 1, 5)   ! (s s sigma)
         M_dSVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 1, 6) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
            M_dVij(j,i,2) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 2, 5)   ! (s p sigma)
            M_dSVij(j,i,2) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 2, 6) ! (s p sigma)
            M_dVij(j,i,3) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 4, 5)   ! (p p sigma)
            M_dSVij(j,i,3) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 4, 6) ! (p p sigma)
            M_dVij(j,i,4) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 5, 5)   ! (p p pi)
            M_dSVij(j,i,4) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 5, 6) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_dVij(j,i,ihop) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, ihop, 5)
               M_dSVij(j,i,ihop) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, ihop, 6)
            enddo
         endselect

      enddo AT20
   enddo AT10
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in Construct_Vij_3TB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_Lag_exp', M_Lag_exp) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_d_Lag_exp', M_d_Lag_exp) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_Vij', M_Vij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_SVij', M_SVij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_dVij', M_dVij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_dSVij', M_dSVij) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
!$omp PARALLEL
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, Laguer, d_Laguer, exp_ad, d_exp_ad, ihop, Fermi, dFermi)
   AT1:do j = 1,nat	! atom #1

      m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! atom #2
         ! Kinds of atoms (elements indices):
         ! [HB 0]
         KOA1 => Scell(NSC)%MDatoms(j)%KOA   ! atom #1
         KOA2 => Scell(NSC)%MDatoms(i)%KOA   ! atom #2
         ! [HB 1]
!          KOA2 => Scell(NSC)%MDatoms(j)%KOA   ! atom #1
!          KOA1 => Scell(NSC)%MDatoms(i)%KOA   ! atom #2

         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! Get Laguerre polynomials and its derivatives to be reused below:
         call get_Laguerres(r, TB(KOA1,KOA2)%rc, Laguer) ! below
         call get_d_Laguerres(r, TB(KOA1,KOA2)%rc, d_Laguer) ! below

         ! Also get the cut-off functions, to include directly into radial dependence:
         Fermi = Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         dFermi = d_Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"

         ! Get the exponential term:
         exp_ad = Laguerre_exponent(r, TB(KOA1,KOA2)%rc) ! below
         d_exp_ad = d_Laguerre_exponent(r, TB(KOA1,KOA2)%rc) ! below

         ! Include the exponenet and the cut-off in the solution:
         ! 1) First, construct the derivative via chain rule:
         d_Laguer(:) = d_Laguer(:)*exp_ad*Fermi + Laguer(:)*(d_exp_ad*Fermi + exp_ad*dFermi)
         ! 2) Now, update the Laguerres themselves:
         Laguer(:) = Laguer(:)*exp_ad*Fermi

         ! Save [Laguerre * exp(-a*r_ij) * cutoff]:
         M_Lag_exp(j,i,:) = Laguer(:)
         M_d_Lag_exp(j,i,:) = d_Laguer(:)

         ! All radial functions for Hamiltonian (functions below):
         M_Vij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, Laguer, 1, 5)   ! (s s sigma)
         M_SVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, Laguer, 1, 6)  ! (s s sigma)
         select case (numpar%basis_size_ind)
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

         ! All derivatives of the radial functions and the overlap matrix:
         M_dVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 1, 5)   ! (s s sigma)
         M_dSVij(j,i,1) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 1, 6) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
            M_dVij(j,i,2) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 2, 5)   ! (s p sigma)
            M_dSVij(j,i,2) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 2, 6) ! (s p sigma)
            M_dVij(j,i,3) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 4, 5)   ! (p p sigma)
            M_dSVij(j,i,3) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 4, 6) ! (p p sigma)
            M_dVij(j,i,4) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, 5, 5)   ! (p p pi)
            M_dSVij(j,i,4) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, 5, 6) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_dVij(j,i,ihop) = radial_function_3TB(r, TB(KOA1,KOA2)%Vrfx, d_Laguer, ihop, 5)
               M_dSVij(j,i,ihop) = radial_function_3TB(r, TB(KOA1,KOA2)%Srfx, d_Laguer, ihop, 6)
            enddo
         endselect

      enddo AT2
   enddo AT1
!$omp end do
!$omp END PARALLEL
#endif

   nullify(rm, nat, m, KOA1, KOA2, KOA3)	! clean up at the end
end subroutine Construct_Vij_3TB




! Tight Binding Hamiltonian within 3TB parametrization:
subroutine construct_TB_H_3TB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Scell, NSC, Err, scc, H_scc_0, H_scc_1)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! material parameters
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   type(Super_cell), dimension(:), intent(inout) :: Scell		! supercell with all the atoms as one object
   integer, intent(in) :: NSC		! number of supercell
   type(Error_handling), intent(inout) :: Err	! error save
   integer, intent(in), optional :: scc   ! flag to do self-consistent charge calculations
   real(8), dimension(:,:), intent(inout), optional :: H_scc_0, H_scc_1 ! zero and first order contributions to Hamiltonian
   character(200) :: Error_descript
   Error_descript = ''

   Scell(NSC)%Ha0 = Scell(NSC)%Ha	! save Hamiltonial from previous time-step
   Scell(NSC)%Ei0 = Scell(NSC)%Ei	! save energy levels for the next timestep
   Scell(NSC)%H_non0 = Scell(NSC)%H_non	! save non-diagonalized Hamiltonian from last time-step

    ! Construct TB Hamiltonian (with DFTB parameters),  orthogonalize it,  and then solve generalized eigenvalue problem:
   if (present(scc) .and. present(H_scc_0) .and. present(H_scc_1)) then
      call Hamil_tot_3TB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Err, scc, H_scc_0, H_scc_1) ! below
   else
      call Hamil_tot_3TB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Err)  ! below
   endif

!    ! Test (comment out for release):
!    call test_nonorthogonal_solution(Scell(NSC)) ! module "TB_NRL"

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_3TB



subroutine Hamil_tot_3TB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_Lag_exp, M_lmn, Mjs, Err, scc, H_scc_0, H_scc_1)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij)
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   type(Error_handling), intent(inout) :: Err	! error save
   integer, intent(in), optional :: scc   ! flag to do self-consistent charge calculations
   real(8), dimension(:,:), intent(inout), optional :: H_scc_0, H_scc_1 ! zero and first order contributions to Hamiltonian
   !-------------------------------------------
   real(8), dimension(:,:), allocatable :: Hij	 ! Hamiltonian
   real(8), dimension(:,:), allocatable :: Sij  ! Overlap
   real(8), dimension(size(Scell(NSC)%Ha,1)) :: Evec, EvecS
   real(8), dimension(:,:), allocatable :: Hij1, Sij1
   integer :: nat, Nsiz, n_orb, do_scc
   integer, target :: j
   integer :: j1, i1, k, l, atom_2, FN, i
   real(8) :: temp, epsylon, Ev, SH_1
   real(8), pointer :: x, y, z
   integer, pointer :: KOA1, KOA2, m
   character(200) :: Error_descript
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part


   Error_descript = ''
   epsylon = 1d-12	! acceptable tolerance : how small an overlap integral can be, before we set it to zero
   ! size of the basis set per atom, below:
   n_orb = identify_DFTB_orbitals_per_atom(numpar%basis_size_ind)  ! module "TB_DFTB"
   nat = Scell(NSC)%Na  ! number of atoms in the supercell
   Nsiz = size(Scell(NSC)%Ha,1) ! size of the total basis set

   if (.not.allocated(Sij)) allocate(Sij(Nsiz,Nsiz))
   if (.not.allocated(Hij)) allocate(Hij(Nsiz,Nsiz))
   Sij = 0.0d0
   Hij = 0.0d0

   ! Identify flag for scc calcilations:
   if (present(scc) .and. present(H_scc_0) .and. present(H_scc_1)) then
      do_scc = scc   ! follow the flag that was passed into here
   else
      do_scc = 0  ! no scc required
   endif

!    print*, 'Hamil_tot_3TB test 0'

   !-----------------------------------

   WNTSCC:if (do_scc == 0) then   ! no scc calculations, construct regular (zero-order) Hamiltonian:

      ! 1) Construct non-orthogonal Hamiltonian H and Overlap matrix S in 2 steps:
#ifdef MPI_USED   ! use the MPI version
      if (.not.allocated(Hij1)) allocate(Hij1(n_orb,n_orb), source = 0.0d0)
      if (.not.allocated(Sij1)) allocate(Sij1(n_orb,n_orb), source = 0.0d0)
      N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
      Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
      Nend = nat
      ! Do the cycle (parallel) calculations:
      do j = Nstart, Nend, N_incr  ! each process does its own part
      !do j = 1,nat	! atom #1
         m => Scell(NSC)%Near_neighbor_size(j)

         do atom_2 = 0,m ! do only for atoms close to that one

            if (atom_2 == 0) then ! the same atom
               i = j    ! atom #2 = atom #1, onsite
            else  ! different atoms
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! atom #2
            endif

            IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
               ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
               ! Contruct a block-hamiltonian:
               call Hamilton_one_3TB(numpar%basis_size_ind, Scell(NSC), j, i, TB_Hamil, Hij1, &
                                    M_Vij, M_Lag_exp, M_lmn, Mjs)   ! below

               ! Construct overlap matrix for this pair of atoms:
               call Get_overlap_S_matrix_3TB(numpar%basis_size_ind, j, i, Sij1, M_SVij, M_lmn)  ! below

               do j1 = 1,n_orb ! all orbitals of atom #1
                  l = (j-1)*n_orb+j1   ! atom #1 (j)
                  do i1 = 1,n_orb ! all orbitals of atom #2
                     k = (i-1)*n_orb+i1   ! atom #2 (i)
                     ! We fill the upper triangle here
                     ! (the order does not matter, just have to be consistent with the Slater-Koster functions):
                     Hij(l,k) = Hij1(j1,i1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian
                     Sij(l,k) = Sij1(j1,i1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices|
                  enddo ! i1
               enddo ! j1
            endif IJ
         enddo ! j
      enddo ! i
      deallocate(Hij1, Sij1)

      ! b) Construct lower triangle - use symmetry:
      do j = Nstart, Nend, N_incr  ! each process does its own part
      !do j = 1,nat	! all atoms
         m => Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 1,m ! do only for atoms close to that one
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            if (i < j) then ! lower triangle
               do j1 = 1,n_orb ! all orbitals of atom #1
                  l = (j-1)*n_orb+j1   ! atom #1
                  do i1 = 1,n_orb ! all orbitals of atom #2
                     k = (i-1)*n_orb+i1   ! atom #2
                     Hij(l,k) = Hij(k,l)
                     Sij(l,k) = Sij(k,l)
                  enddo ! i1
               enddo ! j1
            endif
         enddo ! j
      enddo ! i

      ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
      error_part = 'Error in Hamil_tot_3TB:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{Hij}', Hij) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{Sij}', Sij) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
!$omp parallel private(j, m, atom_2, i, j1, l, i1, k, Hij1, Sij1)
      if (.not.allocated(Hij1)) allocate(Hij1(n_orb,n_orb), source = 0.0d0)
      if (.not.allocated(Sij1)) allocate(Sij1(n_orb,n_orb), source = 0.0d0)
!$omp do
      do j = 1,nat	! atom #1
         m => Scell(NSC)%Near_neighbor_size(j)

         do atom_2 = 0,m ! do only for atoms close to that one

            if (atom_2 == 0) then ! the same atom
               i = j    ! atom #2 = atom #1, onsite
            else  ! different atoms
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! atom #2
            endif

            IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
               ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
               ! Contruct a block-hamiltonian:
               call Hamilton_one_3TB(numpar%basis_size_ind, Scell(NSC), j, i, TB_Hamil, Hij1, &
                                    M_Vij, M_Lag_exp, M_lmn, Mjs)   ! below

               ! Construct overlap matrix for this pair of atoms:
               call Get_overlap_S_matrix_3TB(numpar%basis_size_ind, j, i, Sij1, M_SVij, M_lmn)  ! below

               do j1 = 1,n_orb ! all orbitals of atom #1
                  l = (j-1)*n_orb+j1   ! atom #1 (j)
                  do i1 = 1,n_orb ! all orbitals of atom #2
                     k = (i-1)*n_orb+i1   ! atom #2 (i)
                     ! We fill the upper triangle here
                     ! (the order does not matter, just have to be consistent with the Slater-Koster functions):
                     ! [HA 0]
                     Hij(l,k) = Hij1(j1,i1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian
                     Sij(l,k) = Sij1(j1,i1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices|
                     ! [HA 1]
!                      Hij(l,k) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian
!                      Sij(l,k) = Sij1(i1,j1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
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
            if (i < j) then ! lower triangle
               do j1 = 1,n_orb ! all orbitals of atom #1
                  l = (j-1)*n_orb+j1   ! atom #1
                  do i1 = 1,n_orb ! all orbitals of atom #2
                     k = (i-1)*n_orb+i1   ! atom #2
                     Hij(l,k) = Hij(k,l)
                     Sij(l,k) = Sij(k,l)
                  enddo ! i1
               enddo ! j1
            endif
         enddo ! j
      enddo ! i
!$omp end do
!$omp end parallel
#endif

      ! Check if Hamiltonian is symmetric (for testing purpuses):
!     call check_symmetry(Hij, numpar%MPI_param) ! module "Algebra_tools"

      ! 2)    ! Save the non-orthogonalized Hamiltonian:
      Scell(NSC)%H_non = Hij  ! nondiagonalized Hamiltonian
      Scell(NSC)%Sij = Sij    ! save Overlap matrix

   else WNTSCC ! scc cycle, construct only use second order scc correction:
      Hij = H_scc_0 + H_scc_1 ! zero and second order scc contributions into Hamiltonian
      Scell(NSC)%H_non = Hij  ! save new nondiagonalized Hamiltonian
      Sij = Scell(NSC)%Sij
   endif WNTSCC


   !-----------------------------------
   ! 3) Orthogonalize the Hamiltonian using Lowedin procidure
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization(Nsiz, Sij, Hij, Err, Scell(NSC)%eigen_S) ! module "TB_NRL"

   Scell(NSC)%Hij = Hij ! save orthogonalized but non-diagonalized Hamiltonian

   forall (i = 1:size(Hij,1),  j = 1:size(Hij,2), (ABS(Hij(i,j)) < 1.0d-10))
      Hij(i,j) = 0.0d0
   endforall
   !-----------------------------------
   ! 4) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
   call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript) ! module "Algebra_tools"
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Hamil_tot_DFTB: '//trim(adjustl(Error_descript))
      if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         call Save_error_details(Err, 6, Error_descript)
      endif
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
#ifdef MPI_USED   ! use the MPI version
      do j = Nstart, Nend, N_incr  ! each process does its own part
      !do j = 1,nat	! all atoms
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

       ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
      error_part = 'Error in Hamil_tot_3TB:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Scell(NSC)%PRRx', Scell(NSC)%PRRx) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Scell(NSC)%PRRy', Scell(NSC)%PRRy) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Scell(NSC)%PRRz', Scell(NSC)%PRRz) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
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
#endif

      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0

      Scell(NSC)%PRRx = Scell(NSC)%PRRx * temp
      Scell(NSC)%PRRy = Scell(NSC)%PRRy * temp
      Scell(NSC)%PRRz = Scell(NSC)%PRRz * temp
   endif ! (numpar%optic_model .EQ. 3)

   nullify(KOA1, KOA2, m, x, y, z)
   deallocate(Hij, Sij)
end subroutine Hamil_tot_3TB



subroutine Get_overlap_S_matrix_3TB(basis_ind, i, j, Sij, M_SVij, M_lmn)
   integer, intent(in) :: basis_ind ! index of the basis set: 0=s, 1=sp3, 2=sp3d5
   integer, intent(in) :: i, j  ! atomic indices
   real(8), dimension(:,:), intent(out) :: Sij  ! Overlap matrix
   real(8), dimension(:,:,:), intent(in) :: M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   !---------------------------------------
   real(8), dimension(:), allocatable :: vec_M_SVij12, vec_M_SVij21
   integer :: k
   if (i /= j) then	! it's 2 different atoms:
      ! Construct the overlap matrix including angular part
      select case (basis_ind)
      case (0) ! for s basis set:
         allocate(vec_M_SVij12(1))
         vec_M_SVij12(:) = M_SVij(i,j,:)
         call KS_s(vec_M_SVij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij)	! module "TB_Koster_Slater"
      case (1) ! for sp3 basis set:
         allocate(vec_M_SVij12(4))
         vec_M_SVij12(:) = M_SVij(i,j,:)
         allocate(vec_M_SVij21(4))
         vec_M_SVij21(:) = M_SVij(j,i,:)
         !call KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij)   ! module "TB_Koster_Slater"
         call KS_sp3_hetero_TEST(vec_M_SVij12, vec_M_SVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij)   ! module "TB_Koster_Slater"
      case default ! for sp3d5 basis set:
         allocate(vec_M_SVij12(10))
         vec_M_SVij12(:) = M_SVij(i,j,:)
         allocate(vec_M_SVij21(10))
         vec_M_SVij21(:) = M_SVij(j,i,:)
         !call KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij) ! module "TB_Koster_Slater"
         call KS_sp3d5_hetero_TEST(vec_M_SVij12, vec_M_SVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij) ! module "TB_Koster_Slater"
      end select
   else	! it is the same atom
      Sij = 0.0d0
      forall (k=1:size(Sij,1)) Sij(k,k)=1.0d0
   endif
   if (allocated(vec_M_SVij12)) deallocate(vec_M_SVij12)
   if (allocated(vec_M_SVij21)) deallocate(vec_M_SVij21)
end subroutine Get_overlap_S_matrix_3TB



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
   real(8) :: G_IJK(3,3), H_temp(6), MM(4)
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
      ! This subroutine is used because we fill the LOWER triangle in the Hamiltonian matrix:
      !call KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
      ! This subroutine could be used if we were filling the UPPER triangle first:
      call KS_sp3_hetero_TEST(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   case default ! for sp3d5 basis set:
      allocate(vec_M_SVij12(10))
      vec_M_SVij12(:) = M_Vij12(:)
      allocate(vec_M_SVij21(10))
      vec_M_SVij21(:) = M_Vij21(:)
      !call KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
      ! Same here, the upper triangle vs. the lower one:
      call KS_sp3d5_hetero_TEST(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
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

!!!$omp PARALLEL ! it is already inside parallelized region, do not parallelize again!
!!!$omp do private(atom_3, k, KOA3, at_ind, sh1, sh2)
      AT3:do atom_3 = 1,m ! do only for atoms close to that one
         k = Scell%Near_neighbor_list(i,atom_3) ! this is the list of such close atoms
         ! Make sure the third atoms is not the second atom:
         if (k /= j) then
            KOA3 => Scell%MDatoms(k)%KOA   ! kind of atom #3

            ! Find the combination-of-atoms index:
            at_ind = find_3bdy_ind(KOA1, KOA2, KOA3)  ! module "Dealing_with_3TB"

            ! Factors including LAguerres are independent of shells:
            MM(1) = M_Lag_exp(i,k,1) * M_Lag_exp(j,k,1)
            MM(2) = M_Lag_exp(i,k,1) * M_Lag_exp(j,k,2)
            MM(3) = M_Lag_exp(i,k,2) * M_Lag_exp(j,k,1)
            MM(4) = MM(1) * M_Lag_exp(i,j,1)

            ! Get the radial function for 3-body interaction:
            do sh1 = 1, 1+basis_ind
               do sh2 = 1, 1+basis_ind
                    G_IJK(sh1,sh2) = TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 1) * MM(1) + &
                                     TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 2) * MM(2) + &
                                     TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 3) * MM(3) + &
                                     TB_Hamil(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 4) * MM(4)
               enddo
            enddo
            ! Include angular parts:
            H_3bdy_part(1,1) = H_3bdy_part(1,1) + G_IJK(1,1)  ! s-s * s-s

            if (basis_ind > 0) then ! p3 orbitals:

               H_3bdy_part(1,2:4) = H_3bdy_part(1,2:4) + G_IJK(2,1) * Mjs(j,k,2:4)   ! s-s * p[x,y,z]-s
               H_3bdy_part(2:4,1) = H_3bdy_part(2:4,1) + G_IJK(1,2) * Mjs(i,k,2:4)   ! p[x,y,z]-s * s-s

               ! Calculate repeating part of the K-S matrix elements:
               H_temp(1) = G_IJK(2,2) * Mjs(i,k,2) ! G * px-s
               H_temp(2) = G_IJK(2,2) * Mjs(i,k,3) ! G * px-s
               H_temp(3) = G_IJK(2,2) * Mjs(i,k,4) ! G * px-s
               ! Add it to the Hamiltonian part:
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
                  !H_3bdy_part(5:9,5) = H_3bdy_part(5:9,5) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,6) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dxz-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dxz-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dxz-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dxz-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dxz-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(6,5:9) = H_3bdy_part(6,5:9) + H_temp(1:5)
                  !H_3bdy_part(5:9,6) = H_3bdy_part(5:9,6) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,7) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dyz-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dyz-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dyz-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dyz-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dyz-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(7,5:9) = H_3bdy_part(7,5:9) + H_temp(1:5)
                  !H_3bdy_part(5:9,7) = H_3bdy_part(5:9,7) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,8) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! (dx2-y2)-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! (dx2-y2)-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! (dx2-y2)-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! (dx2-y2)-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! (dx2-y2)-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(8,5:9) = H_3bdy_part(8,5:9) + H_temp(1:5)
                  !H_3bdy_part(5:9,8) = H_3bdy_part(5:9,8) + H_temp(1:5)

                  ! Calculate the K-S matrix elements:
                  H_temp(6) = Mjs(i,k,9) * G_IJK(3,3)
                  H_temp(1) = H_temp(6) * Mjs(j,k,5) ! (d3z2-r2)-s * dxy-s
                  H_temp(2) = H_temp(6) * Mjs(j,k,6) ! (d3z2-r2)-s * dxz-s
                  H_temp(3) = H_temp(6) * Mjs(j,k,7) ! (d3z2-r2)-s * dyz-s
                  H_temp(4) = H_temp(6) * Mjs(j,k,8) ! (d3z2-r2)-s * (dx2-y2)-s
                  H_temp(5) = H_temp(6) * Mjs(j,k,9) ! (d3z2-r2)-s * (d3z2-r2)-s
                  ! Add it into the Hamiltonian part:
                  H_3bdy_part(9,5:9) = H_3bdy_part(9,5:9) + H_temp(1:5)
                  !H_3bdy_part(5:9,9) = H_3bdy_part(5:9,9) + H_temp(1:5)

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

   !-----------------
   ! 2) Average term (testet, correct):
!!!$omp parallel private(atom_2, j, i1, term_s, term_p, term_d) ! in is called from parallelized region, don't parallelize it again
!!!$omp do reduction( + : H_avg)
   do atom_2 = 1, m ! do only for atoms close to that one
      j = Scell%Near_neighbor_list(i, atom_2) ! this is the list of such close atoms
      ! [OS 0]: divergent
!       KOA1 => Scell%MDatoms(i)%KOA
!       KOA2 => Scell%MDatoms(j)%KOA
      ! [OS 1] tested, correct:
      KOA1 => Scell%MDatoms(j)%KOA
      KOA2 => Scell%MDatoms(i)%KOA

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

      !-----------------
      ! 3) Crystal field:
      ! [OS 0]:
!       KOA2 => Scell%MDatoms(j)%KOA
!       KOA1 => Scell%MDatoms(i)%KOA
      ! [OS 1] tested, correct:
      KOA1 => Scell%MDatoms(j)%KOA
      KOA2 => Scell%MDatoms(i)%KOA

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


         ! Diagonal part is excluded (no self-interaction of orbitals):
         if (TB(KOA1,KOA1)%nullify_diag_cf) then
            H_cf(2,2) = 0.0d0 ! px-s * px-s
            H_cf(3,3) = 0.0d0 ! py-s * py-s
            H_cf(4,4) = 0.0d0 ! pz-s * pz-s
         else
            H_cf(2,2) = H_cf(2,2) + H_cf_temp(2)  ! px-s * px-s
            H_cf(3,3) = H_cf(3,3) + H_cf_temp(5) * Mjs(i,j,3)  ! py-s * py-s
            H_cf(4,4) = H_cf(4,4) + matr_spd(2,2) * Mjs(i,j,4) * Mjs(i,j,4) ! pz-s * pz-s
         endif

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
            if (TB(KOA1,KOA1)%nullify_diag_cf) then
               H_cf_temp(1) = 0.0d0                     ! dxy-s * dxy-s
            else
               H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! dxy-s * dxy-s
            endif
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! dxy-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! dxy-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! dxy-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! dxy-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(5,5:9) = H_cf(5,5:9) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,6) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! dxz-s * dxy-s
            if (TB(KOA1,KOA1)%nullify_diag_cf) then
               H_cf_temp(2) = 0.0d0                     ! dxz-s * dxz-s
            else
               H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! dxz-s * dxz-s
            endif
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! dxz-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! dxz-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! dxz-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(6,5:9) = H_cf(6,5:9) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,7) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! dyz-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! dyz-s * dxz-s
            if (TB(KOA1,KOA1)%nullify_diag_cf) then
               H_cf_temp(3) = 0.0d0                     ! dyz-s * dyz-s
            else
               H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! dyz-s * dyz-s
            endif
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! dyz-s * (dx2-y2)-s
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! dyz-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(7,5:9) = H_cf(7,5:9) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,8) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! (dx2-y2)-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! (dx2-y2)-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! (dx2-y2)-s * dyz-s
            if (TB(KOA1,KOA1)%nullify_diag_cf) then
               H_cf_temp(4) = 0.0d0                     ! (dx2-y2)-s * (dx2-y2)-s
            else
               H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! (dx2-y2)-s * (dx2-y2)-s
            endif
            H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! (dx2-y2)-s * (d3z2-r2)-s
            ! Add it into the Hamiltonian part:
            H_cf(8,5:9) = H_cf(8,5:9) + H_cf_temp(1:5)

            ! Calculate the K-S matrix elements:
            H_cf_temp(6) = Mjs(i,j,9) * matr_spd(3,3)
            H_cf_temp(1) = H_cf_temp(6) * Mjs(i,j,5) ! (d3z2-r2)-s * dxy-s
            H_cf_temp(2) = H_cf_temp(6) * Mjs(i,j,6) ! (d3z2-r2)-s * dxz-s
            H_cf_temp(3) = H_cf_temp(6) * Mjs(i,j,7) ! (d3z2-r2)-s * dyz-s
            H_cf_temp(4) = H_cf_temp(6) * Mjs(i,j,8) ! (d3z2-r2)-s * (dx2-y2)-s
            if (TB(KOA1,KOA1)%nullify_diag_cf) then
               H_cf_temp(5) = 0.0d0                     ! (d3z2-r2)-s * (d3z2-r2)-s
            else
               H_cf_temp(5) = H_cf_temp(6) * Mjs(i,j,9) ! (d3z2-r2)-s * (d3z2-r2)-s
            endif
            ! Add it into the Hamiltonian part:
            H_cf(9,5:9) = H_cf(9,5:9) + H_cf_temp(1:5)

         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)
   enddo ! atom_2 = 0, m
!!!$omp enddo
!!!$omp end parallel

   !-----------------
   ! 4) 3-body contributions:
   if (TB(KOA1,KOA1)%include_3body) then  ! only if user defined it to include
      atom_2 = 0  ! to restart
!!!$omp PARALLEL
!!!$omp do private(atom_2, j, KOA2, atom_3, k, KOA3, at_ind, term_s, sh1)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         KOA2 => Scell%MDatoms(j)%KOA   ! index of the second atom

         ! To start summing up
         AT3:do atom_3 = 1,m ! do only for atoms close to that one
            k = Scell%Near_neighbor_list(i,atom_3) ! this is the list of such close atoms
            ! Make sure the third atom is not the second atom:
            if (k /= j) then
               KOA3 => Scell%MDatoms(k)%KOA   ! kind of atom #3

               ! Find the combination-of-atoms index:
               at_ind = find_3bdy_ind(KOA1, KOA2, KOA3)  ! module "Dealing_with_3TB"

               ! Get the 3-body elements for each shells combination:

               term_s = TB(KOA1,KOA2)%Hh3bdy(at_ind, 1) * M_Lag_exp(i,j,1) * M_Lag_exp(j,k,1) * M_Lag_exp(i,k,1) + &
                        TB(KOA1,KOA2)%Hh3bdy(at_ind, 2) * M_Lag_exp(i,j,2) * M_Lag_exp(j,k,1) * M_Lag_exp(i,k,1) + &
                        TB(KOA1,KOA2)%Hh3bdy(at_ind, 3) * M_Lag_exp(i,j,1) * M_Lag_exp(j,k,2) * M_Lag_exp(i,k,1) + &
                        TB(KOA1,KOA2)%Hh3bdy(at_ind, 4) * M_Lag_exp(i,j,1) * M_Lag_exp(j,k,1) * M_Lag_exp(i,k,2) ! Eq.(22) in [1]

               do sh1 = 1, Bsiz  ! for all shells
                  H_3bdy(sh1,sh1) = H_3bdy(sh1,sh1) + term_s   ! no orbital resolution
               enddo ! sh1
            endif ! (k /= i)
         enddo AT3
      enddo AT2
!!!$omp end do
!!!$omp END PARALLEL
   endif ! (TB(KOA1,KOA1)%include_3body)


5001 continue

   ! Collect all the terms into Hamiltonian:
   Hij = E_onsite + H_avg + H_cf + H_3bdy

   deallocate(E_onsite, H_avg, H_cf, H_3bdy)
   nullify(m, KOA1, KOA2, KOA3)
end subroutine Onsite_3TB





!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Derivatives of the attractive part of TB:

! Subroutine for derivative of the Hamiltonian:
subroutine get_dHij_drij_3TB(numpar, Scell, NSC, TB, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei, Mjs, M_Lag_exp, M_d_Lag_exp)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp, M_d_Lag_exp
   !------------------------------------------------------------
   integer :: nat, k
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part
   !------------------------------------------------------------

   nat = size(Scell(NSC)%MDatoms)	! number of atoms
#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = nat
   do k = 1, nat
      Scell(NSC)%MDatoms(k)%forces%att(:) = 0.0d0	! just to start
   enddo
   ! Do the cycle (parallel) calculations:
   ATOMS:do k = Nstart, Nend, N_incr  ! each process does its own part
      ! (tested, correct):
      call get_forces_3TB(k, numpar, Scell, NSC, TB, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, &
                           M_lmn, M_Lag_exp, M_d_Lag_exp, Aij_x_Ei, Mjs) !below

   enddo ATOMS
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in Construct_Vij_BOP:'
   do k = 1, nat
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'MDatoms(k)%forces%att(:)', Scell(NSC)%MDatoms(k)%forces%att(:)) ! module "MPI_subroutines"
   enddo

#else ! use OpenMP instead
   !$omp PARALLEL private(k)
   !$omp do
   ATOMS:do k = 1, nat	! forces for all atoms
      Scell(NSC)%MDatoms(k)%forces%att(:) = 0.0d0	! just to start

      ! (tested, correct):
      call get_forces_3TB(k, numpar, Scell, NSC, TB, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, &
                           M_lmn, M_Lag_exp, M_d_Lag_exp, Aij_x_Ei, Mjs) !below

   enddo ATOMS
   !$omp end do
   !$omp end parallel
#endif
end subroutine get_dHij_drij_3TB



subroutine get_forces_3TB(k, numpar, Scell, NSC, TB, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, M_Lag_exp, M_d_Lag_exp, Aij_x_Ei, Mjs)
   integer, intent(in) :: k	! forces for this atom
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of the supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp, M_d_Lag_exp
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   real(8), dimension(:,:,:), intent(in) :: Mjs ! matrix of overlaps with s-orbital
   !-----------------------------
   integer :: j, m, atom_2, i, j1, i1, nat, Nsiz, n_orb, i4, j4, sh2, sh1
   real(8), dimension(:,:,:), allocatable :: dH, dS
   real(8), dimension(:,:,:), allocatable :: dH1, dS1    ! for various basis sets

   nat = size(Scell(NSC)%MDatoms)	! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1)	! total number of orbitals
   n_orb =  identify_DFTB_orbitals_per_atom(numpar%basis_size_ind)  ! size of the basis set per atom, below

   if (.not.allocated(dH1)) allocate(dH1(3,n_orb,n_orb))
   if (.not.allocated(dS1)) allocate(dS1(3,n_orb,n_orb))
   dH1 = 0.0d0
   dS1 = 0.0d0
   if (.not.allocated(dH)) allocate(dH(3,Nsiz,Nsiz))
   if (.not.allocated(dS)) allocate(dS(3,Nsiz,Nsiz))
   dH = 0.0d0
   dS = 0.0d0

   ! 1) Construct the derivatives of the Hamiltonian and Overlap matrix:
   ! Construct upper triangle - calculate each element:
   ATOM1:do i = 1,nat	! all pairs of atoms contribute to the force (to some degree)
      i4 = (i-1)*n_orb	! orbitals
      m = Scell(NSC)%Near_neighbor_size(i)
      ATOM2:do atom_2 = 0,m ! do only for atoms close to that one
         if (atom_2 == 0) then	! the same atom
            j = i
         else	! two different atoms
            j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         endif
         j4 = (j-1)*n_orb	! orbitals
         IJ:if (j >= i) then ! it's a new pair of atoms, calculate everything
            ! Get the matrix of dH/drij and dS/drij:
            call d_Hamilton_one_3TB(numpar%basis_size_ind, k, Scell(NSC), TB, i, j, atom_2, dH1, M_Vij, M_dVij, &
                  M_lmn, M_Lag_exp, M_d_Lag_exp, dS1, M_SVij, M_dSVij, Mjs) ! this calls the block-hamiltonian

            do j1 = 1,n_orb	! all orbitals
               sh2 = j4+j1
               do i1 = 1,n_orb	! all orbitals
                  sh1 = i4+i1
                  ! [FA 1]
                  dH(:,sh1,sh2) = dH1(:,i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  dS(:,sh1,sh2) = dS1(:,i1,j1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
                  ! [FA 2]
!                   dH(:,sh1,sh2) = dH1(:,j1,i1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
!                   dS(:,sh1,sh2) = dS1(:,j1,i1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices

               enddo ! i1
            enddo ! j1
         endif IJ
      enddo ATOM2
   enddo ATOM1
   ! Clean up the temporary arrays:
   deallocate(dH1, dS1)

   ! b) Construct lower triangle - use symmetry:
   do i = 1,nat	! all pairs of atoms contribute to the force (to some degree)
      i4 = (i-1)*n_orb
      m = Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         j4 = (j-1)*n_orb
         if (j < i) then ! this pair of atoms was already calculated, use the symmetry
            do j1 = 1,n_orb	! all orbitals
               sh2 = j4+j1
               do i1 = 1,n_orb	! all orbitals
                  sh1 = i4+i1
                  dH(:,sh1,sh2) = dH(:,sh2,sh1)
                  dS(:,sh1,sh2) = dS(:,sh2,sh1)
               enddo ! i1
            enddo ! j1
         endif ! (j < i)
      enddo ! i
   enddo ! atom_2

   ! 2) Calculate the forces form the derivatives and the eigervectors:
   call Attract_TB_forces_3TB(Aij, Aij_x_Ei, dH, dS, Scell, NSC, Scell(NSC)%MDatoms(k)%forces%att(:), n_orb)

   deallocate(dH, dS)
end subroutine get_forces_3TB



subroutine Attract_TB_forces_3TB(Aij, Aij_x_Ei, dH, dS, Scell, NSC, Eelectr_s, n_orb)
   real(8), dimension(:,:,:), intent(in) :: dH, dS
   real(8), dimension(:,:), intent(in), target :: Aij, Aij_x_Ei
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(out)  :: Eelectr_s ! part of the forces
   integer, intent(in) :: n_orb ! number of orbitals in the used basis set
   !--------------------------------------
   integer :: j, k, i, ste, n, norb_1, j_norb, n_too
   integer, target :: i2
   integer, pointer :: m, j1

   n = size(Aij,1)	! total number of orbitals
   norb_1 = n_orb - 1
   Eelectr_s = 0.0d0
   n_too = 0 ! to start with

   i2 = 0
   ste = 1
   do i = 1, n	! all orbitals
       if (i .GE. ste) then
          i2 = i2 + 1	! number atom corresponding to this set of orbitals
          ste = ste + n_orb	! counter to switch to the next atom
       endif
       m => Scell(NSC)%Near_neighbor_size(i2)
       do k = 0, m	! all pairs of atoms (to get all orbitals from them)
          if (k == 0) then	! the same atom (contribution from environment)
             j1 => i2
          else	! contribution from a different atom
             j1 => Scell(NSC)%Near_neighbor_list(i2,k) ! this is the list of such close atoms
          endif
          j = (j1-1)*n_orb + 1
          if (j .GT. 0) then
             j_norb = j+norb_1
             Eelectr_s(1) = Eelectr_s(1) + SUM(dH(1,i,j:j_norb)*Aij(i,j:j_norb) - dS(1,i,j:j_norb)*Aij_x_Ei(i,j:j_norb))
             Eelectr_s(2) = Eelectr_s(2) + SUM(dH(2,i,j:j_norb)*Aij(i,j:j_norb) - dS(2,i,j:j_norb)*Aij_x_Ei(i,j:j_norb))
             Eelectr_s(3) = Eelectr_s(3) + SUM(dH(3,i,j:j_norb)*Aij(i,j:j_norb) - dS(3,i,j:j_norb)*Aij_x_Ei(i,j:j_norb))

             ! Mark if there is a potential trouble:
             if (maxval(ABS(Eelectr_s(:))) .GE. 1.0d6) n_too = n_too + 1
!              if (maxval(ABS(Eelectr_s(:))) .GE. 1.0d7) then
!                 write(*,'(a)') 'Trouble in subroutine Attract_TB_forces_NRL, too large attractive force:'
!                 write(*,'(i12,i12)', advance='no') i, j
!                 write(*,'(e25.16,$)') dH(1,i,j:j_norb), Aij(i,j:j_norb)
!                 write(*,'(e25.16,e25.16,e25.16)') Eelectr_s(:)
!              endif

          endif
       enddo
   enddo

   if (n_too > 0) then
      write(*,'(a)') 'Trouble in subroutine Attract_TB_forces_3TB, too large attractive force'
      write(*,'(a,i)') 'For elements: ', n_too
   endif

   nullify(m, j1)
end subroutine Attract_TB_forces_3TB


!ddddddddddddddddddddddddddddddddddddddddddddddddddd
! Derivatives:
subroutine d_Hamilton_one_3TB(basis_ind, k, Scell, TB, i, j, atom_2, dH, M_Vij, M_dVij, M_lmn, M_Lag_exp, M_d_Lag_exp, &
                              dS, M_SVij, M_dSVij, Mjs_in)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(TB_H_3TB), dimension(:,:), intent(in), target :: TB	! parameters of the Hamiltonian of TB
   integer, intent(in) :: i, j, atom_2, k
   real(8), dimension(:,:,:), intent(out) :: dH, dS  ! hamiltonian, all orbitals
   real(8), dimension(:,:,:), intent(in), target :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp, M_d_Lag_exp   ! Laguerres and their derivatives
   real(8), dimension(:,:,:), intent(in), target :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms for S-matrix, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: Mjs_in ! matrix of overlaps with s-orbital
   !---------------------------------------
   real(8), dimension(9) :: M_dlmn	! dl/dx, dl/dy, dl/dz, dm/dx, dm/dy, dm/dz, dn/dx, dn/dy, dn/dz
   real(8), dimension(9) :: M_dlmn_ik, M_dlmn_jk
   real(8), dimension(3, 8) :: dMis, dMjs ! 3p, 5d
   real(8), dimension(9) :: Mjs, Mjks, Miks ! 1s, 3p, 5d
   real(8) :: drij_dsk(3), drikk_dsk(3), drjkk_dsk(3), G_IJK(3,3), dG_IJK(3,3,3), H_temp(6), MM(4), dMM(3,4)
   integer :: n_orb, n_overlap, kk, atom_3, sh1, sh2, at_ind
   real(8), dimension(:), allocatable :: vec_M_Vij12, vec_M_Vij21, vec_M_SVij12, vec_M_SVij21
   real(8), dimension(:), allocatable :: vec_M_dVij12, vec_M_dVij21, vec_M_dSVij12, vec_M_dSVij21
   real(8), dimension(:,:), allocatable :: dH1, dS1
   real(8), dimension(:,:,:), allocatable :: dH_3bdy
   real(8), pointer :: x1, y1, z1, r1
   integer, pointer :: KOA1, KOA2, KOA3, m
   real(8), pointer :: rij, xij, yij, zij
   real(8), pointer :: rik, xik, yik, zik
   real(8), pointer :: rjk, xjk, yjk, zjk

   if (i == j) then ! Onsite contributions (incl. 3-body parts etc.):

      call d_Onsite_3TB(basis_ind, k, i, Scell, TB, M_lmn, M_Lag_exp, M_d_Lag_exp, Mjs_in, dH) ! below
      dS = 0.0d0  ! Here are constants on-site, derivatives = 0

   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:

      ! Depending on the basis set:
      n_overlap = identify_DFTB_basis_size(basis_ind)   ! below
      n_orb = identify_DFTB_orbitals_per_atom(basis_ind)    ! below
      allocate(vec_M_Vij12(n_overlap))
      allocate(vec_M_Vij21(n_overlap))
      allocate(vec_M_SVij12(n_overlap))
      allocate(vec_M_SVij21(n_overlap))
      allocate(vec_M_dVij12(n_overlap))
      allocate(vec_M_dVij21(n_overlap))
      allocate(vec_M_dSVij12(n_overlap))
      allocate(vec_M_dSVij21(n_overlap))
      allocate(dH1(n_orb,n_orb))
      allocate(dS1(n_orb,n_orb))

      x1 => Scell%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
      y1 => Scell%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
      z1 => Scell%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
      r1 => Scell%Near_neighbor_dist(i,atom_2,4)	! at this distance, R

      ! Derivatives of rij by sk:
      ! All functions are from the module "TB_Koster_Slater":
      drij_dsk(1) = drij_dska(i, j, k, x1, y1, z1, r1, Scell%supce, 1, .true.)	! dr_{ij}/ds_{k,x}
      drij_dsk(2) = drij_dska(i, j, k, x1, y1, z1, r1, Scell%supce, 2, .true.)	! dr_{ij}/ds_{k,y}
      drij_dsk(3) = drij_dska(i, j, k, x1, y1, z1, r1, Scell%supce, 3, .true.)	! dr_{ij}/ds_{k,z}

      M_dlmn(1) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 1, drij_dsk(1))	! dl/dsx
      M_dlmn(2) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 2, drij_dsk(2))	! dl/dsy
      M_dlmn(3) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 3, drij_dsk(3))	! dl/dsz
      M_dlmn(4) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 1, drij_dsk(1))	! dm/dsx
      M_dlmn(5) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 2, drij_dsk(2))	! dm/dsy
      M_dlmn(6) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 3, drij_dsk(3))	! dm/dsz
      M_dlmn(7) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 1, drij_dsk(1))	! dn/dsx
      M_dlmn(8) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 2, drij_dsk(2))	! dn/dsy
      M_dlmn(9) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 3, drij_dsk(3))	! dn/dsz

      ! Create vectors from the slices of arrays to pass into subroutines:
      vec_M_Vij21 = M_Vij(i,j,:)
      vec_M_Vij12 = M_Vij(j,i,:)
      vec_M_SVij21 = M_SVij(i,j,:)
      vec_M_SVij12 = M_SVij(j,i,:)

      ! Derivative along dx:
      dH1 = 0.0d0
      dS1 = 0.0d0
      vec_M_dVij21(:) = M_dVij(i,j,:)*drij_dsk(1)
      vec_M_dVij12(:) = M_dVij(j,i,:)*drij_dsk(1)
      vec_M_dSVij21(:) = M_dSVij(i,j,:)*drij_dsk(1)
      vec_M_dSVij12(:) = M_dSVij(j,i,:)*drij_dsk(1)

      select case (basis_ind)
      case (0) ! s
         call d_KS_s(vec_M_Vij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dH1)	! module "TB_Koster_Slater"
         call d_KS_s(vec_M_SVij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dS1)	! module "TB_Koster_Slater"
      case (1) ! sp3
         call d_KS_sp3_hetero(vec_M_Vij12, vec_M_Vij21, vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dS1)	! module "TB_Koster_Slater"
      case (2) ! sp3d5
         call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21,  vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dS1)	! module "TB_Koster_Slater"
      end select
      dH(1,:,:) = dH1(:,:)
      dS(1,:,:) = dS1(:,:)

      ! Derivative along dy:
      dH1 = 0.0d0
      dS1 = 0.0d0
      vec_M_dVij21(:) = M_dVij(i,j,:)*drij_dsk(2)
      vec_M_dVij12(:) = M_dVij(j,i,:)*drij_dsk(2)
      vec_M_dSVij21(:) = M_dSVij(i,j,:)*drij_dsk(2)
      vec_M_dSVij12(:) = M_dSVij(j,i,:)*drij_dsk(2)
      select case (basis_ind)
      case (0) ! s
         call d_KS_s(vec_M_Vij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dH1)	! module "TB_Koster_Slater"
         call d_KS_s(vec_M_SVij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dS1)	! module "TB_Koster_Slater"
      case (1) ! sp3
         call d_KS_sp3_hetero(vec_M_Vij12, vec_M_Vij21, vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dS1)	! module "TB_Koster_Slater"
      case (2) ! sp3d5
         call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21,  vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dS1)	! module "TB_Koster_Slater"
      end select

      dH(2,:,:) = dH1(:,:)
      dS(2,:,:) = dS1(:,:)

      ! Derivative along dz:
      dH1 = 0.0d0
      dS1 = 0.0d0
      vec_M_dVij21(:) = M_dVij(i,j,:)*drij_dsk(3)
      vec_M_dVij12(:) = M_dVij(j,i,:)*drij_dsk(3)
      vec_M_dSVij21(:) = M_dSVij(i,j,:)*drij_dsk(3)
      vec_M_dSVij12(:) = M_dSVij(j,i,:)*drij_dsk(3)
      select case (basis_ind)
      case (0) ! s
         call d_KS_s(vec_M_Vij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dH1)	! module "TB_Koster_Slater"
         call d_KS_s(vec_M_SVij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dS1)	! module "TB_Koster_Slater"
      case (1) ! sp3
         call d_KS_sp3_hetero(vec_M_Vij12, vec_M_Vij21, vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dS1)	! module "TB_Koster_Slater"
      case (2) ! sp3d5
         call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21,  vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dS1)	! module "TB_Koster_Slater"
      end select
      dH(3,:,:) = dH1(:,:)
      dS(3,:,:) = dS1(:,:)

      deallocate(vec_M_Vij12, vec_M_Vij21, vec_M_SVij12, vec_M_SVij21, vec_M_dVij12, vec_M_dVij21, vec_M_dSVij12, &
                  vec_M_dSVij21, dH1, dS1)

      ! 2) 3-body contributions:
      KOA1 => Scell%MDatoms(k)%KOA  ! Kind of atom:
      if (TB(KOA1,KOA1)%include_3body) then  ! only if user defined it to include
         allocate(dH_3bdy(3,n_orb,n_orb))
         ! Number of the nearest neighbors
         m => Scell%Near_neighbor_size(i)
         ! Kind of atom #2
         KOA2 => Scell%MDatoms(j)%KOA   ! index of the second atom
         ! Distances between atom #1 and #2:
         xij => Scell%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
         yij => Scell%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
         zij => Scell%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
         rij => Scell%Near_neighbor_dist(i,atom_2,4)	! at this distance, R
         ! d r_ij / d s_k:
         drij_dsk(1) = drij_dska(i, j, k, xij, yij, zij, rij, Scell%supce, 1, .true.)	! dr_{ij}/ds_{k,x}
         drij_dsk(2) = drij_dska(i, j, k, xij, yij, zij, rij, Scell%supce, 2, .true.)	! dr_{ij}/ds_{k,y}
         drij_dsk(3) = drij_dska(i, j, k, xij, yij, zij, rij, Scell%supce, 3, .true.)	! dr_{ij}/ds_{k,z}

         M_dlmn(1) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 1, drij_dsk(1))	! dl/dsx
         M_dlmn(2) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 2, drij_dsk(2))	! dl/dsy
         M_dlmn(3) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 3, drij_dsk(3))	! dl/dsz
         M_dlmn(4) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 1, drij_dsk(1))	! dm/dsx
         M_dlmn(5) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 2, drij_dsk(2))	! dm/dsy
         M_dlmn(6) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 3, drij_dsk(3))	! dm/dsz
         M_dlmn(7) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 1, drij_dsk(1))	! dn/dsx
         M_dlmn(8) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 2, drij_dsk(2))	! dn/dsy
         M_dlmn(9) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 3, drij_dsk(3))	! dn/dsz


         ! To start summing up:
         dH_3bdy = 0.0d0
         AT3:do atom_3 = 1,m ! do only for atoms close to that one
            kk = Scell%Near_neighbor_list(i,atom_3) ! this is the list of such close atoms
            ! Make sure the third atom is not the second atom:
            if (kk /= j) then
               KOA3 => Scell%MDatoms(kk)%KOA   ! kind of atom #3

               ! Distances between atom #1 and #3:
               xik => Scell%Near_neighbor_dist(i,atom_3,1)	! at this distance, X
               yik => Scell%Near_neighbor_dist(i,atom_3,2)	! at this distance, Y
               zik => Scell%Near_neighbor_dist(i,atom_3,3)	! at this distance, Z
               rik => Scell%Near_neighbor_dist(i,atom_3,4)	! at this distance, R
               ! Distances between atom #2 and #3:
               xjk => Scell%Near_neighbor_dist(j,atom_3,1)	! at this distance, X
               yjk => Scell%Near_neighbor_dist(j,atom_3,2)	! at this distance, Y
               zjk => Scell%Near_neighbor_dist(j,atom_3,3)	! at this distance, Z
               rjk => Scell%Near_neighbor_dist(j,atom_3,4)	! at this distance, R

               ! Derivatives of rij by sk:
               ! All functions are from the module "TB_Koster_Slater":
               drikk_dsk(1) = drij_dska(i, kk, k, xik, yik, zik, rik, Scell%supce, 1, .true.)	! dr_{ikk}/ds_{k,x}
               drikk_dsk(2) = drij_dska(i, kk, k, xik, yik, zik, rik, Scell%supce, 2, .true.)	! dr_{ikk}/ds_{k,y}
               drikk_dsk(3) = drij_dska(i, kk, k, xik, yik, zik, rik, Scell%supce, 3, .true.)	! dr_{ikk}/ds_{k,z}

               drjkk_dsk(1) = drij_dska(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 1, .true.)	! dr_{jkk}/ds_{k,x}
               drjkk_dsk(2) = drij_dska(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 2, .true.)	! dr_{jkk}/ds_{k,y}
               drjkk_dsk(3) = drij_dska(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 3, .true.)	! dr_{jkk}/ds_{k,z}


               ! Find the combination-of-atoms index:
               at_ind = find_3bdy_ind(KOA1, KOA2, KOA3)  ! module "Dealing_with_3TB"

               ! Factors including Laguerres are independent of shells:
               MM(1) = M_Lag_exp(i,kk,1) * M_Lag_exp(j,kk,1)
               MM(2) = M_Lag_exp(i,kk,1) * M_Lag_exp(j,kk,2)
               MM(3) = M_Lag_exp(i,kk,2) * M_Lag_exp(j,kk,1)
               MM(4) = MM(1) * M_Lag_exp(i,j,1)

               dMM(:,1) = M_d_Lag_exp(i,kk,1)*M_Lag_exp(j,kk,1)*drikk_dsk(:) + &
                          M_Lag_exp(i,kk,1)*M_d_Lag_exp(j,kk,1)*drjkk_dsk(:)
               dMM(:,2) = M_d_Lag_exp(i,kk,1)*M_Lag_exp(j,kk,2)*drikk_dsk(:) + &
                          M_Lag_exp(i,kk,1)*M_d_Lag_exp(j,kk,2)*drjkk_dsk(:)
               dMM(:,3) = M_d_Lag_exp(i,kk,2)*M_Lag_exp(j,kk,1)*drikk_dsk(:) + &
                          M_Lag_exp(i,kk,2)*M_d_Lag_exp(j,kk,1)*drjkk_dsk(:)
               dMM(:,4) = dMM(:,1)*M_Lag_exp(i,j,1) + MM(1)*M_d_Lag_exp(i,j,1)*drij_dsk(:)

               ! Get the radial function for 3-body interaction:
               do sh1 = 1, 1+basis_ind
                  do sh2 = 1, 1+basis_ind
                     G_IJK(sh1,sh2) = TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 1) * MM(1) + &
                                     TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 2) * MM(2) + &
                                     TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 3) * MM(3) + &
                                     TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 4) * MM(4)

                     dG_IJK(:,sh1,sh2) = TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 1) * dMM(:,1) + &
                                     TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 2) * dMM(:,2) + &
                                     TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 3) * dMM(:,3) + &
                                     TB(KOA1,KOA2)%V3bdy(at_ind, sh1, sh2, 4) * dMM(:,4)
                  enddo
               enddo
               ! Include angular parts:
               dH_3bdy(:,1,1) = dH_3bdy(:,1,1) + dG_IJK(:,1,1)   ! s-s * s-s / d s_{x,y,z}

               if (basis_ind > 0) then ! p3 orbitals:
                  ! Slater-Koster angular parts:
                  Miks(2) = Mjs_in(i,kk,2)  ! px-s
                  Miks(3) = Mjs_in(i,kk,3)  ! py-s
                  Miks(4) = Mjs_in(i,kk,4)  ! pz-s

                  Mjks(2) = Mjs_in(j,kk,2)  ! px-s
                  Mjks(3) = Mjs_in(j,kk,3)  ! py-s
                  Mjks(4) = Mjs_in(j,kk,4)  ! pz-s


                  M_dlmn_ik(1) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 1, 1, drikk_dsk(1))	! dl/dsx
                  M_dlmn_ik(2) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 1, 2, drikk_dsk(2))	! dl/dsy
                  M_dlmn_ik(3) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 1, 3, drikk_dsk(3))	! dl/dsz
                  M_dlmn_ik(4) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 2, 1, drikk_dsk(1))	! dm/dsx
                  M_dlmn_ik(5) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 2, 2, drikk_dsk(2))	! dm/dsy
                  M_dlmn_ik(6) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 2, 3, drikk_dsk(3))	! dm/dsz
                  M_dlmn_ik(7) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 3, 1, drikk_dsk(1))	! dn/dsx
                  M_dlmn_ik(8) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 3, 2, drikk_dsk(2))	! dn/dsy
                  M_dlmn_ik(9) = ddija_dskb_kd(i, kk, k, xik, yik, zik, rik, Scell%supce, 3, 3, drikk_dsk(3))	! dn/dsz

                  M_dlmn_jk(1) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 1, 1, drjkk_dsk(1))	! dl/dsx
                  M_dlmn_jk(2) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 1, 2, drjkk_dsk(2))	! dl/dsy
                  M_dlmn_jk(3) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 1, 3, drjkk_dsk(3))	! dl/dsz
                  M_dlmn_jk(4) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 2, 1, drjkk_dsk(1))	! dm/dsx
                  M_dlmn_jk(5) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 2, 2, drjkk_dsk(2))	! dm/dsy
                  M_dlmn_jk(6) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 2, 3, drjkk_dsk(3))	! dm/dsz
                  M_dlmn_jk(7) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 3, 1, drjkk_dsk(1))	! dn/dsx
                  M_dlmn_jk(8) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 3, 2, drjkk_dsk(2))	! dn/dsy
                  M_dlmn_jk(9) = ddija_dskb_kd(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 3, 3, drjkk_dsk(3))	! dn/dsz


                  dMis(1,1) = -M_dlmn_ik(1) ! px-s / ds_x
                  dMis(2,1) = -M_dlmn_ik(2) ! px-s / ds_y
                  dMis(3,1) = -M_dlmn_ik(3) ! px-s / ds_z
                  dMis(1,2) = -M_dlmn_ik(4) ! py-s / ds_x
                  dMis(2,2) = -M_dlmn_ik(5) ! py-s / ds_y
                  dMis(3,2) = -M_dlmn_ik(6) ! py-s / ds_z
                  dMis(1,3) = -M_dlmn_ik(7) ! pz-s / ds_x
                  dMis(2,3) = -M_dlmn_ik(8) ! pz-s / ds_y
                  dMis(3,3) = -M_dlmn_ik(9) ! pz-s / ds_z

                  dMjs(1,1) = -M_dlmn_jk(1) ! px-s / ds_x
                  dMjs(2,1) = -M_dlmn_jk(2) ! px-s / ds_y
                  dMjs(3,1) = -M_dlmn_jk(3) ! px-s / ds_z
                  dMjs(1,2) = -M_dlmn_jk(4) ! py-s / ds_x
                  dMjs(2,2) = -M_dlmn_jk(5) ! py-s / ds_y
                  dMjs(3,2) = -M_dlmn_jk(6) ! py-s / ds_z
                  dMjs(1,3) = -M_dlmn_jk(7) ! pz-s / ds_x
                  dMjs(2,3) = -M_dlmn_jk(8) ! pz-s / ds_y
                  dMjs(3,3) = -M_dlmn_jk(9) ! pz-s / ds_z


                  dH_3bdy(:,1,2) = dH_3bdy(:,1,2) + dG_IJK(:,2,1) * Mjks(2) + &
                                    G_IJK(2,1) * dMjs(:,1)   ! s-s * px-s / d s_{x,y,z}
                  dH_3bdy(:,1,3) = dH_3bdy(:,1,3) + dG_IJK(:,2,1) * Mjks(3) + &
                                    G_IJK(2,1) * dMjs(:,2)   ! s-s * py-s / d s_{x,y,z}
                  dH_3bdy(:,1,4) = dH_3bdy(:,1,4) + dG_IJK(:,2,1) * Mjks(4) + &
                                    G_IJK(2,1) * dMjs(:,3)   ! s-s * pz-s / d s_{x,y,z}

                  dH_3bdy(:,2,1) = dH_3bdy(:,2,1) + dG_IJK(:,2,1) * Miks(2) + &
                                    G_IJK(2,1) * dMis(:,1)   ! px-s * s-s / d s_{x,y,z}
                  dH_3bdy(:,3,1) = dH_3bdy(:,3,1) + dG_IJK(:,2,1) * Miks(3) + &
                                    G_IJK(2,1) * dMis(:,2)   ! py-s * s-s / d s_{x,y,z}
                  dH_3bdy(:,4,1) = dH_3bdy(:,4,1) + dG_IJK(:,2,1) * Miks(4) + &
                                    G_IJK(2,1) * dMis(:,3)   ! pz-s * s-s / d s_{x,y,z}


                  dH_3bdy(:,2,2) = dH_3bdy(:,2,2) + dG_IJK(:,2,2) * Miks(2) * Mjks(2) + &
                     G_IJK(2,2)*( dMis(:,1) * Mjks(2) + Miks(2) * dMjs(:,1) ) ! px-s px-s / d s_{x,y,z}
                  dH_3bdy(:,2,3) = dH_3bdy(:,2,3) + dG_IJK(:,2,2) * Miks(2) * Mjks(3) + &
                     G_IJK(2,2)*( dMis(:,1) * Mjks(3) + Miks(2) * dMjs(:,2) ) ! px-s py-s / d s_{x,y,z}
                  dH_3bdy(:,2,4) = dH_3bdy(:,2,4) + dG_IJK(:,2,2) * Miks(2) * Mjks(4) + &
                     G_IJK(2,2)*( dMis(:,1) * Mjks(4) + Miks(2) * dMjs(:,3) ) ! px-s pz-s / d s_{x,y,z}

                  dH_3bdy(:,3,2) = dH_3bdy(:,3,2) + dG_IJK(:,2,2) * Miks(3) * Mjks(2) + &
                     G_IJK(2,2)*( dMis(:,2) * Mjks(2) + Miks(3) * dMjs(:,1) ) ! py-s px-s / d s_{x,y,z}
                  dH_3bdy(:,4,2) = dH_3bdy(:,4,2) + dG_IJK(:,2,2) * Miks(4) * Mjks(2) + &
                     G_IJK(2,2)*( dMis(:,3) * Mjks(2) + Miks(4) * dMjs(:,1) ) ! pz-s px-s / d s_{x,y,z}


                  dH_3bdy(:,3,3) = dH_3bdy(:,3,3) + dG_IJK(:,2,2) * Miks(3) * Mjks(3) + &
                     G_IJK(2,2)*( dMis(:,2) * Mjks(3) + Miks(3) * dMjs(:,2) ) ! py-s py-s / d s_{x,y,z}
                  dH_3bdy(:,3,4) = dH_3bdy(:,3,4) + dG_IJK(:,2,2) * Miks(3) * Mjks(4) + &
                     G_IJK(2,2)*( dMis(:,2) * Mjks(4) + Miks(3) * dMjs(:,3) ) ! py-s pz-s / d s_{x,y,z}
                  dH_3bdy(:,4,3) = dH_3bdy(:,4,3) + dG_IJK(:,2,2) * Miks(4) * Mjks(3) + &
                     G_IJK(2,2)*( dMis(:,3) * Mjks(3) + Miks(4) * dMjs(:,2) ) ! pz-s py-s / d s_{x,y,z}

                  dH_3bdy(:,4,4) = dH_3bdy(:,4,4) + dG_IJK(:,2,2) * Miks(4) * Mjks(4) + &
                     G_IJK(2,2)*( dMis(:,3) * Mjks(4) + Miks(4) * dMjs(:,3) ) ! pz-s pz-s / d s_{x,y,z}

                  if (basis_ind > 1) then ! d5 orbitals:

                     ! Define d-symmetry factors to be reused below:
                     dMjs(1,4) = m_sqrt3 * (M_dlmn(1)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(4)) ! dxy-s / ds_x
                     dMjs(2,4) = m_sqrt3 * (M_dlmn(2)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(5)) ! dxy-s / ds_y
                     dMjs(3,4) = m_sqrt3 * (M_dlmn(3)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(6)) ! dxy-s / ds_z
                     dMjs(1,5) = m_sqrt3 * (M_dlmn(1)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(7)) ! dxz-s / ds_x
                     dMjs(2,5) = m_sqrt3 * (M_dlmn(2)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(8)) ! dxz-s / ds_y
                     dMjs(3,5) = m_sqrt3 * (M_dlmn(3)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(9)) ! dxz-s / ds_z
                     dMjs(1,6) = m_sqrt3 * (M_dlmn(4)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(7)) ! dyz-s / ds_x
                     dMjs(2,6) = m_sqrt3 * (M_dlmn(5)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(8)) ! dyz-s / ds_y
                     dMjs(3,6) = m_sqrt3 * (M_dlmn(6)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(9)) ! dyz-s / ds_z
                     dMjs(1,7) = m_sqrt3 * (M_dlmn(1)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(4)) ! (dx2-y2)-s / ds_x
                     dMjs(2,7) = m_sqrt3 * (M_dlmn(2)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(5)) ! (dx2-y2)-s / ds_y
                     dMjs(3,7) = m_sqrt3 * (M_dlmn(3)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(6)) ! (dx2-y2)-s / ds_z
                     dMjs(1,8) = 2.0d0*M_dlmn(7)*M_lmn(3,i,j) - (M_dlmn(1)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(4)) ! (d3z2-r2)-s / ds_x
                     dMjs(2,8) = 2.0d0*M_dlmn(8)*M_lmn(3,i,j) - (M_dlmn(2)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(5)) ! (d3z2-r2)-s / ds_y
                     dMjs(3,8) = 2.0d0*M_dlmn(9)*M_lmn(3,i,j) - (M_dlmn(3)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(6)) ! (d3z2-r2)-s / ds_z

                     Mjs(5) = Mjs_in(i,j,5)  ! dxy-s
                     Mjs(6) = Mjs_in(i,j,6)  ! dxz-s
                     Mjs(7) = Mjs_in(i,j,7)  ! dyz-s
                     Mjs(8) = Mjs_in(i,j,8)  ! (dx2-y2)-s
                     Mjs(9) = Mjs_in(i,j,9)  ! (d3z2-r2)-s


!                      dH_3bdy(1,5) = dH_3bdy(1,5) + G_IJK(1,3) * Mjs(j,k,5) ! s-s * dxy-s
!                      dH_3bdy(1,6) = dH_3bdy(1,6) + G_IJK(1,3) * Mjs(j,k,6) ! s-s * dxz-s
!                      dH_3bdy(1,7) = dH_3bdy(1,7) + G_IJK(1,3) * Mjs(j,k,7) ! s-s * dxz-s
!                      dH_3bdy(1,8) = dH_3bdy(1,8) + G_IJK(1,3) * Mjs(j,k,8) ! s-s * (dx2-y2)-s
!                      dH_3bdy(1,9) = dH_3bdy(1,9) + G_IJK(1,3) * Mjs(j,k,9) ! s-s * (d3z2-r2)-s
!
!                      dH_3bdy(5,1) = dH_3bdy(5,1) + G_IJK(3,1) * Mjs(i,k,5) ! dxy-s * s-s
!                      dH_3bdy(6,1) = dH_3bdy(6,1) + G_IJK(3,1) * Mjs(i,k,6) ! dxz-s * s-s
!                      dH_3bdy(7,1) = dH_3bdy(7,1) + G_IJK(3,1) * Mjs(i,k,7) ! dxz-s * s-s
!                      dH_3bdy(8,1) = dH_3bdy(8,1) + G_IJK(3,1) * Mjs(i,k,8) ! (dx2-y2)-s * s-s
!                      dH_3bdy(9,1) = dH_3bdy(9,1) + G_IJK(3,1) * Mjs(i,k,9) ! (d3z2-r2)-s * s-s

                     ! Calculate repeating part the K-S matrix elements:
!                      H_temp(1) = Mjs(i,k,2) * Mjs(j,k,5) ! px-s * dxy-s
!                      H_temp(2) = Mjs(i,k,2) * Mjs(j,k,6) ! px-s * dxz-s
!                      H_temp(3) = Mjs(i,k,2) * Mjs(j,k,7) ! px-s * dxz-s
!                      H_temp(4) = Mjs(i,k,2) * Mjs(j,k,8) ! px-s * (dx2-y2)-s
!                      H_temp(5) = Mjs(i,k,2) * Mjs(j,k,9) ! px-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(2,5:9) = dH_3bdy(2,5:9) + H_temp(1:5) * G_IJK(2,3)
!                      dH_3bdy(5:9,2) = dH_3bdy(5:9,2) + H_temp(1:5) * G_IJK(3,2)

                     ! Calculate repeating part the K-S matrix elements:
!                      H_temp(1) = Mjs(i,k,3) * Mjs(j,k,5) ! py-s * dxy-s
!                      H_temp(2) = Mjs(i,k,3) * Mjs(j,k,6) ! py-s * dxz-s
!                      H_temp(3) = Mjs(i,k,3) * Mjs(j,k,7) ! py-s * dxz-s
!                      H_temp(4) = Mjs(i,k,3) * Mjs(j,k,8) ! py-s * (dx2-y2)-s
!                      H_temp(5) = Mjs(i,k,3) * Mjs(j,k,9) ! py-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(3,5:9) = dH_3bdy(3,5:9) + H_temp(1:5) * G_IJK(2,3)
!                      dH_3bdy(5:9,3) = dH_3bdy(5:9,3) + H_temp(1:5) * G_IJK(3,2)

                     ! Calculate repeating part the K-S matrix elements:
!                      H_temp(1) = Mjs(i,k,4) * Mjs(j,k,5) ! dxy-s * pz-s
!                      H_temp(2) = Mjs(i,k,4) * Mjs(j,k,6) ! dxz-s * pz-s
!                      H_temp(3) = Mjs(i,k,4) * Mjs(j,k,7) ! dyz-s * pz-s
!                      H_temp(4) = Mjs(i,k,4) * Mjs(j,k,8) ! (dx2-y2)-s * pz-s
!                      H_temp(5) = Mjs(i,k,4) * Mjs(j,k,9) ! (d3z2-r2)-s * pz-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(4,5:9) = dH_3bdy(4,5:9) + H_temp(1:5) * G_IJK(2,3)
!                      dH_3bdy(5:9,4) = dH_3bdy(5:9,4) + H_temp(1:5) * G_IJK(3,2)

                     ! Calculate the K-S matrix elements:
!                      H_temp(6) = Mjs(i,k,5) * G_IJK(3,3)
!                      H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dxy-s * dxy-s
!                      H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dxy-s * dxz-s
!                      H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dxy-s * dyz-s
!                      H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dxy-s * (dx2-y2)-s
!                      H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dxy-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(5,5:9) = dH_3bdy(5,5:9) + H_temp(1:5)

                     ! Calculate the K-S matrix elements:
!                      H_temp(6) = Mjs(i,k,6) * G_IJK(3,3)
!                      H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dxz-s * dxy-s
!                      H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dxz-s * dxz-s
!                      H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dxz-s * dyz-s
!                      H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dxz-s * (dx2-y2)-s
!                      H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dxz-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(6,5:9) = dH_3bdy(6,5:9) + H_temp(1:5)

                     ! Calculate the K-S matrix elements:
!                      H_temp(6) = Mjs(i,k,7) * G_IJK(3,3)
!                      H_temp(1) = H_temp(6) * Mjs(j,k,5) ! dyz-s * dxy-s
!                      H_temp(2) = H_temp(6) * Mjs(j,k,6) ! dyz-s * dxz-s
!                      H_temp(3) = H_temp(6) * Mjs(j,k,7) ! dyz-s * dyz-s
!                      H_temp(4) = H_temp(6) * Mjs(j,k,8) ! dyz-s * (dx2-y2)-s
!                      H_temp(5) = H_temp(6) * Mjs(j,k,9) ! dyz-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(7,5:9) = dH_3bdy(7,5:9) + H_temp(1:5)

                     ! Calculate the K-S matrix elements:
!                      H_temp(6) = Mjs(i,k,8) * G_IJK(3,3)
!                      H_temp(1) = H_temp(6) * Mjs(j,k,5) ! (dx2-y2)-s * dxy-s
!                      H_temp(2) = H_temp(6) * Mjs(j,k,6) ! (dx2-y2)-s * dxz-s
!                      H_temp(3) = H_temp(6) * Mjs(j,k,7) ! (dx2-y2)-s * dyz-s
!                      H_temp(4) = H_temp(6) * Mjs(j,k,8) ! (dx2-y2)-s * (dx2-y2)-s
!                      H_temp(5) = H_temp(6) * Mjs(j,k,9) ! (dx2-y2)-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(8,5:9) = dH_3bdy(8,5:9) + H_temp(1:5)

                     ! Calculate the K-S matrix elements:
!                      H_temp(6) = Mjs(i,k,9) * G_IJK(3,3)
!                      H_temp(1) = H_temp(6) * Mjs(j,k,5) ! (d3z2-r2)-s * dxy-s
!                      H_temp(2) = H_temp(6) * Mjs(j,k,6) ! (d3z2-r2)-s * dxz-s
!                      H_temp(3) = H_temp(6) * Mjs(j,k,7) ! (d3z2-r2)-s * dyz-s
!                      H_temp(4) = H_temp(6) * Mjs(j,k,8) ! (d3z2-r2)-s * (dx2-y2)-s
!                      H_temp(5) = H_temp(6) * Mjs(j,k,9) ! (d3z2-r2)-s * (d3z2-r2)-s
                     ! Add it into the Hamiltonian part:
!                      dH_3bdy(9,5:9) = dH_3bdy(9,5:9) + H_temp(1:5)

                  endif ! (basis_ind > 1)
               endif ! (basis_ind > 0)

            endif ! (kk /= j)
         enddo AT3
         ! add the contributions into the output:
         dH = dH + dH_3bdy

         deallocate(dH_3bdy)
      endif ! (TB(KOA1,KOA1)%include_3body)

   endif ! (i == j)


   ! Clean up at the end:
   nullify (x1, y1, z1, r1, KOA1, KOA2, KOA3, m)
   nullify(xij, yij, zij, rij)
   nullify(xik, yik, zik, rik)
   nullify(xjk, yjk, zjk, rjk)
end subroutine d_Hamilton_one_3TB





subroutine d_Onsite_3TB(basis_ind, k, i, Scell, TB, M_lmn, M_Lag_exp, M_d_Lag_exp, Mjs_in, dHij)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   integer, intent(in) :: k, i         ! indix of the atom#1, and atoms #2=#3
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   real(8), dimension(:,:,:), intent(in) :: M_lmn  ! matrix of l, m, n
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: M_d_Lag_exp   ! matrix of derivatives of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: Mjs_in ! matrix of overlaps with s-orbital
   real(8), dimension(:,:,:), intent(out) :: dHij	! derivatives of the overlap integerals along x,y,z [eV/A]
   !---------------------
   integer :: i1, Bsiz, atom_2, sh1, sh2, j, at_ind, atom_3, kk
   real(8) :: r, term_s, term_p, term_d, term_3bdy(3)
   real(8) :: matr_spd(3,3), d_matr_spd(3,3), H_cf_temp(9), drij_dsk(3), drikk_dsk(3), drjkk_dsk(3)
   real(8), dimension(9) :: M_dlmn	! dl/dx, dl/dy, dl/dz, dm/dx, dm/dy, dm/dz, dn/dx, dn/dy, dn/dz
   real(8), dimension(3, 8) :: dMjs ! 3p, 5d
   real(8), dimension(9) :: Mjs ! 1s, 3p, 5d
   real(8), dimension(:,:,:), allocatable :: E_onsite, H_avg, H_cf, H_3bdy
   integer, pointer :: m, KOA1, KOA2, KOA3
   real(8), pointer :: r1, x1, y1, z1
   real(8), pointer :: rij, xij, yij, zij
   real(8), pointer :: rik, xik, yik, zik
   real(8), pointer :: rjk, xjk, yjk, zjk



   ! Kind of atom:
   KOA1 => Scell%MDatoms(i)%KOA
   ! Number of the nearest neighbors
   m => Scell%Near_neighbor_size(i)

   ! Size of the basis:
   select case (basis_ind)
   case (0) ! s
      Bsiz = 1
   case (1) ! sp3
      Bsiz = 4
   case(2)  ! sp3d5
      Bsiz = 9
   end select

   ! To start with:
   allocate(E_onsite(3,Bsiz,Bsiz), source=0.0d0)
   allocate(H_avg(3,Bsiz,Bsiz), source=0.0d0)
   allocate(H_cf(3,Bsiz,Bsiz), source=0.0d0)
   allocate(H_3bdy(3,Bsiz,Bsiz), source=0.0d0)

   ! The onsite energies are constructed out of 4 terms [1]:
   ! 1) onsite eigenvalues
   ! 2) average contribution of atoms around
   ! 3) crystall field contribution
   ! 4) 3-body contributions
   ! Let's calculate them all:


   !-----------------
   ! 1) Onsite energies of spin-unpolirized orbital values:
   E_onsite = 0.0d0     ! constants -> all derivatives = 0

   !-----------------
   ! 2,3) Average and crystal field terms (tested, correct):
   do atom_2 = 1, m ! do only for atoms close to that one
      j = Scell%Near_neighbor_list(i, atom_2) ! this is the list of such close atoms
      ! [OS 0]
!       KOA1 => Scell%MDatoms(i)%KOA
!       KOA2 => Scell%MDatoms(j)%KOA
      ! [OS 1]
      KOA1 => Scell%MDatoms(j)%KOA
      KOA2 => Scell%MDatoms(i)%KOA


      x1 => Scell%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
      y1 => Scell%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
      z1 => Scell%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
      r1 => Scell%Near_neighbor_dist(i,atom_2,4)	! at this distance, R

      ! Derivatives of rij by sk:
      ! All functions are from the module "TB_Koster_Slater":
      drij_dsk(1) = drij_dska(i, j, k, x1, y1, z1, r1, Scell%supce, 1, .true.)	! dr_{ij}/ds_{k,x}
      drij_dsk(2) = drij_dska(i, j, k, x1, y1, z1, r1, Scell%supce, 2, .true.)	! dr_{ij}/ds_{k,y}
      drij_dsk(3) = drij_dska(i, j, k, x1, y1, z1, r1, Scell%supce, 3, .true.)	! dr_{ij}/ds_{k,z}


      ! 2) Average terms:
      term_s = SUM( TB(KOA1,KOA2)%Hhavg(1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s orbitals
      H_avg(:,1,1) = H_avg(:,1,1) + term_s * drij_dsk(:)

      if (basis_ind > 0) then ! p3 orbitals:
         term_p = SUM( TB(KOA1,KOA2)%Hhavg(2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p3 orbitals
         do i1 = 2, 4
            H_avg(:,i1,i1) = H_avg(:,i1,i1) + term_p * drij_dsk(:)
         enddo

         if (basis_ind > 1) then ! d5 orbitals:
            term_d = SUM( TB(KOA1,KOA2)%Hhavg(3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d5 orbitals
            do i1 = 5, 9
               H_avg(:,i1,i1) = H_avg(:,i1,i1) + term_d * drij_dsk(:)
            enddo
         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)

      !-----------------
      ! 3) Crystal field:
      ! [OS 0]
!       KOA1 => Scell%MDatoms(i)%KOA
!       KOA2 => Scell%MDatoms(j)%KOA
      ! [OS 1] tested, correct
      KOA1 => Scell%MDatoms(j)%KOA
      KOA2 => Scell%MDatoms(i)%KOA


      M_dlmn(1) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 1, drij_dsk(1))	! dl/dsx
      M_dlmn(2) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 2, drij_dsk(2))	! dl/dsy
      M_dlmn(3) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 1, 3, drij_dsk(3))	! dl/dsz
      M_dlmn(4) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 1, drij_dsk(1))	! dm/dsx
      M_dlmn(5) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 2, drij_dsk(2))	! dm/dsy
      M_dlmn(6) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 2, 3, drij_dsk(3))	! dm/dsz
      M_dlmn(7) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 1, drij_dsk(1))	! dn/dsx
      M_dlmn(8) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 2, drij_dsk(2))	! dn/dsy
      M_dlmn(9) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell%supce, 3, 3, drij_dsk(3))	! dn/dsz

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
      ! and its derivatives:
      d_matr_spd(1,1) = SUM( TB(KOA1,KOA2)%Hhcf(1,1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s-s orbitals
      if (basis_ind > 0) then ! p3 orbitals:
         d_matr_spd(1,2) = SUM( TB(KOA1,KOA2)%Hhcf(1,2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s-p orbitals
         d_matr_spd(2,1) = SUM( TB(KOA1,KOA2)%Hhcf(2,1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p-s orbitals
         d_matr_spd(2,2) = SUM( TB(KOA1,KOA2)%Hhcf(2,2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p-p orbitals
         if (basis_ind > 1) then ! d5 orbitals:
            d_matr_spd(1,3) = SUM( TB(KOA1,KOA2)%Hhcf(1,3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s-d orbitals
            d_matr_spd(3,1) = SUM( TB(KOA1,KOA2)%Hhcf(3,1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d-s orbitals
            d_matr_spd(2,3) = SUM( TB(KOA1,KOA2)%Hhcf(2,3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p-d orbitals
            d_matr_spd(3,2) = SUM( TB(KOA1,KOA2)%Hhcf(3,2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d-p orbitals
            d_matr_spd(3,3) = SUM( TB(KOA1,KOA2)%Hhcf(3,3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d-d orbitals
         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)

      ! Include angular parts:
      ! Mjs(1) = Mjs_in(i,j,1) ! s-s  = 1.0
      H_cf(:,1,1) = H_cf(:,1,1) + d_matr_spd(1,1) * drij_dsk(:)  ! s-s * s-s

      if (basis_ind > 0) then ! p3 orbitals:

         ! Define derivatives of the p-symmetry factors to be reused below:
         ! [F 1]
         dMjs(1,1) = -M_dlmn(1) ! px-s / ds_x
         dMjs(2,1) = -M_dlmn(2) ! px-s / ds_y
         dMjs(3,1) = -M_dlmn(3) ! px-s / ds_z
         dMjs(1,2) = -M_dlmn(4) ! py-s / ds_x
         dMjs(2,2) = -M_dlmn(5) ! py-s / ds_y
         dMjs(3,2) = -M_dlmn(6) ! py-s / ds_z
         dMjs(1,3) = -M_dlmn(7) ! pz-s / ds_x
         dMjs(2,3) = -M_dlmn(8) ! pz-s / ds_y
         dMjs(3,3) = -M_dlmn(9) ! pz-s / ds_z

         ! [E 2], tested
         Mjs(2) = Mjs_in(i,j,2)  ! px-s
         Mjs(3) = Mjs_in(i,j,3)  ! py-s
         Mjs(4) = Mjs_in(i,j,4)  ! pz-s
!          ! [E 1]
!          Mjs(2) = -Mjs_in(i,j,2)  ! px-s
!          Mjs(3) = -Mjs_in(i,j,3)  ! py-s
!          Mjs(4) = -Mjs_in(i,j,4)  ! pz-s


         H_cf_temp(1) = d_matr_spd(2,1) * Mjs(2)
         H_cf(1,1,2) = H_cf(1,1,2) + H_cf_temp(1) * drij_dsk(1) + matr_spd(2,1) * dMjs(1,1)  ! s-s * px-s / ds_x
         H_cf(2,1,2) = H_cf(2,1,2) + H_cf_temp(1) * drij_dsk(2) + matr_spd(2,1) * dMjs(2,1)  ! s-s * px-s / ds_y
         H_cf(3,1,2) = H_cf(3,1,2) + H_cf_temp(1) * drij_dsk(3) + matr_spd(2,1) * dMjs(3,1)  ! s-s * px-s / ds_z
         ! Lower triangle use symmetry:
         H_cf(:,2,1) = H_cf(:,1,2)

         H_cf_temp(1) = d_matr_spd(2,1) * Mjs(3)
         H_cf(1,1,3) = H_cf(1,1,3) + H_cf_temp(1) * drij_dsk(1) + matr_spd(2,1) * dMjs(1,2)  ! s-s * py-s / ds_x
         H_cf(2,1,3) = H_cf(2,1,3) + H_cf_temp(1) * drij_dsk(2) + matr_spd(2,1) * dMjs(2,2)  ! s-s * py-s / ds_y
         H_cf(3,1,3) = H_cf(3,1,3) + H_cf_temp(1) * drij_dsk(3) + matr_spd(2,1) * dMjs(3,2)  ! s-s * py-s / ds_z
         ! Lower triangle use symmetry:
         H_cf(:,3,1) = H_cf(:,1,3)

         H_cf_temp(1) = d_matr_spd(2,1) * Mjs(4)
         H_cf(1,1,4) = H_cf(1,1,4) + H_cf_temp(1) * drij_dsk(1) + matr_spd(2,1) * dMjs(1,3)  ! s-s * pz-s / ds_x
         H_cf(2,1,4) = H_cf(2,1,4) + H_cf_temp(1) * drij_dsk(2) + matr_spd(2,1) * dMjs(2,3)  ! s-s * pz-s / ds_y
         H_cf(3,1,4) = H_cf(3,1,4) + H_cf_temp(1) * drij_dsk(3) + matr_spd(2,1) * dMjs(3,3)  ! s-s * pz-s / ds_z
         ! Lower triangle use symmetry:
         H_cf(:,4,1) = H_cf(:,1,4)

         ! Diagonal part is excluded (no self-interaction of orbitals):
         if (TB(KOA1,KOA1)%nullify_diag_cf) then
            H_cf(:,2,2) = 0.0d0
            H_cf(:,3,3) = 0.0d0
            H_cf(:,4,4) = 0.0d0
         else
            H_cf_temp(1) = d_matr_spd(2,2) * Mjs(2) * Mjs(2)
            H_cf(:,2,2) = H_cf(:,2,2) + H_cf_temp(1) * drij_dsk(:) + &
                        matr_spd(2,2) * 2.0d0* dMjs(:,1) * Mjs(2) ! px-s * px-s /ds_x
            H_cf_temp(1) = d_matr_spd(2,2) * Mjs(3) * Mjs(3)
            H_cf(:,3,3) = H_cf(:,3,3) + H_cf_temp(1) * drij_dsk(:) + &
                        matr_spd(2,2) * 2.0d0* dMjs(:,2) * Mjs(3) ! py-s * py-s /ds_x
            H_cf_temp(1) = d_matr_spd(2,2) * Mjs(4) * Mjs(4)
            H_cf(:,4,4) = H_cf(:,4,4) + H_cf_temp(1) * drij_dsk(:) + &
                        matr_spd(2,2) * 2.0d0* dMjs(:,3) * Mjs(4) ! pz-s * pz-s /ds_x
         endif

         ! Off-diagonal part:
         H_cf_temp(1) = d_matr_spd(2,2) * Mjs(2) * Mjs(3)
         H_cf(1,2,3) = H_cf(1,2,3) + H_cf_temp(1) * drij_dsk(1) + &
                        matr_spd(2,2) * (dMjs(1,1) * Mjs(3) + Mjs(2) * dMjs(1,2) ) ! px-s * py-s /ds_x
         H_cf(2,2,3) = H_cf(2,2,3) + H_cf_temp(1) * drij_dsk(2) + &
                        matr_spd(2,2) * (dMjs(2,1) * Mjs(3) + Mjs(2) * dMjs(2,2) ) ! px-s * py-s /ds_y
         H_cf(3,2,3) = H_cf(3,2,3) + H_cf_temp(1) * drij_dsk(3) + &
                        matr_spd(2,2) * (dMjs(3,1) * Mjs(3) + Mjs(2) * dMjs(3,2) ) ! px-s * py-s /ds_z
         ! For lower triangle use symmetry:
         H_cf(:,3,2) = H_cf(:,2,3) ! py-s * px-s /ds_{x,y,z}

         H_cf_temp(1) = d_matr_spd(2,2) * Mjs(2) * Mjs(4)
         H_cf(1,2,4) = H_cf(1,2,4) + H_cf_temp(1) * drij_dsk(1) + &
                        matr_spd(2,2) * (dMjs(1,1) * Mjs(4) + Mjs(2) * dMjs(1,3) ) ! px-s * pz-s /ds_x
         H_cf(2,2,4) = H_cf(2,2,4) + H_cf_temp(1) * drij_dsk(2) + &
                        matr_spd(2,2) * (dMjs(2,1) * Mjs(4) + Mjs(2) * dMjs(2,3) ) ! px-s * pz-s /ds_y
         H_cf(3,2,4) = H_cf(3,2,4) + H_cf_temp(1) * drij_dsk(3) + &
                        matr_spd(2,2) * (dMjs(3,1) * Mjs(4) + Mjs(2) * dMjs(3,3) ) ! px-s * pz-s /ds_z
         ! For lower triangle use symmetry:
         H_cf(:,4,2) = H_cf(:,2,4) ! pz-s * px-s /ds_{x,y,z}


         H_cf_temp(1) = d_matr_spd(2,2) * Mjs(3) * Mjs(4)
         H_cf(1,3,4) = H_cf(1,3,4) + H_cf_temp(1) * drij_dsk(1) + &
                        matr_spd(2,2) * (dMjs(1,2) * Mjs(4) + Mjs(3) * dMjs(1,3) ) ! py-s * pz-s /ds_x
         H_cf(2,3,4) = H_cf(2,3,4) + H_cf_temp(1) * drij_dsk(2) + &
                        matr_spd(2,2) * (dMjs(2,2) * Mjs(4) + Mjs(3) * dMjs(2,3) ) ! py-s * pz-s /ds_y
         H_cf(3,3,4) = H_cf(3,3,4) + H_cf_temp(1) * drij_dsk(3) + &
                        matr_spd(2,2) * (dMjs(3,2) * Mjs(4) + Mjs(3) * dMjs(3,3) ) ! py-s * pz-s /ds_z
         ! For lower triangle use symmetry:
         H_cf(:,4,3) = H_cf(:,3,4)  ! pz-s * py-s /ds_{x,y,z}


         if (basis_ind > 1) then ! d5 orbitals:

            ! Define d-symmetry factors to be reused below:
            dMjs(1,4) = m_sqrt3 * (M_dlmn(1)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(4)) ! dxy-s / ds_x
            dMjs(2,4) = m_sqrt3 * (M_dlmn(2)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(5)) ! dxy-s / ds_y
            dMjs(3,4) = m_sqrt3 * (M_dlmn(3)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(6)) ! dxy-s / ds_z
            dMjs(1,5) = m_sqrt3 * (M_dlmn(1)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(7)) ! dxz-s / ds_x
            dMjs(2,5) = m_sqrt3 * (M_dlmn(2)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(8)) ! dxz-s / ds_y
            dMjs(3,5) = m_sqrt3 * (M_dlmn(3)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(9)) ! dxz-s / ds_z
            dMjs(1,6) = m_sqrt3 * (M_dlmn(4)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(7)) ! dyz-s / ds_x
            dMjs(2,6) = m_sqrt3 * (M_dlmn(5)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(8)) ! dyz-s / ds_y
            dMjs(3,6) = m_sqrt3 * (M_dlmn(6)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(9)) ! dyz-s / ds_z
            dMjs(1,7) = m_sqrt3 * (M_dlmn(1)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(4)) ! (dx2-y2)-s / ds_x
            dMjs(2,7) = m_sqrt3 * (M_dlmn(2)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(5)) ! (dx2-y2)-s / ds_y
            dMjs(3,7) = m_sqrt3 * (M_dlmn(3)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(6)) ! (dx2-y2)-s / ds_z
            dMjs(1,8) = 2.0d0*M_dlmn(7)*M_lmn(3,i,j) - (M_dlmn(1)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(4)) ! (d3z2-r2)-s / ds_x
            dMjs(2,8) = 2.0d0*M_dlmn(8)*M_lmn(3,i,j) - (M_dlmn(2)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(5)) ! (d3z2-r2)-s / ds_y
            dMjs(3,8) = 2.0d0*M_dlmn(9)*M_lmn(3,i,j) - (M_dlmn(3)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(6)) ! (d3z2-r2)-s / ds_z

            Mjs(5) = Mjs_in(i,j,5)  ! dxy-s
            Mjs(6) = Mjs_in(i,j,6)  ! dxz-s
            Mjs(7) = Mjs_in(i,j,7)  ! dyz-s
            Mjs(8) = Mjs_in(i,j,8)  ! (dx2-y2)-s
            Mjs(9) = Mjs_in(i,j,9)  ! (d3z2-r2)-s


            ! Diagonal terms:
            if (TB(KOA1,KOA1)%nullify_diag_cf) then
               H_cf(:,5,5) = 0.0d0  ! dxy-s * dxy-s
               H_cf(:,6,6) = 0.0d0  ! dxz-s * dxz-s
               H_cf(:,7,7) = 0.0d0  ! dyz-s * dyz-s
               H_cf(:,8,8) = 0.0d0  ! (dx2-y2) * (dx2-y2)
               H_cf(:,9,9) = 0.0d0  ! (d3z2-r2) * (d3z2-r2)
            else
               H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5)
               H_cf(:,5,5) = H_cf(:,5,5) + H_cf_temp(1) * drij_dsk(:) + &
                              matr_spd(3,3) * 2.0d0 * Mjs(5) * dMjs(:,4)  ! dxy-s * dxy-s / ds_{x,y,z}
               H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6)
               H_cf(:,6,6) = H_cf(:,6,6) + H_cf_temp(1) * drij_dsk(:) + &
                              matr_spd(3,3) * 2.0d0 * Mjs(6) * dMjs(:,5)  ! dxz-s * dxz-s / ds_{x,y,z}
               H_cf_temp(1) = d_matr_spd(3,3) * Mjs(7)
               H_cf(:,7,7) = H_cf(:,7,7) + H_cf_temp(1) * drij_dsk(:) + &
                              matr_spd(3,3) * 2.0d0 * Mjs(7) * dMjs(:,6)  ! dyz-s * dyz-s / ds_{x,y,z}
               H_cf_temp(1) = d_matr_spd(3,3) * Mjs(8)
               H_cf(:,8,8) = H_cf(:,8,8) + H_cf_temp(1) * drij_dsk(:) + &
                              matr_spd(3,3) * 2.0d0 * Mjs(8) * dMjs(:,7)  ! (dx2-y2) * (dx2-y2) / ds_{x,y,z}
               H_cf_temp(1) = d_matr_spd(3,3) * Mjs(9)
               H_cf(:,9,9) = H_cf(:,9,9) + H_cf_temp(1) * drij_dsk(:) + &
                              matr_spd(3,3) * 2.0d0 * Mjs(9) * dMjs(:,8)  ! (d3z2-r2) * (d3z2-r2) / ds_{x,y,z}
            endif

            ! Off-diagonal terms:
            H_cf_temp(1) = d_matr_spd(1,3) * Mjs(5)
            H_cf(1,1,5) = H_cf(1,1,5) + H_cf_temp(1) * drij_dsk(1) + matr_spd(1,3) * dMjs(1,4)  ! s-s * dxy-s / ds_x
            H_cf(2,1,5) = H_cf(2,1,5) + H_cf_temp(1) * drij_dsk(2) + matr_spd(1,3) * dMjs(2,4)  ! s-s * dxy-s / ds_y
            H_cf(3,1,5) = H_cf(3,1,5) + H_cf_temp(1) * drij_dsk(3) + matr_spd(1,3) * dMjs(3,4)  ! s-s * dxy-s / ds_z
            ! For lower triangle use symmetry:
            H_cf(:,5,1) = H_cf(:,1,5)  ! dxy-s * s-s/ ds_{x,y,z}

            H_cf_temp(1) = d_matr_spd(1,3) * Mjs(6)
            H_cf(1,1,6) = H_cf(1,1,6) + H_cf_temp(1) * drij_dsk(1) + matr_spd(1,3) * dMjs(1,5)  ! s-s * dxz-s / ds_x
            H_cf(2,1,6) = H_cf(2,1,6) + H_cf_temp(1) * drij_dsk(2) + matr_spd(1,3) * dMjs(2,5)  ! s-s * dxz-s / ds_y
            H_cf(3,1,6) = H_cf(3,1,6) + H_cf_temp(1) * drij_dsk(3) + matr_spd(1,3) * dMjs(3,5)  ! s-s * dxz-s / ds_z
            ! For lower triangle use symmetry:
            H_cf(:,6,1) = H_cf(:,1,6)  ! dxz-s * s-s/ ds_{x,y,z}

            H_cf_temp(1) = d_matr_spd(1,3) * Mjs(7)
            H_cf(1,1,7) = H_cf(1,1,7) + H_cf_temp(1) * drij_dsk(1) + matr_spd(1,3) * dMjs(1,6)  ! s-s * dyz-s / ds_x
            H_cf(2,1,7) = H_cf(2,1,7) + H_cf_temp(1) * drij_dsk(2) + matr_spd(1,3) * dMjs(2,6)  ! s-s * dyz-s / ds_y
            H_cf(3,1,7) = H_cf(3,1,7) + H_cf_temp(1) * drij_dsk(3) + matr_spd(1,3) * dMjs(3,6)  ! s-s * dyz-s / ds_z
            ! For lower triangle use symmetry:
            H_cf(:,7,1) = H_cf(:,1,7)  ! dyz-s * s-s/ ds_{x,y,z}

            H_cf_temp(1) = d_matr_spd(1,3) * Mjs(8)
            H_cf(1,1,8) = H_cf(1,1,8) + H_cf_temp(1) * drij_dsk(1) + matr_spd(1,3) * dMjs(1,7)  ! s-s * (dx2-y2)-s / ds_x
            H_cf(2,1,8) = H_cf(2,1,8) + H_cf_temp(1) * drij_dsk(2) + matr_spd(1,3) * dMjs(2,7)  ! s-s * (dx2-y2)-s / ds_y
            H_cf(3,1,8) = H_cf(3,1,8) + H_cf_temp(1) * drij_dsk(3) + matr_spd(1,3) * dMjs(3,7)  ! s-s * (dx2-y2)-s / ds_z
            ! For lower triangle use symmetry:
            H_cf(:,8,1) = H_cf(:,1,8)  ! (dx2-y2)-s * s-s / ds_{x,y,z}

            H_cf_temp(1) = d_matr_spd(1,3) * Mjs(9)
            H_cf(1,1,9) = H_cf(1,1,9) + H_cf_temp(1) * drij_dsk(1) + matr_spd(1,3) * dMjs(1,8)  ! s-s * (d3z2-r2)-s / ds_x
            H_cf(2,1,9) = H_cf(2,1,9) + H_cf_temp(1) * drij_dsk(2) + matr_spd(1,3) * dMjs(2,8)  ! s-s * (d3z2-r2)-s / ds_y
            H_cf(3,1,9) = H_cf(3,1,9) + H_cf_temp(1) * drij_dsk(3) + matr_spd(1,3) * dMjs(3,8)  ! s-s * (d3z2-r2)-s / ds_z
            ! For lower triangle use symmetry:
            H_cf(:,9,1) = H_cf(:,1,9)  ! (d3z2-r2)-s * s-s / ds_{x,y,z}


            ! Calculate repeating part the K-S matrix elements:
            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(5)
            H_cf(1,2,5) = H_cf(1,2,5) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,1) * Mjs(5) + Mjs(2) * dMjs(1,4) ) ! px-s * dxy-s / ds_x
            H_cf(2,2,5) = H_cf(2,2,5) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,1) * Mjs(5) + Mjs(2) * dMjs(2,4) ) ! px-s * dxy-s / ds_y
            H_cf(3,2,5) = H_cf(3,2,5) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,1) * Mjs(5) + Mjs(2) * dMjs(3,4) ) ! px-s * dxy-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,5,2) = H_cf(:,2,5)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(6)
            H_cf(1,2,6) = H_cf(1,2,6) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,1) * Mjs(6) + Mjs(2) * dMjs(1,5) ) ! px-s * dxz-s / ds_x
            H_cf(2,2,6) = H_cf(2,2,6) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,1) * Mjs(6) + Mjs(2) * dMjs(2,5) ) ! px-s * dxz-s / ds_y
            H_cf(3,2,6) = H_cf(3,2,6) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,1) * Mjs(6) + Mjs(2) * dMjs(3,5) ) ! px-s * dxz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,6,2) = H_cf(:,2,6)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(7)
            H_cf(1,2,7) = H_cf(1,2,7) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,1) * Mjs(7) + Mjs(2) * dMjs(1,6) ) ! px-s * dyz-s / ds_x
            H_cf(2,2,7) = H_cf(2,2,7) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,1) * Mjs(7) + Mjs(2) * dMjs(2,6) ) ! px-s * dyz-s / ds_y
            H_cf(3,2,7) = H_cf(3,2,7) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(2,1) * Mjs(7) + Mjs(2) * dMjs(3,6) ) ! px-s * dyz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,7,2) = H_cf(:,2,7)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(8)
            H_cf(1,2,8) = H_cf(1,2,8) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,1) * Mjs(8) + Mjs(2) * dMjs(1,7) ) ! px-s * (dx2-y2)-s / ds_x
            H_cf(2,2,8) = H_cf(2,2,8) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,1) * Mjs(8) + Mjs(2) * dMjs(2,7) ) ! px-s * (dx2-y2)-s / ds_y
            H_cf(3,2,8) = H_cf(3,2,8) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,1) * Mjs(8) + Mjs(2) * dMjs(3,7) ) ! px-s * (dx2-y2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,8,2) = H_cf(:,2,8)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(9)
            H_cf(1,2,9) = H_cf(1,2,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,1) * Mjs(9) + Mjs(2) * dMjs(1,8) ) ! px-s * (d3z2-r2)-s / ds_x
            H_cf(2,2,9) = H_cf(2,2,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,1) * Mjs(9) + Mjs(2) * dMjs(2,8) ) ! px-s * (d3z2-r2)-s / ds_y
            H_cf(3,2,9) = H_cf(3,2,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,1) * Mjs(9) + Mjs(2) * dMjs(3,8) ) ! px-s * (d3z2-r2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,2) = H_cf(:,2,9)


            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(5)
            H_cf(1,3,5) = H_cf(1,3,5) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,2) * Mjs(5) + Mjs(3) * dMjs(1,4) ) ! py-s * dxy-s / ds_x
            H_cf(2,3,5) = H_cf(2,3,5) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,2) * Mjs(5) + Mjs(3) * dMjs(2,4) ) ! py-s * dxy-s / ds_y
            H_cf(3,3,5) = H_cf(3,3,5) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,2) * Mjs(5) + Mjs(3) * dMjs(3,4) ) ! py-s * dxy-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,5,3) = H_cf(:,3,5)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(6)
            H_cf(1,3,6) = H_cf(1,3,6) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,2) * Mjs(6) + Mjs(3) * dMjs(1,5) ) ! py-s * dxz-s / ds_x
            H_cf(2,3,6) = H_cf(2,3,6) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,2) * Mjs(6) + Mjs(3) * dMjs(2,5) ) ! py-s * dxz-s / ds_y
            H_cf(3,3,6) = H_cf(3,3,6) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,2) * Mjs(6) + Mjs(3) * dMjs(3,5) ) ! py-s * dxz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,6,3) = H_cf(:,3,6)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(7)
            H_cf(1,3,7) = H_cf(1,3,7) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,2) * Mjs(7) + Mjs(3) * dMjs(1,6) ) ! py-s * dyz-s / ds_x
            H_cf(2,3,7) = H_cf(2,3,7) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,2) * Mjs(7) + Mjs(3) * dMjs(2,6) ) ! py-s * dyz-s / ds_y
            H_cf(3,3,7) = H_cf(3,3,7) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,2) * Mjs(7) + Mjs(3) * dMjs(3,6) ) ! py-s * dyz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,7,3) = H_cf(:,3,7)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(8)
            H_cf(1,3,8) = H_cf(1,3,8) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,2) * Mjs(8) + Mjs(3) * dMjs(1,7) ) ! py-s * (dx2-y2) / ds_x
            H_cf(2,3,8) = H_cf(2,3,8) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,2) * Mjs(8) + Mjs(3) * dMjs(2,7) ) ! py-s * (dx2-y2) / ds_y
            H_cf(3,3,8) = H_cf(3,3,8) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,2) * Mjs(8) + Mjs(3) * dMjs(3,7) ) ! py-s * (dx2-y2) / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,8,3) = H_cf(:,3,8)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(9)
            H_cf(1,3,9) = H_cf(1,3,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,2) * Mjs(9) + Mjs(3) * dMjs(1,8) ) ! py-s * (d3z2-r2) / ds_x
            H_cf(2,3,9) = H_cf(2,3,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,2) * Mjs(9) + Mjs(3) * dMjs(2,8) ) ! py-s * (d3z2-r2) / ds_y
            H_cf(3,3,9) = H_cf(3,3,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,2) * Mjs(9) + Mjs(3) * dMjs(3,8) ) ! py-s * (d3z2-r2) / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,3) = H_cf(:,3,9)


            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(5)
            H_cf(1,4,5) = H_cf(1,4,5) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,3) * Mjs(5) + Mjs(4) * dMjs(1,4) ) ! pz-s * dxy / ds_x
            H_cf(2,4,5) = H_cf(2,4,5) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,3) * Mjs(5) + Mjs(4) * dMjs(2,4) ) ! pz-s * dxy / ds_y
            H_cf(3,4,5) = H_cf(3,4,5) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,3) * Mjs(5) + Mjs(4) * dMjs(3,4) ) ! pz-s * dxy / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,5,4) = H_cf(:,4,5)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(6)
            H_cf(1,4,6) = H_cf(1,4,6) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,3) * Mjs(6) + Mjs(4) * dMjs(1,5) ) ! pz-s * dxz / ds_x
            H_cf(2,4,6) = H_cf(2,4,6) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,3) * Mjs(6) + Mjs(4) * dMjs(2,5) ) ! pz-s * dxz / ds_y
            H_cf(3,4,6) = H_cf(3,4,6) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,3) * Mjs(6) + Mjs(4) * dMjs(3,5) ) ! pz-s * dxz / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,6,4) = H_cf(:,4,6)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(7)
            H_cf(1,4,7) = H_cf(1,4,7) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,3) * Mjs(7) + Mjs(4) * dMjs(1,6) ) ! pz-s * dyz / ds_x
            H_cf(2,4,7) = H_cf(2,4,7) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,3) * Mjs(7) + Mjs(4) * dMjs(2,6) ) ! pz-s * dyz / ds_y
            H_cf(3,4,7) = H_cf(3,4,7) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,3) * Mjs(7) + Mjs(4) * dMjs(3,6) ) ! pz-s * dyz / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,7,4) = H_cf(:,4,7)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(8)
            H_cf(1,4,8) = H_cf(1,4,8) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,3) * Mjs(8) + Mjs(4) * dMjs(1,7) ) ! pz-s * (dx2-y2) / ds_x
            H_cf(2,4,8) = H_cf(2,4,8) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,3) * Mjs(8) + Mjs(4) * dMjs(2,7) ) ! pz-s * (dx2-y2) / ds_y
            H_cf(3,4,8) = H_cf(3,4,8) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,3) * Mjs(8) + Mjs(4) * dMjs(3,7) ) ! pz-s * (dx2-y2) / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,8,4) = H_cf(:,4,8)

            H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(9)
            H_cf(1,4,9) = H_cf(1,4,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(2,3) * ( dMjs(1,3) * Mjs(9) + Mjs(4) * dMjs(1,8) ) ! pz-s * (d3z2-r2) / ds_x
            H_cf(2,4,9) = H_cf(2,4,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(2,3) * ( dMjs(2,3) * Mjs(9) + Mjs(4) * dMjs(2,8) ) ! pz-s * (d3z2-r2) / ds_y
            H_cf(3,4,9) = H_cf(3,4,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(2,3) * ( dMjs(3,3) * Mjs(9) + Mjs(4) * dMjs(3,8) ) ! pz-s * (d3z2-r2) / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,4) = H_cf(:,4,9)


            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(6)
            H_cf(1,5,6) = H_cf(1,5,6) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,4) * Mjs(6) + Mjs(5) * dMjs(1,5) ) ! dxy-s * dxz-s / ds_x
            H_cf(2,5,6) = H_cf(2,5,6) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,4) * Mjs(6) + Mjs(5) * dMjs(2,5) ) ! dxy-s * dxz-s / ds_y
            H_cf(3,5,6) = H_cf(3,5,6) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,4) * Mjs(6) + Mjs(5) * dMjs(3,5) ) ! dxy-s * dxz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,6,5) = H_cf(:,5,6)

            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(7)
            H_cf(1,5,7) = H_cf(1,5,7) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,4) * Mjs(7) + Mjs(5) * dMjs(1,6) ) ! dxy-s * dyz-s / ds_x
            H_cf(2,5,7) = H_cf(2,5,7) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,4) * Mjs(7) + Mjs(5) * dMjs(2,6) ) ! dxy-s * dyz-s / ds_y
            H_cf(3,5,7) = H_cf(3,5,7) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,4) * Mjs(7) + Mjs(5) * dMjs(3,6) ) ! dxy-s * dyz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,7,5) = H_cf(:,5,7)

            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(8)
            H_cf(1,5,8) = H_cf(1,5,8) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,4) * Mjs(8) + Mjs(5) * dMjs(1,7) ) ! dxy-s * (dx2-y2) / ds_x
            H_cf(2,5,8) = H_cf(2,5,8) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,4) * Mjs(8) + Mjs(5) * dMjs(2,7) ) ! dxy-s * (dx2-y2) / ds_y
            H_cf(3,5,8) = H_cf(3,5,8) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,4) * Mjs(8) + Mjs(5) * dMjs(3,7) ) ! dxy-s * (dx2-y2) / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,8,5) = H_cf(:,5,8)

            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(9)
            H_cf(1,5,9) = H_cf(1,5,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,4) * Mjs(9) + Mjs(5) * dMjs(1,8) ) ! dxy-s * (d3z2-r2) / ds_x
            H_cf(2,5,9) = H_cf(2,5,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,4) * Mjs(9) + Mjs(5) * dMjs(2,8) ) ! dxy-s * (d3z2-r2) / ds_y
            H_cf(3,5,9) = H_cf(3,5,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,4) * Mjs(9) + Mjs(5) * dMjs(3,8) ) ! dxy-s * (d3z2-r2) / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,5) = H_cf(:,5,9)


            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6) * Mjs(7)
            H_cf(1,6,7) = H_cf(1,6,7) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,5) * Mjs(7) + Mjs(6) * dMjs(1,6) ) ! dxz-s * dyz-s / ds_x
            H_cf(2,6,7) = H_cf(2,6,7) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,5) * Mjs(7) + Mjs(6) * dMjs(2,6) ) ! dxz-s * dyz-s / ds_y
            H_cf(3,6,7) = H_cf(3,6,7) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,5) * Mjs(7) + Mjs(6) * dMjs(3,6) ) ! dxz-s * dyz-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,7,6) = H_cf(:,6,7)

            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6) * Mjs(8)
            H_cf(1,6,8) = H_cf(1,6,8) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,5) * Mjs(8) + Mjs(6) * dMjs(1,7) ) ! dxz-s * (dx2-y2)-s / ds_x
            H_cf(2,6,8) = H_cf(2,6,8) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,5) * Mjs(8) + Mjs(6) * dMjs(2,7) ) ! dxz-s * (dx2-y2)-s / ds_y
            H_cf(3,6,8) = H_cf(3,6,8) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,5) * Mjs(8) + Mjs(6) * dMjs(3,7) ) ! dxz-s * (dx2-y2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,8,6) = H_cf(:,6,8)

            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6) * Mjs(9)
            H_cf(1,6,9) = H_cf(1,6,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,5) * Mjs(9) + Mjs(6) * dMjs(1,8) ) ! dxz-s * (d3z2-r2)-s / ds_x
            H_cf(2,6,9) = H_cf(2,6,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,5) * Mjs(9) + Mjs(6) * dMjs(2,8) ) ! dxz-s * (d3z2-r2)-s / ds_y
            H_cf(3,6,9) = H_cf(3,6,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,5) * Mjs(9) + Mjs(6) * dMjs(3,8) ) ! dxz-s * (d3z2-r2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,6) = H_cf(:,6,9)


            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(7) * Mjs(8)
            H_cf(1,7,8) = H_cf(1,7,8) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,6) * Mjs(8) + Mjs(7) * dMjs(1,7) ) ! dyz-s * (dx2-y2)-s / ds_x
            H_cf(2,7,8) = H_cf(2,7,8) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,6) * Mjs(8) + Mjs(7) * dMjs(2,7) ) ! dyz-s * (dx2-y2)-s / ds_y
            H_cf(3,7,8) = H_cf(3,7,8) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,6) * Mjs(8) + Mjs(7) * dMjs(3,7) ) ! dyz-s * (dx2-y2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,8,7) = H_cf(:,7,8)

            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(7) * Mjs(9)
            H_cf(1,7,9) = H_cf(1,7,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,7) * Mjs(9) + Mjs(7) * dMjs(1,8) ) ! dyz-s * (d3z2-r2)-s / ds_x
            H_cf(2,7,9) = H_cf(2,7,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,7) * Mjs(9) + Mjs(7) * dMjs(2,8) ) ! dyz-s * (d3z2-r2)-s / ds_y
            H_cf(3,7,9) = H_cf(3,7,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,7) * Mjs(9) + Mjs(7) * dMjs(3,8) ) ! dyz-s * (d3z2-r2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,7) = H_cf(:,7,9)


            H_cf_temp(1) = d_matr_spd(3,3) * Mjs(8) * Mjs(9)
            H_cf(1,8,9) = H_cf(1,8,9) + H_cf_temp(1) * drij_dsk(1) + &
               matr_spd(3,3) * ( dMjs(1,7) * Mjs(9) + Mjs(8) * dMjs(1,8) ) ! (dx2-y2)-s * (d3z2-r2)-s / ds_x
            H_cf(2,8,9) = H_cf(2,8,9) + H_cf_temp(1) * drij_dsk(2) + &
               matr_spd(3,3) * ( dMjs(2,7) * Mjs(9) + Mjs(8) * dMjs(2,8) ) ! (dx2-y2)-s * (d3z2-r2)-s / ds_y
            H_cf(3,8,9) = H_cf(3,8,9) + H_cf_temp(1) * drij_dsk(3) + &
               matr_spd(3,3) * ( dMjs(3,7) * Mjs(9) + Mjs(8) * dMjs(3,8) ) ! (dx2-y2)-s * (d3z2-r2)-s / ds_z
            ! Lower triangle use symmetry:
            H_cf(:,9,8) = H_cf(:,8,9)

         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)

5003 continue

   enddo ! atom_2 = 0, m

   !-----------------
   ! 4) 3-body contributions:
   if (TB(KOA1,KOA1)%include_3body) then  ! only if user defined it to include
      atom_2 = 0  ! to restart
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         KOA2 => Scell%MDatoms(j)%KOA   ! index of the second atom
         ! Distances between atom #1 and #2:
         xij => Scell%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
         yij => Scell%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
         zij => Scell%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
         rij => Scell%Near_neighbor_dist(i,atom_2,4)	! at this distance, R
         ! d r_ij / d s_k:
         drij_dsk(1) = drij_dska(i, j, k, xij, yij, zij, rij, Scell%supce, 1, .true.)	! dr_{ij}/ds_{k,x}
         drij_dsk(2) = drij_dska(i, j, k, xij, yij, zij, rij, Scell%supce, 2, .true.)	! dr_{ij}/ds_{k,y}
         drij_dsk(3) = drij_dska(i, j, k, xij, yij, zij, rij, Scell%supce, 3, .true.)	! dr_{ij}/ds_{k,z}

         ! To start summing up
         AT3:do atom_3 = 1,m ! do only for atoms close to that one
            kk = Scell%Near_neighbor_list(i,atom_3) ! this is the list of such close atoms
            ! Make sure the third atom is not the second atom:
            if (kk /= j) then
               KOA3 => Scell%MDatoms(kk)%KOA   ! kind of atom #3

               ! Distances between atom #1 and #3:
               xik => Scell%Near_neighbor_dist(i,atom_3,1)	! at this distance, X
               yik => Scell%Near_neighbor_dist(i,atom_3,2)	! at this distance, Y
               zik => Scell%Near_neighbor_dist(i,atom_3,3)	! at this distance, Z
               rik => Scell%Near_neighbor_dist(i,atom_3,4)	! at this distance, R
               ! Distances between atom #2 and #3:
               xjk => Scell%Near_neighbor_dist(j,atom_3,1)	! at this distance, X
               yjk => Scell%Near_neighbor_dist(j,atom_3,2)	! at this distance, Y
               zjk => Scell%Near_neighbor_dist(j,atom_3,3)	! at this distance, Z
               rjk => Scell%Near_neighbor_dist(j,atom_3,4)	! at this distance, R

               ! Derivatives of rij by sk:
               ! All functions are from the module "TB_Koster_Slater":
               drikk_dsk(1) = drij_dska(i, kk, k, xik, yik, zik, rik, Scell%supce, 1, .true.)	! dr_{ikk}/ds_{k,x}
               drikk_dsk(2) = drij_dska(i, kk, k, xik, yik, zik, rik, Scell%supce, 2, .true.)	! dr_{ikk}/ds_{k,y}
               drikk_dsk(3) = drij_dska(i, kk, k, xik, yik, zik, rik, Scell%supce, 3, .true.)	! dr_{ikk}/ds_{k,z}
               drjkk_dsk(1) = drij_dska(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 1, .true.)	! dr_{jkk}/ds_{k,x}
               drjkk_dsk(2) = drij_dska(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 2, .true.)	! dr_{jkk}/ds_{k,y}
               drjkk_dsk(3) = drij_dska(j, kk, k, xjk, yjk, zjk, rjk, Scell%supce, 3, .true.)	! dr_{jkk}/ds_{k,z}

               ! Find the combination-of-atoms index:
               at_ind = find_3bdy_ind(KOA1, KOA2, KOA3)  ! module "Dealing_with_3TB"

               ! Get the 3-body elements for each shells combination:
               term_3bdy(:) = TB(KOA1,KOA2)%Hh3bdy(at_ind, 1) * ( &
                  M_d_Lag_exp(i,j,1) * M_Lag_exp(j,kk,1) * M_Lag_exp(i,kk,1) * drij_dsk(:)  + &
                  M_Lag_exp(i,j,1) * M_d_Lag_exp(j,kk,1) * M_Lag_exp(i,kk,1) * drjkk_dsk(:) + &
                  M_Lag_exp(i,j,1) * M_Lag_exp(j,kk,1) * M_d_Lag_exp(i,kk,1) * drikk_dsk(:) ) + &
                              TB(KOA1,KOA2)%Hh3bdy(at_ind, 2) * ( &
                  M_d_Lag_exp(i,j,2) * M_Lag_exp(j,kk,1) * M_Lag_exp(i,kk,1) * drij_dsk(:)  + &
                  M_Lag_exp(i,j,2) * M_d_Lag_exp(j,kk,1) * M_Lag_exp(i,kk,1) * drjkk_dsk(:) + &
                  M_Lag_exp(i,j,2) * M_Lag_exp(j,kk,1) * M_d_Lag_exp(i,kk,1) * drikk_dsk(:) ) + &
                              TB(KOA1,KOA2)%Hh3bdy(at_ind, 3) * ( &
                  M_d_Lag_exp(i,j,1) * M_Lag_exp(j,kk,2) * M_Lag_exp(i,kk,1) * drij_dsk(:)  + &
                  M_Lag_exp(i,j,1) * M_d_Lag_exp(j,kk,2) * M_Lag_exp(i,kk,1) * drjkk_dsk(:) + &
                  M_Lag_exp(i,j,1) * M_Lag_exp(j,kk,2) * M_d_Lag_exp(i,kk,1) * drikk_dsk(:) ) + &
                              TB(KOA1,KOA2)%Hh3bdy(at_ind, 4) * ( &
                  M_d_Lag_exp(i,j,1) * M_Lag_exp(j,kk,1) * M_Lag_exp(i,kk,2) * drij_dsk(:)  + &
                  M_Lag_exp(i,j,1) * M_d_Lag_exp(j,kk,1) * M_Lag_exp(i,kk,2) * drjkk_dsk(:) + &
                  M_Lag_exp(i,j,1) * M_Lag_exp(j,kk,1) * M_d_Lag_exp(i,kk,2) * drikk_dsk(:) )

               do sh1 = 1, Bsiz  ! for all shells
                  H_3bdy(:,sh1,sh1) = H_3bdy(:,sh1,sh1) + term_3bdy(:)   ! no orbital resolution
               enddo ! sh1
            endif ! (k /= i)
         enddo AT3
      enddo AT2
   endif ! (TB(KOA1,KOA1)%include_3body)

5002 continue
   ! Collect all the terms into derivative of the Hamiltonian:
   dHij = E_onsite + H_avg + H_cf + H_3bdy

   deallocate(E_onsite, H_avg, H_cf, H_3bdy)
   nullify(m, KOA1, KOA2, KOA3, x1, y1, z1, r1)
   nullify(xij, yij, zij, rij)
   nullify(xik, yik, zik, rik)
   nullify(xjk, yjk, zjk, rjk)
end subroutine d_Onsite_3TB




!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! Attractive forces for supercell from the derivatives of the Hamiltonian:
subroutine Attract_TB_Forces_Press_3TB(Scell, NSC, TB, numpar, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, &
                                       Aij_x_Ei, Mjs_in, M_Lag_exp, M_d_Lag_exp)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: M_d_Lag_exp   ! matrix of derivatives of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: Mjs_in ! matrix of overlaps with s-orbital
   !------------------------------------------------------------------
   real(8), allocatable, dimension(:,:) :: dwr_press, dS_press
   real(8), allocatable, dimension(:,:,:) :: dHij, dSij
   integer i, j, k, n
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   if (numpar%p_const) then	! calculate this for P=const Parrinello-Rahman MD
      n = size(Aij,1)
      allocate(dwr_press(9,n))
      allocate(dS_press(9,n))
      allocate(dHij(9,n,size(Aij,2)))
      allocate(dSij(9,n,size(Aij,2)))
      dHij = 0.0d0
      dSij = 0.0d0
      dwr_press = 0.0d0
      dS_press = 0.0d0

      call dHamil_tot_Press_3TB(Scell, NSC, TB, numpar, M_Vij, M_dVij, M_SVij, M_dSVij, &
                                 M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dHij, dSij)   ! below

#ifdef MPI_USED   ! use the MPI version
      N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
      Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
      Nend = n
      ! Do the cycle (parallel) calculations:
      do i = Nstart, Nend, N_incr  ! each process does its own part
      !do i = 1, n
         do j = 1, 9
            dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(i,:)) ! old, tested, good
            dS_press(j,i) = dS_press(j,i) + SUM(dSij(j,i,:)*Aij_x_Ei(i,:))
         enddo ! i
      enddo ! j
      ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
      error_part = 'Error in Attract_TB_Forces_Press_3TB:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dwr_press', dwr_press) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dS_press', dS_press) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
      !$omp PARALLEL DO private(i,j)
      do j = 1, 9
         do i = 1, n
            dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(i,:)) ! old, tested, good
            dS_press(j,i) = dS_press(j,i) + SUM(dSij(j,i,:)*Aij_x_Ei(i,:))
         enddo ! i
      enddo ! j
      !$OMP END PARALLEL DO
#endif

      Scell(NSC)%SCforce%att = 0.0d0
      do i = 1,3
         do k = 1,3
!             Scell(NSC)%SCforce%att(k,i) = SUM(dwr_press((i-1)*3+k,:)) - SUM(dS_press((i-1)*3+k,:))
            Scell(NSC)%SCforce%att(i,k) = SUM(dwr_press((i-1)*3+k,:)) - SUM(dS_press((i-1)*3+k,:)) ! test
         enddo ! k
      enddo ! i
      deallocate(dwr_press, dS_press, dHij, dSij)
   endif
end subroutine Attract_TB_Forces_Press_3TB



subroutine dHamil_tot_Press_3TB(Scell, NSC, TB, numpar, M_Vij, M_dVij, M_SVij, M_dSVij, &
                                 M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dHij, dSij)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(in) :: Mjs_in
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: M_d_Lag_exp   ! matrix of derivatives of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(inout) :: dHij, dSij
   !-----------------------------------------------------
   real(8), dimension(:,:,:), allocatable :: dHij1, dSij1
   integer :: i, j, j1, i1, atom_2, m, nat, i2, j2
   integer :: i4, j4, norb, n_overlap
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   ! Depending on the basis set:
   n_overlap = identify_DFTB_basis_size(numpar%basis_size_ind)   ! module "TB_DFTB"
   norb = identify_DFTB_orbitals_per_atom(numpar%basis_size_ind)    ! module "TB_DFTB"
   nat = size(Scell(NSC)%MDatoms)	! number of atoms
   dHij = 0.0d0 ! to start with
   dSij = 0.0d0 ! to start with

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = nat
   if (.not.allocated(dHij1)) allocate(dHij1(9,norb,norb))
   if (.not.allocated(dSij1)) allocate(dSij1(9,norb,norb))
   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(i)
      i4 = (i-1)*norb
      do atom_2 = 0,m ! do only for atoms close to that one
         if (atom_2 == 0) then
            j = i
         else
            j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         endif
         if (j .GT. 0) then
            j4 = (j-1)*norb
            call dHamilton_one_Press_3TB(i, atom_2, numpar%basis_size_ind, Scell, NSC, TB, norb, n_overlap, &
                                          M_Vij, M_dVij, M_SVij, M_dSVij, &
                                          M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dHij1, dSij1)  ! below
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis
            do j1 = 1,norb	! all orbitals
               j2 = j4+j1
               do i1 = 1,norb	! all orbitals
                  i2 = i4+i1
                  dHij(:,i2,j2) = dHij1(:,i1,j1)	! construct the total Hamiltonian from
                  dSij(:,i2,j2) = dSij1(:,i1,j1)	! construct the total Overlap Matrix from
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0) then
      enddo ! j
   enddo ! i
   if (allocated(dHij1)) deallocate(dHij1)
   if (allocated(dSij1)) deallocate(dSij1)
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in dHamil_tot_Press_3TB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dHij', dHij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dSij', dSij) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
   !$omp parallel private(i,m,i4,atom_2,j,j4,j1,i1,j2,i2,dHij1,dSij1)
   if (.not.allocated(dHij1)) allocate(dHij1(9,norb,norb))
   if (.not.allocated(dSij1)) allocate(dSij1(9,norb,norb))
   !$omp do
   do i = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(i)
      i4 = (i-1)*norb
      do atom_2 = 0,m ! do only for atoms close to that one
         if (atom_2 == 0) then
            j = i
         else
            j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         endif
         if (j .GT. 0) then
            j4 = (j-1)*norb
            call dHamilton_one_Press_3TB(i, atom_2, numpar%basis_size_ind, Scell, NSC, TB, norb, n_overlap, &
                                          M_Vij, M_dVij, M_SVij, M_dSVij, &
                                          M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dHij1, dSij1)  ! below
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis
            do j1 = 1,norb	! all orbitals
               j2 = j4+j1
               do i1 = 1,norb	! all orbitals
                  i2 = i4+i1
                  dHij(:,i2,j2) = dHij1(:,i1,j1)	! construct the total Hamiltonian from
                  dSij(:,i2,j2) = dSij1(:,i1,j1)	! construct the total Overlap Matrix from
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0) then
      enddo ! j
   enddo ! i
   !$omp end do
   if (allocated(dHij1)) deallocate(dHij1)
   if (allocated(dSij1)) deallocate(dSij1)
   !$omp end parallel
#endif
end subroutine dHamil_tot_Press_3TB


subroutine dHamilton_one_Press_3TB(i, atom_2, basis_size, Scell, NSC, TB, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, &
                                    M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dHij_press, dSij_press)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
   integer, intent(in) :: i, atom_2
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: basis_size, NSC ! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   integer, intent(in) :: norb, n_overlap  ! number of orbitals per atom (depends on the basis set), and number of overlap functions
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! Hopping for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! Overlap for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(in) :: Mjs_in
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: M_d_Lag_exp   ! matrix of derivatives of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(out) :: dHij_press, dSij_press
   !---------------------------------------------------
   integer :: ki, kj
   real(8), dimension(9,norb,norb) ::  dH, dS ! hopping integrals
   dHij_press = 0.0d0
   if (atom_2 == 0) then ! Onsite contributions are constant, so derivatives are always zero:
      dH = 0.0d0  ! to start with
      call d_Onsite_Press_3TB(i, Scell(NSC), TB, basis_size, norb, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dH) ! below
      dS = 0.0d0
      do ki = 1, norb
         do kj = 1, norb
            dHij_press(:,ki,kj) = dH(:,ki,kj)   ! Hopping Integrals
         enddo ! kj
      enddo  ! ki
      dSij_press = 0.0d0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings:
      call dHopping_Press_3TB(i, atom_2, Scell, NSC, TB, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn,  dH, dS)    ! below
      do ki = 1, norb
         do kj = 1, norb
            dHij_press(:,ki,kj) = dH(:,ki,kj)   ! Hopping Integrals ! CORRECT, tested on sp3 basis set
            dSij_press(:,ki,kj) = dS(:,ki,kj)   ! Hopping Integrals
         enddo ! kj
      enddo  ! ki
   endif
end subroutine dHamilton_one_Press_3TB


subroutine d_Onsite_Press_3TB(i, Scell, TB, basis_ind, norb, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Mjs_in, M_Lag_exp, M_d_Lag_exp, dH) ! below
   integer, intent(in) :: i
   type(Super_cell), intent(inout), target :: Scell	! supercell with all the atoms as one object
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   integer, intent(in) :: basis_ind, norb   ! number of orbitals per atom
   real(8), dimension(:,:,:), intent(in), target :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in), target :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(in) :: Mjs_in ! symmetry factors for s-s,p-s,d-s overlap
   real(8), dimension(:,:,:), intent(in) :: M_Lag_exp   ! matrix of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(in) :: M_d_Lag_exp   ! matrix of derivatives of laguerre * exp(-a*r_ij) * cutoff
   real(8), dimension(:,:,:), intent(out) :: dH
   !---------------------------------------------
   integer :: i1, j1, k1, Bsiz, atom_2, sh1, sh2, j, at_ind, atom_3, kk
   real(8) :: r, term_s, term_p, term_d, term_3bdy(3)
   real(8) :: matr_spd(3,3), d_matr_spd(3,3), H_cf_temp(9), drij_dsk(3), drikk_dsk(3), drjkk_dsk(3)
   real(8) :: rij(3), sij(3), drij_dh
   real(8), dimension(3) :: M_dlmn
   real(8), dimension(8) :: dMjs ! 3p, 5d
   real(8), dimension(9) :: Mjs ! 1s, 3p, 5d
   real(8), dimension(:,:,:), allocatable :: E_onsite, H_avg, H_cf, H_3bdy
   integer, pointer :: m, KOA1, KOA2, KOA3
   real(8), pointer :: r1, x1, y1, z1

   dH = 0.0d0   ! to start with

   ! Kind of atom:
   KOA1 => Scell%MDatoms(i)%KOA
   ! Number of the nearest neighbors
   m => Scell%Near_neighbor_size(i)

   ! Size of the basis:
   Bsiz = norb

   ! To start with:
   allocate(E_onsite(9,Bsiz,Bsiz), source=0.0d0)
   allocate(H_avg(9,Bsiz,Bsiz), source=0.0d0)
   allocate(H_cf(9,Bsiz,Bsiz), source=0.0d0)
   allocate(H_3bdy(9,Bsiz,Bsiz), source=0.0d0)

   ! The onsite energies are constructed out of 4 terms [1]:
   ! 1) onsite eigenvalues
   ! 2) average contribution of atoms around
   ! 3) crystall field contribution
   ! 4) 3-body contributions
   ! Let's calculate them all:


   !-----------------
   ! 1) Onsite energies of spin-unpolirized orbital values:
   E_onsite = 0.0d0     ! constants -> all derivatives = 0

   !-----------------
   ! 2,3) Average and crystal field terms (tested, correct):
   do atom_2 = 1, m ! do only for atoms close to that one
      j = Scell%Near_neighbor_list(i, atom_2) ! this is the list of such close atoms
      ! [OS 0]
!       KOA1 => Scell%MDatoms(i)%KOA
!       KOA2 => Scell%MDatoms(j)%KOA
      ! [OS 1] tested, correct
      KOA1 => Scell%MDatoms(j)%KOA
      KOA2 => Scell%MDatoms(i)%KOA

      r1 => Scell%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
      rij(:) = Scell%Near_neighbor_dist(i,atom_2,1:3)  ! at this distance, X
      sij(:) = Scell%Near_neighbor_dist_s(i,atom_2,:) ! at this distance, SX

      ! Radial parts for aferage atom contribution:
      term_s = SUM( TB(KOA1,KOA2)%Hhavg(1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s orbitals
      if (basis_ind > 0) then ! p3 orbitals:
         term_p = SUM( TB(KOA1,KOA2)%Hhavg(2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p3 orbitals
         if (basis_ind > 1) then ! d5 orbitals:
            term_d = SUM( TB(KOA1,KOA2)%Hhavg(3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d5 orbitals
         endif
      endif

      ! Radial parts for crystal field:
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
      ! and its derivatives:
      d_matr_spd(1,1) = SUM( TB(KOA1,KOA2)%Hhcf(1,1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s-s orbitals
      if (basis_ind > 0) then ! p3 orbitals:
         d_matr_spd(1,2) = SUM( TB(KOA1,KOA2)%Hhcf(1,2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s-p orbitals
         d_matr_spd(2,1) = SUM( TB(KOA1,KOA2)%Hhcf(2,1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p-s orbitals
         d_matr_spd(2,2) = SUM( TB(KOA1,KOA2)%Hhcf(2,2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p-p orbitals
         if (basis_ind > 1) then ! d5 orbitals:
            d_matr_spd(1,3) = SUM( TB(KOA1,KOA2)%Hhcf(1,3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! s-d orbitals
            d_matr_spd(3,1) = SUM( TB(KOA1,KOA2)%Hhcf(3,1,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d-s orbitals
            d_matr_spd(2,3) = SUM( TB(KOA1,KOA2)%Hhcf(2,3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! p-d orbitals
            d_matr_spd(3,2) = SUM( TB(KOA1,KOA2)%Hhcf(3,2,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d-p orbitals
            d_matr_spd(3,3) = SUM( TB(KOA1,KOA2)%Hhcf(3,3,1:4) * M_d_Lag_exp(i,j,1:4) )  ! d-d orbitals
         endif ! (basis_ind > 0)
      endif ! (basis_ind > 1)


      do i1 = 1, 3	! gamma
         do j1 = 1, 3	! delta
            k1 = (i1-1)*3 + j1
            ! all the components of the h_alpha_beta(3,3):

            drij_dh = drij_dhab(rij(j1), sij(i1), r1)	! dr_{ij}/dh_{gamma,delta}, module "TB_Koster_Slater"
            !drij_dh = drij_dhab(rij(i1), sij(j1), r1)	! dr_{ij}/dh_{gamma,delta}, module "TB_Koster_Slater"

            ! 2) Average terms:
            H_avg(k1,1,1) = H_avg(k1,1,1) + term_s * drij_dh
            if (basis_ind > 0) then ! p3 orbitals:
               do sh1 = 2, 4
                  H_avg(k1,sh1,sh1) = H_avg(k1,sh1,sh1) + term_p * drij_dh
               enddo
               if (basis_ind > 1) then ! d5 orbitals:
                  do sh1 = 5, 9
                     H_avg(k1,sh1,sh1) = H_avg(k1,sh1,sh1) + term_d * drij_dh
                  enddo
               endif ! (basis_ind > 0)
            endif ! (basis_ind > 1)


            !-----------------
            ! 3) Crystal field:
            H_cf(k1,1,1) = H_cf(k1,1,1) + d_matr_spd(1,1) * drij_dh ! s-s * s-s / dh

            M_dlmn(1) = dda_dhgd(1, j1, rij(1), rij(j1), sij(i1), r1)	! dl/dh{gamma,delta}, module "TB_Koster_Slater"
            M_dlmn(2) = dda_dhgd(2, j1, rij(2), rij(j1), sij(i1), r1)	! dm/dh{gamma,delta}, module "TB_Koster_Slater"
            M_dlmn(3) = dda_dhgd(3, j1, rij(3), rij(j1), sij(i1), r1)	! dn/dh{gamma,delta}, module "TB_Koster_Slater"
!             M_dlmn(1) = dda_dhgd(1, j1, rij(1), rij(i1), sij(j1), r1)	! dl/dh{gamma,delta}, module "TB_Koster_Slater"
!             M_dlmn(2) = dda_dhgd(2, j1, rij(2), rij(i1), sij(j1), r1)	! dm/dh{gamma,delta}, module "TB_Koster_Slater"
!             M_dlmn(3) = dda_dhgd(3, j1, rij(3), rij(i1), sij(j1), r1)	! dn/dh{gamma,delta}, module "TB_Koster_Slater"
            ! Koster-Slater smmetry factors:
            if (basis_ind > 0) then ! p3 orbitals:
               Mjs(2) = Mjs_in(i,j,2)  ! px-s
               Mjs(3) = Mjs_in(i,j,3)  ! py-s
               Mjs(4) = Mjs_in(i,j,4)  ! pz-s

               dMjs(1) = -M_dlmn(1) ! dl/dh{gamma,delta}
               dMjs(2) = -M_dlmn(2) ! dm/dh{gamma,delta}
               dMjs(3) = -M_dlmn(3) ! dn/dh{gamma,delta}

               H_cf_temp(1) = d_matr_spd(2,1) * Mjs(2)
               H_cf(k1,1,2) = H_cf(k1,1,2) + H_cf_temp(1) * drij_dh + matr_spd(2,1) * dMjs(1)  ! s-s * px-s / dh
               ! Lower triangle use symmetry:
               H_cf(k1,2,1) = H_cf(k1,1,2)

               H_cf_temp(1) = d_matr_spd(2,1) * Mjs(3)
               H_cf(k1,1,3) = H_cf(k1,1,3) + H_cf_temp(1) * drij_dh + matr_spd(2,1) * dMjs(2)  ! s-s * py-s / dh
               ! Lower triangle use symmetry:
               H_cf(k1,3,1) = H_cf(k1,1,3)

               H_cf_temp(1) = d_matr_spd(2,1) * Mjs(4)
               H_cf(k1,1,4) = H_cf(k1,1,4) + H_cf_temp(1) * drij_dh + matr_spd(2,1) * dMjs(3)  ! s-s * pz-s / ds_x
               ! Lower triangle use symmetry:
               H_cf(k1,4,1) = H_cf(k1,1,4)

               ! Diagonal part is excluded (no self-interaction of orbitals):
               if (TB(KOA1,KOA2)%nullify_diag_cf) then
                  H_cf(k1,2,2) = 0.0d0
                  H_cf(k1,3,3) = 0.0d0
                  H_cf(k1,4,4) = 0.0d0
               else
                  H_cf_temp(1) = d_matr_spd(2,2) * Mjs(2) * Mjs(2)
                  H_cf(k1,2,2) = H_cf(k1,2,2) + H_cf_temp(1) * drij_dh + &
                        matr_spd(2,2) * 2.0d0* dMjs(1) * Mjs(2) ! px-s * px-s /ds_x
                  H_cf_temp(1) = d_matr_spd(2,2) * Mjs(3) * Mjs(3)
                  H_cf(k1,3,3) = H_cf(k1,3,3) + H_cf_temp(1) * drij_dh + &
                        matr_spd(2,2) * 2.0d0* dMjs(2) * Mjs(3) ! py-s * py-s /ds_x
                  H_cf_temp(1) = d_matr_spd(2,2) * Mjs(4) * Mjs(4)
                  H_cf(k1,4,4) = H_cf(k1,4,4) + H_cf_temp(1) * drij_dh + &
                        matr_spd(2,2) * 2.0d0* dMjs(3) * Mjs(4) ! pz-s * pz-s /ds_x
               endif

               ! Off-diagonal part:
               H_cf_temp(1) = d_matr_spd(2,2) * Mjs(2) * Mjs(3)
               H_cf(k1,2,3) = H_cf(k1,2,3) + H_cf_temp(1) * drij_dh + &
                        matr_spd(2,2) * (dMjs(1) * Mjs(3) + Mjs(2) * dMjs(2) ) ! px-s * py-s /ds_x
               ! For lower triangle use symmetry:
               H_cf(k1,3,2) = H_cf(k1,2,3) ! py-s * px-s /ds_{x,y,z}

               H_cf_temp(1) = d_matr_spd(2,2) * Mjs(2) * Mjs(4)
               H_cf(k1,2,4) = H_cf(k1,2,4) + H_cf_temp(1) * drij_dh + &
                        matr_spd(2,2) * (dMjs(1) * Mjs(4) + Mjs(2) * dMjs(3) ) ! px-s * pz-s /ds_x
               ! For lower triangle use symmetry:
               H_cf(k1,4,2) = H_cf(k1,2,4) ! pz-s * px-s /ds_{x,y,z}


               H_cf_temp(1) = d_matr_spd(2,2) * Mjs(3) * Mjs(4)
               H_cf(k1,3,4) = H_cf(k1,3,4) + H_cf_temp(1) * drij_dh + &
                        matr_spd(2,2) * (dMjs(2) * Mjs(4) + Mjs(3) * dMjs(3) ) ! py-s * pz-s /ds_x
               ! For lower triangle use symmetry:
               H_cf(k1,4,3) = H_cf(k1,3,4)  ! pz-s * py-s /ds_{x,y,z}

               if (basis_ind > 1) then ! d5 orbitals:
                  Mjs(5) = Mjs_in(i,j,5)  ! dxy-s
                  Mjs(6) = Mjs_in(i,j,6)  ! dxz-s
                  Mjs(7) = Mjs_in(i,j,7)  ! dyz-s
                  Mjs(8) = Mjs_in(i,j,8)  ! (dx2-y2)-s
                  Mjs(9) = Mjs_in(i,j,9)  ! (d3z2-r2)-s
                  ! Define d-symmetry factors to be reused below:
                  dMjs(4) = m_sqrt3 * (M_dlmn(1)*M_lmn(2,i,j) + M_lmn(1,i,j)*M_dlmn(2)) ! dxy-s / dh
                  dMjs(5) = m_sqrt3 * (M_dlmn(1)*M_lmn(3,i,j) + M_lmn(1,i,j)*M_dlmn(3)) ! dxz-s / dh
                  dMjs(6) = m_sqrt3 * (M_dlmn(2)*M_lmn(3,i,j) + M_lmn(2,i,j)*M_dlmn(3)) ! dyz-s / dh
                  dMjs(7) = m_sqrt3 * (M_dlmn(1)*M_lmn(1,i,j) - M_lmn(2,i,j)*M_dlmn(2)) ! (dx2-y2)-s / dh
                  dMjs(8) = 2.0d0*M_dlmn(3)*M_lmn(3,i,j) - (M_dlmn(1)*M_lmn(1,i,j) + M_lmn(2,i,j)*M_dlmn(2)) ! (d3z2-r2)-s / dh


                  ! Diagonal terms:
                  if (TB(KOA1,KOA2)%nullify_diag_cf) then
                     H_cf(k1,5,5) = 0.0d0  ! dxy-s * dxy-s
                     H_cf(k1,6,6) = 0.0d0  ! dxz-s * dxz-s
                     H_cf(k1,7,7) = 0.0d0  ! dyz-s * dyz-s
                     H_cf(k1,8,8) = 0.0d0  ! (dx2-y2) * (dx2-y2)
                     H_cf(k1,9,9) = 0.0d0  ! (d3z2-r2) * (d3z2-r2)
                  else
                     H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5)
                     H_cf(k1,5,5) = H_cf(k1,5,5) + H_cf_temp(1) * drij_dh + &
                              matr_spd(3,3) * 2.0d0 * Mjs(5) * dMjs(4)  ! dxy-s * dxy-s / dh
                     H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6)
                     H_cf(k1,6,6) = H_cf(k1,6,6) + H_cf_temp(1) * drij_dh + &
                              matr_spd(3,3) * 2.0d0 * Mjs(6) * dMjs(5)  ! dxz-s * dxz-s / dh
                     H_cf_temp(1) = d_matr_spd(3,3) * Mjs(7)
                     H_cf(k1,7,7) = H_cf(k1,7,7) + H_cf_temp(1) * drij_dh + &
                              matr_spd(3,3) * 2.0d0 * Mjs(7) * dMjs(6)  ! dyz-s * dyz-s / dh
                     H_cf_temp(1) = d_matr_spd(3,3) * Mjs(8)
                     H_cf(k1,8,8) = H_cf(k1,8,8) + H_cf_temp(1) * drij_dh + &
                              matr_spd(3,3) * 2.0d0 * Mjs(8) * dMjs(7)  ! (dx2-y2) * (dx2-y2) / dh
                     H_cf_temp(1) = d_matr_spd(3,3) * Mjs(9)
                     H_cf(k1,9,9) = H_cf(k1,9,9) + H_cf_temp(1) * drij_dh + &
                              matr_spd(3,3) * 2.0d0 * Mjs(9) * dMjs(8)  ! (d3z2-r2) * (d3z2-r2) / dh
                  endif

                  ! Off-diagonal terms:
                  H_cf_temp(1) = d_matr_spd(1,3) * Mjs(5)
                  H_cf(k1,1,5) = H_cf(k1,1,5) + H_cf_temp(1) * drij_dh + matr_spd(1,3) * dMjs(4)  ! s-s * dxy-s / dh
                  ! For lower triangle use symmetry:
                  H_cf(k1,5,1) = H_cf(k1,1,5)  ! dxy-s * s-s/ dh

                  H_cf_temp(1) = d_matr_spd(1,3) * Mjs(6)
                  H_cf(k1,1,6) = H_cf(k1,1,6) + H_cf_temp(1) * drij_dh + matr_spd(1,3) * dMjs(5)  ! s-s * dxz-s / dh
                  ! For lower triangle use symmetry:
                  H_cf(k1,6,1) = H_cf(k1,1,6)  ! dxz-s * s-s/ dh

                  H_cf_temp(1) = d_matr_spd(1,3) * Mjs(7)
                  H_cf(k1,1,7) = H_cf(k1,1,7) + H_cf_temp(1) * drij_dh + matr_spd(1,3) * dMjs(6)  ! s-s * dyz-s / dh
                  ! For lower triangle use symmetry:
                  H_cf(k1,7,1) = H_cf(k1,1,7)  ! dyz-s * s-s/ dh

                  H_cf_temp(1) = d_matr_spd(1,3) * Mjs(8)
                  H_cf(k1,1,8) = H_cf(k1,1,8) + H_cf_temp(1) * drij_dh + matr_spd(1,3) * dMjs(7)  ! s-s * (dx2-y2)-s / dh
                  ! For lower triangle use symmetry:
                  H_cf(k1,8,1) = H_cf(k1,1,8)  ! (dx2-y2)-s * s-s / dh

                  H_cf_temp(1) = d_matr_spd(1,3) * Mjs(9)
                  H_cf(k1,1,9) = H_cf(k1,1,9) + H_cf_temp(1) * drij_dh + matr_spd(1,3) * dMjs(8)  ! s-s * (d3z2-r2)-s / dh
                  ! For lower triangle use symmetry:
                  H_cf(k1,9,1) = H_cf(k1,1,9)  ! (d3z2-r2)-s * s-s / dh


                  ! Calculate repeating part the K-S matrix elements:
                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(5)
                  H_cf(k1,2,5) = H_cf(k1,2,5) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(1) * Mjs(5) + Mjs(2) * dMjs(4) ) ! px-s * dxy-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,5,2) = H_cf(k1,2,5)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(6)
                  H_cf(k1,2,6) = H_cf(k1,2,6) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(1) * Mjs(6) + Mjs(2) * dMjs(5) ) ! px-s * dxz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,6,2) = H_cf(k1,2,6)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(7)
                  H_cf(k1,2,7) = H_cf(k1,2,7) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(1) * Mjs(7) + Mjs(2) * dMjs(6) ) ! px-s * dyz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,7,2) = H_cf(k1,2,7)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(8)
                  H_cf(k1,2,8) = H_cf(k1,2,8) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(1) * Mjs(8) + Mjs(2) * dMjs(7) ) ! px-s * (dx2-y2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,8,2) = H_cf(k1,2,8)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(2) * Mjs(9)
                  H_cf(k1,2,9) = H_cf(k1,2,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(1) * Mjs(9) + Mjs(2) * dMjs(8) ) ! px-s * (d3z2-r2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,2) = H_cf(k1,2,9)


                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(5)
                  H_cf(k1,3,5) = H_cf(k1,3,5) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(2) * Mjs(5) + Mjs(3) * dMjs(4) ) ! py-s * dxy-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,5,3) = H_cf(k1,3,5)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(6)
                  H_cf(k1,3,6) = H_cf(k1,3,6) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(2) * Mjs(6) + Mjs(3) * dMjs(5) ) ! py-s * dxz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,6,3) = H_cf(k1,3,6)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(7)
                  H_cf(k1,3,7) = H_cf(k1,3,7) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(2) * Mjs(7) + Mjs(3) * dMjs(6) ) ! py-s * dyz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,7,3) = H_cf(k1,3,7)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(8)
                  H_cf(k1,3,8) = H_cf(k1,3,8) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(2) * Mjs(8) + Mjs(3) * dMjs(7) ) ! py-s * (dx2-y2) / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,8,3) = H_cf(k1,3,8)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(3) * Mjs(9)
                  H_cf(k1,3,9) = H_cf(k1,3,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(2) * Mjs(9) + Mjs(3) * dMjs(8) ) ! py-s * (d3z2-r2) / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,3) = H_cf(k1,3,9)


                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(5)
                  H_cf(k1,4,5) = H_cf(k1,4,5) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(3) * Mjs(5) + Mjs(4) * dMjs(4) ) ! pz-s * dxy / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,5,4) = H_cf(k1,4,5)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(6)
                  H_cf(k1,4,6) = H_cf(k1,4,6) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(3) * Mjs(6) + Mjs(4) * dMjs(5) ) ! pz-s * dxz / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,6,4) = H_cf(k1,4,6)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(7)
                  H_cf(k1,4,7) = H_cf(k1,4,7) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(3) * Mjs(7) + Mjs(4) * dMjs(6) ) ! pz-s * dyz / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,7,4) = H_cf(k1,4,7)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(8)
                  H_cf(k1,4,8) = H_cf(k1,4,8) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(3) * Mjs(8) + Mjs(4) * dMjs(7) ) ! pz-s * (dx2-y2) / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,8,4) = H_cf(k1,4,8)

                  H_cf_temp(1) = d_matr_spd(2,3) * Mjs(4) * Mjs(9)
                  H_cf(k1,4,9) = H_cf(k1,4,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(2,3) * ( dMjs(3) * Mjs(9) + Mjs(4) * dMjs(8) ) ! pz-s * (d3z2-r2) / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,4) = H_cf(k1,4,9)


                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(6)
                  H_cf(k1,5,6) = H_cf(k1,5,6) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(4) * Mjs(6) + Mjs(5) * dMjs(5) ) ! dxy-s * dxz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,6,5) = H_cf(k1,5,6)

                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(7)
                  H_cf(k1,5,7) = H_cf(k1,5,7) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(4) * Mjs(7) + Mjs(5) * dMjs(6) ) ! dxy-s * dyz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,7,5) = H_cf(k1,5,7)

                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(8)
                  H_cf(k1,5,8) = H_cf(k1,5,8) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(4) * Mjs(8) + Mjs(5) * dMjs(7) ) ! dxy-s * (dx2-y2) / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,8,5) = H_cf(k1,5,8)

                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(5) * Mjs(9)
                  H_cf(k1,5,9) = H_cf(k1,5,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(4) * Mjs(9) + Mjs(5) * dMjs(8) ) ! dxy-s * (d3z2-r2) / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,5) = H_cf(k1,5,9)


                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6) * Mjs(7)
                  H_cf(k1,6,7) = H_cf(k1,6,7) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(5) * Mjs(7) + Mjs(6) * dMjs(6) ) ! dxz-s * dyz-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,7,6) = H_cf(k1,6,7)

                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6) * Mjs(8)
                  H_cf(k1,6,8) = H_cf(k1,6,8) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(5) * Mjs(8) + Mjs(6) * dMjs(7) ) ! dxz-s * (dx2-y2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,8,6) = H_cf(k1,6,8)

                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(6) * Mjs(9)
                  H_cf(k1,6,9) = H_cf(k1,6,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(5) * Mjs(9) + Mjs(6) * dMjs(8) ) ! dxz-s * (d3z2-r2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,6) = H_cf(k1,6,9)


                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(7) * Mjs(8)
                  H_cf(k1,7,8) = H_cf(k1,7,8) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(6) * Mjs(8) + Mjs(7) * dMjs(7) ) ! dyz-s * (dx2-y2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,8,7) = H_cf(k1,7,8)

                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(7) * Mjs(9)
                  H_cf(k1,7,9) = H_cf(k1,7,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(7) * Mjs(9) + Mjs(7) * dMjs(8) ) ! dyz-s * (d3z2-r2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,7) = H_cf(k1,7,9)


                  H_cf_temp(1) = d_matr_spd(3,3) * Mjs(8) * Mjs(9)
                  H_cf(k1,8,9) = H_cf(k1,8,9) + H_cf_temp(1) * drij_dh + &
                     matr_spd(3,3) * ( dMjs(7) * Mjs(9) + Mjs(8) * dMjs(8) ) ! (dx2-y2)-s * (d3z2-r2)-s / dh
                  ! Lower triangle use symmetry:
                  H_cf(k1,9,8) = H_cf(k1,8,9)

               endif ! (basis_ind > 1)

            endif ! (basis_ind > 0)

5004 continue

         enddo ! j1
      enddo ! i1
   enddo ! atom_2 = 0, m

   !-----------------
   ! 4) 3-body contributions:
   if (TB(KOA1,KOA1)%include_3body) then  ! only if user defined it to include
      ! not ready yet
      H_3bdy = 0.0d0
   endif ! (TB(KOA1,KOA1)%include_3body)

5002 continue
   ! Collect all the terms into derivative of the Hamiltonian:
   dH = E_onsite + H_avg + H_cf + H_3bdy

   deallocate(E_onsite, H_avg, H_cf, H_3bdy)
   nullify(m, KOA1, KOA2, KOA3, x1, y1, z1, r1)
end subroutine d_Onsite_Press_3TB




subroutine dHopping_Press_3TB(i, atom_2, Scell, NSC, TB, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dH, dS)
! subroutine making the derivatives of the hopping integrals
   integer, intent(in) :: i, atom_2	! atoms
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB	  ! all tight binding parameters
   integer, intent(in) :: norb, n_overlap
   real(8), dimension(:,:,:), intent(in), target :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in), target :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(out) :: dH, dS
   !---------------------------------------------
   integer :: i1, j1, k1
   real(8) :: drij_dh, rij(3), sij(3)
   real(8), dimension(3) :: M_dlmn
   real(8), dimension(n_overlap) :: vec_M_dVij12, vec_M_dVij21, vec_M_dSVij12, vec_M_dSVij21
   real(8), dimension(n_overlap) ::  vec_M_Vij12, vec_M_Vij21, vec_M_SVij12, vec_M_SVij21
   real(8), dimension(norb,norb) :: dH1, dS1
   real(8), pointer :: r
   integer, pointer :: j, KOA1

   dH = 0.0d0
   dS = 0.0d0

   j => Scell(NSC)%Near_neighbor_list(i,atom_2)   ! it interacts with this atom
   r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
   rij(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
   rij(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
   rij(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
   sij(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
   sij(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
   sij(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

!    vec_M_Vij12 = M_Vij(i,j,:)
!    vec_M_Vij21 = M_Vij(j,i,:)
!    vec_M_SVij12 = M_SVij(i,j,:)
!    vec_M_SVij21 = M_SVij(j,i,:)
   vec_M_Vij21 = M_Vij(i,j,:)   ! correct
   vec_M_Vij12 = M_Vij(j,i,:)
   vec_M_SVij21 = M_SVij(i,j,:)
   vec_M_SVij12 = M_SVij(j,i,:)

   do i1 = 1, 3	! gamma
      do j1 = 1, 3	! delta
         k1 = (i1-1)*3 + j1
         ! all the components of the h_alpha_beta(3,3):
         ! TESTED, CORRECT:
         M_dlmn(1) = dda_dhgd(1, j1, rij(1), rij(j1), sij(i1), r)	! dl/dh{gamma,delta}, module "TB_Koster_Slater"
         M_dlmn(2) = dda_dhgd(2, j1, rij(2), rij(j1), sij(i1), r)	! dm/dh{gamma,delta}, module "TB_Koster_Slater"
         M_dlmn(3) = dda_dhgd(3, j1, rij(3), rij(j1), sij(i1), r)	! dn/dh{gamma,delta}, module "TB_Koster_Slater"

         ! TESTED, CORRECT:
         drij_dh = drij_dhab(rij(j1), sij(i1), r)	! dr_{ij}/dh_{gamma,delta}, module "TB_Koster_Slater"

         vec_M_dVij21 = M_dVij(i,j,:)*drij_dh   ! correct
         vec_M_dVij12 = M_dVij(j,i,:)*drij_dh
         vec_M_dSVij21 = M_dSVij(i,j,:)*drij_dh
         vec_M_dSVij12 = M_dSVij(j,i,:)*drij_dh

         select case (norb)
         case (1) ! s
            call d_KS_s(vec_M_Vij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dH1)	! module "TB_Koster_Slater"
            call d_KS_s(vec_M_SVij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dS1)	! module "TB_Koster_Slater"
         case (4) ! sp3
            call d_KS_sp3_hetero(vec_M_Vij12, vec_M_Vij21, vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3),  dH1)	! module "TB_Koster_Slater" ! CORRECT
            call d_KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3),  dS1)	! module "TB_Koster_Slater"
         case default ! sp3d5
            call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21, vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3), dH1)	! module "TB_Koster_Slater"
            call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3), dS1)	! module "TB_Koster_Slater"
         end select
         dH(k1,:,:) = dH1(:,:)
         dS(k1,:,:) = dS1(:,:)
      enddo
   enddo

   ! Kind of atom:
   KOA1 => Scell(NSC)%MDatoms(i)%KOA
   if (TB(KOA1,KOA1)%include_3body) then  ! only if user defined it to include 3-body terms
      ! Not ready yet
   endif

   nullify(r, j, KOA1)
end subroutine dHopping_Press_3TB




!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Multiplication of the two Koster-Slater angular functions (Mis*Mjs):

subroutine get_Mjs_factors(basis_ind, Scell, M_lmn, Mjs, numpar)
   integer, intent(in) :: basis_ind ! size of the basis set
   type(Super_cell), intent(in), target :: Scell  ! supercell with all the atoms as one object
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   real(8), dimension(:,:,:), allocatable, intent(inout) :: Mjs ! K-S multiplication with dummy arguments for radial function
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   !------------------------
   integer :: atom_2, j, i, Nsiz, nat
   integer, pointer :: m
   real(8), pointer :: x1, y1, z1, r1
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   nat = Scell%Na ! number of atoms
   ! Find the number of orbitals:
   Nsiz = identify_DFTB_orbitals_per_atom(basis_ind)  ! module "TB_DFTB"

   if (.not.allocated(Mjs)) allocate( Mjs(nat,nat,Nsiz) )
   Mjs = 0.0d0 ! to start with

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = nat
   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, nat	! atom #1
      m => Scell%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell%Near_neighbor_list(i,atom_2) ! atom #2
         Mjs(i,j,1) = 1.0d0  ! s-s

         if (basis_ind > 0) then ! p3 orbitals:
            Mjs(i,j,2) = -M_lmn(1,i,j)   ! px-s
            Mjs(i,j,3) = -M_lmn(2,i,j)   ! py-s
            Mjs(i,j,4) = -M_lmn(3,i,j)   ! pz-s

            if (basis_ind > 1) then ! d5 orbitals (functions from module "TB_Koster_Slater"):
               Mjs(i,j,5) = t_s_dab(M_lmn(1,i,j),M_lmn(2,i,j),1.0d0)        ! dxy-s
               Mjs(i,j,6) = t_s_dab(M_lmn(1,i,j),M_lmn(3,i,j),1.0d0)        ! dxz-s
               Mjs(i,j,7) = t_s_dab(M_lmn(2,i,j),M_lmn(3,i,j),1.0d0)        ! dyz-s
               Mjs(i,j,8) = t_s_dx2_y2(M_lmn(1,i,j),M_lmn(2,i,j),1.0d0)     ! (dx2-y2)-s
               Mjs(i,j,9) = t_s_dz2_r2(M_lmn(1,i,j),M_lmn(2,i,j),M_lmn(3,i,j),1.0d0)    ! (d3z2-r2)-s
            endif
         endif
      enddo
   enddo
   error_part = 'Error in 3TB: get_Mjs_factors:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Mjs', Mjs) ! module "MPI_subroutines"

#else ! use OpenMP instead
!$omp parallel
!$omp do private(i, j, m, atom_2)
   do i = 1, nat	! atom #1
      m => Scell%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell%Near_neighbor_list(i,atom_2) ! atom #2

!          x1 => Scell%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
!          y1 => Scell%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
!          z1 => Scell%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
!          r1 => Scell%Near_neighbor_dist(i,atom_2,4)	! at this distance, R

         Mjs(i,j,1) = 1.0d0  ! s-s

         if (basis_ind > 0) then ! p3 orbitals:
            ! [ D 2 ]
            Mjs(i,j,2) = -M_lmn(1,i,j)   ! px-s
            Mjs(i,j,3) = -M_lmn(2,i,j)   ! py-s
            Mjs(i,j,4) = -M_lmn(3,i,j)   ! pz-s

            ! [ D 1 ]
!             Mjs(i,j,2) = M_lmn(1,i,j)   ! px-s  test
!             Mjs(i,j,3) = M_lmn(2,i,j)   ! py-s
!             Mjs(i,j,4) = M_lmn(3,i,j)   ! pz-s

            if (basis_ind > 1) then ! d5 orbitals (functions from module "TB_Koster_Slater"):
               Mjs(i,j,5) = t_s_dab(M_lmn(1,i,j),M_lmn(2,i,j),1.0d0)        ! dxy-s
               Mjs(i,j,6) = t_s_dab(M_lmn(1,i,j),M_lmn(3,i,j),1.0d0)        ! dxz-s
               Mjs(i,j,7) = t_s_dab(M_lmn(2,i,j),M_lmn(3,i,j),1.0d0)        ! dyz-s
               Mjs(i,j,8) = t_s_dx2_y2(M_lmn(1,i,j),M_lmn(2,i,j),1.0d0)     ! (dx2-y2)-s
               Mjs(i,j,9) = t_s_dz2_r2(M_lmn(1,i,j),M_lmn(2,i,j),M_lmn(3,i,j),1.0d0)    ! (d3z2-r2)-s
            endif
         endif

      enddo
   enddo
!$omp enddo
!$omp end parallel
#endif
   nullify(m)
end subroutine get_Mjs_factors



!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Tools for radial function calculation within 3TB model:


subroutine get_Laguerres(r_given, rc, Laguer)
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), intent(in) :: rc     ! rescaling coefficient for the distance
   real(8), dimension(6), intent(out) :: Laguer ! polynomials up to 6
   real(8) :: d, ad
   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au * rc   ! [A] -> [Bohr]
   ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]

   call Laguerre_up_to_6(ad, Laguer)   ! module "Algebra_tools"
end subroutine get_Laguerres


pure function Laguerre_exponent(r_given, rc) result(exp_L)
   real(8) exp_L
   real(8), intent(in) :: r_given   ! distance [A]
   real(8), intent(in) :: rc     ! rescaling coefficient for the distance
   real(8) :: d, ad

   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au * rc  ! [A] -> [Bohr]
   ad = m_a * d ! scaling parameter [2]

   ! Get the exponential term:
   exp_L = exp(-0.5d0 * ad)      ! [2]
end function Laguerre_exponent




subroutine get_d_Laguerres(r_given, rc, Laguer)
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), intent(in) :: rc     ! rescaling coefficient for the distance
   real(8), dimension(6), intent(out) :: Laguer ! derivatives of the polynomials up to 6
   real(8) :: d, ad
   ! normalized distance, since coefficients are fitted in this units [1]:
   d = m_a * g_A2au * rc  ! [A] -> [Bohr]
   ad = r_given * d ! scaling parameter [2]
   ! dL/dd:
   call d_Laguerre_up_to_6(ad, Laguer)   ! module "Algebra_tools"
   ! dL/dx = dL/dd * dd/dx:
   Laguer = d * Laguer
end subroutine get_d_Laguerres


pure function d_Laguerre_exponent(r_given, rc) result(exp_L)
   real(8) exp_L
   real(8), intent(in) :: r_given   ! distance [A]
   real(8), intent(in) :: rc     ! rescaling coefficient for the distance
   real(8) :: d, ad

   ! normalized distance, since coefficients are fitted in this units [1]:
   d = m_a * g_A2au * rc  ! [A] -> [Bohr]
   ad = r_given * d ! scaling parameter [2]

   ! Get the exponential term:
   exp_L = -(0.5d0 * d) * exp(-0.5d0 * ad)      ! [2]
end function d_Laguerre_exponent



function radial_function_3TB(r_given, fx_coefs, Laguer, sh_ind, num_lag, rc, no_exp) result(f_out)
   real(8) f_out    ! value found in the array and interpolated
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(:,:), intent(in) :: fx_coefs    ! coefficients before Leguerres in 2-body Hamiltonian
   real(8), dimension(:), intent(in) :: Laguer   ! precalculated Laguerre polinomials
   integer, intent(in) :: sh_ind, num_lag    ! which shell, how many polynomials
   real(8), intent(in), optional :: rc     ! rescaling coefficient for the distance
   logical, intent(in), optional :: no_exp   ! indicate if Laguer = Laguerre * exp(-a*r_IJ), or not
   !---------------------
   real(8) :: d, ad, exp_ad
   integer :: i
   logical :: there_is_exp

   there_is_exp = .true.  ! by default, there is exp. provided inside of Laguer
   if (present(no_exp) .and. present(rc)) then
      if (no_exp) then
         there_is_exp = .false.   ! if it is not provided, recalculate it
      endif
   endif

   ! Construct the function:
   f_out = SUM( fx_coefs(sh_ind, 1:num_lag) * Laguer(1:num_lag) )

   ! In case the Laguerre polynomials did not contain exponent, calculate it:
   if (.not.there_is_exp) then
      ! Get the exponential term:
      exp_ad = Laguerre_exponent(r_given, rc) ! below

      ! Include it into the solution:
      f_out = exp_ad * f_out
   endif
end function radial_function_3TB



! function d_radial_function_3TB(r_given, M_Vij, fx_coefs, d_Laguer, sh_ind, num_lag) result(f_out)
!    real(8) f_out    ! value found in the array and interpolated
!    real(8), intent(in) :: r_given   ! [A] given radial coordinate
!    real(8), intent(in) :: M_Vij     ! radial function
!    real(8), dimension(:,:), intent(in) :: fx_coefs    ! coefficients before Leguerres in 2-body Hamiltonian
!    real(8), dimension(:), intent(in) :: d_Laguer   ! precalculated derivatives of the Laguerre polinomials
!    integer, intent(in) :: sh_ind, num_lag    ! which shell, how many polynomials
!    !---------------------
!    real(8) :: d, ad, exp_ad
!    integer :: i
!
!    ! normalized distance, since coefficients are fitted in this units [1]:
!    d = r_given * g_A2au   ! [A] -> [Bohr]
!    ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]
!
!    ! Get the exponential term:
!    exp_ad = exp(-0.5d0 * ad)
!
!    ! Construct the function:
!    f_out = 0.0d0  ! to start with
!    ! Skip i=1, since dL_0/dr = 0
!    do i = 2, num_lag ! for all polynomials required
!       f_out = f_out + fx_coefs(sh_ind, i) * d_Laguer(i)
!    enddo
!    f_out = exp_ad * f_out - (m_a*0.5d0/g_a0)*M_Vij
! end function d_radial_function_3TB





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
