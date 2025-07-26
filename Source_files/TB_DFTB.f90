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
! This module contains subroutines to deal with TB hamiltonian in the DFTB parametrization: http://www.dftb.org/
! The Koster-Slater files with numerical sets of radial Hamiltonian and Overlap integrals: http://www.dftb.org/parameters/
! [1] skaformat.pdf file from http://www.dftb.org/

MODULE TB_DFTB

use Universal_constants
use Objects
use TB_Koster_Slater
use Little_subroutines, only : linear_interpolation, Fermi_function, d_Fermi_function, Find_in_array_monoton
use Electron_tools, only : find_band_gap
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, test_orthogonalization_c, Loewdin_Orthogonalization, &
                     Loewdin_Orthogonalization_c
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity
use Atomic_tools, only : Reciproc_rel_to_abs

#ifdef MPI_USED
   use MPI_subroutines, only : do_MPI_Allreduce
#endif


implicit none
PRIVATE

public :: Construct_Vij_DFTB, construct_TB_H_DFTB, get_Erep_s_DFTB, get_dHij_drij_DFTB, &
          Attract_TB_Forces_Press_DFTB, dErdr_s_DFTB, dErdr_Pressure_s_DFTB, Complex_Hamil_DFTB, &
          identify_DFTB_orbitals_per_atom, identify_DFTB_basis_size, Get_overlap_S_matrix_DFTB, &
          Hopping_DFTB, get_Erep_s_DFTB_no, dErdr_s_DFTB_no, dErdr_Pressure_s_DFTB_no

 contains
 

subroutine Construct_Vij_DFTB(numpar, TB, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(TB_H_DFTB), dimension(:,:), intent(in) :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable :: M_Vij	! matrix of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dVij	! matrix of derivatives of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_SVij	! matrix of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dSVij	! matrix of derivatives of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   !----------------------------------
   real(8) :: x, y, z, r, sx, sy, sz, r1, fcut, d_fcut, Fermi, dFermi
   integer :: i, j, atom_2, ki, N_bs, ihop
   real(8), pointer :: rm
   integer, pointer :: nat, m, KOA1, KOA2
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part


   nat => Scell(NSC)%Na	! number of atoms in the supercell
   ! number of hopping integrals for this basis set in DFTB:
   N_bs = identify_DFTB_basis_size(numpar%basis_size_ind)  ! below
   
   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,N_bs))	! each pair of atoms, all  V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,N_bs))	! each pair of atoms, all  dV functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,N_bs))	! each pair of atoms, all S functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,N_bs))	! each pair of atoms, all dS functions
    
   M_Vij = 0.0d0
   M_dVij = 0.0d0
   M_SVij = 0.0d0
   M_dSVij = 0.0d0
   
   ! Construct matrix of all the radial functions for each pair of atoms:
#ifdef MPI_USED   ! use the MPI version [tested]
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = nat
   ! Do the cycle (parallel) calculations:
   AT1:do j = Nstart, Nend, N_incr  ! each process does its own part
   m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         KOA2 => Scell(NSC)%MDatoms(j)%KOA  ! Correct order, checked by cohesive energy minimum
         KOA1 => Scell(NSC)%MDatoms(i)%KOA
         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! All radial functions for Hamiltonian:
         M_Vij(j,i,1) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 1)   ! (s s sigma)
         M_SVij(j,i,1) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 1) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
            M_Vij(j,i,2) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 2)   ! (s p sigma)
            M_SVij(j,i,2) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 2) ! (s p sigma)
            M_Vij(j,i,3) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 4)   ! (p p sigma)
            M_SVij(j,i,3) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 4) ! (p p sigma)
            M_Vij(j,i,4) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 5)   ! (p p pi)
            M_SVij(j,i,4) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 5) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_Vij(j,i,ihop) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, ihop)
               M_SVij(j,i,ihop) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, ihop)
            enddo
         endselect

         KOA1 => Scell(NSC)%MDatoms(i)%KOA    ! testing
         KOA2 => Scell(NSC)%MDatoms(j)%KOA
         ! All derivatives of the radial functions and the radial functions for Overlap matrix:
         M_dVij(j,i,1) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 1)   ! (s s sigma)
         M_dSVij(j,i,1) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 1) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
            M_dVij(j,i,2) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 2)   ! (s p sigma)
            M_dSVij(j,i,2) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 2) ! (s p sigma)
            M_dVij(j,i,3) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 4)   ! (p p sigma)
            M_dSVij(j,i,3) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 4) ! (p p sigma)
            M_dVij(j,i,4) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 5)   ! (p p pi)
            M_dSVij(j,i,4) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 5) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_dVij(j,i,ihop) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, ihop)
               M_dSVij(j,i,ihop) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, ihop)
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
   error_part = 'Error in Construct_Vij_DFTB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_Vij', M_Vij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_dVij', M_dVij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_SVij', M_SVij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'M_dSVij', M_dSVij) ! module "MPI_subroutines"

#else
!$omp PARALLEL
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, ihop, Fermi, dFermi)
   AT1:do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one  
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms

!          KOA1 => Scell(NSC)%MDatoms(j)%KOA    ! WRONG ORDER!
!          KOA2 => Scell(NSC)%MDatoms(i)%KOA
         KOA2 => Scell(NSC)%MDatoms(j)%KOA  ! Correct order, checked by cohesive energy minimum
         KOA1 => Scell(NSC)%MDatoms(i)%KOA
         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! All radial functions for Hamiltonian:
         M_Vij(j,i,1) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 1)   ! (s s sigma)
         M_SVij(j,i,1) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 1) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
            M_Vij(j,i,2) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 2)   ! (s p sigma)
            M_SVij(j,i,2) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 2) ! (s p sigma)
            M_Vij(j,i,3) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 4)   ! (p p sigma)
            M_SVij(j,i,3) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 4) ! (p p sigma)
            M_Vij(j,i,4) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 5)   ! (p p pi)
            M_SVij(j,i,4) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 5) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_Vij(j,i,ihop) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, ihop)
               M_SVij(j,i,ihop) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, ihop)
            enddo
         endselect

         KOA1 => Scell(NSC)%MDatoms(i)%KOA    ! testing
         KOA2 => Scell(NSC)%MDatoms(j)%KOA
         ! All derivatives of the radial functions and the radial functions for Overlap matrix:
         M_dVij(j,i,1) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 1)   ! (s s sigma)
         M_dSVij(j,i,1) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 1) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
            M_dVij(j,i,2) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 2)   ! (s p sigma)
            M_dSVij(j,i,2) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 2) ! (s p sigma)
            M_dVij(j,i,3) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 4)   ! (p p sigma)
            M_dSVij(j,i,3) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 4) ! (p p sigma)
            M_dVij(j,i,4) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 5)   ! (p p pi)
            M_dSVij(j,i,4) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 5) ! (p p pi)
         case (2)    ! sp3d5
            do ihop = 2, 10
               M_dVij(j,i,ihop) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, ihop)
               M_dSVij(j,i,ihop) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, ihop)
            enddo
         endselect
         
         ! Add Fermi-like smoothing approach to zero:
         Fermi = Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         dFermi = d_Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         M_dVij(j,i,:) = M_dVij(j,i,:)*Fermi + M_Vij(j,i,:)*dFermi
         M_dSVij(j,i,:) = M_dSVij(j,i,:)*Fermi + M_SVij(j,i,:)*dFermi
         M_Vij(j,i,:) = M_Vij(j,i,:)*Fermi
         M_SVij(j,i,:) = M_SVij(j,i,:)*Fermi

!          ! Testing:
!          M_dVij = 0.0d0
!          M_dSVij = 0.0d0
!          M_Vij = 0.0d0
!          M_SVij = 0.d00

      enddo AT2
   enddo AT1
!$omp end do
   nullify(m, KOA1, KOA2)	! clean up for each thread
!$omp END PARALLEL
#endif
   nullify(rm, nat, m, KOA1, KOA2)	! clean up at the end
end subroutine Construct_Vij_DFTB


function DFTB_radial_function(r_given, r_grid, param_array, ind_array, ind_interpl) result(f_out)
    real(8) :: f_out    ! value found in the array and interpolated
    real(8), intent(in) :: r_given  ! [A] given radial coordinate
    real(8), dimension(:), intent(in), target :: r_grid  ! [A] radial grid
    real(8), dimension(:,:), intent(in) :: param_array    ! [eV] array of parameters to extract a value from
    integer, intent(in) :: ind_array    ! which element we need
    integer, intent(in), optional :: ind_interpl    ! which interpolation scheme to use
    !-------------------------------------
    integer :: i_array, Nsiz
    real(8), pointer :: dr
    dr => r_grid(1) ! grid step in a uniformly distributed grid
    Nsiz = size(r_grid)
    i_array = CEILING(r_given/dr)   ! index of the nearest radial grid point
    if (i_array > Nsiz) then    ! return the last value
       f_out = param_array(Nsiz, ind_array)
    else
       !call linear_interpolation(r_grid, param_array(:,ind_array), r_given, f_out, i_array, x0=0.0d0, y0=0.0d0)   ! module "Little_subroutines"
       call linear_interpolation(r_grid, param_array(:,ind_array), r_given, f_out, i_array)   ! module "Little_subroutines"

       if ( (f_out < min(param_array(i_array-1,ind_array), param_array(i_array,ind_array))) .or. &
            (f_out > max(param_array(i_array-1,ind_array), param_array(i_array,ind_array))) ) then
          print*, 'DFTB_radial_function', r_given, r_grid(i_array-1), r_grid(i_array), ':', f_out, param_array(i_array-1,ind_array), param_array(i_array,ind_array)
       endif
    endif
    nullify(dr)
end function DFTB_radial_function
 
 
function d_DFTB_radial_function(r_given, r_grid, param_array, ind_array, ind_interpl) result(f_out)
    real(8) :: f_out    ! value found in the array and interpolated
    real(8), intent(in) :: r_given  ! [A] given radial coordinate
    real(8), dimension(:), intent(in), target :: r_grid  ! [A] radial grid
    real(8), dimension(:,:), intent(in) :: param_array    ! [eV] array of parameters to extract a value from
    integer, intent(in) :: ind_array    ! which element we need
    integer, intent(in), optional :: ind_interpl    ! which interpolation scheme to use
    !-------------------------------------
    integer :: i_array, Nsiz
    real(8), pointer :: dr
    dr => r_grid(1) ! grid step in a uniformly distributed grid
    Nsiz = size(r_grid)
    i_array = CEILING(r_given/dr)   ! index of the nearest radial grid point
    if (i_array > Nsiz) then    ! return the last value
       f_out = 0.0d0
    elseif (i_array == 1) then   ! too close atoms
       f_out = (param_array(i_array+1,ind_array) - param_array(i_array,ind_array))/dr   ! extrapolate with the same derivative as between first and second points
    else
       f_out = (param_array(i_array,ind_array) - param_array(i_array-1,ind_array))/dr   ! numerical derivative, linear
    endif
    nullify(dr)
 end function d_DFTB_radial_function

 
! Tight Binding Hamiltonian within DFTB parametrization:
subroutine construct_TB_H_DFTB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err, scc, H_scc_0, H_scc_1)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! materil parameters
   type(TB_H_DFTB), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
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
      call Hamil_tot_DFTB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err, scc, H_scc_0, H_scc_1)  ! see below
   else
      call Hamil_tot_DFTB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err)  ! see below
   endif


!    ! Test (comment out for release):
!    call test_nonorthogonal_solution(Scell(NSC)) ! module "TB_NRL"

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_DFTB



subroutine Hamil_tot_DFTB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err, scc, H_scc_0, H_scc_1)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_DFTB), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
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
   n_orb =  identify_DFTB_orbitals_per_atom(numpar%basis_size_ind)  ! size of the basis set per atom, below
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

   !-----------------------------------

   WNTSCC:if (do_scc == 0) then   ! no scc calculations, construct regular (zero-order) Hamiltonian:

#ifdef MPI_USED   ! use the MPI version
      if (.not.allocated(Hij1)) allocate(Hij1(n_orb,n_orb))
      if (.not.allocated(Sij1)) allocate(Sij1(n_orb,n_orb))
      N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
      Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
      Nend = nat
      ! Do the cycle (parallel) calculations:
      do j = Nstart, Nend, N_incr  ! each process does its own part
      !do j = 1,nat	! all atoms
         m => Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            endif

            IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
               KOA1 => Scell(NSC)%MDatoms(j)%KOA
               KOA2 => Scell(NSC)%MDatoms(i)%KOA
               ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
               ! Contruct a block-hamiltonian:
               call Hamilton_one_DFTB(numpar%basis_size_ind, j, i, TB_Hamil(KOA1,KOA2), Hij1, M_Vij, M_lmn)   ! below

               ! Construct overlap matrix for this pair of atoms:
               call Get_overlap_S_matrix_DFTB(numpar%basis_size_ind, j, i, Sij1, M_SVij, M_lmn)  ! below

               do j1 = 1,n_orb ! all orbitals
                  l = (j-1)*n_orb+j1
                  do i1 = 1,n_orb ! all orbitals
                     k = (i-1)*n_orb+i1
                     Hij(k,l) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                     Sij(k,l) = Sij1(i1,j1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
                     if (ABS(Sij(k,l)) <= epsylon) Sij(k,l) = 0.0d0
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
            if (i < j) then ! it's a new pair of atoms, calculate everything
               do j1 = 1,n_orb ! all orbitals
                  l = (j-1)*n_orb+j1
                  do i1 = 1,n_orb ! all orbitals
                     k = (i-1)*n_orb+i1
                     Hij(k,l) = Hij(l,k)
                     Sij(k,l) = Sij(l,k)
                  enddo ! i1
               enddo ! j1
            endif
         enddo ! j
      enddo ! i

      ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
      error_part = 'Error in Hamil_tot_DFTB:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Hij', Hij) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Sij', Sij) ! module "MPI_subroutines"


#else    ! OpenMP to use instead
      ! 1) Construct non-orthogonal Hamiltonian H and Overlap matrix S in 2 steps:
!$omp parallel private(j, m, atom_2, i, KOA1, KOA2, j1, l, i1, k, Hij1, Sij1)
      if (.not.allocated(Hij1)) allocate(Hij1(n_orb,n_orb))
      if (.not.allocated(Sij1)) allocate(Sij1(n_orb,n_orb))
!$omp do
      do j = 1,nat	! all atoms
         m => Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            endif
         
            IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
               KOA1 => Scell(NSC)%MDatoms(j)%KOA
               KOA2 => Scell(NSC)%MDatoms(i)%KOA
               ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
               ! Contruct a block-hamiltonian:
               call Hamilton_one_DFTB(numpar%basis_size_ind, j, i, TB_Hamil(KOA1,KOA2), Hij1, M_Vij, M_lmn)   ! below
            
               ! Construct overlap matrix for this pair of atoms:
               call Get_overlap_S_matrix_DFTB(numpar%basis_size_ind, j, i, Sij1, M_SVij, M_lmn)  ! below

               do j1 = 1,n_orb ! all orbitals
                  l = (j-1)*n_orb+j1
                  do i1 = 1,n_orb ! all orbitals
                     k = (i-1)*n_orb+i1
                     Hij(k,l) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                     Sij(k,l) = Sij1(i1,j1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
!                    Hij(k,l) = Hij1(j1,i1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
!                    Sij(k,l) = Sij1(j1,i1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
                     if (ABS(Sij(k,l)) <= epsylon) Sij(k,l) = 0.0d0
                  enddo ! i1
               enddo ! j1
            endif IJ
         enddo ! j
      enddo ! i
!$omp end do
      deallocate(Hij1, Sij1)
      nullify(KOA1, KOA2, m)
!$omp end parallel

      ! b) Construct lower triangle - use symmetry:
!$omp parallel
!$omp do  private(j, m, atom_2, i, j1, l, i1, k)
      do j = 1,nat	! all atoms
         m => Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 1,m ! do only for atoms close to that one
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            if (i < j) then ! it's a new pair of atoms, calculate everything
               do j1 = 1,n_orb ! all orbitals
                  l = (j-1)*n_orb+j1
                  do i1 = 1,n_orb ! all orbitals
                     k = (i-1)*n_orb+i1
                     Hij(k,l) = Hij(l,k)
                     Sij(k,l) = Sij(l,k)
                  enddo ! i1
               enddo ! j1
            endif
         enddo ! j
      enddo ! i
!$omp end do
      nullify(m)
!$omp end parallel
#endif


   ! 2)    ! Save the non-orthogonalized Hamiltonian:
      Scell(NSC)%H_non = Hij		! nondiagonalized Hamiltonian
      Scell(NSC)%Sij = Sij		! save Overlap matrix

   else WNTSCC ! scc cycle, construct only use second order scc correction:
      Hij = H_scc_0 + H_scc_1 ! zero and second order scc contributions into Hamiltonian
      Scell(NSC)%H_non = Hij  ! save new nondiagonalized Hamiltonian
      Sij = Scell(NSC)%Sij
   endif WNTSCC


   !-----------------------------------
   ! 3) Orthogonalize the Hamiltonian using Lowedin procidure
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization(numpar, Nsiz, Sij, Hij, Err, Scell(NSC)%eigen_S)    ! module "TB_NRL"
   
   Scell(NSC)%Hij = Hij ! save orthogonalized but non-diagonalized Hamiltonian
   
   forall (i = 1:size(Hij,1),  j = 1:size(Hij,2), (ABS(Hij(i,j)) < 1.0d-10))
      Hij(i,j) = 0.0d0
   endforall
   !-----------------------------------
   ! 4) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, numpar%MPI_param, check_M=.true.)
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, numpar%MPI_param, use_DSYEV=.false.)
   call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, numpar%MPI_param, use_DSYEV=.true.) ! module "Algebra_tools"
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
   ! Normalize eigenvectors to | <n|n> |^2 = 1:
!    do j = 1, Nsiz
! !       temp = DSQRT(SUM( Scell(NSC)%Ha(:,j) * Scell(NSC)%Ha(:,j) ))
! !       Scell(NSC)%Ha(:,j) =  Scell(NSC)%Ha(:,j) / temp
!       temp = DSQRT(SUM( Scell(NSC)%Ha(j,:) * Scell(NSC)%Ha(j,:) )) 
!       Scell(NSC)%Ha(j,:) =  Scell(NSC)%Ha(j,:) / temp
!    enddo
   !-----------------------------------
!    ! An example of the solution of the linear eigenproblem with LAPACK subroutines:
!    call dpotrf('U', Nsiz, Sij, Nsiz, j)
!    call dsygst(1, 'U', Nsiz, Hij, Nsiz, Sij, Nsiz,  j)
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, numpar%MPI_param)
!    call mkl_matrix_mult('N', 'N', Sij, Hij, Scell(NSC)%Ha)	! module "Algebra_tools"

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
#ifdef MPI_USED   ! use the MPI version [tested]
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
!                      if (abs(i) > nat*n_orb) print*, 'TEST i', i, j, l, k
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
!                         write(*,'(i4,i4,e,e,e)') k,l, matter%PRRx(k,l), matter%PRRy(k,l), matter%PRRz(k,l)
                  enddo ! i1
               enddo ! j1
            endif ! (i > 0)
         enddo ! j
      enddo ! i
      ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
      error_part = 'Error in Hamil_tot_DFTB:'
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
!                      if (abs(i) > nat*n_orb) print*, 'TEST i', i, j, l, k
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
!                         write(*,'(i4,i4,e,e,e)') k,l, matter%PRRx(k,l), matter%PRRy(k,l), matter%PRRz(k,l)
                  enddo ! i1
               enddo ! j1
            endif ! (i > 0)
         enddo ! j
      enddo ! i
      !$omp end do
      nullify(m, x, y, z)
      !$omp end parallel
#endif

      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0
!       temp = 1.0d0 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
      
      Scell(NSC)%PRRx = Scell(NSC)%PRRx * temp
      Scell(NSC)%PRRy = Scell(NSC)%PRRy * temp
      Scell(NSC)%PRRz = Scell(NSC)%PRRz * temp
   endif ! (numpar%optic_model .EQ. 3)
   
   nullify(KOA1, KOA2, m, x, y, z)
   deallocate(Hij, Sij)
end subroutine Hamil_tot_DFTB



subroutine Hamilton_one_DFTB(basis_ind, i, j, TB, Hij, M_Vij, M_lmn)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   integer, intent(in) :: i, j
   type(TB_H_DFTB), intent(in), target :: TB	! all tight binding parameters
   real(8), dimension(:,:), intent(out) :: Hij  ! hamiltonian, all orbitals in sp3d5 basis set
   real(8), dimension(:,:,:), intent(in) :: M_Vij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   !---------------------------------------/
   integer(4) i1
   if (i == j) then ! Onsite contributions
      Hij = 0.0d0   ! Nondiagonals are zeros
      Hij(1,1) = TB%Es ! s orbital
      if (basis_ind > 0) then ! p3 orbitals:
         do i1 = 2, 4
            Hij(i1,i1) = TB%Ep
         enddo
      endif
      if (basis_ind > 1) then ! d5 orbitals:
         do i1 = 5, 9
            Hij(i1,i1) = TB%Ed
         enddo
      endif
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:
      call Hopping_DFTB(basis_ind, M_Vij(i,j,:), M_Vij(j,i,:), M_lmn(:,i,j), Hij)	! subroutine below
      !print*, '(i,j):', i, j
      !print*, M_Vij(i,j,:)
      !print*, '(j,i):', i, j
      !print*, M_Vij(j,i,:)
      !pause 'Hamilton_one_DFTB'
   endif
   
end subroutine Hamilton_one_DFTB


subroutine Hopping_DFTB(basis_ind, M_Vij12, M_Vij21, M_lmn, ts)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   real(8), dimension(:), intent(in) :: M_Vij12, M_Vij21	! Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:), intent(in) :: M_lmn	! cosine directions
   real(8), dimension(:,:), intent(out) :: ts	! overlap integerals [eV]
   !=============================================
   real(8), dimension(:), allocatable :: vec_M_SVij12, vec_M_SVij21
   ! Construct the overlap integrals including angular part 
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
!       call KS_sp3_hetero_TEST(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   case default ! for sp3d5 basis set:
      allocate(vec_M_SVij12(10))
      vec_M_SVij12(:) = M_Vij12(:)
      allocate(vec_M_SVij21(10))
      vec_M_SVij21(:) = M_Vij21(:)
      call KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   end select
   if (allocated(vec_M_SVij12)) deallocate(vec_M_SVij12)
   if (allocated(vec_M_SVij21)) deallocate(vec_M_SVij21)
end subroutine Hopping_DFTB


subroutine Get_overlap_S_matrix_DFTB(basis_ind, i, j, Sij, M_SVij, M_lmn)
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
         call KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij)	! module "TB_Koster_Slater"
      case default ! for sp3d5 basis set:
         allocate(vec_M_SVij12(10))
         vec_M_SVij12(:) = M_SVij(i,j,:)
         allocate(vec_M_SVij21(10))
         vec_M_SVij21(:) = M_SVij(j,i,:)
         call KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij)	! module "TB_Koster_Slater"
      end select
   else	! it is the same atom
      Sij = 0.0d0
      forall (k=1:size(Sij,1)) Sij(k,k)=1.0d0
   endif
   if (allocated(vec_M_SVij12)) deallocate(vec_M_SVij12)
   if (allocated(vec_M_SVij21)) deallocate(vec_M_SVij21)
end subroutine Get_overlap_S_matrix_DFTB


!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Derivatives of the attractive part of TB:

! Subroutine for derivative of the Hamiltonian:
subroutine get_dHij_drij_DFTB(numpar, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
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
      call get_forces_DFTB(k, numpar, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)	! see below
   enddo ATOMS

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_dHij_drij_DFTB:'
   do k = 1, nat
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'MDatoms(k)%forces%att(:)', Scell(NSC)%MDatoms(k)%forces%att(:)) ! module "MPI_subroutines"
   enddo

#else ! use OpenMP instead
   !$omp PARALLEL private(k) 
   !$omp do
   ATOMS:do k = 1, nat	! forces for all atoms
      Scell(NSC)%MDatoms(k)%forces%att(:) = 0.0d0	! just to start
      call get_forces_DFTB(k, numpar, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)	! see below
   enddo ATOMS
   !$omp end do 
   !$omp end parallel
#endif
end subroutine get_dHij_drij_DFTB



subroutine get_forces_DFTB(k, numpar, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)
   integer, intent(in) :: k	! forces for this atom
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of the supercell
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   !-----------------------------
   integer :: j, m, atom_2, i, j1, i1, nat, Nsiz, n_orb, i4, j4, j4j1, i4i1
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
   
   ! 1) Construct the derivatives of the Hamiltonian and Overlap matrix in 2 steps:
   ! a) Construct upper triangle - calculate each element:
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
            call d_Hamilton_one_DFTB(numpar%basis_size_ind, k, Scell, NSC, i, j, atom_2, dH1, M_Vij, M_dVij, M_lmn, dS1, M_SVij, M_dSVij) ! this calls the block-hamiltonian

            do j1 = 1,n_orb	! all orbitals
               j4j1 = j4+j1
               do i1 = 1,n_orb	! all orbitals
                  i4i1 = i4+i1
                  dH(:,i4i1,j4j1) = dH1(:,i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  dS(:,i4i1,j4j1) = dS1(:,i1,j1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
!                   dH(:,i4i1,j4j1) = dH1(:,j1,i1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40) ! WRONG
!                   dS(:,i4i1,j4j1) = dS1(:,j1,i1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
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
               j4j1 = j4+j1
               do i1 = 1,n_orb	! all orbitals
                  i4i1 = i4+i1
                  dH(:,i4i1,j4j1) = dH(:,j4j1,i4i1)
                  dS(:,i4i1,j4j1) = dS(:,j4j1,i4i1)
               enddo ! i1
            enddo ! j1
         endif ! (j < i)
      enddo ! i
   enddo ! atom_2

   ! 2) Calculate the forces form the derivatives and the eigervectors:
   call Attract_TB_forces_DFTB(Aij, Aij_x_Ei, dH, dS, Scell, NSC, Scell(NSC)%MDatoms(k)%forces%att(:), n_orb)

   deallocate(dH, dS)
end subroutine get_forces_DFTB



subroutine Attract_TB_forces_DFTB(Aij, Aij_x_Ei, dH, dS, Scell, NSC, Eelectr_s, n_orb)
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
   n_too = 0   ! to srtart with

   i2 = 0
   ste = 1
   do i = 1, n	! all orbitals
       if (i .GE. ste) then
          i2 = i2 + 1	! number atom corresponding to these set of orbitals
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

             if (maxval(ABS(Eelectr_s(:))) .GE. 1.0d6) n_too = n_too + 1
!              if (maxval(ABS(Eelectr_s(:))) .GE. 1.0d5) then
!                 write(*,'(a)') 'Trouble in subroutine Attract_TB_forces_DFTB, too large attractive force:'
!                 write(*,*) i, j
!                 write(*,*) dH(1,i,j:j_norb)
!                 write(*,*) Aij(i,j:j_norb)
!                 write(*,'(e25.16,e25.16,e25.16)') Eelectr_s(:)
!              endif
          endif
       enddo
   enddo

   if (n_too > 0) then
      write(*,'(a)') 'Trouble in subroutine Attract_TB_forces_DFTB, too large attractive force'
      write(*,'(a,i)') 'For elements: ', n_too
   endif

   nullify(m, j1)
end subroutine Attract_TB_forces_DFTB


!ddddddddddddddddddddddddddddddddddddddddddddddddddd
! Derivatives:
subroutine d_Hamilton_one_DFTB(basis_ind, k, Scell, NSC, i, j, atom_2, dH, M_Vij, M_dVij, M_lmn, dS, M_SVij, M_dSVij)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC, i, j, atom_2, k
   real(8), dimension(:,:,:), intent(out) :: dH, dS  ! hamiltonian, all orbitals in sp3d5 basis set
   real(8), dimension(:,:,:), intent(in), target :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines and derivatives
   real(8), dimension(:,:,:), intent(in), target :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms for S-matrix, all orbitals, and derivatives
   !---------------------------------------
   real(8) :: drij_dsk(3)
   real(8), dimension(9) :: M_dlmn	! dl/dx, dl/dy, dl/dz, dm/dx, dm/dy, dm/dz, dn/dx, dn/dy, dn/dz
   real(8), pointer :: x1, y1, z1, r1
   integer :: n_orb, n_overlap
   real(8), dimension(:), allocatable :: vec_M_Vij12, vec_M_Vij21, vec_M_SVij12, vec_M_SVij21
   real(8), dimension(:), allocatable :: vec_M_dVij12, vec_M_dVij21, vec_M_dSVij12, vec_M_dSVij21
   real(8), dimension(:,:), allocatable :: dH1, dS1

   if (i == j) then ! Onsite contributions are consts, so derivatives are always 0:
      dH = 0.0d0
      dS = 0.0d0
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:   
      x1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
      y1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
      z1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
      r1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)		! at this distance, R
      
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
      
      ! CHECKED CORRECT:
      drij_dsk(1) = drij_dska(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 1, .true.)	! dr_{ij}/ds_{k,x}, module "TB_Koster_Slater"
      drij_dsk(2) = drij_dska(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 2, .true.)	! dr_{ij}/ds_{k,y}, module "TB_Koster_Slater"
      drij_dsk(3) = drij_dska(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 3, .true.)	! dr_{ij}/ds_{k,z}, module "TB_Koster_Slater"

      ! CHECKED CORRECT:
      M_dlmn(1) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 1, 1, drij_dsk(1))	! dl/dsx, module "TB_Koster_Slater"
      M_dlmn(2) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 1, 2, drij_dsk(2))	! dl/dsy, module "TB_Koster_Slater"
      M_dlmn(3) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 1, 3, drij_dsk(3))	! dl/dsz, module "TB_Koster_Slater"
      M_dlmn(4) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 2, 1, drij_dsk(1))	! dm/dsx, module "TB_Koster_Slater"
      M_dlmn(5) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 2, 2, drij_dsk(2))	! dm/dsy, module "TB_Koster_Slater"
      M_dlmn(6) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 2, 3, drij_dsk(3))	! dm/dsz, module "TB_Koster_Slater"
      M_dlmn(7) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 3, 1, drij_dsk(1))	! dn/dsx, module "TB_Koster_Slater"
      M_dlmn(8) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 3, 2, drij_dsk(2))	! dn/dsy, module "TB_Koster_Slater"
      M_dlmn(9) = ddija_dskb_kd(i, j, k, x1, y1, z1, r1, Scell(NSC)%supce, 3, 3, drij_dsk(3))	! dn/dsz, module "TB_Koster_Slater"

      ! Create vectors from the slices of arrays to pass into subroutines:
!       vec_M_Vij12 = M_Vij(i,j,:)
!       vec_M_Vij21 = M_Vij(j,i,:)
!       vec_M_SVij12 = M_SVij(i,j,:)
!       vec_M_SVij21 = M_SVij(j,i,:)

      vec_M_Vij21 = M_Vij(i,j,:)    ! correct
      vec_M_Vij12 = M_Vij(j,i,:)
      vec_M_SVij21 = M_SVij(i,j,:)
      vec_M_SVij12 = M_SVij(j,i,:)

      ! Derivative along dx:
      dH1 = 0.0d0
      dS1 = 0.0d0
!       vec_M_dVij12(:) = M_dVij(i,j,:)*drij_dsk(1)
!       vec_M_dVij21(:) = M_dVij(j,i,:)*drij_dsk(1)
!       vec_M_dSVij12(:) = M_dSVij(i,j,:)*drij_dsk(1)
!       vec_M_dSVij21(:) = M_dSVij(j,i,:)*drij_dsk(1)
      vec_M_dVij21(:) = M_dVij(i,j,:)*drij_dsk(1) ! correct
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
      case default ! sp3d5
         call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21,  vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dS1)	! module "TB_Koster_Slater"
      end select
      dH(1,:,:) = dH1(:,:)
      dS(1,:,:) = dS1(:,:)
         
      ! Derivative along dy:
      dH1 = 0.0d0
      dS1 = 0.0d0
!       vec_M_dVij12(:) = M_dVij(i,j,:)*drij_dsk(2)
!       vec_M_dVij21(:) = M_dVij(j,i,:)*drij_dsk(2)
!       vec_M_dSVij12(:) = M_dSVij(i,j,:)*drij_dsk(2)
!       vec_M_dSVij21(:) = M_dSVij(j,i,:)*drij_dsk(2)

      vec_M_dVij21(:) = M_dVij(i,j,:)*drij_dsk(2)   ! correct
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
      case default ! sp3d5
         call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21,  vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dS1)	! module "TB_Koster_Slater"
      end select
      
      dH(2,:,:) = dH1(:,:)
      dS(2,:,:) = dS1(:,:)

      ! Derivative along dz:
      dH1 = 0.0d0
      dS1 = 0.0d0
!       vec_M_dVij12(:) = M_dVij(i,j,:)*drij_dsk(3)
!       vec_M_dVij21(:) = M_dVij(j,i,:)*drij_dsk(3)
!       vec_M_dSVij12(:) = M_dSVij(i,j,:)*drij_dsk(3)
!       vec_M_dSVij21(:) = M_dSVij(j,i,:)*drij_dsk(3)

      vec_M_dVij21(:) = M_dVij(i,j,:)*drij_dsk(3)   ! correct
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
      case default ! sp3d5
         call d_KS_sp3d5_hetero(vec_M_Vij12, vec_M_Vij21,  vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dS1)	! module "TB_Koster_Slater"
      end select
      dH(3,:,:) = dH1(:,:)
      dS(3,:,:) = dS1(:,:)
      
      deallocate(vec_M_Vij12, vec_M_Vij21, vec_M_SVij12, vec_M_SVij21, vec_M_dVij12, vec_M_dVij21, vec_M_dSVij12, vec_M_dSVij21, dH1, dS1)
   endif ! (i == j)
   
   ! Clean up at the end:
   nullify (x1, y1, z1, r1)
end subroutine d_Hamilton_one_DFTB




!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! Attractive forces for supercell from the derivatives of the Hamiltonian:
subroutine Attract_TB_Forces_Press_DFTB(Scell, NSC, numpar, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
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
      
      call dHamil_tot_Press_DFTB(Scell, NSC, numpar, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij, dSij)   ! below

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
      error_part = 'Error in Attract_TB_Forces_Press_DFTB:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dwr_press', dwr_press) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dS_press', dS_press) ! module "MPI_subroutines"

#else ! use OpenMP instead
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
            Scell(NSC)%SCforce%att(k,i) = SUM(dwr_press((i-1)*3+k,:)) - SUM(dS_press((i-1)*3+k,:))
!             Scell(NSC)%SCforce%att(i,k) = SUM(dwr_press((i-1)*3+k,:)) - SUM(dS_press((i-1)*3+k,:))
         enddo ! k
      enddo ! i
      deallocate(dwr_press, dS_press, dHij, dSij)
   endif
end subroutine Attract_TB_Forces_Press_DFTB



subroutine dHamil_tot_Press_DFTB(Scell, NSC, numpar, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij, dSij)
! construct the whole Hamilton matrix:
! (with respect to which Rk we take the derivatives, Appendix F of H.Jeschke PhD Thesis)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(inout) :: dHij, dSij
!    real(8), dimension(9,9,9) :: dHij1, dSij1
   real(8), dimension(:,:,:), allocatable :: dHij1, dSij1
   integer :: i, j, j1, i1, atom_2, m, nat, i2, j2
   integer :: i4, j4, norb, n_overlap
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   ! Depending on the basis set:
   n_overlap = identify_DFTB_basis_size(numpar%basis_size_ind)   ! below
   norb = identify_DFTB_orbitals_per_atom(numpar%basis_size_ind)    ! below
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
            call dHamilton_one_Press_DFTB(i, atom_2, Scell, NSC, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij1, dSij1)
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis.
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
   error_part = 'Error in dHamil_tot_Press_DFTB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dHij', dHij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dSij', dSij) ! module "MPI_subroutines"

#else ! use OpenMP instead
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
            call dHamilton_one_Press_DFTB(i, atom_2, Scell, NSC, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij1, dSij1)
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis.
            do j1 = 1,norb	! all orbitals
               j2 = j4+j1
               do i1 = 1,norb	! all orbitals
                  i2 = i4+i1
                  dHij(:,i2,j2) = dHij1(:,i1,j1)	! construct the total Hamiltonian from
                  dSij(:,i2,j2) = dSij1(:,i1,j1)	! construct the total Overlap Matrix from
!                   dHij(:,i2,j2) = dHij1(:,j1,i1)	! construct the total Hamiltonian from   ! WRONG
!                   dSij(:,i2,j2) = dSij1(:,j1,I1)	! construct the total Overlap Matrix from
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
end subroutine dHamil_tot_Press_DFTB ! CHECKED


subroutine dHamilton_one_Press_DFTB(i, atom_2, Scell, NSC, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij_press, dSij_press)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   integer, intent(in) :: i, atom_2
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, intent(in) :: norb, n_overlap  ! number of orbitals per atom (depends on the basis set), and number of overlap functions
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(out) :: dHij_press, dSij_press
   !---------------------------------------------------
   integer :: ki, kj
   real(8), dimension(9,norb,norb) ::  dH, dS ! hopping integrals
   dHij_press = 0.0d0
   if (atom_2 == 0) then ! Onsite contributions are constant, so derivatives are always zero:
      dHij_press = 0.0d0
      dSij_press = 0.0d0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_Press_DFTB(i, atom_2, Scell, NSC, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn,  dH, dS)    ! below
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1, norb
         do kj = 1, norb
!             dHij_press(:,kj,ki) = dH(:,ki,kj)   ! Hopping Integrals ! WRONG
!             dSij_press(:,kj,ki) = dS(:,ki,kj)   ! Hopping Integrals
            dHij_press(:,ki,kj) = dH(:,ki,kj)   ! Hopping Integrals ! CORRECT, tested on sp3 basis set
            dSij_press(:,ki,kj) = dS(:,ki,kj)   ! Hopping Integrals
         enddo ! kj
      enddo  ! ki
   endif
end subroutine dHamilton_one_Press_DFTB



subroutine dHopping_Press_DFTB(i, atom_2, Scell, NSC, norb, n_overlap, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dH, dS)
! subroutine making the derivatives of the hopping integrals
   integer, intent(in) :: i, atom_2	! atoms
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
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
   integer, pointer :: j

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

!          vec_M_dVij12 = M_dVij(i,j,:)*drij_dh
!          vec_M_dVij21 = M_dVij(j,i,:)*drij_dh
!          vec_M_dSVij12 = M_dSVij(i,j,:)*drij_dh
!          vec_M_dSVij21 = M_dSVij(j,i,:)*drij_dh
         vec_M_dVij21 = M_dVij(i,j,:)*drij_dh   ! correct
         vec_M_dVij12 = M_dVij(j,i,:)*drij_dh
         vec_M_dSVij21 = M_dSVij(i,j,:)*drij_dh
         vec_M_dSVij12 = M_dSVij(j,i,:)*drij_dh
         
         select case (norb)
         case (1) ! s
            call d_KS_s(vec_M_Vij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dH1)	! module "TB_Koster_Slater"
            call d_KS_s(vec_M_SVij12, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), dS1)	! module "TB_Koster_Slater"
         case (4) ! sp3
!             call d_KS_sp3_hetero_TEST(vec_M_Vij12, vec_M_Vij21, vec_M_dVij12, vec_M_dVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3),  dH1)	! module "TB_Koster_Slater" ! WRONG
!             call d_KS_sp3_hetero_TEST(vec_M_SVij12, vec_M_SVij21, vec_M_dSVij12, vec_M_dSVij21, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3),  dS1)	! module "TB_Koster_Slater"
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
   
   nullify(r, j)
end subroutine dHopping_Press_DFTB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complex Hamiltonian:
subroutine Complex_Hamil_DFTB(numpar, Scell, NSC, CHij, CSij, Ei, ksx, ksy, ksz, Err)
! This subroutine is unused for CDF calculations! A newer one is in the module "TB"
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(inout):: Ei	! [eV] energy levels for this complex Hamiltonian
   real(8), intent(in) :: ksx, ksy, ksz ! k-point to get Hamiltonian [relative]
   complex, DIMENSION(:,:), allocatable, INTENT(inout) :: CHij, CSij ! Complex Hamiltonian at k-point and overlap matrix
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------------------------------------------------------
   complex, dimension(:,:), allocatable :: CHij_temp, CHij_non, CSij_save
   integer :: Nsiz, j, nat, m, atom_2, i, j1, l, i1, k, norb
   real(8) :: kx, ky, kz, temp
   real(8), target :: nol
   real(8), pointer :: x1, y1, z1
   complex(8) :: expfac, SH_1
   character(200) :: Error_descript
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part


   Error_descript = ''
   nol = 0.0d0
!    temp = g_me*g_e/g_h*1d-10
   nat = size(Scell(NSC)%MDatoms) ! number of atoms
   norb =  identify_DFTB_orbitals_per_atom(numpar%basis_size_ind)  ! size of the basis set per atom, below
   Nsiz = size(Scell(NSC)%Ha,1)

   ! Allocate complex parameters:
   if (.not.allocated(Scell(NSC)%cPRRx)) allocate(Scell(NSC)%cPRRx(Nsiz,Nsiz))
   if (.not.allocated(Scell(NSC)%cPRRy)) allocate(Scell(NSC)%cPRRy(Nsiz,Nsiz))
   if (.not.allocated(Scell(NSC)%cPRRz)) allocate(Scell(NSC)%cPRRz(Nsiz,Nsiz))
   Scell(NSC)%cPRRx = dcmplx(0.0d0,0.0d0)	! to start with
   Scell(NSC)%cPRRy = dcmplx(0.0d0,0.0d0)	! to start with
   Scell(NSC)%cPRRz = dcmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij)) allocate(CHij(Nsiz,Nsiz))
   if (.not.allocated(CSij)) allocate(CSij(Nsiz,Nsiz))
   if (.not.allocated(CSij_save)) allocate(CSij_save(Nsiz,Nsiz))
   CHij = dcmplx(0.0d0,0.0d0)	! to start with
   CSij = dcmplx(0.0d0,0.0d0)	! to start with
   CSij_save = dcmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij_temp)) allocate(CHij_temp(Nsiz,Nsiz))
   if (.not.allocated(CHij_non)) allocate(CHij_non(Nsiz,Nsiz))
   CHij_non = dcmplx(0.0d0,0.0d0)	! to start with
   CHij_temp = dcmplx(0.0d0,0.0d0)	! to start with
   
   call Reciproc(Scell(NSC)%supce, Scell(NSC)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
   call Reciproc_rel_to_abs(ksx, ksy, ksz, Scell, NSC, kx, ky, kz) ! get absolute k-values, molue "Atomic_tools"

   ! 1) Construct complex Hamiltonian and overlap:
#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = nat
   ! Do the cycle (parallel) calculations:
   do j = Nstart, Nend, N_incr  ! each process does its own part
   !do j = 1,nat	! all atoms
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

               CHij_temp(k,l) = DCMPLX(Scell(NSC)%H_non(k,l),0.0d0)*expfac
               CSij(k,l) = DCMPLX(Scell(NSC)%Sij(k,l),0.0d0)*expfac

               if ((isnan(real(CHij_temp(k,l)))) .OR. isnan(aimag(CHij_temp(k,l)))) then
                  print*, i, j, k, l, CHij_temp(k,l)
                  pause 'CHij_temp ISNAN in Complex_Hamil_DFTB'
               endif
               if ((isnan(real(CSij(k,l)))) .OR. isnan(aimag(CSij(k,l)))) then
                  print*, i, j, k, l, CSij(k,l)
                  pause 'CSij ISNAN in Complex_Hamil_DFTB'
               endif

            enddo ! i1
         enddo ! j1
      enddo ! atom_2
   enddo ! j

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in Complex_Hamil_DFTB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'CHij_temp', CHij_temp) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'CSij', CSij) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
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
         
!          if (i > 0) then
            if ((abs(kx) < 1.0d-14) .AND. (abs(ky) < 1.0d-14) .AND. (abs(kz) < 1.0d-14)) then
               expfac = dcmplx(1.0d0,0.0d0)
            else
               expfac = exp(g_CI*dcmplx(kx*x1 + ky*y1 + kz*z1,0.0d0))
!                expfac = exp(g_CI*(kx*x1 + ky*y1 + kz*z1))
            endif
!             write(*,'(i5,i5,es,es,es,es,es)') j, i, x1, y1, z1, expfac

            do j1 = 1,norb ! all orbitals
               l = (j-1)*norb+j1
               do i1 = 1,norb ! all orbitals
                  k = (i-1)*norb+i1

                  CHij_temp(k,l) = DCMPLX(Scell(NSC)%H_non(k,l),0.0d0)*expfac
                  CSij(k,l) = DCMPLX(Scell(NSC)%Sij(k,l),0.0d0)*expfac
!                   CHij_temp(k,l) = DCMPLX(Scell(NSC)%H_non(l,k),0.0d0)*expfac   ! testing
!                   CSij(k,l) = DCMPLX(Scell(NSC)%Sij(l,k),0.0d0)*expfac
!                   if ((k == 188) .and. (l == 224)) then
!                      write(*,'(i5,i5,es,es,es,es,es,es,es,es,es)') k, l, CSij(k,l), CSij(l,k), expfac,  x1, y1, z1
!                   endif
!                   if ((l == 188) .and. (k == 224)) then
!                      write(*,'(i5,i5,es,es,es,es,es,es,es,es,es)') k, l, CSij(k,l), CSij(l,k), expfac,  x1, y1, z1
!                   endif
                  
                  if ((isnan(real(CHij_temp(k,l)))) .OR. isnan(aimag(CHij_temp(k,l)))) then
                     print*, i, j, k, l, CHij_temp(k,l)
                     pause 'CHij_temp ISNAN in Complex_Hamil_DFTB'
                  endif
                  if ((isnan(real(CSij(k,l)))) .OR. isnan(aimag(CSij(k,l)))) then
                     print*, i, j, k, l, CSij(k,l)
                     pause 'CSij ISNAN in Complex_Hamil_DFTB'
                  endif
                  
               enddo ! i1
            enddo ! j1
!          endif ! (i > 0)

      enddo ! atom_2
   enddo ! j
   !$omp end do 
   nullify(x1, y1, z1)
   !$omp end parallel
#endif

!    pause 'Construction ended'
   ! Temporarily save nonorthogonal Hamiltonian and overlap matrix:
   CHij_non = CHij_temp
   CSij_save = CSij
!    print*, 'Checking hermiticity of CH'
!    call check_hermiticity(CHij_non)   ! module "Algebra_tools"
!    print*, 'Checking hermiticity of CS'
!    call check_hermiticity(CSij_save)   ! module "Algebra_tools"
!    pause 'Both checked'

   !call sym_diagonalize(CHij, Ei, Err%Err_descript, numpar%MPI_param) ! modeule "Algebra_tools"
   ! Solve linear eigenproblem:
   ! 1) Orthogonalize the Hamiltonian using Loewdin procidure:
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization_c(numpar, Nsiz, CSij, CHij_temp, Err)	! module "TB_NRL"
   ! 2) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
   call sym_diagonalize(CHij_temp, Ei, Error_descript, numpar%MPI_param)
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Complex_Hamil_DFTB: '//trim(adjustl(Error_descript))
      if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         call Save_error_details(Err, 6, Error_descript)
      endif
      print*, trim(adjustl(Error_descript))
   endif
   ! 3) Convert the eigenvectors back into the non-orthogonal basis:
   call mkl_matrix_mult('N', 'N', CSij, CHij_temp, CHij)	! module "Algebra_tools"
   ! We need to renormalize the wave functions, as they are not normalized to 1 after this procidure:
!    do j = 1, Nsiz
!       CHij(:,j) = CHij(:,j) / DSQRT(dble(SUM( dconjg(CHij(:,j)) * CHij(:,j) )))
!    enddo

!     do j = 1, size(Ei)
!        write(*,'(i5,es,es,es,es,es)') j, Ei(j), CHij_temp(j,1), Scell(NSC)%Ei(j)
!     enddo
!     PAUSE 'Ei pause'
   
   ! 4) Calculate momentum operators:
   ! Optical matrix elements for non-orthogonal TB are taken from:
   ! arXiv:1805.08918v1 -- https://128.84.21.199/abs/1805.08918
#ifdef MPI_USED   ! use the MPI version
   ! Do the cycle (parallel) calculations:
   do j = Nstart, Nend, N_incr  ! each process does its own part
   !do j = 1,nat	! all atoms
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
                     SH_1 = DCMPLX(0.27d0,0.0d0) * (CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l))
                     Scell(NSC)%cPRRx(k,l) = SH_1
                     Scell(NSC)%cPRRy(k,l) = SH_1
                     Scell(NSC)%cPRRz(k,l) = SH_1
                  else	! different atoms at distance {x,y,z}:
                     SH_1 = CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l)
                     Scell(NSC)%cPRRx(k,l) = DCMPLX(x1,0.0d0)*SH_1
                     Scell(NSC)%cPRRy(k,l) = DCMPLX(y1,0.0d0)*SH_1
                     Scell(NSC)%cPRRz(k,l) = DCMPLX(z1,0.0d0)*SH_1
                  endif
                  if (real(Scell(NSC)%cPRRx(k,l)) .GT. 1d10) write(*,'(i5,i5,es,es, es,es, es,es, es, es)') i, j, Scell(NSC)%cPRRx(k,l),  CHij_non(k,l), CSij_save(k,l),  Ei(k), x1
                  if (real(Scell(NSC)%cPRRy(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRy(k,l)
                  if (real(Scell(NSC)%cPRRz(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRz(k,l)
               enddo ! i1
            enddo ! j1
         endif ! (i > 0)
      enddo ! atom_2
   enddo ! j

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in Complex_Hamil_DFTB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Scell(NSC)%cPRRx', Scell(NSC)%cPRRx) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Scell(NSC)%cPRRy', Scell(NSC)%cPRRy) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Scell(NSC)%cPRRz', Scell(NSC)%cPRRz) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
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
                     SH_1 = DCMPLX(0.27d0,0.0d0) * (CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l))
                     Scell(NSC)%cPRRx(k,l) = SH_1
                     Scell(NSC)%cPRRy(k,l) = SH_1
                     Scell(NSC)%cPRRz(k,l) = SH_1
                  else	! different atoms at distance {x,y,z}:
                     SH_1 = CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l)
                     Scell(NSC)%cPRRx(k,l) = DCMPLX(x1,0.0d0)*SH_1
                     Scell(NSC)%cPRRy(k,l) = DCMPLX(y1,0.0d0)*SH_1
                     Scell(NSC)%cPRRz(k,l) = DCMPLX(z1,0.0d0)*SH_1
                  endif
                  if (real(Scell(NSC)%cPRRx(k,l)) .GT. 1d10) write(*,'(i5,i5,es,es, es,es, es,es, es, es)') i, j, Scell(NSC)%cPRRx(k,l),  CHij_non(k,l), CSij_save(k,l),  Ei(k), x1
                  if (real(Scell(NSC)%cPRRy(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRy(k,l)
                  if (real(Scell(NSC)%cPRRz(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRz(k,l)
               enddo ! i1
            enddo ! j1
         endif ! (i > 0)
      enddo ! atom_2
   enddo ! j
   !$omp end do
   nullify(x1, y1, z1)
   !$omp end parallel
#endif

   ! Convert to SI units used later:
   !temp = g_me*g_e/g_h*1d-10 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
   ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
   !temp = 1.0d0
   temp = 1.0d0 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
   Scell(NSC)%cPRRx = Scell(NSC)%cPRRx * temp
   Scell(NSC)%cPRRy = Scell(NSC)%cPRRy * temp
   Scell(NSC)%cPRRz = Scell(NSC)%cPRRz * temp

   deallocate(CHij_non, CHij_temp, CSij_save)
   nullify(x1, y1, z1)
end subroutine Complex_Hamil_DFTB



!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
! Repulsive part:


subroutine get_Erep_s_DFTB_no(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_DFTB_no), dimension(:,:), intent(in) :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a    ! [eV] total repulsive energy
   !=====================================================

   a = 0.0d0   ! no repulsive
end subroutine get_Erep_s_DFTB_no



subroutine get_Erep_s_DFTB(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_DFTB), dimension(:,:), intent(in) :: TB_Repuls
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a    ! [eV] total repulsive energy
   !=====================================================
   real(8) :: E_pot, E_rep_one, E_pot_array(Scell(NSC)%Na)
   integer :: i1, m, atom_2, j1
   integer, pointer :: KOA1, KOA2
   real(8), pointer :: r
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part
   
   a = 0.0d0

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Scell(NSC)%Na
   E_pot_array = 0.0d0
   ! Do the cycle (parallel) calculations:
   do i1 = Nstart, Nend, N_incr  ! each process does its own part
   !do i1 = 1, Scell(NSC)%Na
      !E_pot = 0.0d0
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if (j1 /= i1) then
            KOA1 => Scell(NSC)%MDatoms(i1)%KOA
            KOA2 => Scell(NSC)%MDatoms(j1)%KOA
            r => Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            E_rep_one = DFTB_repulsive_one(TB_Repuls(KOA1,KOA2), r)    ! below
            a = a + E_rep_one
            E_pot_array(i1) = E_pot_array(i1) + E_rep_one
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_Erep_s_DFTB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'a', a) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'E_pot_array', E_pot_array) ! module "MPI_subroutines"

   do i1 = 1, Scell(NSC)%Na
      ! And save for each atom:
      Scell(NSC)%MDAtoms(i1)%Epot = E_pot_array(i1)*0.5d0 ! to exclude double-counting
   enddo

#else    ! OpenMP to use instead
!$omp parallel private(i1, m, atom_2, j1, KOA1, KOA2, r, E_rep_one, E_pot)
!$omp do reduction( + : a)
   do i1 = 1, Scell(NSC)%Na
      E_pot = 0.0d0
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if (j1 /= i1) then
            KOA1 => Scell(NSC)%MDatoms(i1)%KOA
            KOA2 => Scell(NSC)%MDatoms(j1)%KOA
            r => Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            E_rep_one = DFTB_repulsive_one(TB_Repuls(KOA1,KOA2), r)    ! below
            a = a + E_rep_one
            E_pot = E_pot + E_rep_one
         endif ! (j1 .NE. i1)
      enddo ! j1
      ! And save for each atom:
      Scell(NSC)%MDAtoms(i1)%Epot = E_pot*0.5d0 ! to exclude double-counting
      !print*, i1, E_pot, Scell(NSC)%MDAtoms(i1)%Epot
   enddo ! i1
!$omp end do
   nullify(KOA1, KOA2, r)
!$omp end parallel
#endif
   a = a/2.0d0 ! it was doubled
   nullify(KOA1, KOA2, r)
end subroutine get_Erep_s_DFTB


 function DFTB_repulsive_one(TB_Repuls, r) result (F)
   real(8) F    ! Eqs.(1) vs (6) from [1]
   type(TB_Rep_DFTB), intent(in)   :: TB_Repuls
   real(8), intent(in) :: r ! [A] distance between the atoms
   select case (TB_Repuls%ToP) ! which repulsive potential is used: 0=polinomial, 1=spline
   case (0)
      F = DFTB_polinomial(TB_Repuls%c, TB_Repuls%rcut, r)   ! below
   case default ! spline
      F = DFTB_spline(TB_Repuls%a, TB_Repuls%V_rep, TB_Repuls%R, TB_Repuls%rcut_spline, r)  ! below
   end select
end function DFTB_repulsive_one


pure function DFTB_polinomial(C, rcut, r) result (F)
   real(8) F
   real(8), dimension(:), intent(in) :: C   ! polinomial coefficients
   real(8), intent(in) :: rcut   ! cutoff radius [A]
   real(8), intent(in) :: r ! [A] distance between the atoms
   !------------------------------------
   real(8) :: rcutr, dblei
   integer :: i
   F = 0.0d0 ! to start with
   if (r <= rcut) then ! if too far, no interaction
      do i = 1,8
         dblei = dble(i + 1)
         rcutr = rcut - r
         F = F + C(i)*rcutr**dblei
      enddo
   endif
end function DFTB_polinomial


 function DFTB_spline(a, c, R, rcut, r_dist) result (F)
   real(8) F    ! [eV]
   real(8), dimension(:), intent(in) :: a   ! exp cpefficients at short distances
   real(8), dimension(:,:), intent(in) :: c   ! spline coefficients
   real(8), dimension(:), intent(in) :: R   ! [A]
   real(8), intent(in) :: rcut  ! [A]
   real(8), intent(in) :: r_dist ! [A] distance between the atoms
   !----------------------------------------
   integer :: i, Nsiz, i_array
   real(8) :: rr0, rr02, rr03, dr
   Nsiz = size(R)
   
   !print*, 'DFTB_spline:', r_dist, rcut
   
   if (r_dist >= rcut) then
      F = 0.0d0 ! too far, no interaction
   elseif (r_dist > R(Nsiz)) then    ! smooth approach to zero
      rr0 = r_dist - R(Nsiz)
      rr02 = rr0 * rr0
      rr03 = rr02 * rr0
      F = c(Nsiz,1) + c(Nsiz,2)*rr0 + c(Nsiz,3)*rr02 + c(Nsiz,4)*rr03 + c(Nsiz,5)*rr02*rr02 + c(Nsiz,6)*rr03*rr02
   elseif (r_dist <= R(1)) then  ! "exponential wall"
      F = exp(-a(1)*r_dist + a(2))*g_au2ev + a(3)    ! [eV]
   else ! spline
      ! This works only for equidistant arrays:
      !i_array = CEILING((r_dist-R(1))/(R(2)-R(1)))   ! index of the nearest radial grid point for the equally-spaced grid
      ! This works for any monotonous array:
      call Find_in_array_monoton(R, r_dist, i_array)  ! module "Little_subroutines"
      i_array = i_array - 1   ! we need floor instead of ceiling
      !print*, 'DFTB_spline:', i_array, r_dist, R(i_array), R(i_array+1)
      rr0 = r_dist - R(i_array)
      rr02 = rr0 * rr0
      rr03 = rr02 * rr0
      F = c(i_array,1) + c(i_array,2)*rr0 + c(i_array,3)*rr02 + c(i_array,4)*rr03
      if (i_array > Nsiz) then
         write(6,'(a,f,f)') 'Error 1 in DFTB_spline:', i_array, Nsiz
      elseif ((r_dist < R(i_array)) .or. (R(i_array+1) < r_dist) ) then
         write(6,'(a,f,f,f)') 'Error 2 in DFTB_spline:', R(i_array), R(i_array+1), r_dist
         write(6,'(i,f,f,f)') i_array, R(i_array), R(1), R(Nsiz)
      endif
   endif
end function DFTB_spline


!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
! Forces for the repulsive part:

subroutine dErdr_s_DFTB_no(TB_Repuls, Scell, NSC) ! derivatives of the repulsive energy by s
   type(TB_Rep_DFTB_no), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !---------------------------------------
   integer :: ian, n
   n = Scell(NSC)%Na ! number of atoms
   !$omp PARALLEL private(ian)
   !$omp DO
   do ian = 1, n  ! Forces for all atoms
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
   enddo ! ian
   !$omp end do
   !$omp end parallel
END subroutine dErdr_s_DFTB_no


subroutine dErdr_s_DFTB(TB_Repuls, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   type(TB_Rep_DFTB), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(inout) :: numpar   ! all numerical parameters
   !---------------------------------------
   real(8), dimension(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, ian, dik, djk, n, atom_2
   real(8), dimension(:,:), allocatable :: Erx_s
   integer, pointer :: KOA1, KOA2, m, j1
   real(8), pointer ::  x, y, z
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   n = Scell(NSC)%Na ! number of atoms
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = n
   ! Do the cycle (parallel) calculations:
   do ian = Nstart, Nend, N_incr  ! each process does its own part
   !do ian = 1, n	! Forces for all atoms
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
      do i1 = 1, n	! contribution from all atoms
         if (ian == i1) then	! Kroniker delta
            dik = 1
         else
            dik = 0
         endif
         dpsi = 0.0d0
         m => Scell(NSC)%Near_neighbor_size(i1)
         KOA1 => Scell(NSC)%MDatoms(i1)%KOA
         do atom_2 = 1,m		! do only for atoms close to that one
            j1 => Scell(NSC)%Near_neighbor_list(i1, atom_2)	! this is the list of such close atoms
            if (j1 > 0) then
               if (ian == j1) then	! Kroniker delta
                  djk = 1
               else
                  djk = 0
               endif
               cos_if:if ((dik-djk) /= 0) then ! without it, it gives ERROR
                  KOA2 => Scell(NSC)%MDatoms(j1)%KOA
                  x => Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X, Y, Z
                  y => Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z => Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z

                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3) ! correct
                  x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                  x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

                  a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R
                  b =d_DFTB_repulsive_one(TB_Repuls(KOA1, KOA2), a_r) ! below

                  ddlta = dble(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(:) = dpsi(:) + b_delta*x1(:)
               endif cos_if
            endif ! j1 > 0
         enddo ! j1

         Erx_s(:,ian) = Erx_s(:,ian) + dpsi(:) ! potential part in X-coordinate
      enddo ! i1
   enddo ! ian

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in dErdr_s_DFTB:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Erx_s', Erx_s) ! module "MPI_subroutines"

   do ian = 1, n	! Forces for all atoms
      ! Add exponential wall force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)*0.5d0	! factor 0.5 to compensate for double-counting
   enddo

#else    ! OpenMP to use instead
   !$omp PARALLEL private(ian, i1, dik, dpsi, m, KOA1, atom_2, j1, djk, KOA2, x,y,z, x1, b, a_r,ddlta,b_delta)
   !$omp DO
   do ian = 1, n	! Forces for all atoms
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
      do i1 = 1, n	! contribution from all atoms
         if (ian == i1) then	! Kroniker delta
            dik = 1
         else
            dik = 0
         endif
         dpsi = 0.0d0
         m => Scell(NSC)%Near_neighbor_size(i1)
         KOA1 => Scell(NSC)%MDatoms(i1)%KOA
         do atom_2 = 1,m		! do only for atoms close to that one  
            j1 => Scell(NSC)%Near_neighbor_list(i1, atom_2)	! this is the list of such close atoms
            if (j1 > 0) then
               if (ian == j1) then	! Kroniker delta
                  djk = 1
               else
                  djk = 0
               endif
               cos_if:if ((dik-djk) /= 0) then ! without it, it gives ERROR
                  KOA2 => Scell(NSC)%MDatoms(j1)%KOA
                  x => Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X, Y, Z
                  y => Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z => Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z
                  
                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3) ! correct
                  x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                  x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

                  !x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)   ! incorrect
                  !x1(2) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
                  !x1(3) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)
                  
                  a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R
                  b =d_DFTB_repulsive_one(TB_Repuls(KOA1, KOA2), a_r) ! below

                  ddlta = dble(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(:) = dpsi(:) + b_delta*x1(:)
               endif cos_if
            endif ! j1 > 0
         enddo ! j1
         
         Erx_s(:,ian) = Erx_s(:,ian) + dpsi(:) ! potential part in X-coordinate
      enddo ! i1
      ! Add exponential wall force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)*0.5d0	! factor 0.5 to compensate for double-counting
   enddo ! ian
   !$omp end do
   nullify(j1, m, KOA1, KOA2, x, y, z)
   !$omp end parallel
#endif

   deallocate(Erx_s)
   nullify(j1, m, KOA1, KOA2, x, y, z)
END subroutine dErdr_s_DFTB



function d_DFTB_repulsive_one(TB_Repuls, r) result (F)
   real(8) F    ! Eqs.(1) vs (6) from [1]
   type(TB_Rep_DFTB), intent(in)   :: TB_Repuls
   real(8), intent(in) :: r ! [A] distance between the atoms
   select case (TB_Repuls%ToP) ! which repulsive potential is used: 0=polinomial, 1=spline
   case (0)
      F = d_DFTB_polinomial(TB_Repuls%c, TB_Repuls%rcut, r)   ! below
   case default ! spline
      F = d_DFTB_spline(TB_Repuls%a, TB_Repuls%V_rep, TB_Repuls%R, TB_Repuls%rcut_spline, r)  ! below
   end select
end function d_DFTB_repulsive_one


pure function d_DFTB_polinomial(C, rcut, r) result (F)
   real(8) F
   real(8), dimension(:), intent(in) :: C   ! polinomial coefficients
   real(8), intent(in) :: rcut   ! cutoff radius [A]
   real(8), intent(in) :: r ! [A] distance between the atoms
   !------------------------------------
   real(8) :: rcutr, dblei
   integer :: i
   F = 0.0d0 ! to start with
   if (r <= rcut) then ! if too far, no interaction
      do i = 1,8
         dblei = dble(i)
         rcutr = rcut - r
         F = F - (dblei+1.0d0)*C(i)*rcutr**dblei
      enddo
   endif
end function d_DFTB_polinomial





function d_DFTB_spline(a, c, R, rcut, r_dist) result (F)
   real(8) F    ! [eV]
   real(8), dimension(:), intent(in) :: a   ! exp cpefficients at short distances
   real(8), dimension(:,:), intent(in) :: c   ! spline coefficients
   real(8), dimension(:), intent(in) :: R   ! [A]
   real(8), intent(in) :: rcut  ! [A]
   real(8), intent(in) :: r_dist ! [A] distance between the atoms
   !----------------------------------------
   integer :: i, Nsiz, i_array
   real(8) :: rr0, rr02, rr03, dr
   Nsiz = size(R)
   
   if (r_dist >= rcut) then
      F = 0.0d0 ! too far, no interaction
   elseif (r_dist > R(Nsiz)) then    ! smooth approach to zero
      rr0 = r_dist - R(Nsiz)
      rr02 = rr0 * rr0
      rr03 = rr02 * rr0
      F = c(Nsiz,2) + 2.0d0*c(Nsiz,3)*rr0 + 3.0d0*c(Nsiz,4)*rr02 + 4.0d0*c(Nsiz,5)*rr03 + 5.0d0*c(Nsiz,6)*rr02*rr02
   elseif (r_dist <= R(1)) then  ! "exponential wall"
      F = -a(1)*exp(-a(1)*r_dist + a(2))    ! [eV/A]
   else ! spline
      ! This only works for equidistant arrays:
      !i_array = CEILING((r_dist-R(1))/(R(2)-R(1)))   ! index of the nearest radial grid point for the equally-spaced grid
      ! This works for any monotonous array:
      call Find_in_array_monoton(R, r_dist, i_array)  ! module "Little_subroutines"
      i_array = i_array - 1   ! we need floor instead of ceiling

      if (i_array > Nsiz) then
         write(6,'(a,f,f)') 'Error 1 in d_DFTB_spline:', i_array, Nsiz
      elseif ((r_dist < R(i_array)) .or. (R(i_array+1) < r_dist)) then
         write(6,'(a,f,f,f)') 'Error 2 in d_DFTB_spline:', R(i_array), R(i_array+1), r_dist
         write(6,'(i,f,f,f)') i_array, R(i_array), R(1), R(Nsiz)
      endif

      rr0 = r_dist - R(i_array)
      rr02 = rr0 * rr0
      F = c(i_array,2) + 2.0d0*c(i_array,3)*rr0 + 3.0d0*c(i_array,4)*rr02
   endif
end function d_DFTB_spline



subroutine dErdr_Pressure_s_DFTB(TB_Repuls, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_DFTB), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   real(8), dimension(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i, k, l, n, atom_2
   integer, pointer :: KOA1, KOA2, m, j
   real(8) r, rcur(3), scur(3), PForce(3,3)
   real(8) df_psy, psi, dpsy

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = Scell(NSC)%Na ! number of atoms
      Scell(NSC)%SCforce%rep = 0.0d0    ! to start with

      PForce = 0.0d0 ! to start with
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         m => Scell(NSC)%Near_neighbor_size(i)
         KOA1 => Scell(NSC)%MDatoms(i)%KOA
         do atom_2 = 1,m		! do only for atoms close to that one  
            j => Scell(NSC)%Near_neighbor_list(i, atom_2)	! this is the list of such close atoms
            !if (j > 0) then
            if ((i /=  j) .AND. (j  >  0)) then
               KOA2 => Scell(NSC)%MDatoms(j)%KOA
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
               dpsy = d_DFTB_repulsive_one(TB_Repuls(KOA1, KOA2), r) ! function above

               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! (j > 0)
         enddo ! atom_2

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)*0.5d0	! factor 0.5 to compensate for double-counting
            enddo ! l
         enddo ! k
      enddo ! i
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! repulsive part
   endif
   nullify(KOA1, KOA2, m, j)
end subroutine dErdr_Pressure_s_DFTB


subroutine dErdr_Pressure_s_DFTB_no(TB_Repuls, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_DFTB_no), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      Scell(NSC)%SCforce%rep = 0.0d0    ! to start with
   endif
end subroutine dErdr_Pressure_s_DFTB_no




!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! Utilities:

pure function identify_DFTB_basis_size(ind) result (N)
   integer, intent(in) :: ind   ! index of the basis set type
   integer :: N ! number of radial overlap functions
   select case (ind)
   case (0) ! s
      N = 1
   case (1) ! sp3
      N = 4
   case default ! sp3d5
      N = 10
   end select
end function identify_DFTB_basis_size

pure function identify_DFTB_orbitals_per_atom(ind) result(N)
  integer, intent(in) :: ind   ! index of the basis set type
   integer :: N ! number of orbitals per atom
   select case (ind)
   case (0) ! s
      N = 1
   case (1) ! sp3
      N = 4
   case default ! sp3d5
      N = 9
   end select
end function identify_DFTB_orbitals_per_atom

END MODULE TB_DFTB
