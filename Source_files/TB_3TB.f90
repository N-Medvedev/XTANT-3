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
use Atomic_tools, only : Reciproc_rel_to_abs
use Little_subroutines, only : linear_interpolation, Fermi_function, d_Fermi_function, Find_in_array_monoton
use Electron_tools, only : find_band_gap
use TB_Koster_Slater
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, test_orthogonalization_c, Loewdin_Orthogonalization, Loewdin_Orthogonalization_c
use TB_DFTB, only : identify_DFTB_basis_size


implicit none


real(8), parameter :: m_a = 2.0d0   ! exponential decay parameter according to [1]



 contains



subroutine Construct_Vij_3TB(numpar, TB, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)
   type(Numerics_param), intent(in), target :: numpar 	! all numerical parameters
   type(TB_H_3TB), dimension(:,:), intent(in), target :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable :: M_Vij	! matrix of hoppings for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dVij	! matrix of derivatives of hoppings for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_SVij	! matrix of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dSVij  ! matrix of derivatives of Overlap for S-matrix for all pairs of atoms, all orbitals
   !----------------------------------
   real(8) :: x, y, z, r, sx, sy, sz, r1, fcut, d_fcut, Fermi, dFermi, Laguer(6), d_Laguer(6)
   integer :: i, j, atom_2, ki, N_bs, ihop
   real(8), pointer :: rm
   integer, pointer :: nat, m, KOA1, KOA2

   nat => Scell(NSC)%Na	! number of atoms in the supercell
   ! number of hopping integrals for this basis set in 3TB:
   N_bs = identify_DFTB_basis_size(numpar%N_basis_size)  ! module "TB_DFTB"

   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,N_bs))	    ! each pair of atoms, all  V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,N_bs))	! each pair of atoms, all  dV functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,N_bs))	! each pair of atoms, all S functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,N_bs))	! each pair of atoms, all dS functions

   !$OMP WORKSHARE
   M_Vij = 0.0d0
   M_dVij = 0.0d0
   M_SVij = 0.0d0
   M_dSVij = 0.0d0
   !$OMP END WORKSHARE

   ! Construct matrix of all the radial functions for each pair of atoms:
!$omp PARALLEL
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, Laguer, d_Laguer, ihop, Fermi, dFermi)
   AT1:do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms

         KOA1 => Scell(NSC)%MDatoms(j)%KOA
         KOA2 => Scell(NSC)%MDatoms(i)%KOA

         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! Get Laguerre polynomials and its derivatives to be reused below:
         call get_Laguerres(r, Laguer) ! below
         call get_d_Laguerres(r, d_Laguer) ! below

         ! 2-body interaction terms:
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
         case default    ! sp3d5
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
         case default    ! sp3d5
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


         ! 3-body interaction terms:




!          ! Testing:
!          M_dVij = 0.0d0
!          M_dSVij = 0.0d0
!          M_Vij = 0.0d0
!          M_SVij = 0.d00

      enddo AT2
   enddo AT1
!$omp end do
!$omp END PARALLEL

   nullify(rm, nat, m, KOA1, KOA2)	! clean up at the end
end subroutine Construct_Vij_3TB


subroutine get_Laguerres(r_given, Laguer)
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(6), intent(out) :: Laguer ! polynomials up to 6
   real(8) :: d, ad
   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au   ! [A] -> [Bohr]
   ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]

   call Laguerre_up_to_6(ad, Laguer)   ! module "Algebra_tools"
end subroutine get_Laguerres


subroutine get_d_Laguerres(r_given, Laguer)
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(6), intent(out) :: Laguer ! derivatives of the polynomials up to 6
   real(8) :: d, ad
   ! normalized distance, since coefficients are fitted in this units [1]:
   d = r_given * g_A2au   ! [A] -> [Bohr]
   ad = m_a * d ! for some reason, it is multiplied with the scaling parameter [2]

   call d_Laguerre_up_to_6(ad, Laguer)   ! module "Algebra_tools"
end subroutine get_d_Laguerres



subroutine radial_function_3TB(r_given, fx_coefs, Laguer, sh_ind, num_lag) result(f_out)
   real(8) :: f_out    ! value found in the array and interpolated
   real(8), intent(in) :: r_given  ! [A] given radial coordinate
   real(8), dimension(:,:), intent(in) :: fx_coefs    ! coefficients before Leguerres in 2-body Hamiltonian
   real(8), dimension(:), intent(in) :: Laguer   ! precalculated Laguerre polinomials
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
   do i = 1, num_lag ! for all polynomials required
      f_out = f_out + fx_coefs(sh_ind, i) * Laguer(i)
   enddo
   f_out = exp_ad * f_out

end subroutine radial_function_3TB



subroutine d_radial_function_3TB(r_given, M_Vij, fx_coefs, d_Laguer, sh_ind, num_lag) result(f_out)
   real(8) :: f_out    ! value found in the array and interpolated
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
end subroutine d_radial_function_3TB



! Tight Binding Hamiltonian within 3TB parametrization:
subroutine construct_TB_H_3TB(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! materil parameters
   type(TB_H_3TB), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
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
   call Hamil_tot_3TB(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err)	! see below

!    ! Test (comment out for release):
!    call test_nonorthogonal_solution(Scell(NSC)) ! module "TB_NRL"

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_3TB






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


subroutine dErdr_s_3TB(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   type(TB_Rep_3TB), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer ian, n
   n = size(atoms)
   !$omp PARALLEL DO private(ian)
   do ian = 1, n  ! Forces for all atoms
     Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! no added repulsion
   enddo ! ian
   !$OMP END PARALLEL DO
END subroutine dErdr_s_3TB


subroutine dErdr_Pressure_s_3TB(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_3TB), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms ! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      Scell(NSC)%SCforce%rep = 0.0d0   ! no added repulsion
   endif
end subroutine dErdr_Pressure_s_3TB


END MODULE TB_3TB
