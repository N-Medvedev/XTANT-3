! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2021 Nikita Medvedev
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
! This module contains subroutines to deal with TB hamiltonian in the BOP parametrization
! This parameterization is described in:
! [1] https://arxiv.org/abs/1909.04561

MODULE TB_BOP

use Universal_constants
use TB_Koster_Slater
use Objects
use Electron_tools, only : find_band_gap
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, test_orthogonalization_c, Loewdin_Orthogonalization, Loewdin_Orthogonalization_c
use TB_DFTB, only : identify_DFTB_basis_size, identify_DFTB_orbitals_per_atom, Hopping_DFTB, Get_overlap_S_matrix_DFTB
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity
use Atomic_tools, only : Reciproc_rel_to_abs
use Little_subroutines, only : linear_interpolation, Find_in_array_monoton

implicit none

 
 contains


!hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
subroutine construct_TB_H_BOP(numpar, TB_Hamil, matter, M_Vij, M_SVij, M_E0ij, M_lmn, Scell, NSC, Err)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(TB_H_BOP), dimension(:,:), intent(in) :: TB_Hamil	! All parameters of the Hamiltonian of TB
   type(solid), intent(inout) :: matter	! materil parameters
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij, M_E0ij   ! matrix of Overlap functions for all pairs of atoms, all orbitals
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

    ! Construct TB Hamiltonian (with BOP parameters),  orthogonalize it,  and then solve generalized eigenvalue problem:
   call Hamil_tot_BOP(numpar, Scell, NSC, M_Vij, M_SVij, M_E0ij, M_lmn, Err)	! see below

!    ! Test (comment out for release):
!    call test_nonorthogonal_solution(Scell(NSC)) ! module "TB_NRL"

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_BOP


subroutine Hamil_tot_BOP(numpar, Scell, NSC, M_Vij, M_SVij, M_E0ij, M_lmn, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
!    type(TB_H_BOP), dimension(:,:), intent(in) :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij, M_E0ij   ! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------------------
   real(8), dimension(:,:), allocatable :: Hij  ! Hamiltonian
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
   n_orb =  identify_DFTB_orbitals_per_atom(numpar%N_basis_size)  ! size of the basis set per atom, below
   nat = Scell(NSC)%Na  ! number of atoms in the supercell
   Nsiz = size(Scell(NSC)%Ha,1) ! size of the total basis set

   if (.not.allocated(Sij)) allocate(Sij(Nsiz,Nsiz))
   if (.not.allocated(Hij)) allocate(Hij(Nsiz,Nsiz))
   Sij = 0.0d0
   Hij = 0.0d0

   !-----------------------------------
   ! 1) Construct non-orthogonal Hamiltonian H and Overlap matrix S in 2 steps:
!$omp parallel private(j, m, atom_2, i, KOA1, KOA2, j1, l, i1, k, Hij1, Sij1)
if (.not.allocated(Hij1)) allocate(Hij1(n_orb,n_orb), source = 0.0d0)
if (.not.allocated(Sij1)) allocate(Sij1(n_orb,n_orb), source = 0.0d0)
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
            call Hamilton_one_BOP(numpar%N_basis_size, Scell, NSC, M_Vij(j,i,:), M_E0ij, M_lmn, j, i, Hij1)   ! below
            
            ! Construct overlap matrix for this pair of atoms:
            call Get_overlap_S_matrix_BOP(numpar%N_basis_size, M_SVij(j,i,:), M_lmn(:,j,i), j, i, Sij1)  ! below

            do j1 = 1,n_orb ! all orbitals
               l = (j-1)*n_orb+j1
               do i1 = 1,n_orb ! all orbitals
                  k = (i-1)*n_orb+i1
                  Hij(k,l) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  Sij(k,l) = Sij1(i1,j1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
!                   Hij(k,l) = Hij1(j1,i1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
!                   Sij(k,l) = Sij1(j1,i1) ! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
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
         if (i < j) then ! it's a new pair of atoms, calculate everything
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
    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, check_M=.true.)    ! testing
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, use_DSYEV=.false.)
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript) ! module "Algebra_tools"
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Hamil_tot_BOP: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
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
!                      if (abs(i) > nat*n_orb) print*, 'TEST i', i, j, l, k
                     if (i == j) then ! for the diagonal elements according to Trani:
!                         SH_1 = 0.270d0*(Scell(NSC)%Ei(k) - Scell(NSC)%Ei(l))
                        SH_1 =  0.0d0 
!                         SH_1 =  0.270d0*(Scell(NSC)%H_non(k,l) - Scell(NSC)%Ei(k) * Scell(NSC)%Sij(k,l))
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
      !$omp end parallel
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0      
      Scell(NSC)%PRRx = Scell(NSC)%PRRx * temp
      Scell(NSC)%PRRy = Scell(NSC)%PRRy * temp
      Scell(NSC)%PRRz = Scell(NSC)%PRRz * temp
   endif ! (numpar%optic_model .EQ. 3) 
   
   nullify(KOA1, KOA2, m, x, y, z)
   deallocate(Hij, Sij)

end subroutine Hamil_tot_BOP


subroutine Get_overlap_S_matrix_BOP(basis_ind, M_SVij, M_lmn, j, i, Sij)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   real(8), dimension(:), intent(in) :: M_SVij    ! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:), intent(in) :: M_lmn   ! matrix of directional cosines
   integer, intent(in) :: j, i  ! atoms indices
   real(8), dimension(:,:), intent(out) :: Sij  ! overlap matrix, all orbitals in sp3d5 basis set
   !---------------------------------------   
   integer :: k
   if (i == j) then ! Onsite contributions
      Sij = 0.0d0
      forall (k=1:size(Sij,1)) Sij(k,k)=1.0d0
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:
      call Hopping_BOP(basis_ind, M_SVij, M_lmn, Sij)    ! below
   endif
end subroutine Get_overlap_S_matrix_BOP



subroutine Hamilton_one_BOP(basis_ind, Scell, NSC, M_Vij, M_E0ij, M_lmn, j, i, Hij) 
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: M_Vij    ! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_E0ij  ! ons-cite elements
   real(8), dimension(:,:,:), intent(in) :: M_lmn   ! matrix of directional cosines
   integer, intent(in) :: j, i  ! atoms indices
   real(8), dimension(:,:), intent(out) :: Hij  ! hamiltonian, all orbitals in sp3d5 basis set
   !---------------------------------------
   if (i == j) then ! Onsite contributions
      call On_site_BOP(basis_ind, Scell, NSC, M_E0ij, M_lmn, j, Hij)     ! below
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:
      call Hopping_BOP(basis_ind, M_Vij(:), M_lmn(:,j,i), Hij)    ! below
   endif
end subroutine Hamilton_one_BOP


subroutine On_site_BOP(basis_ind, Scell, NSC, M_E0ij, M_lmn, j, ts)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), intent(in) :: M_E0ij    ! Onsite functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn   ! cosine directions
   integer, intent(in) :: j     ! index
   real(8), dimension(:,:), intent(out) :: ts   ! overlap integerals [eV]
   !=============================================
   integer :: atom_2
   integer, pointer :: m, k
   ts = 0.0d0   ! to start with
   
   m => Scell(NSC)%Near_neighbor_size(j)
   do atom_2 = 1, m  ! all neighbours of the atom "j"
      k => Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
      ts(1,1) = ts(1,1) + t_s_s(M_E0ij(j,k, 1)) ! s,  module "TB_Koster_Slater"
      if (basis_ind > 0) then  ! p
         ts(2,2) = ts(2,2) + t_pa_pa(M_lmn(1, j,k), M_E0ij(j,k,2), M_E0ij(j,k,3))  ! px,  module "TB_Koster_Slater"
         ts(3,3) = ts(3,3) + t_pa_pa(M_lmn(2, j,k), M_E0ij(j,k,2), M_E0ij(j,k,3))  ! py,  module "TB_Koster_Slater"
         ts(4,4) = ts(4,4) + t_pa_pa(M_lmn(3, j,k), M_E0ij(j,k,2), M_E0ij(j,k,3))  ! pz,  module "TB_Koster_Slater"
      endif
      if (basis_ind > 1) then  ! d
         ts(5,5) = ts(5,5) + t_dab_dab(M_lmn(1, j,k), M_lmn(2, j,k), M_lmn(3, j,k), M_E0ij(j,k,4), M_E0ij(j,k,5), M_E0ij(j,k,6))  ! dxy,  module "TB_Koster_Slater"
         ts(6,6) = ts(6,6) + t_dab_dab(M_lmn(1, j,k), M_lmn(3, j,k), M_lmn(2, j,k), M_E0ij(j,k,4), M_E0ij(j,k,5), M_E0ij(j,k,6))  ! dxz,  module "TB_Koster_Slater"
         ts(7,7) = ts(7,7) + t_dab_dab(M_lmn(2, j,k), M_lmn(3, j,k), M_lmn(1, j,k), M_E0ij(j,k,4), M_E0ij(j,k,5), M_E0ij(j,k,6))  ! dyz,  module "TB_Koster_Slater"
         ts(8,8) = ts(8,8) + t_dx2_y2_dx2_y2(M_lmn(1, j,k), M_lmn(2, j,k), M_lmn(3, j,k), M_E0ij(j,k,4), M_E0ij(j,k,5), M_E0ij(j,k,6))    ! dx2-y2,  module "TB_Koster_Slater"
         ts(9,9) = ts(9,9) + t_d3z2_r2_d3z2_r2(M_lmn(1, j,k), M_lmn(2, j,k), M_lmn(3, j,k), M_E0ij(j,k,4), M_E0ij(j,k,5), M_E0ij(j,k,6))    ! d3z2-r2,  module "TB_Koster_Slater"
      endif
   enddo
   
   nullify(m, k)
end subroutine On_site_BOP


subroutine Hopping_BOP(basis_ind, M_Vij, M_lmn, ts)
   integer, intent(in) :: basis_ind ! index of the basis set used: 0=s, 1=sp3, 2=sp3d5
   real(8), dimension(:), intent(in) :: M_Vij    ! Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:), intent(in) :: M_lmn   ! cosine directions
   real(8), dimension(:,:), intent(out) :: ts   ! overlap integerals [eV]
   !=============================================
    real(8), dimension(:), allocatable :: vec_M_SVij12, vec_M_SVij21
   ! Construct the overlap integrals including angular part 
   select case (basis_ind)
   case (0) ! for s basis set:
      allocate(vec_M_SVij12(1))
      call BOP_matrices_to_vec(basis_ind, M_Vij, vec_M_SVij12, vec_M_SVij21)   ! below
      call KS_s(vec_M_SVij12, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   case (1) ! for sp3 basis set:
      allocate(vec_M_SVij12(4))
      allocate(vec_M_SVij21(4))
      call BOP_matrices_to_vec(basis_ind, M_Vij, vec_M_SVij12, vec_M_SVij21)   ! below
      call KS_sp3_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   case default ! for sp3d5 basis set:
      allocate(vec_M_SVij12(10))
      allocate(vec_M_SVij21(10))
      call BOP_matrices_to_vec(basis_ind, M_Vij, vec_M_SVij12, vec_M_SVij21)   ! below
      call KS_sp3d5_hetero(vec_M_SVij12, vec_M_SVij21, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module "TB_Koster_Slater"
   end select
   if (allocated(vec_M_SVij12)) deallocate(vec_M_SVij12)
   if (allocated(vec_M_SVij21)) deallocate(vec_M_SVij21)
end subroutine Hopping_BOP


subroutine BOP_matrices_to_vec(basis_ind, M_Vij, vec_M_SVij12, vec_M_SVij21)
   integer, intent(in) :: basis_ind
   real(8), dimension(:), intent(in) :: M_Vij
   real(8), dimension(:), intent(out) :: vec_M_SVij12, vec_M_SVij21
   select case (basis_ind)
   case (0) ! for s basis set:
      vec_M_SVij12(1) = M_Vij(1)    ! (s s sigma)
   case (1) ! for sp3 basis set:
      vec_M_SVij12(1) = M_Vij(1)    ! (s s sigma)
      vec_M_SVij21(1) = M_Vij(1)    ! (s s sigma)
      vec_M_SVij12(2) = M_Vij(2)    ! (s p sigma)
      vec_M_SVij21(2) = -M_Vij(3)    ! (p s sigma)
      vec_M_SVij12(3) = M_Vij(4)    ! (p p sigma)
      vec_M_SVij21(3) = M_Vij(4)    ! (p p sigma)
      vec_M_SVij12(4) = M_Vij(5)    ! (p p pi)
      vec_M_SVij21(4) = M_Vij(5)    ! (p p pi)
   case (2) ! for sp3 basis set:
      vec_M_SVij12(1) = M_Vij(1)    ! (s s sigma)
      vec_M_SVij21(1) = M_Vij(1)    ! (s s sigma)
      vec_M_SVij12(2) = M_Vij(2)    ! (s p sigma)
      vec_M_SVij21(2) = -M_Vij(3)    ! (p s sigma)
      vec_M_SVij12(3) = M_Vij(4)    ! (s d sigma)
      vec_M_SVij21(3) = M_Vij(5)    ! (d s sigma)
      vec_M_SVij12(4) = M_Vij(6)    ! (p p sigma)
      vec_M_SVij21(4) = M_Vij(6)    ! (p p sigma)
      vec_M_SVij12(5) = M_Vij(7)    ! (p p pi)
      vec_M_SVij21(5) = M_Vij(7)    ! (p p pi)
      vec_M_SVij12(6) = M_Vij(8)    ! (p d sigma)
      vec_M_SVij21(6) = -M_Vij(9)    ! (d p sigma)
      vec_M_SVij12(7) = M_Vij(10)    ! (p d pi)
      vec_M_SVij21(7) = -M_Vij(11)    ! (d p pi)
      vec_M_SVij12(8) = M_Vij(12)    ! (d d sigma)
      vec_M_SVij21(8) = M_Vij(12)    ! (d d sigma) 
      vec_M_SVij12(9) = M_Vij(13)    ! (d d pi)
      vec_M_SVij21(9) = M_Vij(13)    ! (d d pi)
      vec_M_SVij12(10) = M_Vij(14)    ! (d d delta)
      vec_M_SVij21(10) = M_Vij(14)    ! (d d delta)
   end select
   ! Reminder (Vector):
   !Vr(1) = (s s sigma)
   !Vr(2) = (s p sigma)         
   !Vr(3) = (p s sigma)
   !Vr(4) = (s d sigma)
   !Vr(5) = (d s sigma)
   !Vr(6) = (p p sigma)
   !Vr(7) = (p p pi)
   !Vr(8) = (p d sigma)         
   !Vr(9) = (d p sigma)
   !Vr(10) = (p d pi)               
   !Vr(11) = (d p pi)
   !Vr(12) = (d d sigma)         
   !Vr(13) = (d d pi)
   !Vr(14) = (d d delta)
   ! Reminder (Matrices):
   !Vr(ir,1) = (s s sigma)
   !Vr(ir,2) = (s p sigma)
   !Vr(ir,3) = (s d sigma)
   !Vr(ir,4) = (p p sigma)
   !Vr(ir,5) = (p p pi)
   !Vr(ir,6) = (p d sigma)
   !Vr(ir,7) = (p d pi)
   !Vr(ir,8) = (d d sigma)
   !Vr(ir,9) = (d d pi)
   !Vr(ir,10) = (d d delta)
end subroutine BOP_matrices_to_vec



!cccccccccccccccccccccccccccccccccccccccccccccccccc
! Construct matrices that enter Hamiltonian:
subroutine Construct_Vij_BOP(numpar, TB_Hamil, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_E0ij, M_dE0ij)
   type(Numerics_param), intent(in), target :: numpar 	! all numerical parameters
   type(TB_H_BOP), dimension(:,:), intent(in), target :: TB_Hamil	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_Vij	! matrix of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dVij	! matrix of derivatives of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_SVij	! matrix of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dSVij	! matrix of derivatives of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_E0ij     ! matrix of functions for on-site energies for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable, intent(inout) :: M_dE0ij   ! matrix of derivatives of functions for on-site energies for all pairs of atoms, all orbitals
   !----------------------------------
   integer :: j, atom_2, Vr_ind, N_bs, N_onsite
   integer, pointer :: nat, KOA1, KOA2, m, i
   real(8), pointer :: r
   
   ! Total number of atoms in the supercell:
   nat => Scell(NSC)%Na
   ! number of hopping integrals for this basis set in BOP:
   N_bs = identify_BOP_params_size(numpar%N_basis_size)  ! below
   N_onsite = identify_BOP_onsite_size(numpar%N_basis_size)  ! below
   
   ! Allocate the arrays for the overlap parameters:
   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,N_bs), source = 0.0d0)     ! each pair of atoms, all  V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,N_bs), source = 0.0d0)   ! each pair of atoms, all  dV functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,N_bs), source = 0.0d0)   ! each pair of atoms, all S functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,N_bs), source = 0.0d0) ! each pair of atoms, all dS functions
   if (.not.allocated(M_E0ij)) allocate(M_E0ij(nat,nat,N_onsite), source = 0.0d0)   ! all on-site parameters
   if (.not.allocated(M_dE0ij)) allocate(M_dE0ij(nat,nat,N_onsite), source = 0.0d0)   ! all on-site parameters

   ! Set all the parameters:
!$omp parallel
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, Vr_ind)
   AT1:do j = 1, nat
      m => Scell(NSC)%Near_neighbor_size(j) ! how many atoms are close enough to the j-th one
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         i => Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         KOA1 => Scell(NSC)%MDatoms(j)%KOA
         KOA2 => Scell(NSC)%MDatoms(i)%KOA
         r => Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]
         
         ! All on-site parameters:
         Vr_ind = 1 ! s
         M_E0ij(j,i,Vr_ind) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .true., &
                                      TB_Hamil(KOA1,KOA2)%E_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_ni(Vr_ind,:)) ! below
         if (numpar%N_basis_size > 0) then  ! sp3
            do Vr_ind = 2,3 ! p
               M_E0ij(j,i,Vr_ind) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .true., &
                                      TB_Hamil(KOA1,KOA2)%E_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_ni(Vr_ind,:)) ! below
            enddo
         endif   
         if (numpar%N_basis_size > 1) then  ! sp3d5
            do Vr_ind = 4,6 ! d
               M_E0ij(j,i,Vr_ind) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .true., &
                                      TB_Hamil(KOA1,KOA2)%E_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_ni(Vr_ind,:)) ! below
            enddo
         endif
         
         ! All radial functions for Hamiltonian:
         Vr_ind = 1 ! s
         M_Vij(j,i,1) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s s sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
         M_SVij(j,i,1) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s s sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
         select case (numpar%N_basis_size)
         case (1)    ! sp3
            Vr_ind = 2
            M_Vij(j,i,2) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s p sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_SVij(j,i,2) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s s sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            Vr_ind = 3
            M_Vij(j,i,3) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p s sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_SVij(j,i,3) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p s sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            Vr_ind = 6    
            M_Vij(j,i,4) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &    ! (p p sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_SVij(j,i,4) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p p sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            Vr_ind = 7
            M_Vij(j,i,5) =  BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &    ! (p p pi)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_SVij(j,i,5) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p p pi)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
         case default    ! sp3d5
            do Vr_ind = 2, 14
               M_Vij(j,i,Vr_ind) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
               M_SVij(j,i,Vr_ind) = BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            enddo
         endselect

         !ddddddddddddddddddddddddddddddddddd
         ! Derivatives of all on-site parameters:
         Vr_ind = 1 ! s
         M_dE0ij(j,i,Vr_ind) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .true., &
                                      TB_Hamil(KOA1,KOA2)%E_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_ni(Vr_ind,:)) ! below
         if (numpar%N_basis_size > 0) then  ! sp3
            do Vr_ind = 2,3 ! p
               M_dE0ij(j,i,Vr_ind) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .true., &
                                      TB_Hamil(KOA1,KOA2)%E_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_ni(Vr_ind,:)) ! below
            enddo
         endif   
         if (numpar%N_basis_size > 1) then  ! sp3d5
            do Vr_ind = 4,6 ! d
               M_dE0ij(j,i,Vr_ind) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .true., &
                                      TB_Hamil(KOA1,KOA2)%E_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%E_ni(Vr_ind,:)) ! below
            enddo
         endif
         
         ! Derivatives of all radial functions for Hamiltonian:
         Vr_ind = 1 ! s
         M_dVij(j,i,1) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s s sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
         M_dSVij(j,i,1) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s s sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
         select case (numpar%N_basis_size)
         case (1)    ! sp3
            Vr_ind = 2
            M_dVij(j,i,2) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s p sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_dSVij(j,i,2) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (s s sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            Vr_ind = 3
            M_dVij(j,i,3) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p s sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_dSVij(j,i,3) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p s sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            Vr_ind = 6    
            M_dVij(j,i,4) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &    ! (p p sigma)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_dSVij(j,i,4) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p p sigma)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            Vr_ind = 7
            M_dVij(j,i,5) =  d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &    ! (p p pi)
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
            M_dSVij(j,i,5) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &   ! (p p pi)
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
         case default    ! sp3d5
            do Vr_ind = 2, 14
               M_dVij(j,i,Vr_ind) =  d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &
                                      TB_Hamil(KOA1,KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%H_ni(Vr_ind,:)) ! below
               M_dSVij(j,i,Vr_ind) = d_BOP_radial_function(r, TB_Hamil(KOA1,KOA2)%rcut, TB_Hamil(KOA1,KOA2)%dcut, .false., &
                                      TB_Hamil(KOA1,KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1,KOA2)%S_ni(Vr_ind,:)) ! below
            enddo
         endselect
         
      enddo AT2
   enddo AT1
!$omp end do
!$omp end parallel
   
   nullify(nat, KOA1, KOA2, m, i, r)
end subroutine Construct_Vij_BOP


!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
! Radial functions:

function BOP_radial_function(r, rcut, dcut, on_site, ci, li, ni) result(BOP_func)
   real(8) BOP_func
   real(8), intent(in) :: r, rcut, dcut
   logical, intent(in) :: on_site
   real(8), dimension(:), intent(in) :: ci, li, ni
   real(8) :: Rmax, pure_f, cut_f

   Rmax = rcut + dcut  ! [A]
   if (R >= Rmax) then ! too far apart, don't interact
      pure_f = 0.0d0
      cut_f = 0.0d0
   else if (R < rcut) then ! no smooth cut-off
      pure_f = BOP_pure_f(r, on_site, ci, li, ni) ! below
      cut_f = 1.0d0
   else    ! smooth cut-off region
      pure_f = BOP_pure_f(r, on_site, ci, li, ni) ! below
      cut_f = BOP_cut_off(rcut, dcut, r)    ! below
   endif

   BOP_func = pure_f * cut_f
end function BOP_radial_function



function BOP_pure_f(r, on_site, ci, li, ni) result(pure_f)
   real(8) pure_f
   real(8), intent(in) :: r
   logical, intent(in) :: on_site
   real(8), dimension(:), intent(in) :: ci, li, ni
   real(8) :: eps, arg, fnctn
   integer :: i
   eps = 1.0d-10    ! precision
   pure_f = 0.0d0   ! to start with
   if (abs(ci(1)) < eps) then ! no such level exists, shift it away
      if (on_site) then ! account for energylevel shifts
         pure_f = 100.0d0   ! [eV] chosen arbitrarily, similar to NRL
      endif
   else  ! real level
      do i = 1, 7    ! all BOP fitting functions
         if ( abs(ci(i)) < eps ) exit    ! exclude zeros
         if ( (abs(R) < eps) .or. (abs(li(i)) < eps ) ) then
            fnctn = ci(i)
         else
            arg = li(i) * R**ni(i)
            fnctn = ci(i) * exp(-arg)
         endif
         pure_f = pure_f + fnctn
      enddo
   endif
end function BOP_pure_f



function d_BOP_radial_function(r, rcut, dcut, on_site, ci, li, ni) result(BOP_func)
   real(8) BOP_func
   real(8), intent(in) :: r, rcut, dcut
   logical, intent(in) :: on_site
   real(8), dimension(:), intent(in) :: ci, li, ni
   real(8) :: Rmax, pure_f, cut_f, d_pure_f, d_cut_f

   Rmax = rcut + dcut  ! [A]
   if (R >= Rmax) then ! too far apart, don't interact
      pure_f = 0.0d0
      d_pure_f = 0.0d0
      cut_f = 0.0d0
      d_cut_f = 0.0d0
   else if (R < rcut) then ! no smooth cut-off
      pure_f = BOP_pure_f(r, on_site, ci, li, ni) ! below
      d_pure_f = d_BOP_pure_f(r, on_site, ci, li, ni) ! below
      cut_f = 1.0d0
      d_cut_f = 0.0d0
   else    ! smooth cut-off region
      pure_f = BOP_pure_f(r, on_site, ci, li, ni) ! below
      d_pure_f = d_BOP_pure_f(r, on_site, ci, li, ni) ! below
      cut_f = BOP_cut_off(rcut, dcut, r)    ! below
      d_cut_f = d_BOP_cut_off(rcut, dcut, r)    ! below
   endif

   BOP_func = pure_f * d_cut_f + d_pure_f * cut_f
end function d_BOP_radial_function



function d_BOP_pure_f(r, on_site, ci, li, ni) result(pure_f)
   real(8) pure_f
   real(8), intent(in) :: r
   logical, intent(in) :: on_site
   real(8), dimension(:), intent(in) :: ci, li, ni
   real(8) :: eps, arg, fnctn
   integer :: i
   eps = 1.0d-10    ! precision
   pure_f = 0.0d0   ! to start with
   if (abs(ci(1)) > eps) then ! no such level exists, shift it away
      do i = 1, 7    ! all BOP fitting functions
         if ( abs(ci(i)) < eps ) exit    ! exclude zeros
         if ( (abs(R) < eps) .or. (abs(li(i)) < eps ) ) then
            fnctn = 0.0d0
         else
            arg = li(i) * R**(ni(i)-1.0d0)
            fnctn = - ci(i) * arg * ni(i) * exp(-arg * R)
         endif
         pure_f = pure_f + fnctn
      enddo
   endif
end function d_BOP_pure_f



pure function BOP_cut_off(rcut, dcut, R) result(F_cut)
   real(8) F_cut
   real(8), intent(in) :: rcut, dcut   ! coeffs. of Eq.(19) [1]
   real(8), intent(in) ::  R ! interatomic distance [A]
   real(8) :: arg
   arg = (R - rcut - dcut)/dcut
   F_cut = 0.5d0 * ( 1.0d0 - cos( g_Pi * arg) )
end function BOP_cut_off


pure function d_BOP_cut_off(rcut, dcut, R) result(F_cut)
   real(8) F_cut
   real(8), intent(in) :: rcut, dcut   ! coeffs. of Eq.(19) [1]
   real(8), intent(in) ::  R ! interatomic distance [A]
   real(8) :: arg
   arg = (R - rcut - dcut)/dcut
   F_cut = 0.5d0*g_Pi/dcut * sin( g_Pi * arg )
end function d_BOP_cut_off


!uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
! Utilities:
pure function identify_BOP_onsite_size(ind) result (N)
   integer, intent(in) :: ind   ! index of the basis set type
   integer :: N ! number of on-site energies
   select case (ind)
   case (0) ! s
      N = 1
   case (1) ! sp3
      N = 3
   case default ! sp3d5
      N = 6
   end select
end function identify_BOP_onsite_size


pure function identify_BOP_params_size(ind) result (N)
   integer, intent(in) :: ind   ! index of the basis set type
   integer :: N ! number of radial overlap functions
   select case (ind)
   case (0) ! s
      N = 1
   case (1) ! sp3
      N = 5
   case default ! sp3d5
      N = 14
   end select
end function identify_BOP_params_size


!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Repulsive part of TB:

subroutine get_Erep_s_BOP(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_BOP), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   integer :: i1, m, atom_2, j1
   integer, pointer :: KOA1, KOA2
   real(8), pointer :: r

   a = 0.0d0    ! to start with
!$omp parallel private(i1, m, atom_2, j1, KOA1, KOA2, r)
!$omp do reduction( + : a)
   do i1 = 1, Scell(NSC)%Na
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if (j1 /= i1) then
            KOA1 => Scell(NSC)%MDatoms(i1)%KOA
            KOA2 => Scell(NSC)%MDatoms(j1)%KOA
            r => Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            a = a + BOP_repulsive_one(TB_Repuls(KOA1,KOA2), r)    ! below
!             print*, 'get_Erep_s_BOP', r, a
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1
!$omp end do
!$omp end parallel
   a = a/2.0d0 ! it was doubled
   nullify(KOA1, KOA2, r)
end subroutine get_Erep_s_BOP


function BOP_repulsive_one(TB_Repuls, r) result(Erep)
   real(8) Erep ! [eV]
   type(TB_Rep_BOP), intent(in) :: TB_Repuls
   real(8), intent(in) :: r
   integer :: Nsiz, icur
   Nsiz = size(TB_Repuls%R)
   if (r >= TB_Repuls%R(Nsiz)) then ! too far, no interaction
      Erep = 0.0d0
   elseif (r < TB_Repuls%R(1)) then ! too close, extrapolate the first point data
      icur = 1
      call linear_interpolation(TB_Repuls%R, TB_Repuls%V_rep, r, Erep, icur, x0=0.0d0, y0=1.0d5) ! module "Little_subroutines"
   else ! within the range, interpolate numerical data
      call Find_in_array_monoton(TB_Repuls%R, r, icur) ! module "Little_subroutines"
      call linear_interpolation(TB_Repuls%R, TB_Repuls%V_rep, r, Erep, icur) ! module "Little_subroutines"
   endif
end function BOP_repulsive_one


subroutine dErdr_s_BOP(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   type(TB_Rep_BOP), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer ian, n
   n = size(atoms)
   !$omp PARALLEL DO private(ian)
   do ian = 1, n  ! Forces for all atoms
     Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
   enddo ! ian
   !$OMP END PARALLEL DO
END subroutine dErdr_s_BOP


subroutine dErdr_Pressure_s_BOP(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_BOP), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms ! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      Scell(NSC)%SCforce%rep = 0.0d0
   endif
end subroutine dErdr_Pressure_s_BOP

END MODULE TB_BOP
