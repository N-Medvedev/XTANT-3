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
! This module contains subroutines to deal with TB hamiltonian in the NRL parametrization
! This parameterization is described in:
! [1] D.A. Papaconstantopoulos and M.J. Mehl, J. Phys.: Condens. Matter 15 (2003) R413–R440

MODULE TB_NRL
use Universal_constants
use TB_Koster_Slater
use Objects
!use Variables
! use Algebra_tools
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity
use Little_subroutines
use Atomic_tools
use Electron_tools
use Nonadiabatic

implicit none


! ! this interface finds by itself which of the two subroutine to use depending on the array passed:
! interface test_orthogonalization
!    module procedure test_orthogonalization_r	! real version
!    module procedure test_orthogonalization_c	! complex version
! end interface test_orthogonalization


! public :: test_orthogonalization

! Modular parameters:
real(8) :: m_one_third, m_two_third, m_four_third

parameter (m_one_third = 1.0d0/3.0d0)
parameter (m_two_third = 2.0d0/3.0d0)
parameter (m_four_third = 2.0d0*m_two_third)

 contains

!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Attractive part of TB:

! Tight Binding Hamiltonian within NRL parametrization:
subroutine construct_TB_H_NRL(numpar, matter, TB_Hamil, M_Vij, M_SVij, M_lmn, Scell, NSC, Err)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! materil parameters
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   type(Super_cell), dimension(:), intent(inout) :: Scell		! supercell with all the atoms as one object
   integer, intent(in) :: NSC		! number of supercell
   type(Error_handling), intent(inout) :: Err	! error save
   character(200) :: Error_descript
   real(8) :: R_cut
   integer i, Nx, Ny, Nz
   Error_descript = ''

!$OMP WORKSHARE
   Scell(NSC)%Ha0 = Scell(NSC)%Ha	! save Hamiltonial from previous time-step
   Scell(NSC)%Ei0 = Scell(NSC)%Ei	! save energy levels for the next timestep
   Scell(NSC)%H_non0 = Scell(NSC)%H_non	! save non-diagonalized Hamiltonian from last time-step
!$OMP END WORKSHARE

    ! Construct TB Hamiltonian (with Mehl parameters),  orthogonalize it,  and then solve generalized eigenvalue problem:
   call Hamil_tot_NRL(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err)	! see below

   ! Check whether the supercell is large enough:
!    call get_near_neighbours(Scell, dm=R_cut) ! module "Atomic_tools"
!    if ( (R_cut <= Scell(NSC)%Supce(1,1)/2.0d0) .and. (R_cut <= Scell(NSC)%Supce(2,2)/2.0d0) .and. (R_cut <= Scell(NSC)%Supce(3,3)/2.0d0) ) then
!       call Hamil_tot_NRL(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err)	! see below
! !       call Hamil_tot_NRL_unitcell(numpar, Scell, NSC, TB_Hamil, Err)
!    else
!       call Hamil_tot_NRL_unitcell(numpar, Scell, NSC, TB_Hamil, Err)
! !       call Hamil_tot_NRL(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err)	! see below
!    endif

   ! test (comment out for release):
!    call test_nonorthogonal_solution(Scell(NSC))

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_NRL


subroutine test_nonorthogonal_solution(Scell)
! The subroutine tests if solution of the linear eigenproblem (Schoedinger eq. with overlap matrix) was correct.
! Note that the subroutine is very slow, so use it only for testing, not for actual calculations.
   type(Super_cell), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   real(8) Ev, Evec(size(Scell%Ei)), Svec(size(Scell%Ei)), Sv, epsylon, SM(size(Scell%Ei),size(Scell%Ei)), EvM(size(Scell%Ei),size(Scell%Ei)), OvM(size(Scell%Ei),size(Scell%Ei))
   integer i, j, k, M
!    epsylon = 1.0d-12
   epsylon = 1.0d-11

   call print_time_step('test_nonorthogonal_solution', 0.0d0, msec=.true.)
   
   
   M = size(Evec)
   Ev = 0.0d0
   do i = 1,M
      do k = 1,M
         Evec(k) = SUM(Scell%H_non(k,:)*Scell%Ha(:,i)) ! correct
         Svec(k) = SUM(Scell%Sij(k,:)*Scell%Ha(:,i)) ! correct
      enddo
      
      Ev = SUM(Evec(:)*Scell%Ha(:,i))
      Sv = SUM(Svec(:)*Scell%Ha(:,i))
      if (ABS(Sv-1.0d0) > epsylon) then ! the normalization broke down:
         write(*,'(a,i5,f,f,f)') 'Sv: ', i, Sv, Ev, Scell%Ei(i)
      endif
      
      Sv = Sv*Scell%Ei(i)
      ! Solution of the linear eigenproblem must be zero:
      if (ABS(Ev - Sv)/min(ABS(Ev),ABS(Sv)) .GT. epsylon) then ! if it's not, something went wrong:
         print*, 'Ev:', i, Ev, Sv
      endif
   enddo

   call print_time_step('test_nonorthogonal_solution', 1.0d0, msec=.true.)
   
!   Matrix multiplication method:
!   Get the < i | H | j >
   CALL dgemm ('T','N', M, M, M, 1.0d0, Scell%H_non, M, Scell%Ha, M, 0.0d0, SM, M)	
   CALL dgemm ('T','N', M, M, M, 1.0d0, SM, M, Scell%Ha, M, 0.0d0, EvM, M)
!   Get the < i | S | j >
   CALL dgemm ('T','N', M, M, M, 1.0d0, Scell%Sij, M, Scell%Ha, M, 0.0d0, SM, M)	
   CALL dgemm ('T','N', M, M, M, 1.0d0, SM, M, Scell%Ha, M, 0.0d0, OvM, M)
!   Check generalized eigenvalue problem:
   do k = 1,M
      do i = 1, M
         if (ABS(EvM(k,i) - Scell%Ei(k)*OvM(k,i)) > epsylon) then
            write(*,'(i5,i5, es,es,es)') k, i, EvM(k,i), Scell%Ei(k)*OvM(k,i)
         endif
      enddo
   enddo

   call print_time_step('test_nonorthogonal_solution', 2.0d0, msec=.true.)
   
!    PAUSE 'test_nonorthogonal_solution'
end subroutine test_nonorthogonal_solution


subroutine Hamil_tot_NRL(numpar, Scell, NSC, TB_Hamil, M_Vij, M_SVij, M_lmn, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   type(Error_handling), intent(inout) :: Err	! error save
   !--------------------------------------
   real(8), dimension(:,:), allocatable :: Hij	 ! Hamiltonian
   real(8), dimension(:,:), allocatable :: Sij, SH, HS
   real(8), dimension(size(Scell(NSC)%Ha,1)) :: Evec, EvecS
   real(8), dimension(9,9) :: Hij1, Sij1
   real(8), dimension(10) :: Vij1
   integer :: nat, Nsiz, n_orb
   integer, target :: j
   integer :: j1, i1, k, l, atom_2, FN, i
   real(8) :: temp, epsylon, Ev, SH_1, SH_2
   real(8), pointer :: x, y, z
   integer, pointer :: KOA1, KOA2, m
   character(200) :: Error_descript
   
!    print*, 'Hamil_tot_NRL     is called'
   
   Error_descript = ''
   epsylon = 1d-12	! acceptable tolerance : how small an overlap integral can be, before we set it to zero
   
   n_orb = 9	! number of orbital for sp3d5 basis
   nat = Scell(NSC)%Na
   Nsiz = size(Scell(NSC)%Ha,1)
   
   if (.not.allocated(Sij)) allocate(Sij(Nsiz,Nsiz))
   if (.not.allocated(Hij)) allocate(Hij(Nsiz,Nsiz))
   Sij = 0.0d0
   Hij = 0.0d0
   !-----------------------------------
   ! 1) Construct non-orthogonal Hamiltonian H and Overlap matrix S in 2 steps:
   
   ! a) Construct upper triangle - calculate each element:
!$omp parallel
!$omp do  private(j, m, atom_2, i, KOA1, KOA2, j1, l, i1, k, Hij1, Sij1)
   do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one  
         if (atom_2 .EQ. 0) then
            i = j
         else
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         endif
         
         IJ:if (i >= j) then ! it's a new pair of atoms, calculate everything
            ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
            call Hamilton_one_NRL(Scell, NSC, i, j, atom_2, TB_Hamil, Hij1, M_Vij, M_lmn) ! this calles the block-hamiltonian

            KOA1 => Scell(NSC)%MDatoms(j)%KOA
            KOA2 => Scell(NSC)%MDatoms(i)%KOA
            ! Construct overlap matrix for this pair of atoms:
            call Get_overlap_S_matrix(Scell, NSC, i, atom_2, j, TB_Hamil(KOA1,KOA2), Sij1, M_SVij, M_lmn)

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


   ! 2) Convert the Hamiltonian on-site energies from Rydbergs into eVs:
   !$OMP WORKSHARE
   Hij = Hij*g_Ry	! [Ry] into [eV]
   ! Save the non-orthogonalized Hamiltonian:
   Scell(NSC)%H_non = Hij		! nondiagonalized Hamiltonian
   Scell(NSC)%Sij = Sij		! save Overlap matrix
   !$OMP END WORKSHARE
   
   !-----------------------------------
   ! 3) Orthogonalize the Hamiltonian using Lowedin procidure
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization(Nsiz, Sij, Hij, Err, Scell(NSC)%eigen_S)	! below
   
   !$OMP WORKSHARE
   Scell(NSC)%Hij = Hij ! save orthogonalized but non-diagonalized Hamiltonian
   !$OMP END WORKSHARE
   
   !-----------------------------------
   ! 4) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, check_M=.true.)
   call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript) ! module "Algebra_tools"
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Hamil_tot_NRL: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif
   Scell(NSC)%Hij_sol = Hij		! eigenvectors of nondiagonalized Hamiltonian
   
   !-----------------------------------
   ! 5) Convert the eigenvectors back into the non-orthogonal basis:
   call mkl_matrix_mult('N', 'N', Sij, Hij, Scell(NSC)%Ha)	! module "Algebra_tools"
   ! Normalize eigenvectors to | <n|n> |^2 = 1:
!    do j = 1, Nsiz
! !       print*, 'Eigen:', j, SUM( Scell(NSC)%Ha(:,j) * Scell(NSC)%Ha(:,j) )
!       temp = DSQRT(SUM( Scell(NSC)%Ha(:,j) * Scell(NSC)%Ha(:,j) ))
!       Scell(NSC)%Ha(:,j) =  Scell(NSC)%Ha(:,j) / temp
!    enddo
   
   !-----------------------------------
!    ! An example of the solution of the linear eigenproblem with LAPACK subroutines:
!    call dpotrf('U', Nsiz, Sij, Nsiz, j)
!    call dsygst(1, 'U', Nsiz, Hij, Nsiz, Sij, Nsiz,  j)
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript)
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
!                         SH_1 =  0.0d0 
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
      !$omp end parallel
!       deallocate(SH, HS)
      ! Convert to SI units used later:
      !temp = g_me*g_e/g_h*1d-10 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
      
      Scell(NSC)%PRRx = Scell(NSC)%PRRx * temp
      Scell(NSC)%PRRy = Scell(NSC)%PRRy * temp
      Scell(NSC)%PRRz = Scell(NSC)%PRRz * temp
   endif ! (numpar%optic_model .EQ. 3) 
   nullify(KOA1, KOA2, m, x, y, z)

   deallocate(Hij, Sij)
end subroutine Hamil_tot_NRL



subroutine Hamil_tot_NRL_unitcell(numpar, Scell, NSC, TB_Hamil, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Error_handling), intent(inout) :: Err	! error save
   !--------------------------------------
   real(8), dimension(:,:), allocatable :: Hij, Sij	 ! Hamiltonian
!    real(8), dimension(:,:), allocatable :: Hij_p, Sij_p
   real(8), dimension(size(Scell(NSC)%Ha,1)) :: Evec, EvecS
   real(8), dimension(9,9) :: Hij1, Sij1
   real(8), dimension(10) :: Vij1
   integer :: nat, Nsiz, n_orb, zb(3)
   integer, target :: j
   integer :: j1, i1, k, l, atom_2, FN, i
   integer :: Nx, Ny, Nz, x_cell, y_cell, z_cell
   real(8) :: temp, epsylon, Ev, SH_1, SH_2, Ncells
   real(8), pointer :: x, y, z
   integer, pointer :: KOA1, KOA2, m
   character(200) :: Error_descript
   logical :: origin_cell
   
   Error_descript = ''
   epsylon = 1d-12	! acceptable tolerance : how small an overlap integral can be, before we set it to zero
   
   n_orb = 9	! number of orbital for sp3d5 basis
   nat = Scell(NSC)%Na
   Nsiz = size(Scell(NSC)%Ha,1)
   
   if (.not.allocated(Sij)) allocate(Sij(Nsiz,Nsiz))
   if (.not.allocated(Hij)) allocate(Hij(Nsiz,Nsiz))
   Sij = 0.0d0
   Hij = 0.0d0

   !-----------------------------------
   ! 0) Get the number of cells to be included:
   call get_mirror_cell_num_NRL(Scell, NSC, numpar, Nx, Ny, Nz)	! see below

   ! 1) Construct non-orthogonal Hamiltonian H and Overlap matrix S in 2 steps:
   
   ! a) Construct upper triangle - calculate each element:
   call Cycles_to_get_H1_and_S1(Scell, NSC, nat, n_orb, Nx, Ny, Nz, TB_Hamil, numpar, Hij1, Sij1, Hij, Sij)	! below

   ! !$omp parallel  shared(Hij, Sij)
! print*, 'TEST PARALLEL', allocated(Hij), allocated(Sij)
!    if (.not.allocated(Hij_p)) allocate(Hij_p(Nsiz,Nsiz))
!    if (.not.allocated(Sij_p)) allocate(Sij_p(Nsiz,Nsiz))
!    Sij_p = 0.0d0
!    Hij_p = 0.0d0
! print*, 'TEST _p', allocated(Hij_p), allocated(Sij_p)
! !$omp do reduction( + : Hij_p, Sij_p)  private(j, i, x_cell, y_cell, z_cell, zb, j1, l, i1, k, Hij1, Sij1)
!    do j = 1,nat	! all atoms
!    print*, 'TEST DO PARALLEL', allocated(Hij), allocated(Sij)
!       do i = j,nat	! all pairs of atoms
!          print*, 'TEST DO PARALLEL2', allocated(Hij), allocated(Sij)
!          XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
!             YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
!                ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
!                   zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
! 
!                   print*, j,i, 'Before Hij1' , allocated(Hij), allocated(Sij)
!                   
!                   ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
!                   call Hamilton_one_NRL_cell(Scell, NSC, j, i, zb, TB_Hamil, Hij1)	! this calles the block-hamiltonian
!                   print*, j,i, 'Hij1'
!                   ! Construct overlap matrix for this pair of atoms:
!                   call Get_overlap_S_matrix_cell(Scell, NSC, j, i, zb, TB_Hamil, Sij1)
!                   print*, j,i, 'Sij1'
! 
!                   do j1 = 1,n_orb ! all orbitals
!                      l = (j-1)*n_orb+j1
!                      do i1 = 1,n_orb ! all orbitals
!                         k = (i-1)*n_orb+i1
!                         Hij_p(k,l) = Hij_p(k,l) + Hij1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
!                         if (ABS(Sij1(i1,j1)) >= epsylon) then
!                            Sij_p(k,l) = Sij_p(k,l) + Sij1(i1,j1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
!                         endif
!                      enddo ! i1
!                   enddo ! j1
!                enddo ZC
!             enddo YC
!          enddo XC
!       enddo ! i
!    enddo ! j
! !$omp end do
!    Hij = Hij_p
!    Sij = Sij_p
!    deallocate(Hij_p, Sij_p)
! !$omp end parallel

   ! b) Construct lower triangle - use symmetry:
!$omp parallel
!$omp do  private(j, i, j1, l, i1, k)
   do j = 1,nat	! all atoms
      do i = 1,j	! do only for atoms close to that one  
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

   ! 2) Convert the Hamiltonian on-site energies from Rydbergs into eVs:
   !$OMP WORKSHARE
   Hij = Hij*g_Ry	! [Ry] into [eV]
   ! Save the non-orthogonalized Hamiltonian:
   Scell(NSC)%H_non = Hij		! nondiagonalized Hamiltonian
   Scell(NSC)%Sij = Sij		! save Overlap matrix
   !$OMP END WORKSHARE
   
   !-----------------------------------
   ! 3) Orthogonalize the Hamiltonian using Lowedin procidure
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization(Nsiz, Sij, Hij, Err, Scell(NSC)%eigen_S) ! below
   
   !$OMP WORKSHARE
   where (ABS(Hij) < epsylon) Hij = 0.0d0
   Scell(NSC)%Hij = Hij ! save orthogonalized but non-diagonalized Hamiltonian
   !$OMP END WORKSHARE
   
   !-----------------------------------
   ! 4) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, check_M=.true.)
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript)
   call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript, use_DSYEV=.true.)
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Hamil_tot_NRL: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif
   Scell(NSC)%Hij_sol = Hij		! eigenvectors of nondiagonalized Hamiltonian
   
   !-----------------------------------
   ! 5) Convert the eigenvectors back into the non-orthogonal basis:
   call mkl_matrix_mult('N', 'N', Sij, Hij, Scell(NSC)%Ha)	! module "Algebra_tools"
   
   !-----------------------------------
!    ! An example of the solution of the linear eigenproblem with LAPACK subroutines:
!    call dpotrf('U', Nsiz, Sij, Nsiz, j)
!    call dsygst(1, 'U', Nsiz, Hij, Nsiz, Sij, Nsiz,  j)
!    call sym_diagonalize(Hij, Scell(NSC)%Ei, Error_descript)
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
                        if (j1 == i1) then
                           ! Skip the same orbital, no overlap
                        else
                           ! [1] Nonorthagonal expression:
                           !SH_1 =  0.270d0*(Scell(NSC)%H_non(k,l) - Scell(NSC)%Ei(k) * Scell(NSC)%Sij(k,l))
                           ! Testing alternative orthogonalized Hamiltonian:
                           SH_1 = 1.10d0*Scell(NSC)%Hij(k,l)
                           !Testing orthogonal expression (Trani):
                           !SH_1 = 0.270d0*Scell(NSC)%H_non(k,l)

                           Scell(NSC)%PRRx(k,l) = SH_1
                           Scell(NSC)%PRRy(k,l) = SH_1
                           Scell(NSC)%PRRz(k,l) = SH_1
                        endif
                     else
                        ! [1] Nonorthagonal expression:
                        SH_1 = (Scell(NSC)%H_non(k,l) - Scell(NSC)%Ei(k) * Scell(NSC)%Sij(k,l))
                        ! Testing alternative orthogonalized Hamiltonian:
                        !SH_1 = Scell(NSC)%Hij(k,l)
                        !Testing orthogonal expression (Trani):
                        !SH_1 = Scell(NSC)%H_non(k,l)

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
!       deallocate(SH, HS)
      ! Convert to SI units used later:
      !temp = g_me*g_e/g_h*1d-10 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0 !/ 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
      
      Scell(NSC)%PRRx = Scell(NSC)%PRRx * temp
      Scell(NSC)%PRRy = Scell(NSC)%PRRy * temp
      Scell(NSC)%PRRz = Scell(NSC)%PRRz * temp
   endif ! (numpar%optic_model .EQ. 3) 
   nullify(KOA1, KOA2, m, x, y, z)

   deallocate(Hij, Sij)
end subroutine Hamil_tot_NRL_unitcell


subroutine Cycles_to_get_H1_and_S1(Scell, NSC, nat, n_orb, Nx, Ny, Nz, TB_Hamil, numpar, Hij1, Sij1, Hij_p, Sij_p)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   integer, intent(in) :: nat, n_orb, Nx, Ny, Nz
   real(8), dimension(:,:), intent(inout) :: Hij1, Sij1, Hij_p, Sij_p
   !-------------------------
   integer :: x_cell, y_cell, z_cell, zb(3)
   integer :: j, i, j1, l, i1, k

   !$omp parallel private(j, i, x_cell, y_cell, z_cell, zb, j1, l, i1, k, Hij1, Sij1)
   !$omp do reduction( + : Hij_p, Sij_p)
   XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
      YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
         ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
            zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
            do j = 1,nat	! all atoms
               do i = j,nat	! all pairs of atoms


                  ! First, for the non-orthagonal Hamiltonian for this pair of atoms:
                  call Hamilton_one_NRL_cell(Scell, NSC, j, i, zb, TB_Hamil, numpar, Hij1)	! this calles the block-hamiltonian
                  ! Construct overlap matrix for this pair of atoms:
                  call Get_overlap_S_matrix_cell(Scell, NSC, j, i, zb, TB_Hamil, Sij1)

                  do j1 = 1,n_orb ! all orbitals
                     l = (j-1)*n_orb+j1
                     do i1 = 1,n_orb ! all orbitals
                        k = (i-1)*n_orb+i1
                        Hij_p(k,l) = Hij_p(k,l) + Hij1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                        if (ABS(Sij1(i1,j1)) >= 1.0d-10) Sij_p(k,l) = Sij_p(k,l) + Sij1(i1,j1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
                     enddo ! i1
                  enddo ! j1
               enddo ! i
            enddo ! j
         enddo ZC
      enddo YC
   enddo XC
   !$omp end do
   !$omp end parallel
end subroutine Cycles_to_get_H1_and_S1



subroutine get_mirror_cell_num_NRL(Scell, NSC, numpar, Nx, Ny, Nz)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   integer, intent(out) :: Nx, Ny, Nz ! number of super-cells images that contribute to total energy
   real(8), dimension(3) :: zb ! super cell indices
   real(8) :: R_cut

   ! Get the cut-off distance for the potential:
   call get_near_neighbours(Scell, numpar, dm=R_cut) ! module "Atomic_tools"
!    R_cut = R_cut*2.0d0	! to make sure atoms fit into half of the supercell
   
   call get_number_of_image_cells(Scell, NSC, Scell(NSC)%MDAtoms, R_cut, Nx, Ny, Nz) ! module "Atomic_tools"
end subroutine get_mirror_cell_num_NRL



subroutine Loewdin_Orthogonalization(Nsiz, Sij, Hij, Err, eigen_S) ! below
   integer, intent(in) :: Nsiz
   real(8), dimension(:,:), intent(inout) :: Sij, Hij
   type(Error_handling), intent(inout) :: Err	! error save
   real(8), dimension(:), intent(out), optional :: eigen_S   ! eigenvalues of S
   !------------------------------
   real(8), dimension(:,:), allocatable :: Xij, Sij_temp !, Sij_save
   real(8), dimension(:,:), allocatable :: s_mat
   real(8), dimension(:), allocatable :: Ev
   integer :: N_neg, N_zero, i1
   real(8) :: epsylon
   character(200) :: Error_descript
   Error_descript = ''
   epsylon = 1d-8
   
   ! Save the overlap matrix to test the orthogonalization later:
!     allocate(Sij_save(Nsiz,Nsiz))
!     Sij_save = Sij
   allocate(Xij(Nsiz,Nsiz))
   allocate(Sij_temp(Nsiz,Nsiz))
   
   ! Orthogonalize the Hamiltonian, based on  [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144], in a few steps:
   ! 1) diagonalize S matrix:
   allocate(Ev(Nsiz))	! array with eigenvalues of nonorthogonal matrix
   Ev = 0.0d0	! to start from
!    call sym_diagonalize(Sij, Ev, Error_descript) ! module "Algebra_tools"
   call sym_diagonalize(Sij, Ev, Error_descript, use_DSYEV=.true.) ! module "Algebra_tools"
   
   ! If requested, output the eigenvalues of the overlap matrix:
   if (present(eigen_S)) eigen_S = Ev

   ! Now Sij is the collection of eigenvectors, Ev contains eigenvalues
   ! Moreover, Sij is now unitary matrix that can be used as U from Eq.(3.166) [Szabo "Modern Quantum Chemistry" 1986, p. 143]
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Loewdin_Orthogonalization: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif
   ! Test if diagonalization of S matrix worked well:
   !$OMP WORKSHARE
   N_neg = COUNT(Ev < 0.0d0)
   N_zero = COUNT(ABS(Ev) < epsylon)
   !$OMP END WORKSHARE
   if ( (N_zero > 0) .or. (N_neg > 0) ) then
      print*, 'Subroutine Loewdin_Orthogonalization has problems in the overlap matrix:'
      print*, 'Negative: ', N_neg, 'Zeros: ', N_zero
!       print*, 'Above one:', COUNT(Ev > 1.0d0)
!       do i1 = 1, size(Ev)
!          if (Ev(i1) <= epsylon) print*, 'NEG:', Ev(i1)
!          if (Ev(i1) > 1.0d0 + epsylon) print*, 'BIG:', Ev(i1)
!       enddo
      where(Ev(:)<0.0d0) Ev(:) = ABS(Ev(:))
   endif

   ! 2) construct S = [1/sqrt(Ev)]:
   allocate(s_mat(Nsiz,Nsiz)) ! temporary array with matrix of 1/sqrt(Ev)
   !$OMP WORKSHARE
   s_mat = 0.0d0 ! to start from
   ! Construct S matrix if eigenstates are non-zero:
   forall (i1=1:size(Ev), ABS(Ev(i1)) > epsylon) s_mat(i1,i1) = 1.0d0/dsqrt(Ev(i1))
   !$OMP END WORKSHARE

   ! 3) construct X = S*s_mat [Szabo "Modern Quantum Chemistry" 1986, p. 144, Eq.(3.169)]
   call mkl_matrix_mult('N', 'N', Sij, s_mat, Xij) ! module "Algebra_tools"
   
   ! 4) orthoginalize the Hij matrix [Szabo "Modern Quantum Chemistry" 1986, p. 144, Eq.(3.177)]:
   call dsymm('L', 'U', Nsiz, Nsiz, 1.0d0, Hij, Nsiz, Xij, Nsiz, 0.0d0, Sij_temp, Nsiz)	! from BLAS library
   
   call mkl_matrix_mult('T', 'N', Xij, Sij_temp, Hij)	! module "Algebra_tools"
   
!     ! test (comment out for release):
!       call test_orthogonalization_r(Sij_save, Sij, Xij) ! Check  if it is orthogonalized correctly
   
   Sij = Xij	! save X matrix
end subroutine Loewdin_Orthogonalization



subroutine Complex_Hamil_NRL(numpar, Scell, NSC, CHij, CSij, Ei, ksx, ksy, ksz, Err)
! This subroutine is unused for CDF calculations! A newer one is in the module "TB"
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(inout):: Ei	! [eV] energy levels for this complex Hamiltonian
   real(8), intent(in) :: ksx, ksy, ksz ! k-point to get Hamiltonian [relative]
   complex, DIMENSION(:,:), allocatable, INTENT(inout) :: CHij, CSij ! Complex Hamiltonian at k-point and overlap matrix
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------------------------------------------------------
   complex, dimension(:,:), allocatable :: CHij_temp, CHij_non, CSij_save, CHij_orth
   integer :: Nsiz, j, nat, m, atom_2, i, j1, l, i1, k, norb
   real(8) :: kx, ky, kz, temp
   real(8), target :: nol
   real(8), pointer :: x1, y1, z1
   complex(8) :: expfac, SH_1
   character(200) :: Error_descript
   Error_descript = ''
   nol = 0.0d0
!    temp = g_me*g_e/g_h*1d-10
   nat = size(Scell(NSC)%MDatoms) ! number of atoms
   norb = 9	! sp3d5 number of orbitals per atom
   Nsiz = size(Scell(NSC)%Ha,1)

   ! Allocate complex parameters:
   if (.not.allocated(Scell(NSC)%cPRRx)) allocate(Scell(NSC)%cPRRx(Nsiz,Nsiz))
   if (.not.allocated(Scell(NSC)%cPRRy)) allocate(Scell(NSC)%cPRRy(Nsiz,Nsiz))
   if (.not.allocated(Scell(NSC)%cPRRz)) allocate(Scell(NSC)%cPRRz(Nsiz,Nsiz))
   Scell(NSC)%cPRRx = 0.0d0
   Scell(NSC)%cPRRy = 0.0d0
   Scell(NSC)%cPRRz = 0.0d0
   if (.not.allocated(CHij)) allocate(CHij(Nsiz,Nsiz))
   if (.not.allocated(CSij)) allocate(CSij(Nsiz,Nsiz))
   if (.not.allocated(CSij_save)) allocate(CSij_save(Nsiz,Nsiz))
   CHij = dcmplx(0.0d0,0.0d0)	! to start with
   CSij = dcmplx(0.0d0,0.0d0)	! to start with
   CSij_save = dcmplx(0.0d0,0.0d0)	! to start with
   if (.not.allocated(CHij_temp)) allocate(CHij_temp(Nsiz,Nsiz))  ! nonorthogonalized Hamiltonian
   if (.not.allocated(CHij_non)) allocate(CHij_non(Nsiz,Nsiz))
   if (.not.allocated(CHij_orth)) allocate(CHij_orth(Nsiz,Nsiz))  ! orthogonalized Hamiltonian
   CHij_non = dcmplx(0.0d0,0.0d0)	! to start with
   CHij_temp = dcmplx(0.0d0,0.0d0)	! to start with
   CHij_orth = dcmplx(0.0d0,0.0d0)	! to start with
   
   call Reciproc(Scell(NSC)%supce, Scell(NSC)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
   call Reciproc_rel_to_abs(ksx, ksy, ksz, Scell, NSC, kx, ky, kz) ! get absolute k-values

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
         
         if (i > 0) then
            if ((abs(kx) < 1.0d-14) .AND. (abs(ky) < 1.0d-14) .AND. (abs(kz) < 1.0d-14)) then
               expfac = dcmplx(1.0d0,0.0d0)
            else
               expfac = exp(g_CI*dcmplx(kx*x1 + ky*y1 + kz*z1,0.0d0))
            endif

            do j1 = 1,norb ! all orbitals for sp3d5
               l = (j-1)*norb+j1
               do i1 = 1,norb ! all orbitals for sp3d5
                  k = (i-1)*norb+i1

                  CHij_temp(k,l) = DCMPLX(Scell(NSC)%H_non(k,l),0.0d0)*expfac
                  CSij(k,l) = DCMPLX(Scell(NSC)%Sij(k,l),0.0d0)*expfac

!                   write(*,'(i5,i5,es,es,es)') k, l, CSij(k,l), Scell(NSC)%Sij(k,l)
                  
                  if ((isnan(real(CHij_temp(k,l)))) .OR. isnan(aimag(CHij_temp(k,l)))) then
                     print*, i, j, k, l, CHij_temp(k,l)
                     pause 'CHij_temp ISNAN in Complex_Hamil_NRL'
                  endif
                  if ((isnan(real(CSij(k,l)))) .OR. isnan(aimag(CSij(k,l)))) then
                     print*, i, j, k, l, CSij(k,l)
                     pause 'CSij ISNAN in Complex_Hamil_NRL'
                  endif
                  
               enddo ! i1
            enddo ! j1
         endif ! (i > 0)
      enddo ! atom_2
   enddo ! j
   !$omp end do 
   !$omp end parallel

   ! Temporarily save nonorthogonal Hamiltonian and overlap matrix:
   CHij_non = CHij_temp
   CSij_save = CSij

   !call sym_diagonalize(CHij, Ei, Err%Err_descript) ! modeule "Algebra_tools"
   ! Solve linear eigenproblem:
   ! 1) Orthogonalize the Hamiltonian using Loewdin procidure:
   ! according to [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144]:
   call Loewdin_Orthogonalization_c(Nsiz, CSij, CHij_temp, Err)	! below
   ! Save orthogonalized Hamiltonian (for optical coefficients below)
   CHij_orth = CHij_temp

   ! 2) Diagonalize the orthogonalized Hamiltonian to get electron energy levels (eigenvalues of H):
   call sym_diagonalize(CHij_temp, Ei, Error_descript)
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Complex_Hamil_NRL: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif
   ! 3) Convert the eigenvectors back into the non-orthogonal basis:
   call mkl_matrix_mult('N', 'N', CSij, CHij_temp, CHij)	! module "Algebra_tools"
   ! We need to renormalize the wave functions, as they are not normalized to 1 after this procidure:
!    do j = 1, Nsiz
!       CHij(:,j) = CHij(:,j) / DSQRT(dble(SUM( dconjg(CHij(:,j)) * CHij(:,j) )))
!    enddo

   
   
!    do j = 1, size(Ei)
!       write(*,'(i5,es,es,es,es,es)') j, Ei(j), CHij_temp(j,1)
!    enddo
!    PAUSE 'Ei pause'
   
   ! 4) Calculate momentum operators:
   ! Optical matrix elements for non-orthogonal TB are taken from:
   ! [1] arXiv:1805.08918v1 -- https://128.84.21.199/abs/1805.08918
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

                     write(*,'(a,i5,i5,i5,i5)') 'SH_1: ', k, l, i, j

                     if (j1 == i1) then
                        ! Skip the same orbital, no overlap
                     else
                        ! [1] Nonorthagonal expression:
                        !SH_1 = DCMPLX(0.270d0,0.0d0) * (CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l))
                        ! Testing alternative orthogonalized Hamiltonian:
                        SH_1 = DCMPLX(1.1d0,0.0d0) * CHij_orth(k,l)
                        ! Testing orthogonal expression (Trani):
                        !SH_1 = DCMPLX(temp*0.270d0,0.0d0)*CHij_non(k,l)
                        print*, k, l, SH_1

                        Scell(NSC)%cPRRx(k,l) = SH_1
                        Scell(NSC)%cPRRy(k,l) = SH_1
                        Scell(NSC)%cPRRz(k,l) = SH_1
                     endif
                  else	! different atoms at distance {x,y,z}:
                     ! [1] Nonorthagonal expression:
                     SH_1 = CHij_non(k,l) - DCMPLX(Ei(k),0.0d0)*CSij_save(k,l)
                     ! Testing alternative orthogonalized Hamiltonian:
                     !SH_1 = CHij_orth(k,l)
                     ! Testing orthogonal expression (Trani):
                     !SH_1 = CHij_non(k,l)

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
   !$omp end parallel
   
   ! Convert to SI units used later:
   !temp = g_me*g_e/g_h*1d-10 / 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
   ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
   temp = 1.0d0 !/ 2.0d0	! UNCLEAR WHERE THE 1/2 COMES FROM ???
!    Scell(NSC)%cPRRx = Scell(NSC)%cPRRx * temp
!    Scell(NSC)%cPRRy = Scell(NSC)%cPRRy * temp
!    Scell(NSC)%cPRRz = Scell(NSC)%cPRRz * temp
   
   deallocate(CHij_non, CHij_temp, CSij_save)
   nullify(x1, y1, z1)
end subroutine Complex_Hamil_NRL



subroutine Loewdin_Orthogonalization_c(Nsiz, Sij, Hij, Err)	! below
   integer, intent(in) :: Nsiz
   complex, dimension(:,:), intent(inout) :: Sij, Hij
   type(Error_handling), intent(inout) :: Err	! error save
   !------------------------------
   complex, dimension(:,:), allocatable :: Xij, Sij_temp
!    complex, dimension(:,:), allocatable :: Sijsave
   real(8), dimension(:,:), allocatable :: s_mat
   real(8), dimension(:), allocatable :: Ev
   integer :: N_neg, N_zero, i1
   real(8) :: epsylon
   character(200) :: Error_descript
   Error_descript = ''
   epsylon = 1d-12
   
   ! Save the overlap matrix to test the orthogonalization later:
!    allocate(Sijsave(Nsiz,Nsiz))
!    Sijsave = Sij
   allocate(Xij(Nsiz,Nsiz))
   allocate(Sij_temp(Nsiz,Nsiz))
   Xij = dcmplx(0.0d0,0.0d0)
   Sij_temp = dcmplx(0.0d0,0.0d0)
   
   ! Orthogonalize the Hamiltonian, based on  [Szabo "Modern Quantum Chemistry" 1986, pp. 142-144], in a few steps:
   ! 1) diagonalize S matrix:
   allocate(Ev(Nsiz))	! array with eigenvalues of nonorthogonal matrix
   Ev = 0.0d0	! to start from
   call sym_diagonalize(Sij, Ev, Error_descript) ! module "Algebra_tools"
   
   ! Now Sij is the collection of eigenvectors, Ev contains eigenvalues
   ! Moreover, Sij is now unitary matrix that can be used as U from Eq.(3.166) [Szabo "Modern Quantum Chemistry" 1986, p. 143]
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine Loewdin_Orthogonalization_c: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif
   ! Test if diagonalization of S matrix worked well:
   !$OMP WORKSHARE
   N_neg = COUNT(Ev < 0.0d0)
   N_zero = COUNT(ABS(Ev) < epsylon)
   !$OMP END WORKSHARE
   if ( (N_zero > 0) .or. (N_neg > 0) ) then
      print*, 'Subroutine Loewdin_Orthogonalization_c has zero or negative eigenvalues in the COMPLEX overlap matrix!'
      print*, 'Negative: ', N_neg, 'Zero: ', N_zero
!       do i1 = 1, size(Ev)
!          print*, 'NEG:', i1, Ev(i1)
!       enddo
      where(Ev(:)<0.0d0) Ev(:) = ABS(Ev(:))
   endif

   ! 2) construct S = [1/sqrt(Ev)]:
   allocate(s_mat(Nsiz,Nsiz)) ! temporary array with matrix of 1/sqrt(Ev)
   !$OMP WORKSHARE
   s_mat = 0.0d0 ! to start from
   ! Construct S matrix if eigenstates are non-zero:
   forall (i1=1:size(Ev), ABS(Ev(i1)) > epsylon) s_mat(i1,i1) = 1.0d0/dsqrt(Ev(i1))
   !$OMP END WORKSHARE

   ! 3) construct X = S*s_mat [Szabo "Modern Quantum Chemistry" 1986, p. 144, Eq.(3.169)]
   !    Xij = matmul(Sij, dcmplx(s_mat,0.0d0))  ! same result with standart routine
   call mkl_matrix_mult('N', 'N', Sij, cmplx(s_mat,0.0d0), Xij) ! module "Algebra_tools"
   
   ! 4) orthoginalize the Cij matrix [Szabo "Modern Quantum Chemistry" 1986, p. 144, Eq.(3.177)]:
   ! Note that the matrix is Hermitian, not Symmetric, thus the following subroutine cannot be used!:
!    call zsymm('L', 'U', Nsiz, Nsiz, dcmplx(1.0d0,0d0), Hij, Nsiz, Xij, Nsiz, dcmplx(0.0d0,0.0d0), Sij_temp, Nsiz)	! from BLAS library
   ! Use the standard fortran routine for matrices multiplication:
   Sij_temp = matmul(Hij, Xij)
   
!    Hij = matmul(dconjg(transpose(Xij)), Sij_temp) ! same result with standart routine:
   call mkl_matrix_mult('C', 'N', Xij, Sij_temp, Hij)	! module "Algebra_tools"
   
!     ! test (comment out for release):
!     call test_orthogonalization_c(Sijsave, Xij) ! Checked, correct!
   
    Sij = Xij	! save X matrix

    if (allocated(Xij)) deallocate(Xij)
    if (allocated(Sij_temp)) deallocate(Sij_temp)
!     if (allocated(Sijsave)) deallocate(Sijsave)
    if (allocated(s_mat)) deallocate(s_mat)
    if (allocated(Ev)) deallocate(Ev)
end subroutine Loewdin_Orthogonalization_c



subroutine test_orthogonalization_c(Sij_undiag, OBX) ! Check if it is orthogonalized correctly
   complex, dimension(:,:), intent(in) :: Sij_undiag	! undiagonalized overlap matrix
   complex, dimension(:,:), intent(in) :: OBX		! constructed X matrix
   complex, dimension(size(OBX,1),size(OBX,2)) :: Sij_temp, Sij_test	! test overlap of orthogonalized matrix
   integer :: N, i, j, i1, j1
   real(8) :: tolerance
   complex :: Sij_single
   tolerance = 1.0d-10		! acceptable inaccuracy in the basis orthogonalization
   N = size(OBX,1)
   Sij_temp = matmul(Sij_undiag, OBX)
   call mkl_matrix_mult('C', 'N',OBX, (Sij_temp), Sij_test) ! module "Algebra_tools"
   
!$omp PARALLEL
! !$omp do  private(i, j)
!$omp do  private(i, j, Sij_single)
   do i = 1, N
      do j = 1, N
         if (i == j) then ! all diagonal elements must be equal to 1
            if (  ABS(dble(Sij_test(i,j)) - 1.0d0) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Re ORTO:', i, j, Sij_test(i,j), Sij_undiag(i,j), OBX(i,j)
            if (  ABS(aimag(Sij_test(i,j))) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Im ORTO:', i, j, Sij_test(i,j), Sij_undiag(i,j), OBX(i,j)
         else ! all nondiagonal elements must be zeros
            if (  ABS(dble(Sij_test(i,j))) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Re ORTO:',  i, j, Sij_test(i,j), Sij_undiag(i,j), OBX(i,j)
            if (  ABS(aimag(Sij_test(i,j))) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Im ORTO:', i, j, Sij_test(i,j), Sij_undiag(i,j), OBX(i,j)
         endif
!--------------------------------------------
! DO NOT DELETE THE FOLLOWING LINES -- IT IS AN EXAMPLE OF HOW TO CORRECTLY WORK WITH INDICES OF THE MATRICES!
!           Sij_single = dcmplx(0.0d0,0.0d0)
!           do i1 = 1, N
!              do j1 = 1, N
!                 Sij_single = Sij_single + DCONJG(OBX(i1,i))*OBX(j1,j)*Sij_undiag(i1,j1) ! set symmetric elements
!              enddo
!           enddo
!            if (i == j) then ! all diagonal elements must be equal to 1
!               if (  ABS(dble(Sij_single) - 1.0d0) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e,e,e)') 'Re ORTO:', i, j, Sij_single, Sij_undiag(i,j), OBX(i,j), Sij_test(i,j)
!               if (  ABS(aimag(Sij_single)) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Im ORTO:', i, j, Sij_single, Sij_undiag(i,j), OBX(i,j)
!            else ! all nondiagonal elements must be zeros
!               if (  ABS(dble(Sij_single)) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Re ORTO:',  i, j, Sij_single, Sij_undiag(i,j), OBX(i,j)
!               if (  ABS(aimag(Sij_single)) > tolerance ) write(*,'(a,i5,i5,e,e,e,e,e,e)') 'Im ORTO:', i, j, Sij_single, Sij_undiag(i,j), OBX(i,j)
!            endif
!--------------------------------------------
      enddo
   enddo   
!$omp end do
!$omp END PARALLEL
! PAUSE 'test_orthogonalization DONE'
end subroutine test_orthogonalization_c



subroutine test_orthogonalization_r(Sij_undiag, Sij, OBX) ! Check if it is orthogonalized correctly
   real(8), dimension(:,:), intent(in) :: Sij_undiag	! undiagonalized overlap matrix
   real(8), dimension(:,:), intent(in) :: Sij		! eigenvectors of Sij_undiag, to use as unitary matrix
   real(8), dimension(:,:), intent(in) :: OBX		! constructed X matrix
   real(8), dimension(size(OBX,1),size(OBX,2)) :: Sij_temp, Sij_test	! test overlap of orthogonalized matrix
   integer :: N, i, j, i1, j1
   real(8) :: Sij_single, tolerance
   N = size(OBX,1)
   call dsymm('L', 'U', N, N, 1.0d0, Sij_undiag, N,  OBX, N, 0.0d0, Sij_temp, N) ! from BLAS library
   !call  dgemm ('T', 'N', N, N, N, 1.0d0, OBX, N, Sij_temp, N, 0.0d0, Sij_test, N) ! BLAS non-symmetric matrix multiplication
   call mkl_matrix_mult('T', 'N', OBX, Sij_temp, Sij_test) ! module "Algebra_tools"
   
   tolerance = 1.0d-10		! acceptable inaccuracy in the basis orthogonalization
   
!$omp PARALLEL
!$omp do  private(i, j)
   do i = 1, N
      do j = 1, N
         if (i == j) then ! all diagonal elements must be equal to 1
            if (  ABS(Sij_test(i,j) - 1.0d0) > tolerance ) write(*,'(a,i5,i5,f,f)') 'Failed test_orthogonalization_r:', i, j, Sij_test(i,j) !, Sij_single
         else ! all nondiagonal elements must be zeros
            if (  ABS(Sij_test(i,j)) > tolerance ) write(*,'(a,i5,i5,f,f)') 'Failed test_orthogonalization_r:', i, j, Sij_test(i,j) !,Sij_single
         endif
!--------------------------------------------
! DO NOT DELETE THE FOLLOWING LINES -- IT IS AN EXAMPLE OF HOW TO CORRECTLY WORK WITH INDICES OF THE MATRICES!
!           Sij_single = 0.0d0
!           do i1 = 1, N
!              do j1 = 1, N
!                 Sij_single = Sij_single + OBX(i1,i)*OBX(j1,j)*Sij_undiag(i1,j1) ! set symmetric elements
!              enddo
!           enddo
!--------------------------------------------
      enddo
   enddo   
!$omp end do
!$omp END PARALLEL
! PAUSE 'test_orthogonalization DONE'
end subroutine test_orthogonalization_r


subroutine Get_overlap_S_matrix(Scell, NSC, i, atom_2, j, TB, Sij, M_SVij, M_lmn)
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, INTENT(IN) :: i, j, atom_2
   type(TB_H_NRL),  intent(in), target :: TB	! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: Sij  ! hamiltonian
   real(8), dimension(:,:,:), intent(in), target :: M_SVij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   !---------------------------------------
   real(8) :: r
   real(8), dimension(10) :: vec_M_SVij
   integer :: k
   if (i /= j) then	! it's 2 different atoms:
      vec_M_SVij(:) = M_SVij(i,j,:)
      ! Construct the overlap matrix including angular part for sp3d5 basis set:
      call  KS_sp3d5(vec_M_SVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), Sij)	! module TB_Koster_Slater
   else	! it is the same atom
      Sij = 0.0d0
      forall (k=1:9) Sij(k,k)=1.0d0
   endif
end subroutine Get_overlap_S_matrix


subroutine Get_overlap_S_matrix_cell(Scell, NSC, i, j, zb, TB, Sij)
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, INTENT(IN) :: i, j
   integer, dimension(:), intent(in) :: zb	! vector to the mirror cell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: Sij  ! hamiltonian
   !---------------------------------------
   real(8) :: r, r1, x, y, z
   real(8), dimension(10) :: vec_M_SVij
   integer :: k
   logical :: origin_cell
   real(8), pointer :: rm
   integer, pointer :: KOA1, KOA2
   
   
   origin_cell = ALL(zb==0) ! if it is the origin cell
   if ((i /= j) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
!    if (i /= j) then	! it's 2 different atoms:
      call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i, j, r1, x=x, y=y, z=z) ! module "Atomic_tools"
      
      KOA1 => Scell(NSC)%MDatoms(j)%KOA
      KOA2 => Scell(NSC)%MDatoms(i)%KOA
      
      r = r1*g_A2au	! convert into [Bohr]
      rm => TB(KOA1,KOA2)%Rc ! cut-off radius
   
      ! All radial functions for Overlap matrix:
      vec_M_SVij(1) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(1), TB(KOA1,KOA2)%qllm(1), TB(KOA1,KOA2)%rllm(1), TB(KOA1,KOA2)%sllm(1), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(1) = (s s sigma)
      vec_M_SVij(2) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(2), TB(KOA1,KOA2)%qllm(2), TB(KOA1,KOA2)%rllm(2), TB(KOA1,KOA2)%sllm(2), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(2) = (s p sigma)
      vec_M_SVij(3) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(5), TB(KOA1,KOA2)%qllm(5), TB(KOA1,KOA2)%rllm(5), TB(KOA1,KOA2)%sllm(5), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(3) = (s d sigma)
      vec_M_SVij(4) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(3), TB(KOA1,KOA2)%qllm(3), TB(KOA1,KOA2)%rllm(3), TB(KOA1,KOA2)%sllm(3), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(4) = (p p sigma)
      vec_M_SVij(5) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(4), TB(KOA1,KOA2)%qllm(4), TB(KOA1,KOA2)%rllm(4), TB(KOA1,KOA2)%sllm(4), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(5) = (p p pi)
      vec_M_SVij(6) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(6), TB(KOA1,KOA2)%qllm(6), TB(KOA1,KOA2)%rllm(6), TB(KOA1,KOA2)%sllm(6), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(6) = (p d sigma)
      vec_M_SVij(7) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(7), TB(KOA1,KOA2)%qllm(7), TB(KOA1,KOA2)%rllm(7), TB(KOA1,KOA2)%sllm(7), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(7) = (p d pi)
      vec_M_SVij(8) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(8), TB(KOA1,KOA2)%qllm(8), TB(KOA1,KOA2)%rllm(8), TB(KOA1,KOA2)%sllm(8), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(8) = (d d sigma)
      vec_M_SVij(9) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(9), TB(KOA1,KOA2)%qllm(9), TB(KOA1,KOA2)%rllm(9), TB(KOA1,KOA2)%sllm(9), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(9) = (d d pi)
      vec_M_SVij(10) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(10), TB(KOA1,KOA2)%qllm(10), TB(KOA1,KOA2)%rllm(10), TB(KOA1,KOA2)%sllm(10), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(10) = (d d delta)

      ! Construct the overlap matrix including angular part for sp3d5 basis set:
      call  KS_sp3d5(vec_M_SVij, x/r1, y/r1, z/r1, Sij)	! module TB_Koster_Slater
   else	! it is the same atom
      Sij = 0.0d0
      forall (k=1:9) Sij(k,k)=1.0d0
   endif
   nullify (KOA1, KOA2, rm)
end subroutine Get_overlap_S_matrix_cell




subroutine Hamilton_one_NRL(Scell, NSC, i, j, atom_2, TB, Hij, M_Vij, M_lmn)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, intent(in) :: i, j, atom_2
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   real(8), dimension(:,:), intent(out) :: Hij  ! hamiltonian, all orbitals in sp3d5 basis set
   real(8), dimension(:,:,:), intent(in) :: M_Vij	! matrix of Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines
   !---------------------------------------/
   integer(4) ki, kj, i1, j1, k1, ik, n_orb
   real(8) x,y,z,r,r1, x0, y0, z0, sx, sy, sz !  interatomic distance projections, and the total distance
   integer :: cell_x1, cell_y1, cell_z1 ! cell numbers
   real(8), dimension(3) :: zb
   real(8) :: rho	! an embedded-atom-like ‘density’
   real(8), dimension(:), pointer :: al, bl, cl, dl	! set of parameters in Eq.(17)
   real(8), pointer :: rm
   integer, pointer :: KOA1, KOA2

   n_orb = 9	! number of arbitals for sp3d5 basis set

   KOA1 => Scell(NSC)%MDatoms(j)%KOA
   KOA2 => Scell(NSC)%MDatoms(i)%KOA

   if (i == j) then ! Onsite contributions
      Hij = 0.0d0   ! Nondiagonals are zeros

      al => TB(KOA1,KOA2)%al
      bl => TB(KOA1,KOA2)%bl
      cl => TB(KOA1,KOA2)%cl
      dl => TB(KOA1,KOA2)%dl
      
      ! s orbital of sp3d5 basis set:
      call rho_NRL(Scell, NSC, TB, j, rho)  ! to get rho
      Hij(1,1) = On_site_NRL(rho, al(1), bl(1), cl(1), dl(1)) 	! function below
      ! p3 orbitals of sp3d5 basis set:
      do i1 = 2, 4
         Hij(i1,i1) = On_site_NRL(rho, al(2), bl(2), cl(2), dl(2)) 	! function below
      enddo
      ! d5 orbitals of sp3d5 basis set:
      do i1 = 5, 7
         Hij(i1,i1) = On_site_NRL(rho, al(3), bl(3), cl(3), dl(3)) 	! function below
      enddo
      do i1 = 8, 9
         if (TB(KOA1,KOA2)%ind_split == 1) then	! splitting between t2g and e2
            Hij(i1,i1) = On_site_NRL(rho, al(4), bl(4), cl(4), dl(4)) 	! function below
         else	! no splitting between t2g and e2
            Hij(i1,i1) = On_site_NRL(rho, al(3), bl(3), cl(3), dl(3)) 	! function below
         endif
      enddo
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:
      r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
      r = r*g_A2au	! convert into [Bohr]
      rm => TB(KOA1,KOA2)%Rc ! cut-off radius
!       Hij = 0.0d0   ! Nondiagonals are zeros
      call Hopping_NRL(M_Vij(i,j,:), M_lmn(:,i,j), r, rm, Hij)	! subroutine below
   endif
   
   nullify (al, bl, cl, dl, KOA1, KOA2, rm)
end subroutine Hamilton_one_NRL



subroutine Hamilton_one_NRL_cell(Scell, NSC, i, j, zb, TB, numpar, Hij) ! this calles the block-hamiltonian
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, intent(in) :: i, j
   integer, dimension(:), intent(in) :: zb	! vector to the mirror cell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   real(8), dimension(:,:), intent(inout) :: Hij  ! hamiltonian, all orbitals in sp3d5 basis set
   !---------------------------------------/
   integer(4) ki, kj, i1, j1, k1, ik, n_orb
   real(8) x,y,z,r,r1, x0, y0, z0, sx, sy, sz !  interatomic distance projections, and the total distance
   integer :: cell_x1, cell_y1, cell_z1 ! cell numbers
   real(8) :: rho	! an embedded-atom-like ‘density’
   real(8), dimension(:), pointer :: al, bl, cl, dl	! set of parameters in Eq.(17)
   real(8), dimension(10) :: Vsr
   real(8), pointer :: rm
   integer, pointer :: KOA1, KOA2
   logical :: origin_cell
   
   n_orb = 9	! number of arbitals for sp3d5 basis set
   
   KOA1 => Scell(NSC)%MDatoms(j)%KOA
   KOA2 => Scell(NSC)%MDatoms(i)%KOA

!    if (i == j) then ! Onsite contributions   
   origin_cell = ALL(zb==0) ! if it is the origin cell
   if ((i == j) .and. (origin_cell)) then ! exclude self-interaction only within original super cell
      Hij = 0.0d0   ! Nondiagonals are zeros

      al => TB(KOA1,KOA2)%al
      bl => TB(KOA1,KOA2)%bl
      cl => TB(KOA1,KOA2)%cl
      dl => TB(KOA1,KOA2)%dl
      
      ! s orbital of sp3d5 basis set:
      call rho_NRL_cell(Scell, NSC, TB, numpar, j, rho)	! to get rho      
      
      Hij(1,1) = On_site_NRL(rho, al(1), bl(1), cl(1), dl(1)) 	! function below
      ! p3 orbitals of sp3d5 basis set:
      do i1 = 2, 4
         Hij(i1,i1) = On_site_NRL(rho, al(2), bl(2), cl(2), dl(2)) 	! function below
      enddo
      ! d5 orbitals of sp3d5 basis set:
      do i1 = 5, 7
         Hij(i1,i1) = On_site_NRL(rho, al(3), bl(3), cl(3), dl(3)) 	! function below
      enddo
      do i1 = 8, 9
         if (TB(KOA1,KOA2)%ind_split == 1) then	! splitting between t2g and e2
            Hij(i1,i1) = On_site_NRL(rho, al(4), bl(4), cl(4), dl(4)) 	! function below
         else	! no splitting between t2g and e2
            Hij(i1,i1) = On_site_NRL(rho, al(3), bl(3), cl(3), dl(3)) 	! function below
         endif
      enddo
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:
      call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i, j, r1, x=x, y=y, z=z) ! module "Atomic_tools"

      r = r1*g_A2au	! convert into [Bohr]
      rm => TB(KOA1,KOA2)%Rc ! cut-off radius
      
!       Hij = 0.0d0   ! Nondiagonals are zeros
      
      ! All radial functions for Hamiltonian:
      Vsr(1) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(1), TB(KOA1,KOA2)%fllm(1), TB(KOA1,KOA2)%gllm(1), TB(KOA1,KOA2)%hllm(1), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(1) = (s s sigma)
      Vsr(2) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(2), TB(KOA1,KOA2)%fllm(2), TB(KOA1,KOA2)%gllm(2), TB(KOA1,KOA2)%hllm(2), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(2) = (s p sigma)
      Vsr(3) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(5), TB(KOA1,KOA2)%fllm(5), TB(KOA1,KOA2)%gllm(5), TB(KOA1,KOA2)%hllm(5), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(3) = (s d sigma)
      Vsr(4) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(3), TB(KOA1,KOA2)%fllm(3), TB(KOA1,KOA2)%gllm(3), TB(KOA1,KOA2)%hllm(3), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(4) = (p p sigma)
      Vsr(5) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(4), TB(KOA1,KOA2)%fllm(4), TB(KOA1,KOA2)%gllm(4), TB(KOA1,KOA2)%hllm(4), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(5) = (p p pi)
      Vsr(6) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(6), TB(KOA1,KOA2)%fllm(6), TB(KOA1,KOA2)%gllm(6), TB(KOA1,KOA2)%hllm(6), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(6) = (p d sigma)
      Vsr(7) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(7), TB(KOA1,KOA2)%fllm(7), TB(KOA1,KOA2)%gllm(7), TB(KOA1,KOA2)%hllm(7), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(7) = (p d pi)
      Vsr(8) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(8), TB(KOA1,KOA2)%fllm(8), TB(KOA1,KOA2)%gllm(8), TB(KOA1,KOA2)%hllm(8), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(8) = (d d sigma)
      Vsr(9) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(9), TB(KOA1,KOA2)%fllm(9), TB(KOA1,KOA2)%gllm(9), TB(KOA1,KOA2)%hllm(9), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(9) = (d d pi)
      Vsr(10) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(10), TB(KOA1,KOA2)%fllm(10), TB(KOA1,KOA2)%gllm(10), TB(KOA1,KOA2)%hllm(10), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(10) = (d d delta)

!       call Hopping_NRL(M_Vij(i,j,:), M_lmn(:,i,j), r, rm, Hij)	! subroutine below
      call KS_sp3d5(Vsr, x/r1, y/r1, z/r1, Hij)		! module TB_Koster_Slater
   endif
   
   nullify (al, bl, cl, dl, KOA1, KOA2, rm)
end subroutine Hamilton_one_NRL_cell



subroutine Hopping_NRL(M_Vij, M_lmn, r, rm, ts)
   real(8), dimension(:), intent(in) :: M_Vij	! Overlap functions for all pairs of atoms, all orbitals
   real(8), dimension(:), intent(in) :: M_lmn	! cosine directions
   real(8), intent(in) :: r		! relative distances (projections) between the atoms
   real(8), intent(in):: rm 	! cut-off radius
   real(8), dimension(9,9), intent(out) :: ts	! overlap integerals [a.u.]
   !=============================================
   real(8), dimension(10) :: vec_M_Vij
   ! Construct the overlap integrals including angular part for sp3d5 basis set:
   vec_M_Vij = M_Vij
   call  KS_sp3d5(vec_M_Vij, M_lmn(1), M_lmn(2), M_lmn(3), ts)	! module TB_Koster_Slater
end subroutine Hopping_NRL


subroutine Construct_Vij_NRL(numpar, TB, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)	! module "TB_NRL"
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable :: M_Vij	! matrix of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dVij	! matrix of derivatives of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_SVij	! matrix of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dSVij	! matrix of derivatives of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   !----------------------------------
   real(8) :: x, y, z, r, sx, sy, sz, r1
   integer :: i, j, atom_2, ki
   real(8), pointer :: rm
   integer, pointer :: nat, m, KOA1, KOA2
   
   nat => Scell(NSC)%Na	! number of atoms in the supercell
   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,10))	! each pair of atoms, all 10 V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,10))	! each pair of atoms, all 10 V functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,10))	! each pair of atoms, all 10 V functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,10))	! each pair of atoms, all 10 V functions
   
   !$OMP WORKSHARE
   M_Vij = 0.0d0
   M_dVij = 0.0d0
   M_SVij = 0.0d0
   M_dSVij = 0.0d0
   !$OMP END WORKSHARE
   
   ! Construct matrix of all the radial functions for each pair of atoms, in 2 steps:
   
! 1) Upper triangle of the matrix (i,j) - calculate every element:
!$omp PARALLEL
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, rm)
   AT1:do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one  
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         
         IJ:if (i > j) then ! it's a new pair of atoms, calculate everything
            KOA1 => Scell(NSC)%MDatoms(j)%KOA
            KOA2 => Scell(NSC)%MDatoms(i)%KOA
            r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R
            ! In the formula, r must be in Bohr, so convert into it from Angstroms:
            r = r*g_A2au

            rm => TB(KOA1,KOA2)%Rc ! cut-off radius [Bohr]

            ! All radial functions for Hamiltonian:
            M_Vij(j,i,1) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(1), TB(KOA1,KOA2)%fllm(1), TB(KOA1,KOA2)%gllm(1), TB(KOA1,KOA2)%hllm(1), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(1) = (s s sigma)
            M_Vij(j,i,2) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(2), TB(KOA1,KOA2)%fllm(2), TB(KOA1,KOA2)%gllm(2), TB(KOA1,KOA2)%hllm(2), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(2) = (s p sigma)
            M_Vij(j,i,3) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(5), TB(KOA1,KOA2)%fllm(5), TB(KOA1,KOA2)%gllm(5), TB(KOA1,KOA2)%hllm(5), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(3) = (s d sigma)
            M_Vij(j,i,4) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(3), TB(KOA1,KOA2)%fllm(3), TB(KOA1,KOA2)%gllm(3), TB(KOA1,KOA2)%hllm(3), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(4) = (p p sigma)
            M_Vij(j,i,5) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(4), TB(KOA1,KOA2)%fllm(4), TB(KOA1,KOA2)%gllm(4), TB(KOA1,KOA2)%hllm(4), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(5) = (p p pi)
            M_Vij(j,i,6) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(6), TB(KOA1,KOA2)%fllm(6), TB(KOA1,KOA2)%gllm(6), TB(KOA1,KOA2)%hllm(6), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(6) = (p d sigma)
            M_Vij(j,i,7) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(7), TB(KOA1,KOA2)%fllm(7), TB(KOA1,KOA2)%gllm(7), TB(KOA1,KOA2)%hllm(7), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(7) = (p d pi)
            M_Vij(j,i,8) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(8), TB(KOA1,KOA2)%fllm(8), TB(KOA1,KOA2)%gllm(8), TB(KOA1,KOA2)%hllm(8), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(8) = (d d sigma)
            M_Vij(j,i,9) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(9), TB(KOA1,KOA2)%fllm(9), TB(KOA1,KOA2)%gllm(9), TB(KOA1,KOA2)%hllm(9), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(9) = (d d pi)
            M_Vij(j,i,10) = Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(10), TB(KOA1,KOA2)%fllm(10), TB(KOA1,KOA2)%gllm(10), TB(KOA1,KOA2)%hllm(10), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(10) = (d d delta)

!             M_Vij(:,:,1:10) = 0.0d0
            
            ! All derivatives of the radial functions:
            M_dVij(j,i,1) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(1), TB(KOA1,KOA2)%fllm(1), TB(KOA1,KOA2)%gllm(1), TB(KOA1,KOA2)%hllm(1), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(1) = (s s sigma)
            M_dVij(j,i,2) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(2), TB(KOA1,KOA2)%fllm(2), TB(KOA1,KOA2)%gllm(2), TB(KOA1,KOA2)%hllm(2), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(2) = (s p sigma)
            M_dVij(j,i,3) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(5), TB(KOA1,KOA2)%fllm(5), TB(KOA1,KOA2)%gllm(5), TB(KOA1,KOA2)%hllm(5), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(3) = (s d sigma)
            M_dVij(j,i,4) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(3), TB(KOA1,KOA2)%fllm(3), TB(KOA1,KOA2)%gllm(3), TB(KOA1,KOA2)%hllm(3), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(4) = (p p sigma)
            M_dVij(j,i,5) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(4), TB(KOA1,KOA2)%fllm(4), TB(KOA1,KOA2)%gllm(4), TB(KOA1,KOA2)%hllm(4), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(5) = (p p pi)
            M_dVij(j,i,6) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(6), TB(KOA1,KOA2)%fllm(6), TB(KOA1,KOA2)%gllm(6), TB(KOA1,KOA2)%hllm(6), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(6) = (p d sigma)
            M_dVij(j,i,7) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(7), TB(KOA1,KOA2)%fllm(7), TB(KOA1,KOA2)%gllm(7), TB(KOA1,KOA2)%hllm(7), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(7) = (p d pi)
            M_dVij(j,i,8) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(8), TB(KOA1,KOA2)%fllm(8), TB(KOA1,KOA2)%gllm(8), TB(KOA1,KOA2)%hllm(8), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(8) = (d d sigma)
            M_dVij(j,i,9) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(9), TB(KOA1,KOA2)%fllm(9), TB(KOA1,KOA2)%gllm(9), TB(KOA1,KOA2)%hllm(9), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(9) = (d d pi)
            M_dVij(j,i,10) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%ellm(10), TB(KOA1,KOA2)%fllm(10), TB(KOA1,KOA2)%gllm(10), TB(KOA1,KOA2)%hllm(10), rm, TB(KOA1,KOA2)%lden, 0, 0.0d0)	!    V(10) = (d d delta)

!             M_dVij(:,:,1:10) = 0.0d0
            
            ! All radial functions for Overlap matrix:
            M_SVij(j,i,1) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(1), TB(KOA1,KOA2)%qllm(1), TB(KOA1,KOA2)%rllm(1), TB(KOA1,KOA2)%sllm(1), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(1) = (s s sigma)
            M_SVij(j,i,2) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(2), TB(KOA1,KOA2)%qllm(2), TB(KOA1,KOA2)%rllm(2), TB(KOA1,KOA2)%sllm(2), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(2) = (s p sigma)
            M_SVij(j,i,3) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(5), TB(KOA1,KOA2)%qllm(5), TB(KOA1,KOA2)%rllm(5), TB(KOA1,KOA2)%sllm(5), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(3) = (s d sigma)
            M_SVij(j,i,4) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(3), TB(KOA1,KOA2)%qllm(3), TB(KOA1,KOA2)%rllm(3), TB(KOA1,KOA2)%sllm(3), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(4) = (p p sigma)
            M_SVij(j,i,5) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(4), TB(KOA1,KOA2)%qllm(4), TB(KOA1,KOA2)%rllm(4), TB(KOA1,KOA2)%sllm(4), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(5) = (p p pi)
            M_SVij(j,i,6) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(6), TB(KOA1,KOA2)%qllm(6), TB(KOA1,KOA2)%rllm(6), TB(KOA1,KOA2)%sllm(6), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(6) = (p d sigma)
            M_SVij(j,i,7) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(7), TB(KOA1,KOA2)%qllm(7), TB(KOA1,KOA2)%rllm(7), TB(KOA1,KOA2)%sllm(7), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(7) = (p d pi)
            M_SVij(j,i,8) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(8), TB(KOA1,KOA2)%qllm(8), TB(KOA1,KOA2)%rllm(8), TB(KOA1,KOA2)%sllm(8), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(8) = (d d sigma)
            M_SVij(j,i,9) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(9), TB(KOA1,KOA2)%qllm(9), TB(KOA1,KOA2)%rllm(9), TB(KOA1,KOA2)%sllm(9), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(9) = (d d pi)
            M_SVij(j,i,10) = Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(10), TB(KOA1,KOA2)%qllm(10), TB(KOA1,KOA2)%rllm(10), TB(KOA1,KOA2)%sllm(10), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(10) = (d d delta)

!             M_SVij(:,:,1:10) = 0.0d0
            
            ! All derivatives of the radial functions for Overlap matrix:
            M_dSVij(j,i,1) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(1), TB(KOA1,KOA2)%qllm(1), TB(KOA1,KOA2)%rllm(1), TB(KOA1,KOA2)%sllm(1), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(1) = (s s sigma)
            M_dSVij(j,i,2) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(2), TB(KOA1,KOA2)%qllm(2), TB(KOA1,KOA2)%rllm(2), TB(KOA1,KOA2)%sllm(2), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(2) = (s p sigma)
            M_dSVij(j,i,3) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(5), TB(KOA1,KOA2)%qllm(5), TB(KOA1,KOA2)%rllm(5), TB(KOA1,KOA2)%sllm(5), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(3) = (s d sigma)
            M_dSVij(j,i,4) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(3), TB(KOA1,KOA2)%qllm(3), TB(KOA1,KOA2)%rllm(3), TB(KOA1,KOA2)%sllm(3), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(4) = (p p sigma)
            M_dSVij(j,i,5) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(4), TB(KOA1,KOA2)%qllm(4), TB(KOA1,KOA2)%rllm(4), TB(KOA1,KOA2)%sllm(4), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(5) = (p p pi)
            M_dSVij(j,i,6) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(6), TB(KOA1,KOA2)%qllm(6), TB(KOA1,KOA2)%rllm(6), TB(KOA1,KOA2)%sllm(6), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(6) = (p d sigma)
            M_dSVij(j,i,7) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(7), TB(KOA1,KOA2)%qllm(7), TB(KOA1,KOA2)%rllm(7), TB(KOA1,KOA2)%sllm(7), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 0.0d0)	!    V(7) = (p d pi)
            M_dSVij(j,i,8) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(8), TB(KOA1,KOA2)%qllm(8), TB(KOA1,KOA2)%rllm(8), TB(KOA1,KOA2)%sllm(8), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(8) = (d d sigma)
            M_dSVij(j,i,9) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(9), TB(KOA1,KOA2)%qllm(9), TB(KOA1,KOA2)%rllm(9), TB(KOA1,KOA2)%sllm(9), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(9) = (d d pi)
            M_dSVij(j,i,10) = d_Overlap_function_NRL(r, TB(KOA1,KOA2)%pllm(10), TB(KOA1,KOA2)%qllm(10), TB(KOA1,KOA2)%rllm(10), TB(KOA1,KOA2)%sllm(10), rm, TB(KOA1,KOA2)%lden, TB(KOA1,KOA2)%ind_overlap, 1.0d0)	!    V(10) = (d d delta)
            
!             M_dSVij(:,:,1:10) = 0.0d0
            
         endif IJ
      enddo AT2
   enddo AT1
!$omp end do
!$omp END PARALLEL

! 2) Lower triangle of the matrix (i,j) - use symmetry:
!$omp PARALLEL
!$omp do private(j, m, atom_2, i)
   do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 1,m ! do only for atoms close to that one  
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         if (i < j) then ! this pair of atoms was already calculated, use the symmetry of the Hamiltonial and S matrix:
            M_Vij(j,i,:) = M_Vij(i,j,:)
            M_dVij(j,i,:) = M_dVij(i,j,:)
            M_SVij(j,i,:) = M_SVij(i,j,:)
            M_dSVij(j,i,:) = M_dSVij(i,j,:)
         endif
      enddo
   enddo
!$omp end do
!$omp END PARALLEL

   nullify(rm, nat, m, KOA1, KOA2)	! clean up at the end
end subroutine Construct_Vij_NRL


! Density of on-site density function:
subroutine rho_NRL(Scell, NSC, TB_Hamil, j, rho)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB_Hamil   ! parameters of the Hamiltonian of TB
   integer, intent(in) :: j
   REAL(8), INTENT(out) :: rho	 ! Hamiltonian
   !-----------------
   integer :: nat
   integer :: atom_2, m
   real(8) :: r
   real(8), pointer :: r1
   real(8), pointer :: lambd, rc, lden
   integer, pointer :: KOA1, KOA2, i

   nat = Scell(NSC)%Na

   rho = 0.0d0
   m = Scell(NSC)%Near_neighbor_size(j)
!    !$omp parallel private(atom_2,i,KOA1,KOA2,lambd,rc,lden,r1,r)
!    !$omp do reduction( + : rho)
   do atom_2 = 1,m ! do only for atoms close to that one  
      i => Scell(NSC)%Near_neighbor_list(j,atom_2)		! this is the list of such close atoms
      KOA1 => Scell(NSC)%MDatoms(j)%KOA
      KOA2 => Scell(NSC)%MDatoms(i)%KOA
      lambd => TB_Hamil(KOA1,KOA2)%lambd
      rc => TB_Hamil(KOA1,KOA2)%Rc
      lden => TB_Hamil(KOA1,KOA2)%lden
!       call shortest_distance(Scell, NSC, Scell(NSC)%MDAtoms, i, j, r)	! module "Atomic_tools"
      r1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R
      
      ! In this formula, r must be in Bohr, so convert into it from Angstroms:
      r = r1*g_A2au

      rho = rho + exp_F_NRL(r, lambd, rc, lden)	! function below
   enddo ! atom_2
!    !$omp end do
!    !$omp end parallel
   
   nullify(lambd, rc, lden, i, KOA1, KOA2, r1)
end subroutine rho_NRL



subroutine rho_NRL_cell(Scell, NSC, TB_Hamil, numpar, j, rho)	! to get rho
 type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   integer, intent(in) :: j
   REAL(8), INTENT(out) :: rho	 ! Hamiltonian
   !-----------------
   integer :: nat, Nx, Ny, Nz
   integer :: i, m, x_cell, y_cell, z_cell
   real(8) :: r, r1
   integer, dimension(3) :: zb	! vector to the mirror-cell
   real(8), pointer :: lambd, rc, lden
   integer, pointer :: KOA1, KOA2
   logical :: origin_cell
   
   ! 0) Get the number of cells to be included:
   call get_mirror_cell_num_NRL(Scell, NSC, numpar, Nx, Ny, Nz)	! see below

   nat = Scell(NSC)%Na
   rho = 0.0d0
   do i = 1,nat
      XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
         YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
            ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
               zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell_x
               
               origin_cell = ALL(zb==0) ! if it is the origin cell
               if ((i /= j) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                  KOA1 => Scell(NSC)%MDatoms(j)%KOA
                  KOA2 => Scell(NSC)%MDatoms(i)%KOA
                  lambd => TB_Hamil(KOA1,KOA2)%lambd
                  rc => TB_Hamil(KOA1,KOA2)%Rc
                  lden => TB_Hamil(KOA1,KOA2)%lden

                  call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), j, i, r1) ! module "Atomic_tools"
                  r = r1*g_A2au	! convert into [Bohr]

                  rho = rho + exp_F_NRL(r, lambd, rc, lden)	! function below
               endif
               
            enddo ZC
         enddo YC
      enddo XC
   enddo

   nullify(lambd, rc, lden, KOA1, KOA2)
end subroutine rho_NRL_cell






pure function On_site_NRL(rho, al, bl, cl, dl) result(hil)
   real(8) :: hil
   real(8), intent(in) :: rho, al, bl, cl, dl
   real(8) :: rho_23
   rho_23 = rho**m_two_third
   hil = al + bl*rho_23 + cl*rho_23*rho_23 + dl*rho*rho	! Eq.(17) from [1]
end function On_site_NRL


pure function Overlap_function_NRL(r, pllm, qllm, rllm, lambd, rc, dl, ind, ll) result(Hllm)	! Eq.(19), (20) or (23) from [1]
   real(8) :: Hllm
   real(8), intent(in) :: r, pllm, qllm, rllm, lambd, rc, dl
   integer, intent(in) :: ind	! 0=Eq.(19) or (20);  1=Eq.(23)
   real(8), intent(in) :: ll	! delta (l l')
   !-------------------------------
   real(8) :: r2
   select case (ind)
   case (1)
      r2 = r*r
      Hllm = (ll + pllm*r + qllm*r2 + rllm*r*r2)*exp_F_NRL(r, lambd, rc, dl)	! function below
   case default
      Hllm = (pllm + qllm*r + rllm*r*r)*exp_F_NRL(r, lambd, rc, dl)	! function below
   endselect
end function Overlap_function_NRL


pure function exp_F_NRL(r, lambd, rc, dl) result(exp_F)	! part of Eq.(15), (19), (20), and (23) from [1]
   real(8) :: exp_F
   real(8), intent(in) :: r, rc, dl, lambd
   if (r >= rc) then
      exp_F = 0.0d0
   else
      exp_F = dexp(-lambd*lambd*r)*F_NRL(r, rc, dl)	! function below
   endif
end function exp_F_NRL


pure function F_NRL(r, rc, dl) result(F) ! Eq.(16) from [1]
   real(8) :: F
   real(8), intent(in) :: r, rc, dl
   if (r >= rc) then
      F = 0.0d0
   else
      F = 1.0d0/(1.0d0 + exp( (r - rc)/dl + 5.0d0 ))
!      F = 1.0d0/(1.0d0 + dexp( (r - rc + 5.0d0)/dl ))
   endif
end function F_NRL

pure function d_F_NRL(r, rc, dl) result(F) ! Eq.(16) from [1], dF/drij
   real(8) :: F
   real(8), intent(in) :: r, rc, dl
   real(8) :: exp_r, denom
   if (r >= rc) then
      F = 0.0d0
   else
      exp_r = exp( (r - rc)/dl + 5.0d0 )
!       exp_r = dexp( (r - rc + 5.0d0)/dl )
      denom = 1.0d0 + exp_r
      F = -exp_r/(dl*denom*denom)
   endif
end function d_F_NRL


!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Derivatives of the attractive part of TB:

! Subroutine for derivative of the Hamiltonian:
subroutine get_dHij_drij_NRL(TB_Hamil, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei) ! module "TB_NRL"
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   !------------------------------------------------------------
   integer :: nat, k
   !------------------------------------------------------------

   nat = size(Scell(NSC)%MDatoms)	! number of atoms
   
   !$omp PARALLEL private(k) 
   !$omp do
   ATOMS:do k = 1, nat	! forces for all atoms
      Scell(NSC)%MDatoms(k)%forces%att(:) = 0.0d0	! just to start
      call get_forces(k, TB_Hamil, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)	! see below
   enddo ATOMS
   !$omp end do 
   !$omp end parallel
end subroutine get_dHij_drij_NRL



subroutine get_forces(k, TB, Scell, NSC, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)
   integer, intent(in) :: k	! forces for this atom
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB	! parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of the supercell
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   !-----------------------------
   integer :: j, m, atom_2, i, j1, i1, nat, Nsiz, n_orb, i4, j4, j4j1, i4i1, INFO
   real(8), dimension(:,:,:), allocatable :: dH, dS
   real(8), dimension(3,9,9) :: dH1, dS1
   
   nat = size(Scell(NSC)%MDatoms)	! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1)	! total number of orbitals
   n_orb = 9		! number of orbitals per atom for sp3d5 basis set

!    print*, 'get_forces: 0'
   
   if (.not.allocated(dH)) allocate(dH(3,Nsiz,Nsiz), source=0.0d0, stat=INFO)
   if (INFO /= 0) print*, 'Problem allocating dH in get_forces, module TB_NRL'
   if (.not.allocated(dS)) allocate(dS(3,Nsiz,Nsiz), source=0.0d0, stat=INFO)
   if (INFO /= 0) print*, 'Problem allocating dS in get_forces, module TB_NRL'
!    dH = 0.0d0
!    dS = 0.0d0
   
!    print*, 'get_forces: 1'
   
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
            call d_Hamilton_one_NRL(k, Scell, NSC, i, j, atom_2, TB, dH1, M_Vij, M_dVij, M_lmn, dS1, M_SVij, M_dSVij) ! this calls the block-hamiltonian
            do j1 = 1,n_orb	! all orbitals
               j4j1 = j4+j1
               do i1 = 1,n_orb	! all orbitals
                  i4i1 = i4+i1
                  dH(:,i4i1,j4j1) = dH1(:,i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  dS(:,i4i1,j4j1) = dS1(:,i1,j1)	! construct the total Overlap Matrix from the blocks of one-atom overlap matrices
               enddo ! i1
            enddo ! j1
         endif IJ
      enddo ATOM2
   enddo ATOM1
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

   ! 2) Convert the derivative of the Hamiltonian into eVs from Ry:
   dH = dH*g_Ry	! [Ry] into [eV]

!    print*, 'get_forces: 2'

   ! 3) Calculate the forces form the derivatives and the eigervectors:
   call Attract_TB_forces_NRL(Aij, Aij_x_Ei, dH, dS, Scell, NSC, Scell(NSC)%MDatoms(k)%forces%att(:))

!    print*, 'get_forces: 3'
   
   if (allocated(dH)) deallocate(dH, STAT=INFO)
   if (INFO /= 0) print*, 'Problem deallocating dH in get_forces from TB_NRL'
   
!    print*, 'get_forces: 4'
   
   if (allocated(dS)) deallocate(dS, STAT=INFO)
   if (INFO /= 0) print*, 'Problem deallocating dS in get_forces from TB_NRL'
   
!    print*, 'get_forces: 5'
end subroutine get_forces



subroutine Attract_TB_forces_NRL(Aij, Aij_x_Ei, dH, dS, Scell, NSC, Eelectr_s)
   real(8), dimension(:,:,:), intent(in) :: dH, dS
   real(8), dimension(:,:), intent(in), target :: Aij, Aij_x_Ei
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(out)  :: Eelectr_s ! part of the forces

   integer :: j, k, i, ste, n, n_orb, norb_1, j_norb, n_too
   integer, target :: i2
   integer, pointer :: m, j1

   n = size(Aij,1)	! total number of orbitals
   n_orb = 9	! sp3d5 basis set
   norb_1 = n_orb - 1
   Eelectr_s = 0.0d0
   n_too = 0   ! to start with

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

             if (maxval(ABS(Eelectr_s(:))) .GE. 1.0d7) n_too = n_too + 1
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
      write(*,'(a)') 'Trouble in subroutine Attract_TB_forces_NRL, too large attractive force'
      write(*,'(a,i)') 'For elements: ', n_too
   endif

   nullify(m, j1)
end subroutine Attract_TB_forces_NRL

!ddddddddddddddddddddddddddddddddddddddddddddddddddd
! Derivatives:
subroutine d_Hamilton_one_NRL(k, Scell, NSC, i, j, atom_2, TB, dH, M_Vij, M_dVij, M_lmn, dS, M_SVij, M_dSVij)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC, i, j, atom_2, k
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   real(8), dimension(:,:,:), intent(out) :: dH, dS  ! hamiltonian, all orbitals in sp3d5 basis set
   real(8), dimension(:,:,:), intent(in), target :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines and derivatives
   real(8), dimension(:,:,:), intent(in), target :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms for S-matrix, all orbitals, and derivatives
   !---------------------------------------
   integer :: ki, i1, n_orb
   real(8) :: drij_dsk(3)
   real(8), dimension(9) :: M_dlmn	! dl/dx, dl/dy, dl/dz, dm/dx, dm/dy, dm/dz, dn/dx, dn/dy, dn/dz
   real(8) :: rho, d_rho(3)	! an embedded-atom-like ‘density’ and its derivative
   real(8) :: dH_on
   real(8), pointer :: rm, x1, y1, z1, r1
   real(8), dimension(:), pointer :: al, bl, cl, dl
   integer, pointer :: KOA1, KOA2
   real(8), dimension(10) :: vec_M_Vij, vec_M_SVij
   real(8), dimension(10) :: vec_M_dVij, vec_M_dSVij
   real(8), dimension(9,9) :: dH1, dS1

   if (i == j) then ! Onsite contributions
      dH(:,:,:)  = 0.0d0
      KOA1 => Scell(NSC)%MDatoms(j)%KOA
      KOA2 => Scell(NSC)%MDatoms(i)%KOA
      al => TB(KOA1,KOA2)%al
      bl => TB(KOA1,KOA2)%bl
      cl => TB(KOA1,KOA2)%cl
      dl => TB(KOA1,KOA2)%dl
      ! s orbital of sp3d5 basis set:
      call rho_NRL(Scell, NSC, TB, j, rho)	! to get rho
      ! Get the derivatives of rho:
      call d_On_site_NRL(Scell, NSC, TB, rho, j, k, al(1), bl(1), cl(1), dl(1), d_rho)	! function below
      dH(:,1,1) = d_rho(:)
      ! p3 orbitals of sp3d5 basis set:
      do i1 = 2, 4
         call d_On_site_NRL(Scell, NSC, TB, rho, j, k, al(2), bl(2), cl(2), dl(2), d_rho)	! function below
         dH(:,i1,i1) = d_rho(:)
      enddo
      ! d5 orbitals of sp3d5 basis set:
      do i1 = 5, 7
         call d_On_site_NRL(Scell, NSC, TB, rho, j, k, al(3), bl(3), cl(3), dl(3), d_rho)	! function below
         dH(:,i1,i1) = d_rho(:)
      enddo
      do i1 = 8, 9
         if (TB(KOA1,KOA2)%ind_split == 1) then	! splitting between t2g and e2
            call d_On_site_NRL(Scell, NSC, TB, rho, j, k, al(4), bl(4), cl(4), dl(4), d_rho)	! function below
         else	! no splitting between t2g and e2
            call d_On_site_NRL(Scell, NSC, TB, rho, j, k, al(3), bl(3), cl(3), dl(3), d_rho)	! function below
         endif
         dH(:,i1,i1) = d_rho(:)
      enddo
      dS = 0.0d0	! const gives 0 always
   else	! For pairs of atoms, fill the hamiltonain with Hopping Integrals:   
      KOA1 => Scell(NSC)%MDatoms(i)%KOA
      KOA2 => Scell(NSC)%MDatoms(j)%KOA
      rm => TB(KOA1,KOA2)%Rc ! cut-off radius
      
      x1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,1)	! at this distance, X
      y1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,2)	! at this distance, Y
      z1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,3)	! at this distance, Z
      r1 => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)		! at this distance, R
      
      vec_M_Vij = M_Vij(i,j,:)
      vec_M_SVij = M_SVij(i,j,:)

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
         
      ! Convert [A] -> [Bohr]:
      drij_dsk(:) = drij_dsk(:)*g_A2au
         
      ! Derivative along dx:
      dH1 = 0.0d0
      dS1 = 0.0d0
      vec_M_dVij(:) = M_dVij(i,j,:)*drij_dsk(1)
      vec_M_dSVij(:) = M_dSVij(i,j,:)*drij_dsk(1)
      call d_KS_sp3d5(vec_M_Vij, vec_M_dVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dH1)	! module "TB_Koster_Slater"
      call d_KS_sp3d5(vec_M_SVij, vec_M_dSVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(4), M_dlmn(7),  dS1)	! module "TB_Koster_Slater"
      dH(1,:,:) = dH1(:,:)
      dS(1,:,:) = dS1(:,:)
         
      ! Derivative along dy:
      dH1 = 0.0d0
      dS1 = 0.0d0
      vec_M_dVij(:) = M_dVij(i,j,:)*drij_dsk(2)
      vec_M_dSVij(:) = M_dSVij(i,j,:)*drij_dsk(2)
      call d_KS_sp3d5(vec_M_Vij, vec_M_dVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dH1)	! module "TB_Koster_Slater"
      call d_KS_sp3d5(vec_M_SVij, vec_M_dSVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(2), M_dlmn(5), M_dlmn(8),  dS1)	! module "TB_Koster_Slater"
      dH(2,:,:) = dH1(:,:)
      dS(2,:,:) = dS1(:,:)
         
      ! Derivative along dz:
      dH1 = 0.0d0
      dS1 = 0.0d0
      vec_M_dVij(:) = M_dVij(i,j,:)*drij_dsk(3)
      vec_M_dSVij(:) = M_dSVij(i,j,:)*drij_dsk(3)
      call d_KS_sp3d5(vec_M_Vij, vec_M_dVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dH1)	! module "TB_Koster_Slater"
      call d_KS_sp3d5(vec_M_SVij, vec_M_dSVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(3), M_dlmn(6), M_dlmn(9),  dS1)	! module "TB_Koster_Slater"
      dH(3,:,:) = dH1(:,:)
      dS(3,:,:) = dS1(:,:)
   endif ! (i == j)
   
   nullify (KOA1, KOA2, rm, x1, y1, z1, r1, al, bl, cl, dl)
end subroutine d_Hamilton_one_NRL



!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! Subroutine to get attractive forces for supercell from the derivatives of the Hamiltonian:
subroutine Attract_TB_Forces_Press_NRL(TB_Hamil, Scell, NSC, numpar, Aij, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, Aij_x_Ei)
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:), intent(in) :: Aij, Aij_x_Ei
   !------------------------------------------------------------------
   real(8), allocatable, dimension(:,:) :: dwr_press, dS_press
   real(8), allocatable, dimension(:,:,:) :: dHij, dSij
   integer i, j, k, n
   if (numpar%p_const) then	! calculate this for P=const Parrinello-Rahman MD
      n = size(Aij,1)
      allocate(dwr_press(9,n))
      allocate(dS_press(9,n))
      allocate(dHij(9,n,size(Aij,2)))
      allocate(dSij(9,n,size(Aij,2)))
      !$OMP WORKSHARE
      dHij = 0.0d0
      dSij = 0.0d0
      dwr_press = 0.0d0
      dS_press = 0.0d0
      !$OMP END WORKSHARE
      
      call dHamil_tot_Press_NRL(Scell(NSC)%MDatoms, Scell, NSC, numpar, TB_Hamil, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij, dSij)
      
      ! 2) Convert the derivative of the Hamiltonian into eVs from Ry:
      dHij = dHij*g_Ry*g_A2au	! [Ry] into [eV], and [1/a0] -> [1/A]
      dSij = dSij*g_A2au	! [1/a0] -> [1/A]

      !$omp PARALLEL DO private(i,j)
      do j = 1, 9
         do i = 1, n 
            dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(i,:)) ! old, tested, good
            dS_press(j,i) = dS_press(j,i) + SUM(dSij(j,i,:)*Aij_x_Ei(i,:))
         enddo ! i
      enddo ! j
      !$OMP END PARALLEL DO
      
      Scell(NSC)%SCforce%att = 0.0d0
      do i = 1,3
         do k = 1,3
            Scell(NSC)%SCforce%att(k,i) = SUM(dwr_press((i-1)*3+k,:)) - SUM(dS_press((i-1)*3+k,:))
!             Scell(NSC)%SCforce%att(i,k) = SUM(dwr_press((i-1)*3+k,:)) - SUM(dS_press((i-1)*3+k,:))
         enddo ! k
      enddo ! i
      deallocate(dwr_press, dS_press, dHij, dSij)
   endif
end subroutine Attract_TB_Forces_Press_NRL


subroutine dHamil_tot_Press_NRL(atoms, Scell, NSC, numpar, TB_Hamil, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij, dSij)
! construct the whole Hamilton matrix:
! (with respect to which Rk we take the derivatives, Appendix F of H.Jeschke PhD Thesis)
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(inout) :: dHij, dSij
   real(8), dimension(9,9,9) :: dHij1, dSij1
   integer :: i, j, j1, i1, ki, atom_2, m, nat, k, i2, j2, NumTB
   integer i4, j4, norb
   norb = 9	! sp3d5 basis set
   nat = size(atoms)
   dHij = 0.0d0
   dSij = 0.0d0
   !$omp PARALLEL DO private(i,m,i4,atom_2,j,j4,j1,i1,j2,i2,dHij1,dSij1)
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
            call dHamilton_one_Press_NRL(i, atom_2, Scell(NSC)%MDatoms, Scell, NSC, TB_Hamil, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij1, dSij1)
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis.
            do j1 = 1,norb	! all orbitals
               j2 = j4+j1
               do i1 = 1,norb	! all orbitals
                  i2 = i4+i1
                  dHij(:,i2,j2) = dHij1(:,i1,j1)	! construct the total Hamiltonian from
                  dSij(:,i2,j2) = dSij1(:,i1,j1)	! construct the total Overlap Matrix from
!                   dHij(:,i2,j2) = dHij1(:,j1,i1)	! construct the total Hamiltonian from
!                   dSij(:,i2,j2) = dSij1(:,j1,I1)	! construct the total Overlap Matrix from
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0) then
      enddo ! j
   enddo ! i
   !$OMP END PARALLEL DO
end subroutine dHamil_tot_Press_NRL ! CHECKED


subroutine dHamilton_one_Press_NRL(i, atom_2, atoms, Scell, NSC, TB, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dHij_press, dSij_press)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   integer(4), intent(IN) :: i, atom_2
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_NRL), dimension(:,:),intent(in), target ::  TB	! all tight binding parameters
   real(8), dimension(:,:,:), intent(in) :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(out) :: dHij_press, dSij_press
   integer :: ki, kj, j, norb
   integer :: i1
   real(8) :: rho, d_rho(9)
   real(8), dimension(9,9,9) ::  dH, dS ! hopping integrals
   real(8), dimension(:), pointer :: al, bl, cl, dl
   integer, pointer :: KOA1, KOA2
   norb = 9	! sp3d5 basis set
   dHij_press = 0.0d0
   if (atom_2 .EQ. 0) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      j = i
      KOA1 => Scell(NSC)%MDatoms(j)%KOA
      KOA2 => Scell(NSC)%MDatoms(i)%KOA

      al => TB(KOA1,KOA2)%al
      bl => TB(KOA1,KOA2)%bl
      cl => TB(KOA1,KOA2)%cl
      dl => TB(KOA1,KOA2)%dl
      ! s orbital of sp3d5 basis set:
      call rho_NRL(Scell, NSC, TB, j, rho)	! to get rho
      call d_On_site_NRL_Press(Scell, NSC, TB, rho, j, al(1), bl(1), cl(1), dl(1), d_rho)	! function below
      dHij_press(:,1,1) = d_rho(:)
      ! p3 orbitals of sp3d5 basis set:
      do i1 = 2, 4
         call d_On_site_NRL_Press(Scell, NSC, TB, rho, j, al(2), bl(2), cl(2), dl(2), d_rho)	! function below
         dHij_press(:,i1,i1) = d_rho(:)
      enddo
      ! d5 orbitals of sp3d5 basis set:
      do i1 = 5, 7
         call d_On_site_NRL_Press(Scell, NSC, TB, rho, j, al(3), bl(3), cl(3), dl(3), d_rho)	! function below
         dHij_press(:,i1,i1) = d_rho(:)
      enddo
      do i1 = 8, 9
         if (TB(KOA1,KOA2)%ind_split == 1) then	! splitting between t2g and e2
            call d_On_site_NRL_Press(Scell, NSC, TB, rho, j, al(4), bl(4), cl(4), dl(4), d_rho)	! function below
         else	! no splitting between t2g and e2
            call d_On_site_NRL_Press(Scell, NSC, TB, rho, j, al(3), bl(3), cl(3), dl(3), d_rho)	! function below
         endif
         dHij_press(:,i1,i1) = d_rho(:)
      enddo
      dSij_press = 0.0d0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_Press_NRL(i, atom_2, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn,  dH, dS)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1, norb
         do kj = 1, norb
            dHij_press(:,kj,ki) = dH(:,ki,kj)   ! Hopping Integrals
            dSij_press(:,kj,ki) = dS(:,ki,kj)   ! Hopping Integrals
!             dHij_press(:,ki,kj) = dH(:,ki,kj)   ! Hopping Integrals
!             dSij_press(:,ki,kj) = dS(:,ki,kj)   ! Hopping Integrals
         enddo ! kj
      enddo  ! ki
   endif
   
   nullify(al, bl, cl, dl, KOA1, KOA2)
end subroutine dHamilton_one_Press_NRL



subroutine dHopping_Press_NRL(i, atom_2, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij, M_lmn, dH, dS)
! subroutine making the derivatives of the hopping integrals
   integer, intent(in) :: i, atom_2	! atoms
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   real(8), dimension(:,:,:), intent(in), target :: M_Vij, M_dVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in), target :: M_SVij, M_dSVij	! matrix of Overlap functions for all pairs of atoms, all orbitals, and derivatives
   real(8), dimension(:,:,:), intent(in) :: M_lmn	! matrix of directional cosines l, m, n; and derivatives
   real(8), dimension(:,:,:), intent(out) :: dH, dS
   !---------------------------------------------
   integer :: i1, j1, k1
   real(8) :: drij_dh, rij(3), sij(3)
   real(8), dimension(3) :: M_dlmn
   real(8), dimension(10) :: vec_M_dVij, vec_M_dSVij
   real(8), dimension(10) ::  vec_M_Vij, vec_M_SVij
   real(8), dimension(9,9) :: dH1, dS1
   real(8), pointer :: r
   integer, pointer :: j

   dH = 0.0d0
   dS = 0.0d0

   !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
   j => Scell(NSC)%Near_neighbor_list(i,atom_2)   ! it interacts with this atom
   r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
   rij(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
   rij(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
   rij(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
   sij(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
   sij(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
   sij(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

   vec_M_Vij = M_Vij(i,j,:)
   vec_M_SVij = M_SVij(i,j,:)
   
   do i1 = 1, 3	! gamma
      do j1 = 1, 3	! delta
         k1 = (i1-1)*3 + j1
         ! all the components of the h_alpha_beta(3,3):
         ! TESTED, CORRECT:
         M_dlmn(1) = dda_dhgd(1, j1, rij(1), rij(j1), sij(i1), r)	! dl/dh{gamma,delta}, module "TB_Koster_Slater"
         M_dlmn(2) = dda_dhgd(2, j1, rij(2), rij(j1), sij(i1), r)	! dm/dh{gamma,delta}, module "TB_Koster_Slater"
         M_dlmn(3) = dda_dhgd(3, j1, rij(3), rij(j1), sij(i1), r)	! dn/dh{gamma,delta}, module "TB_Koster_Slater"
!          M_dlmn = -M_dlmn*g_au2A
         M_dlmn = M_dlmn*g_au2A
         
         ! TESTED, CORRECT:
         drij_dh = drij_dhab(rij(j1), sij(i1), r)	! dr_{ij}/dh_{gamma,delta}, module "TB_Koster_Slater"

         vec_M_dVij = M_dVij(i,j,:)*drij_dh
         vec_M_dSVij = M_dSVij(i,j,:)*drij_dh
         
!          call d_KS_sp3d5(vec_M_Vij, vec_M_dVij, -M_lmn(1,i,j), -M_lmn(2,i,j), -M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3), dH1)	! module "TB_Koster_Slater"
!          call d_KS_sp3d5(vec_M_SVij, vec_M_dSVij, -M_lmn(1,i,j), -M_lmn(2,i,j), -M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3), dS1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_TEST(vec_M_Vij, vec_M_dVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3), dH1)	! module "TB_Koster_Slater"
         call d_KS_sp3d5_TEST(vec_M_SVij, vec_M_dSVij, M_lmn(1,i,j), M_lmn(2,i,j), M_lmn(3,i,j), M_dlmn(1), M_dlmn(2), M_dlmn(3), dS1)	! module "TB_Koster_Slater"
!          dH(k1,:,:) = TRANSPOSE(dH1(:,:))
!          dS(k1,:,:) = TRANSPOSE(dS1(:,:))
         dH(k1,:,:) = dH1(:,:)
         dS(k1,:,:) = dS1(:,:)
      enddo
   enddo
   
   nullify(r, j)
end subroutine dHopping_Press_NRL




!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Derivative of the density of on-site density function:
subroutine d_rho_NRL(Scell, NSC, TB_Hamil, j, k, d_rho)	! d rho/drij
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB_Hamil   ! parameters of the Hamiltonian of TB
   integer, intent(in) :: j,k 
   real(8), dimension(3), intent(out) :: d_rho	 ! "Density"
   !-----------------
   integer :: i, atom_2, m
   real(8) :: r, drij_dsk(3), d_exp
   real(8), pointer :: lambd, rc, lden, r1, x1, y1, z1
   integer, pointer :: KOA1, KOA2

   d_rho = 0.0d0
   m = Scell(NSC)%Near_neighbor_size(j)
!    !$omp parallel private(i,atom_2,r,lambd,rc,lden,KOA1,KOA2, r1,x1,y1,z1,drij_dsk,d_exp)
!    !$omp do reduction( + : d_rho)
   do atom_2 = 1,m ! do only for atoms close to that one  
      i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
      KOA1 => Scell(NSC)%MDatoms(j)%KOA
      KOA2 => Scell(NSC)%MDatoms(i)%KOA
      lambd => TB_Hamil(KOA1,KOA2)%lambd
      rc => TB_Hamil(KOA1,KOA2)%Rc
      lden => TB_Hamil(KOA1,KOA2)%lden
      x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1)	! at this distance, X
      y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2)	! at this distance, Y
      z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3)	! at this distance, Z
      r1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,4)		! at this distance, R

      ! Derivatives of rij by sk:
      drij_dsk(1) = drij_dska(j, i, k, x1, y1, z1, r1, Scell(NSC)%supce, 1, .true.)	! dr_{ij}/ds_{k,x}, module "TB_Koster_Slater"
      drij_dsk(2) = drij_dska(j, i, k, x1, y1, z1, r1, Scell(NSC)%supce, 2, .true.)	! dr_{ij}/ds_{k,y}, module "TB_Koster_Slater"
      drij_dsk(3) = drij_dska(j, i, k, x1, y1, z1, r1, Scell(NSC)%supce, 3, .true.)	! dr_{ij}/ds_{k,z}, module "TB_Koster_Slater"

      ! In this formula, r must be in Bohr, so convert into it from Angstroms:
      r = r1*g_A2au

      ! Convert [A] -> [Bohr]:
      drij_dsk(:) = drij_dsk(:)*g_A2au

      d_exp = d_exp_F_NRL(r, lambd, rc, lden)	! function below
      d_rho(:) = d_rho(:) + d_exp*drij_dsk(:)
   enddo ! atom_2
!    !$omp end do
!    !$omp end parallel
   
   nullify(lambd, rc, lden, KOA1, KOA2, r1, x1, y1, z1)
end subroutine d_rho_NRL


subroutine d_On_site_NRL(Scell, NSC, TB_Hamil, rho, j, k, al, bl, cl, dl, dhil)! dh/dsk
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC, j, k ! number of supercell, number of atom
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(3) :: dhil
   real(8), intent(in) :: rho, al, bl, cl, dl
   real(8) :: dpho_drij(3), rho13, temp
   if (rho > 0.0d0) then	! there are neighbours
      rho13 = rho**m_one_third
      temp = m_two_third*bl/rho13+ m_four_third*cl*rho13 + 2.0d0*dl*rho		! Eq.(17) from [1]
      call d_rho_NRL(Scell, NSC, TB_Hamil, j, k, dpho_drij)	! d rho/drij
   else	! there are no neighbours => no forces
      temp = 0.0d0
      dpho_drij(:) = 0.0d0
   endif
   dhil(:) = temp*dpho_drij(:)
end subroutine d_On_site_NRL


subroutine d_On_site_NRL_Press(Scell, NSC, TB_Hamil, rho, j, al, bl, cl, dl, dhil)! dh/dsk
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC, j	! number of supercell, number of atom
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(9) :: dhil
   real(8), intent(in) :: rho, al, bl, cl, dl
   real(8) :: dpho_dhab(9), rho13, temp
   if (rho > 0.0d0) then	! there are neighbours
      rho13 = rho**m_one_third
      temp = m_two_third*bl/rho13+ m_four_third*cl*rho13 + 2.0d0*dl*rho		! Eq.(17) from [1]
      call d_rho_NRL_Press(Scell, NSC, TB_Hamil, j, dpho_dhab)	! d rho/drij
   else	! there are no neighbours => no forces
      temp = 0.0d0
      dpho_dhab(:) = 0.0d0
   endif
   dhil(:) = temp*dpho_dhab(:)
end subroutine d_On_site_NRL_Press


subroutine d_rho_NRL_Press(Scell, NSC, TB_Hamil, j, d_rho)	! d rho/drij
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB_Hamil   ! parameters of the Hamiltonian of TB
   integer, intent(in) :: j
   real(8), dimension(9), intent(out) :: d_rho	 ! "Density" derivative by h_{alpha,beta}
   !-----------------
   integer :: i, atom_2, m
   real(8) :: r, drij_dhab_c(9), d_exp
   real(8), pointer :: lambd, rc, lden, r1, x1, y1, z1, sx1, sy1, sz1
   integer, pointer :: KOA1, KOA2

   d_rho = 0.0d0
   m = Scell(NSC)%Near_neighbor_size(j)
   do atom_2 = 1,m ! do only for atoms close to that one
      i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
      KOA1 => Scell(NSC)%MDatoms(j)%KOA
      KOA2 => Scell(NSC)%MDatoms(i)%KOA
      lambd => TB_Hamil(KOA1,KOA2)%lambd
      rc => TB_Hamil(KOA1,KOA2)%Rc
      lden => TB_Hamil(KOA1,KOA2)%lden
      x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1)	! at this distance, X
      y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2)	! at this distance, Y
      z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3)	! at this distance, Z
      r1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,4)		! at this distance, R
      sx1 => Scell(NSC)%Near_neighbor_dist_s(j,atom_2,1)	! at this distance, SX
      sy1 => Scell(NSC)%Near_neighbor_dist_s(j,atom_2,2)	! at this distance, SY
      sz1 => Scell(NSC)%Near_neighbor_dist_s(j,atom_2,3)	! at this distance, SZ

      ! Derivatives of rij by hab:
      drij_dhab_c(1) =  drij_dhab(x1, sx1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(2) =  drij_dhab(y1, sx1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(3) =  drij_dhab(z1, sx1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(4) =  drij_dhab(x1, sy1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(5) =  drij_dhab(y1, sy1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(6) =  drij_dhab(z1, sy1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(7) =  drij_dhab(x1, sz1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(8) =  drij_dhab(y1, sz1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      drij_dhab_c(9) =  drij_dhab(z1, sz1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
      
!       drij_dhab_c(1) =  drij_dhab(x1, sx1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(2) =  drij_dhab(x1, sy1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(3) =  drij_dhab(x1, sz1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(4) =  drij_dhab(y1, sx1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(5) =  drij_dhab(y1, sy1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(6) =  drij_dhab(y1, sz1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(7) =  drij_dhab(z1, sx1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(8) =  drij_dhab(z1, sy1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"
!       drij_dhab_c(9) =  drij_dhab(z1, sz1, r1)	! dr_{ij}/dh_{a,b}, module "TB_Koster_Slater"

      ! In this formula, r must be in Bohr, so convert into it from Angstroms:
      r = r1*g_A2au

      d_exp = d_exp_F_NRL(r, lambd, rc, lden)	! function below
      d_rho(:) = d_rho(:) + d_exp*drij_dhab_c(:)
   enddo ! atom_2

   nullify(lambd, rc, lden, KOA1, KOA2, r1, x1, y1, z1, sx1, sy1, sz1)
end subroutine d_rho_NRL_Press


pure function d_Overlap_function_NRL(r, pllm, qllm, rllm, lambd, rc, dl, ind, ll) result(Hllm)	! Eq.(19), (20) or (23) from [1], dH/drij
   real(8) :: Hllm
   real(8), intent(in) :: r, pllm, qllm, rllm, lambd, rc, dl
   integer, intent(in) :: ind	! 0=Eq.(19) or (20);  1=Eq.(23)
   real(8), intent(in) :: ll	! delta (l l')
   !-------------------------------
   real(8) :: r2
   select case (ind)
   case (1)
      r2 = r*r
      Hllm = (pllm + 2.0d0*qllm*r + 3.0d0*rllm*r2)*exp_F_NRL(r, lambd, rc, dl)	+ (ll + pllm*r + qllm*r2 + rllm*r*r2)*d_exp_F_NRL(r, lambd, rc, dl)
   case default
      Hllm = (qllm + 2.0d0*rllm*r)*exp_F_NRL(r, lambd, rc, dl) + (pllm + qllm*r + rllm*r*r)*d_exp_F_NRL(r, lambd, rc, dl)
   endselect
end function d_Overlap_function_NRL


pure function d_exp_F_NRL(r, lambd, rc, dl) result(exp_F)	! part of Eq.(15), (19), (20), and (23) from [1], d(exp*F)/drij
   real(8) :: exp_F
   real(8), intent(in) :: r, rc, dl, lambd
   real(8) :: exp_r, d_exp_r, l2
   if (r >= rc) then
      exp_F = 0.0d0
   else
      l2 = lambd*lambd
      exp_r = exp(-l2*r)
      exp_F = exp_r*(d_F_NRL(r, rc, dl) - l2*F_NRL(r, rc, dl))	! functions below
   endif
end function d_exp_F_NRL



!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Repulsive part of TB:

subroutine get_Erep_s_NRL(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_NRL), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   a = 0.0d0	! There is no repulsive term in NRL parameterization
end subroutine get_Erep_s_NRL


subroutine dErdr_s_NRL(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   type(TB_Rep_NRL), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
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
END subroutine dErdr_s_NRL


subroutine dErdr_Pressure_s_NRL(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_NRL), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms ! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      Scell(NSC)%SCforce%rep = 0.0d0
   endif
end subroutine dErdr_Pressure_s_NRL

!oooooooooooooooooooooooooooooooooooooooo
! Obsolete subroutines:

subroutine Hamilton_one_NRL_unoptimized(Scell, NSC, i, j, atoms, TB, Hij, x1, y1, z1, sx1, sy1, sz1, cell_x, cell_y, cell_z)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, intent(IN) :: i, j
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_NRL), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   real(8), dimension(:,:), intent(out) :: Hij  ! hamiltonian, all orbitals in sp3d5 basis set
   real(8), intent(out), optional :: x1, y1, z1    ! shortest distances
   real(8), intent(out), optional :: sx1, sy1, sz1 ! shortest distance in relative coordinates
   integer, intent(out), optional :: cell_x, cell_y, cell_z ! cell numbers
   !---------------------------------------
   integer(4) ki, kj, i1, j1, k1, ik, n_orb
   real(8) x,y,z,r,r1, x0, y0, z0, sx, sy, sz !  interatomic distance projections, and the total distance
   integer :: cell_x1, cell_y1, cell_z1 ! cell numbers
   real(8), dimension(3) :: zb
   real(8) :: rho	! an embedded-atom-like ‘density’
   real(8), dimension(:), pointer :: al, bl, cl, dl	! set of parameters in Eq.(17)
   integer, pointer :: KOA1, KOA2

   n_orb = 9	! number of arbitals for sp3d5 basis set

   KOA1 => Scell(NSC)%MDatoms(j)%KOA
   KOA2 => Scell(NSC)%MDatoms(i)%KOA

   if (i == j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      Hij = 0.0d0   ! Nondiagonals are zeros

      al => TB(KOA1,KOA2)%al
      bl => TB(KOA1,KOA2)%bl
      cl => TB(KOA1,KOA2)%cl
      dl => TB(KOA1,KOA2)%dl
      
      ! s orbital of sp3d5 basis set:
      call rho_NRL(Scell, NSC, TB, j, rho)	! to get rho
      Hij(1,1) = On_site_NRL(rho, al(1), bl(1), cl(1), dl(1)) 	! function below
      ! p3 orbitals of sp3d5 basis set:
      do i1 = 2, 4
         Hij(i1,i1) = On_site_NRL(rho, al(2), bl(2), cl(2), dl(2)) 	! function below
      enddo
      ! d5 orbitals of sp3d5 basis set:
      do i1 = 5, 7
         Hij(i1,i1) = On_site_NRL(rho, al(3), bl(3), cl(3), dl(3)) 	! function below
      enddo
      do i1 = 8, 9
         if (TB(KOA1,KOA2)%ind_split == 1) then	! splitting between t2g and e2
            Hij(i1,i1) = On_site_NRL(rho, al(4), bl(4), cl(4), dl(4)) 	! function below
         else	! no splitting between t2g and e2
            Hij(i1,i1) = On_site_NRL(rho, al(3), bl(3), cl(3), dl(3)) 	! function below
         endif
      enddo
      
      if (present(x1)) x1 = 0.0d0
      if (present(y1)) y1 = 0.0d0
      if (present(z1)) z1 = 0.0d0
      if (present(sx1)) sx1 = 0.0d0
      if (present(sy1)) sy1 = 0.0d0
      if (present(sz1)) sz1 = 0.0d0
      if (present(cell_x)) cell_x = 0
      if (present(cell_y)) cell_y = 0
      if (present(cell_z)) cell_z = 0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call shortest_distance(Scell, NSC, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz, cell_x=cell_x1, cell_y=cell_y1, cell_z=cell_z1)
      if (present(x1)) x1 = x
      if (present(y1)) y1 = y
      if (present(z1)) z1 = z
      if (present(sx1)) sx1 = sx
      if (present(sy1)) sy1 = sy
      if (present(sz1)) sz1 = sz
      if (present(cell_x)) cell_x = cell_x1
      if (present(cell_y)) cell_y = cell_y1
      if (present(cell_z)) cell_z = cell_z1
      
      ! In this formula, r must be in Bohr, so convert into it from Angstroms:
      x = x*g_A2au
      y = y*g_A2au
      z = z*g_A2au
      r = r*g_A2au
      
      call Hopping_NRL_unoptimized(TB(KOA1,KOA2), Hij, x, y, z, r)	! subroutine below
      
!       if ((i == 23) .and. (j == 1)) then
!         print*, 'r', x*g_au2A, y*g_au2A, z*g_au2A, r*g_au2A, sx, sy, sz, cell_x1, cell_y1, cell_z1
!          do ki = 1, n_orb
!            do kj = 1, n_orb
!                  print*, 'AA', i, j, ki, kj, Hij(kj,ki)
! !              !Hij(kj,ki) = ts(kj,ki)   ! Hopping Integrals
!            enddo ! kj
!          enddo  ! ki
!          pause 'AFTER Hopping_NRL'
!        endif
!        
!       if ((j == 23) .and. (i == 1)) then
!         print*, 'r', x*g_au2A, y*g_au2A, z*g_au2A, r*g_au2A, sx, sy, sz, cell_x1, cell_y1, cell_z1
!         do ki = 1, n_orb
!           do kj = 1, n_orb
!                 print*, 'BB', i, j, ki, kj, Hij(kj,ki)*g_Ry
!              !Hij(kj,ki) = ts(kj,ki)   ! Hopping Integrals
!           enddo ! kj
!         enddo  ! ki
!         pause 'AFTER Hopping_NRL'
!        endif
       
   endif
   nullify (al, bl, cl, dl, KOA1, KOA2)
end subroutine Hamilton_one_NRL_unoptimized


subroutine Hopping_NRL_unoptimized(TB, ts, x, y, z, r0)
! subroutine making the hopping integrals
! ts --- the Hopping Integrals themselves
! x,y,z --- the distances between naighbour atoms we are analizing
   type(TB_H_NRL), intent(in), target :: TB	! all tight binding parameters
   REAL(8), DIMENSION(9,9), INTENT(out) :: ts
   real(8), intent(in) :: x,y,z	! relative distances (projections) between the atoms
   real(8), intent(in), optional :: r0
   !=============================================
   real(8), pointer :: rm
   real(8) r, r1, xr, yr, zr, Vsr(10)
   integer :: i

   if (present(r0)) then
      r = r0 ! it is given, no need to recalculate
   else
      r = SQRT(x*x + y*y + z*z)   ! total distance between the two atoms within the supercell
   endif
   
   rm => TB%Rc ! cut-off radius
   
   if (r <= rm) then   ! these atoms do interact:
      ! Get the parameterized radial part of the overlap integrals:
!       do i = 1, size(Vsr)	! for all coefficients
!          Vsr(i) = Overlap_function_NRL(r, TB%ellm(i), TB%fllm(i), TB%gllm(i), TB%hllm(i), rm, TB%lden, 0, 0.0d0)	! function below
!       enddo
      
      Vsr(1) = Overlap_function_NRL(r, TB%ellm(1), TB%fllm(1), TB%gllm(1), TB%hllm(1), rm, TB%lden, 0, 0.0d0)	!    V(1) = (s s sigma)
      Vsr(2) = Overlap_function_NRL(r, TB%ellm(2), TB%fllm(2), TB%gllm(2), TB%hllm(2), rm, TB%lden, 0, 0.0d0)	!    V(2) = (s p sigma)
      Vsr(3) = Overlap_function_NRL(r, TB%ellm(5), TB%fllm(5), TB%gllm(5), TB%hllm(5), rm, TB%lden, 0, 0.0d0)	!    V(3) = (s d sigma)
      Vsr(4) = Overlap_function_NRL(r, TB%ellm(3), TB%fllm(3), TB%gllm(3), TB%hllm(3), rm, TB%lden, 0, 0.0d0)	!    V(4) = (p p sigma)
      Vsr(5) = Overlap_function_NRL(r, TB%ellm(4), TB%fllm(4), TB%gllm(4), TB%hllm(4), rm, TB%lden, 0, 0.0d0)	!    V(5) = (p p pi)
      Vsr(6) = Overlap_function_NRL(r, TB%ellm(6), TB%fllm(6), TB%gllm(6), TB%hllm(6), rm, TB%lden, 0, 0.0d0)	!    V(6) = (p d sigma)
      Vsr(7) = Overlap_function_NRL(r, TB%ellm(7), TB%fllm(7), TB%gllm(7), TB%hllm(7), rm, TB%lden, 0, 0.0d0)	!    V(7) = (p d pi)
      Vsr(8) = Overlap_function_NRL(r, TB%ellm(8), TB%fllm(8), TB%gllm(8), TB%hllm(8), rm, TB%lden, 0, 0.0d0)	!    V(8) = (d d sigma)
      Vsr(9) = Overlap_function_NRL(r, TB%ellm(9), TB%fllm(9), TB%gllm(9), TB%hllm(9), rm, TB%lden, 0, 0.0d0)	!    V(9) = (d d pi)
      Vsr(10) = Overlap_function_NRL(r, TB%ellm(10), TB%fllm(10), TB%gllm(10), TB%hllm(10), rm, TB%lden, 0, 0.0d0)	!    V(10) = (d d delta)

!       print* , 'Vsr'
!       write(*,'(f,$)') Vsr(:)
!       pause ' -- TEST Vsr'
   
      xr = x/r
      yr = y/r
      zr = z/r
      ! Construct the overlap integrals including angular part for sp3d5 basis set:
      call  KS_sp3d5(Vsr, xr, yr, zr, ts)	! module TB_Koster_Slater
      
   else ! these atoms too far to interact via covalent bonding:
      ts = 0.0d0
   endif

   nullify(rm)
end subroutine Hopping_NRL_unoptimized



subroutine Get_overlap_S_matrix_unoptimized(Scell, NSC, i, atom_2, j, atoms, TB, Sij)
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer, INTENT(IN) :: i, j, atom_2
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_NRL),  intent(in), target :: TB	! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: Sij  ! hamiltonian
   !---------------------------------------
   real(8), dimension(10) :: Vsr
   real(8) x,y,z,r, xr, yr, zr
   real(8), pointer :: rm, x1, y1, z1, r1
   integer :: k, k2

   rm => TB%Rc ! cut-off radius
   
!    call shortest_distance(Scell, NSC, atoms, i, j, r, x1=x, y1=y, z1=z)	! module "Atomic_tools"
   if (atom_2 == 0) then	! it's the same atom
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
      r = 0.0d0
   else
      x1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,1) ! at this distance, X
      y1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,2) ! at this distance, Y
      z1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,3) ! at this distance, Z
      r1 => Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R
      ! In this formula, r must be in Bohr, so convert into it from Angstroms:
      x = x1*g_A2au
      y = y1*g_A2au
      z = z1*g_A2au
      r = r1*g_A2au
   endif

   if (r <= rm) then   ! these atoms do interact:
 
      ! Get the parameterized radial part of the overlap matrix:
!       do k = 1, 10
!          write(*,'(i2,f,f,f,f)') k, TB%pllm(k), TB%qllm(k), TB%rllm(k), TB%sllm(k)
!       enddo
!       PAUSE 'Overlap'
      
      ! (The order of Vsr does not coincide with the order of pllm,qllm,rllm,sllm due to the order of parameters in the input files from NRL database)
      Vsr(1) = Overlap_function_NRL(r, TB%pllm(1), TB%qllm(1), TB%rllm(1), TB%sllm(1), rm, TB%lden, TB%ind_overlap, 1.0d0)	!    V(1) = (s s sigma)
      Vsr(2) = Overlap_function_NRL(r, TB%pllm(2), TB%qllm(2), TB%rllm(2), TB%sllm(2), rm, TB%lden, TB%ind_overlap, 0.0d0)	!    V(2) = (s p sigma)
      Vsr(3) = Overlap_function_NRL(r, TB%pllm(5), TB%qllm(5), TB%rllm(5), TB%sllm(5), rm, TB%lden, TB%ind_overlap, 0.0d0)	!    V(3) = (s d sigma)
      Vsr(4) = Overlap_function_NRL(r, TB%pllm(3), TB%qllm(3), TB%rllm(3), TB%sllm(3), rm, TB%lden, TB%ind_overlap, 1.0d0)	!    V(4) = (p p sigma)
      Vsr(5) = Overlap_function_NRL(r, TB%pllm(4), TB%qllm(4), TB%rllm(4), TB%sllm(4), rm, TB%lden, TB%ind_overlap, 1.0d0)	!    V(5) = (p p pi)
      Vsr(6) = Overlap_function_NRL(r, TB%pllm(6), TB%qllm(6), TB%rllm(6), TB%sllm(6), rm, TB%lden, TB%ind_overlap, 0.0d0)	!    V(6) = (p d sigma)
      Vsr(7) = Overlap_function_NRL(r, TB%pllm(7), TB%qllm(7), TB%rllm(7), TB%sllm(7), rm, TB%lden, TB%ind_overlap, 0.0d0)	!    V(7) = (p d pi)
      Vsr(8) = Overlap_function_NRL(r, TB%pllm(8), TB%qllm(8), TB%rllm(8), TB%sllm(8), rm, TB%lden, TB%ind_overlap, 1.0d0)	!    V(8) = (d d sigma)
      Vsr(9) = Overlap_function_NRL(r, TB%pllm(9), TB%qllm(9), TB%rllm(9), TB%sllm(9), rm, TB%lden, TB%ind_overlap, 1.0d0)	!    V(9) = (d d pi)
      Vsr(10) = Overlap_function_NRL(r, TB%pllm(10), TB%qllm(10), TB%rllm(10), TB%sllm(10), rm, TB%lden, TB%ind_overlap, 1.0d0)	!    V(10) = (d d delta)
      
      if (i /= j) then	! it's 2 different atoms:
         xr = x/r
         yr = y/r
         zr = z/r
         ! Construct the overlap matrix including angular part for sp3d5 basis set:
         call  KS_sp3d5(Vsr, xr, yr, zr, Sij)	! module TB_Koster_Slater
      else	! it is the same atom
         ! no angular part:
!          call  KS_sp3d5(Vsr, 0.0d0, 0.0d0, 0.0d0, Sij)	! module TB_Koster_Slater
!          do k = 1, 9
!             do k2 = 1, 9
!                Sij(k,k2) = 0.0d0
!             enddo
!          enddo
         Sij = 0.0d0
         forall (k=1:9) Sij(k,k)=1.0d0
      endif
   else	! too far atoms do not overlap:
      Sij = 0.0d0   ! atmos too far apart do not overlap
   endif

   nullify(rm, x1, y1, z1, r1)
end subroutine Get_overlap_S_matrix_unoptimized


END MODULE TB_NRL
