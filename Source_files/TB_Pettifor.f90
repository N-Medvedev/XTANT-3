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
! This module contains subroutines to deal with TB hamiltonian in the Pettifor parametrization

MODULE TB_Pettifor
use Universal_constants
use Objects
use TB_Koster_Slater
use Algebra_tools, only : Kronecker_delta, sym_diagonalize, Reciproc
use Atomic_tools, only : get_fraction_of_given_sort, shortest_distance, Reciproc_rel_to_abs
use Electron_tools, only : find_band_gap

implicit none
PRIVATE

public :: Construct_M_Vs, Complex_Hamil_tot, dHij_s, Attract_TB_Forces_Press, dErdr_s, dErdr_Pressure_s, &
         construct_TB_H_Pettifor, get_Erep_s, dHij_r, dE2rep_dr2

 contains


!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Attractive part of TB:


! Second derivative of the TB Hamiltonian by r:
subroutine dHij_r(TB_Hamil, atoms, Scell, numpar, M_Vs, M_dVs, M_d2Vs, M_cos, F, dF)
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   real(8), dimension(:,:,:), intent(in) :: M_d2Vs ! matrix of functions d2Vs
   real(8), dimension(:,:,:), intent(in) :: M_cos	! matrix of directional cosines
   real(8), dimension(:,:), allocatable, intent(inout) :: F, dF	! forces and derivatives by r, [eV/A] and [eV/A^2]
   !------------------------------
   REAL(8), DIMENSION(3) :: Eelectr_r, Eelectr_2r
   real(8), dimension(:,:), allocatable :: dHijx_r_all ! allocatable
   real(8), dimension(:,:), allocatable :: dHijy_r_all ! allocatable
   real(8), dimension(:,:), allocatable :: dHijz_r_all ! allocatable
   real(8), dimension(:,:), allocatable :: d2Hijx_r2_all, d2Hijy_r2_all, d2Hijz_r2_all	! For derivatives of the forces
   integer :: N, Nat, k
   Nat = size(Scell(1)%MDatoms)	! number of atoms
   N = size(Scell(1)%Ei)	! number of the energy levels
   if (.not.allocated(F)) allocate(F(3,N))
   F = 0.0d0
   if (.not.allocated(dF)) allocate(dF(3,N))
   dF = 0.0d0
   
!$omp PARALLEL private(k, Eelectr_r, Eelectr_2r, dHijx_r_all, dHijy_r_all, dHijz_r_all, d2Hijx_r2_all, d2Hijy_r2_all, d2Hijz_r2_all) 
   if (.not.allocated(dHijx_r_all)) allocate(dHijx_r_all(N,N))
   if (.not.allocated(dHijy_r_all)) allocate(dHijy_r_all(N,N))
   if (.not.allocated(dHijz_r_all)) allocate(dHijz_r_all(N,N))
   if (.not.allocated(d2Hijx_r2_all)) allocate(d2Hijx_r2_all(N,N))
   if (.not.allocated(d2Hijy_r2_all)) allocate(d2Hijy_r2_all(N,N))
   if (.not.allocated(d2Hijz_r2_all)) allocate(d2Hijz_r2_all(N,N))   
!$omp do
   do k = 1, Nat	! forces and derivatives for atoms:
      call dHamil_tot_r(dHijx_r_all, dHijy_r_all, dHijz_r_all, d2Hijx_r2_all, d2Hijy_r2_all, d2Hijz_r2_all, TB_Hamil, Scell, 1, numpar, k, M_Vs, M_dVs, M_d2Vs, M_cos) ! see below
      call Attract_TB_H3_near(Scell(1)%Aij, dHijx_r_all, dHijy_r_all, dHijz_r_all, Scell, 1, Eelectr_r)	! see below
      call Attract_TB_H3_near(Scell(1)%Aij, d2Hijx_r2_all, d2Hijy_r2_all, d2Hijz_r2_all, Scell, 1, Eelectr_2r)	! see below
      F(:,k) = Eelectr_r(:)	! save attractive forces [eV/A]
      dF(:,k) = Eelectr_2r(:)	! save derivatives of attractive forces [eV/A^2]
   enddo ! k
!$omp end do 
   if (allocated(dHijx_r_all)) deallocate(dHijx_r_all)
   if (allocated(dHijy_r_all)) deallocate(dHijy_r_all)
   if (allocated(dHijz_r_all)) deallocate(dHijz_r_all)
   if (allocated(d2Hijx_r2_all)) deallocate(d2Hijx_r2_all)
   if (allocated(d2Hijy_r2_all)) deallocate(d2Hijy_r2_all)
   if (allocated(d2Hijz_r2_all)) deallocate(d2Hijz_r2_all)
!$omp end parallel
end subroutine dHij_r


!11111111111111111111111111111111111111111111111111111111111111111111111
! Derivative of the atomic TB Hamiltonian:
subroutine dHij_s(TB_Hamil, atoms, Scell, NSC, numpar, Aij, M_x1, M_xrr) ! attractive forces for atoms, module "TB_Hamiltonian"
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors   
   REAL(8), DIMENSION(:,:), INTENT(in) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis
   REAL(8), DIMENSION(:,:,:), INTENT(in) :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in) :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   !========================================= 
   REAL(8), DIMENSION(3) :: Eelectr_s
   real(8), dimension(:,:), allocatable :: dHijx_s_all ! allocatable
   real(8), dimension(:,:), allocatable :: dHijy_s_all ! allocatable
   real(8), dimension(:,:), allocatable :: dHijz_s_all ! allocatable
   real(8), dimension(:,:,:), allocatable :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), allocatable :: M_dVs ! matrix of functions dVs
   integer k, nat, nat4, my_id, OMP_GET_THREAD_NUM, num_th, OMP_SET_NUM_THREADS
   nat = size(atoms)
   nat4 = nat*4 ! sp^3 parametrization has 4 orbitals per atom

   ! Construct array of functions Vs and dVs for all pairs of atoms to use for forces:
   call Construct_M_Vs(Scell, NSC, TB_Hamil, M_Vs, M_dVs) ! subroitine below


!$omp PARALLEL private(k, Eelectr_s, dHijx_s_all, dHijy_s_all, dHijz_s_all) 
   if (.not.allocated(dHijx_s_all)) allocate(dHijx_s_all(nat4,nat4))
   if (.not.allocated(dHijy_s_all)) allocate(dHijy_s_all(nat4,nat4))
   if (.not.allocated(dHijz_s_all)) allocate(dHijz_s_all(nat4,nat4))   
!$omp do
   do k = 1,nat ! initial conditions for atoms:
!       my_id = OMP_GET_THREAD_NUM() ! identify which thread it is
!       num_th = OMP_GET_NUM_THREADS() ! get nmber of threads available
      
      Scell(NSC)%MDatoms(k)%forces%att(:) = 0.0d0
      call dHamil_tot_s(dHijx_s_all, dHijy_s_all, dHijz_s_all, TB_Hamil, Scell, NSC, numpar, k, M_x1, M_xrr, M_Vs, M_dVs) ! see below

      call Attract_TB_H3_near(Aij, dHijx_s_all, dHijy_s_all, dHijz_s_all, Scell, NSC, Eelectr_s) ! see below
      Scell(NSC)%MDatoms(k)%forces%att(:) = Eelectr_s(:) ! save attractive forces

   enddo ! k
!$omp end do 
   if (allocated(dHijx_s_all)) deallocate(dHijx_s_all)
   if (allocated(dHijy_s_all)) deallocate(dHijy_s_all)
   if (allocated(dHijz_s_all)) deallocate(dHijz_s_all)
!$omp end parallel

   if (allocated(M_Vs)) deallocate(M_Vs)
   if (allocated(M_dVs)) deallocate(M_dVs)
end subroutine dHij_s



subroutine Construct_M_Vs(Scell, NSC, TB_Hamil, M_Vs, M_dVs, M_d2Vs, M_cos)
   type(Super_cell), dimension(:), intent(in), target :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_Vs	! matrix of functions Vs
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_dVs	! matrix of functions dVs, first derivatives
   real(8), dimension(:,:,:), allocatable, intent(out), optional :: M_d2Vs	! matrix of functions d2Vs, second rerivaties
   real(8), dimension(:,:,:), allocatable, intent(out), optional :: M_cos	! matrix of cosine directions d^{alpha}
   !---------------------------
   real(8), pointer :: r, x, y, z
   integer, pointer :: nat,  m, j
   integer i, atom_2, N, k
   nat => Scell(NSC)%Na ! number of atoms
   N = size(TB_Hamil(Scell(NSC)%MDatoms(1)%KOA, Scell(NSC)%MDatoms(1)%KOA)%V0) ! number of orbitals per atom
   
   if (.not.allocated(M_Vs)) allocate(M_Vs(N,nat,nat))	! for all orbitals and pairs of atoms
   M_Vs = 0.0d0
   if (.not.allocated(M_dVs)) allocate(M_dVs(N,nat,nat))	! for all orbitals and pairs of atoms
   M_dVs = 0.0d0
   if (present(M_d2Vs)) then
      if (.not.allocated(M_d2Vs)) allocate(M_d2Vs(N,nat,nat))	! for all orbitals and pairs of atoms
      M_d2Vs = 0.0d0
   endif
   if (present(M_cos)) then
      if (.not.allocated(M_cos)) allocate(M_cos(nat,nat,3))	! for 3 directions (x, y, z)  and pairs of atoms
      M_cos = 0.0d0
   endif
   
   !$omp PARALLEL DO private(i, m, atom_2, j, k, x, y, z, r)
   do i = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R

            ! Vs and its derivatives:
            do k = 1, N	! for all electron orbitals (s, px, py, pz):
               M_Vs(k,i,j) = Vs(TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), k, r) ! function below
               M_dVs(k,i,j) = dVs(TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), k, r, M_Vs(k,i,j)) ! function below
            enddo
            
            ! Second derivatives of Vs: 
            if (present(M_d2Vs)) then
               do k = 1, N	! all orbitals
                  M_d2Vs(k,i,j) = d2Vs( TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), k, r, M_Vs(k,i,j), M_dVs(k,i,j) ) ! function below
               enddo ! k = 1, N
            endif ! (present(M_d2Vs)) 
            
            ! Directional cosines:
            if (present(M_cos)) then
               x => Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
               y => Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
               z => Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
               M_cos(i,j,1) = x/r
               M_cos(i,j,2) = y/r
               M_cos(i,j,3) = z/r
            endif ! (present(M_cos))
            
         endif !  (j .GT. 0) 
      enddo ! atom_2 = 1,m 
   enddo ! do i = 1,nat
   !$omp END PARALLEL DO 
   nullify(nat, m, j, r, x, y, z)
end subroutine Construct_M_Vs



subroutine Attract_TB_H3_near(Aij, dHijx_s, dHijy_s, dHijz_s, Scell, NSC, Eelectr_s)
   REAL(8), DIMENSION(:,:), intent(in) :: dHijx_s, dHijy_s, dHijz_s
   REAL(8), DIMENSION(:,:), INTENT(in), target :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis 
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   REAL(8), DIMENSION(:), intent(out)  :: Eelectr_s ! part of the forces

   integer i, j, k, i2, ste, n
   integer, pointer :: m, j1
   real(8), pointer :: Aijij, Aijij1, Aijij2, Aijij3

   n = size(Aij,1)
   Eelectr_s = 0.0d0

   i2 = 0
   ste = 1
   do i = 1, n
       if (i .GE. ste) then
          i2 = i2 + 1
          ste = ste + 4
       endif
       m => Scell(NSC)%Near_neighbor_size(i2)
       do k = 1, m
          j1 => Scell(NSC)%Near_neighbor_list(i2,k) ! this is the list of such close atoms
          j = (j1-1)*4 + 1
          if (j .GT. 0) then
              Aijij => Aij(i,j)
              Aijij1 => Aij(i,j+1)
              Aijij2 => Aij(i,j+2)
              Aijij3 => Aij(i,j+3)
              Eelectr_s(1) = Eelectr_s(1) + (dHijx_s(i,j)*Aijij) + (dHijx_s(i,j+1)*Aijij1) + (dHijx_s(i,j+2)*Aijij2) + (dHijx_s(i,j+3)*Aijij3)
              Eelectr_s(2) = Eelectr_s(2) + (dHijy_s(i,j)*Aijij) + (dHijy_s(i,j+1)*Aijij1) + (dHijy_s(i,j+2)*Aijij2) + (dHijy_s(i,j+3)*Aijij3)
              Eelectr_s(3) = Eelectr_s(3) + (dHijz_s(i,j)*Aijij) + (dHijz_s(i,j+1)*Aijij1) + (dHijz_s(i,j+2)*Aijij2) + (dHijz_s(i,j+3)*Aijij3)
             if (minval(ABS(Eelectr_s(:))) .GE. 1.0d6) then
                write(*,'(a)') 'Trouble in subroutine Attract_TB_H3_near, too large attractive force:'
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j), dHijy_s(i,j), dHijz_s(i,j), Aijij
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+1), dHijy_s(i,j+1), dHijz_s(i,j+1), Aijij1
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+2), dHijy_s(i,j+2), dHijz_s(i,j+2), Aijij2
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+3), dHijy_s(i,j+3), dHijz_s(i,j+3), Aijij3
                write(*,'(e25.16,e25.16,e25.16)') Eelectr_s(:)
             endif
          endif
       enddo
   enddo
   nullify(m, j1, Aijij, Aijij1, Aijij2, Aijij3)
end subroutine Attract_TB_H3_near



subroutine dHamil_tot_r(dHijx, dHijy, dHijz,  d2Hijx, d2Hijy, d2Hijz, TB_Hamil, Scell, NSC, numpar, k, M_Vs, M_dVs, M_d2Vs, M_cos)
   integer, INTENT(IN) :: k ! number of atom
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   REAL(8), DIMENSION(:,:), INTENT(out) :: dHijx, dHijy, dHijz,   d2Hijx, d2Hijy, d2Hijz
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   real(8), dimension(:,:,:), intent(in) :: M_d2Vs ! matrix of functions d2Vs
   real(8), dimension(:,:,:), intent(in) :: M_cos	! matrix of directional cosines
   !-----------------------------------------------
   REAL(8), DIMENSION(:,:), allocatable :: dHijx1, dHijy1, dHijz1,   d2Hijx1, d2Hijy1, d2Hijz1
   integer :: nat, i, N, i4, j4, atom_2, j1, i1, j4j1, i4i1
   integer, pointer :: m, j
   nat = size(Scell(NSC)%MDAtoms)	! number of atoms
   N = size(TB_Hamil(Scell(NSC)%MDatoms(1)%KOA, Scell(NSC)%MDatoms(1)%KOA)%V0) ! that's how many orbitals per atom
   allocate(dHijx1(N,N))
   allocate(dHijy1(N,N))
   allocate(dHijz1(N,N))
   allocate(d2Hijx1(N,N))
   allocate(d2Hijy1(N,N))
   allocate(d2Hijz1(N,N))
   dHijx = 0.0d0
   dHijy = 0.0d0
   dHijz = 0.0d0
   d2Hijx1 = 0.0d0
   d2Hijy1 = 0.0d0
   d2Hijz1 = 0.0d0

   do i = 1, nat	! all atoms
      i4 = (i-1)*N
      m => Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            j4 = (j-1)*N
            call dHamilton_one_r(TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), i, j, k, dHijx1, dHijy1, dHijz1,  d2Hijx1, d2Hijy1, d2Hijz1,  Scell, NSC, atom_2, M_Vs, M_dVs, M_d2Vs, M_cos)
            do j1 = 1,N ! all orbitals
               j4j1 = j4+j1
               do i1 = 1,N ! all orbitals
                  i4i1 = i4+i1
                  dHijx(i4i1, j4j1) = dHijx1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  dHijy(i4i1, j4j1) = dHijy1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  dHijz(i4i1, j4j1) = dHijz1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  d2Hijx(i4i1, j4j1) = d2Hijx1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  d2Hijy(i4i1, j4j1) = d2Hijy1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  d2Hijz(i4i1, j4j1) = d2Hijz1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i
   
   nullify(m, j)
   deallocate(dHijx1,dHijy1,dHijz1, d2Hijx1, d2Hijy1, d2Hijz1)
end subroutine dHamil_tot_r



subroutine dHamil_tot_s(dHijx, dHijy, dHijz, TB_Hamil, Scell, NSC, numpar, k, M_x1, M_xrr, M_Vs, M_dVs)
! construct the whole Hamilton matrix:
! it uses input of Hij, the Hamiltonain block,
! input Es and Ep - the onsite atomic energies
! X - the array of coordinates of all atoms
! nat - the number of atoms
! supce - the size of the supercell, used for periodical boundary conditions
! k - is the atom for which we calculate the forces 
! (with respect to which Rk we take the derivatives, Appendix F of H.Jeschke PhD Thesis)
   integer, INTENT(IN) :: k ! number of atom
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   REAL(8), DIMENSION(:,:), INTENT(out) :: dHijx, dHijy, dHijz
   REAL(8), DIMENSION(:,:,:), INTENT(in) :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in) :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   !-----------------------------------------------
   REAL(8), DIMENSION(4,4) :: dHijx1, dHijy1, dHijz1
   integer i, j1, i1, ki, atom_2
   integer i4, j4, i4i1, j4j1
   integer, pointer :: nat, m, j

   nat => Scell(NSC)%Na
   dHijx = 0.0d0
   dHijy = 0.0d0
   dHijz = 0.0d0

   do i = 1,nat	! all atoms
      !n1 = size(TB_Hamil(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)%V0) ! that's how many orbitals per atom
      i4 = (i-1)*4
      m => Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            j4 = (j-1)*4
            ! Which kinds of atoms are these (which TB parameters to use):
            !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA)
            !call dHamilton_one_s(TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), i, j, k, dHijx, dHijy, dHijz, Scell, NSC, atom_2)
            call dHamilton_one_s(TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), i, j, k, dHijx1, dHijy1, dHijz1, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
            do j1 = 1,4 ! all orbitals
               j4j1 = j4+j1
               do i1 = 1,4 ! all orbitals
                  i4i1 = i4+i1
!                  dHijx((i-1)*4+i1, (j-1)*4+j1) = dHijx(i1,j1) ! construct from the blocks of one-atom Hamiltonian
!                  dHijy((i-1)*4+i1, (j-1)*4+j1) = dHijy(i1,j1) ! construct from the blocks of one-atom Hamiltonian
!                  dHijz((i-1)*4+i1, (j-1)*4+j1) = dHijz(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  dHijx(i4i1, j4j1) = dHijx1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  dHijy(i4i1, j4j1) = dHijy1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
                  dHijz(i4i1, j4j1) = dHijz1(i1,j1) ! construct from the blocks of one-atom Hamiltonian
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i
   nullify(nat, m, j)
end subroutine dHamil_tot_s




subroutine dHamilton_one_r(TB, i, j, k, dHijx, dHijy, dHijz, d2Hijx, d2Hijy, d2Hijz, Scell, NSC, atom_2, M_Vs, M_dVs, M_d2Vs, M_cos)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   type(TB_H_Pettifor), intent(in) :: TB ! all tight binding parameters
   integer, INTENT(IN) :: i, j, k, atom_2
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   REAL(8), DIMENSION(:,:), INTENT(out) :: dHijx, dHijy, dHijz,  d2Hijx, d2Hijy, d2Hijz
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   real(8), dimension(:,:,:), intent(in) :: M_d2Vs ! matrix of functions d2Vs
   real(8), dimension(:,:,:), intent(in) :: M_cos	! matrix of directional cosines
   !----------------------------
   integer(4) ki, kj, i1, j1, k1, N
   REAL(8), DIMENSION(4,4) :: dtsx, dtsy, dtsz	! first derivatives of hopping integrals
   REAL(8), DIMENSION(4,4) :: d2tsx, d2tsy, d2tsz	! second derivatives of hopping integrals
   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      dHijx = 0.0d0   ! diagonals are zeros
      dHijy = 0.0d0   ! diagonals are zeros
      dHijz = 0.0d0   ! diagonals are zeros
      d2Hijx = 0.0d0
      d2Hijy = 0.0d0
      d2Hijz = 0.0d0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_r(dtsx, dtsy, dtsz, d2tsx, d2tsy, d2tsz, i, j, k, TB, Scell, NSC, atom_2, M_Vs, M_dVs, M_d2Vs, M_cos)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      N = size(TB%V0) ! that's how many orbitals per atom
      do ki = 1,N
         dHijx(:,ki) = dtsx(ki,:)   ! derivatie of the hopping Integrals by x
         dHijy(:,ki) = dtsy(ki,:)   ! derivatie of the hopping Integrals by y
         dHijz(:,ki) = dtsz(ki,:)   ! derivatie of the hopping Integrals by z
         d2Hijx(:,ki) = d2tsx(ki,:)  ! second derivatie of the hopping Integrals by x
         d2Hijy(:,ki) = d2tsy(ki,:)  ! second derivatie of the hopping Integrals by y
         d2Hijz(:,ki) = d2tsz(ki,:)  ! second derivatie of the hopping Integrals by z
      enddo  ! ki
   endif
end subroutine dHamilton_one_r




subroutine dHamilton_one_s(TB, i, j, k, dHijx, dHijy, dHijz, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   type(TB_H_Pettifor), intent(in) :: TB ! all tight binding parameters
   integer, INTENT(IN) :: i, j, k, atom_2
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   REAL(8), DIMENSION(:,:), INTENT(out) :: dHijx, dHijy, dHijz
   REAL(8), DIMENSION(:,:,:), INTENT(in) :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in) :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   !----------------------------
   integer(4) ki, kj, i1, j1, k1
   REAL(8), DIMENSION(4,4) :: dtsx, dtsy, dtsz ! hopping integrals
   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      dHijx = 0.0d0   ! diagonals are zeros
      dHijy = 0.0d0   ! diagonals are zeros
      dHijz = 0.0d0   ! diagonals are zeros
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_s(dtsx, dtsy, dtsz, i, j, k, TB, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1,4
         dHijx(:,ki) = dtsx(ki,:)   ! Hopping Integrals
         dHijy(:,ki) = dtsy(ki,:)   ! Hopping Integrals
         dHijz(:,ki) = dtsz(ki,:)   ! Hopping Integrals
!         
!          do kj = 1,4
!              dHijx(kj,ki) = dtsx(ki,kj)   ! Hopping Integrals
!              dHijy(kj,ki) = dtsy(ki,kj)   ! Hopping Integrals
!              dHijz(kj,ki) = dtsz(ki,kj)   ! Hopping Integrals
!              if (ABS(dHijx(kj,ki)) .LE. 1d-14) then
!                 dHijx(kj,ki) = 0.0d0
!              endif
!              if (ABS(dHijy(kj,ki)) .LE. 1d-14) then
!                 dHijy(kj,ki) = 0.0d0
!              endif
!              if (ABS(dHijz(kj,ki)) .LE. 1d-14) then
!                 dHijz(kj,ki) = 0.0d0
!              endif
!          enddo ! kj
      enddo  ! ki
   endif
end subroutine dHamilton_one_s



subroutine dHopping_r(dtsx, dtsy, dtsz,  d2tsx, d2tsy, d2tsz,  i, j, k, TB, Scell, NSC, atom_2, M_Vs, M_dVs, M_d2Vs, M_cos)
! subroutine making the derivatives of the hopping integrals
! dtsz, dtsy, dtsx --- the derivatives Hopping Integrals in Z,Y,X axis
! Xco --- atomic coordinates of all the atoms, 
! supce --- the supercell size
! i, j - a pair of atoms, with respect to which we calculate a force acting on 
! the atom "k".
   type(TB_H_Pettifor), intent(in) :: TB ! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: dtsx, dtsy, dtsz,  d2tsx, d2tsy, d2tsz
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   INTEGER(4), intent(in) :: i,j,k, atom_2
   real(8), dimension(:,:,:), intent(in), target :: M_Vs	! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in), target :: M_dVs	! matrix of functions dVs
   real(8), dimension(:,:,:), intent(in), target :: M_d2Vs	! matrix of functions d2Vs
   real(8), dimension(:,:,:), intent(in), target :: M_cos	! matrix of directional cosines
   !---------------------------------------
   real(8) :: x0, y0, z0, sx, sy, sz	! relative distances (projections) between the atoms 
   real(8) r1, temp, temp1, temp2, temp3, z_temp !, d2Vs1, d2Vs2, d2Vs3, d2Vs4
   INTEGER(4) :: ik, i1, j1, k1, dik, djk, dij
   REAL(8), DIMENSION(3) :: dlds, dmds, dnds, zb, d2lds, d2mds, d2nds
   real(8) drdrx, drdry, drdrz, b, d2rdr2x, d2rdr2y, d2rdr2z
   real(8) :: drdrx2, drdry2, drdrz2
   real(8), pointer :: dVs1, dVs2, dVs3, dVs4, Vs1, Vs2, Vs3, Vs4, d2Vs1, d2Vs2, d2Vs3, d2Vs4
   real(8), pointer :: x,y,z,r, xrr, yrr, zrr,  xr, yr, zr
   REAL(8), DIMENSION(:), pointer :: x1

   dik = Kronecker_delta(i,k)	! module "Atomic_tools"
   djk = Kronecker_delta(j,k)	! module "Atomic_tools"

   if (ABS(dik - djk) > 1.0d-12) then ! only then it is non-zero
      x => Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
      y => Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
      z => Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
      r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R

      ! Directional cosines:
      xr => M_cos(i,j,1)
      yr => M_cos(i,j,2)
      zr => M_cos(i,j,3)
!       xr = x/r
!       yr = y/r
!       zr = z/r
      
      ! Vksi(r) functions and their derivatives:
      Vs1 => M_Vs(1,i,j)	! V s s sigma
      Vs2 => M_Vs(2,i,j)	! V s p sigma
      Vs3 => M_Vs(3,i,j)	! V p p sigma
      Vs4 => M_Vs(4,i,j)	! V p p pi
      dVs1 => M_dVs(1,i,j)
      dVs2 => M_dVs(2,i,j)
      dVs3 => M_dVs(3,i,j)
      dVs4 => M_dVs(4,i,j)
      d2Vs1 => M_d2Vs(1,i,j)
      d2Vs2 => M_d2Vs(2,i,j)
      d2Vs3 => M_d2Vs(3,i,j)
      d2Vs4 => M_d2Vs(4,i,j)
!       d2Vs1 = d2Vs(TB, 1, r, Vs1, dVs1)	! second derevative of the attractive potential contribution
!       d2Vs2 = d2Vs(TB, 2, r, Vs2, dVs2)	! second derevative of the attractive potential contribution
!       d2Vs3 = d2Vs(TB, 3, r, Vs3, dVs3)	! second derevative of the attractive potential contribution
!       d2Vs4 = d2Vs(TB, 4, r, Vs4, dVs4)	! second derevative of the attractive potential contribution
      
      ! Derivative d r_{ij} / d r_{k,alpha}
      drdrx = drij_drka(i, j, k, x, r)	! module "TB_Koster_Slater"
      drdry = drij_drka(i, j, k, y, r)	! module "TB_Koster_Slater"
      drdrz = drij_drka(i, j, k, z, r)	! module "TB_Koster_Slater"
      
      ! Squares of derivatives: 
      drdrx2 = drdrx*drdrx 
      drdry2 = drdry*drdry
      drdrz2 = drdrz*drdrz
      
      ! Second derivatives d2 r_{ij} / d r2_{k,alpha} = d d_{alpha} / d r_{k,alpha}:
      d2rdr2x = ddija_drkb(i, j, k, 1, 1, x, x, r)	! module "TB_Koster_Slater"
      d2rdr2y = ddija_drkb(i, j, k, 2, 2, y, y, r)	! module "TB_Koster_Slater"
      d2rdr2z = ddija_drkb(i, j, k, 3, 3, z, z, r)	! module "TB_Koster_Slater"
      
      ! Derivatives of cosine directions:
      dlds(1) = ddija_drkb(i, j, k, 1, 1, x, x, r)	! module "TB_Koster_Slater"
      dlds(2) = ddija_drkb(i, j, k, 1, 2, x, y, r)	! module "TB_Koster_Slater"
      dlds(3) = ddija_drkb(i, j, k, 1, 3, x, z, r)	! module "TB_Koster_Slater"
      dmds(1) = ddija_drkb(i, j, k, 2, 1, y, x, r)	! module "TB_Koster_Slater"
      dmds(2) = ddija_drkb(i, j, k, 2, 2, y, y, r)	! module "TB_Koster_Slater"
      dmds(3) = ddija_drkb(i, j, k, 2, 3, y, z, r)	! module "TB_Koster_Slater"
      dnds(1) = ddija_drkb(i, j, k, 3, 1, z, x, r)	! module "TB_Koster_Slater"
      dnds(2) = ddija_drkb(i, j, k, 3, 2, z, y, r)	! module "TB_Koster_Slater"
      dnds(3) = ddija_drkb(i, j, k, 3, 3, z, z, r)	! module "TB_Koster_Slater"
      
      !Second derivatives of cosine directions:
      d2lds(1) = d2dija_drkb2(i, j, k, 1, 1, x, x, r)	! module "TB_Koster_Slater"
      d2lds(2) = d2dija_drkb2(i, j, k, 1, 2, x, y, r)	! module "TB_Koster_Slater"
      d2lds(3) = d2dija_drkb2(i, j, k, 1, 3, x, z, r)	! module "TB_Koster_Slater"
      d2mds(1) = d2dija_drkb2(i, j, k, 2, 1, y, x, r)	! module "TB_Koster_Slater"
      d2mds(2) = d2dija_drkb2(i, j, k, 2, 2, y, y, r)	! module "TB_Koster_Slater"
      d2mds(3) = d2dija_drkb2(i, j, k, 2, 3, y, z, r)	! module "TB_Koster_Slater"
      d2nds(1) = d2dija_drkb2(i, j, k, 3, 1, z, x, r)	! module "TB_Koster_Slater"
      d2nds(2) = d2dija_drkb2(i, j, k, 3, 2, z, y, r)	! module "TB_Koster_Slater"
      d2nds(3) = d2dija_drkb2(i, j, k, 3, 3, z, z, r)	! module "TB_Koster_Slater"
      
      ! Derivatives of the overlap integrals:
      ! First:
      dtsx(1,1) = dVs1*drdrx ! ss
      dtsy(1,1) = dVs1*drdry ! ss
      dtsz(1,1) = dVs1*drdrz ! ss
      ! Second:
      d2tsx(1,1) = d2Vs1*drdrx2 + dVs1*d2rdr2x
      d2tsy(1,1) = d2Vs1*drdry2 + dVs1*d2rdr2y
      d2tsz(1,1) = d2Vs1*drdrz2 + dVs1*d2rdr2z
      
      ! First:
      temp = xr*dVs2
      dtsx(2,1) = -(dlds(1)*Vs2 + temp*drdrx) ! -spx
      dtsy(2,1) = -(dlds(2)*Vs2 + temp*drdry) ! -spx
      dtsz(2,1) = -(dlds(3)*Vs2 + temp*drdrz) ! -spx
      ! Second:
      d2tsx(2,1) = -( d2lds(1)*Vs2 + dVs2*(2.0d0*dlds(1)*drdrx + d2rdr2x*xr) + d2Vs2*drdrx2*xr )
      d2tsy(2,1) = -( d2lds(2)*Vs2 + dVs2*(2.0d0*dlds(2)*drdry + d2rdr2y*xr) + d2Vs2*drdry2*xr )
      d2tsz(2,1) = -( d2lds(3)*Vs2 + dVs2*(2.0d0*dlds(3)*drdrz + d2rdr2z*xr) + d2Vs2*drdrz2*xr )
      
      ! First:
      dtsx(1,2) = -dtsx(2,1)
      dtsy(1,2) = -dtsy(2,1)
      dtsz(1,2) = -dtsz(2,1)
      ! Second:
      d2tsx(1,2) = -d2tsx(2,1)
      d2tsy(1,2) = -d2tsy(2,1)
      d2tsz(1,2) = -d2tsz(2,1)

      ! First:
      temp = yr*dVs2
      dtsx(3,1) = -(dmds(1)*Vs2 + temp*drdrx) ! -spy
      dtsy(3,1) = -(dmds(2)*Vs2 + temp*drdry) ! -spy
      dtsz(3,1) = -(dmds(3)*Vs2 + temp*drdrz) ! -spy
      ! Second:
      d2tsx(3,1) = -( d2mds(1)*Vs2 + dVs2*(2.0d0*dmds(1)*drdrx + d2rdr2x*yr) + d2Vs2*drdrx2*yr )
      d2tsy(3,1) = -( d2mds(2)*Vs2 + dVs2*(2.0d0*dmds(2)*drdry + d2rdr2y*yr) + d2Vs2*drdry2*yr )
      d2tsz(3,1) = -( d2mds(3)*Vs2 + dVs2*(2.0d0*dmds(3)*drdrz + d2rdr2z*yr) + d2Vs2*drdrz2*yr )
      
      ! First:
      dtsx(1,3) = -dtsx(3,1)
      dtsy(1,3) = -dtsy(3,1)
      dtsz(1,3) = -dtsz(3,1)
      ! Second:
      d2tsx(1,3) = -d2tsx(3,1)
      d2tsy(1,3) = -d2tsy(3,1)
      d2tsz(1,3) = -d2tsz(3,1)
      
      ! First:
      temp = zr*dVs2
      dtsx(4,1) = -(dnds(1)*Vs2 + temp*drdrx)  ! -spz
      dtsy(4,1) = -(dnds(2)*Vs2 + temp*drdry)  ! -spz
      dtsz(4,1) = -(dnds(3)*Vs2 + temp*drdrz)  ! -spz
      ! Second:
      d2tsx(4,1) = -( d2nds(1)*Vs2 + dVs2*(2.0d0*dnds(1)*drdrx + d2rdr2x*zr) + d2Vs2*drdrx2*zr )
      d2tsy(4,1) = -( d2nds(2)*Vs2 + dVs2*(2.0d0*dnds(2)*drdry + d2rdr2y*zr) + d2Vs2*drdry2*zr )
      d2tsz(4,1) = -( d2nds(3)*Vs2 + dVs2*(2.0d0*dnds(3)*drdrz + d2rdr2z*zr) + d2Vs2*drdrz2*zr )

      ! First:
      dtsx(1,4) = -dtsx(4,1)
      dtsy(1,4) = -dtsy(4,1)
      dtsz(1,4) = -dtsz(4,1)
      ! Second:
      d2tsx(1,4) = -d2tsx(4,1)
      d2tsy(1,4) = -d2tsy(4,1)
      d2tsz(1,4) = -d2tsz(4,1)
      
      ! First:
      temp = 2.0d0*xr*(Vs3 - Vs4)
      temp1 = (xr*xr*(dVs3 - dVs4) + dVs4)
      dtsx(2,2) = temp*dlds(1) + temp1*drdrx ! pxpx
      dtsy(2,2) = temp*dlds(2) + temp1*drdry ! pxpx
      dtsz(2,2) = temp*dlds(3) + temp1*drdrz ! pxpx
      ! Second:
      temp = Vs3 - Vs4
      temp3 =  (d2Vs3 - d2Vs4)*xr*xr + d2Vs4
      d2tsx(2,2) = 2.0d0*temp*(dlds(1)*dlds(1) + xr*d2lds(1)) + (4.0d0*(dVs3 - dVs4)*xr*dlds(1) + temp3*drdrx)*drdrx + temp1*d2rdr2x
      d2tsy(2,2) = 2.0d0*temp*(dlds(2)*dlds(2) + xr*d2lds(2)) + (4.0d0*(dVs3 - dVs4)*xr*dlds(2) + temp3*drdry)*drdry + temp1*d2rdr2y
      d2tsz(2,2) = 2.0d0*temp*(dlds(3)*dlds(3) + xr*d2lds(3)) + (4.0d0*(dVs3 - dVs4)*xr*dlds(3) + temp3*drdrz)*drdrz + temp1*d2rdr2z

      ! First:
      !temp = Vs3 - Vs4
      temp1 = xr*yr*(dVs3 - dVs4)
      dtsx(3,2) = (xr*dmds(1) + yr*dlds(1))*temp + temp1*drdrx ! pxpy
      dtsy(3,2) = (xr*dmds(2) + yr*dlds(2))*temp + temp1*drdry ! pxpy
      dtsz(3,2) = (xr*dmds(3) + yr*dlds(3))*temp + temp1*drdrz ! pxpy
      ! Second:
      d2tsx(3,2) = temp*(d2lds(1)*yr + d2mds(1)*xr + 2.0d0*(dlds(1)*dmds(1))) + (dVs3 - dVs4)*(2.0d0*(dlds(1)*yr + dmds(1)*xr)*drdrx + xr*yr*d2rdr2x) + (d2Vs3 - d2Vs4)*xr*yr*drdrx2
      d2tsy(3,2) = temp*(d2lds(2)*yr + d2mds(2)*xr + 2.0d0*(dlds(2)*dmds(2))) + (dVs3 - dVs4)*(2.0d0*(dlds(2)*yr + dmds(2)*xr)*drdry + xr*yr*d2rdr2y) + (d2Vs3 - d2Vs4)*xr*yr*drdry2
      d2tsz(3,2) = temp*(d2lds(3)*yr + d2mds(3)*xr + 2.0d0*(dlds(3)*dmds(3))) + (dVs3 - dVs4)*(2.0d0*(dlds(3)*yr + dmds(3)*xr)*drdrz + xr*yr*d2rdr2z) + (d2Vs3 - d2Vs4)*xr*yr*drdrz2

      ! First:
      dtsx(2,3) = dtsx(3,2)
      dtsy(2,3) = dtsy(3,2)
      dtsz(2,3) = dtsz(3,2)
      ! Second:
      d2tsx(2,3) = d2tsx(3,2)
      d2tsy(2,3) = d2tsy(3,2)
      d2tsz(2,3) = d2tsz(3,2)

      ! First:
      temp1 = xr*zr*(dVs3 - dVs4)
      dtsx(4,2) = (xr*dnds(1) + zr*dlds(1))*temp + temp1*drdrx ! pxpz
      dtsy(4,2) = (xr*dnds(2) + zr*dlds(2))*temp + temp1*drdry ! pxpz
      dtsz(4,2) = (xr*dnds(3) + zr*dlds(3))*temp + temp1*drdrz ! pxpz
      ! Second:
      d2tsx(4,2) = temp*(d2lds(1)*zr + d2nds(1)*xr + 2.0d0*(dlds(1)*dnds(1))) + (dVs3 - dVs4)*(2.0d0*(dlds(1)*zr + dnds(1)*xr)*drdrx + xr*zr*d2rdr2x) + (d2Vs3 - d2Vs4)*xr*zr*drdrx2
      d2tsy(4,2) = temp*(d2lds(2)*zr + d2nds(2)*xr + 2.0d0*(dlds(2)*dnds(2))) + (dVs3 - dVs4)*(2.0d0*(dlds(2)*zr + dnds(2)*xr)*drdry + xr*zr*d2rdr2y) + (d2Vs3 - d2Vs4)*xr*zr*drdry2
      d2tsz(4,2) = temp*(d2lds(3)*zr + d2nds(3)*xr + 2.0d0*(dlds(3)*dnds(3))) + (dVs3 - dVs4)*(2.0d0*(dlds(3)*zr + dnds(3)*xr)*drdrz + xr*zr*d2rdr2z) + (d2Vs3 - d2Vs4)*xr*zr*drdrz2

      ! First:
      dtsx(2,4) = dtsx(4,2)
      dtsy(2,4) = dtsy(4,2)
      dtsz(2,4) = dtsz(4,2)
      ! Second:
      d2tsx(2,4) = d2tsx(4,2)
      d2tsy(2,4) = d2tsy(4,2)
      d2tsz(2,4) = d2tsz(4,2)

      ! First:
      temp2 = 2.0d0*yr*(Vs3 - Vs4)
      temp1 = yr*yr*(dVs3 - dVs4) + dVs4
      dtsx(3,3) = temp2*dmds(1) + temp1*drdrx ! pypy
      dtsy(3,3) = temp2*dmds(2) + temp1*drdry ! pypy
      dtsz(3,3) = temp2*dmds(3) + temp1*drdrz ! pypy
      ! Second:
      temp3 =  (d2Vs3 - d2Vs4)*yr*yr + d2Vs4
      !temp = Vs3 - Vs4
      d2tsx(3,3) = 2.0d0*temp*(dmds(1)*dmds(1) + yr*d2mds(1)) + (4.0d0*(dVs3 - dVs4)*yr*dmds(1) + temp3*drdrx)*drdrx + temp1*d2rdr2x
      d2tsy(3,3) = 2.0d0*temp*(dmds(2)*dmds(2) + yr*d2mds(2)) + (4.0d0*(dVs3 - dVs4)*yr*dmds(2) + temp3*drdry)*drdry + temp1*d2rdr2y
      d2tsz(3,3) = 2.0d0*temp*(dmds(3)*dmds(3) + yr*d2mds(3)) + (4.0d0*(dVs3 - dVs4)*yr*dmds(3) + temp3*drdrz)*drdrz + temp1*d2rdr2z

      ! First:
      temp1 = yr*zr*(dVs3 - dVs4)
      dtsx(4,3) = (yr*dnds(1) + zr*dmds(1))*temp + temp1*drdrx ! pypz
      dtsy(4,3) = (yr*dnds(2) + zr*dmds(2))*temp + temp1*drdry ! pypz
      dtsz(4,3) = (yr*dnds(3) + zr*dmds(3))*temp + temp1*drdrz ! pypz
      ! Second:
      d2tsx(4,3) = temp*(d2mds(1)*zr + d2nds(1)*yr + 2.0d0*(dmds(1)*dnds(1))) + (dVs3 - dVs4)*(2.0d0*(dmds(1)*zr + dnds(1)*yr)*drdrx + yr*zr*d2rdr2x) + (d2Vs3 - d2Vs4)*yr*zr*drdrx2
      d2tsy(4,3) = temp*(d2mds(2)*zr + d2nds(2)*yr + 2.0d0*(dmds(2)*dnds(2))) + (dVs3 - dVs4)*(2.0d0*(dmds(2)*zr + dnds(2)*yr)*drdry + yr*zr*d2rdr2y) + (d2Vs3 - d2Vs4)*yr*zr*drdry2
      d2tsz(4,3) = temp*(d2mds(3)*zr + d2nds(3)*yr + 2.0d0*(dmds(3)*dnds(3))) + (dVs3 - dVs4)*(2.0d0*(dmds(3)*zr + dnds(3)*yr)*drdrz + yr*zr*d2rdr2z) + (d2Vs3 - d2Vs4)*yr*zr*drdrz2

      ! First:
      dtsx(3,4) = dtsx(4,3)
      dtsy(3,4) = dtsy(4,3)
      dtsz(3,4) = dtsz(4,3)
      ! Second:
      d2tsx(3,4) = d2tsx(4,3)
      d2tsy(3,4) = d2tsy(4,3)
      d2tsz(3,4) = d2tsz(4,3)

      ! First:
      temp1 = zr*zr*(dVs3 - dVs4) + dVs4
      z_temp = 2.0d0*zr*temp
      dtsx(4,4) = z_temp*dnds(1) + temp1*drdrx ! pzpz
      dtsy(4,4) = z_temp*dnds(2) + temp1*drdry ! pzpz
      dtsz(4,4) = z_temp*dnds(3) + temp1*drdrz ! pzpz
      ! Second:
      temp3 =  (d2Vs3 - d2Vs4)*zr*zr + d2Vs4
      !temp = Vs3 - Vs4
      d2tsx(4,4) = 2.0d0*temp*(dnds(1)*dnds(1) + zr*d2nds(1)) + (4.0d0*(dVs3 - dVs4)*zr*dnds(1) + temp3*drdrx)*drdrx + temp1*d2rdr2x
      d2tsy(4,4) = 2.0d0*temp*(dnds(2)*dnds(2) + zr*d2nds(2)) + (4.0d0*(dVs3 - dVs4)*zr*dnds(2) + temp3*drdry)*drdry + temp1*d2rdr2y
      d2tsz(4,4) = 2.0d0*temp*(dnds(3)*dnds(3) + zr*d2nds(3)) + (4.0d0*(dVs3 - dVs4)*zr*dnds(3) + temp3*drdrz)*drdrz + temp1*d2rdr2z

      nullify(d2Vs1, d2Vs2, d2Vs3, d2Vs4, dVs1, dVs2, dVs3, dVs4, Vs1, Vs2, Vs3, Vs4, x, y, z, r, x1, xrr, yrr, zrr, xr, yr, zr)
   else ! it's all zeros:
      ! First:
      dtsx = 0.0d0
      dtsy = 0.0d0
      dtsz = 0.0d0
      ! Second:
      d2tsx = 0.0d0
      d2tsy = 0.0d0
      d2tsz = 0.0d0
   endif
end subroutine dHopping_r



subroutine dHopping_s(dtsx, dtsy, dtsz, i, j, k, TB, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
! subroutine making the derivatives of the hopping integrals
! dtsz, dtsy, dtsx --- the derivatives Hopping Integrals in Z,Y,X axis
! Xco --- atomic coordinates of all the atoms, 
! supce --- the supercell size
! i, j - a pair of atoms, with respect to which we calculate a force acting on 
! the atom "k".
   type(TB_H_Pettifor), intent(in) :: TB ! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: dtsx, dtsy, dtsz
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   INTEGER(4), intent(in) :: i,j,k, atom_2
   REAL(8), DIMENSION(:,:,:), INTENT(in), target :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in), target :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), intent(in), target :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in), target :: M_dVs ! matrix of functions dVs
   !---------------------------------------
   real(8) :: x0, y0, z0, sx, sy, sz	! relative distances (projections) between the atoms 
   real(8) r1, temp, temp1, temp2, z_temp
   INTEGER(4) :: ik, i1, j1, k1, dik, djk, dij, i_test, j_test
   REAL(8), DIMENSION(3) :: dlds, dmds, dnds, zb
   real(8) drdsx, drdsy, drdsz, b, xr, yr, zr
   real(8), pointer :: dVs1, dVs2, dVs3, dVs4, Vs1, Vs2, Vs3, Vs4
   real(8), pointer :: x,y,z,r, xrr, yrr, zrr
   REAL(8), DIMENSION(:), pointer :: x1

   if (i == k) then
      dik = 1
   else
      dik = 0
   endif
   if (k == j) then
      djk = 1
   else
      djk = 0
   endif

   if ((dik - djk) /= 0) then ! only then it is non-zero
      x => Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
      y => Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
      z => Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
      r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R

      b = dble(dik - djk)/r
      
      x1 => M_x1(:,i,j)
!       x1(1) => M_x1(1,i,j) !x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
!       x1(2) => M_x1(2,i,j) !x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
!       x1(3) => M_x1(3,i,j) !x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

      drdsx = b*x1(1)
      drdsy = b*x1(2)
      drdsz = b*x1(3)
      
!        print*, 'drdsx', drdsx, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 1, .true.)
!        print*, 'drdsy', drdsy, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 2, .true.)
!        print*, 'drdsz', drdsz, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 3, .true.)
!        write(*,'(a,i3,i3,i3,es,es,es,es)') 'x', i, j, k, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 1, with_delta=.true.), drdsx, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 1), x1(1)/r
!        write(*,'(a,i3,i3,i3,es,es,es,es)') 'y', i, j, k, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 2, with_delta=.true.), drdsy, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 2), x1(2)/r
!        write(*,'(a,i3,i3,i3,es,es,es,es)') 'z', i, j, k, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 3, with_delta=.true.), drdsz, drij_dska(i, j, k, x, y, z, r, Scell(NSC)%supce, 3), x1(3)/r
      
      xrr => M_xrr(1,i,j) !x/(r*r)
      yrr => M_xrr(2,i,j) !y/(r*r)
      zrr => M_xrr(3,i,j) !z/(r*r)
      
      ! CORRECT:
      dlds(1) = b*Scell(NSC)%supce(1,1) - xrr*drdsx
      dlds(2) = b*Scell(NSC)%supce(2,1) - xrr*drdsy
      dlds(3) = b*Scell(NSC)%supce(3,1) - xrr*drdsz
      dmds(1) = b*Scell(NSC)%supce(1,2) - yrr*drdsx
      dmds(2) = b*Scell(NSC)%supce(2,2) - yrr*drdsy
      dmds(3) = b*Scell(NSC)%supce(3,2) - yrr*drdsz
      dnds(1) = b*Scell(NSC)%supce(1,3) - zrr*drdsx
      dnds(2) = b*Scell(NSC)%supce(2,3) - zrr*drdsy
      dnds(3) = b*Scell(NSC)%supce(3,3) - zrr*drdsz
      ! NOT CORRECT:
!       dlds(1) = b*Scell(NSC)%supce(1,1) - xrr*drdsx
!       dlds(2) = b*Scell(NSC)%supce(1,2) - xrr*drdsy
!       dlds(3) = b*Scell(NSC)%supce(1,3) - xrr*drdsz
!       dmds(1) = b*Scell(NSC)%supce(2,1) - yrr*drdsx
!       dmds(2) = b*Scell(NSC)%supce(2,2) - yrr*drdsy
!       dmds(3) = b*Scell(NSC)%supce(2,3) - yrr*drdsz
!       dnds(1) = b*Scell(NSC)%supce(3,1) - zrr*drdsx
!       dnds(2) = b*Scell(NSC)%supce(3,2) - zrr*drdsy
!       dnds(3) = b*Scell(NSC)%supce(3,3) - zrr*drdsz
 
!     print*, dlds(1), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 1, 1, drdsx) 
!     print*, dlds(2), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 1, 2, drdsy) 
!     print*, dlds(3), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 1, 3, drdsz) 
!     print*, dmds(1), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 2, 1, drdsx) 
!     print*, dmds(2), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 2, 2, drdsy) 
!     print*, dmds(3), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 2, 3, drdsz) 
!     print*, dnds(1), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 3, 1, drdsx) 
!     print*, dnds(2), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 3, 2, drdsy) 
!     print*, dnds(3), ddija_dskb_kd(i, j, k, x, y, z, r, Scell(NSC)%supce, 3, 3, drdsz) 
!     pause 'TEST dnds'
 
 
      ! All of the hopping integrals ts for each pair of atoms i and j,
      ! according to the Eqs.(2.23)-(2.29) from PhD Thesis of H.Jeschke, Page 38
      xr = (x/r)
      yr = (y/r)
      zr = (z/r)
      dVs1 => M_dVs(1,i,j) !dVs(TB,1,r)
      dVs2 => M_dVs(2,i,j) !dVs(TB,2,r)
      dVs3 => M_dVs(3,i,j) !dVs(TB,3,r)
      dVs4 => M_dVs(4,i,j) !dVs(TB,4,r)
      Vs1 => M_Vs(1,i,j) !Vs(TB,1,r)
      Vs2 => M_Vs(2,i,j) !Vs(TB,2,r)
      Vs3 => M_Vs(3,i,j) !Vs(TB,3,r)
      Vs4 => M_Vs(4,i,j) !Vs(TB,4,r)

!       ! testing:
!       call d_KS_sp3(M_Vs(:,i,j), M_dVs(:,i,j)*drdsx, xr, yr, zr, dlds(1), dmds(1), dnds(1), dtsx)	! module "TB_Koster_Slater"
!       call d_KS_sp3(M_Vs(:,i,j), M_dVs(:,i,j)*drdsy, xr, yr, zr, dlds(2), dmds(2), dnds(2), dtsy)	! module "TB_Koster_Slater"
!       call d_KS_sp3(M_Vs(:,i,j), M_dVs(:,i,j)*drdsz, xr, yr, zr, dlds(3), dmds(3), dnds(3), dtsz)	! module "TB_Koster_Slater"
!       write(*,*) 'dtsx-1'
!       do i_test = 1,4
!          do j_test = 1,4
!             print*, 'x', i_test, j_test, dtsx(i_test,j_test)
!             print*, 'z', i_test, j_test, dtsy(i_test,j_test)
!             print*, 'y', i_test, j_test, dtsz(i_test,j_test)
!          enddo
!       enddo

      dtsx(1,1) = dVs1*drdsx ! ss
      dtsy(1,1) = dVs1*drdsy ! ss
      dtsz(1,1) = dVs1*drdsz ! ss

      temp = xr*dVs2
      dtsx(2,1) = -(dlds(1)*Vs2 + temp*drdsx) ! -spx
      dtsy(2,1) = -(dlds(2)*Vs2 + temp*drdsy) ! -spx
      dtsz(2,1) = -(dlds(3)*Vs2 + temp*drdsz) ! -spx

      dtsx(1,2) = -dtsx(2,1)
      dtsy(1,2) = -dtsy(2,1)
      dtsz(1,2) = -dtsz(2,1)

      temp = yr*dVs2
      dtsx(3,1) = -(dmds(1)*Vs2 + temp*drdsx)  ! -spy
      dtsy(3,1) = -(dmds(2)*Vs2 + temp*drdsy) ! -spy
      dtsz(3,1) = -(dmds(3)*Vs2 + temp*drdsz) ! -spy

      dtsx(1,3) = -dtsx(3,1)
      dtsy(1,3) = -dtsy(3,1)
      dtsz(1,3) = -dtsz(3,1)

      temp = zr*dVs2
      dtsx(4,1) = -(dnds(1)*Vs2 + temp*drdsx)  ! -spz
      dtsy(4,1) = -(dnds(2)*Vs2 + temp*drdsy)  ! -spz
      dtsz(4,1) = -(dnds(3)*Vs2 + temp*drdsz)  ! -spz

      dtsx(1,4) = -dtsx(4,1)
      dtsy(1,4) = -dtsy(4,1)
      dtsz(1,4) = -dtsz(4,1)

      temp = 2.0d0*xr*(Vs3 - Vs4)
      temp1 = (xr*xr*(dVs3 - dVs4) + dVs4)
      dtsx(2,2) = temp*dlds(1) + temp1*drdsx ! pxpx
      dtsy(2,2) = temp*dlds(2) + temp1*drdsy ! pxpx
      dtsz(2,2) = temp*dlds(3) + temp1*drdsz ! pxpx

      temp = Vs3 - Vs4
      temp1 = xr*yr*(dVs3 - dVs4)
      dtsx(3,2) = (xr*dmds(1) + yr*dlds(1))*temp + temp1*drdsx ! pxpy
      dtsy(3,2) = (xr*dmds(2) + yr*dlds(2))*temp + temp1*drdsy ! pxpy
      dtsz(3,2) = (xr*dmds(3) + yr*dlds(3))*temp + temp1*drdsz ! pxpy

      dtsx(2,3) = dtsx(3,2)
      dtsy(2,3) = dtsy(3,2)
      dtsz(2,3) = dtsz(3,2)

      temp1 = xr*zr*(dVs3 - dVs4)
      dtsx(4,2) = (xr*dnds(1) + zr*dlds(1))*temp + temp1*drdsx ! pxpz
      dtsy(4,2) = (xr*dnds(2) + zr*dlds(2))*temp + temp1*drdsy ! pxpz
      dtsz(4,2) = (xr*dnds(3) + zr*dlds(3))*temp + temp1*drdsz ! pxpz

      dtsx(2,4) = dtsx(4,2)
      dtsy(2,4) = dtsy(4,2)
      dtsz(2,4) = dtsz(4,2)

      temp2 = 2.0d0*yr*(Vs3 - Vs4)
      temp1 = yr*yr*(dVs3 - dVs4) + dVs4
      dtsx(3,3) = temp2*dmds(1) + temp1*drdsx ! pypy
      dtsy(3,3) = temp2*dmds(2) + temp1*drdsy ! pypy
      dtsz(3,3) = temp2*dmds(3) + temp1*drdsz ! pypy

      temp1 = yr*zr*(dVs3 - dVs4)
      dtsx(4,3) = (yr*dnds(1) + zr*dmds(1))*temp + temp1*drdsx ! pypz
      dtsy(4,3) = (yr*dnds(2) + zr*dmds(2))*temp + temp1*drdsy ! pypz
      dtsz(4,3) = (yr*dnds(3) + zr*dmds(3))*temp + temp1*drdsz ! pypz

      dtsx(3,4) = dtsx(4,3)
      dtsy(3,4) = dtsy(4,3)
      dtsz(3,4) = dtsz(4,3)

      temp1 = zr*zr*(dVs3 - dVs4) + dVs4
      z_temp = 2.0d0*zr*temp
      dtsx(4,4) = z_temp*dnds(1) + temp1*drdsx ! pzpz
      dtsy(4,4) = z_temp*dnds(2) + temp1*drdsy ! pzpz
      dtsz(4,4) = z_temp*dnds(3) + temp1*drdsz ! pzpz
      !dtsx(4,4) = 2.0d0*zr*dnds(1)*temp + temp1*drdsx ! pzpz
      !dtsy(4,4) = 2.0d0*zr*dnds(2)*temp + temp1*drdsy ! pzpz
      !dtsz(4,4) = 2.0d0*zr*dnds(3)*temp + temp1*drdsz ! pzpz
      
!       write(*,*) 'dtsx-2'
!       do i_test = 1,4
!          do j_test = 1,4
!             print*, 'x', i_test, j_test, dtsx(i_test,j_test)
!             print*, 'z', i_test, j_test, dtsy(i_test,j_test)
!             print*, 'y', i_test, j_test, dtsz(i_test,j_test)
!          enddo
!       enddo
!       pause 'dHopping_s'
      
      nullify(dVs1, dVs2, dVs3, dVs4, Vs1, Vs2, Vs3, Vs4, x, y, z, r, x1, xrr, yrr, zrr)
   else ! it's all zeros:
      dtsx = 0.0d0
      dtsy = 0.0d0
      dtsz = 0.0d0
   endif
end subroutine dHopping_s



subroutine Attract_TB_Forces_Press(TB_Hamil, atoms, Scell, NSC, numpar, Aij)
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   REAL(8), DIMENSION(:,:), INTENT(in) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD
   !REAL(8), DIMENSION(9,size(Aij,1)) :: dwr_press
   !REAL(8), DIMENSION(9,size(Aij,1),size(Aij,2)) :: dHij
   REAL(8), allocatable, DIMENSION(:,:) :: dwr_press
   REAL(8), allocatable, DIMENSION(:,:,:) :: dHij
   integer i, j, k, n
   if (numpar%p_const) then	! calculate this for P=const Parrinello-Rahman MD
      n = size(Aij,1)
      allocate(dwr_press(9,n))
      allocate(dHij(9,n,size(Aij,2)))
      dHij = 0.0d0
      dwr_press = 0.0d0
      call dHamil_tot_Press(Scell(NSC)%MDatoms, Scell, NSC, numpar, TB_Hamil, dHij)
      !$omp PARALLEL DO private(i,j)
      do i = 1, n
          do j = 1,9
            dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(i,:)) ! old, tested, good
            !dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(:,i)) ! wrong!!
         enddo ! j
      enddo ! i
      !$OMP END PARALLEL DO
      Scell(NSC)%SCforce%att = 0.0d0
      do i = 1,3
         do k = 1,3
            !print*, Scell(NSC)%SCforce%att(k, i), SUM(dwr_press((i-1)*3+k,:))
            !Scell(NSC)%SCforce%att(k, i) = Scell(NSC)%SCforce%att(k, i) + SUM(dwr_press((i-1)*3+k,:))
            Scell(NSC)%SCforce%att(k, i) = SUM(dwr_press((i-1)*3+k,:))
         enddo ! k
      enddo ! i
!       pause 'Attract_TB_Forces_Press'
      deallocate(dwr_press)
      deallocate(dHij)
   endif
end subroutine Attract_TB_Forces_Press


subroutine dHamil_tot_Press(atoms, Scell, NSC, numpar, TB_Hamil, dHij)
! construct the whole Hamilton matrix:
! (with respect to which Rk we take the derivatives, Appendix F of H.Jeschke PhD Thesis)
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   REAL(8), DIMENSION(:,:,:), INTENT(inout) :: dHij
   REAL(8), DIMENSION(9,4,4) :: dHij1
   integer :: i, j, j1, i1, ki, atom_2, m, nat, k, i2, j2, NumTB
   integer i4, j4
   !nat = matter%Na
   nat = size(atoms)
   dHij = 0.0d0
   do i = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(i)
      i4 = (i-1)*4
      do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            j4 = (j-1)*4
            !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA)
            call dHamilton_one_Press(i, atom_2, Scell(NSC)%MDatoms, Scell, NSC, TB_Hamil(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA), dHij1)
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis.
            do j1 = 1,4 ! all orbitals
               j2 = j4+j1
               do i1 = 1,4 ! all orbitals
                  !do k = 1,9
                     !dHij(k, (i-1)*4+i1, (j-1)*4+j1) = dHij1(k,i1,j1)	! construct the total Hamiltonian from 
                  !enddo ! k
                  i2 = i4+i1
                  dHij(:,i2,j2) = dHij1(:,i1,j1)	! construct the total Hamiltonian from
!                   dHij(1, i2, j2) = dHij1(1,i1,j1)	! construct the total Hamiltonian from
!                   dHij(2, i2, j2) = dHij1(4,i1,j1)	! construct the total Hamiltonian from
!                   dHij(3, i2, j2) = dHij1(7,i1,j1)	! construct the total Hamiltonian from
!                   dHij(4, i2, j2) = dHij1(2,i1,j1)	! construct the total Hamiltonian from
!                   dHij(5, i2, j2) = dHij1(5,i1,j1)	! construct the total Hamiltonian from
!                   dHij(6, i2, j2) = dHij1(8,i1,j1)	! construct the total Hamiltonian from
!                   dHij(7, i2, j2) = dHij1(3,i1,j1)	! construct the total Hamiltonian from
!                   dHij(8, i2, j2) = dHij1(6,i1,j1)	! construct the total Hamiltonian from
!                   dHij(9, i2, j2) = dHij1(9,i1,j1)	! construct the total Hamiltonian from
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0) then
      enddo ! j
   enddo ! i
end subroutine dHamil_tot_Press ! CHECKED


subroutine dHamilton_one_Press(i, atom_2, atoms, Scell, NSC, TB, dHij_press)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   integer(4), INTENT(IN) :: i, atom_2
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Pettifor), intent(in) ::  TB	! all tight binding parameters   
   REAL(8), DIMENSION(:,:,:), INTENT(out) :: dHij_press
   integer :: ki, kj, j
   REAL(8), DIMENSION(9,4,4) :: dts_press ! hopping integrals
   dHij_press = 0.0d0
   j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      dHij_press = 0.0d0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_Press_h(i, atom_2, atoms, Scell, NSC, TB, dts_press)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1,4
         do kj = 1,4
            dHij_press(:,kj,ki) = dts_press(:,ki,kj)   ! Hopping Integrals ! OLD, TESTED, GOOD
            !dHij_press(:,ki,kj) = dts_press(:,ki,kj)   ! Hopping Integrals  ! WRONG!! IMMEDIATE DIVERGENCE
         enddo ! kj
      enddo  ! ki
   endif
end subroutine dHamilton_one_Press


subroutine dHopping_Press_h(i, atom_2, atoms, Scell, NSC, TB, dts_press)
! subroutine making the derivatives of the hopping integrals
   INTEGER, intent(in) :: i,atom_2 ! atoms
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_Pettifor), intent(in) :: TB	! all tight binding parameters
   REAL(8), DIMENSION(:,:,:), INTENT(out) :: dts_press
   REAL(8), DIMENSION(9,4,4) :: dtsx
   integer :: k, ik
   real(8) :: x,y,z, x0, y0, z0, sx, sy, sz	! relative distances (projections) between the atoms 
   real(8) r, r1, diff_dircos_l, diff_dircos_m, diff_dircos_n, drijdhgd
   INTEGER(4) :: j,i1,j1,k1, dik, djk, dij
   real(8), dimension(3) :: zb
   real(8) sxr, syr, szr, Vs2r, dVs2r, Vs3r, dVs3r, Vs4r, dVs4r, x_r, y_r, z_r
   real(8), dimension(4) :: Vsr, dVsr

   dts_press = 0.0d0

   !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
   j = Scell(NSC)%Near_neighbor_list(i,atom_2)   ! it interacts with this atom
   x = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
   y = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
   z = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
   r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
   sx = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
   sy = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
   sz = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

   sxr = sx/r
   syr = sy/r
   szr = sz/r

   x_r = x/r
   y_r = y/r
   z_r = z/r

   Vs2r = Vs(TB,2,r)
   dVs2r = dVs(TB,2,r)
   Vs3r = Vs(TB,3,r)
   dVs3r = dVs(TB,3,r)
   Vs4r = Vs(TB,4,r)
   dVs4r = dVs(TB,4,r)

   Vsr(1) = 0.0d0
   Vsr(2) = Vs2r
   Vsr(3) = Vs3r
   Vsr(4) = Vs4r
   dVsr(1) = dVs(TB,1,r)
   dVsr(2) = dVs2r
   dVsr(3) = dVs3r
   dVsr(4) = dVs4r

   ! All of the hopping integrals ts for each pair of atoms i and j,
   ! according to the Eqs.(2.23)-(2.29) from PhD Thesis of H.Jeschke, Page 38
  do k1 = 1,9 ! all the components of the h_alpha_beta(3,3):
   diff_dircos_l = 0.0d0
   diff_dircos_m = 0.0d0
   diff_dircos_n = 0.0d0
   if (k1 .EQ. 1) then ! (1,1) ax
      drijdhgd = x*sxr
      diff_dircos_l = sxr
   else if (k1 .EQ. 2) then ! (2,1) ay
      drijdhgd = y*sxr
      diff_dircos_m = sxr
   else if (k1 .EQ. 3) then ! (3,1) az
      drijdhgd = z*sxr
      diff_dircos_n = sxr
   else if (k1 .EQ. 4) then ! (1,2) bx 
      drijdhgd = x*syr
      diff_dircos_l = syr
   else if (k1 .EQ. 5) then ! (2,2) by
      drijdhgd = y*syr
      diff_dircos_m = syr
   else if (k1 .EQ. 6) then ! (3,2) bz
      drijdhgd = z*syr
      diff_dircos_n = syr
   else if (k1 .EQ. 7) then ! (1,3) cx
      drijdhgd = x*szr
      diff_dircos_l = szr
   else if (k1 .EQ. 8) then ! (2,3) cy
      drijdhgd = y*szr
      diff_dircos_m = szr
   else if (k1 .EQ. 9) then ! (3,3) cz
      drijdhgd = z*szr
      diff_dircos_n = szr
   endif

   diff_dircos_l = diff_dircos_l - x*drijdhgd/(r*r)
   diff_dircos_m = diff_dircos_m - y*drijdhgd/(r*r)
   diff_dircos_n = diff_dircos_n - z*drijdhgd/(r*r)
   
!    print*, k1, drijdhgd, drij_dhab(x, sx, r)
!    print*, k1, drijdhgd, drij_dhab(y, sx, r)
!    print*, k1, drijdhgd, drij_dhab(z, sx, r)
!    print*, k1, drijdhgd, drij_dhab(x, sy, r)
!    print*, k1, drijdhgd, drij_dhab(y, sy, r)
!    print*, k1, drijdhgd, drij_dhab(z, sy, r)
!    print*, k1, drijdhgd, drij_dhab(x, sz, r)
!    print*, k1, drijdhgd, drij_dhab(y, sz, r)
!    print*, k1, drijdhgd, drij_dhab(z, sz, r)
!    pause
   
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 1, x, x, sx, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 2, x, y, sx, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 3, x, z, sx, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 1, x, x, sy, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 2, x, y, sy, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 3, x, z, sy, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 1, x, x, sz, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 2, x, y, sz, r)
!    print*, 'diff_dircos_l', k1, diff_dircos_l, dda_dhgd(1, 3, x, z, sz, r)
!    pause

!    call d_KS_sp3(Vsr, dVsr*drijdhgd, x_r, y_r, z_r, diff_dircos_l, diff_dircos_m, diff_dircos_n, dtsx(k1,:,:))
!    print*, k1, dtsx(k1,:,:)
!    pause 'TEST 1'

   dtsx(k1,1,1) = dVs(TB,1,r)*drijdhgd ! ss
   dtsx(k1,2,1) = -(diff_dircos_l*Vs2r + x_r*dVs2r*drijdhgd)  ! -spx
   dtsx(k1,3,1) = -(diff_dircos_m*Vs2r + y_r*dVs2r*drijdhgd)  ! -spy
   dtsx(k1,4,1) = -(diff_dircos_n*Vs2r + z_r*dVs2r*drijdhgd)  ! -spz
   dtsx(k1,1,2) = -dtsx(k1,2,1)
   dtsx(k1,2,2) = 2.0d0*x_r*(diff_dircos_l)*(Vs3r - Vs4r) + (x_r*x_r*(dVs3r - dVs4r) + dVs4r)*drijdhgd ! pxpx
   dtsx(k1,3,2)=(x_r*(diff_dircos_m)+y_r*(diff_dircos_l))*(Vs3r-Vs4r)+x_r*y_r*(dVs3r-dVs4r)*drijdhgd ! pxpy
   dtsx(k1,4,2)=(x_r*(diff_dircos_n)+z_r*(diff_dircos_l))*(Vs3r-Vs4r)+x_r*z_r*(dVs3r-dVs4r)*drijdhgd ! pxpz
   dtsx(k1,1,3) = -dtsx(k1,3,1)
   dtsx(k1,2,3)=dtsx(k1,3,2)
   dtsx(k1,3,3) = 2.0d0*y_r*(diff_dircos_m)*(Vs3r-Vs4r)+(y_r*y_r*(dVs3r-dVs4r)+dVs4r)*drijdhgd ! pypy
   dtsx(k1,4,3)=(y_r*(diff_dircos_n)+z_r*(diff_dircos_m))*(Vs3r-Vs4r)+y_r*z_r*(dVs3r-dVs4r)*drijdhgd ! pypz
   dtsx(k1,1,4) = -dtsx(k1,4,1)
   dtsx(k1,2,4)=dtsx(k1,4,2)
   dtsx(k1,3,4)=dtsx(k1,4,3)
   dtsx(k1,4,4)=2.0d0*z_r*(diff_dircos_n)*(Vs3r-Vs4r)+(z_r*z_r*(dVs3r-dVs4r)+dVs4r)*drijdhgd ! pzpz
   
!    print*, k1, dtsx(k1,:,:)
!    pause 'TEST 2'
   
  enddo ! k1 = 1,9 ! all the components of the h_alpha_beta(3,3):
  dts_press = dtsx
end subroutine dHopping_Press_h


! Tight Binding Hamiltonian within Pettifor parametrization:
subroutine construct_TB_H_Pettifor(numpar, matter, TB_Hamil, Scell, NSC, Ha, Err)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! materil parameters
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), intent(out) :: Ha ! Constructed Hamiltonian
   type(Error_handling), intent(inout) :: Err	! error save
   character(200) :: Error_descript
   integer i
   Error_descript = ''

   ! Get lists of nearest neighbors:
!    call get_near_neighbours(Scell) ! see "Atomic_tools"

   Scell(NSC)%Ha0 = Scell(NSC)%Ha ! save Hamiltonial from previous time-step
   Scell(NSC)%Ei0 = Scell(NSC)%Ei ! save energy levels for the next timestep
   Scell(NSC)%H_non0 = Scell(NSC)%H_non	! save non-diagonalized Hamiltonian from last time-step

!    call print_time('BEFORE', ind=1)

   ! Construct TB Hamiltonian (with Pettifor parameters):
   call Hamil_tot(numpar, Scell, NSC, TB_Hamil, Ha) ! see below
   Scell(NSC)%H_non = Ha ! save non-diagonalized Hamiltonian

!    call print_time('AFTER', ind=1)

   ! Diagonalize the Hamiltonian to get electron energy levels:
   call sym_diagonalize(Ha, Scell(NSC)%Ei, Error_descript, check_M=.true.)
   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine construct_TB_H_Pettifor: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_Pettifor



subroutine Complex_Hamil_tot(numpar, Scell, NSC, atoms, TB, Hij, CHij, ksx, ksy, ksz)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB   ! parameters of the Hamiltonian of TB
   REAL(8), DIMENSION(:,:), allocatable, INTENT(inout), optional :: Hij	 ! Real Hamiltonian at Gamma-point
   real(8), intent(in) :: ksx, ksy, ksz ! k-point to get Hamiltonian [relative]
   complex, DIMENSION(:,:), allocatable, INTENT(inout), optional :: CHij ! Complex Hamiltonian at k-point
   integer nat, Ne, i, j, j1, i1, k, l, m, atom_2, NumTB, Nsiz
   integer :: cell_x, cell_y, cell_z ! cell numbers
   real(8) :: x, y, z, r, sx, sy, sz, temp
   real(8) :: kx, ky, kz ! k-point to get Hamiltonian at, [1/A]
   REAL(8), DIMENSION(4,4) :: Hij1
   complex(8) :: C_I, expfac

   nat = size(atoms) ! number of atoms
   Ne = Scell(NSC)%Ne  ! initial number of VB electrons
   if (present(Hij)) then
      Nsiz = size(Hij,1)
   else
      Nsiz = size(TB(1,1)%V0)*Scell(1)%Na
   endif

   if (numpar%optic_model .EQ. 2) then
!       temp = g_me*g_e/g_h*1d-10
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0
      if (.not.allocated(Scell(NSC)%cPRRx)) allocate(Scell(NSC)%cPRRx(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%cPRRy)) allocate(Scell(NSC)%cPRRy(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%cPRRz)) allocate(Scell(NSC)%cPRRz(Nsiz,Nsiz))
      Scell(NSC)%cPRRx = (0.0d0,0.0d0)
      Scell(NSC)%cPRRy = (0.0d0,0.0d0)
      Scell(NSC)%cPRRz = (0.0d0,0.0d0)
   endif

   call Reciproc(Scell(NSC)%supce, Scell(NSC)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
   call Reciproc_rel_to_abs(ksx, ksy, ksz, Scell, NSC, kx, ky, kz) ! get absolute k-values
!     print*, matter%k_supce
!     print*, 'k', kx, ky, kz
!     pause 'PAUSE'

   if (.not.present(Hij) .AND. .not.present(CHij)) then
      print*, 'Neither real, nor complex Hamiltonian is present,'
      print*, 'so, what do you expect this subroutine to do?'
      print*, 'The subroutine "Complex_Hamil_tot" from the module TB_Hamiltonian'
      print*, 'requires at least one of them to be present.'
   endif

   if (present(Hij) .AND. .not.present(CHij)) then ! just get real Hamiltonian:
      call Hamil_tot(numpar, Scell, NSC, TB, Hij) ! see below
   endif

   if (.not.present(Hij) .AND. present(CHij)) then ! get complex Hamiltonian and the real one as well:
      C_I = dcmplx(0.0d0,1.0d0)	! complex unity
      if (.not.allocated(CHij)) then
         allocate(CHij(Ne,Ne))
         !CHij = cmplx(0.0d0,0.0d0) ! old!!!
      endif
      CHij = dcmplx(0.0d0,0.0d0)
      !if (.not.allocated(CHij1)) allocate(CHij1(Ne,Ne))
      do j = 1,nat	! all atoms
         m = Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one  
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            endif ! (atom_2 .EQ. 0)
            if (i .GT. 0) then
               !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)
               call Hamilton_one(Scell, NSC, i, j, atoms, TB(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA), Hij1, x1=x, y1=y, z1=z, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z)
               !expfac = exp(2.0d0*g_Pi*C_I*(kx*real(cell_x) + ky*real(cell_y) + kz*real(cell_z)))
               if ((abs(kx) .LT. 1.0d-14) .AND. (abs(ky) .LT. 1.0d-14) .AND. (abs(kz) .LT. 1.0d-14)) then
                  expfac = cmplx(1.0d0,0.0d0)
               else
                  expfac = exp(C_I*(kx*x + ky*y + kz*z))
!                   expfac = cexp(C_I*(kx*x + ky*y + kz*z))
               endif
               
               do j1 = 1,4 ! all orbitals
                  l = (j-1)*4+j1
                  do i1 = 1,4 ! all orbitals
                     k = (i-1)*4+i1
                     !if ((j1 .EQ. 1) .AND. (i1 .EQ. 1)) then
                     CHij(k,l) = DCMPLX(Hij1(j1,i1),0.0d0)*expfac ! make complex number of two real ones
!                      print*, 'a', Hij1(j1,i1), expfac

                     !CHij(k,l) = CMPLX(Hij1(i1,j1),0.0d0)*expfac ! make complex number of two real ones

!                      if (i .EQ. j) then
!                         CHij(k,l) = CMPLX(Hij1(j1,i1),0.0d0) ! make complex number of two real ones
!                         !CHij(k,l) = CMPLX(Hij1(i1,j1),0.0d0) ! make complex number of two real ones
!                      else
!                         CHij(k,l) = CMPLX(Hij1(j1,i1),0.0d0)*expfac ! make complex number of two real ones
!                         !CHij(k,l) = CMPLX(Hij1(i1,j1),0.0d0)*expfac ! make complex number of two real ones
!                      endif

                     if (numpar%optic_model .EQ. 2) then ! create matrix element:
                        if (i .EQ. j) then ! for the dielectric function, according to Trani:
                           Scell(NSC)%cPRRx(k,l) = DCMPLX(temp*0.270d0,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = DCMPLX(temp*0.270d0,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = DCMPLX(temp*0.270d0,0.0d0)*CHij(k,l)
                        else
                           Scell(NSC)%cPRRx(k,l) = DCMPLX(temp*x,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = DCMPLX(temp*y,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = DCMPLX(temp*z,0.0d0)*CHij(k,l)
                        endif
                        if (real(Scell(NSC)%cPRRx(k,l)) .GT. 1d10) print*, 'cPRRx', i, j, Scell(NSC)%cPRRx(k,l)
                        if (real(Scell(NSC)%cPRRy(k,l)) .GT. 1d10) print*, 'cPRRy',i, j, Scell(NSC)%cPRRy(k,l)
                        if (real(Scell(NSC)%cPRRz(k,l)) .GT. 1d10) print*, 'cPRRz',i, j, Scell(NSC)%cPRRz(k,l)
                     endif

                  enddo ! i1
               enddo !j1
            endif ! (i .GT. 0)
         enddo ! atom_2
      enddo ! j
   endif

   if (present(Hij) .AND. present(CHij)) then ! get complex Hamiltonian from the real one:
      C_I = dcmplx(0.0d0,1.0d0)	! complex unity
      if (.not.allocated(CHij)) then
         allocate(CHij(Ne,Ne))
      endif
      CHij = dcmplx(0.0d0,0.0d0)
      do j = 1,nat	! all atoms
         m = Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one  
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            endif ! (atom_2 .EQ. 0)
            if (i .GT. 0) then
               call shortest_distance(Scell, NSC, atoms, i, j, r, x1=x, y1=y, z1=z, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z)
               !expfac = exp(2.0d0*g_Pi*C_I*(kx*real(cell_x) + ky*real(cell_y) + kz*real(cell_z)))
               if ((abs(kx) .LT. 1.0d-14) .AND. (abs(ky) .LT. 1.0d-14) .AND. (abs(kz) .LT. 1.0d-14)) then
                  expfac = dcmplx(1.0d0,0.0d0)
               else
!                   expfac = cexp(C_I*(kx*x + ky*y + kz*z))
                  expfac = exp(C_I*(kx*x + ky*y + kz*z))
               endif
               do j1 = 1,4 ! all orbitals
                  l = (j-1)*4+j1
                  do i1 = 1,4 ! all orbitals
                     k = (i-1)*4+i1
                     !if ((j1 .EQ. 1) .AND. (i1 .EQ. 1)) then
                     !if ((l .EQ. k) .or. ((j1 .EQ. 1) .AND. (i1 .EQ. 1))) then
                     CHij(k,l) = DCMPLX(Hij(k,l),0.0d0)*expfac
!                       if (i .EQ. j) then
!                          CHij(k,l) = CMPLX(Hij(k,l),0.0d0)
!                          !CHij(k,l) = Hij(l,k)
!                       else
!                          CHij(k,l) = CMPLX(Hij(k,l),0.0d0)*expfac
!                           !CHij(k,l) = Hij(l,k)*expfac
!                       endif

                     if ((isnan(real(CHij(k,l)))) .OR. isnan(aimag(CHij(k,l)))) then
                        print*, i, j, k, l, CHij(k,l)
                        pause 'ISNAN'
                     endif

                     if (numpar%optic_model .EQ. 2) then ! create matrix element:
                        if (i .EQ. j) then ! for the dielectric function, according to Trani:
                           Scell(NSC)%cPRRx(k,l) = temp*0.270d0*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = temp*0.270d0*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = temp*0.270d0*CHij(k,l)
                        else
                           Scell(NSC)%cPRRx(k,l) = temp*x*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = temp*y*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = temp*z*CHij(k,l)
                        endif
                        if (real(Scell(NSC)%cPRRx(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRx(k,l)
                        if (real(Scell(NSC)%cPRRy(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRy(k,l)
                        if (real(Scell(NSC)%cPRRz(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRz(k,l)
                     endif
                  enddo ! i1
               enddo !j1
            endif ! (i .GT. 0)
         enddo ! atom_2
      enddo ! j
!       CHij = CONJG(CHij) ! complex conjugate
   endif

end subroutine Complex_Hamil_tot



! hamiltonian for atoms:
subroutine Hamil_tot(numpar, Scell, NSC, TB_Hamil, Hij)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_Pettifor), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   REAL(8), DIMENSION(:,:), INTENT(inout) :: Hij	 ! Hamiltonian
   REAL(8), DIMENSION(4,4) :: Hij1
   REAL(8), DIMENSION(4) :: Vij1
   integer :: nat, Ne, NumTB, Nsiz
   integer :: i, j, j1, i1, k, l, atom_2, m, FN
   real(8) :: r, x, y, z, temp
   nat = Scell(NSC)%Na
   Ne = Scell(NSC)%Ne
   Nsiz = size(Hij,1)
   if (numpar%optic_model .EQ. 3) then
!       temp = g_me*g_e/g_h*1d-10
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0
      if (.not.allocated(Scell(NSC)%PRRx)) allocate(Scell(NSC)%PRRx(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%PRRy)) allocate(Scell(NSC)%PRRy(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%PRRz)) allocate(Scell(NSC)%PRRz(Nsiz,Nsiz))
      Scell(NSC)%PRRx = 0.0d0
      Scell(NSC)%PRRy = 0.0d0
      Scell(NSC)%PRRz = 0.0d0
   endif

   Hij = 0.0d0

   ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis:
   do j = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one  
         if (atom_2 .EQ. 0) then
            i = j
         else
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         endif
         if (i .GT. 0) then ! if there really is a nearest neighbor

            ! Which kinds of atoms are these (which TB parameters to use):
            !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)
            if (numpar%optic_model .EQ. 3) then ! create matrix element:
               !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z) ! module "Atomic_tools"
               call Hamilton_one(Scell, NSC, i, j, Scell(NSC)%MDatoms, TB_Hamil(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA), Hij1, x1=x, y1=y, z1=z) ! block-hamiltonian and the shortest distances
            else
               call Hamilton_one(Scell, NSC, i, j, Scell(NSC)%MDatoms, TB_Hamil(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA), Hij1) ! this calles the block-hamiltonian
            endif
            do j1 = 1,4 ! all orbitals
               l = (j-1)*4+j1
               do i1 = 1,4 ! all orbitals
                  k = (i-1)*4+i1
                  !Hij((i-1)*4+i1, (j-1)*4+j1) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  Hij(k,l) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)

                  if (numpar%optic_model .EQ. 3) then ! create matrix element:
                     if (i .EQ. j) then ! for the dielectric function, according to Trani:
                        Scell(NSC)%PRRx(k,l) = temp*0.270d0*Hij(k,l)
                        Scell(NSC)%PRRy(k,l) = temp*0.270d0*Hij(k,l)
                        Scell(NSC)%PRRz(k,l) = temp*0.270d0*Hij(k,l)
                     else
                        Scell(NSC)%PRRx(k,l) = temp*x*Hij(k,l)
                        Scell(NSC)%PRRy(k,l) = temp*y*Hij(k,l)
                        Scell(NSC)%PRRz(k,l) = temp*z*Hij(k,l)
                     endif
                     if (Scell(NSC)%PRRx(k,l) .GT. 1d10) print*, i, j, Scell(NSC)%PRRx(k,l)
                     if (Scell(NSC)%PRRy(k,l) .GT. 1d10) print*, i, j, Scell(NSC)%PRRy(k,l)
                     if (Scell(NSC)%PRRz(k,l) .GT. 1d10) print*, i, j, Scell(NSC)%PRRz(k,l)
!                      write(*,'(i4,i4,e,e,e)') k,l, matter%PRRx(k,l), matter%PRRy(k,l), matter%PRRz(k,l)
                  endif
               enddo ! i1
            enddo ! j1
         endif ! (i > 0)
      enddo ! j
   enddo ! i

!   open(NEWUNIT=FN, FILE = 'OUTPUT_Hamiltonain_Si.dat')
!      do i = 1, size(Hij,1)
!         do j = 1, size(Hij,2)
!            write(FN, '(i,i,e)') i, j, Hij(i,j)
!         enddo
!      enddo
!   close(FN)
!    pause 'Hamil_tot'

end subroutine Hamil_tot


subroutine Hamilton_one(Scell, NSC, i, j, atoms, TB, Hij,  x1, y1, z1, sx1, sy1, sz1, cell_x, cell_y, cell_z)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer(4), INTENT(IN) :: i, j
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Pettifor), intent(in) :: TB	! all tight binding parameters
   REAL(8), DIMENSION(4,4), INTENT(out) :: Hij  ! hamiltonian
   real(8), intent(out), optional :: x1, y1, z1    ! shortest distances
   real(8), intent(out), optional :: sx1, sy1, sz1 ! shortest distance in relative coordinates
   integer, intent(out), optional :: cell_x, cell_y, cell_z ! cell numbers
   integer(4) ki, kj, i1, j1, k1, ik
   REAL(8), DIMENSION(4,4) :: ts ! hopping integrals
   real(8) x,y,z,r,r1, x0, y0, z0, sx, sy, sz !  interatomic distance projections, and the total distance
   integer :: cell_x1, cell_y1, cell_z1 ! cell numbers
   real(8), dimension(3) :: zb

   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      Hij = 0.0d0   ! Nondiagonals are zeros
      Hij(1,1) = TB%Es ! diagonals are equal to onsite energies
      do ki = 2,4
         Hij(ki,ki) = TB%Ep ! diagonals are equal to onsite energies 
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
      !call Hopping_old(TB, ts, -x, -y, -z)
      call Hopping(TB, ts, x, y, z, r)
      !call Hopping_NEW(TB, ts, x, y, z, r) ! this subroutine is ~25% slower than the old one
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1, 4
         do kj = 1,4
            Hij(kj,ki) = ts(kj,ki)   ! Hopping Integrals
         enddo ! kj
      enddo  ! ki
   endif
end subroutine Hamilton_one



subroutine Hopping(TB, ts, x, y, z, r0)
! subroutine making the hopping integrals
! ts --- the Hopping Integrals themselves
! x,y,z --- the distances between naighbour atoms we are analizing
   type(TB_H_Pettifor), intent(in), target :: TB	! all tight binding parameters
   REAL(8), DIMENSION(4,4), INTENT(out) :: ts
   real(8), intent(in) :: x,y,z	! relative distances (projections) between the atoms
   real(8), intent(in), optional :: r0
   !=============================================
   real(8), pointer :: rm
   INTEGER(4) :: i,j, k
   INTEGER(4) :: i1,j1
   real(8) r, r1, xr, yr, zr, Vsr(4)
   if (present(r0)) then
      r = r0 ! it is given, no need to recalculate
   else
      r = SQRT(x*x + y*y + z*z)   ! total distance between the two atoms within the supercell
   endif
   
   rm => TB%rm ! cut-off radius
   
   if (r <= rm) then   ! these atoms do interact:
      ! All of the hopping integrals ts for each pair of atoms i and j,
      ! according to the Eqs.(2.23)-(2.29) from PhD Thesis of H.Jeschke, Page 38
      xr = x/r
      yr = y/r
      zr = z/r
      Vsr(1) = Vs(TB,1,r)	! calling the function given below
      Vsr(2) = Vs(TB,2,r)	! calling the function given below
      Vsr(3) = Vs(TB,3,r)	! calling the function given below
      Vsr(4) = Vs(TB,4,r)	! calling the function given below
      ts(1,1) = Vsr(1)
      ts(1,2) = -xr*Vsr(2)
      ts(1,3) = -yr*Vsr(2)
      ts(1,4) = -zr*Vsr(2)
      ts(2,1) = xr*Vsr(2)
      ts(2,2) = xr*xr*(Vsr(3) - Vsr(4)) + Vsr(4)
      ts(2,3) = xr*yr*(Vsr(3) - Vsr(4))
      ts(2,4) = xr*zr*(Vsr(3) - Vsr(4))
      ts(3,1) = yr*Vsr(2)
      ts(3,2) = xr*yr*(Vsr(3) - Vsr(4))
      ts(3,3) = yr*yr*(Vsr(3) - Vsr(4)) + Vsr(4)
      ts(3,4) = yr*zr*(Vsr(3) - Vsr(4))
      ts(4,1) = zr*Vsr(2)
      ts(4,2) = xr*zr*(Vsr(3) - Vsr(4))
      ts(4,3) = yr*zr*(Vsr(3) - Vsr(4))
      ts(4,4) = zr*zr*(Vsr(3) - Vsr(4)) + Vsr(4)
      
!         print*, 'ts-1', ts
!         call KS_sp3_hetero(Vsr, Vsr, xr, yr, zr, ts)
!         print*, 'ts-2', ts
!         pause

   else ! these atoms too far to interact via covalent bonding:
      ts = 0.0d0
   endif

   nullify(rm)
end subroutine Hopping




function Vs(TB,i,r) ! functions Vsssigma, /p sigma, /p pi
! for the Hopping Integrals according to Eqs.(2.23)-(2.29), (2.32), 
! and the Appendix E with all parameters,
! from H.Jeschke PhD thesis.
! i   ---  stands for index of V, i=1 is "s s sigma", i=2 is "s p sigma", i=3 is "p p sigma", i=4 is "p p pi"
! x, y, z --- are the coordinates of interatomic distance (positive or negative!) [A]
! r  --- absolute interatomic distance [A]
   type(TB_H_Pettifor), intent(in), target :: TB	! all tight binding parameters
   INTEGER(4), intent(in) :: i ! index of V
   real(8), intent(in) :: r  ! absolute distance between the atoms	
   real(8) Vs ! function itself (Eqs.(2.23)-(2.27))
   !------------------------
   real(8), dimension(:), pointer :: V0 ! multiplicator in V, Eq.(2.31)
   ! Parameters from Appenix E of Haralds PhD thesis:
   real(8), pointer :: r0, n, r1, rm
   real(8), dimension(:), pointer :: nc, rc, c0, c1, c2, c3
   real(8) r_r1, r_r12
   integer j

   V0 => TB%V0
   r0 => TB%r0
   n => TB%n
   r1 => TB%r1
   rm => TB%rm
   nc => TB%nc
   rc => TB%rc
   c0 => TB%c0
   c1 => TB%c1
   c2 => TB%c2
   c3 => TB%c3
   if (r .GT. rm) then ! cutt-off distance
      Vs = 0.0d0
   elseif (r .GE. r1) then
      r_r1 = r-r1
      r_r12 = r_r1*r_r1 
      Vs = V0(i)*(c0(i) + c1(i)*r_r1 + c2(i)*r_r12 + c3(i)*r_r12*r_r1) ! smooth zeroise in between
      !Vs = (c0(i) + c1(i)*r_r1 + c2(i)*r_r12 + c3(i)*r_r12*r_r1) ! smooth zeroise in between
   else 
      Vs = V0(i)*((r0/r)**n)*dexp(n*(-(r/rc(i))**nc(i) + (r0/rc(i))**nc(i))) ! close atoms
   endif

   nullify(V0, r0, n, r1, rm, nc, rc, c0, c1, c2, c3)
end function Vs


function dVs(TB,i,r,Vs_precalc) ! functions Vsssigma, /p sigma, /p pi
! for the Hopping Integrals according to Eqs.(2.23)-(2.29), (2.32), 
! and the Appendix E with all parameters,
! from H.Jeschke PhD thesis.
! i   ---  stands for index of V, i=1 is "s s sigma", i=2 is "s p sigma", i=3 is "p p sigma", i=4 is "p p pi"
! x, y, z --- are the coordinates of interatomic distance (positive or negative!)
! r  --- absolute interatomic distance
   type(TB_H_Pettifor), intent(in), target :: TB	! all tight binding parameters
   INTEGER(4), intent(in) :: i ! index of V
   real(8), intent(in) :: r  ! absolute distance between the atoms	
   real(8), intent(in), optional :: Vs_precalc	! precalculated Vs function
   real(8) dVs ! function itself (Eqs.(2.23)-(2.27))
   !------------
   real(8), dimension(:), pointer :: V0 ! multiplicator in V, Eq.(2.31)
   real(8), pointer :: r0, n, r1, rm
   real(8), dimension(:), pointer :: nc, rc, c0, c1, c2, c3
   real(8) r_r1
   V0 => TB%V0
   r0 => TB%r0
   n => TB%n
   r1 => TB%r1
   rm => TB%rm
   nc => TB%nc
   rc => TB%rc
   c0 => TB%c0
   c1 => TB%c1
   c2 => TB%c2
   c3 => TB%c3
   ! Different functional dependencies for different distances:
   if (r .GT. rm) then ! cut-off distance
      dVs = 0.0d0
   elseif (r .GT. r1) then
      r_r1 = r-r1
      dVs = V0(i)*(c1(i) + c2(i)*r_r1*2.0d0 + c3(i)*r_r1*r_r1*3.0d0) ! smooth zeroise in between
      !dVs = (c1(i) + c2(i)*r_r1*2.0d0 + c3(i)*r_r1*r_r1*3.0d0) ! smooth zeroise in between
   else
      !dVs = -n*Vs(TB,i,r)*(1.0d0/r + nc(i)/rc(i)*(r/rc(i))**(nc(i)-1.0d0))  ! close to atoms
      if (present(Vs_precalc)) then	! no need to recalculate it, it is provided:
         dVs = -n*Vs_precalc/r*(1.0d0 + nc(i)*(r/rc(i))**(nc(i)))  ! close to atoms
      else	! calculate Vs function:
         dVs = -n*Vs(TB,i,r)/r*(1.0d0 + nc(i)*(r/rc(i))**(nc(i)))  ! close to atoms
      endif
   endif
   nullify(V0, r0, n, r1, rm, nc, rc, c0, c1, c2, c3)
end function dVs



function d2Vs(TB, i, r, Vs_r, dVs_r)	! second derevative of the attractive potential contribution
   ! Hopping Integrals according to Eqs.(2.23)-(2.29), (2.32), Jeschke PhD thesis
   type(TB_H_Pettifor), intent(in), target :: TB	! all tight binding parameters
   INTEGER(4) i ! index of V
   real(8), intent(in) :: r, Vs_r, dVs_r
   real(8) d2Vs ! function itself (Eqs.(2.23)-(2.27))
   !------------
   real(8), dimension(:), pointer :: V0 ! multiplicator in V, Eq.(2.31)
   real(8), pointer :: r0, n, r1, rm
   real(8), dimension(:), pointer :: nc, rc, c0, c1, c2, c3
   real(8) r_r1, rrc, rrcnc
   V0 => TB%V0
   r0 => TB%r0
   n => TB%n
   r1 => TB%r1
   rm => TB%rm
   nc => TB%nc
   rc => TB%rc
   c0 => TB%c0
   c1 => TB%c1
   c2 => TB%c2
   c3 => TB%c3
   ! Different functional dependencies for different distances:
   if ((r .GT. rm) .or. (Vs_r <= 0.0d0)) then ! cut-off distance
      d2Vs = 0.0d0
   elseif (r .GT. r1) then
      r_r1 = r-r1
      d2Vs = V0(i)*(c2(i)*2.0d0 + c3(i)*r_r1*6.0d0) ! smooth zeroise in between
      !d2Vs = (c2(i)*2.0d0 + c3(i)*r_r1*6.0d0) ! smooth zeroise in between
   else
      d2Vs = (-1.0d0/r + dVs_r/Vs_r)*dVs_r - n/r*Vs_r*nc(i)*nc(i)/rc(i)*(r/rc(i))**(nc(i)-1)
   endif
   nullify(V0, r0, n, r1, rm, nc, rc, c0, c1, c2, c3)
end function d2Vs



!11111111111111111111111111111111111111111111111111111111111111111111111



!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
!Repulsive part of TB:

subroutine get_Erep_s(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Pettifor"
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_Pettifor), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   real(8) a_r, b, E_pot
   INTEGER(4) i1, j1, m, atom_2, NumTB
   real(8), DIMENSION(3) :: zb
   a = 0.0d0
   do i1 = 1, Scell(NSC)%Na   ! for all atoms
      b = 0.0d0
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if (j1 .NE. i1) then
            !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA)
            !call shortest_distance(matter, atoms, i1, j1, a_r) ! from module Atomic_tools
            a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            b = b + phi(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r)
         endif ! (j1 .NE. i1)
      enddo ! j1
      !a = a + fx(TB,b)
      ! Potential energy:
      E_pot = fx(TB_Repuls(Scell(NSC)%MDatoms(i1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),b) + TB_Repuls(Scell(NSC)%MDatoms(i1)%KOA,Scell(NSC)%MDatoms(i1)%KOA)%E0_TB
      ! Add to the total one:
      a = a + E_pot
      ! And save for each atom:
      Scell(NSC)%MDAtoms(i1)%Epot = E_pot
   enddo ! i1
end subroutine get_Erep_s


function fx(TB_Repuls,a_x) ! polinomial form, for repulsive potential
   ! uses Eq. 2.36, Page 40 in H.Jeschke PhD thesis.
   type(TB_Rep_Pettifor), intent(in), target :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_x
   real(8) fx
   real(8), pointer :: a0, a1, a2, a3, a4
   real(8) :: a_x2
   a0 => TB_Repuls%a0(1)
   a1 => TB_Repuls%a0(2)
   a2 => TB_Repuls%a0(3)
   a3 => TB_Repuls%a0(4)
   a4 => TB_Repuls%a0(5)
   a_x2 = a_x*a_x
   fx = a0 + a1*a_x + a2*a_x2 + a3*a_x2*a_x + a4*a_x2*a_x2
   nullify(a0, a1, a2, a3, a4) 
end function fx


function phi(TB_Repuls,a_x)	! repulsive potential itself
   ! it corresponds to the Eq. 2.37, Page 40 in H.Jeschke PhD thesis.
   type(TB_Rep_Pettifor), intent(in), target :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8) a_x
   real(8) phi
   !---------------
   real(8), pointer :: phi0, d0, m, mc, dc
   real(8), pointer :: c0, c1, c2, c3, d1, dm
   real(8) :: x_d1, x_d2
   phi0 => TB_Repuls%phi0
   d0 => TB_Repuls%d0
   m => TB_Repuls%m
   mc => TB_Repuls%mc
   dc => TB_Repuls%dc
   c0 => TB_Repuls%c0(1)
   c1 => TB_Repuls%c0(2)
   c2 => TB_Repuls%c0(3)
   c3 => TB_Repuls%c0(4)
   d1 => TB_Repuls%d1
   dm => TB_Repuls%dm
   if (a_x .GT. dm) then
      phi = 0.0d0
   elseif (a_x .GT. d1) then
      x_d1 = a_x-d1
      x_d2 = x_d1*x_d1
      phi = c0 + c1*x_d1 + c2*x_d2 + c3*x_d1*x_d2
   elseif (a_x .LE. 0.0d0) then
      phi = 0.0d0
   else
      phi = phi0*((d0/a_x)**m)*exp(m*(-(a_x/dc)**mc + (d0/dc)**mc))
   endif
   nullify(phi0, d0, m, mc, dc, c0, c1, c2, c3, d1, dm)
end function phi



subroutine dErdr_s(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   type(TB_Rep_Pettifor), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   REAL(8), DIMENSION(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, x, y, z, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, j1, ian, dik, djk, n, m, atom_2, NumTB
   real(8), DIMENSION(3) :: zb
   real(8), dimension(:,:), allocatable :: Erx_s
   n = size(atoms)
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

   !$omp PARALLEL DO private(ian,i1,dpsi,psi,m,j1,x,y,z,a_r,dik,djk,x1,b,ddlta,b_delta,a,atom_2,NumTB)
   do ian = 1, n  ! Forces for all atoms
     Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
     do i1 = 1, n
         dpsi = 0.0d0
         psi = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i1)
         do atom_2 = 1,m ! do only for atoms close to that one
            j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
            if (j1 .GT. 0) then
                !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA)
                if (ian .EQ. i1) then
                   dik = 1
                else
                   dik = 0
                endif
                if (ian .EQ. j1) then
                   djk = 1
                else
                   djk = 0
                endif

                if ((j1 .NE. i1) .OR. ((dik-djk) .NE. 0)) then ! without it, it gives ERROR
                  x = Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X
                  y = Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z = Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z
                  a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R

                  psi = psi + phi(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r) ! Eq.(F21), p.147 H.Jeschke PhD Thesis

                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
                  x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                  x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
!                   x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
!                   x1(2) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
!                   x1(3) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)
                  b = dphi(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r)
                  ddlta = real(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
               endif
            endif
         enddo ! j1

         a = dfx(TB_Repuls(Scell(NSC)%MDatoms(i1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),psi)
         Erx_s(1,ian) = Erx_s(1,ian) + a*dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + a*dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + a*dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Erx_s(:,ian) ! all repulsive forces
   enddo ! ian
   !$OMP END PARALLEL DO
   deallocate(Erx_s)
END subroutine dErdr_s



subroutine dE2rep_dr2(TB_Repuls, atoms, Scell, numpar, F, dF) ! second derivatives of the repulsive energy by r
   ! CHECKED, WORK CORRECTLY: dE/dr = 1/h*dE/ds
   type(TB_Rep_Pettifor), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:), allocatable, intent(out) :: F, dF	! force and its derivative
   !---------------------------------------------------
   real(8) dpsi(3), dpsi2(3),  ddlta, dik, djk
   real(8) drdrx, drdry, drdrz, d2rdr2x, d2rdr2y, d2rdr2z
   integer :: NSC ! number of supercell
   integer i1, j1, ian, n, m, atom_2
   real(8), dimension(:), allocatable, target :: psi_vec, dfx_vec, d2fx_vec
   real(8), dimension(:,:), allocatable :: Erx_s, Erx_s2
   real(8), dimension(:,:), allocatable, target :: M_phi, M_dPhi, M_d2Phi
   real(8), pointer :: phi_point, dphi_point, d2phi_point, x, y ,z, a_r
   NSC = 1	! we only have one supercell so far
   n = size(atoms)
   allocate(psi_vec(n))	! all psi's
   psi_vec = 0.0d0
   allocate(dfx_vec(n))	! all dfx's
   dfx_vec = 0.0d0
   allocate(d2fx_vec(n))	! all d2fx's
   d2fx_vec = 0.0d0
   
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0
   allocate(Erx_s2(3,n)) ! x,y,z-forces for each atoms
   Erx_s2 = 0.0d0
   
   if (.not.allocated(F)) allocate(F(3,n))
   if (.not.allocated(dF)) allocate(dF(3,n))
   F(:,:) = 0.0d0	! just to start with, forces
   dF(:,:) = 0.0d0	! just to start with, derivatives of forces
   
   ! Get matrices of phi(i,j) functions for all pairs of atoms (i,j):
   call construct_phi_matrix(Scell, TB_Repuls, M_Phi, M_dPhi, M_d2Phi)	! see below
   ! Get all psi(phi) functions for all atoms (i):
   call construct_psi_matrix(M_Phi, psi_vec)	! see below
   ! Get all fx(psi(phi)) functions for all atoms (i):
   call construct_dfx_matrix(Scell, TB_Repuls, psi_vec, dfx_vec, d2fx_vec)	! see below
   
   !$omp PARALLEL DO private(ian,i1,dpsi,dpsi2,m,j1,x,y,z,a_r,dik,djk,dphi_point,d2phi_point,ddlta, drdrx, drdry, drdrz, d2rdr2x, d2rdr2y, d2rdr2z, atom_2)
   do ian = 1, n	! Forces for all atoms
     do i1 = 1, n	! contribution from all other atoms
         dpsi = 0.0d0
         dpsi2 = 0.0d0
         dik = Kronecker_delta(ian, i1)	! module "Atomic_tools"
         m = Scell(NSC)%Near_neighbor_size(i1)
         do atom_2 = 1,m ! do only for atoms close to that one
            j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
            if (j1 .GT. 0) then
               djk = Kronecker_delta(ian, j1)	! module "Atomic_tools"
               ddlta = dik - djk
               if ((j1 .NE. i1) .OR. (ABS(ddlta) .GE. 1.0d-12)) then ! without it, it gives ERROR
                  x => Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X
                  y => Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z => Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z
                  a_r => Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R
                  
                  ! Derivatives d r_{i,j} / d r_{k,alpha}
                  drdrx = drij_drka(i1, j1, ian, x, a_r)	! module "TB_Koster_Slater"
                  drdry = drij_drka(i1, j1, ian, y, a_r)	! module "TB_Koster_Slater"
                  drdrz = drij_drka(i1, j1, ian, z, a_r)	! module "TB_Koster_Slater"
                  
                  ! Second derivatives d2 r_{ij} / d r2_{k,alpha}
                  d2rdr2x = ddija_drkb(i1, j1, ian, 1, 1, x, x, a_r)	! module "TB_Koster_Slater"
                  d2rdr2y = ddija_drkb(i1, j1, ian, 2, 2, y, y, a_r)	! module "TB_Koster_Slater"
                  d2rdr2z = ddija_drkb(i1, j1, ian, 3, 3, z, z, a_r)	! module "TB_Koster_Slater"

                  ! Precalculated patrices of Phi-functions and d Phi / d r_{i,j}
                  dphi_point => M_dPhi(j1, i1)
                  d2phi_point => M_d2Phi(j1, i1)
                  
                  ! d Psi / d r_{k,alpha}
                  dpsi(1)= dpsi(1) + dphi_point*drdrx
                  dpsi(2)= dpsi(2) + dphi_point*drdry
                  dpsi(3)= dpsi(3) + dphi_point*drdrz

                  ! d2 Psi / d r2_{k,alpha}
                  dpsi2(1) = dpsi2(1) + d2phi_point*drdrx*drdrx + dphi_point*d2rdr2x
                  dpsi2(2) = dpsi2(2) + d2phi_point*drdry*drdry + dphi_point*d2rdr2y
                  dpsi2(3) = dpsi2(3) + d2phi_point*drdrz*drdrz + dphi_point*d2rdr2z
               endif
            endif
         enddo ! j1
         !write(*,'(a,i3,i3,es,es,es)') 'test a2', ian, i1, d2fx_vec(i1), a2, ABS(d2fx_vec(i1)-a2)
         
         Erx_s(:,ian) = Erx_s(:,ian) + dfx_vec(i1)*dpsi(:)	! repulsive part in X, Y, Z-coordinate
         Erx_s2(:,ian) = Erx_s2(:,ian) + d2fx_vec(i1)*dpsi(:)*dpsi(:) + dfx_vec(i1)*dpsi2(:)	! derivaties of the repulsive part
      enddo ! i1
      F(:,ian) = Erx_s(:,ian)	! all repulsive forces for this atom
      dF(:,ian) = Erx_s2(:,ian)	! all derivatives of repulsive forces for this atom
   enddo ! ian
   !$OMP END PARALLEL DO
   
   deallocate(Erx_s, Erx_s2, M_Phi, M_dPhi, M_d2Phi,  psi_vec, dfx_vec, d2fx_vec)
   nullify(phi_point, dphi_point, d2phi_point, x, y ,z, a_r)
end subroutine dE2rep_dr2


subroutine construct_dfx_matrix(Scell, TB_Repuls, psi_vec, dfx_vec, d2fx_vec)
   type(TB_Rep_Pettifor), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   real(8), dimension(:), intent(in) :: psi_vec	! vector of Psy functions
   real(8), dimension(:), allocatable, intent(out) :: dfx_vec	! vector of a functions
   real(8), dimension(:), allocatable, intent(out) :: d2fx_vec	! vector of a2 functions
   integer :: i1
   if (.not.allocated(dfx_vec)) allocate(dfx_vec(size(psi_vec)))
   if (.not.allocated(d2fx_vec)) allocate(d2fx_vec(size(psi_vec)))
   !$omp PARALLEL DO private(i1)
   do i1 = 1, size(psi_vec)
      dfx_vec(i1) = dfx(TB_Repuls(Scell(1)%MDatoms(i1)%KOA,Scell(1)%MDatoms(i1)%KOA),psi_vec(i1))	! function below
      d2fx_vec(i1) = df2x(TB_Repuls(Scell(1)%MDatoms(i1)%KOA,Scell(1)%MDatoms(i1)%KOA),psi_vec(i1))	! function below
   enddo
   !$OMP END PARALLEL DO
end subroutine construct_dfx_matrix


subroutine construct_psi_matrix(M_Phi, psi_vec)
   real(8), dimension(:,:), intent(in) :: M_Phi	! matrix of Phi functions
   real(8), dimension(:), allocatable, intent(out) :: psi_vec	! vector of Psy functions
   integer :: i1
   if (.not.allocated(psi_vec)) allocate(psi_vec(size(M_Phi,1)))
   !$omp PARALLEL DO private(i1)
   do i1 = 1, size(psi_vec)
      psi_vec(i1) = SUM(M_Phi(:, i1))
   enddo
   !$OMP END PARALLEL DO
end subroutine construct_psi_matrix


subroutine construct_phi_matrix(Scell, TB_Repuls, M_Phi, M_dPhi, M_d2Phi)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   type(TB_Rep_Pettifor), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   real(8), dimension(:,:), allocatable, intent(out) :: M_Phi	! matrix of Phi functions
   real(8), dimension(:,:), allocatable, intent(out) :: M_dPhi	! matrix of first derivatives of Phi functions
   real(8), dimension(:,:), allocatable, intent(out), optional :: M_d2Phi	! matrix of second derivatives of of Phi functions
   !---------------------
   real(8) :: a_r
   integer :: N, i1, m, atom_2, j1
   N = size(Scell(1)%MDAtoms)	! number of atoms
   ! Start from allocating and setting to zero:
   if (.not.allocated(M_Phi)) allocate(M_Phi(N,N))
   M_Phi = 0.0d0
   if (.not.allocated(M_dPhi)) allocate(M_dPhi(N,N))
   M_dPhi = 0.0d0
   if (present(M_d2Phi)) then
      if (.not.allocated(M_d2Phi)) allocate(M_d2Phi(N,N))
      M_d2Phi = 0.0d0
   endif
   ! Calculate the nonzero elements:
   !$omp PARALLEL DO private(i1, m, atom_2, j1, a_r)
   do i1 = 1, N
      m = Scell(1)%Near_neighbor_size(i1)
      do atom_2 = 1,m ! do only for atoms close to that one
         j1 = Scell(1)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if (j1 .GT. 0) then
            a_r = Scell(1)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R
            ! Phi function:
            M_Phi(i1,j1) = phi(TB_Repuls(Scell(1)%MDatoms(j1)%KOA, Scell(1)%MDatoms(i1)%KOA),a_r) ! function above (or below)
            ! Derivative of Phi function:
            M_dPhi(i1,j1) = dphi(TB_Repuls(Scell(1)%MDatoms(j1)%KOA, Scell(1)%MDatoms(i1)%KOA), a_r, M_Phi(i1,j1) )	! function above (or below)
            ! Second derivative of Phi function:
            if (present(M_d2Phi)) then
               M_d2Phi(i1,j1) = d2phi(TB_Repuls(Scell(1)%MDatoms(j1)%KOA, Scell(1)%MDatoms(i1)%KOA), a_r, M_Phi(i1,j1), M_dPhi(i1,j1))	! function above (or below)
            endif ! present(M_d2Phi)
         endif ! (j1 .GT. 0)
      enddo ! atom_2
   enddo ! i1
   !$OMP END PARALLEL DO
end subroutine construct_phi_matrix



subroutine dErdr_Pressure_s(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by h
   type(TB_Rep_Pettifor), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms ! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   REAL(8), DIMENSION(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i,j,k,l, i1, i2, i3, ik, m, atom_2, n, NumTB
   real(8) x,y,z,sx,sy,sz,r, rcur(3), scur(3)
   real(8) df_psy, psi, dpsy

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = size(atoms)
      Scell(NSC)%SCforce%rep = 0.0d0
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         psi = 0.0d0
         dpsy = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i)
         do atom_2 = 1,m ! do only for atoms close to that one
            j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
            if ((i .NE. j) .AND. (j .GT. 0)) then
               !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)
               !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

               psi = psi + phi(TB_Repuls(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA),r) ! Eq.(F21), p.147 H.Jeschke PhD Thesis
               dpsy = dphi(TB_Repuls(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA),r)
               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! i=j
         enddo ! j
         df_psy = dfx(TB_Repuls(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(i)%KOA),psi)
         do l = 1,3 ! supce indices
            do k = 1,3  ! supce indices
               Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + df_psy*Rep_Pr(l,k)  ! checked
            enddo ! l
         enddo ! k
      enddo ! i
      
!       do k = 1,3 ! supce indices
!          do l = 1,3  ! supce indices
!             print*, l,k, Scell(NSC)%SCforce%rep(l,k)
!          enddo ! l
!       enddo ! k
   endif
end subroutine dErdr_Pressure_s


! The derivative of the function by its argument. To be multiplied later by the derivative of the argument by the variable:
! d Phi / d r_{ij} * d r_{ij} / d r_{k,beta}
function dphi(TB_Repuls, a_x, Phi_precalc)	! derevative of the repulsive potential itself
   ! it corresponds to the Eq. 2.37, Page 40 in H.Jeschke PhD thesis.
   type(TB_Rep_Pettifor), intent(in), target :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_x
   real(8), intent(in), optional :: Phi_precalc
   real(8) dphi
   !-----------------
   real(8) :: x_d1
   real(8), pointer :: phi0, d0, m, mc, dc
   real(8), pointer :: c0, c1, c2, c3, d1, dm
   phi0 => TB_Repuls%phi0
   d0 => TB_Repuls%d0
   m => TB_Repuls%m
   mc => TB_Repuls%mc
   dc => TB_Repuls%dc
   c0 => TB_Repuls%c0(1)
   c1 => TB_Repuls%c0(2)
   c2 => TB_Repuls%c0(3)
   c3 => TB_Repuls%c0(4)
   d1 => TB_Repuls%d1
   dm => TB_Repuls%dm
   x_d1 = a_x-d1

   if (a_x .GT. dm) then
      dphi = 0.0d0
   elseif (a_x .GT. d1) then
      dphi = c1 + c2*x_d1*2.0d0 + c3*x_d1*x_d1*3.0d0
   elseif (a_x .LE. 0.0d0) then
      dphi = 0.0d0
   else
      if (present(Phi_precalc)) then ! Phi function is already known:
         dphi = -m*Phi_precalc/a_x*(1.0d0 + mc*(a_x/dc)**mc)
      else ! call Phi function:
         dphi = -m*phi(TB_Repuls,a_x)/a_x*(1.0d0 + mc*(a_x/dc)**mc)
      endif
   endif
   nullify(phi0, d0, m, mc, dc, c0, c1, c2, c3, d1, dm)
end function dphi


! The derivative of the function by its argument. To be multiplied later by the derivative or the argument by the variable:
! d2 Phi / d r2_{ij},  thus, total derivative by the variable is: 
! d2 Phi / d r2_{k,beta} = d2 Phi / d r2_{ij} * (d r_{ij} / d r_{k,beta})^2 + d Phi / d r_{ij} * d2 r_{ij} / d r2_{k,beta})
function d2phi(TB_Repuls, a_x, phi_r, dphi)	! second derevative of the repulsive potential itself
   type(TB_Rep_Pettifor), intent(in), target :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_x, phi_r, dphi 
   real(8) :: d2phi
   !-----------------
   real(8) :: x_d1
   real(8), pointer :: phi0, d0, m, mc, dc
   real(8), pointer :: c0, c1, c2, c3, d1, dm
   phi0 => TB_Repuls%phi0
   d0 => TB_Repuls%d0
   m => TB_Repuls%m
   mc => TB_Repuls%mc
   dc => TB_Repuls%dc
   c0 => TB_Repuls%c0(1)
   c1 => TB_Repuls%c0(2)
   c2 => TB_Repuls%c0(3)
   c3 => TB_Repuls%c0(4)
   d1 => TB_Repuls%d1
   dm => TB_Repuls%dm
   
   if ((a_x .GT. dm) .or. (phi_r <= 0.0d0)) then
      d2phi = 0.0d0
   elseif (a_x .GT. d1) then
      x_d1 = a_x-d1
      d2phi = c2*2.0d0 + c3*x_d1*6.0d0
   elseif (a_x .LE. 0.0d0) then
      d2phi = 0.0d0
   else
      d2phi = (-1.0d0/a_x + dphi/phi_r)*dphi - m/a_x*phi_r*mc*mc/dc*(a_x/dc)**(mc-1)
   endif
   nullify(phi0, d0, m, mc, dc, c0, c1, c2, c3, d1, dm)
end function d2phi


function dfx(TB_Repuls,x) ! derivative polinomial form, for repulsive potential
   ! uses Eq. 2.36, Page 40 in H.Jeschke PhD thesis.
   type(TB_Rep_Pettifor), intent(in), target :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8) x
   real(8) dfx
   !----------
   real(8), pointer :: a0, a1, a2, a3, a4
   real(8) :: x2
   a0 => TB_Repuls%a0(1)
   a1 => TB_Repuls%a0(2)
   a2 => TB_Repuls%a0(3)
   a3 => TB_Repuls%a0(4)
   a4 => TB_Repuls%a0(5)
   x2 = x*x
   dfx = a1 + a2*x*2.0d0 + a3*x2*3.0d0 + a4*x*x2*4.0d0
   nullify(a0, a1, a2, a3, a4)
end function dfx


function df2x(TB_Repuls,x) ! derivative polinomial form, for repulsive potential
   ! Uses a derivative of Eq. 2.36, Page 40 in H.Jeschke PhD thesis
   type(TB_Rep_Pettifor), intent(in), target :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8) x
   real(8) df2x
   !----------
   real(8), pointer :: a2, a3, a4
   real(8) :: x2
   a2 => TB_Repuls%a0(3)
   a3 => TB_Repuls%a0(4)
   a4 => TB_Repuls%a0(5)
   x2 = x*x
   df2x = a2*2.0d0 + a3*x*6.0d0 + a4*x2*12.0d0
   nullify(a2, a3, a4)
end function df2x



END MODULE TB_Pettifor
