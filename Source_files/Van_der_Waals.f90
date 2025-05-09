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
! This module contains subroutines to deal with van der Waals forces within TB

MODULE Van_der_Waals
use Universal_constants
use Objects
use Little_subroutines, only : fast_pow, count_3d, Fermi_function, d_Fermi_function, d2_Fermi_function
use Atomic_tools, only : get_interplane_indices, shortest_distance, get_near_neighbours, get_number_of_image_cells, &
            distance_to_given_cell
#ifdef MPI_USED
   use MPI_subroutines, only : do_MPI_Allreduce
#endif

implicit none
PRIVATE

public :: Construct_B, get_vdW_s, get_vdW_s_D, get_vdW_interlayer, d_vdW_forces

 contains


subroutine test_vdW(TB_Waals)
   class(TB_vdW), dimension(:,:), allocatable, intent(inout):: TB_Waals ! van der Waals parameters within TB
   if (allocated(TB_Waals)) then
      select type(TB_Waals)
      type is (TB_vdW_Girifalco)
         print*, 'type is TB_vdW_Girifalco'
      type is (TB_vdW_Dumitrica)
         print*, 'type is TB_vdW_Dumitrica'
      type is (TB_vdW_LJ_cut)
         print*, 'type is TB_vdW_LJ_cut'
      end select
      pause 'test_vdW'
   else
      print*, 'For this material vdW class is undefined!'
      pause 'test_vdW'
   endif
end subroutine test_vdW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Girifalco type of van der Waals forces:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine get_vdW_interlayer(TB_Waals, Scell, NSC, matter, numpar, a)   ! vdW energy
! This subroutine is only used for comparison of the interlayer vdW energy with other works
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_vdW_Girifalco), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter	! all material parameters
   real(8), intent(out) :: a
   !=====================================================
   real(8), dimension(Scell(NSC)%Na) :: indices ! working array of indices
   real(8) :: sum_a, a_r
   INTEGER(4) i1, j1, m, atom_2
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   sum_a = 0.0d0
   
   ! Get in which plane each atom recides:
   call get_interplane_indices(Scell, NSC, numpar, Scell(NSC)%MDatoms, matter, indices) ! module "Atomic_tools"
!    print*, 'indices', indices
!    pause 'indices'

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Scell(NSC)%Na
   ! Do the cycle (parallel) calculations:
   do i1 = Nstart, Nend, N_incr  ! each process does its own part
   !do i1 = 1, Scell(NSC)%Na ! all atoms
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if ( (j1 /= i1) .and. (indices(i1) /= indices(j1)) ) then ! count only interplane energy:
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
            sum_a = sum_a + vdW_Girifalco(Scell, NSC, TB_Waals, i1, j1, a_r)    ! function below
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_vdW_interlayer:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'sum_a', sum_a) ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp PARALLEL private(i1,j1,m,atom_2,a_r)
   !$omp do reduction( + : sum_a)
   do i1 = 1, Scell(NSC)%Na ! all atoms
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
      !do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
         if ( (j1 /= i1) .and. (indices(i1) /= indices(j1)) ) then ! count only interplane energy:
         !if ((j1 .NE. i1) .and. (Scell(NSC)%MDatoms(i1)%S(3) /= Scell(NSC)%MDatoms(j1)%S(3))) then ! for testing only! Plotting interlayer energy.
            !a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
            sum_a = sum_a + vdW_Girifalco(Scell, NSC, TB_Waals, i1, j1, a_r)    ! function below
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1
   !$omp end do
   !$omp end parallel
#endif
   a = sum_a
end subroutine get_vdW_interlayer


subroutine get_mirror_cell_num(Scell, NSC, numpar, atoms, Nx, Ny, Nz)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !type(TB_vdW_Girifalco), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   integer, intent(out) :: Nx, Ny, Nz ! number of super-cells images that contribute to total energy
   real(8), dimension(3) :: zb ! super cell indices
   real(8) :: R_cut

   ! Get the cut-off distance for vdW potential:
   call get_near_neighbours(Scell, numpar, include_vdW=.true., dm=R_cut) ! module "Atomic_tools"
   
   call get_number_of_image_cells(Scell, NSC, atoms, R_cut, Nx, Ny, Nz) ! module "Atomic_tools"

   ! It has to be at least one image cell:
   if (Nx < 1) Nx = 1
   if (Ny < 1) Ny = 1
   if (Nz < 1) Nz = 1
!    !Test
!    Nx = 0
!    Ny = 0
!    Nz = 0
end subroutine get_mirror_cell_num


subroutine get_vdW_s(TB_Waals, Scell, NSC, numpar, a)   ! vdW energy
! This subroutine calcualtes the full vdW energy
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !type(TB_vdW_Girifalco), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   class(TB_vdW), dimension(:,:), allocatable, intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   real(8) :: sum_a, a_r, V_vdW, E_pot, E_pot_one
   real(8), dimension(Scell(1)%Na) :: E_pot_array
   integer :: Nx, Ny, Nz, zb(3)
   INTEGER(4) i1, j1, m, atom_2, x_cell, y_cell, z_cell
   logical :: origin_cell
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part
   
   ! Find how many image cells along each direction we need to include:
   call get_mirror_cell_num(Scell, NSC, numpar, Scell(NSC)%MDatoms, Nx, Ny, Nz) ! subroutine above
   
   sum_a = 0.0d0
   ! Construct matrix of all the radial functions for each pair of atoms:
#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Scell(NSC)%Na
   E_pot_array = 0.0d0  ! initializing
   ! Do the cycle (parallel) calculations:
   do i1 = Nstart, Nend, N_incr  ! each process does its own part
   !do i1 = 1, Scell(NSC)%Na ! all atoms
      XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
         YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
            ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
               zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
               origin_cell = ALL(zb==0) ! if it is the origin cell

               do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                     !call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, a_r) ! module "Atomic_tools"
                     !sum_a = sum_a + vdW_Girifalco(Scell, NSC, TB_Waals, i1, j1, a_r)    ! function below
                     E_pot_one = vdW_energy(Scell, NSC, TB_Waals, i1, j1, a_r)    ! function below
                     sum_a = sum_a + E_pot_one
                     E_pot_array(i1) = E_pot_array(i1) + E_pot_one
                  endif ! (j1 .NE. i1)
               enddo ! j1
            enddo ZC2
         enddo YC2
      enddo XC2
   enddo ! i1
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_vdW_s:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'sum_a', sum_a) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'E_pot_array', E_pot_array) ! module "MPI_subroutines"

   ! And save for each atom:
   Scell(NSC)%MDAtoms(:)%Epot = Scell(NSC)%MDAtoms(:)%Epot + E_pot_array(:)*0.5d0 ! exclude double-counting

#else ! use OpenMP instead
   !$omp PARALLEL private(i1,j1,a_r,x_cell,y_cell,z_cell,zb,origin_cell, E_pot, E_pot_one) shared(NSC)
   !$omp do reduction( + : sum_a)
   XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
      YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
         ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
            zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
            origin_cell = ALL(zb==0) ! if it is the origin cell
            do i1 = 1, Scell(NSC)%Na ! all atoms
               E_pot = 0.0d0  ! to restart
               do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                     !call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, a_r) ! module "Atomic_tools"
                     !sum_a = sum_a + vdW_Girifalco(Scell, NSC, TB_Waals, i1, j1, a_r)    ! function below
                     E_pot_one = vdW_energy(Scell, NSC, TB_Waals, i1, j1, a_r)    ! function below
                     sum_a = sum_a + E_pot_one
                     E_pot = E_pot + E_pot_one
                  endif ! (j1 .NE. i1)
               enddo ! j1
               ! And save for each atom:
               Scell(NSC)%MDAtoms(i1)%Epot = Scell(NSC)%MDAtoms(i1)%Epot + E_pot*0.5d0 ! exclude double-counting
            enddo ! i1
         enddo ZC
      enddo YC
   enddo XC
   !$omp end do
   !$omp end parallel
#endif
   a = sum_a * 0.5d0 ! compensate for the double counting
   !pause 'get_vdW_s pause'
end subroutine get_vdW_s


! Construct matrices of multipliers often used later to calculate derivatives of vdW:
subroutine Construct_B(TB_Waals, Scell, NSC, numpar, atoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, &
                       Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz)
   !type(TB_vdW_Girifalco), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   class(TB_vdW), dimension(:,:), allocatable, intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   real(8), dimension(:,:,:), allocatable, intent(out) :: Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij ! intermediate calculations matrices
   integer, intent(out) :: Nx, Ny, Nz ! number of super-cell images to consider
   !-----------------------------
   real(8) :: x, y, z, sx, sy, sz
   integer :: x_cell, y_cell, z_cell, coun_cell
   integer :: i1, j1, n

   n = size(atoms) ! total number of atoms

   ! Find how many image cells along each direction we need to include:
   call get_mirror_cell_num(Scell, NSC, numpar, atoms, Nx, Ny, Nz) ! subroutine above

   coun_cell = (2*Nx+1)*(2*Ny+1)*(2*Nz+1) ! total number of image cells

   if (.not.allocated(Bij))   allocate(Bij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(A_rij)) allocate(A_rij(coun_cell,n,n), source=1.0d30)
   if (.not.allocated(Xij))   allocate(Xij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(Yij))   allocate(Yij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(Zij))   allocate(Zij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(SXij))   allocate(SXij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(SYij))   allocate(SYij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(SZij))   allocate(SZij(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(XijSupce))   allocate(XijSupce(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(YijSupce))   allocate(YijSupce(coun_cell,n,n), source=0.0d0)
   if (.not.allocated(ZijSupce))   allocate(ZijSupce(coun_cell,n,n), source=0.0d0)

   ! Fill the arrays that were already allocated:
   call cell_cycle_B(Scell, NSC, TB_Waals, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz, n) ! below

end subroutine Construct_B


subroutine cell_cycle_B(Scell, NSC, TB_Waals, numpar, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz, n)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_vdW), dimension(:,:), allocatable, intent(in) :: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(inout) :: numpar   ! all numerical parameters
   real(8), dimension(:,:,:), intent(inout) :: Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij
   integer, intent(in) :: Nx, Ny, Nz ! number of super-cell images to consider
   integer, intent(in) :: n ! number of atoms
   !-------------------------
   real(8) :: x, y, z, sx, sy, sz
   integer :: x_cell, y_cell, z_cell, coun_cell
   integer :: i1, j1
   integer, dimension(3) :: zb
   logical :: origin_cell
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

#ifdef MPI_USED   ! use the MPI version [tested]
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = n
   ! Initializing:
   Xij = 0.0d0
   Yij = 0.0d0
   Zij = 0.0d0
   SXij = 0.0d0
   SYij = 0.0d0
   SZij = 0.0d0
   XijSupce = 0.0d0
   YijSupce = 0.0d0
   ZijSupce = 0.0d0
   Bij = 0.0d0

   do i1 = Nstart, Nend, N_incr  ! each process does its own part
   !do i1 = 1,n
      XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
         YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
            ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
               coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"
               zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
               origin_cell = ALL(zb==0) ! if it is the origin cell
               do j1 = 1,n
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, A_rij(coun_cell,i1,j1), &
                                                x=x, y=y, z=z, sx=sx, sy=sy, sz=sz) ! module "Atomic_tools"

                     ! Save necessary elements:
                     Xij(coun_cell,i1,j1) = x
                     Yij(coun_cell,i1,j1) = y
                     Zij(coun_cell,i1,j1) = z
                     SXij(coun_cell,i1,j1) = sx
                     SYij(coun_cell,i1,j1) = sy
                     SZij(coun_cell,i1,j1) = sz
                     XijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
                     YijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                     ZijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

                     Bij(coun_cell,i1,j1) = dvdW(TB_Waals, Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA, A_rij(coun_cell,i1,j1)) ! below
                  else ! No self-interaction
                     A_rij(coun_cell,i1,j1) = 1.0d30
                     Xij(coun_cell,i1,j1) = 0.0d0
                     Yij(coun_cell,i1,j1) = 0.0d0
                     Zij(coun_cell,i1,j1) = 0.0d0
                     Bij(coun_cell,i1,j1) = 0.0d0
                  endif
               enddo ! j1
            enddo ZC
         enddo YC
      enddo XC
   enddo ! i1

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in cell_cycle_B:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Xij', Xij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Yij', Yij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Zij', Zij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'SXij', SXij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'SYij', SYij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'SZij', SZij) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'XijSupce', XijSupce) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'YijSupce', YijSupce) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'ZijSupce', ZijSupce) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Bij', Bij) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
   !$omp PARALLEL private(i1,j1,x,y,z,sx,sy,sz,x_cell,y_cell,z_cell,zb,origin_cell,coun_cell)
   !$omp do
   XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
      YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
         ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
            coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"
            zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
            origin_cell = ALL(zb==0) ! if it is the origin cell
            do i1 = 1,n
               do j1 = 1,n
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, A_rij(coun_cell,i1,j1), &
                                                x=x, y=y, z=z, sx=sx, sy=sy, sz=sz) ! module "Atomic_tools"

                     ! Save necessary elements:
                     Xij(coun_cell,i1,j1) = x
                     Yij(coun_cell,i1,j1) = y
                     Zij(coun_cell,i1,j1) = z
                     SXij(coun_cell,i1,j1) = sx
                     SYij(coun_cell,i1,j1) = sy
                     SZij(coun_cell,i1,j1) = sz
                     XijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
                     YijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                     ZijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

                     Bij(coun_cell,i1,j1) = dvdW(TB_Waals, Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA, A_rij(coun_cell,i1,j1)) ! below
                  else ! No self-interaction
                     A_rij(coun_cell,i1,j1) = 1.0d30
                     Xij(coun_cell,i1,j1) = 0.0d0
                     Yij(coun_cell,i1,j1) = 0.0d0
                     Zij(coun_cell,i1,j1) = 0.0d0
                     Bij(coun_cell,i1,j1) = 0.0d0
                  endif
               enddo ! j1
            enddo ! i1
         enddo ZC
      enddo YC
   enddo XC
   !$omp end do
   !$omp end parallel
#endif
end subroutine cell_cycle_B




subroutine d_vdW_forces(Scell, NSC, numpar, F_vdW, dF_vdW)   ! vdW force and second derivative
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), dimension(:,:), allocatable, intent(out) :: F_vdW, dF_vdW	! force and its derivative
   !=====================================================
   real(8) :: a_r, F(3), dF(3), F_r, dF_r, x, y, z, drdrx, drdry, drdrz, d2rdr2x, d2rdr2y, d2rdr2z
   integer :: Nx, Ny, Nz, zb(3)
   INTEGER(4) :: n, i1, j1, m, atom_2, x_cell, y_cell, z_cell
   logical :: origin_cell
   integer, pointer :: KOA1, KOA2
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   n = Scell(NSC)%Na   ! Number of atoms
   ! Make sure the forces are allocated:
   if (.not.allocated(F_vdW)) allocate(F_vdW(3,n))
   if (.not.allocated(dF_vdW)) allocate(dF_vdW(3,n))
   F_vdW(:,:) = 0.0d0	! just to start with, forces
   dF_vdW(:,:) = 0.0d0	! just to start with, derivatives of forces


   ! Check if there is any vdW forces in this parameterization:
   if (.not.allocated(Scell(NSC)%TB_Waals)) return ! nothing to do, if not

   ! Find how many image cells along each direction we need to include:
   call get_mirror_cell_num(Scell, NSC, numpar, Scell(NSC)%MDatoms, Nx, Ny, Nz) ! subroutine above
   !print*, 'd_vdW_forces:', Nx, Ny, Nz


#ifdef MPI_USED   ! use the MPI version [tested]
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Scell(NSC)%Na
   do i1= Nstart, Nend, N_incr  ! each process does its own part
   !do i1 = 1, Scell(NSC)%Na ! all atoms
      KOA1 => Scell(NSC)%MDatoms(i1)%KOA   ! atom #1
      XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
         YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
            ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
               zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
               origin_cell = ALL(zb==0) ! if it is the origin cell

               F = 0.0d0   ! to restart
               dF = 0.0d0  ! to restart
               do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell

                     KOA2 => Scell(NSC)%MDatoms(j1)%KOA   ! atom #2

                     !call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, a_r, x, y, z) ! module "Atomic_tools"

                     ! Derivatives d r_{i,j} / d r_{i,alpha}:
                     drdrx = x/a_r
                     drdry = y/a_r
                     drdrz = z/a_r

                     ! Second derivatives d2 r_{ij} / d r2_{i,alpha}:
                     d2rdr2x = ddija_dria(x, a_r)  ! below
                     d2rdr2y = ddija_dria(y, a_r)  ! below
                     d2rdr2z = ddija_dria(z, a_r)  ! below

                     ! Get the derivatives of the potential by |r|:
                     F_r = dvdW(Scell(NSC)%TB_Waals, KOA1, KOA2, a_r)   ! below
                     dF_r = d2vdW(Scell(NSC)%TB_Waals, KOA1, KOA2, a_r) ! below

                     ! Construct the force and derivative:
                     F(1) = F(1) + F_r*drdrx
                     F(2) = F(2) + F_r*drdry
                     F(3) = F(3) + F_r*drdrz
                     dF(1) = dF(1) + dF_r*drdrx + F_r*d2rdr2x
                     dF(2) = dF(2) + dF_r*drdry + F_r*d2rdr2y
                     dF(3) = dF(3) + dF_r*drdrz + F_r*d2rdr2z

                  endif ! (j1 .NE. i1)
               enddo ! j1
               ! And save for each atom:
               F_vdW(:,i1) = F_vdW(:,i1) + F
               dF_vdW(:,i1) = dF_vdW(:,i1) + dF
            enddo ZC
         enddo YC
      enddo XC
   enddo ! i1
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in d_vdW_forces:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'F_vdW', F_vdW) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'dF_vdW', dF_vdW) ! module "MPI_subroutines"

#else    ! OpenMP to use instead
   !$omp PARALLEL private(i1, j1, a_r, x_cell, y_cell, z_cell, zb, origin_cell, x, y, z, drdrx, drdry, drdrz, d2rdr2x, d2rdr2y, d2rdr2z, KOA1, KOA2, F_r, dF_r, F, dF)
   !$omp do
   XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
      YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
         ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
            zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
            origin_cell = ALL(zb==0) ! if it is the origin cell
            do i1 = 1, Scell(NSC)%Na ! all atoms
               KOA1 => Scell(NSC)%MDatoms(i1)%KOA   ! atom #1
               F = 0.0d0   ! to restart
               dF = 0.0d0  ! to restart
               do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell

                     KOA2 => Scell(NSC)%MDatoms(j1)%KOA   ! atom #2

                     !call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, a_r, x, y, z) ! module "Atomic_tools"

                     ! Derivatives d r_{i,j} / d r_{i,alpha}:
                     drdrx = x/a_r
                     drdry = y/a_r
                     drdrz = z/a_r

                     ! Second derivatives d2 r_{ij} / d r2_{i,alpha}:
                     d2rdr2x = ddija_dria(x, a_r)  ! below
                     d2rdr2y = ddija_dria(y, a_r)  ! below
                     d2rdr2z = ddija_dria(z, a_r)  ! below

                     ! Get the derivatives of the potential by |r|:
                     F_r = dvdW(Scell(NSC)%TB_Waals, KOA1, KOA2, a_r)   ! below
                     dF_r = d2vdW(Scell(NSC)%TB_Waals, KOA1, KOA2, a_r) ! below

                     ! Construct the force and derivative:
                     F(1) = F(1) + F_r*drdrx
                     F(2) = F(2) + F_r*drdry
                     F(3) = F(3) + F_r*drdrz
                     dF(1) = dF(1) + dF_r*drdrx + F_r*d2rdr2x
                     dF(2) = dF(2) + dF_r*drdry + F_r*d2rdr2y
                     dF(3) = dF(3) + dF_r*drdrz + F_r*d2rdr2z

                  endif ! (j1 .NE. i1)
               enddo ! j1
               ! And save for each atom:
               F_vdW(:,i1) = F_vdW(:,i1) + F
               dF_vdW(:,i1) = dF_vdW(:,i1) + dF
            enddo ! i1
         enddo ZC
      enddo YC
   enddo XC
   !$omp end do
   nullify(KOA1, KOA2)
   !$omp end parallel
#endif
   nullify(KOA1, KOA2)
end subroutine d_vdW_forces


pure function ddija_dria(r_a, r) result(dd)
   real(8) dd  ! second derivative of the distance
   real(8), intent(in) :: r_a, r
   dd = 1.0d0/r * (1.0d0 - (r_a/r)**2)
end function ddija_dria



function vdW_energy(Scell, NSC, TB_Waals, i1, j1, a_r) result(V)   ! vdW potential in various forms
   real(8) V ! van der Waals potential energy [eV]
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   class(TB_vdW), dimension(:,:), allocatable, intent(in) :: TB_Waals ! van der Waals parameters within TB
   integer, intent(in) :: i1, j1 ! atoms
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !====================================
   integer, pointer :: KOAi, KOAj

   select type (TB_Waals)
   type is (TB_vdW_Girifalco)
      V = vdW_Girifalco(Scell, NSC, TB_Waals, i1, j1, a_r)  ! below

   type is (TB_vdW_LJ_cut)
      KOAi => Scell(NSC)%MDatoms(i1)%KOA ! kind of first atom
      KOAj => Scell(NSC)%MDatoms(j1)%KOA ! kind of second atom
      V = double_truncated_LJ(TB_Waals(KOAi,KOAj), a_r) ! below

   type is (TB_vdW_Dumitrica) ! UNFINISHED, DO NOT USE
      V = 0.0d0

   end select

   nullify(KOAi, KOAj)
end function vdW_energy



function dvdW(TB_Waals, KOA1, KOA2, a_r) result(dV)   ! wrapper around select_type; ASSOCIATE does not work inside OMP region
   real(8) :: dV  ! derivative of the LJ potential with cut-offs included
   class(TB_vdW), intent(in), dimension(:,:), allocatable :: TB_Waals   ! parameters of the repulsive part of TB-H
   integer, intent(in) :: KOA1, KOA2   ! kind of atom #1 and #2
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !-------------
   ! Select the vdW potential:
   select type (TB_Waals)
   type is (TB_vdW_Girifalco)
      dV = d_vdW_Girifalco(TB_Waals(KOA1,KOA2), a_r)  ! below

   type is (TB_vdW_LJ_cut)
      dV = d_double_truncated_LJ(TB_Waals(KOA1,KOA2), a_r) ! below

   type is (TB_vdW_ILJ_cut)
      dV = d_ILJ(TB_Waals(KOA1,KOA2), a_r) ! below

   type is (TB_vdW_Dumitrica) ! UNFINISHED, DO NOT USE
      dV = 0.0d0
   end select
end function dvdW



function d2vdW(TB_Waals, KOA1, KOA2, a_r) result(dV)  ! wrapper around select_type; ASSOCIATE does not work inside OMP region
   real(8) :: dV  ! derivative of the LJ potential with cut-offs included
   class(TB_vdW), intent(in), dimension(:,:), allocatable :: TB_Waals   ! parameters of the repulsive part of TB-H
   integer, intent(in) :: KOA1, KOA2   ! kind of atom #1 and #2
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !-------------
   ! Select the vdW potential:
   select type (TB_Waals)
   type is (TB_vdW_Girifalco)
      dV = d2_vdW_Girifalco(TB_Waals(KOA1,KOA2), a_r)  ! below

   type is (TB_vdW_LJ_cut)
      dV = d2_double_truncated_LJ(TB_Waals(KOA1,KOA2), a_r) ! below

   type is (TB_vdW_ILJ_cut)
      dV = d2_ILJ(TB_Waals(KOA1,KOA2), a_r) ! below

   type is (TB_vdW_Dumitrica) ! UNFINISHED, DO NOT USE
      dV = 0.0d0
   end select
end function d2vdW



pure function double_truncated_LJ(TB_Waals, a_r) result(V)
   real(8) V   ! [eV] potential
   type(TB_vdW_LJ_cut), intent(in) :: TB_Waals   ! LJ parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !----------------------------
   real(8) :: V_LJ, f_short, f_long
   ! Pure LJ potential:
   V_LJ = General_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, a_r) ! below
   ! Long-range cutoff:
   f_long = large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   ! Short-range cutoff:
   f_short = short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   ! collect terms:
   V = V_LJ * f_short * f_long   ! [eV]
   !if ((f_short < 1.0d0) .and. (f_long < 1.0d0)) print*, a_r, V_LJ, f_short, f_long
end function double_truncated_LJ


pure function d_double_truncated_LJ(TB_Waals, a_r) result(V)
   real(8) V   ! [eV] potential
   type(TB_vdW_LJ_cut), intent(in) :: TB_Waals   ! LJ parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !----------------------------
   real(8) :: V_LJ, f_short, f_long, d_V_LJ, d_f_short, d_f_long
   ! Pure LJ potential and its derivative:
   V_LJ = General_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, a_r) ! below
   d_V_LJ = d_General_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, a_r) ! below
   ! Long-range cutoff and its derivative:
   f_long = large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   d_f_long = d_large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   ! Short-range cutoff and its derivative:
   f_short = short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   d_f_short = d_short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   ! collect terms:
   V = V_LJ*(d_f_short*f_long + f_short*d_f_long) + d_V_LJ*(f_short*f_long)   ! [eV/A]
end function d_double_truncated_LJ


pure function d2_double_truncated_LJ(TB_Waals, a_r) result(V)
   real(8) V   ! [eV] potential
   type(TB_vdW_LJ_cut), intent(in) :: TB_Waals   ! LJ parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !----------------------------
   real(8) :: V_LJ, f_short, f_long,  d_V_LJ, d_f_short, d_f_long,  d2_V_LJ, d2_f_short, d2_f_long
   ! Pure LJ potential and its derivative:
   V_LJ = General_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, a_r) ! below
   d_V_LJ = d_General_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, a_r) ! below
   d2_V_LJ = d2_General_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, a_r) ! below

   ! Long-range cutoff and its derivative:
   f_long = large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   d_f_long = d_large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   d2_f_long = d2_large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below

   ! Short-range cutoff and its derivative:
   f_short = short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   d_f_short = d_short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   d2_f_short = d2_short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below

   ! collect terms:
   V = d2_V_LJ * (f_short*f_long) + &
       2.0d0*d_V_LJ * (d_f_short*f_long + f_short*d_f_long) + &
       V_LJ * (d2_f_short*f_long + 2.0d0*d_f_short*d_f_long + f_short*d2_f_long)   ! [eV/A^2]
end function d2_double_truncated_LJ




pure function ILJ(TB_Waals, a_r) result(V)
   real(8) V   ! [eV] potential
   type(TB_vdW_ILJ_cut), intent(in) :: TB_Waals   ! LJ parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !----------------------------
   real(8) :: V_LJ, f_short, f_long
   ! Pure LJ potential:
   V_LJ = Improved_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, TB_Waals%m, a_r) ! below
   ! Long-range cutoff:
   f_long = large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   ! Short-range cutoff:
   f_short = short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   ! collect terms:
   V = V_LJ * f_short * f_long   ! [eV]
   !if ((f_short < 1.0d0) .and. (f_long < 1.0d0)) print*, a_r, V_LJ, f_short, f_long
end function ILJ


pure function d_ILJ(TB_Waals, a_r) result(V)
   real(8) V   ! [eV] derivative of the potential
   type(TB_vdW_ILJ_cut), intent(in) :: TB_Waals   ! LJ parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !----------------------------
   real(8) :: V_LJ, f_short, f_long, d_V_LJ, d_f_short, d_f_long
   ! Pure LJ potential and its derivative:
   V_LJ = Improved_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, TB_Waals%m, a_r) ! below
   d_V_LJ = d_Improved_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, TB_Waals%m, a_r) ! below
   ! Long-range cutoff and its derivative:
   f_long = large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   d_f_long = d_large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   ! Short-range cutoff and its derivative:
   f_short = short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   d_f_short = d_short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   ! collect terms:
   V = V_LJ*(d_f_short*f_long + f_short*d_f_long) + d_V_LJ*(f_short*f_long)   ! [eV/A]
end function d_ILJ


pure function d2_ILJ(TB_Waals, a_r) result(V)
   real(8) V   ! [eV] second derivative of the potential
   type(TB_vdW_ILJ_cut), intent(in) :: TB_Waals   ! LJ parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !----------------------------
   real(8) :: V_LJ, f_short, f_long,  d_V_LJ, d_f_short, d_f_long,  d2_V_LJ, d2_f_short, d2_f_long
   ! Pure LJ potential and its derivatives:
   V_LJ = Improved_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, TB_Waals%m, a_r) ! below
   d_V_LJ = d_Improved_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, TB_Waals%m, a_r) ! below
   d2_V_LJ = d2_Improved_LJ_potential(TB_Waals%eps, TB_Waals%r0, TB_Waals%n, TB_Waals%m, a_r) ! below

   ! Long-range cutoff and its derivatives:
   f_long = large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   d_f_long = d_large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below
   d2_f_long = d2_large_range_cutoff(TB_Waals%d0_cut, TB_Waals%dd_cut, a_r) ! below

   ! Short-range cutoff and its derivatives:
   f_short = short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   d_f_short = d_short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below
   d2_f_short = d2_short_range_cutoff(TB_Waals%d0_short, TB_Waals%dd_short, a_r) ! below

   ! collect terms:
   V = d2_V_LJ * (f_short*f_long) + &
       2.0d0*d_V_LJ * (d_f_short*f_long + f_short*d_f_long) + &
       V_LJ * (d2_f_short*f_long + 2.0d0*d_f_short*d_f_long + f_short*d2_f_long)   ! [eV/A^2]
end function d2_ILJ



pure function General_LJ_potential(eps, r0, n, r) result(V)
   real(8) :: V   ! LJ potential in Mie form
   real(8), intent(in) :: eps, r0, n   ! [eV], {A], [-] : potential minimum, position of minimum, power
   real(8), intent(in) :: r   ! [A] interatomic distance
   !------------------
   real(8) :: r0r, r0rn
   if (r > 0.0d0) then
      r0r = r0/r
      r0rn = r0r**n
      V = eps * (r0rn - 2.0d0)*r0rn
   else
      V = 0.0d0
   endif
end function General_LJ_potential


pure function d_General_LJ_potential(eps, r0, n, r) result(V)
   real(8) :: V   ! LJ potential in Mie form
   real(8), intent(in) :: eps, r0, n   ! [eV], {A], [-] : potential minimum, position of minimum, power
   real(8), intent(in) :: r   ! [A] interatomic distance
   !------------------
   real(8) :: r0r, r0rn
   if (r > 0.0d0) then
      r0r = r0/r
      r0rn = r0r**n
      V = -eps * 2.0d0*n/r * r0rn * (r0rn - 1.0d0)
   else
      V = 0.0d0
   endif
end function d_General_LJ_potential


pure function d2_General_LJ_potential(eps, r0, n, r) result(V)
   real(8) :: V   ! LJ potential in Mie form
   real(8), intent(in) :: eps, r0, n   ! [eV], {A], [-] : potential minimum, position of minimum, power
   real(8), intent(in) :: r   ! [A] interatomic distance
   !------------------
   real(8) :: r0r, r0rn
   if (r > 0.0d0) then
      r0r = r0/r
      r0rn = r0r**n
      V = eps * 2.0d0*n/r**2 * ( (n+1.0d0)*r0rn * (r0rn - 1.0d0) + n*r0rn**2 )
   else
      V = 0.0d0
   endif
end function d2_General_LJ_potential



pure function Improved_LJ_potential(eps, r0, n, m, r) result(V)
   ! Improved Lennard-Jones: https://www.mdpi.com/1420-3049/26/13/3906
   ! V=eps*( m/(n-m)*(r0/r)^(n) - n/(n-m)*(r0/r)^m )
   ! reducing to LJ for n=12, m=6
   real(8) :: V   ! LJ potential in Mie form
   real(8), intent(in) :: eps, r0, n, m   ! [eV], {A], [-] : potential minimum, position of minimum, powers
   real(8), intent(in) :: r   ! [A] interatomic distance
   !------------------
   real(8) :: r0r, n_m
   if (r > 0.0d0) then
      r0r = r0/r
      n_m = n - m
      V = eps * (m/n_m*r0r**n - n/n_m*r0r**m)
   else
      V = 0.0d0
   endif
end function Improved_LJ_potential


pure function d_Improved_LJ_potential(eps, r0, n, m, r) result(V)
   real(8) :: V   ! LJ potential in Mie form
   real(8), intent(in) :: eps, r0, n, m  ! [eV], {A], [-] : potential minimum, position of minimum, powers
   real(8), intent(in) :: r   ! [A] interatomic distance
   !------------------
   real(8) :: r0r, n_m
   if (r > 0.0d0) then
      r0r = r0/r
      n_m = n - m
      V = -eps*(r0r**n - r0r**m)*n*m/(n_m*r)
   else
      V = 0.0d0
   endif
end function d_Improved_LJ_potential


pure function d2_Improved_LJ_potential(eps, r0, n, m, r) result(V)
   real(8) :: V   ! LJ potential in Mie form
   real(8), intent(in) :: eps, r0, n, m  ! [eV], {A], [-] : potential minimum, position of minimum, powers
   real(8), intent(in) :: r   ! [A] interatomic distance
   !------------------
   real(8) :: r0r, n_m
   if (r > 0.0d0) then
      r0r = r0/r
      n_m = n - m
      V = eps*n*m/(n_m * r**2)*((n+1.0d0)*r0r**n - (m+1.0d0)*r0r**m)
   else
      V = 0.0d0
   endif
end function d2_Improved_LJ_potential




pure function large_range_cutoff(d0, dd, r) result(fl)
   real(8) fl
   real(8), intent(in) :: d0, dd, r ! [A], cutoff distance, its width, interatomic distance
   fl = Fermi_function(d0, dd, r)   ! module "Little_sobroutine"
end function large_range_cutoff

pure function d_large_range_cutoff(d0, dd, r) result(fl)
   real(8) fl
   real(8), intent(in) :: d0, dd, r ! [A], cutoff distance, its width, interatomic distance
   fl = d_Fermi_function(d0, dd, r)   ! module "Little_sobroutine"
end function d_large_range_cutoff

pure function d2_large_range_cutoff(d0, dd, r) result(fl)
   real(8) fl
   real(8), intent(in) :: d0, dd, r ! [A], cutoff distance, its width, interatomic distance
   fl = d2_Fermi_function(d0, dd, r)   ! module "Little_sobroutine"
end function d2_large_range_cutoff


pure function short_range_cutoff(d0, dd, r) result(fsh)
   real(8) fsh
   real(8), intent(in) :: d0, dd, r ! [A], cutoff distance, its width, interatomic distance
   fsh = 1.0d0 - Fermi_function(d0, dd, r)   ! module "Little_sobroutine"
end function short_range_cutoff

pure function d_short_range_cutoff(d0, dd, r) result(fsh)
   real(8) fsh
   real(8), intent(in) :: d0, dd, r ! [A], cutoff distance, its width, interatomic distance
   fsh = -d_Fermi_function(d0, dd, r)   ! module "Little_sobroutine"
end function d_short_range_cutoff

pure function d2_short_range_cutoff(d0, dd, r) result(fsh)
   real(8) fsh
   real(8), intent(in) :: d0, dd, r ! [A], cutoff distance, its width, interatomic distance
   fsh = -d2_Fermi_function(d0, dd, r)   ! module "Little_sobroutine"
end function d2_short_range_cutoff





function vdW_Girifalco(Scell, NSC, TB_Waals, i1, j1, a_r) ! function according to Girifalco, based on lenard-Jones potential + cut-offs
   real(8) :: vdW_Girifalco ! van der Waals potential energy [eV]
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_vdW_Girifalco), dimension(:,:), intent(in), target :: TB_Waals ! van der Waals parameters within TB
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   integer, intent(in) :: i1, j1 ! atoms
   !====================================
   real(8), pointer :: C12	! [eV*A^12] Lenard-Jones C12
   real(8), pointer :: C6	! [eV*A^6] Lenard-Jones C6
   real(8), pointer :: dm 	! [A] radius where to switch to polinomial
   real(8), pointer :: d_cut 	! [A] cut-off radius at large r
   real(8), pointer :: a	! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8), pointer :: b 
   real(8), pointer :: c
   real(8), pointer :: d
   real(8), pointer :: dsm	! [A] radius where to switch to polinomial
   real(8), pointer :: ds_cut	! [A] cut-off radius at small r
   real(8), pointer :: as	! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8), pointer :: bs
   real(8), pointer :: cs
   real(8), pointer :: ds
   real(8), pointer :: es
   real(8), pointer :: fs
   integer, pointer :: KOAi, KOAj
   !====================================   
   real(8) :: a_r2, a_r3, a_r4, r
   !real(8) :: f_cut_small, small_f, f_cut_large
   
   KOAi => Scell(NSC)%MDatoms(i1)%KOA ! kind of first atom
   KOAj => Scell(NSC)%MDatoms(j1)%KOA ! kind of second atom
   C6 => TB_Waals(KOAi, KOAj)%C6
   C12 => TB_Waals(KOAi, KOAj)%C12
   dm => TB_Waals(KOAi, KOAj)%dm
   d_cut => TB_Waals(KOAi, KOAj)%d_cut
   a => TB_Waals(KOAi, KOAj)%a
   b => TB_Waals(KOAi, KOAj)%b
   c => TB_Waals(KOAi, KOAj)%c
   d => TB_Waals(KOAi, KOAj)%d
   dsm => TB_Waals(KOAi, KOAj)%dsm
   ds_cut => TB_Waals(KOAi, KOAj)%ds_cut
   as => TB_Waals(KOAi, KOAj)%as
   bs => TB_Waals(KOAi, KOAj)%bs
   cs => TB_Waals(KOAi, KOAj)%cs
   ds => TB_Waals(KOAi, KOAj)%ds
   es => TB_Waals(KOAi, KOAj)%es
   fs => TB_Waals(KOAi, KOAj)%fs
   
   if ((a_r > d_cut) .or. (a_r < ds_cut)) then ! anything beyond cut-offs is zero:
      vdW_Girifalco = 0.0d0
   elseif (a_r > dm) then ! in this region we use polinomial function:
      a_r2 = a_r*a_r
      vdW_Girifalco = a*a_r2*a_r + b*a_r2 + c*a_r + d
   elseif (a_r < dsm) then ! in this region we use polinomial function:
      a_r2 = a_r*a_r
      a_r4 = a_r2*a_r2
      vdW_Girifalco = as*a_r4*a_r + bs*a_r4 + cs*a_r2*a_r + ds*a_r2 + es*a_r + fs
   else ! at shorter distances we use proper potential:
      r = a_r !+ r_LJ
      vdW_Girifalco = LJ_potential(C12, C6, r)
   endif

   nullify(KOAi, KOAj, C6, C12, dm, d_cut, a, b, c, d, dsm, ds_cut, as, bs, cs, ds, es, fs)
end function vdW_Girifalco


function d_vdW_Girifalco(TB_Waals, a_r) result(dvdW)
   real(8) :: dvdW ! derivative of the LJ potential with cut-offs included
   type(TB_vdW_Girifalco), intent(in), target :: TB_Waals   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   real(8), pointer :: C12	! [eV*A^12] Lenard-Jones C12
   real(8), pointer :: C6	! [eV*A^6] Lenard-Jones C6
   real(8), pointer :: dm 	! [A] radius where to switch to polinomial
   real(8), pointer :: d_cut 	! [A] cut-off radius
   real(8), pointer :: a		! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8), pointer :: b 
   real(8), pointer :: c
   real(8), pointer :: d
   real(8), pointer :: dsm	! [A] radius where to switch to polinomial
   real(8), pointer :: ds_cut	! [A] cut-off radius at small r
   real(8), pointer :: as	! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8), pointer :: bs
   real(8), pointer :: cs
   real(8), pointer :: ds
   real(8), pointer :: es
   real(8), pointer :: fs
   !====================================
   real(8) :: r, E_vdW, f_cut_small, f_cut_large, d_f_small, d_f_large, d_E_vdW
   real(8) :: a_r2, a_r4
   
   C6 => TB_Waals%C6
   C12 => TB_Waals%C12
   dm => TB_Waals%dm
   d_cut => TB_Waals%d_cut
   a => TB_Waals%a
   b => TB_Waals%b
   c => TB_Waals%c
   d => TB_Waals%d
   dsm => TB_Waals%dsm
   ds_cut => TB_Waals%ds_cut
   as => TB_Waals%as
   bs => TB_Waals%bs
   cs => TB_Waals%cs
   ds => TB_Waals%ds
   es => TB_Waals%es
   fs => TB_Waals%fs
   
   if ((a_r > d_cut) .or. (a_r < ds_cut)) then ! anything beyond cut-offs is zero:
      dvdW = 0.0d0
   elseif (a_r > dm) then ! in this region we use polinomial function:
      dvdW = 3.0d0*a*a_r*a_r + 2.0d0*b*a_r + c
   elseif (a_r < dsm) then ! in this region we use polinomial function:
      a_r2 = a_r*a_r
      a_r4 = a_r2*a_r2
      dvdW = 5.0d0*as*a_r4 + 4.0d0*bs*a_r2*a_r + 3.0d0*cs*a_r2 + 2.0d0*ds*a_r + es
   else ! at shorter distances we use proper potential:
      r = a_r !+ r_LJ
      d_E_vdW = d_LJ_potential(C12, C6, r) ! derivative of  LJ part of the potential
      dvdW = d_E_vdW
   endif
   nullify(C12, C6, dm, d_cut, a, b, c, d, dsm, ds_cut, as, bs, cs, ds, es, fs)
end function d_vdW_Girifalco



function d2_vdW_Girifalco(TB_Waals, a_r) result(d2vdW)
   real(8) :: d2vdW ! second derivative of the LJ potential with cut-offs included
   type(TB_vdW_Girifalco), intent(in), target :: TB_Waals   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   real(8), pointer :: C12	! [eV*A^12] Lenard-Jones C12
   real(8), pointer :: C6	! [eV*A^6] Lenard-Jones C6
   real(8), pointer :: dm 	! [A] radius where to switch to polinomial
   real(8), pointer :: d_cut 	! [A] cut-off radius
   real(8), pointer :: a		! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8), pointer :: b
   real(8), pointer :: c
   real(8), pointer :: d
   real(8), pointer :: dsm	! [A] radius where to switch to polinomial
   real(8), pointer :: ds_cut	! [A] cut-off radius at small r
   real(8), pointer :: as	! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8), pointer :: bs
   real(8), pointer :: cs
   real(8), pointer :: ds
   real(8), pointer :: es
   real(8), pointer :: fs
   !====================================
   real(8) :: r, E_vdW, f_cut_small, f_cut_large, d_f_small, d_f_large
   real(8) :: a_r2, a_r3

   C6 => TB_Waals%C6
   C12 => TB_Waals%C12
   dm => TB_Waals%dm
   d_cut => TB_Waals%d_cut
   a => TB_Waals%a
   b => TB_Waals%b
   c => TB_Waals%c
   d => TB_Waals%d
   dsm => TB_Waals%dsm
   ds_cut => TB_Waals%ds_cut
   as => TB_Waals%as
   bs => TB_Waals%bs
   cs => TB_Waals%cs
   ds => TB_Waals%ds
   es => TB_Waals%es
   fs => TB_Waals%fs

   if ((a_r > d_cut) .or. (a_r < ds_cut)) then ! anything beyond cut-offs is zero:
      d2vdW = 0.0d0
   elseif (a_r > dm) then ! in this region we use polinomial function:
      d2vdW = 6.0d0*a*a_r + 2.0d0*b
   elseif (a_r < dsm) then ! in this region we use polinomial function:
      a_r2 = a_r*a_r
      a_r3 = a_r2*a_r
      d2vdW = 20.0d0*as*a_r3 + 12.0d0*bs*a_r2 + 6.0d0*cs*a_r + 2.0d0*ds
   else ! at shorter distances we use proper potential:
      r = a_r !+ r_LJ
      d2vdW = d2_LJ_potential(C12, C6, r) ! below
   endif
   nullify(C12, C6, dm, d_cut, a, b, c, d, dsm, ds_cut, as, bs, cs, ds, es, fs)
end function d2_vdW_Girifalco



function LJ_potential(C12, C6, r) ! Lenard-Jones potential
   real(8) :: LJ_potential
   real(8), intent(in) :: C12, C6, r
   real(8) :: a_r6, a_r12
   a_r6 = fast_pow(r,6) ! function from module "Little_subroutines"
   a_r12 = a_r6*a_r6
   LJ_potential = C12/a_r12 - C6/a_r6
end function LJ_potential


function d_LJ_potential(C12, C6, r) ! derivative of Lenard-Jones potential
   real(8) :: d_LJ_potential
   real(8), intent(in) :: C12, C6, r
   real(8) :: a_r6, a_r7, a_r13
   a_r6 = fast_pow(r,6) ! function from module "Little_subroutines"
   a_r7 = a_r6*r
   a_r13 = a_r6*a_r7
   d_LJ_potential = -12.0d0*C12/a_r13 + 6.0d0*C6/a_r7
end function d_LJ_potential


function d2_LJ_potential(C12, C6, r) result(d2F) ! second derivative of Lenard-Jones potential
   real(8) :: d2F
   real(8), intent(in) :: C12, C6, r
   real(8) :: a_r8, a_r14
   a_r8 = r**8
   a_r14 = r**14
   d2F = 156.0d0*C12/a_r14 - 42.0d0*C6/a_r8
end function d2_LJ_potential




pure function f_cut_L(a_r, r_L, d_L) ! cut-off function at large distances
   real(8) :: f_cut_L
   real(8), intent(in) :: a_r, r_L, d_L
   f_cut_L = 1.0d0/(1.0d0 + dexp((a_r - r_L)/d_L))
end function f_cut_L


pure function d_f_cut_L(a_r, r_L, d_L) ! derivative of cut-off function at large distances
   real(8) :: d_f_cut_L
   real(8), intent(in) :: a_r, r_L, d_L
   !------------------
   real(8) :: exp_r, exp_r2, arg
   arg = (a_r - r_L)/d_L
   if (arg >= log(HUGE(a_r))) then
      d_f_cut_L = 0.0d0
   else
      exp_r = dexp(arg)
      exp_r2 = 1.0d0 + exp_r
      d_f_cut_L = -exp_r/(d_L*exp_r2*exp_r2)
   endif
end function d_f_cut_L


pure function d2_f_cut_L(a_r, r_L, d_L) result(d2F) ! second derivative of cut-off function at large distances
   real(8) :: d2F
   real(8), intent(in) :: a_r, r_L, d_L
   !------------------
   real(8) :: exp_r, exp_r2, arg
   arg = (a_r - r_L)/d_L
   if (arg >= log(HUGE(a_r))) then
      d2F = 0.0d0
   else
      exp_r = dexp(arg)
      exp_r2 = 1.0d0 + exp_r
      d2F = -exp_r/(d_L**2 * exp_r2**2) * (2.0d0*exp_r/exp_r2 - 1.0d0)
   endif
end function d2_f_cut_L




pure function f_cut_S(a_r, r_S, d_S) ! cut-off function at small distances
   real(8) :: f_cut_S
   real(8), intent(in) :: a_r, r_S, d_S
   real(8) :: small_f
   small_f = 1.0d0 + dexp(-(a_r - r_S)/d_S)
   f_cut_S = 1.0d0/(small_f*small_f)
end function f_cut_S


pure function d_f_cut_S(a_r, r_S, d_S) ! derivative of cut-off function at small distances
   real(8) :: d_f_cut_S
   real(8), intent(in) :: a_r, r_S, d_S
   real(8) :: small_f, exp_r, arg
   arg = -(a_r - r_S)/d_S
   if (arg >= log(HUGE(a_r))) then
      d_f_cut_S = 0.0d0
   else
      exp_r = dexp(arg)
      small_f = 1.0d0 + exp_r
      d_f_cut_S = 2.0d0*exp_r/(d_S*small_f*small_f*small_f)
   endif
end function d_f_cut_S


pure function d2_f_cut_S(a_r, r_S, d_S) result(d2F) ! second derivative of cut-off function at small distances
   real(8) :: d2F
   real(8), intent(in) :: a_r, r_S, d_S
   real(8) :: small_f, exp_r, arg
   arg = -(a_r - r_S)/d_S
   if (arg >= log(HUGE(a_r))) then
      d2F = 0.0d0
   else
      exp_r = dexp(arg)
      small_f = 1.0d0 + exp_r
      d2F = 2.0d0/(d_S**2) * exp_r/(small_f**3) * (3.0d0*exp_r/small_f - 1.0d0)
   endif
end function d2_f_cut_S





!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! Obsolete and outdated subroutines

! Derivatives of the vdW energy by s:
subroutine dvdWdr_s(atoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, Nx, Ny, Nz)
   !type(TB_vdW), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: Bij, A_rij, Xij, Yij, Zij ! intermediate calculations matrices
   integer, intent(in) :: Nx, Ny, Nz ! number of super-cells to go through
   !---------------------
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   REAL(8), DIMENSION(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, x, y, z, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, j1, ian, dik, djk, n, m, atom_2
   integer :: x_cell, y_cell, z_cell, coun_cell
   integer, DIMENSION(3) :: zb
   real(8), dimension(:,:), allocatable :: Erx_s
   logical :: origin_cell
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   n = size(atoms) ! total number of atoms
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0


#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = n
   ! Do the cycle (parallel) calculations:
   do ian = Nstart, Nend, N_incr  ! each process does its own part
   !do ian = 1, n  ! Forces for all atoms
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
                           b = Bij(coun_cell,i1,j1) ! dvdW(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))
                           a_r = A_rij(coun_cell,i1,j1) ! call shortest_distance(Scell, NSC, atoms, i1, j1, A_rij(coun_cell,i1,j1), x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)

                           ddlta = dble(dik - djk)/a_r
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
   enddo ! ian
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in dvdWdr_s:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{Erx_s}', Erx_s) ! module "MPI_subroutines"

   do ian = 1, n  ! Forces for all atoms
      ! Add van der Waals force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)
   enddo

#else ! use OpenMP instead
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
                           b = Bij(coun_cell,i1,j1) ! dvdW(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))
                           a_r = A_rij(coun_cell,i1,j1) ! call shortest_distance(Scell, NSC, atoms, i1, j1, A_rij(coun_cell,i1,j1), x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz) 

                           ddlta = dble(dik - djk)/a_r
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
!       print*, 'OLD', ian, Erx_s(:,ian)
   enddo ! ian
   !$omp end do
   !$omp end parallel
#endif

   deallocate(Erx_s)
!    pause 'dvdWdr_s'
END subroutine dvdWdr_s



! Derivatives of the vdW energy by h:
subroutine dvdWdr_Pressure_s(atoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz)
   !type(TB_vdW_Girifalco), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
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
                        dpsy = Bij(coun_cell,i,j)    ! dvdW(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))

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
!       print*, 'OLD P', PForce
   endif
!    pause 'dvWdr_Pressure_s'
end subroutine dvdWdr_Pressure_s


subroutine dvdWdr_s_OLD(TB_Waals, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   !type(TB_Rep_Molteni), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   class(TB_vdW), dimension(:,:), allocatable, intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   !-----------------------
   REAL(8), DIMENSION(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, x, y, z, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, j1, ian, dik, djk, n, m, atom_2, NumTB
   real(8), DIMENSION(3) :: zb
   real(8), dimension(:,:), allocatable :: Erx_s
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   n = size(atoms)
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = n
   ! Do the cycle (parallel) calculations:
   do ian = Nstart, Nend, N_incr  ! each process does its own part
   !do ian = 1, n  ! Forces for all atoms
     !Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
     do i1 = 1, n
         dpsi = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i1)
         do atom_2 = 1,m ! do only for atoms close to that one
            j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
            if (j1 .GT. 0) then
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

!                   x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
!                   x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
!                   x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
                  x1(2) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
                  x1(3) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)

!                   b = dphi_M(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r,1)
                  b = dvdW(TB_Waals, Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA, a_r)

                  ddlta = real(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
                endif
            endif
         enddo ! j1

         Erx_s(1,ian) = Erx_s(1,ian) + dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      ! Add van der Waals force to already calculated other forces:
      !Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian) !*0.5d0
      print*, 'OLD',ian, Erx_s(:,ian)
   enddo ! ian

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in dvdWdr_s_OLD:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{Erx_s}', Erx_s) ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp PARALLEL DO private(ian,i1,dpsi,m,j1,x,y,z,a_r,dik,djk,x1,b,ddlta,b_delta,atom_2)
   do ian = 1, n  ! Forces for all atoms
     !Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
     do i1 = 1, n
         dpsi = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i1)
         do atom_2 = 1,m ! do only for atoms close to that one
            j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
            if (j1 .GT. 0) then
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

!                   x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
!                   x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
!                   x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
                  x1(2) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
                  x1(3) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)

!                   b = dphi_M(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r,1)
                  b = dvdW(TB_Waals, Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA, a_r)

                  ddlta = real(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
                endif
            endif
         enddo ! j1

         Erx_s(1,ian) = Erx_s(1,ian) + dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      ! Add van der Waals force to already calculated other forces:
      !Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian) !*0.5d0
      print*, 'OLD',ian, Erx_s(:,ian)
   enddo ! ian
   !$OMP END PARALLEL DO
#endif
   deallocate(Erx_s)
END subroutine dvdWdr_s_OLD


subroutine dvdWdr_Pressure_s_OLD(TB_Waals, atoms, Scell, NSC, numpar)! derivatives of the repulsive energy by h
!    type(TB_Rep_Molteni), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   class(TB_vdW), dimension(:,:), allocatable, intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms ! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   REAL(8), DIMENSION(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i,j,k,l, i1, i2, i3, ik, m, atom_2, n, NumTB
   real(8) x,y,z,sx,sy,sz,r, rcur(3), scur(3), PForce(3,3)
   real(8) df_psy, psi, dpsy

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = size(atoms)
!       Scell(NSC)%SCforce%rep = 0.0d0
      PForce = 0.0d0 ! to start with
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i)
         do atom_2 = 1,m ! do only for atoms close to that one
            j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
            if ((i .NE. j) .AND. (j .GT. 0)) then
               !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

               dpsy = dvdW(TB_Waals, Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA, r)
               
               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! i=j
         enddo ! j

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               !Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + Rep_Pr(l,k) !*0.5d0
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)
!                print*, 'OLD P', l, k, PForce(l,k)
            enddo ! l
         enddo ! k
      enddo ! i
      !Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add vdW part to existing TB part
      print*, 'OLD P', PForce
   endif
end subroutine dvdWdr_Pressure_s_OLD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UNFINISHED PART FOR ((TB_vdW_Dumitrica))
!van der Waals force for atoms:
subroutine d_vdW_s_D(TB_Waals, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Pettifor"
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_vdW_Dumitrica), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   real(8) a_r, a_r2, a_r3, C6, sum_a, dik, djk, X, Y, Z, x1(3), dpsi(3), b, ddlta, b_delta
   INTEGER(4) i1, j1, m, atom_2, ian

   do ian = 1, Scell(NSC)%Na ! for all atoms get forces
      dpsi = 0.0d0
      do i1 = 1, Scell(NSC)%Na  ! all atoms
         !m = Scell(NSC)%Near_neighbor_size(i1)
         !do atom_2 = 1, m ! do only for atoms close to that one
         !   j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
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
            if ((j1 .NE. i1)  .OR. ((dik-djk) .NE. 0)) then
               !a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
               call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r, x1=X, y1=Y, z1=Z) ! module "Atomic_tools"
!                x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
!                x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
!                x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
               x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
               x1(2) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
               x1(3) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)
               
               a_r2 = a_r*a_r
               a_r3 = a_r2*a_r
               C6 = TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA)%C6
               b = df_damp_D(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA), a_r)*C6/(a_r3*a_r3) + f_damp_D(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA), a_r)*(-6.0d0*C6)/(a_r*a_r3*a_r3) ! function below
               ddlta = real(dik - djk)/a_r
               b_delta = b*ddlta
               dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
               dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
               dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
            endif ! (j1 .NE. i1)
         enddo ! j1
      enddo ! i1
      a = -0.5d0*sum_a
   enddo ! ian
end subroutine d_vdW_s_D


!van der Waals potential for atoms:
subroutine get_vdW_s_D(TB_Waals, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Pettifor"
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_vdW_Dumitrica), dimension(:,:), intent(in)   :: TB_Waals ! van der Waals parameters within TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   real(8) :: a_r, a_r2, a_r3, C6, sum_a, E_rep, E_rep_one, E_rep_array(Scell(NSC)%Na)
   INTEGER(4) i1, j1, m, atom_2
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   sum_a = 0.0d0

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Scell(NSC)%Na
   E_rep_array = 0.0d0
   do i1 = Nstart, Nend, N_incr  ! each process does its own part
   !do i1 = 1, Scell(NSC)%Na ! all atoms
      do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
         if (j1 .NE. i1) then
            !a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
            a_r2 = a_r*a_r
            a_r3 = a_r2*a_r
            C6 = TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA)%C6
            E_rep_one = df_damp_D(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA), a_r)*C6/(a_r3*a_r3) ! function below
            sum_a = sum_a + E_rep_one
            E_rep_array(i1) = E_rep_array(i1) + E_rep_one
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_Coulomb_scc_energy:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{sum_a}', sum_a) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{E_rep_array}', E_rep_array) ! module "MPI_subroutines"

   do i1 = 1, Scell(NSC)%Na ! all atoms
      ! And save for each atom:
      Scell(NSC)%MDAtoms(i1)%Epot = Scell(NSC)%MDAtoms(i1)%Epot - E_rep_array(i1)*0.5d0 ! exclude double-counting
   enddo
#else ! use OpenMP instead
   !$omp PARALLEL private(i1,j1,a_r,a_r2,a_r3,C6, E_rep, E_rep_one)
   !$omp do reduction( + : sum_a)
   do i1 = 1, Scell(NSC)%Na ! all atoms
      E_rep = 0.0d0
      !m = Scell(NSC)%Near_neighbor_size(i1)
      !do atom_2 = 1, m ! do only for atoms close to that one
      !   j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
      do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
         if (j1 .NE. i1) then
            !a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
            a_r2 = a_r*a_r
            a_r3 = a_r2*a_r
            C6 = TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA)%C6
            E_rep_one = df_damp_D(TB_Waals(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA), a_r)*C6/(a_r3*a_r3) ! function below
            sum_a = sum_a + E_rep_one
            E_rep = E_rep + E_rep_one
         endif ! (j1 .NE. i1)
      enddo ! j1
      ! And save for each atom:
      Scell(NSC)%MDAtoms(i1)%Epot = Scell(NSC)%MDAtoms(i1)%Epot - E_rep*0.5d0 ! exclude double-counting
   enddo ! i1
   !$omp end do
   !$omp end parallel
#endif
   a = -0.5d0*sum_a
end subroutine get_vdW_s_D


! Derivative of the damping function:
function df_damp_D(TB_Waals, r)
   real(8) :: df_damp_D ! function, dimensionless
   type(TB_vdW_Dumitrica), intent(in) :: TB_Waals ! van der Waals parameters within TB
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: a, C6, f, ar, ar2, ar3, invr
   real(8) :: k(6) ! factorials
   C6 = TB_Waals%C6
   a = TB_Waals%alpha
   k(1) = 1.0d0
   k(2) = 2.0d0
   k(3) = 6.0d0
   k(4) = 24.0d0
   k(5) = 120.0d0
   k(6) = 720.0d0
   f = a
   invr = 1.0d0/r
   ar = a*r
   ar2 = ar*ar
   ar3 = ar2*ar
   f = f + (a-invr)*ar
   f = f + (a-2.0d0*invr)*ar2/k(2)
   f = f + (a-3.0d0*invr)*ar3/k(3)
   f = f + (a-4.0d0*invr)*ar2*ar2/k(4)
   f = f + (a-5.0d0*invr)*ar3*ar2/k(5)
   f = f + (a-6.0d0*invr)*ar3*ar3/k(6)
   f = f*dexp(-ar)
   df_damp_D = f
end function df_damp_D


! Damping function:
function f_damp_D(TB_Waals, r)
   real(8) :: f_damp_D ! function, dimensionless
   type(TB_vdW_Dumitrica), intent(in) :: TB_Waals ! van der Waals parameters within TB
   real(8), intent(in) :: r ! [A] interatomic distance
   real(8) :: a, C6, f, ar, ar2, ar3
   real(8) :: k(6) ! factorials
   C6 = TB_Waals%C6
   a = TB_Waals%alpha
!    k(1) = 1.0d0
   k(2) = 2.0d0
   k(3) = 6.0d0
   k(4) = 24.0d0
   k(5) = 120.0d0
   k(6) = 720.0d0
   f = 1.0d0
   ar = a*r
   ar2 = ar*ar
   ar3 = ar2*ar
   f = f + ar
   f = f + ar2/k(2)
   f = f + ar3/k(3)
   f = f + ar2*ar2/k(4)
   f = f + ar3*ar2/k(5)
   f = f + ar3*ar3/k(6)
   f = 1.0d0 - f*dexp(-ar)
   f_damp_D = f
end function f_damp_D

END MODULE Van_der_Waals
