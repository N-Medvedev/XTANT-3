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
! This module contains subroutines to read the basis sets, and normalize the functions properly
! The basis sets can be downloaded from: https://bse.pnl.gov/bse/portal
! Now they are saved in the folder BASIS_SETS

module BS_Basis_sets

use Universal_constants
use Objects
use BS_Cartesian_Gaussians, only: Overlap_of_high_order_gaussians
use BS_Spherical_Gaussians, only: find_TM_size
use Algebra_tools, only : double_factorial

! For OpenMP external library
#ifdef OMP_inside
   USE OMP_LIB, only : omp_get_wtime
#endif



implicit none
PRIVATE

character(50) :: m_basis_sets_folder

! Where do we keep all the basis sets:
parameter(m_basis_sets_folder = '!BASIS_SETS')

public :: set_xTB_AO

 
 contains

!------------------------------------------
!  Construct overlap integrals of nonorthogonal basis set (Cartesian):
subroutine construct_S_matrix(MDAtoms, ChemEl, Sij, Sijx, Sijy, Sijz, C_transform)
   type(Atom), dimension(:), intent(in), target :: MDAtoms	! all atoms in the supercell contribute to the total basis
   type(At_data), dimension(:), intent(in), target :: ChemEl	! data for each element involved (basis sets etc.)
   real(8), dimension(:,:), intent(inout) :: Sij	! overlap matrix of nonorthogonal basis set
   real(8), dimension(:,:), intent(inout) :: Sijx	! overlap matrix of nonorthogonal basis set in X direction (UNNORMALIZED)
   real(8), dimension(:,:), intent(inout) :: Sijy	! overlap matrix of nonorthogonal basis set in Y direction (UNNORMALIZED)
   real(8), dimension(:,:), intent(inout) :: Sijz	! overlap matrix of nonorthogonal basis set in Z direction (UNNORMALIZED)
   real(8), dimension(:,:), allocatable, intent(inout) :: C_transform    ! transformation matrix from Cartesian Gausians to Spherical (Pure) Gaussians
   !------------------------------------------
   integer :: N_atoms, count1, count2
   integer, dimension(size(MDAtoms)) :: N_basis
   integer :: i1, j1, k1, N_basis1, N_GTO1
   integer :: i2, j2, k2, N_basis2, N_GTO2
   integer :: TM_siz1, TM_siz2
   real(8) :: Sij_single, Sij_single_x, Sij_single_y, Sij_single_z, C1C2
   logical :: do_transform_matrix
   !---------------------------------------------
   integer, pointer :: KOA1, KOA2, AMa(:), AMb(:)
   real(8), pointer ::  R1(:), R2(:)
   real(8), pointer :: alpha1, C1, Norm1 
   real(8), pointer :: alpha2, C2, Norm2
   !------------------------------------------

   N_atoms = size(MDAtoms)	! how many atoms in the supercell
   N_basis = 0 ! to start
   
   ! Create an array of indices to use later in parallelized region:
   count1 = 0
   do i1 = 2, N_atoms
      KOA1 => MDAtoms(i1-1)%KOA
      count1 = count1 + size(ChemEl(KOA1)%Cart_Basis)
      N_basis(i1) =  count1
   enddo
   
   ! Do we need to construct the transformation matrix?
   if (.not.allocated(C_transform)) then  ! if it does not exist yet, construct it 
      ! Define the size of transformation matrix:
      TM_siz1 =  find_TM_size(ChemEl(KOA1)%Cart_Basis(:))   ! module "BS_Spherical_Gaussians"
      KOA1 => MDAtoms(N_atoms)%KOA
      TM_siz2 = N_basis(N_atoms) + size(ChemEl(KOA1)%Cart_Basis)
      ! Allocate the transformation matrix:
      allocate(C_transform(TM_siz1,TM_siz2), source=0.0d0)
      do_transform_matrix = .true.
   else ! if it exists, reuse it
      do_transform_matrix = .false.
   endif
   

!$omp PARALLEL
!$omp do  private(i1, KOA1, R1, N_basis1, j1, count1, AMa, N_GTO1, i2, KOA2, R2, N_basis2, j2, count2, AMb, &
!$omp  N_GTO2, Sij_single, Sij_single_x, Sij_single_y, Sij_single_z, k1, alpha1, C1, Norm1, k2, alpha2, C2, Norm2, C1C2)
   do i1 = 1, N_atoms ! all atoms
      KOA1 => MDAtoms(i1)%KOA	! this atoms i is of this type KOA
      R1 => MDAtoms(i1)%R		! coordinates of this atoms
      N_basis1 = size(ChemEl(KOA1)%Cart_Basis)	! corresponding basis type for this type of atom
      
      do j1 = 1, N_basis1 ! basis functions for each atom
         !count1 = count1 + 1 ! count functions to use as indices
         count1 = N_basis(i1) + j1	 ! count functions to use as indices
!          print*, 'count1', count1, N_basis(i1) + j1
         AMa => ChemEl(KOA1)%Cart_Basis(j1)%AM ! angular momenta for this basis function

         ! contracted gaussians:
         N_GTO1 = size(ChemEl(KOA1)%Cart_Basis(j1)%GTO) ! how many primitives are used in this contracted GTO
         ! Overlap integrals with all the functions:
         do i2 = i1, N_atoms ! for only nonrepeating pairs of atoms (since the matrix is symmetric)
            KOA2 => MDAtoms(i2)%KOA	! this atoms i is if this type KOA
            R2 => MDAtoms(i2)%R			! coordinates of this atoms
            N_basis2 = size(ChemEl(KOA2)%Cart_Basis)	! corresponding basis type for this type of atom
            
            do j2 = 1, N_basis2 ! basis functions for each atom we are overlapping with
               !count2 = count2 + 1 ! count functions to use as indices
               count2 =  N_basis(i2) + j2	! count functions to use as indices
!                print*, 'count2', count2, N_basis(i2) + j2
               AMb => ChemEl(KOA2)%Cart_Basis(j2)%AM ! angular momenta for this basis function

               ! contracted gaussians:
               N_GTO2 = size(ChemEl(KOA2)%Cart_Basis(j2)%GTO) ! how many primitives are used in this contracted GTO
               Sij_single = 0.0d0 ! to start from
               !-------------------
               ! Overlap for each contracted gaussian:
               do k1 = 1, N_GTO1 ! for all contracted GTOs:
                  alpha1 => ChemEl(KOA1)%Cart_Basis(j1)%GTO(k1)%alpha
                  C1 => ChemEl(KOA1)%Cart_Basis(j1)%GTO(k1)%C
                  Norm1 => ChemEl(KOA1)%Cart_Basis(j1)%GTO(k1)%Norm
                  
                  do k2 = 1, N_GTO2 ! for all contracted GTOs:
                     alpha2 => ChemEl(KOA2)%Cart_Basis(j2)%GTO(k2)%alpha
                     C2 => ChemEl(KOA2)%Cart_Basis(j2)%GTO(k2)%C
                     Norm2 => ChemEl(KOA2)%Cart_Basis(j2)%GTO(k2)%Norm
                     
                     ! Get the overlap matrix:
                     call Overlap_of_high_order_gaussians(R1(1), R1(2), R1(3), alpha1, AMa, Norm1, R2(1), R2(2), R2(3), alpha2, AMb, Norm2, &
                                            Sij_single, Sij_single_x, Sij_single_y, Sij_single_z) ! module "BS_Cartesian_Gaussians"
                     C1C2 = C1*C2
                     Sij(count1, count2) = Sij(count1, count2) + C1C2*Sij_single	! overlap matrix of nonorthogonal basis set 
                     Sijx(count1, count2) = Sijx(count1, count2) + C1C2*Sij_single_x	! overlap matrix of nonorthogonal basis set in X direction
                     Sijy(count1, count2) = Sijy(count1, count2) + C1C2*Sij_single_y	! overlap matrix of nonorthogonal basis set in Y direction
                     Sijz(count1, count2) = Sijz(count1, count2) + C1C2*Sij_single_z	! overlap matrix of nonorthogonal basis set in Z direction
                     
                  enddo ! k2 = 1, N_GTO2
               enddo ! k = 1, N_GTO
!                write(*,'(a,i3,i3,e,i3)') 'coef', count1, count2, Sij(count1, count2), OMP_get_thread_num()
               !-------------------
               !Sij(count2, count1) = Sij(count1, count2)
            enddo ! j2 = 1, N_basis2
         enddo ! i2 = i1, N_atoms
      enddo ! j1 = 1, N_basis1
   enddo ! i1 = 1, N_atoms
!$omp end do
!$omp do  private(i1, j1)
   do i1 = 1, size(Sij,1)
      do j1 = i1, size(Sij,2)
        Sij(j1,i1) = Sij(i1,j1) ! set symmetric elements for lower triangle
        Sijx(j1,i1) = Sijx(i1,j1) ! set symmetric elements for lower triangle
        Sijy(j1,i1) = Sijy(i1,j1) ! set symmetric elements for lower triangle
        Sijz(j1,i1) = Sijz(i1,j1) ! set symmetric elements for lower triangle
      enddo
   enddo
!$omp end do
!$omp end parallel

   ! free the memory:
   nullify(KOA1, alpha1, C1, Norm1, KOA2, alpha2, C2, Norm2, AMa, AMb, R1, R2)
end subroutine construct_S_matrix 
 
 
!------------------------------------------------------------------------------
! Reading basis set from files:
subroutine find_basis_set_for_element(El_Name, Basis_name, Basis, path_sep, Error_message)
   character(*), intent(in) :: El_Name 		! chemical element name
   character(*), intent(in), optional :: Basis_name	! basis set name night be provided
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   character(1), intent(in) :: path_sep
   character(*), intent(out) :: Error_message
   !--------------------------------------------------
   character(500) :: command, folder, char1, File_name
   character(1000) :: string
   character(100), dimension(:), allocatable :: BASIS_SET_NAMES
   integer FN, Reason, N, i, N_basis_set
   logical :: file_exist, basis_exists, basis_found
   !folder = 'BASIS_SETS'
   folder = trim(adjustl(m_basis_sets_folder))
   File_name = 'temp.txt'
   !command = 'dir /b/o '//trim(adjustl(folder))//' > '//trim(adjustl(File_name))
   ! Find out which OS it is:
   ! Get the names of all data files in the folder using system commands:
   if (path_sep .EQ. '\') then	! if it is Windows
      command = 'dir '//trim(adjustl(folder))//'\*.gbs /b/o >'//trim(adjustl(File_name))
   else
      command = "ls -t"//trim(adjustl(folder))//" | grep '.gbs' >"//trim(adjustl(File_name))
   endif
   
   print*, '--------------------------------------------'
   
   ! Read all file names with basis sets into the file:
   call system(trim(adjustl(command)))
   
   FN = 9998
   inquire(file=trim(adjustl(File_name)),exist=file_exist)
   if (file_exist) then
      open(unit=FN, file='temp.txt')
   else
      print*, 'Could not read the list of basis sets:'
      print*, 'File ', trim(adjustl(File_name)), ' not found!'
   endif
   
   ! Count how many basis sets are there:
   call Count_lines_in_file(FN, N)
   allocate(BASIS_SET_NAMES(N))
   
   ! Read from this file:
   do i = 1, N
      read(FN,*) BASIS_SET_NAMES(i)
   enddo
   ! Close and delete temporary file:
   if (file_exist) close(FN, status='delete')

   ! Check if the requested basis set is present:
   basis_exists = .false.
   if (present(Basis_name)) then
      do i = 1, N
         if (trim(adjustl(Basis_name)) ==  trim(adjustl(BASIS_SET_NAMES(i)))) then
            print*, 'The requested basis set file ', trim(adjustl(Basis_name)), ' exists in our database'
            basis_exists = .true.
            N_basis_set = i ! the number of basis set
         endif
      enddo
       if (.not.basis_exists) then
          print*, 'The given basis set name ', trim(adjustl(Basis_name)), ' is not found in our database:'
          print*, 'Directory ', trim(adjustl(folder))
          print*, 'Using default basis set instead...'
       endif
   endif

   ! Find basis set and read parameters:
   if (.not.basis_exists) then ! we don't have this basis set
      print*, 'Searching for default basis set in our database...'
      call read_basis_for_element(trim(adjustl(folder))//path_sep//BASIS_SET_NAMES, El_Name, Basis, basis_found)
   else
      print*, 'Searching for selected element in our database...'
      call read_basis_for_element(trim(adjustl(folder))//path_sep//BASIS_SET_NAMES, El_Name, Basis, basis_found, N_basis_set)
   endif

   if (.not.basis_found) then ! We don't have basis set for the chosen element:
      write(Error_message, '(a,a)') 'No basis set found in the database for the element ', El_Name
      print*, trim(adjustl(Error_message))
   else ! we have basis set; now normalize it properly
      call normalize_basis_set(Basis) ! see below
   endif

   print*, '--------------------------------------------'

!     do i = 1, size(Basis)
!        print*, Basis(i)%MO_type	! type of orbital: S, P, D ...
!        print*, 'a', Basis(i)%GTO(:)%alpha	! contracted Gaussian Type Orbital, defined above
!        print*, 'C', Basis(i)%GTO(:)%C
!        print*, 'N', Basis(i)%GTO(:)%Norm
!        print*, Basis(i)%AM(:)	! Cartesian angular momentum quantum number in X,Y,Z directions
!        print*, Basis(i)%zeta	! zeta for STO adjustment: [Szabo, Ostlund, "Modern Quantum Chemistry", p. 185]
!        print*, Basis(i)%index
!     enddo
end subroutine find_basis_set_for_element

 
subroutine read_basis_for_element(BASIS_SET_NAMES, El_Name, Basis, basis_found, N_basis_set)
   character(100), dimension(:), intent(in) :: BASIS_SET_NAMES	! names of all basis sets
   character(*), intent(in) :: El_Name 		! chemical element name
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(in), optional :: N_basis_set	! basis set number, if known
   logical, intent(out) :: basis_found
   !------------------------------------------
   integer :: N, N_all, i

   N_all = size(BASIS_SET_NAMES) ! that's how many basis sets we have in all
   
   if (present(N_basis_set)) then
      N = N_basis_set ! start from this one
   else
      N = 1 ! start from the start
   endif
   
   basis_found = .false.
   
   call find_basis_parameters(BASIS_SET_NAMES, N, El_Name, Basis, basis_found)
   
   if (.not.basis_found) then ! search for alterntive
      print*, 'Trying other basis sets:'
      do i = 1, N_all
         call find_basis_parameters(BASIS_SET_NAMES, i, El_Name, Basis, basis_found)
         if (basis_found) then
            print*, 'Basis set for element ', El_Name, ' from the file ', trim(adjustl(BASIS_SET_NAMES(i)))
            print*, 'will be used for further calculations'
            exit
         endif
      enddo ! i
   endif

end subroutine read_basis_for_element


subroutine find_basis_parameters(BASIS_SET_NAMES, N, El_Name, Basis, basis_found)
   character(100), dimension(:), intent(in) :: BASIS_SET_NAMES	! names of all basis sets
   integer, intent(in) :: N	! number of the basis set file
   character(*), intent(in) :: El_Name 		! chemical element name
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   logical, intent(inout) :: basis_found ! did we find the basis set for this element in this database?
   !--------------------------------------------
   character(100) :: mark_line, read_line
   character(3) :: read_element
   integer :: FN, Reason, i, N_lines
   mark_line = '****' ! that is separator used in the basis sets files
   
   FN = 9997
   open(unit=FN, file=trim(adjustl(BASIS_SET_NAMES(N))))
   call Count_lines_in_file(FN, N_lines) ! how many line in total
     
   i = 0
   CH_BASIS:do
      i = i + 1
      read(FN,*,IOSTAT=Reason) read_line
      if (Reason .LT. 0) then ! ... end of file reached ...
         basis_found = .false. ! no basis set in this file for this element
      else
         if ( trim(adjustl(read_line)) == trim(adjustl(mark_line)) ) then ! here new element starts:
            read(FN,*,IOSTAT=Reason) read_element
            i = i + 1
            if (trim(adjustl(read_element)) == El_Name) then
               print*, 'Basis set for element ', El_Name, ' is found in the file ', trim(adjustl(BASIS_SET_NAMES(N)))
               call read_basis_parameters(FN, Basis, i, basis_found)
               exit CH_BASIS
            endif
         endif
      endif
!       print*, i, trim(adjustl(read_line))
      if (i >= N_lines) then
         print*, 'Basis set for element ', El_Name, ' is not found in the file ', trim(adjustl(BASIS_SET_NAMES(N)))
         exit
      endif
   enddo CH_BASIS
   close(FN)
end subroutine find_basis_parameters


subroutine read_basis_parameters(FN, Basis, i_line, basis_found)
   integer, intent(in) :: FN	! file number to read the parameters from (must be already opened and read up to the line with parametrs)
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(inout) :: i_line	! number of line that we are reading from the file
   logical, intent(inout) :: basis_found	! did we find acceptable basis set for the element given?
   !------------------------------------
   character(100) :: read_line
   character(2) :: MO_type
   integer :: Reason, i, num_func, i_basis, MO_index
   real(8) :: temp_z
   !------------------------------------
   i_line = i_line - 1 ! to start counting from the correct number, since we add 1 at the beginning of the cycle below
   i_basis = 0
   MO_index = 0
   ! Read all the parameters from the file:
   READ_CICLE:do
      i_line = i_line + 1 ! we are reading the next line
      read(FN,'(a)',IOSTAT=Reason) read_line    
!       print*, 'LINE:', trim(adjustl(read_line))
            
      if (Reason < 0) then !... end of file reached ...
         exit READ_CICLE
      else if (Reason > 0) then !... something went wrong
         print*, 'Problem durin reading BASIS SET file #', FN, ' in line #', i_line
         exit READ_CICLE
      else ! normal reading line:
         ! If this line contains a symbol, this is a starting line for a new molecular orbital (MO) parameters:
!          print*, 'TEST 0: ', read_line
         if ( IT_IS_A_LETTER(read_line(1:1)) ) then
            read(read_line,*) MO_type, num_func, temp_z
            MO_index = MO_index + 1 ! next set of functions for the next MO starts in the file:
!             print*, 'TYPE: ',MO_type, num_func, temp_z
            ! Read from the file all the parameters for this orbital:
            call get_basis_set_parameters(Basis, i_basis, num_func=num_func, MO_type=MO_type, zeta=temp_z, FN=FN, i_line=i_line, MO_index=MO_index) ! see below            
         else if (trim(adjustl(read_line(1:4))) == '****') then ! this line contains the end of element basis set definition
            exit READ_CICLE
         else ! this line contains parameters of the basis set function
            ! read the function parameters in the following lines:
            print*, 'Unknown format of the line encounterred, #', i_line
         endif
      endif ! (Reason < 0)      
   enddo READ_CICLE
   basis_found = .true. ! we found the basis set we can use
   print*, 'Basis set parameters are read from the file.'
end subroutine read_basis_parameters


subroutine get_basis_set_parameters(Basis, i_basis, i_primitive_basis, num_func, MO_type, alpha, zeta, C, AM, FN, i_line, MO_index)
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(inout) :: i_basis ! which basis function to write parameters into
   integer, intent(inout), optional :: i_primitive_basis ! which primitive basis function to write parameters into
   integer, intent(inout), optional :: num_func	! how many primitive gaussian basis functions are to be used
   character(*), intent(in), optional :: MO_type	! type of orbital: S, P, D ...
   real(8), intent(in), optional :: alpha	! alpha {1/r^2} => [1/au^2]
   real(8), intent(in), optional :: zeta	! zeta
   real(8), intent(in), optional :: C		! coefficients before exponents
   real(8), dimension(3), intent(in), optional :: AM	! Cartesian angular momentum quantum number in X,Y,Z directions
   integer, intent(in), optional :: FN	! file number to read basis set parameters from (must be already opened!)
   integer, intent(inout), optional :: i_line	! from which line we continue reading the file
   integer, intent(in), optional :: MO_index	! index of the MO
   !-------------------------
   integer :: i, Size_old, Reason, leng, j
   real(8), dimension(:,:), allocatable :: Coef
   
   ! Set angular momenta (and other parameters) for orbitals:
   if (present(MO_type)) then
      leng = LEN(trim(adjustl(MO_type))) ! how many functions are given in the file depends on the mark "S", "SP" and so on
      
      ! Read GTO primitives parameteres from the file, if desired by the user:
      if (present(FN) .and. present(i_line) .and. present(num_func)) then
         allocate( Coef(num_func,1+leng) )
         do i = 1, num_func
            i_line = i_line + 1
            read(FN,*,IOSTAT=Reason) Coef(i,:) ! read alphas and Cs from a file with basis set parameters (3-21G.gbs and some such)
         enddo
      endif
      
      ! Construct angular momenta:
      TY:do i = 1, leng ! for all parameters
         Size_old = size(Basis) ! how many functions there were before adding new ones
         
         select case (trim(adjustl(MO_type(i:i)))) ! orbital type
         case('s', 'S') ! this adds only one function to the basis set:
            ! Extend the array of basis set functions by 1 (s):
            call extend_basis(Basis, increase_by=1, set_GTO_size=num_func) ! see below
            ! The subroutine set_AM is given below; it sets cartesian angular momenta into the basis set:
            call set_AM(i_basis, Basis, 0, 0, 0) ! s
         case('p', 'P') ! this adds 3 functions to the basis set:
            ! Extend the array of basis set functions by 3 (px, py, pz):
            call extend_basis(Basis, increase_by=3, set_GTO_size=num_func) ! see below
            ! Set cartesian angular momenta into the basis set for these new basis functions:
            call set_AM(i_basis, Basis, 1, 0, 0) ! px
            call set_AM(i_basis, Basis, 0, 1, 0) ! py
            call set_AM(i_basis, Basis, 0, 0, 1) ! pz
         case('d', 'D') ! this adds 5 functions to the basis set:
            ! Extend the array of basis set functions by 6:
            call extend_basis(Basis, increase_by=6, set_GTO_size=num_func) ! see below
            ! Set cartesian angular momenta into the basis set for these new basis functions:
            call set_AM(i_basis, Basis, 2, 0, 0) ! dx2
            call set_AM(i_basis, Basis, 0, 2, 0) ! dy2
            call set_AM(i_basis, Basis, 0, 0, 2) ! dz2
            call set_AM(i_basis, Basis, 1, 1, 0) ! dxy
            call set_AM(i_basis, Basis, 0, 1, 1) ! dyz
            call set_AM(i_basis, Basis, 1, 0, 1) ! dxz
         case('f', 'F') ! this adds 10 functions to the basis set:
            ! Extend the array of basis set functions by 10:
            call extend_basis(Basis, increase_by=10, set_GTO_size=num_func) ! see below
            ! Set cartesian angular momenta into the basis set for these new basis functions:
            call set_AM(i_basis, Basis, 3, 0, 0) ! fx3
            call set_AM(i_basis, Basis, 0, 3, 0) ! fy3
            call set_AM(i_basis, Basis, 0, 0, 3) ! fz3
            call set_AM(i_basis, Basis, 2, 1, 0) ! fx2y1
            call set_AM(i_basis, Basis, 2, 0, 1) ! fx2z1
            call set_AM(i_basis, Basis, 1, 2, 0) ! fx1y2
            call set_AM(i_basis, Basis, 1, 0, 2) ! fx1z2
            call set_AM(i_basis, Basis, 0, 2, 1) ! fy2z1
            call set_AM(i_basis, Basis, 0, 1, 2) ! fz2y1
            call set_AM(i_basis, Basis, 1, 1, 1) ! fx1y1z1
         case('g', 'G') ! this adds 15 functions to the basis set:
            ! Extend the array of basis set functions by 15:
            call extend_basis(Basis, increase_by=15, set_GTO_size=num_func) ! see below
            ! Set cartesian angular momenta into the basis set for these new basis functions:
            call set_AM(i_basis, Basis, 4, 0, 0) ! gx4
            call set_AM(i_basis, Basis, 0, 4, 0) ! gy4
            call set_AM(i_basis, Basis, 0, 0, 4) ! gz4
            call set_AM(i_basis, Basis, 3, 1, 0) ! gx3y1
            call set_AM(i_basis, Basis, 3, 0, 1) ! gx3z1
            call set_AM(i_basis, Basis, 1, 3, 0) ! gy3x1
            call set_AM(i_basis, Basis, 0, 3, 1) ! gy3z1
            call set_AM(i_basis, Basis, 1, 0, 3) ! gz3x1
            call set_AM(i_basis, Basis, 0, 1, 3) ! gz3y1
            call set_AM(i_basis, Basis, 2, 2, 0) ! gx2y2
            call set_AM(i_basis, Basis, 2, 0, 2) ! gx2z2
            call set_AM(i_basis, Basis, 0, 2, 2) ! gy2z2
            call set_AM(i_basis, Basis, 2, 1, 1) ! gx2y1z1
            call set_AM(i_basis, Basis, 1, 2, 1) ! gx1y2z1
            call set_AM(i_basis, Basis, 1, 1, 2) ! gx1y1z2
         case default
            print*, 'Unknown orbital type encounterred'
         endselect

         ! Save index of MO (since many atoms will use the same basis set parameters, mark it with an index to find it later):
         if (present(MO_index)) then
            Basis(Size_old+1:size(Basis))%index = MO_index
         endif
         
         do j = Size_old+1, size(Basis) ! for all new basis functions:
            ! Save basis function type (S, P, D, F, G):
            Basis(j)%MO_type = MO_type(i:i)

            ! Save the coefficients:
            Basis(j)%GTO(:)%alpha = Coef(:,1)
            Basis(j)%GTO(:)%C = Coef(:,i+1)
            
            ! STO zeta to be adjusted to:
            if (present(zeta)) then
               Basis(j)%zeta = zeta ! save what was read from the file
               Basis(j)%GTO(:)%alpha = Basis(j)%GTO(:)%alpha*(zeta*zeta) ! renormalize accordingly (Szabo, "Quantum Modern Chemistry", p.185)
            endif
         enddo ! i = Size_old, size(Basis)
         
      enddo TY
      
      if (allocated(Coef)) deallocate(Coef) ! after the use is finished, free the memory
   else if (present(zeta)) then ! Only STO zeta to be adjusted to:
      Basis(i_basis)%zeta = zeta
      Basis(i_basis)%GTO(:)%alpha = Basis(i_basis)%GTO(:)%alpha*(zeta*zeta)
   endif ! (present(MO_type))
   
   if (present(num_func)) then ! that's how many primitive GTO functions are to be used:
      if (.not.allocated(Basis(i_basis)%GTO)) allocate(Basis(i_basis)%GTO(num_func))
   endif
   
   if (present(alpha).and. present(i_primitive_basis)) then ! alpha parameters for the given primitive GTO:
      Basis(i_basis)%GTO(i_primitive_basis)%alpha = alpha
   endif
   
   if (present(C) .and. present(i_primitive_basis)) then ! coefficients for the given primitive GTO:
      Basis(i_basis)%GTO(i_primitive_basis)%C = C
   endif
   
   if (present(AM)) then ! in case we want to overwrite angular momenta:
      Basis(i_basis)%AM = AM
   endif
end subroutine get_basis_set_parameters


pure subroutine set_AM(i_basis, Basis, k, l, m, MO_type)
   integer, intent(inout) :: i_basis ! index of the basis function
   type(Basis_set), dimension(:), intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(in) :: k, l, m ! cartesian angular momenta
   character(1), intent(in), optional :: MO_type    ! type of orbital
   i_basis = i_basis + 1 ! next basis function
   Basis(i_basis)%AM(1) = k
   Basis(i_basis)%AM(2) = l
   Basis(i_basis)%AM(3) = m
   if (present(MO_type)) then ! save type of orbital
      Basis(i_basis)%MO_type = MO_type
   endif
end subroutine set_AM





subroutine set_xTB_AO(MO_type, Basis)
   character(1), intent(in) :: MO_type  ! orbital type (for xTB only S, P, D are used; F and G are for completeness)
   type(Basis_set_STO), intent(inout) :: Basis ! basis set functions coefficients
   !-----------------------
   integer :: i_basis  ! index of the basis set; number of GTO primitives; principal quantum number

   i_basis = 0  ! to start with

   select case (trim(adjustl(MO_type))) ! orbital type
   case('s', 'S') ! this adds only one function to the basis set:
      ! The subroutine set_AM is given below; it sets cartesian angular momenta into the basis set:
      allocate(Basis%AM_set(1)) ! s
      call set_AM_single(i_basis, Basis, 0, 0, 0) ! s
   case('p', 'P') ! this adds 3 functions to the basis set:
      ! Set cartesian angular momenta into the basis set for these new basis functions:
      allocate(Basis%AM_set(3)) ! p
      call set_AM_single(i_basis, Basis, 1, 0, 0) ! px
      call set_AM_single(i_basis, Basis, 0, 1, 0) ! py
      call set_AM_single(i_basis, Basis, 0, 0, 1) ! pz
   case('d', 'D') ! this adds 5 functions to the basis set:
      ! Set cartesian angular momenta into the basis set for these new basis functions:
      allocate(Basis%AM_set(6)) ! d
      call set_AM_single(i_basis, Basis, 2, 0, 0) ! dx2
      call set_AM_single(i_basis, Basis, 0, 2, 0) ! dy2
      call set_AM_single(i_basis, Basis, 0, 0, 2) ! dz2
      call set_AM_single(i_basis, Basis, 1, 1, 0) ! dxy
      call set_AM_single(i_basis, Basis, 0, 1, 1) ! dyz
      call set_AM_single(i_basis, Basis, 1, 0, 1) ! dxz
   case('f', 'F') ! this adds 10 functions to the basis set:
      ! Set cartesian angular momenta into the basis set for these new basis functions:
      allocate(Basis%AM_set(10)) ! f
      call set_AM_single(i_basis, Basis, 3, 0, 0) ! fx3
      call set_AM_single(i_basis, Basis, 0, 3, 0) ! fy3
      call set_AM_single(i_basis, Basis, 0, 0, 3) ! fz3
      call set_AM_single(i_basis, Basis, 2, 1, 0) ! fx2y1
      call set_AM_single(i_basis, Basis, 2, 0, 1) ! fx2z1
      call set_AM_single(i_basis, Basis, 1, 2, 0) ! fx1y2
      call set_AM_single(i_basis, Basis, 1, 0, 2) ! fx1z2
      call set_AM_single(i_basis, Basis, 0, 2, 1) ! fy2z1
      call set_AM_single(i_basis, Basis, 0, 1, 2) ! fz2y1
      call set_AM_single(i_basis, Basis, 1, 1, 1) ! fx1y1z1
   case('g', 'G') ! this adds 15 functions to the basis set:
      ! Set cartesian angular momenta into the basis set for these new basis functions:
      allocate(Basis%AM_set(15)) ! g
      call set_AM_single(i_basis, Basis, 4, 0, 0) ! gx4
      call set_AM_single(i_basis, Basis, 0, 4, 0) ! gy4
      call set_AM_single(i_basis, Basis, 0, 0, 4) ! gz4
      call set_AM_single(i_basis, Basis, 3, 1, 0) ! gx3y1
      call set_AM_single(i_basis, Basis, 3, 0, 1) ! gx3z1
      call set_AM_single(i_basis, Basis, 1, 3, 0) ! gy3x1
      call set_AM_single(i_basis, Basis, 0, 3, 1) ! gy3z1
      call set_AM_single(i_basis, Basis, 1, 0, 3) ! gz3x1
      call set_AM_single(i_basis, Basis, 0, 1, 3) ! gz3y1
      call set_AM_single(i_basis, Basis, 2, 2, 0) ! gx2y2
      call set_AM_single(i_basis, Basis, 2, 0, 2) ! gx2z2
      call set_AM_single(i_basis, Basis, 0, 2, 2) ! gy2z2
      call set_AM_single(i_basis, Basis, 2, 1, 1) ! gx2y1z1
      call set_AM_single(i_basis, Basis, 1, 2, 1) ! gx1y2z1
      call set_AM_single(i_basis, Basis, 1, 1, 2) ! gx1y1z2
   case default
      print*, 'Unknown orbital type encounterred'
   endselect
end subroutine set_xTB_AO



pure subroutine set_AM_single(i_basis, Basis, k, l, m)
   integer, intent(inout) :: i_basis ! index of the basis function
   type(Basis_set_STO), intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(in) :: k, l, m ! cartesian angular momenta
   i_basis = i_basis + 1 ! next basis function
   Basis%AM_set(i_basis)%AM(1) = k
   Basis%AM_set(i_basis)%AM(2) = l
   Basis%AM_set(i_basis)%AM(3) = m
end subroutine set_AM_single





pure subroutine normalize_basis_set(Basis)
! The used GTOs are normalized as follows:
! f(x,y,z) = SUM( Norm[i]*C[i]*exp(-a[i]*r^2) ; i=1..size(Basis%GTO(:)) )
! int( f(x,y,z)^2, dx dy dz = -infinity..infinity ) = 1.0e0
   type(Basis_set), dimension(:), intent(inout), target :: Basis ! basis set functions coefficients
   real(8), pointer :: alpha
   integer, pointer :: l, m, n
   real(8) :: Norm, factorials, pows, temp
   integer :: i, j, N_basis, N_GTO, L_sum
   N_basis = size(Basis)
   pows = 0.75d0     ! parameter to use later
   do i = 1, N_basis ! for all basis functions
      ! Cartesian angular momenta:
      l => Basis(i)%AM(1)
      m => Basis(i)%AM(2)
      n => Basis(i)%AM(3)
      L_sum = l+m+n

      factorials = GTO_normalization_factorials(l, m, n)    ! below

      N_GTO = size(Basis(i)%GTO) ! number of contracted gaussian primitives
      do j = 1, N_GTO
         alpha => Basis(i)%GTO(j)%alpha
         Basis(i)%GTO(j)%Norm = GTO_normalization_coefficient(factorials, alpha, dble(L_sum), pows) ! below
      enddo
   enddo
   nullify(l, m, n, alpha)
end subroutine normalize_basis_set


pure function GTO_normalization_coefficient(factorials, alpha, L_sum, pows) result(Norm)
   real(8) :: Norm	! output : the normalization coefficient
   real(8), intent(in) :: factorials, alpha, L_sum
   real(8), intent(in), optional :: pows
   !--------------------
   real(8) :: pow, temp
   if (present(pows)) then
      pow = pows
   else
      pow = 0.75d0      ! 3.0d0/4.0d0
   endif
   temp = (4.0d0*alpha)**(L_sum/2.0d0)
   Norm = ((2.0d0*alpha/g_Pi)**pow)*temp/DSQRT(factorials)
end function GTO_normalization_coefficient


pure function GTO_normalization_factorials(l, m, n) result(fac)
   integer, intent(in) :: l, m, n	! angular momenta
   real(8) :: fac	! factorials entering normalization coefficient for GTO
   !-----------------
   real(8) :: A, B
   A = double_factorial(2*l-1)		! module "Algebra_tools"
   B = A*double_factorial(2*m-1)		! module "Algebra_tools"
   fac = B*double_factorial(2*n-1)	! module "Algebra_tools"
end function GTO_normalization_factorials


subroutine extend_basis(Basis, increase_by, set_GTO_size)
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(in), optional :: increase_by ! by how many functions should the basis set be increased
   integer, intent(in), optional :: set_GTO_size ! number of GTO contruction functions (contraction length)
   !-----------------------
   type(Basis_set), dimension(:), allocatable :: Temp_Basis ! basis set functions coefficients for temporary storage
   integer :: N_size, N_new_size
   ! In case the basis set was not used before, allocate it:
   if (.not.allocated(Basis)) then ! the basis set already exists
      if (present(increase_by)) then ! this is the size of the basis set:
         allocate(Basis(increase_by))
      else ! use default value:
         allocate(Basis(0))
      endif
      ! Set default values for the basis set (mostly zeros)
      call set_default_values_basis(Basis) ! see below
      N_size = 0
   else
      ! Define the size of the basis set parameters:
      N_size = size(Basis)
   endif
   
   ! And the new size:
   if (present(increase_by)) then ! this is the size of the basis set:
      N_new_size = N_size + increase_by
   else ! use default value of 1 function:
      N_new_size = N_size + 1
   endif
   
   ! Allocate temporary storage for existing parameters:
   allocate(Temp_Basis(N_size))
   
   ! Temporary store already existing basis set parameters:
   Temp_Basis = Basis ! all parameters together
   
   ! Extend the basis set:
   ! 1) deallocate it:
   deallocate(Basis)
   ! 2) reallocate with the new size:
   allocate(Basis(N_new_size))
   ! 3) set all the parameters to default values:
   if (present(set_GTO_size)) then
      call set_default_values_basis(Basis, set_GTO_size=set_GTO_size) ! see below 
   else
      call set_default_values_basis(Basis) ! see below 
   endif
   ! 4) restore preexisting parameters into the old part of the basis set:
   Basis(1:N_size) = Temp_Basis(1:N_size)
end subroutine extend_basis


subroutine set_default_values_basis(Basis, set_size, set_GTO_size)
   type(Basis_set), dimension(:), allocatable, intent(inout) :: Basis ! basis set functions coefficients
   integer, intent(in), optional :: set_size ! the basis set size to be defined
   integer, intent(in), optional :: set_GTO_size ! number of GTO contruction functions (contraction length)
   !-----------------------------
   integer i, N
   ! In a surprizing case the basis was not yet allocated, do it:
   if (.not.allocated(Basis)) then ! the basis set already exists
      if (present(set_size)) then ! this is the size of the basis set:
         allocate(Basis(set_size))
      else ! use default value of 1 function:
         allocate(Basis(1))
      endif
   endif
   N = size(Basis)
   
   ! For all elements, set default values:
   do i = 1, N
      if (.not.allocated(Basis(i)%GTO) .and. present(set_GTO_size)) then
         allocate(Basis(i)%GTO(set_GTO_size))
      else
         print*, 'Subroutine set_default_values_basis is invoked without given set_GTO_size'
      endif
      Basis(i)%GTO(:)%alpha = 0.0d0
      Basis(i)%GTO(:)%C = 0.0d0
      Basis(i)%GTO(:)%Norm = 1.0d0
      Basis(i)%zeta = 0.0d0
      Basis(i)%AM(:) = 0
      Basis(i)%MO_type = ''
      Basis(i)%index = 0
   enddo
end subroutine set_default_values_basis



pure function IT_IS_A_LETTER(z) ! cehcks whether the chatacter belongs to letters or not
   logical :: IT_IS_A_LETTER
   character(*), intent(in) :: z
   if ( (z>='A' .and. z<='Z') .or. (z>='a' .and. z<='z') ) then
      IT_IS_A_LETTER = .true.
   else
      IT_IS_A_LETTER = .false.
   endif
end function IT_IS_A_LETTER



!-----------------------------------------------------
subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=604) 
       enddo
       604 continue
    endif
    i = 0
    do
        read(File_num,*, end=603)
        i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines_in_file


subroutine print_elapsed_time(start, text)
   real(8), intent(inout) :: start ! starting time
   character(*), intent(in), optional :: text	! text to print out
   character(200) :: string
   real(8) :: finish

#ifdef OMP_inside
   finish = omp_get_wtime ( ) ! OMP provided subroutine
   if (present(text)) then
      string = trim(adjustl(text))
   else
      string = 'Elapsed time: '
   endif
   write(*,'(a,es12.2,a)')  trim(adjustl(string))//'	', finish - start, ' [sec]'
   start = finish
#endif
end subroutine print_elapsed_time


 
end module BS_Basis_sets
