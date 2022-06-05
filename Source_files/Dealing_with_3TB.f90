! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2022 Nikita Medvedev
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
! This module contains subroutines for dealing with the parameters of 3TB parameterization
! The parameters are taken from the ThreeBodyTB code: https://github.com/usnistgov/ThreeBodyTB.jl

MODULE Dealing_with_3TB
use Universal_Constants   ! let it use universal constants
use Objects   ! since it uses derived types, it must know about them from module 'Objects'
use Dealing_with_files, only : read_file, Count_lines_in_file
use Dealing_with_DFTB, only: same_or_different_atom_types

implicit none

! Modular parameters:
character(14) :: m_3TB_directory
character(20) :: m_3TB_onsite_data
character(3) :: m_3TB_elements_data
character(6) :: m_3TB_binary_data

parameter (m_3TB_directory = '3TB_PARAMETERS')
parameter (m_3TB_onsite_data = 'Onsite_energies.dat')
parameter (m_3TB_elements_data = 'els')
parameter (m_3TB_binary_data = 'binary')

 contains


subroutine read_3TB_onsite_file(FN, Element1, TB_Hamil, N_basis_siz, error_message)
   integer, intent(in) :: FN    ! onsite parameters of 3TB file to read from (must be already opened)
   character(*), intent(in) :: Element1   ! element from the periodic table
   type(TB_H_3TB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   integer, intent(out) :: N_basis_siz ! basis size for this element
   character(*), intent(inout) :: error_message
   !-------------------------------------------
   integer :: Reason, count_lines, i, Nsiz, BSsiz
   real(8) :: Es, Ep, Ed
   character(3) :: El_name
   logical :: read_well, found_element

   ! To start with:
   N_basis_siz = 0
   count_lines = 0
   error_message = ''
   found_element = .false.

   ! Find out how many elements we have in the data file:
   call Count_lines_in_file(FN, Nsiz, skip_lines=1)  ! module "Dealing_with_files"

   ! Skip the first comment line:
   read(FN,*, IOSTAT=Reason)
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not. read_well) goto 4000 ! exit

   ! Find the element we need:
   ELEM:do i = 1, Nsiz ! check all lines until find the right one
      read(FN,*, IOSTAT=Reason) El_name, BSsiz, Es, Ep, Ed
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not. read_well) then
         goto 4000 ! exit
      endif
      ! Check if this is the element we are looking for:
      if ( trim(adjustl(El_name)) == trim(adjustl(Element1)) ) then
         found_element = .true.
         ! Save the onsite data into the Hamiltonian parameters:
         N_basis_siz = BSsiz
         TB_Hamil%Es = Es
         TB_Hamil%Ep = Ep
         TB_Hamil%Ed = Ed
         exit ELEM  ! no need to continue, found it
      endif
   enddo ELEM

4000 continue
   if (.not.found_element) then ! cannot continue calculations without parameters
      error_message = 'Onsite 3TB parameters for element '//trim(adjustl(Element1))//' not found'
   endif
end subroutine read_3TB_onsite_file


subroutine read_3TB_2bdy_file(FN, Element1, Element2, TB_Hamil, error_message)
   integer, intent(in) :: FN    ! 2-bosy parameters of 3TB file to read from (must be already opened)
   character(*), intent(in) :: Element1, Element2   ! elements from the periodic table
   type(TB_H_3TB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   character(*), intent(inout) :: error_message
   !-------------------------------------------
   integer :: Reason, i, count_lines, Nsiz, sizeH, sizeS
   real(8) :: Es, Ep, Ed
   character(3) :: El_name1, El_name2
   character(30) :: found_tag
   character(100) :: temp_ch
   logical :: read_well, found_element
   real(8), dimension(:), allocatable :: Hdata, Sdata   ! text with data containing H and S
   integer, dimension(:), allocatable :: inds   ! indices to intepret these data

   ! To start with:
   count_lines = 0
   error_message = ''
   found_element = .false.

   ! Find how many lines are in the file:
   call Count_lines_in_file(FN, Nsiz)  ! module "Dealing_with_files"

   ! Read file with 2-body parameters
   ! It is done in a few steps:
   ! 1) Find the dimensions of the data strings
   RD2bd: do i = 1, Nsiz
      read(FN,*, IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

      ! Find the tag in this string, if any:
      call extract_tag(temp_ch, found_tag)  ! below

      ! Find the sizes of the arrays with the H and S parameters:
      call react_to_tag(FN, found_tag, count_lines, read_well, error_message, sizeH=sizeH, sizeS=sizeS)   ! below
   enddo RD2bd
   ! Rewind the file to start reading it again:
   count_lines = 0
   REWIND(FN)

   ! 2) Read the data from the strings into the arrays:
   ! Knowing the sizes of arrays, allocate them:
   allocate(Hdata(sizeH), source=0.0d0)
   allocate(Sdata(sizeS), source=0.0d0)
   RD2bd2: do i = 1, Nsiz
      read(FN,*, IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

      ! Find the tag in this string, if any:
      call extract_tag(temp_ch, found_tag)  ! below

      ! Find the sizes of the arrays with the H and S parameters:
      call react_to_tag(FN, found_tag, count_lines, read_well, error_message, Hdata=Hdata, Sdata=Sdata)   ! below
   enddo RD2bd2

end subroutine read_3TB_2bdy_file



subroutine react_to_tag(FN, found_tag, count_lines, read_well, error_message, sizeH, sizeS, Hdata, Sdata, inds)
   integer, intent(in) :: FN  ! file number to read from
   character(*), intent(in) :: found_tag    ! tag identifies how to interprete this line in 3TB file
   integer, intent(inout) :: count_lines  ! which line we are reading from the file
   logical, intent(inout) :: read_well    ! did it read well
   character(*), intent(inout) :: error_message ! error message if id did not read well
   integer, intent(out), optional :: sizeH, sizeS  ! sizes of data arrays containing H and S coefficients
   real(8), dimension(:), intent(out), optional :: Hdata, Sdata   ! text with data containing H and S
   integer, dimension(:), intent(out), optional :: inds   ! indices to intepret these data
   !-----------------------
   integer :: N_ver, Reason

   print*, 'Tag: ', trim(adjustl(found_tag))

   ! Check if the tag is a key word that we need:
   select case(trim(adjustl(found_tag)))
   case ('datH')
      if (present(Hdata)) then ! read the H coefficients into the array
         backspace(FN) ! to read this line again, extracting the parameter we need
         count_lines = count_lines - 1
         read(FN,'(10X,(es))', IOSTAT=Reason) Hdata
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         print*, 'Hdata', Hdata
      endif
   case ('datS')
      if (present(Sdata)) then ! read the S coefficients into the array
         backspace(FN) ! to read this line again, extracting the parameter we need
         count_lines = count_lines - 1
         read(FN,'(10X,(es))', IOSTAT=Reason) Sdata
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         print*, 'Sdata', Sdata
      endif
   case ('sizeH')
      if (present(sizeH)) then
         backspace(FN) ! to read this line again, extracting the parameter we need
         count_lines = count_lines - 1
         read(FN,'(11X,i2)', IOSTAT=Reason) sizeH
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         print*, 'sizeH:', sizeH
      endif
   case ('sizeS')
      if (present(sizeS)) then
         backspace(FN) ! to read this line again, extracting the parameter we need
         count_lines = count_lines - 1
         read(FN,'(11X,i2)', IOSTAT=Reason) sizeS
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         print*, 'sizeS:', sizeS
      endif
   case ('version')
      backspace(FN) ! to read this line again, extracting the parameter we need
      count_lines = count_lines - 1
      read(FN,'(13X,i1)', IOSTAT=Reason) N_ver
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (N_ver /= 2) then
         print*, '*=============================================*'
         print*, ' Attention: 3TB parameterization version is wrong!'
         print*, ' Version:', N_ver
         print*, ' It may cause some problems...'
         print*, ' *=============================================*'
      endif
   end select

end subroutine react_to_tag



subroutine extract_tag(string_read, found_tag)
   character(*), intent(in) :: string_read   ! string we read from xml file
   character(*), intent(out) :: found_tag    ! tag found in this string, if any
   !-----------
   integer :: pos_start, pos_end

   ! Find opening and closing tags: "<" and ">":
   pos_start = INDEX(trim(adjustl(string_read)), '<')  ! intrinsic
   pos_end = INDEX(trim(adjustl(string_read)), '>')    ! intrinsic

   if (pos_end > 0) then ! there is a tag
      found_tag = trim(adjustl(string_read(pos_start+1:pos_end-1)))
   else ! there is no tag in this string
      found_tag = ''
   endif
end subroutine extract_tag




subroutine construct_3TB_filenames(Folder_name, Element1, Element2, path_sep, Filename_2body, Filename_3body, INFO)
   character(*), intent(in) :: Folder_name   ! directory with input data
   character(*), intent(in) :: Element1, Element2   ! elements names for which we require parameterization
   character(1), intent(in) :: path_sep   ! path separator in Windows vs Linux
   character(*), intent(inout) :: Filename_2body, Filename_3body  ! corresponding filenames with 3TB parameters
   integer, intent(inout) :: INFO  ! info on files
   !---------------------------
   integer :: ToA
   character(150) :: folder
   character(13) :: name_2body, name_3body
   character(30) :: file_name2bdy, file_name3bdy
   logical :: file_exists, file_exists2

   INFO = 0 ! no error at the start

   ! Part of the file name used in all 3TB parameters files:
   name_2body = 'coef.el.2bdy.'
   name_3body = 'coef.el.3bdy.'

   ! Identify if it is the same element or not:
   ToA = same_or_different_atom_types(Element1, Element2)   ! module "Dealing_with_DFTB"

   ! Prepare the parts of the path and file name:
   if (ToA == 0) then ! same element
      ! Add the directory path:
      write(folder,'(a)') trim(adjustl(Folder_name))//trim(adjustl(m_3TB_elements_data))
      write(file_name2bdy,'(a)') trim(adjustl(name_2body))//trim(adjustl(Element1))//'.xml'
      write(file_name3bdy,'(a)') trim(adjustl(name_3body))//trim(adjustl(Element1))//'.xml'

      ! Collect the path and name of the files:
      write(Filename_2body,'(a)') trim(adjustl(folder))//path_sep//trim(adjustl(file_name2bdy))
      write(Filename_3body,'(a)') trim(adjustl(folder))//path_sep//trim(adjustl(file_name3bdy))

      ! Check if file with the requested parameters exists:
      inquire(file=trim(adjustl(Filename_2body)),exist=file_exists)
      if (.not.file_exists) then
         INFO = 1 ! file with 3-body parameters of an element not found
      endif

      ! Check if file with the requested parameters exists:
      inquire(file=trim(adjustl(Filename_3body)),exist=file_exists)
      if (.not.file_exists) then
         INFO = 2 ! file with 3-body parameters of an element not found
      endif

   else ! two different elements
      write(folder,'(a)') trim(adjustl(Folder_name))//trim(adjustl(m_3TB_binary_data))
      write(file_name2bdy,'(a)') trim(adjustl(name_2body))//trim(adjustl(Element1))//'.'//trim(adjustl(Element2))//'.xml'

      ! Collect the path and name of the files:
      write(Filename_2body,'(a)') trim(adjustl(folder))//path_sep//trim(adjustl(file_name2bdy))

      ! Check if such a file exists:
      inquire(file=trim(adjustl(Filename_2body)),exist=file_exists)
      if (.not.file_exists) then ! check if the order of elements is different
         write(file_name2bdy,'(a)') trim(adjustl(name_2body))//trim(adjustl(Element2))//'.'//trim(adjustl(Element1))//'.xml'
         write(Filename_2body,'(a)') trim(adjustl(folder))//path_sep//trim(adjustl(file_name2bdy))
         inquire(file=trim(adjustl(Filename_2body)),exist=file_exists2)
         ! Check if file with the requested parameters exists:
         if (.not.file_exists2) then
            INFO = 3 ! file with 2-body parameters of a binary compound not found
         endif
      endif

      write(file_name3bdy,'(a)') trim(adjustl(name_3body))//trim(adjustl(Element1))//'.'//trim(adjustl(Element2))//'.xml'

      ! Collect the path and name of the files:
      write(Filename_3body,'(a)') trim(adjustl(folder))//path_sep//trim(adjustl(file_name3bdy))

      ! Check if such a file exists:
      inquire(file=trim(adjustl(file_name3bdy)),exist=file_exists)
      if (.not.file_exists) then ! check if the order of elements is different
         write(file_name3bdy,'(a)') trim(adjustl(name_3body))//trim(adjustl(Element2))//'.'//trim(adjustl(Element1))//'.xml'
         write(Filename_3body,'(a)') trim(adjustl(folder))//path_sep//trim(adjustl(file_name3bdy))
         inquire(file=trim(adjustl(Filename_3body)),exist=file_exists2)
         ! Check if file with the requested parameters exists:
         if (.not.file_exists2) then
            INFO = 4 ! file with 3-body parameters of a binary compound not found
         endif
      endif
   endif

end subroutine construct_3TB_filenames



pure subroutine idnetify_basis_size_3TB(TB_Hamil, Nsiz)
   integer, intent(out) :: Nsiz  ! index: 0=s, 1=sp3, 2=sp3d5
   type(TB_H_3TB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   real(8) :: eps
   integer :: Nd, Np, Ns
   integer :: i,j 
   eps = 1.0d-6   ! [eV] precision
   Ns = 0   ! to start with
   Np = 1   ! to start with
   Nd = 1   ! to start with
   ! Check if d-band is used for any element:   
   if ( maxval(ABS(TB_Hamil(:,:)%Ed - 100.0d0)) < eps) Nd = 0
   ! Check if p band is used for any element:
   if ( maxval(ABS(TB_Hamil(:,:)%Ep - 100.0d0)) < eps) Np = 0
   ! Set basis set size:
   Nsiz = Ns + Np + Nd
end subroutine idnetify_basis_size_3TB


END MODULE Dealing_with_3TB
