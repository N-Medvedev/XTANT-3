! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2022-2023 Nikita Medvedev
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
PRIVATE


interface find_3bdy_ind
   module procedure find_3bdy_ind_ch	! find by element names
   module procedure find_3bdy_ind_int	! find by kind of atoms
end interface find_3bdy_ind


! Modular parameters:
character(14) :: m_3TB_directory
character(20) :: m_3TB_onsite_data
character(3) :: m_3TB_elements_data
character(6) :: m_3TB_binary_data

parameter (m_3TB_directory = '3TB_PARAMETERS')
parameter (m_3TB_onsite_data = 'Onsite_energies.dat')
parameter (m_3TB_elements_data = 'els')
parameter (m_3TB_binary_data = 'binary')


public :: find_3bdy_ind,  m_3TB_directory, m_3TB_onsite_data, read_3TB_onsite_file, construct_3TB_filenames, &
                            read_3TB_2bdy_file, read_3TB_3bdy_file


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



subroutine read_3TB_3bdy_file(FN, Element1, Element2, TB_Hamil, error_message)
   integer, intent(in) :: FN    ! 2-bosy parameters of 3TB file to read from (must be already opened)
   character(*), intent(in) :: Element1, Element2   ! elements from the periodic table
   type(TB_H_3TB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   character(*), intent(inout) :: error_message
   !-------------------------------------------
   integer :: Reason, i, count_lines, Nsiz, sizeH, sizeS
   logical :: read_well
   character(30) :: found_tag
   character(100) :: temp_ch
   character(5000) :: string_ind
   real(8), dimension(:), allocatable :: Hdata   ! text with data containing 3-body parameters for H

   ! To start with:
   count_lines = 0
   error_message = ''

   ! Find how many lines are in the file:
   call Count_lines_in_file(FN, Nsiz)  ! module "Dealing_with_files"

   ! Read file with 3-body parameters
   ! It is done in a few steps:
   ! 1) Find the dimensions of the data strings
   RD3bd: do i = 1, Nsiz
      read(FN,*, IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) exit RD3bd

      ! Find the tag in this string, if any:
      call extract_tag(temp_ch, found_tag)  ! below

      ! Find the sizes of the arrays with the H and S parameters:
      call react_to_tag(FN, found_tag, count_lines, read_well, error_message, sizeH=sizeH)   ! below
      if (.not.read_well) exit RD3bd
   enddo RD3bd
   if (.not.read_well) goto 5002

   ! Rewind the file to start reading it again:
   count_lines = 0
   REWIND(FN)

   ! 2) Read the data from the strings into the arrays:
   ! Knowing the size of array, allocate it:
   allocate(Hdata(sizeH), source=0.0d0)

   RD2bd2: do i = 1, Nsiz
      read(FN,*, IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) exit RD2bd2

      ! Find the tag in this string, if any:
      call extract_tag(temp_ch, found_tag)  ! below

      ! Find the sizes of the arrays with the H and S parameters:
      call react_to_tag(FN, found_tag, count_lines, read_well, error_message, Hdata=Hdata, inds=string_ind)   ! below
      if (.not.read_well) exit RD2bd2
   enddo RD2bd2
   if (.not.read_well) goto 5002

   ! Assign the read parameters to the variables to be used in the code:
   call sort_3bdy_parameters(TB_Hamil, Element1, Element2, Hdata, string_ind)   ! below

   ! Convert Hamiltonian parameters from Ry to eV:
   TB_Hamil%V3bdy = TB_Hamil%V3bdy * g_Ry    ! [Ry] -> [eV]
   TB_Hamil%Hh3bdy = TB_Hamil%Hh3bdy * g_Ry  ! [Ry] -> [eV]

5002 continue
end subroutine read_3TB_3bdy_file



subroutine sort_3bdy_parameters(TB_Hamil, Element1, Element2, Hdata, string_ind)
   type(TB_H_3TB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   character(*), intent(in) :: Element1, Element2   ! elements from the periodic table
   real(8), dimension(:), intent(in) :: Hdata   ! text with data containing H
   character(*), intent(in) :: string_ind   ! string with all the indices to intepret these data
   !-------------------------
   integer :: i, block_start, block_end, current_block, str_len, H_ind
   integer :: block2_start, block2_end
   integer :: colon_pos, ind3(3), ind4(4), ind5(5), ind6(6), ind8(8), ind10(10), ind12(12), ind15(15), ind18(18)
   integer :: sh1, sh2
   character(1) :: ind_var
   character(3) :: ch_ind_atom1, ch_ind_atom2, ch_ind_atom3, ch_ind_sh1, ch_ind_sh2

   ! To start counting brakets in the data
   current_block = 0
   str_len = LEN(trim(string_ind))  ! last symbol in the string of data

   ! Go through all the indices data and sort all of them out:
   strind:do

      ! 1) Find opening and closing tags for a data block: "[" and "]":
      block_start = INDEX(string_ind(current_block+1:str_len), '[')  ! intrinsic
      block_end = INDEX(string_ind(current_block+1:str_len), ']')    ! intrinsic

      ! 2) Read indices:
      ! first, find what kind of variables are saved in this block: H or O
      colon_pos = INDEX(string_ind(current_block+block_start:current_block+block_end), ':', BACK=.true.) ! intrinsic
      ind_var = string_ind(current_block+block_start+colon_pos:current_block+block_start+colon_pos)   ! this symbol is the lable

      ! second, find the positions for the second block that contains the indices:
      block2_start = INDEX(string_ind(current_block+block_end+1:str_len), '[')  ! intrinsic
      block2_end = INDEX(string_ind(current_block+block_end+1:str_len), ']')    ! intrinsic

      ! Having the lables and the indices, sort them out:
      select case (ind_var)   ! choose the variable: H or O
      case('H')   ! there are 5 more indices: atom 1, shell 1, atom 2, shell 2, atom 3
         ! Read these indices:
         read(string_ind(current_block+block_start+1:current_block+block_end),*) &
                  ch_ind_atom1, ch_ind_sh1, ch_ind_atom2, ch_ind_sh2, ch_ind_atom3

         ! Find if the indices belong to this Hamiltonian parameter set:
         AT1: if (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom1(2:))) ) then

            ! find to which index this combination of 3 atoms corresponds to:
            H_ind = find_3bdy_ind(trim(adjustl(ch_ind_atom1(2:))), trim(adjustl(ch_ind_atom2(2:))), trim(adjustl(ch_ind_atom3(2:)))) ! below

            ! find indices of the shells:
            sh1 = ind_of_shell( trim(adjustl(ch_ind_sh1(2:2))) )   ! below
            sh2 = ind_of_shell( trim(adjustl(ch_ind_sh2(2:2))) )   ! below

            ! check if the atoms are of the same or different type:
            if ( trim(adjustl(ch_ind_atom1(2:))) == trim(adjustl(ch_ind_atom2(2:))) ) then   ! the same type, 3 indices
               ! Read indices of the overlapping atoms from the second block:
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind3(:)

               ! Now, save the parameters into Hamiltonian:
               TB_Hamil%V3bdy(H_ind, sh1, sh2, 1) = Hdata( ind3(1) )
               TB_Hamil%V3bdy(H_ind, sh1, sh2, 2) = Hdata( ind3(2) ) * 0.5d0
               TB_Hamil%V3bdy(H_ind, sh1, sh2, 3) = Hdata( ind3(2) ) * 0.5d0 ! due to permutation symmetry, g2=g3
               TB_Hamil%V3bdy(H_ind, sh1, sh2, 4) = Hdata( ind3(3) )

            else  ! different types, 4 indices
               ! Read indices of the overlapping atoms from the second block:
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)

               ! Now, save the parameters into Hamiltonian:
               do i = 1, 4
                  TB_Hamil%V3bdy(H_ind, sh1, sh2, i) = Hdata( ind4(i) )
               enddo
            endif
         endif AT1

      case('O')   ! there are 3 more indices: atom 1, atom 2, atom 3
         read(string_ind(current_block+block_start+1:current_block+block_end),*) ch_ind_atom1, ch_ind_atom2, ch_ind_atom3

         ! Find if the indices belong to this Hamiltonian parameter set:
         AT3: if (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom1(2:))) ) then
            ! Read indices of the overlapping atoms from the second block:
            read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)

            ! find to which index this combination of 3 atoms corresponds to:
            H_ind = find_3bdy_ind(trim(adjustl(ch_ind_atom1(2:))), trim(adjustl(ch_ind_atom2(2:))), trim(adjustl(ch_ind_atom3(2:)))) ! below

            ! Knowing indices of the parameters and or the H, save them:
            do i = 1, 4
               TB_Hamil%Hh3bdy(H_ind,i) = Hdata( ind4(i) )
            enddo
         endif AT3

      end select

      current_block = current_block + block_end + block2_end  ! mark these 2 blocks as done, next starts after this one

      if ((block_end == 0) .or. (block2_end == 0) .or. (current_block >= str_len)) exit strind
   enddo strind
end subroutine sort_3bdy_parameters


function ind_of_shell(sh_name) result(sh)
   integer sh  ! index of the shell
   character(*), intent(in) :: sh_name ! shell name: s, p, d
   select case (trim(adjustl(sh_name)))
   case ('s', 'S')
      sh = 1
   case ('p', 'P')
      sh = 2
   case ('d', 'D')
      sh = 3
   end select
end function ind_of_shell



subroutine read_3TB_2bdy_file(FN, Element1, Element2, TB_Hamil, error_message)
   integer, intent(in) :: FN    ! 2-bosy parameters of 3TB file to read from (must be already opened)
   character(*), intent(in) :: Element1, Element2   ! elements from the periodic table
   type(TB_H_3TB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   character(*), intent(inout) :: error_message
   !-------------------------------------------
   integer :: Reason, i, count_lines, Nsiz, sizeH, sizeS
   logical :: read_well
   character(30) :: found_tag
   character(100) :: temp_ch
   character(5000) :: string_ind
   real(8), dimension(:), allocatable :: Hdata, Sdata   ! text with data containing H and S

   ! To start with:
   count_lines = 0
   error_message = ''
   ! Default values, no overlap coefficients to start with:
   TB_Hamil%Hhavg = 0.0d0
   TB_Hamil%Hhcf = 0.0d0
   TB_Hamil%Vrfx = 0.0d0
   TB_Hamil%Srfx = 0.0d0

   ! Find how many lines are in the file:
   call Count_lines_in_file(FN, Nsiz)  ! module "Dealing_with_files"

   ! Read file with 2-body parameters
   ! It is done in a few steps:
   ! 1) Find the dimensions of the data strings
   RD2bd: do i = 1, Nsiz
      read(FN,*, IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) exit RD2bd

      ! Find the tag in this string, if any:
      call extract_tag(temp_ch, found_tag)  ! below

      ! Find the sizes of the arrays with the H and S parameters:
      call react_to_tag(FN, found_tag, count_lines, read_well, error_message, sizeH=sizeH, sizeS=sizeS)   ! below
      if (.not.read_well) exit RD2bd
   enddo RD2bd
   if (.not.read_well) goto 5001

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
      if (.not.read_well) exit RD2bd2

      ! Find the tag in this string, if any:
      call extract_tag(temp_ch, found_tag)  ! below

      ! Find the sizes of the arrays with the H and S parameters:
      call react_to_tag(FN, found_tag, count_lines, read_well, error_message, Hdata=Hdata, Sdata=Sdata, inds=string_ind)   ! below
      if (.not.read_well) exit RD2bd2
   enddo RD2bd2
   if (.not.read_well) goto 5001

   ! Assign the read parameters to the variables to be used in the code:
   call sort_2bdy_parameters(TB_Hamil, Element1, Element2, Hdata, Sdata, string_ind)   ! below

   ! Convert Hamiltonian parameters from Ry to eV:
   TB_Hamil%Hhavg = TB_Hamil%Hhavg * g_Ry ! [Ry] -> [eV]
   TB_Hamil%Hhcf = TB_Hamil%Hhcf * g_Ry ! [Ry] -> [eV]
   TB_Hamil%Vrfx = TB_Hamil%Vrfx * g_Ry ! [Ry] -> [eV]

5001 continue
end subroutine read_3TB_2bdy_file



subroutine sort_2bdy_parameters(TB_Hamil, Element1, Element2, Hdata, Sdata, string_ind)
   type(TB_H_3TB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   character(*), intent(in) :: Element1, Element2   ! elements from the periodic table
   real(8), dimension(:), intent(in) :: Hdata, Sdata   ! text with data containing H and S
   character(*), intent(in) :: string_ind   ! string with all the indices to intepret these data
   !-------------------------
   integer :: i, block_start, block_end, current_block, str_len
   integer :: block2_start, block2_end
   integer :: colon_pos, ind4(4), ind5(5), ind6(6), ind8(8), ind10(10), ind12(12), ind15(15), ind18(18)
   character(1) :: ind_var
   character(3) :: ch_ind_atom1, ch_ind_atom2, ch_ind_sh1, ch_ind_sh2

   ! To start counting brakets in the data
   current_block = 0
   str_len = LEN(trim(string_ind))  ! last symbol in the string of data

   ! Go through all the indices data and sort all of them out:
   strind:do

      ! 1) Find opening and closing tags for a data block: "[" and "]":
      block_start = INDEX(string_ind(current_block+1:str_len), '[')  ! intrinsic
      block_end = INDEX(string_ind(current_block+1:str_len), ']')    ! intrinsic

      ! 2) Read indices:
      ! first, find what kind of variables are saved in this block: H, S or O
      colon_pos = INDEX(string_ind(current_block+block_start:current_block+block_end), ':', BACK=.true.) ! intrinsic
      ind_var = string_ind(current_block+block_start+colon_pos:current_block+block_start+colon_pos)   ! this symbol is the lable

      ! second, find the positions for the second block that contains the indices:
      block2_start = INDEX(string_ind(current_block+block_end+1:str_len), '[')  ! intrinsic
      block2_end = INDEX(string_ind(current_block+block_end+1:str_len), ']')    ! intrinsic

      !print*, block_start, block_end, string_ind(current_block+block_start:current_block+block_end)
      !print*, block2_start, block2_end, string_ind(current_block+block_end+block2_start:current_block+block_end+block2_end)

      ! Having the lables and the indices, sort them out:
      select case (ind_var)   ! choose the variable: H, S or O
      case('H')   ! there are 4 more indices: atom 1, shell 1, atom 2, shell 2
         ! Read these indices:
         read(string_ind(current_block+block_start+1:current_block+block_end),*) ch_ind_atom1, ch_ind_sh1, ch_ind_atom2, ch_ind_sh2

         ! Find if the indices belong to this Hamiltonian parameter set:
!          AT1: if ( (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom1(2:))) ) .and. &
!               (trim(adjustl(Element2)) == trim(adjustl(ch_ind_atom2(2:))) ) ) then         ! [A 0]
         AT1: if ( (trim(adjustl(Element2)) == trim(adjustl(ch_ind_atom1(2:))) ) .and. &
              (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom2(2:))) ) ) then           ! [A 1]

            ! Read indices of the overlapping shells from the second block:
            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind5(:)
               !print*, ind5(:)

               ! Now, save the (s s sigma) parameters into Hamiltonian:
               do i = 1, 5
                  TB_Hamil%Vrfx(1,i) = Hdata( ind5(i) ) ! (s s sigma)
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) .or. &
                 ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind5(:)
               !print*, ind5(:)

               ! Now, save the (s p sigma) parameters into Hamiltonian:
               do i = 1, 5
                  TB_Hamil%Vrfx(2,i) = Hdata( ind5(i) ) ! (s p sigma)
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) .or. &
                 ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind5(:)
               !print*, ind5(:)

               ! Now, save the (s d sigma) parameters into Hamiltonian:
               do i = 1, 5
                  TB_Hamil%Vrfx(3,i) = Hdata( ind5(i) )  ! (s d sigma)
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind10(:)
               !print*, ind10(:)
               ! Now, save the parameters into Hamiltonian:
               do i = 1, 5
                  TB_Hamil%Vrfx(4,i) = Hdata( ind10(i) )    ! (p p sigma)
                  TB_Hamil%Vrfx(5,i) = Hdata( ind10(i+5) )  ! (p p pi)
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) .or. &
                 ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind10(:)
               !print*, ind10(:)
               ! Now, save the parameters into Hamiltonian:
               do i = 1, 5
                  TB_Hamil%Vrfx(6,i) = Hdata( ind10(i) )    ! (p d sigma)
                  TB_Hamil%Vrfx(7,i) = Hdata( ind10(i+5) )  ! (p d pi)
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind15(:)
               !print*, ind15(:)
               ! Now, save the parameters into Hamiltonian:
               do i = 1, 5
                  TB_Hamil%Vrfx(8,i) = Hdata( ind15(i) )       ! (d d sigma)
                  TB_Hamil%Vrfx(9,i) = Hdata( ind15(i+5) )     ! (d d pi)
                  TB_Hamil%Vrfx(10,i) = Hdata( ind15(i+10) )   ! (d d delta)
               enddo
            endif

         endif AT1


         !print*, ind_var, ch_ind_atom1, ch_ind_sh1, ch_ind_atom2, ch_ind_sh2
      case('S')   ! there are 4 more indices: atom 1, shell 1, atom 2, shell2
         read(string_ind(current_block+block_start+1:current_block+block_end),*) ch_ind_atom1, ch_ind_sh1, ch_ind_atom2, ch_ind_sh2

         ! Find if the indices belong to this Hamiltonian parameter set:
!          AT2: if ( (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom1(2:))) ) .and. &
!               (trim(adjustl(Element2)) == trim(adjustl(ch_ind_atom2(2:))) ) ) then           ! [B 0]
         AT2: if ( (trim(adjustl(Element2)) == trim(adjustl(ch_ind_atom1(2:))) ) .and. &
              (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom2(2:))) ) ) then            ! [B 1]

            ! Read indices of the overlapping shells from the second block:
            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind6(:)
               ! Now, save the (s s sigma) parameters into Overlap:
               do i = 1, 6
                  TB_Hamil%Srfx(1,i) = Sdata( ind6(i) )
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) .or. &
                 ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind6(:)
               ! Now, save the (s p sigma) parameters into Overlap:
               do i = 1, 6
                  TB_Hamil%Srfx(2,i) = Sdata( ind6(i) )
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) .or. &
                 ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind6(:)
               ! Now, save the (s d sigma) parameters into Overlap:
               do i = 1, 6
                  TB_Hamil%Srfx(3,i) = Sdata( ind6(i) )
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind12(:)
               ! Now, save the parameters into Hamiltonian:
               do i = 1, 6
                  TB_Hamil%Srfx(4,i) = Sdata( ind12(i) )    ! (p p sigma)
                  TB_Hamil%Srfx(5,i) = Sdata( ind12(i+6) )  ! (p p pi)
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) .or. &
                 ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind12(:)
               ! Now, save the parameters into Hamiltonian:
               do i = 1, 6
                  TB_Hamil%Srfx(6,i) = Sdata( ind12(i) )    ! (p d sigma)
                  TB_Hamil%Srfx(7,i) = Sdata( ind12(i+6) )  ! (p d pi)
               enddo
            endif

            if ( ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind18(:)
               ! Now, save the parameters into Hamiltonian:
               do i = 1, 6
                  TB_Hamil%Srfx(8,i) = Sdata( ind18(i) )       ! (d d sigma)
                  TB_Hamil%Srfx(9,i) = Sdata( ind18(i+6) )     ! (d d pi)
                  TB_Hamil%Srfx(10,i) = Sdata( ind18(i+12) )   ! (d d delta)
               enddo
            endif

         endif AT2

      case('O')   ! there are 3 more indices: atom 1, shell 1, shell2
         read(string_ind(current_block+block_start+1:current_block+block_end),*) ch_ind_atom1, ch_ind_sh1, ch_ind_sh2

         ! Find if the indices belong to this Hamiltonian parameter set:
         AT3: if (trim(adjustl(Element1)) == trim(adjustl(ch_ind_atom1(2:))) ) then

            ! Read indices of the overlapping shells from the second block:
            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4 ! Diagonal elements go into average atom part
                  !TB_Hamil%Hhcf(1,1,i) = Hdata( ind4(i) )
                  TB_Hamil%Hhavg(1,i) = Hdata( ind4(i) )
!                   print*, 's-s O:', TB_Hamil%Hhavg(1,i)
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4
                  TB_Hamil%Hhcf(1,2,i) = Hdata( ind4(i) )
                  !print*, 'cf s-p', TB_Hamil%Hhcf(1,2,i)
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 's') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4
                  TB_Hamil%Hhcf(1,3,i) = Hdata( ind4(i) )
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4
                  TB_Hamil%Hhcf(2,1,i) = Hdata( ind4(i) )
                  !print*, 'cf p-s', TB_Hamil%Hhcf(2,1,i)
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4
                  TB_Hamil%Hhcf(2,3,i) = Hdata( ind4(i) )
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 's') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4
                  TB_Hamil%Hhcf(3,1,i) = Hdata( ind4(i) )
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind4(:)
               !print*, ind4(:)
               do i = 1, 4
                  TB_Hamil%Hhcf(3,2,i) = Hdata( ind4(i) )
               enddo
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'p') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'p') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind8(:)
               !print*, ind8(:)
               do i = 1, 4
!                   TB_Hamil%Hhcf(2,2,i) = Hdata( ind8(i) )
!                   TB_Hamil%Hhavg(2,i) = Hdata( ind8(i+4) )
                  TB_Hamil%Hhavg(2,i) = Hdata( ind8(i) )
                  TB_Hamil%Hhcf(2,2,i) = Hdata( ind8(i+4) ) !/ 1.5d0
!                   print*, 'p-p O:', ind8(i), TB_Hamil%Hhavg(2,i), TB_Hamil%Hhcf(2,2,i)
               enddo
!                pause '3TB reading'
            endif

            if ( (trim(adjustl(ch_ind_sh1(2:2))) == 'd') .and. (trim(adjustl(ch_ind_sh2(2:2))) == 'd') ) then
               read(string_ind(current_block+block_end+block2_start+1:current_block+block_end+block2_end-1),*) ind8(:)
               !print*, ind8(:)
               do i = 1, 4
!                   TB_Hamil%Hhcf(3,3,i) = Hdata( ind8(i) )
!                   TB_Hamil%Hhavg(3,i) = Hdata( ind8(i+4) )
                  TB_Hamil%Hhavg(3,i) = Hdata( ind8(i) )
                  TB_Hamil%Hhcf(3,3,i) = Hdata( ind8(i+4) ) !/ 1.5d0
!                   print*, 'd-d O:', TB_Hamil%Hhavg(3,i)
               enddo
            endif

         endif AT3

      end select

      current_block = current_block + block_end + block2_end  ! mark these 2 blocks as done, next starts after this one

      if ((block_end == 0) .or. (block2_end == 0) .or. (current_block >= str_len)) exit strind
   enddo strind
end subroutine sort_2bdy_parameters




subroutine react_to_tag(FN, found_tag, count_lines, read_well, error_message, sizeH, sizeS, Hdata, Sdata, inds)
   integer, intent(in) :: FN  ! file number to read from
   character(*), intent(in) :: found_tag    ! tag identifies how to interprete this line in 3TB file
   integer, intent(inout) :: count_lines  ! which line we are reading from the file
   logical, intent(inout) :: read_well    ! did it read well
   character(*), intent(inout) :: error_message ! error message if id did not read well
   integer, intent(out), optional :: sizeH, sizeS  ! sizes of data arrays containing H and S coefficients
   real(8), dimension(:), intent(out), optional :: Hdata, Sdata   ! text with data containing H and S
   character(*), intent(out), optional :: inds   ! string with all the indices to intepret these data
   !-----------------------
   integer :: N_ver, Reason, str_len
   character(30) :: temp_ch
   character(5000) :: temp_string

   ! Check if the tag is a key word that we need:
   select case(trim(adjustl(found_tag)))
   case ('inds') ! dictionary of indices, attributing array of data to proper parameters
      if (present(inds)) then ! read the list of indices
         backspace(FN) ! to read this line again, extracting the parameters we need
         count_lines = count_lines - 1
         !First read the data as text:
         read(FN, '(a)', IOSTAT=Reason) temp_string
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         ! Now, interpret the data from the string (excluding the tag at the end with "/"):
         str_len = LEN(trim(temp_string)) ! length of the string
         inds = temp_string(55:str_len-7)
      endif

   case ('datH') ! Hamiltonian 2-body parameters
      if (present(Hdata)) then ! read the H coefficients into the array
         backspace(FN) ! to read this line again, extracting the parameters we need
         count_lines = count_lines - 1
         !First read the data as text:
         read(FN, '(a)', IOSTAT=Reason) temp_string
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         ! Now, interpret the data from the string (excluding the tag at the end with "/"):
         str_len = LEN(trim(temp_string)) ! length of the string
         read(temp_string(11:str_len-7), *, IOSTAT=Reason) Hdata

         if (Hdata(size(Hdata)) == 0.0d0) then ! probably reading of the string failed
            print*, 'Hdata', Hdata
            error_message = 'Reading from 3TB H parameters failed'
            read_well = .false.
         endif
      endif

   case ('datS') ! Overlap matrix 2-body parameters
      if (present(Sdata)) then ! read the S coefficients into the array
         backspace(FN) ! to read this line again, extracting the parameters we need
         count_lines = count_lines - 1

         !First read the data as text:
         read(FN, '(a)', IOSTAT=Reason) temp_string
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         ! Now, interpret the data from the string (excluding the tag at the end with "/"):
         str_len = LEN(trim(temp_string)) ! length of the string
         read(temp_string(11:str_len-7), *, IOSTAT=Reason) Sdata

         if (Sdata(size(Sdata)) == 0.0d0) then ! probably reading of the string failed
            print*, 'Sdata', Sdata
            error_message = 'Reading from 3TB S parameters failed'
            read_well = .false.
         endif
      endif

   case ('sizeH') ! size of the array with H parameters
      if (present(sizeH)) then
         backspace(FN) ! to read this line again, extracting the parameter we need
         count_lines = count_lines - 1

         !First read the data as text:
         !read(FN,'(11X,i3)', IOSTAT=Reason) sizeH
         read(FN, '(a)', IOSTAT=Reason) temp_ch
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         ! Now, interpret the data from the string (excluding the tag at the end with "/"):
         str_len = LEN(trim(temp_ch)) ! length of the string
         read(temp_ch(12:str_len-8), *, IOSTAT=Reason) sizeH
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         if (sizeH < 1) then ! probably reading of the string failed
            print*, 'sizeH', Sdata
            error_message = 'Reading from 3TB sizeH parameter failed'
            read_well = .false.
         endif
      endif

   case ('sizeS') ! size of the array with S parameters
      if (present(sizeS)) then
         backspace(FN) ! to read this line again, extracting the parameter we need
         count_lines = count_lines - 1

         !First read the data as text:
         !read(FN,'(11X,i2)', IOSTAT=Reason) sizeS
         read(FN, '(a)', IOSTAT=Reason) temp_ch
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"

         ! Now, interpret the data from the string (excluding the tag at the end with "/"):
         str_len = LEN(trim(temp_ch)) ! length of the string
         read(temp_ch(12:str_len-8), *, IOSTAT=Reason) sizeS
      endif

   case ('version') ! index of the parameterization type
      backspace(FN) ! to read this line again, extracting the parameter we need
      count_lines = count_lines - 1
      read(FN,'(13X,i1)', IOSTAT=Reason) N_ver
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (N_ver /= 2) then
         print*, '*=============================================*'
         print*, ' Attention: 3TB parameterization version is not supported!'
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
      inquire(file=trim(adjustl(Filename_3body)),exist=file_exists)
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



pure function find_3bdy_ind_ch(El1, El2, El3) result(ind)
   integer :: ind ! index of the parameters to use in setting 3-body Hamiltonian
   character(*), intent(in) :: El1, El2, El3 ! elements of these 3 atoms

   if ((trim(adjustl(El1)) == trim(adjustl(El2))) .and. (trim(adjustl(El2)) == trim(adjustl(El3)))) then ! same elements (1 1 1) or (2 2 2) etc.
      ind = 1
   else if (trim(adjustl(El1)) == trim(adjustl(El2))) then ! (1 1 2) or (2 2 1) etc.
      ind = 1
   else if (trim(adjustl(El1)) == trim(adjustl(El3))) then ! (1 2 1) or (2 1 2) etc.
      ind = 2
   else if (trim(adjustl(El2)) == trim(adjustl(El3))) then ! (1 2 2) or (2 1 1) etc.
      ind = 3
   endif
end function find_3bdy_ind_ch

pure function find_3bdy_ind_int(El1, El2, El3) result(ind)
   integer :: ind ! index of the parameters to use in setting 3-body Hamiltonian
   integer, intent(in) :: El1, El2, El3 ! kinds of elements of these 3 atoms

   if ((El1 == El2) .and. (El2 == El3)) then ! same elements (1 1 1) or (2 2 2) etc.
      ind = 1
   else if (El1 == El2) then ! (1 1 2) or (2 2 1) etc.
      ind = 1
   else if (El1 == El3) then ! (1 2 1) or (2 1 2) etc.
      ind = 2
   else if (El2 == El3) then ! (1 2 2) or (2 1 1) etc.
      ind = 3
   endif
end function find_3bdy_ind_int



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
