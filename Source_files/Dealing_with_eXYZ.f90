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
! This module contains subroutines for dealing with the extended XYZ format

MODULE Dealing_with_eXYZ
use Dealing_with_files, only : read_file, Count_lines_in_file

implicit none
PRIVATE

!==============================================
! For substitution of atoms
type Substitute_data
   logical :: required  ! is it required?
   character(3) :: Atoms_to_substitute
   character(3) :: Atoms_substituting
   real(8) :: percentage   ! percentage of atoms to be substituted
endtype Substitute_data
!==============================================

public :: interpret_XYZ_comment_line, Substitute_data

 contains

subroutine interpret_XYZ_comment_line(line_2, Supce, ind_S, ind_R, ind_V, ind_atoms, SC_X, SC_Y, &
                                      it_is_mixture, substitution_data, Error_descript)
   character(*), intent(in) :: line_2  ! line #2 read from XYZ file
   integer, intent(out) :: ind_S, ind_R, ind_V, ind_atoms  ! indices: Element; Coordinates; Velocities; Atoms_setting
   real(8), dimension(3,3), intent(inout) :: Supce ! supercell vectors
   real(8), intent(inout) :: SC_X, SC_Y   ! sepuercell sizes along X and Y [A]
   logical, intent(inout) :: it_is_mixture  ! flag for alloy/mixture
   type(Substitute_data), intent(inout) :: substitution_data  ! flag for substitution
   character(*), intent(inout) :: Error_descript   ! error message, if any
   !----------------------
   integer :: block_start, block_end, block2_end, colon_pos, colon_pos2, current_block, str_len, eq_pos
   integer :: Reason, count_lines, i, block_max, block_min
   character(50) :: text_read, text_read2
   logical :: read_well

   ! Defaults:
   Error_descript = ''
   ind_S = 0 ! Element index
   ind_R = 0 ! Relative coordinates
   ind_V = 0 ! no velocities, only coordinates
   ind_atoms = -1  ! Atoms are set randomly
   SC_X = 0.0d0   ! ][A
   SC_Y = 0.0d0   ! [A]
   it_is_mixture = .false.  ! to start with, assume no mixture
   substitution_data%required = .false.  ! to start with, assume no substitution

   ! Start reading the line:
   count_lines = 0   ! to start with
   current_block = 0 ! to start with
   str_len = LEN(trim(line_2))  ! last symbol in the string of data

   ! Go through all the indices data and sort all of them out:
   strind:do
      ! Equal sign position:
      eq_pos = INDEX(line_2(current_block+1:str_len), '=')  ! intrinsic
      ! Read the variable:
      if (eq_pos > 0) then ! if there is "=" sign
         text_read = line_2(current_block+1:current_block+eq_pos-1)
      else  ! if there is no "=" sign, just check the entire line
         text_read = line_2(current_block+1:str_len)
      endif

      ! Block with data inside:
      block_start = INDEX(line_2(current_block+1:str_len), '"')  ! intrinsic
      block_end   = INDEX(line_2(current_block+block_start+1:str_len), '"')    ! intrinsic

      !print*, 'interpret_XYZ_comment_line 1:', trim(adjustl(text_read)), block_start, block_end, eq_pos
      !print*, 'interpret_XYZ_comment_line 2:', line_2(current_block+1:str_len)
      !print*, 'interpret_XYZ_comment_line 3:', line_2(current_block+block_start+1:current_block+block_start+block_end-1)

      ! Interprete the block with the data:
      select case (trim(adjustl(text_read)))
      !-----------------
      case ('Lattice', 'lattice', 'Supercell', 'supercell')
         ! Supercell vectors:
         !read(line_2(current_block+block_start+1:current_block+block_start+block_end-1),*,IOSTAT=Reason) Supce(:,:)
         read(line_2(current_block+block_start+1:current_block+block_start+block_end-1),*,IOSTAT=Reason) &
               Supce(1,1), Supce(1,2), Supce(1,3), &
               Supce(2,1), Supce(2,2), Supce(2,3), &
               Supce(3,1), Supce(3,2), Supce(3,3)  ! the order of reading: vector-wise instead of column-wise
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a)') 'Could not read block ', count_lines, ' in line #2 in XYZ file'
            print*, trim(adjustl(Error_descript))
            goto 2012
         endif

      !-----------------
      case ('Properties', 'properties')
         ind_atoms = 1  ! atoms are set in this file
         ! Read properties, if given:
         block_min = max(current_block+eq_pos+1, current_block+block_start+1)
         block_max = max(current_block+block_start+block_end-1 , str_len)
         call interprete_properties_line( line_2(block_min:block_max), ind_S, ind_R, ind_V ) ! below

      !-----------------
      case ('Random', 'random')
         ind_atoms = 0  ! atoms to be set at random
         if (eq_pos > 0) then ! X and Y sizes of supercell are specified
            call interpret_random_line(line_2(current_block+block_start+1:current_block+block_start+block_end-1) , &
                  SC_X, SC_Y) ! below
         endif

      !-----------------
      case ('Alloy', 'alloy', 'Mixture', 'mixture')
         ind_atoms = -1  ! atomic positions are set in this file
         it_is_mixture = .true.
         !print*, "We've got an alloy, everyone!"

      !-----------------
      case ('Substitution', 'substitution', 'Substitute', 'substitute')
         substitution_data%required = .true.
         ind_atoms = 1  ! atoms are set in this file (some to be substituted later)
         ! Read properties, if given:
         block_min = max(current_block+eq_pos+1, current_block+block_start+1)
         block_max = max(current_block+block_start+block_end-1 , str_len)
         call interprete_substitute_line( line_2(current_block+block_start+1:current_block+block_start+block_end-1), substitution_data) ! below

      !-----------------
      case default
         print*, 'Unknown flag in 2d line in xyz-file:', trim(adjustl(text_read))

      end select

      ! Next block, if exists:
      block2_end = INDEX(line_2(current_block+block_start+block_end+1:str_len), '=')    ! intrinsic
      current_block = current_block + block_start + block_end  ! mark these 2 blocks as done, next starts after this one
      if ((block_end == 0) .or. (block2_end == 0) .or. (current_block >= str_len)) exit strind
   enddo strind

2012 continue
end subroutine interpret_XYZ_comment_line



subroutine interprete_substitute_line(prop_block, substitution_data)
   character(*), intent(in) :: prop_block
   type(Substitute_data), intent(inout) :: substitution_data  ! flag for substitution
   !------------------------
   character(3) :: text_read1, text_read2
   real(8) :: real_read
   integer :: Reason, count_lines
   logical :: read_well

   ! Defaults to start with:
   substitution_data%Atoms_to_substitute = ""
   substitution_data%Atoms_substituting = ""
   substitution_data%percentage = 0.0d0

   ! Read the data:
   count_lines = 0
   read(prop_block,*,IOSTAT=Reason) text_read1, text_read2, real_read
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      print*, 'Problem in soubroutine interprete_substitute_line:'
      print*, 'Could not interprete block data in line #2: ', prop_block
      print*, 'Icorrect format: atoms to substitute, substituting atoms, and percentage'
      print*, 'No atomic substitution will be made!'
      substitution_data%required =.false.
   else  ! Interpret it:
      substitution_data%Atoms_to_substitute = trim(adjustl(text_read1))
      substitution_data%Atoms_substituting = trim(adjustl(text_read2))
      substitution_data%percentage = real_read
   endif

   !print*, substitution_data%required, substitution_data%Atoms_to_substitute, substitution_data%Atoms_substituting, substitution_data%percentage
end subroutine interprete_substitute_line


subroutine interpret_random_line(prop_block, SC_X, SC_Y)
   character(*), intent(in) :: prop_block
   real(8), intent(inout) :: SC_X, SC_Y
   !------------------------
   character(1) :: text_read1, text_read2
   real(8) :: read_val1, read_val2, eps
   integer :: Reason, count_lines
   logical :: read_well

   count_lines = 0
   read(prop_block,*,IOSTAT=Reason) text_read1, read_val1, text_read2, read_val2
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      print*, 'Problem in soubroutine interpret_random_line:'
      print*, 'Could not interprete block data in line #2: ', prop_block
      print*, 'Icorrect format in the data with X and Y sizes of simulation box'
      print*, 'Using cubic simulation box instead'
   else  ! Interpret it:
      select case (trim(adjustl(text_read1)))
      case ('X','x')
         SC_X = read_val1
      case ('Y','y')
         SC_Y = read_val1
      end select

      select case (trim(adjustl(text_read2)))
      case ('X','x')
         SC_X = read_val2
      case ('Y','y')
         SC_Y = read_val2
      end select
   endif

   ! Check for consistency:
   eps = 1.0d-10
   if ( ( abs(SC_Y) < eps ) .and. ( abs(SC_X) > eps ) ) then ! only one dimension (X) was set
      SC_Y = SC_X ! make it square
   endif
   if ( ( abs(SC_X) < eps ) .and. ( abs(SC_Y) > eps ) ) then ! only one dimension (Y) was set
      SC_X = SC_Y ! make it square
   endif

end subroutine interpret_random_line


subroutine interprete_properties_line(prop_block, ind_S, ind_R, ind_V)
   ! Properties allowed in OVITO (not all supported in XTANT yet):
   ! https://www.ovito.org/manual/reference/file_formats/input/xyz.html
   character(*), intent(in) :: prop_block
   integer, intent(inout) :: ind_S, ind_R, ind_V
   !------------------------
   character(20) :: text_read2, read_char
   integer :: colon_pos, current_block, block_len

   ! Defaults:
   ind_S = -1
   ind_R = -1
   ind_V = -1

   current_block = 1 ! to start
   colon_pos = INDEX(prop_block, ':') ! intrinsic
   block_len = LEN(prop_block)

   !print*, 'interprete_properties_line:', colon_pos, prop_block

   do while (colon_pos > 0)
      ! Read the marker:
      text_read2 = prop_block(current_block:current_block+colon_pos-2)
      read_char = prop_block(current_block+colon_pos:current_block+colon_pos)

      !print*, 'interprete_properties_line:', colon_pos, text_read2, read_char

      ! Interpret it:
      select case (trim(adjustl(text_read2)))
      case ('Species','species')
         if ( (trim(adjustl(read_char)) == 'S' ) .or. &
              (trim(adjustl(read_char)) == 's' ) ) then
            ind_S = 1   ! Element name is set
         else
            ind_S = 0   ! KOA index is set
         endif
      case ('Pos','pos')
         if ( (trim(adjustl(read_char)) == 'S' ) .or. &
              (trim(adjustl(read_char)) == 's' ) ) then
            ind_R = 0   ! Relative coordinate S is set
         else
            ind_R = 1   ! Absolute coordinate R is set
         endif
      case ('Vel','vel', 'Velo', 'velo')
         if ( (trim(adjustl(read_char)) == 'S' ) .or. &
              (trim(adjustl(read_char)) == 's' ) ) then
            ind_V = 0   ! Relative velocity VS is set
         else
            ind_V = 1   ! Absolute velocity V is set
         endif
      case ('Mass', 'mass')
         ! To read user-specified mass (unused)
      case ('Charge', 'charge')
         ! To read user-specified charge (unused)
      case ('kinetic_energy')
         ! To read user-specified charge (unused)
      end select

      ! Next colon:
      current_block = current_block + colon_pos  ! next part
      colon_pos = INDEX(prop_block(current_block:block_len), ':') ! intrinsi
   enddo
   !print*, ind_S, ind_R, ind_V
end subroutine interprete_properties_line


END MODULE Dealing_with_eXYZ
