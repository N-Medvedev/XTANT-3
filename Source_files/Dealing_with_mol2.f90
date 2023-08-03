! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2023 Nikita Medvedev
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
! This module contains subroutines for dealing with the mol2 format
! https://www.structbio.vanderbilt.edu/archives/amber-archive/2007/att-1568/01-mol2_2pg_113.pdf

MODULE Dealing_with_mol2
use Objects, only: Super_cell, Solid, Numerics_param, Error_handling, Save_error_details
use Universal_Constants, only : g_Pi
use Dealing_with_files, only : read_file, number_of_columns
use Atomic_tools, only :  Coordinates_abs_to_rel, check_periodic_boundaries, cell_vectors_defined_by_angles
use Dealing_with_POSCAR, only : get_KOA_from_element

implicit none
PRIVATE

public :: read_mol2

 contains



subroutine read_mol2(FN, File_name, Scell, SCN, matter, numpar, Err) ! basic mol2 format
   integer, intent(in) :: FN  ! file with mol2 extension to read from (must be opened)
   character(*), intent(in) :: File_name  ! name of this file
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: SCN ! number of the supercell (always =1)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------
   real(8) :: a, b, c, alpha, beta, gamm, a_vec(3), b_vec(3), c_vec(3)
   integer :: count_lines, Reason, Nat, i, temp_i, KOA
   logical :: read_well, read_atoms, read_supercell
   character(1000) :: read_line, Error_descript
   character(30) :: cut_line, read_marker
   character(10) :: El_name_read
   character(3) :: El_name
   character(*), parameter :: comment_flag = '#'   ! comments in mol2 are starting with it
   character(*), parameter :: marker_flag = '@' ! functional markers in mol2 are starting with it

   count_lines = 0      ! to start with
   read_atoms = .false. ! didn't read yet
   read_supercell  = .false. ! didn't read yet

   ! First line in mol2, comment:
   readmol2:do ! read until the end of file
      read(FN,'(a)',IOSTAT=Reason) read_line
      call read_file(Reason, count_lines, read_well)
      if (Reason .LT. 0) exit readmol2
      if (.not. read_well) then
         write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         exit readmol2
      endif

      select case ( trim(adjustl(read_line(1:1))) )
      case (marker_flag)
         ! Here is something useful:
         cut_line = trim(adjustl(read_line(10:)))
         read(cut_line,*,IOSTAT=Reason) read_marker
         count_lines = count_lines - 1    ! reading the same line
         call read_file(Reason, count_lines, read_well)
         !print*, trim(adjustl(read_marker)) ! test reading line

         select case ( trim(adjustl(read_marker)) )
         !----------------------------
         case ('MOLECULE', 'molecule', 'Molecule') ! to get number of atoms
            read(FN, '(a)', IOSTAT=Reason) read_line  ! skip the name of the molecule
            if (.not. read_well) then
               write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               exit readmol2
            endif

            read(FN,*,IOSTAT=Reason) Nat ! number of atoms
            if (.not. read_well) then
               write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               exit readmol2
            endif

            Scell(SCN)%Na = Nat  ! knowng the number of atoms, save it and use it
            allocate(Scell(SCN)%MDAtoms(Nat))

         !----------------------------
         case ('ATOM', 'Atom', 'atom') ! to get atomic species and coordinates
            ! Read atomis parameters (element name and coordinates):
            do i = 1, Scell(SCN)%Na
               read(FN,*,IOSTAT=Reason) temp_i, El_name_read, Scell(SCN)%MDAtoms(i)%R(:)
               if (.not. read_well) then
                  write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
                  call Save_error_details(Err, 3, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  exit readmol2
               endif
               ! Make sure it is just the element name:
               call trim_element_name(El_name_read, El_name)   ! below

               ! Get the internal element number:
               call get_KOA_from_element(El_name, matter, KOA) ! module "Dealing_with_POSCAR"
               if (KOA <= 0) then
                  write(Error_descript,'(a,i3,a,$)') 'In the target, there is no element ', trim(adjustl(El_name)), &
                                                      ' from file '//trim(adjustl(File_name))
                  call Save_error_details(Err, 3, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  exit readmol2
               endif
               Scell(SCN)%MDAtoms(i)%KOA = KOA

               !print*, i, El_name, KOA, Scell(SCN)%MDAtoms(i)%R(:)
            enddo
            read_atoms = .true.  ! we read atomic data

         !----------------------------
         case ('CRYSIN', 'crysin', 'Crysin') ! to get supercell sizes
            read(FN,*,IOSTAT=Reason) a, b, c, alpha, beta, gamm
            if (.not. read_well) then
               write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               exit readmol2
            endif

            ! Convert angles from deg to rad:
            alpha = alpha*g_Pi/180.0d0
            beta = beta*g_Pi/180.0d0
            gamm = gamm*g_Pi/180.0d0
            ! Construct supercell vectors from the absolute values and angles:
            call cell_vectors_defined_by_angles(a, b, c, alpha, beta, gamm, a_vec, b_vec, c_vec, Reason)   ! module "Atomic_tools"
            if (Reason /= 0) then
               write(Error_descript,'(a)') 'Supercell parameters are wrong in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               exit readmol2
            endif
            ! Save supercell vectors:
            Scell(SCN)%Supce(:,1) = a_vec
            Scell(SCN)%Supce(:,2) = b_vec
            Scell(SCN)%Supce(:,3) = c_vec

            read_supercell = .true. ! we read supercell data
         end select

      case (comment_flag)
         ! skip comment line
      case default
         ! also skip this line
      end select
   enddo readmol2

   ! Make sure we read all we need:
   if (.not.read_atoms) then
      write(Error_descript,'(a)') 'Could not find atomic coordinates in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 4, Error_descript)
      print*, trim(adjustl(Error_descript))
      return
   endif

   if (.not.read_supercell) then
      write(Error_descript,'(a)') 'Could not find supercell vectors in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 4, Error_descript)
      print*, trim(adjustl(Error_descript))
      return
   endif

   ! Previous step:
   Scell(SCN)%Supce0 = Scell(SCN)%Supce   ! initial

   ! Get the relative coordinates from the absolute ones provided:
   call Coordinates_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"

   ! Make sure all atoms are inside the box
   call check_periodic_boundaries(matter, Scell, SCN) ! module "Atomic_tools"

!    ! Test:
!    print*, 'Supce:', Scell(SCN)%Supce
!    do i = 1, Scell(SCN)%Na
!        print*, 'S', i, Scell(SCN)%MDAtoms(i)%KOA, Scell(SCN)%MDAtoms(i)%S(:)
!        print*, 'R', i, Scell(SCN)%MDAtoms(i)%KOA, Scell(SCN)%MDAtoms(i)%R(:)
!    enddo
!    pause 'Dealing_with_mol2'
end subroutine read_mol2


subroutine trim_element_name(El_name_read, El_name)
   character(*), intent(in) :: El_name_read  ! read name from file
   character(*), intent(out) :: El_name   ! make sure it is just the chemical element name
   !----------------------------
   integer :: i, name_length, counter
   character(*), parameter :: numbers = '0123456789'
   El_name = ''   ! to start with
   name_length = LEN(trim(adjustl(El_name_read)))
   counter = 0
   do i = 1, name_length
      if (verify(trim(adjustl(El_name_read(i:i))), trim(adjustl(numbers))) /= 0) then ! it is not a number
         counter = counter + 1
         El_name(counter:counter) = trim(adjustl(El_name_read(i:i)))
      else ! it is a number
         ! skip it
      endif
   enddo
end subroutine trim_element_name



END MODULE Dealing_with_mol2
