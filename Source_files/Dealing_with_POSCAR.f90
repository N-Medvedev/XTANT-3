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
! This module contains subroutines for dealing with the POSCAR format
! https://www.vasp.at/wiki/index.php/POSCAR

MODULE Dealing_with_POSCAR
use Objects
use Dealing_with_files, only : read_file, number_of_columns
use Atomic_tools, only :  Coordinates_abs_to_rel, velocities_abs_to_rel

implicit none
PRIVATE

public :: read_POSCAR, get_KOA_from_element

 contains



subroutine read_POSCAR(FN_POSCAR, File_name_POSCAR, Scell, SCN, matter, numpar, Err) ! poscar format
   ! Description of the format: https://www.vasp.at/wiki/index.php/POSCAR
   integer, intent(in) :: FN_POSCAR ! extended XYZ file number (must be already open)
   character(*), intent(in) :: File_name_POSCAR ! extended XYZ file name
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: SCN ! number of the supercell (always =1)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------
   real(8) :: Scalfac(3), coord(3)
   integer :: count_lines, Reason, Nat, Ncol, i, kind_count, atoms_count, KOA
   integer, dimension(:), allocatable :: i_temp
   logical :: read_well
   character(1000) :: read_line, Error_descript
   character(10) :: ch_temp
   character(3), dimension(:), allocatable :: El_name
   character(*), parameter :: numbers = '0123456789'

   Scalfac(:) = 0.0d0   ! to start with
   count_lines = 0      ! to start with

   ! First line in POSCAR, comment:
   read(FN_POSCAR,*,IOSTAT=Reason) ! skip
   call read_file(Reason, count_lines, read_well)

   ! Second line in POSCAR, Scaling factor(s):
   read(FN_POSCAR,*,IOSTAT=Reason) read_line
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 4400
   endif
   Ncol = number_of_columns(read_line) ! module "Dealing_with_files"
   select case(Ncol)
   case (3)
      read(read_line,*,IOSTAT=Reason) Scalfac   ! three numbers
   case default
      read(read_line,*,IOSTAT=Reason) Scalfac(1)   ! one number
      Scalfac(2:3) = Scalfac(1)
   end select

   ! Lines 3-5 in POSCAR, Lattice:
   read(FN_POSCAR,*,IOSTAT=Reason) coord(:)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 4400
   endif
   Scell(SCN)%Supce(1,:) = coord(:)  ! three numbers

   read(FN_POSCAR,*,IOSTAT=Reason) coord(:)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 4400
   endif
   Scell(SCN)%Supce(2,:) = coord(:)  ! three numbers

   read(FN_POSCAR,*,IOSTAT=Reason) coord(:)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 4400
   endif
   Scell(SCN)%Supce(3,:) = coord(:)  ! three numbers

   ! Scaling factor:
   Scell(SCN)%Supce(:,1) = Scell(SCN)%Supce(:,1) * Scalfac(1)
   Scell(SCN)%Supce(:,2) = Scell(SCN)%Supce(:,2) * Scalfac(2)
   Scell(SCN)%Supce(:,3) = Scell(SCN)%Supce(:,3) * Scalfac(3)

   ! Previous step:
   Scell(SCN)%Supce0 = Scell(SCN)%Supce   ! initial

   ! Line 6 in POSCAR, Species names (optional):
   read_line = '' ! to restart
   read(FN_POSCAR, '(a200)', IOSTAT=Reason) read_line ! read names or numbers of species
   !print*, Reason, 'read_line', trim(adjustl(read_line))
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 4400
   endif
   ! Check if the line contains optional names, or the mandatory numbers:
   !if (verify(trim(adjustl(read_line(1:1))), trim(adjustl(numbers))) /= 0) then ! it was optional line with names
   read(read_line,*,IOSTAT=Reason) ch_temp
   if (verify(trim(adjustl(ch_temp(1:1))), trim(adjustl(numbers))) /= 0) then ! it was optional line with names
      Ncol = number_of_columns(trim(adjustl(read_line))) ! module "Dealing_with_files"
      allocate(El_name(Ncol))
      read(read_line,*,IOSTAT=Reason) El_name(:) ! read kinds of species

      ! Next line, numbers of atoms:
      read(FN_POSCAR,'(a200)',IOSTAT=Reason) read_line ! read numbers of species
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 4400
      endif
   endif
   ! Having the numbers of species, get the total number of atoms:
   Ncol = number_of_columns(trim(adjustl(read_line))) ! module "Dealing_with_files"
   allocate(i_temp(Ncol), source = 0)
   read(read_line,*,IOSTAT=Reason) i_temp(:) ! read numbers of species
   !print*, 'Ncol', Ncol
   Nat = SUM(i_temp) ! count all the atoms
   !print*, 'Nat=', Nat, size(i_temp), i_temp, ': ', read_line
   !pause 'poscar'
   Scell(SCN)%Na = Nat
   allocate(Scell(SCN)%MDAtoms(Nat))

   ! Selective dynamics (optional) OR specification of coordinate format:
   read(FN_POSCAR,*,IOSTAT=Reason) read_line ! read names or numbers of species
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 4400
   endif
   ! Check if it was optional line or not
   select case(read_line(1:1))
   case('S', 's') ! read the next mandatory line with specification of coordinate format
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 4400
      endif
   endselect
   ! Check what format coordinates are given in (relative or absolute):

   ! Next lines contain at least the atomic type and coordinates:
   kind_count = 1
   atoms_count = i_temp(kind_count)
   do i = 1, Scell(SCN)%Na
      if (i > atoms_count) then
         kind_count = kind_count + 1   ! next index of kinds of atoms
         atoms_count = atoms_count + i_temp(kind_count)  ! next boundary for this kind of atoms
      endif

      ! Define kind of atom:
      if (allocated(El_name)) then ! element names are given
         ! Find the index from the element name:
         call get_KOA_from_element(El_name(kind_count), matter, KOA) ! below
         if (KOA <= 0) then
            write(Error_descript,'(a,i3,a,$)') 'In the target, there is no element ', trim(adjustl(El_name(kind_count))), &
                                               ' from file '//trim(adjustl(File_name_POSCAR))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 4400
         endif
         Scell(SCN)%MDAtoms(i)%KOA = KOA
      else ! assume elements follow from the chemical formula
         Scell(SCN)%MDAtoms(i)%KOA = kind_count
      endif

      ! Define atomic coordinates:
      read(FN_POSCAR,*,IOSTAT=Reason) coord(:)
      if (.not. read_well) then
         write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_POSCAR))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 4400
      endif

      select case(read_line(1:1))
      case ('C', 'c', 'K', 'k') ! Absolute
         Scell(SCN)%MDAtoms(i)%R(:) = coord(:) * Scalfac(:)
         Scell(SCN)%MDAtoms(i)%R0(:) = Scell(SCN)%MDAtoms(i)%R(:)
         ! Get the relative coordinates from the absolute ones provided:
         call Coordinates_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
      case default ! Relative
         Scell(SCN)%MDAtoms(i)%S(:) = coord(:)
         Scell(SCN)%MDAtoms(i)%S0(:) = Scell(SCN)%MDAtoms(i)%S(:)
      endselect
      !print*, i, Scell(SCN)%MDAtoms(i)%KOA, Scell(SCN)%MDAtoms(i)%S(:)
   enddo

   ! Optional lines from here (incl. velosities)
   read_well = .true.   ! to start with
   RDID: do while (read_well)
      read(FN_POSCAR,*,IOSTAT=Reason) read_line ! read names or numbers of species
      call read_file(Reason, count_lines, read_well)

      ! Sort the atomic velocities that were read:
      select case(trim(adjustl(read_line)))
      case ('CARTESIAN', 'Cartesian', 'cartesian') ! Absolute
         do i = 1, Scell(SCN)%Na
            read(FN_POSCAR,*,IOSTAT=Reason) Scell(SCN)%MDAtoms(i)%V(:)
            call read_file(Reason, count_lines, read_well)
            Scell(SCN)%MDAtoms(i)%V0(:) = Scell(SCN)%MDAtoms(i)%V(:)
            ! Get the relative velocities from the absolute ones provided:
            call velocities_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
         enddo
      case('DIRECT', 'Direct', 'direct')  ! relative
         do i = 1, Scell(SCN)%Na
            read(FN_POSCAR,*,IOSTAT=Reason) Scell(SCN)%MDAtoms(i)%SV(:)
            call read_file(Reason, count_lines, read_well)
            Scell(SCN)%MDAtoms(i)%SV0(:) = Scell(SCN)%MDAtoms(i)%SV(:)
         enddo
      end select
   enddo RDID

   if (allocated(El_name)) deallocate(El_name)
   deallocate(i_temp)
4400 continue
!    pause 'read_POSCAR'
end subroutine read_POSCAR



subroutine get_KOA_from_element(El_name, matter, KOA)
   character(3), intent(in) :: El_name ! name of the element
   type(Solid), intent(inout) :: matter   ! all material parameters
   integer ,intent(out) :: KOA   ! number of element from the given array
   !--------------------------
   integer :: i
   KOA = 0  ! to start with
   do i = 1, size(matter%Atoms)
      if ( trim(adjustl(matter%Atoms(i)%Name)) == trim(adjustl(El_name)) ) then
         KOA = i  ! that's the element in our array
         exit
      endif
   enddo
end subroutine get_KOA_from_element


END MODULE Dealing_with_POSCAR
