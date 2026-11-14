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
! This module contains subroutines for dealing with the POSCAR format
! https://www.vasp.at/wiki/index.php/POSCAR

MODULE Dealing_with_POSCAR
use Objects
use Dealing_with_files, only : read_file, number_of_columns
use Atomic_tools, only :  Coordinates_abs_to_rel, velocities_abs_to_rel
use Algebra_tools, only : Det_3x3

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
   real(8) :: Scalfac(3), coord(3), cell_volume, scale_factor
   integer :: count_lines, Reason, Nat, Ncol, i, kind_count, atoms_count, KOA, velocity_mode
   integer, dimension(:), allocatable :: i_temp
   logical :: read_well, cartesian_coordinates, allocated_atoms
   character(1000) :: read_line, Error_descript
   character(10) :: ch_temp
   character(3), dimension(:), allocatable :: El_name
   character(*), parameter :: numbers = '0123456789'

   Scalfac = 1.0d0            ! to start with
   count_lines = 0            ! to start with
   allocated_atoms = .false.  ! to start with

   !print*, "read_POSCAR"
   !pause 'Starting to read_POSCAR'

   ! First line in POSCAR, comment:
   !read(FN_POSCAR,*,IOSTAT=Reason) ! skip
   read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line   ! skip
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read the comment line')
      goto 4400
   endif
   !print*, "read_POSCAR", Reason, read_well

   ! Second line in POSCAR, Scaling factor(s):
   read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read the scaling factors')
      goto 4400
   endif
   read_line = adjustl(read_line)

   Ncol = number_of_columns(read_line) ! module "Dealing_with_files"
   if (Ncol /= 1 .and. Ncol /= 3) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Expected one or three scaling factors')
      goto 4400
   endif
   if (Ncol == 1) then
      read(read_line,*,IOSTAT=Reason) Scalfac(1)
   else
      read(read_line,*,IOSTAT=Reason) Scalfac
   endif
   if (Reason /= 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid scaling factor')
      goto 4400
   endif


   ! Lines 3-5 in POSCAR, Lattice:
   do i = 1, 3
      read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read a lattice vector')
         goto 4400
      endif
      read(read_line,*,IOSTAT=Reason) coord
      if (Reason /= 0) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid lattice vector')
         goto 4400
      endif
      Scell(SCN)%Supce(i,:) = coord
   enddo

   ! Get the volume:
   call Det_3x3(Scell(SCN)%Supce, Scell(SCN)%V)
   cell_volume = abs(Scell(SCN)%V)
   if (cell_volume <= tiny(cell_volume)) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Lattice vectors define a zero-volume cell')
      goto 4400
   endif

   ! Apply scaling:
   if (Ncol == 1) then  ! Single scaling factor
      if (Scalfac(1) == 0.0d0) then
         call report_poscar_error(Err, File_name_POSCAR, 2, 'The scalar scaling factor must not be zero')
         goto 4400
      else if (Scalfac(1) < 0.0d0) then
         scale_factor = (abs(Scalfac(1))/cell_volume)**(1.0d0/3.0d0)
         Scalfac = scale_factor
      else
         Scalfac(2:3) = Scalfac(1)
      endif
   else if (any(Scalfac <= 0.0d0)) then
      call report_poscar_error(Err, File_name_POSCAR, 2, 'Three scaling factors must be positive')
      goto 4400
   endif

   ! Scaling factor:
   !Scell(SCN)%Supce(:,1) = Scell(SCN)%Supce(:,1) * Scalfac(1)
   !Scell(SCN)%Supce(:,2) = Scell(SCN)%Supce(:,2) * Scalfac(2)
   !Scell(SCN)%Supce(:,3) = Scell(SCN)%Supce(:,3) * Scalfac(3)
   do i = 1, 3
      Scell(SCN)%Supce(i,:) = Scell(SCN)%Supce(i,:) * Scalfac(i)
   enddo
   ! Update volume:
   call Det_3x3(Scell(SCN)%Supce, Scell(SCN)%V)
   ! Previous step:
   Scell(SCN)%Supce0 = Scell(SCN)%Supce   ! initial

   ! Line 6 in POSCAR, Species names (optional):
   read_line = '' ! to restart
   read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line ! read names or numbers of species
   !print*, Reason, 'read_line', trim(adjustl(read_line))
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read species names or atom counts')
      goto 4400
   endif
   read_line = adjustl(read_line)

   Ncol = number_of_columns(trim(adjustl(read_line))) ! module "Dealing_with_files"
   if (Ncol == 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Species or atom-count line is empty')
      goto 4400
   endif

   ! Check if the line contains optional names, or the mandatory numbers:
   !if (verify(trim(adjustl(read_line(1:1))), trim(adjustl(numbers))) /= 0) then ! it was optional line with names
   read(read_line,*,IOSTAT=Reason) ch_temp
   if (Reason /= 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid species or atom-count line')
      goto 4400
   endif

   if (verify(ch_temp(1:1), numbers) /= 0 .and. ch_temp(1:1) /= '+' .and. ch_temp(1:1) /= '-') then
      allocate(El_name(Ncol), stat=Reason)
      if (Reason /= 0) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not allocate species names')
         goto 4400
      endif
      read(read_line,*,IOSTAT=Reason) El_name   ! read kinds of species
      if (Reason /= 0) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid species-name line')
         goto 4400
      endif

      ! Next line, numbers of atoms:
      read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read atom counts')
         goto 4400
      endif
      read_line = adjustl(read_line)
   endif


   ! Having the numbers of species, get the total number of atoms:
   Ncol = number_of_columns(trim(adjustl(read_line))) ! module "Dealing_with_files"
   if (Ncol == 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Atom-count line is empty')
      goto 4400
   endif
   if (allocated(El_name)) then
      if (size(El_name) /= Ncol) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, &
                                  'Species-name and atom-count columns do not match')
         goto 4400
      endif
!    else if (Ncol > size(matter%Atoms)) then
!       call report_poscar_error(Err, File_name_POSCAR, count_lines, &
!                                'Atom-count line has more species than the configured material')
!       goto 4400
   endif


   !allocate(i_temp(Ncol), source = 0)
   allocate(i_temp(Ncol), stat=Reason)
   if (Reason /= 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not allocate atom counts')
      goto 4400
   endif

   read(read_line,*,IOSTAT=Reason) i_temp(:) ! read numbers of species
   if (Reason /= 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid atom-count line')
      goto 4400
   endif
   if (any(i_temp < 0)) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Atom counts must not be negative')
      goto 4400
   endif

   !print*, 'Ncol', Ncol
   !Nat = SUM(i_temp) ! count all the atoms
   Nat = 0
   do i = 1, Ncol
      if (i_temp(i) > huge(Nat) - Nat) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Total atom count is too large')
         goto 4400
      endif
      Nat = Nat + i_temp(i)
   enddo
   if (Nat <= 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'POSCAR must contain at least one atom')
      goto 4400
   endif

   !print*, 'Nat=', Nat, size(i_temp), i_temp, ': ', read_line
   !pause 'poscar'
   allocate(Scell(SCN)%MDAtoms(Nat), stat=Reason)
   if (Reason /= 0) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not allocate atom data')
      goto 4400
   endif
   allocated_atoms = .true.
   Scell(SCN)%Na = Nat


   ! Selective dynamics (optional) OR specification of coordinate format:
   read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line ! read names or numbers of species
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read coordinate mode')
      goto 4400
   endif
   read_line = adjustl(read_line)
   ! Check if it was optional line or not
   select case(read_line(1:1))
   case ('S', 's')   ! needs to read next line with the coordinates
      read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read coordinate mode after selective dynamics')
         goto 4400
      endif
      read_line = adjustl(read_line)
   end select
   !definition of coordinates format:
   select case(read_line(1:1))
   case ('C', 'c', 'K', 'k')
      cartesian_coordinates = .true.
   case ('D', 'd')
      cartesian_coordinates = .false.
   case default
      call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Unknown POSCAR coordinate mode')
      goto 4400
   end select


   ! Check what format coordinates are given in (relative or absolute):

   ! Next lines contain at least the atomic type and coordinates:
   kind_count = 1
   atoms_count = i_temp(kind_count)
   do i = 1, Scell(SCN)%Na
      !if (i > atoms_count) then
      !   kind_count = kind_count + 1   ! next index of kinds of atoms
      !   atoms_count = atoms_count + i_temp(kind_count)  ! next boundary for this kind of atoms
      !endif
      do while (i > atoms_count .and. kind_count < Ncol)
         kind_count = kind_count + 1   ! next index of kinds of atoms
         atoms_count = atoms_count + i_temp(kind_count)  ! next boundary for this kind of atoms
      enddo

      ! Define kind of atom:
      if (allocated(El_name)) then ! element names are given
         ! Find the index from the element name:
         call get_KOA_from_element(El_name(kind_count), matter, KOA) ! below
         if (KOA <= 0) then
            call report_poscar_error(Err, File_name_POSCAR, count_lines, &
                                     'POSCAR species is not present in the configured material')
            goto 4400
         endif
         Scell(SCN)%MDAtoms(i)%KOA = KOA
      else ! assume elements follow from the chemical formula
         Scell(SCN)%MDAtoms(i)%KOA = kind_count
      endif

      ! Define atomic coordinates:
      !read(FN_POSCAR,*,IOSTAT=Reason) coord(:)
      read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read atomic coordinates')
         goto 4400
      endif

      read(read_line,*,IOSTAT=Reason) coord
      if (Reason /= 0) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid atomic coordinates')
         goto 4400
      endif

      if (cartesian_coordinates) then ! Absolute
         Scell(SCN)%MDAtoms(i)%R(:) = coord(:) * Scalfac(:)
         Scell(SCN)%MDAtoms(i)%R0(:) = Scell(SCN)%MDAtoms(i)%R(:)
         ! Get the relative coordinates from the absolute ones provided:
         call Coordinates_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
      else ! Relative
         Scell(SCN)%MDAtoms(i)%S(:) = coord(:)
         Scell(SCN)%MDAtoms(i)%S0(:) = Scell(SCN)%MDAtoms(i)%S(:)
      endif
      !print*, i, Scell(SCN)%MDAtoms(i)%KOA, Scell(SCN)%MDAtoms(i)%S(:)
   enddo


   ! Optional lines from here (incl. velosities)
   velocity_mode = 0
   do
      read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
      if (Reason < 0) exit
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read optional velocity data')
         goto 4400
      endif
      read_line = adjustl(read_line)
      if (len_trim(read_line) == 0) cycle
      select case(read_line(1:1))
      case ('C', 'c')
         velocity_mode = 1
         exit
      case ('D', 'd')
         velocity_mode = 2
         exit
      end select
   enddo

   if (velocity_mode == 1 .or. velocity_mode == 2) then
      do i = 1, Scell(SCN)%Na
         read(FN_POSCAR,'(a)',IOSTAT=Reason) read_line
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Could not read atomic velocities')
            goto 4400
         endif
         read(read_line,*,IOSTAT=Reason) coord
         if (Reason /= 0) then
            call report_poscar_error(Err, File_name_POSCAR, count_lines, 'Invalid atomic velocities')
            goto 4400
         endif
         if (velocity_mode == 1) then
            Scell(SCN)%MDAtoms(i)%V = coord
            Scell(SCN)%MDAtoms(i)%V0 = coord
         else
            Scell(SCN)%MDAtoms(i)%SV = coord
            Scell(SCN)%MDAtoms(i)%SV0 = coord
         endif
      enddo
      ! Get the relative velocities from the absolute ones provided:
      if (velocity_mode == 1) call velocities_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
   endif

   if (allocated(El_name)) deallocate(El_name)
   if (allocated(i_temp)) deallocate(i_temp)

   !pause 'read_POSCAR'

   return   ! if sucsessful, exit

4400 continue  ! if not sucsessful, deallocated MDAtoms and reset Scell(SCN)%Na
   if (allocated(El_name)) deallocate(El_name)
   if (allocated(i_temp)) deallocate(i_temp)
   if (allocated_atoms) then
      deallocate(Scell(SCN)%MDAtoms)
      Scell(SCN)%Na = 0
   endif

end subroutine read_POSCAR



subroutine get_KOA_from_element(El_name, matter, KOA)
   character(3), intent(in) :: El_name ! name of the element
   type(Solid), intent(in) :: matter   ! all material parameters
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




subroutine report_poscar_error(Err, File_name_POSCAR, line_number, message)
   type(Error_handling), intent(inout) :: Err
   character(*), intent(in) :: File_name_POSCAR, message
   integer, intent(in) :: line_number
   character(1000) :: Error_descript

   write(Error_descript,'(a,i0,a,a)') trim(message)//' on line ', line_number, ' in file ', trim(File_name_POSCAR)
   call Save_error_details(Err, 3, trim(Error_descript))
   print*, trim(Error_descript)
end subroutine report_poscar_error


END MODULE Dealing_with_POSCAR
