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
! This module contains subroutines for dealing with CDF-format file

MODULE Dealing_with_CDF
use Universal_Constants   ! let it use universal constants
use Objects, only : Save_error_details, Solid, Numerics_param, Error_handling, At_data
use Periodic_table, only : Decompose_compound
use Dealing_with_files, only : read_file
use Dealing_with_EADL, only : READ_EADL_TYPE_FILE_real, READ_EADL_TYPE_FILE_int, select_imin_imax

implicit none
PRIVATE

public :: read_CDF_file, write_CDF_file

 contains


subroutine read_CDF_file(FN, matter, numpar, Err, File_name, file_atomic_parameters, file_EADL_file, INFO)
   integer, intent(in) :: FN  ! file number (must be already opened)
   type(Solid), intent(inout) :: matter   ! all material parameters
   type(Numerics_param), intent(inout) :: numpar   ! all numerical parameters
   type(Error_handling), intent(inout) :: Err   ! error save
   character(*), intent(in) :: File_name, file_atomic_parameters, file_EADL_file ! files with CDF, atomic parameters, and EADL
   integer, intent(out) :: INFO  ! error flag
   !===============================================
   integer, dimension(:), allocatable :: at_numbers
   real(8), dimension(:), allocatable :: at_percentage
   character(3), dimension(:), allocatable :: at_short_names ! name of the element
   character(25), dimension(:), allocatable :: at_names ! full name of the element
   real(8), dimension(:), allocatable :: at_masses ! mass of each element [Mp]
   integer, dimension(:), allocatable :: at_NVE    ! number of valence electrons
   integer FN1, i, j, k, Z
   character(200) :: error_message, Error_descript
   character(200) :: Folder_name, Folder_name2
   logical file_exist, file_opened, read_well
   integer Reason, count_lines, temp
   real(8) retemp
   character(11) :: ch_nul
   parameter (ch_nul = '')

   INFO = 0 ! no errors to start with
   count_lines = 0   ! start counting
   read(FN,*,IOSTAT=Reason)	! skip first line with the name of the material
   call read_file(Reason, count_lines, read_well)  ! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 3, Error_descript) ! modlue "Objects"
      print*, trim(adjustl(Error_descript))
      INFO = 1   ! error in reading file
      return ! exit the subroutine
   endif

   read(FN,*,IOSTAT=Reason) matter%Chem ! chemical formula of the material
   call read_file(Reason, count_lines, read_well) ! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 3, Error_descript) ! modlue "Objects"
      print*, trim(adjustl(Error_descript))
      INFO = 1   ! error in reading file
      return ! exit the subroutine
   endif

   Folder_name2 = trim(adjustl(numpar%input_path))//trim(adjustl(file_atomic_parameters))  ! 'Atomic_parameters'
   call Decompose_compound(Folder_name2, matter%Chem, numpar%path_sep, INFO, error_message, matter%N_KAO, at_numbers, &
               at_percentage, at_short_names, at_names, at_masses, at_NVE) ! molude 'Periodic_table'
   if (INFO .NE. 0) then
      call Save_error_details(Err, INFO, error_message) ! modlue "Objects"
      print*, trim(adjustl(error_message))
      INFO = 1   ! error in reading file
      return ! exit the subroutine
   endif
   if (.not.allocated(matter%Atoms)) allocate(matter%Atoms(matter%N_KAO))

   do i = 1, matter%N_KAO ! for all sorts of atoms
      matter%Atoms(i)%Z = at_numbers(i)
      matter%Atoms(i)%Name = at_short_names(i)
      matter%Atoms(i)%Ma = at_masses(i)*g_Mp ! [kg]
      matter%Atoms(i)%percentage = at_percentage(i)
      matter%Atoms(i)%NVB = at_NVE(i)
      !print*, matter%Atoms(i)%Name, matter%Atoms(i)%Ma, matter%Atoms(i)%NVB
      !pause 'get_CDF_data'
   enddo

   read(FN,*,IOSTAT=Reason) retemp ! skip this line - density is given elsewhere
   call read_file(Reason, count_lines, read_well) ! module "Dealing_with_files"
   if (.not. read_well) then
      write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 3, Error_descript) ! modlue "Objects"
      print*, trim(adjustl(Error_descript))
      INFO = 1   ! error in reading file
      return ! exit the subroutine
   endif
   if (matter%dens .LE. 0.0d0) matter%dens = retemp ! in case density is not given elsewhere

   AT_NUM:do i = 1, matter%N_KAO ! for each kind of atoms:
      read(FN,*,IOSTAT=Reason) matter%Atoms(i)%sh  ! number of shells in this element
      call read_file(Reason, count_lines, read_well) ! module "Dealing_with_files"
      if (.not. read_well) then
         write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         call Save_error_details(Err, 3, Error_descript) ! modlue "Objects"
         print*, trim(adjustl(Error_descript))
         INFO = 1   ! error in reading file
         return ! exit the subroutine
      endif

      if (.not. allocated(matter%Atoms(i)%Ne_shell)) allocate(matter%Atoms(i)%Ne_shell(matter%Atoms(i)%sh), source=-0.5d0) ! Ne
      if (.not. allocated(matter%Atoms(i)%Shell_name)) allocate(matter%Atoms(i)%Shell_name(matter%Atoms(i)%sh), source=ch_nul) ! shell names
      if (.not. allocated(matter%Atoms(i)%Shl_dsgnr)) allocate(matter%Atoms(i)%Shl_dsgnr(matter%Atoms(i)%sh),source=-1) ! shell disignator
      if (.not. allocated(matter%Atoms(i)%N_CDF)) allocate(matter%Atoms(i)%N_CDF(matter%Atoms(i)%sh), source=0) ! #CDF functions
      if (.not. allocated(matter%Atoms(i)%Ip)) allocate(matter%Atoms(i)%Ip(matter%Atoms(i)%sh), source=-1.0d0) ! ionization potentials
      if (.not. allocated(matter%Atoms(i)%Ek)) allocate(matter%Atoms(i)%Ek(matter%Atoms(i)%sh), source=0.0d0) ! Ekin of shells
      if (.not. allocated(matter%Atoms(i)%TOCS)) allocate(matter%Atoms(i)%TOCS(matter%Atoms(i)%sh), source=0) ! type of cross-section
      if (.not. allocated(matter%Atoms(i)%TOCSph)) allocate(matter%Atoms(i)%TOCSph(matter%Atoms(i)%sh), source=0) ! type of photon cross-section
      if (.not. allocated(matter%Atoms(i)%El_MFP)) allocate(matter%Atoms(i)%El_MFP(matter%Atoms(i)%sh)) ! electron MFPs
      if (.not. allocated(matter%Atoms(i)%Ph_MFP)) allocate(matter%Atoms(i)%Ph_MFP(matter%Atoms(i)%sh)) ! photon MFPs
      matter%Atoms(i)%Ek = 0.0d0 ! starting value
      !if (.not. allocated(matter%Atoms(i)%Radiat)) then ! currently unused, maybe implemented later analogous to TREKIS
      !   allocate(matter%Atoms(i)%Radiat(matter%Atoms%sh)) ! allocate Radiative-times
      !   matter%Atoms(i)%Radiat = 1d24 ! to start with
      !endif
      if (.not. allocated(matter%Atoms(i)%Auger)) then
         allocate(matter%Atoms(i)%Auger(matter%Atoms(i)%sh), source=-1.2d24) ! Auger times
      endif
      if (.not. allocated(matter%Atoms(i)%CDF)) allocate(matter%Atoms(i)%CDF(matter%Atoms(i)%sh)) ! CDF functions

      do j = 1, matter%Atoms(i)%sh	! for all shells:
         ! Number of CDF functions, shell-designator, ionization potential, number of electrons, Auger-time:
         read(FN,*,IOSTAT=Reason) matter%Atoms(i)%N_CDF(j), matter%Atoms(i)%Shl_dsgnr(j), matter%Atoms(i)%Ip(j), &
                                    matter%Atoms(i)%Ne_shell(j), matter%Atoms(i)%Auger(j)
         call read_file(Reason, count_lines, read_well) ! module "Dealing_with_files"
         if (.not. read_well) then
            write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript) ! modlue "Objects"
            print*, trim(adjustl(Error_descript))
            INFO = 1   ! error in reading file
            return ! exit the subroutine
         endif
         matter%Atoms(i)%Shell_name(j) = ''  ! no name yet, to be defined in the next subroutine

         ! Check if any atomic data are undefined or unphysical, use those from the EPICS-database:
         call check_atomic_data(matter, numpar, Err, i, j, matter%Atoms(i)%sh, file_atomic_parameters, file_EADL_file)  ! below

         DOCDF:if (matter%Atoms(i)%N_CDF(j) .GT. 0) then ! do this shell with provided CDF coefficients
            matter%Atoms(i)%TOCS(j) = 1 ! CDF cross-section
            matter%Atoms(i)%TOCSph(j) = 1 ! CDF cross-section
            allocate(matter%Atoms(i)%CDF(j)%A(matter%Atoms(i)%N_CDF(j)))
            allocate(matter%Atoms(i)%CDF(j)%E0(matter%Atoms(i)%N_CDF(j)))
            allocate(matter%Atoms(i)%CDF(j)%G(matter%Atoms(i)%N_CDF(j)))

            do k = 1, matter%Atoms(i)%N_CDF(j)  ! for all CDF-functions for this shell
               read(FN,*,IOSTAT=Reason) matter%Atoms(i)%CDF(j)%E0(k), matter%Atoms(i)%CDF(j)%A(k), matter%Atoms(i)%CDF(j)%G(k)
!                 write(*,*) matter%Atoms(i)%CDF(j)%E0(k), matter%Atoms(i)%CDF(j)%A(k), matter%Atoms(i)%CDF(j)%G(k)
               if (.not. read_well) then
                  write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
                  call Save_error_details(Err, 3, Error_descript) ! modlue "Objects"
                  print*, trim(adjustl(Error_descript))
                  INFO = 1   ! error in reading file
                  return ! exit the subroutine
               endif
            enddo

         else DOCDF
            matter%Atoms(i)%TOCS(j) = 0   ! BEB cross-section
            matter%Atoms(i)%TOCSph(j) = 0 ! EPDL cross-section
         endif DOCDF
      enddo
   enddo AT_NUM
end subroutine read_CDF_file


subroutine write_CDF_file(FN, Material, Chem, dens, Efermi, c_sound, Atoms)
   integer, intent(in) :: FN  ! file to write CDF data into
   character(*), intent(in) :: Material   ! material name
   character(*), intent(in) :: Chem       ! chemical formula
   real(8), intent(in) :: dens   ! density of the target [g/cm^3]
   real(8), intent(in) :: Efermi ! [eV] fermi energy
   real(8), intent(in) :: c_sound  ! speed of sound [m/s]
   type(At_data), dimension(:), intent(in) :: Atoms   ! all kinds of atoms of the compound
   !-------------------------
   integer :: N_KOA, i, j, k
   character(200) :: line_to_write

   N_KOA = size(Atoms)  ! kinds of atoms

   write(FN,'(a)') trim(adjustl(Material)) ! material name
   write(FN,'(a)') trim(adjustl(Chem))     ! chemical formula
   write(line_to_write,'(f15.6, f15.6, f15.6, a)') dens, c_sound, Efermi, '   ! density [g/cm^3], speed of sound [m/s], Fermi level [eV]'
   write(FN,'(a)') trim(adjustl(line_to_write))

   AT_NUM:do i = 1, N_KOA ! for each kind of atoms:
      write(line_to_write,'(i4,a)') Atoms(i)%sh, '   ! number of shells in element: '//trim(adjustl(Atoms(i)%Name))
      write(FN,'(a)') trim(adjustl(line_to_write))

      if (allocated(Atoms(i)%N_CDF)) then
         SH_NUM:do j = 1, Atoms(i)%sh  ! for all shells
            write(line_to_write,'(i3, i4, f15.6, f15.6, es15.6E3, a)') Atoms(i)%N_CDF(j), Atoms(i)%Shl_dsgnr(j), &
               Atoms(i)%Ip(j), Atoms(i)%Ne_shell(j), Atoms(i)%Auger(j), &
               '  ! Oscillators, designator:'//trim(adjustl(Atoms(i)%Shell_name(j)))//', Ip [eV], Ne, Auger [fs]'
            write(FN,'(a)') trim(adjustl(line_to_write))

            CDF_NUM:do k = 1, Atoms(i)%N_CDF(j)  ! for all CDF-functions for this shell
               write(line_to_write,'(f15.6, f15.6, f15.6, a)') Atoms(i)%CDF(j)%E0(k), Atoms(i)%CDF(j)%A(k), Atoms(i)%CDF(j)%G(k),&
                                                                  '   ! E0, A, G'
               write(FN,'(a)') trim(adjustl(line_to_write))
            enddo CDF_NUM
         enddo SH_NUM
      else
         write(FN,'(a)') '************************************************'
         write(FN,'(a)') 'The Ritchie-Howie CDF-coefficients are undefined, cannot print them out!'
         write(FN,'(a)') 'Probably, the cdf-file was not specified or specified incorrectly in INPUT.txt'
      endif
   enddo AT_NUM
end subroutine write_CDF_file



subroutine check_atomic_data(matter, numpar, Err, i, cur_shl, shl_tot, file_atomic_parameters, file_EADL_file)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   integer, intent(in) :: i, cur_shl, shl_tot	! No. of atom, current number of shell, total number of shells for this atom
   character(*), intent(in) :: file_atomic_parameters, file_EADL_file
   !===============================================
   character(200) :: Folder_name, File_name, Error_descript
   integer FN, INFO, j, Z, imax, imin, N_shl
   real(8), dimension(:), allocatable :: Nel   ! number of electrons in each shell
   integer, dimension(:), allocatable :: Shl_num ! shell designator
   character(11), dimension(:), allocatable :: Shell_name
   logical :: file_exist, file_opened

   ! Open eadl.all database:
   !Folder_name = 'INPUT_DATA'//trim(adjustl(numpar%path_sep))//'Atomic_parameters'
   Folder_name = trim(adjustl(numpar%input_path))//trim(adjustl(file_atomic_parameters))   !'Atomic_parameters'
   !File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//'eadl.all'
   File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(file_EADL_file))

   !call open_file('readonly', File_name, FN, INFO, Error_descript)
   INFO = 0
   inquire(file=trim(adjustl(File_name)),exist=file_exist)
   if (.not.file_exist) then
      INFO = 1
      Error_descript = 'File '//trim(adjustl(File_name))//' does not exist, the program terminates'
      call Save_error_details(Err, INFO, Error_descript)
      print*, trim(adjustl(Error_descript))
   else
      !open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), status = 'old')
      FN=104
      open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      inquire(file=trim(adjustl(File_name)),opened=file_opened)
      if (.not.file_opened) then
         INFO = 2
         Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
         call Save_error_details(Err, INFO, Error_descript)
         print*, trim(adjustl(Error_descript))
      endif
   endif

   select case (INFO)
   case (0)
      Z =  matter%Atoms(i)%Z ! atomic number

      ! Shell name:
      if (LEN(trim(adjustl(matter%Atoms(i)%Shell_name(cur_shl)))) < 1) then ! take if from EADL-database
         if ( matter%Atoms(i)%Shl_dsgnr(cur_shl) >= 63 ) then
            matter%Atoms(i)%Shell_name(cur_shl) = 'Valence'
         else ! core shell
            call READ_EADL_TYPE_FILE_int(FN, File_name, Z, 912, INFO, error_message=Error_descript, N_shl=N_shl, &
                  Nel=Nel, Shell_name=Shell_name)    ! module "Dealing_with_EADL"
            matter%Atoms(i)%Shell_name(cur_shl) = Shell_name(cur_shl)
         endif
      endif
      ! Number of electrons in this shell:
      if (matter%Atoms(i)%Ne_shell(cur_shl) .LT. 0) then ! take if from EADL-database
         call READ_EADL_TYPE_FILE_int(FN, File_name, Z, 912, INFO, error_message=Error_descript, N_shl=N_shl, &
                  Nel=Nel, Shl_num=Shl_num) ! module "Dealing_with_EADL"
         if (INFO .NE. 0) then
            call Save_error_details(Err, INFO, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 9999
         endif
         imax = 0
         imin = 0
         call select_imin_imax(imin, imax, matter%Atoms(i)%Shl_dsgnr(cur_shl)) ! sum subshells, module "Dealing_with_EADL"
         if (imax .GT. 0) then
            matter%Atoms(i)%Ne_shell(cur_shl) = SUM(Nel, MASK=((Shl_num .GE. imin) .AND. (Shl_num .LE. imax)))
         else
            matter%Atoms(i)%Ne_shell(cur_shl) = Nel(cur_shl)
         endif
      endif
      ! Ionization potential:
      if (matter%Atoms(i)%Ip(cur_shl) .LT. 0.0d0) then ! take if from EADL-database
         ! Read auger-times:
         call READ_EADL_TYPE_FILE_real(FN, File_name, Z, 913, matter%Atoms(i)%Ip, cur_shl=cur_shl, shl_tot=shl_tot, &
               Shl_dsgtr=matter%Atoms(i)%Shl_dsgnr(cur_shl), INFO=INFO, error_message=Error_descript) ! module "Dealing_with_EADL"
         if (INFO .NE. 0) then
            call Save_error_details(Err, INFO, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 9999
         endif
         ! Read kinetic energies:
         call READ_EADL_TYPE_FILE_real(FN, File_name, Z, 914, matter%Atoms(i)%Ek, cur_shl=cur_shl, shl_tot=shl_tot, &
               Shl_dsgtr=matter%Atoms(i)%Shl_dsgnr(cur_shl), INFO=INFO, error_message=Error_descript) ! module "Dealing_with_EADL"
         if (INFO .NE. 0) then
            call Save_error_details(Err, INFO, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 9999
         endif
         !matter%Atoms(i)%TOCS(cur_shl) = 0 ! BEB cross-section
      endif
      ! Auger-decay time:
      if (matter%Atoms(i)%Auger(cur_shl) .LT. 0.0d0) then ! take if from EADL-database
        ! Read auger-times:
        call READ_EADL_TYPE_FILE_real(FN, File_name, Z, 922, matter%Atoms(i)%Auger, cur_shl=cur_shl, shl_tot=shl_tot, &
               Shl_dsgtr=matter%Atoms(i)%Shl_dsgnr(cur_shl), INFO=INFO, error_message=Error_descript) ! module "Dealing_with_EADL"
        if (INFO .NE. 0) then
            call Save_error_details(Err, INFO, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 9999
         endif
      endif
      close (FN)
   case default
      call Save_error_details(Err, INFO, Error_descript)
      print*, trim(adjustl(Error_descript))
   end select
9999 continue
end subroutine check_atomic_data



END MODULE Dealing_with_CDF
