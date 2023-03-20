! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2023 Nikita Medvedev
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
! This module contains subroutines for dealing with EADL and EPDL97 databases
MODULE Dealing_with_EADL
use Universal_Constants   ! let it use universal constants
use Objects   ! since it uses derived types, it must know about them from module 'Objects'
!use Dealing_with_files ! count lines in file etc.
use Little_subroutines, only : Find_monotonous_LE

implicit none
PRIVATE

 character(12) :: m_EADL_file, m_EPDL_file
 ! Outdated databases:
!  parameter(m_EADL_file = 'eadl.all')
!  parameter(m_EPDL_file = 'epdl97.all')
 ! Updated databases:
 parameter(m_EADL_file = 'EADL2017.all')
 parameter(m_EPDL_file = 'EPDL2017.all')

public :: m_EADL_file, m_EPDL_file, READ_EADL_TYPE_FILE_int, READ_EADL_TYPE_FILE_real, select_imin_imax, &
         READ_EPDL_TYPE_FILE_real, next_designator, Read_EPDL_data, define_PQN

 contains

! Reading EADL data for case of integer arrays:
subroutine READ_EADL_TYPE_FILE_int(FN, File_name, Z_needed, I_needed, INFO, error_message, N_shl, Nel, Shell_name, PQN, Shl_num, Ip, Ek, Radiat, Auger, REDO)
   integer, intent (inout) :: FN
   character(100), intent(in) :: File_name
   integer, intent(in) :: Z_needed, I_needed
   integer, intent(inout) :: INFO  ! 0=file read well; 1=no file; 2=couldn't open; 3=error while reading
   character(100), intent(inout) :: error_message
   ! output varuables:
   integer, intent(out) :: N_shl ! how many shells we have in this atom
   real(8), dimension(:), allocatable, intent(inout) :: Nel   ! number of electrons in each shell
   character(11), dimension(:), allocatable, intent(out), optional :: Shell_name   ! names of the shells
   integer, dimension(:), allocatable, intent(out), optional :: PQN   ! principal quantum number
   integer, dimension(:), allocatable, intent(out), optional :: Shl_num ! shell designator
   real(8), dimension(:), allocatable, intent(out), optional :: Ip ! ionization potential
   real(8), dimension(:), allocatable, intent(out), optional :: Ek ! mean kinetic energy of the shell
   real(8), dimension(:), allocatable, intent(out), optional :: Radiat ! radiative decay time of the shell
   real(8), dimension(:), allocatable, intent(out), optional :: Auger  ! Auger decay time of the shell
   logical, intent(in), optional :: REDO ! redo the calculations of N_shl, or not?
   !=====================================================
   integer Z, A, Yi, Yo, Date, C, I, S, EndCheck, counter
   real(8) AW, X1
   real(8) READ1, READ2, READ3, READ4
   integer run_i, Reason, imin, imax
   logical File_opened

   inquire(file=trim(adjustl(File_name)),opened=File_opened,number=FN)
   if (.not. File_opened) then
      FN = 336
      open(unit=FN, FILE = trim(adjustl(File_name)), status = 'old', action = 'read')
      !open(newunit=FN, FILE = trim(adjustl(File_name)), status = 'old', action = 'read')
   endif

   Z = 0
   do while (Z .NE. Z_needed)
      call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
      call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
      if (Z .GE. 100) then !print*, 'The element is not found in the EADL database...'
         write(error_message,*) '(INT) Element #', Z_needed, ' was not found in the database: eadl.all'
         INFO = 4
      endif
      if (Z .GE. 100) exit
   enddo

   if (Z .EQ. Z_needed) then ! this is the element we were looking for
      do while (I .NE. I_needed) ! find the value that we need
         call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
         call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
         if (Z .GT. Z_needed) then !print*, 'INT: The I-value is not found in the EADL database...'
            write(error_message,*) '(INT) The I-value ', Z, Z_needed, I_needed, ' is not found in the EADL database: eadl.all'
            INFO = 4
         endif
         if (Z .GT. Z_needed) exit
      enddo
      do run_i = 1, counter+1
         backspace(FN)
      enddo

      ! for all shells:
      if (present(REDO)) then
         if (REDO) then ! recalculate it: 
            N_shl = counter ! that's how many shells we have in this atom
         endif
         ! otherwise, don't, use the given value only
      else ! do it by defaut:
         N_shl = counter ! that's how many shells we have in this atom
      endif
      if (.not. allocated(Nel)) allocate(Nel(N_shl)) ! allocate number of electrons

      if (present(Shell_name)) then
         if (.not. allocated(Shell_name)) allocate(Shell_name(N_shl)) ! allocate shell names
      endif
      if (present(Shl_num)) then
         if (.not. allocated(Shl_num)) allocate(Shl_num(N_shl)) ! allocate shell disignator for each shell
      endif
      if (present(Ip)) then
         if (.not. allocated(Ip)) allocate(Ip(N_shl)) ! allocate ionization potentials
      endif
      if (present(Ek)) then
         if (.not. allocated(Ek)) allocate(Ek(N_shl)) ! allocate mean kinetic energies of the shells
      endif
      if (present(Radiat)) then
         if (.not. allocated(Radiat)) then
            allocate(Radiat(N_shl)) ! allocate auger-times
            Radiat = 1d-24 ! to be inversed later
         endif
      endif
      if (present(Auger)) then
         if (.not. allocated(Auger)) then
            allocate(Auger(N_shl)) ! allocate radiative times
            Auger = 1d-24 ! to be inversed later
         endif
      endif

      if (present(PQN) .and. present(Shell_name)) then
         if (.not. allocated(PQN)) allocate(PQN(N_shl)) ! allocate principal quantum number
         do run_i = 1, N_shl
            !read(FN,9991) READ1, READ2 ! Shell designator; No. of electrons in the shell
            read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
            if (Reason /= 0) then	! possible old format is used, try it:
               backspace(FN)	! to try the same line
               read(FN,9991) READ1, READ2
            endif 
            Nel(run_i) = READ2 ! number of electrons in this shell
            if (present(Shl_num)) Shl_num(run_i) = READ1 ! shell designator 
            call define_PQN(INT(READ1), Shell_name(run_i), PQN(run_i))
         enddo
      else if (present(PQN)) then
         do run_i = 1, N_shl
            !read(FN,9991) READ1, READ2 ! Shell designator; No. of electrons in the shell
            read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
            if (Reason /= 0) then	! possible old format is used, try it:
               backspace(FN)	! to try the same line
               read(FN,9991) READ1, READ2
            endif 
            Nel(run_i) = READ2 ! number of electrons in this shell
            if (present(Shl_num)) Shl_num(run_i) = READ1 ! shell designator 
            call define_PQN(INT(READ1), PQN=PQN(run_i))
         enddo
      else if (present(Shell_name)) then
         do run_i = 1, N_shl
            !read(FN,9991) READ1, READ2 ! Shell designator; No. of electrons in the shell
            read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
            if (Reason /= 0) then	! possible old format is used, try it:
               backspace(FN)	! to try the same line
               read(FN,9991) READ1, READ2
            endif 
            Nel(run_i) = READ2 ! number of electrons in this shell
            if (present(Shl_num)) Shl_num(run_i) = READ1 ! shell designator 
            call define_PQN(INT(READ1), Shell_name=Shell_name(run_i))
         enddo
      else
         do run_i = 1, N_shl
            !read(FN,9991) READ1, READ2 ! Shell designator; No. of electrons in the shell
            read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
            if (Reason /= 0) then	! possible old format is used, try it:
               backspace(FN)	! to try the same line
               read(FN,9991) READ1, READ2
            endif 
            Nel(run_i) = READ2 ! number of electrons in this shell
            if (present(Shl_num)) Shl_num(run_i) = READ1 ! shell designator
         enddo
      endif
   endif

!    flush(FN)
   rewind(FN)
9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
end subroutine READ_EADL_TYPE_FILE_int


! Reading EADL data for case of real arrays:
subroutine READ_EADL_TYPE_FILE_real(FN, File_name, Z_needed, I_needed, Array1, cur_shl, diff_cur_shl, shl_tot, Shl_dsgtr, INFO, error_message, cont)
!    call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(N_at)%Auger, N_at, cur_shl, shl) ! read auger-times
   integer, intent(inout) :: FN
   character(100), intent(in) :: File_name
   integer, intent(in) :: Z_needed, I_needed
   real(8), dimension(:), allocatable, intent(inout) :: Array1 ! array with variables that we need
   integer, intent(in), optional :: cur_shl, diff_cur_shl, shl_tot, Shl_dsgtr  ! shell number, in which shell to write the data, number of shells, shell designator
   integer, intent(inout) :: INFO  ! 0=file read well; 1=no file; 2=couldn't open; 3=error while reading
   character(100), intent(inout) :: error_message
   logical, intent(in), optional :: cont    ! continue reading file, or start from the brginning?

   integer Z, A, Yi, Yo, Date, C, I, S, EndCheck, counter
   real(8) AW, X1
   real(8) READ1, READ2, READ3, READ4
   integer run_i, imax, imin, Reason, icont
   logical File_opened, Exist_val

   Exist_val = .true.

   inquire(file=trim(adjustl(File_name)),opened=File_opened,number=FN)
   if (.not. File_opened) then
      !open(newunit=FN, FILE = trim(adjustl(File_name)), status = 'old', action = 'read') ! EADL
      FN = 100
      open(unit=FN, FILE = trim(adjustl(File_name)), status = 'old', action = 'read') ! EADL
   endif
   Z = 0 ! to start with
   
   do while (Z .NE. Z_needed)
      call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
      call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
      if (Z .GE. 100) then 
         write(error_message,*) 'Element #', Z_needed, ' was not found in the database: eadl.all'
         INFO = 4
      endif
      if (Z .GE. 100) exit
   enddo

   if (Z .EQ. Z_needed) then ! this is the element we were looking for
      do while (I .NE. I_needed) ! find the value that we need
         call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
         call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
         if (Z .GT. Z_needed) then
            write(error_message,*) 'The I-value ', Z, Z_needed, I_needed, ' is not found in the EADL database: eadl.all'
!             INFO = 4
            Exist_val = .false.
         endif
         if (Z .GT. Z_needed) exit
      enddo
      do run_i = 1, counter+1
         backspace(FN)
      enddo

      if (present(cur_shl)) then ! this shell only:
         if (.not. allocated(Array1)) then
            allocate(Array1(shl_tot)) ! now we know the dimensions...
            Array1 = 0.0d0
         endif
         prsnt:if (present(Shl_dsgtr)) then
          if (Shl_dsgtr .LT. 63) then
            READ1 = 0
            do while (READ1 .NE. Shl_dsgtr)
               !read(FN,9991, IOSTAT=Reason) READ1, READ2
               read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
               if (Reason /= 0) then	! possible old format is used, try it:
                  backspace(FN)	! to try the same line
                  read(FN,9991) READ1, READ2
               endif 
               if (READ1 .EQ. Shl_dsgtr) exit
               IF (Reason .LT. 0) THEN ! end of file reached, nothing found
                  rewind(FN)
                  exit
               endif
            enddo
            if (Reason .LT. 0) then ! If there is no such shell in the database, try to sum up sub-shells:
                Z = 0
                do while (Z .NE. Z_needed)
                   call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
                   call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
                   if (Z .GE. 100) then
                      write(error_message,*) 'The element ', Z_needed, ' is not found in the EADL database: eadl.all'
                      INFO = 4
                   endif
                   if (Z .GE. 100) exit
                enddo
                do while (I .NE. I_needed) ! find the value that we need
                   call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
                   call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
                   if (Z .GT. Z_needed) then
                      write(error_message,*) 'INT: The I-value is not found in the EADL database: eadl.all'
                      INFO = 4
                   endif
                   if (Z .GT. Z_needed) exit
                enddo
                do run_i = 1, counter+1
                   backspace(FN)
                enddo
               READ3 = 0.0d0
               call select_imin_imax(imin, imax, Shl_dsgtr) ! sum subshells:

               READ3 = 0.0d0
               icont = 0
               !do run_i = imin, imax
               do run_i = 1, imax
!                   read(FN,9991,IOSTAT=Reason) READ1, READ2 ! Shell designator, No. of electrons in the shell
                  read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
                  if (Reason /= 0) then	! possible old format is used, try it:
                     backspace(FN)	! to try the same line
                     read(FN,9991) READ1, READ2
                  endif 
!                   print*, READ1, READ2*1d6
                  if (READ1 .GE. imin) then
                     if (isnan(READ2)) then
                     else
                        READ3 = READ3 + READ2 ! number of electrons in this shell
                     endif
                     icont = icont + 1
                  endif
                  if (READ1 .EQ. imax) exit
               enddo
               if (present(diff_cur_shl)) then
                  Array1(diff_cur_shl) = (READ3/real(icont))*1d6
               else
                  Array1(cur_shl) = (READ3/real(icont))*1d6
               endif
            else
               if (present(diff_cur_shl)) then
                  Array1(diff_cur_shl) = READ2*1d6
               else
                  Array1(cur_shl) = READ2*1d6
               endif
            endif
          else
            if (present(diff_cur_shl)) then
               Array1(diff_cur_shl) = 1d22
            else
               Array1(cur_shl) = 1d22
            endif
          endif
         else prsnt
            print*, 'In this case Shl_dsgtr MUST be present!'
            print*, 'in the subroutine  READ_EADL_TYPE_FILE_real'
         endif prsnt
         if (present(diff_cur_shl)) then
            if (isnan(Array1(diff_cur_shl))) Array1(diff_cur_shl) = 0.0d0 ! just in case...
         else
            if (isnan(Array1(cur_shl))) Array1(cur_shl) = 0.0d0 ! just in case...
         endif

      else ! for all shells:
         if (.not. allocated(Array1)) then
            allocate(Array1(counter)) ! now we know the dimensions...
            Array1 = 0.0d0
         endif
         if (Exist_val) then	! if the value is found in the database, use it, if not leave default
            do run_i = 1, size(Array1)
               !read(FN,9991) READ1, READ2
                read(FN,*, IOSTAT=Reason) READ1, READ2 
                if (Reason /= 0) then	! possibly old format is used, try it:
                   backspace(FN)	! to try the same line
                   read(FN,9991) READ1, READ2 !
                endif
               Array1(run_i) = READ2*1d6
            enddo
         endif
      endif
   endif
!    flush(FN)
   rewind(FN)
9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
end subroutine READ_EADL_TYPE_FILE_real


subroutine select_imin_imax(imin, imax, shl_dsgntr)
   integer, intent(in) :: shl_dsgntr
   integer, intent(out) :: imin, imax ! according to EADL database, those are corresponding shells:
   select case (shl_dsgntr) ! sum subshells:
   case (2)
      imin = 3
      imax = 6
   case (4)
      imin = 5
      imax = 6
   case (7)
      imin = 8
      imax = 14
   case (9)
      imin = 10
      imax = 11
   case (12)
      imin = 13
      imax = 14
   case (15)
      imin = 16
      imax = 25
   case (26)
      imin = 27
      imax = 39
   case (40)
      imin = 41
      imax = 56
   case (57)
      imin = 58
      imax = 61
   case default
      imin = 0
      imax = 0
   endselect
end subroutine select_imin_imax


subroutine next_designator(Shl_dsgtr, Last_shl) ! find the designator for the VB (next shell after last one)
   integer, intent(in) :: Last_shl
   integer, intent(out) :: Shl_dsgtr
   select case (Last_shl) ! sum subshells:
   case (:1)
      Shl_dsgtr = 2
   case (2:6)
      Shl_dsgtr = 7
   case (7:14)
      Shl_dsgtr = 15
   case (15:25)
      Shl_dsgtr = 26
   case (26:39)
      Shl_dsgtr = 40
   case (40:)
      Shl_dsgtr = 57
   endselect
end subroutine next_designator


!NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
! New subroutines to deal with EPDL (much faster version):

subroutine Read_EPDL_data(path_sep, FN, File_name, INFO, Z, Shl_dsgnr, Phot_abs_CS_tot, Phot_abs_CS_shl)
   character(1), intent(in) :: path_sep	! path separator
   integer, intent(in) :: FN	! file number of ENDL database
   character(100), intent(in) :: File_name	! file name of ENDL database
   integer, intent(inout) :: INFO	! info whether file read well
   integer, intent(in) :: Z	! atomic number
   integer, dimension(:), intent(in) :: Shl_dsgnr	! array of EADL shell designators
   real(8), dimension(:,:), allocatable, intent(out) :: Phot_abs_CS_tot	! Total photoabsorption cross section:  Photon energy [eV], Cross section [A^2]
   real(8), dimension(:,:,:), allocatable, intent(out) :: Phot_abs_CS_shl	! per shell cross sections of photoabsorption (to be extracted from EPDL): shell, Cross section [A^2]
   !-----------------------------
   integer :: EndCheck, i, counter, j, Nsiz, shl_num_1, Nshl
   integer :: Z0, A0, Yi0, Yo0, Date0, C0, I0, S0, X10, Iflag0
   real(8) :: AW0, temp
   
   Nshl = size(Shl_dsgnr)	! number of shells in this atom
   
   ! Start by reading the header of the first block
   call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10, Iflag0)	! see below
   if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine

   EndCheck = 0	! to start with
   do while (Z0 < Z)	! skip lines until you find the element you need
      do while (EndCheck /= 1)	! until end of the block
         read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
         if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
      enddo
      EndCheck = 0	! restart checker
      ! Read the header of the next block:
      call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10, Iflag0)	! see below
      if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
!       print*, 'EPDL:', Z0
   enddo
   
   ! We found the element we need:
   Z0Z:do while (Z0 == Z)	! read all the data for this element
      ! Read data from the blocks of this elements:
      EndCheck = 0
      do while (EndCheck /= 1)	! read the entire block
         read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
         if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
         if (EndCheck /= 1) then	! if it is not the end of a block, read the line
            backspace(FN)	! get back to read this line again
            select case (C0)	! which process is it?
            case (71)	! Coherent scattering
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (72)	! Incoherent scattering
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (73)	! Photoabsorption
               select case (I0)	! which kind of data is it?
               case (0)	! cross section
                  select case (S0)
                  case (0)	! total cross section
!                      print*, 'Iflag0 tot', Iflag0
                     ! Allocate array if needed:
                     ! First, count how many data points we have here:
                     if (.not. allocated(Phot_abs_CS_tot)) then
                        counter = 0
                        do while (EndCheck /= 1)	! until end of the block
                           counter = counter + 1		! count lines
                           read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
                           if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        enddo
                        ! Then, rewind the data back to read:
                        do j = 1, counter	! rewind back to the start of the block
                           backspace(FN)	! get back to read this line again
                        enddo
                        ! Finaly, allocate arrays to the given number:
                        Nsiz = counter-1
                        allocate(Phot_abs_CS_tot(2,Nsiz))
                        Phot_abs_CS_tot = 0.0d0
                        allocate(Phot_abs_CS_shl(Nshl,2,Nsiz))
                        Phot_abs_CS_shl(:,1,:) = 1d24
                        Phot_abs_CS_shl(:,2,:) = -1.0d-12
                     else	! in case we know the size, for some reason...
                         Nsiz = size(Phot_abs_CS_tot,2)
                     endif
                     
                     EndCheck = 0
                     counter = 0
                     ! Now read the total cross section into this array:
                     do while (EndCheck /= 1)	! until end of the block
                        counter = counter + 1
                        read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
                        if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        if (EndCheck == 1) then
                           exit	! done with this block
                        else
                           backspace(FN)
                           read(FN,*,IOSTAT=INFO) Phot_abs_CS_tot(1,counter), Phot_abs_CS_tot(2,counter)
                           if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        endif
                     enddo ! while (EndCheck /= 1)
                  case (91)	! cross section per shell
                     counter = 0
                     do while (EndCheck /= 1)	! until end of the block
                        counter = counter + 1		! count lines
                        read(FN,9997,IOSTAT=INFO) EndCheck	! check if it is the end of a block
                        if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        if (EndCheck == 1) then
                           exit	! done with this block
                        else
                           backspace(FN)
                           call find_shell_by_designator(X10, Shl_dsgnr, shl_num_1)
                           read(FN,*,IOSTAT=INFO) Phot_abs_CS_shl(shl_num_1,1,counter), Phot_abs_CS_shl(shl_num_1,2,counter)
                           if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
                        endif
                     enddo 
                  case default ! skip line
                     read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
                  endselect
               case (10)	! average energy of emitted particle
                  read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
               case (11)	! average energy of atom
                  read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
               case default ! skip line
                  read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
               end select
            case (74)	! Pair production
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case (75)	! Triplet production
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            case default ! skip line
               read(FN,*,IOSTAT=INFO) ! We don't use this in the code, so just skip such lines
            end select
         endif	! it is the end of the block, go to the next iteration of the cycle
         
      enddo ! while (EndCheck /= 1)
      ! After the end of the block, read the header of the next one:
      call Read_ENDL_header(FN, File_name, INFO, Z0, A0, Yi0, Yo0, AW0, Date0, C0, I0, S0, X10, Iflag0)	! see below
      if (INFO /= 0) goto 9998	! if there was any error, or end of file is reached, exit the subroutine
   enddo Z0Z ! (Z0 == Z)

9998 continue
9997   format(71X,I1)			! block separator
end subroutine Read_EPDL_data



subroutine Read_ENDL_header(FN, File_name, INFO, Z, A, Yi, Yo, AW, Date, C, I, S, X1, Iflag)
   integer, intent(in) :: FN	! file number of ENDL database
   character(100), intent(in) :: File_name	! file name of ENDL database
   integer, intent(inout) :: INFO	! info whether file read well
   integer, intent(out), optional :: Z 	! atomic number
   integer, intent(out), optional :: A	! mass number (in all cases=0, for elemental data)
   integer, intent(out), optional :: Yi	! incident particle designator (see Table II)
   integer, intent(out), optional :: Yo	! outgoing particle designator (see Table II)
   real(8), intent(out), optional :: AW	! atomic mass (amu)
   integer, intent(out), optional :: Date	! date of evaluation (YYMMDD)
   integer, intent(out), optional :: Iflag	! interpolation flag (used in EPDL, not EADL)
   integer, intent(out), optional :: C	! reaction descriptor (see Table II)
   integer, intent(out), optional :: I		! reaction property (see Table II)
   integer, intent(out), optional :: S	! reaction modifier (see Table II)
   integer, intent(out), optional :: X1	! subshell designator (see Table VI)
   !----------------------------------------------
   integer :: Z0, A0, Yi0, Yo0, Date0, C0, I0, S0, Iflag0
   real(8) :: AW0, X10
!    character(72) :: test
   ! First header line:
   
!    read(FN,'(a)') test
!    print*, ':::: ', test
!    backspace(FN)
   if (present(Iflag)) then ! EPDL database:
      read(FN,'(i3, i3, 1X, i2, 1X, i2, 1X, E11.4, 1x, i6, i1)',IOSTAT=INFO) Z0, A0, Yi0, Yo0, AW0, Date0, Iflag0
   else ! EADL database:
      read(FN,'(i3, i3, 1X, i2, 1X, i2, 1X, E11.4, 1x, i6)',IOSTAT=INFO) Z0, A0, Yi0, Yo0, AW0, Date0
   endif
!    print*, '1) Read_ENDL_header', INFO
!    print*, Z0, A0, Yi0, Yo0, AW0, Date0
   if (INFO /= 0) then	! something went wrong while reading this line
      goto 9998	! exit the subroutine
   endif
   
   ! Second header line:
   read(FN,'(i2, i3, i3, 14X, E11.4)',IOSTAT=INFO) C0, I0, S0, X10
!    print*, '2) Read_ENDL_header', INFO
!    print*, C0, I0, S0, X10
   if (INFO /= 0) then	! something went wrong while reading this line
      goto 9998	! exit the subroutine
   endif
   
   ! Save all required output:
   if (present(Z)) Z = Z0
   if (present(A)) A = A0
   if (present(Yi)) Yi = Yi0
   if (present(Yo)) Yo = Yo0
   if (present(AW)) AW = AW0
   if (present(Date)) Date = Date0
   if (present(Iflag)) Iflag = Iflag0
   if (present(C)) C = C0
   if (present(I)) I = I0
   if (present(S)) S = S0
   if (present(X1)) X1 = NINT(X10)
   
9998 continue
end subroutine Read_ENDL_header



subroutine find_shell_by_designator(shl_dsgntr, sch_array, shl_num)
   integer, intent(in) :: shl_dsgntr	! shell designator
   integer, dimension(:), intent(in) :: sch_array	! array of all designators
   integer, intent(out) :: shl_num	! shell number in arrays
   call Find_monotonous_LE(dble(sch_array), dble(shl_dsgntr), shl_num)	! module "Little_subroutines"
end subroutine find_shell_by_designator



!OOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! Old subroutines to deal with EPDL:

! Reading data from EPDL file: photoabsrotbtion cross-section for given photon energy
subroutine READ_EPDL_TYPE_FILE_real(FN2, File_name, Z_needed, C_needed, I_needed, S_needed, Photoabs_sigma, Eph, Shl_design, Last_dsgtr, Photoabs_sigma_tot)
   integer, intent (inout) :: FN2
   character(100), intent(in) :: File_name
   integer, intent(in) :: Z_needed, C_needed, I_needed, S_needed
   real(8), intent(in) :: Eph ! photon energy [eV]
   !type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   !integer, intent(in) :: Nat, Shl ! number of atom and a shell
   integer, intent(in) :: Shl_design ! shell designator - for this shell we are looking the cross section
   integer, intent(in), optional :: Last_dsgtr ! last shell which wasn't VB
   real(8), intent(inout) :: Photoabs_sigma ! Cross section for this shell [A^2]
   real(8), intent(inout), optional ::  Photoabs_sigma_tot  ! Total cross section [A^2]
   integer Z, A, Yi, Yo, Date, Iflag, C, I, S, EndCheck, counter, run_i
   real(8) AW, X1
   integer imin, imax, Shl_dsgtr
   real(8) READ1, READ2, READ3, READ4, E_photon
   logical File_opened

   inquire(file=trim(adjustl(File_name)),opened=File_opened)
   if (.not. File_opened) then
      !FN2 = 338
      open(FN2, FILE = trim(adjustl(File_name)), status = 'old', action = 'read') ! EPDL97
   endif

   E_photon = Eph/1d6 ! [MeV] units in the EPDL-database
   Photoabs_sigma = 0.0d0
   Shl_dsgtr = Shl_design
   if (present(Last_dsgtr)) then
       if (Shl_design .GE. 63) then
          call next_designator(Shl_dsgtr, Last_dsgtr) ! find the designator for the VB (next shell after last one)
       else
          Shl_dsgtr = Shl_design
       endif
   endif
   call select_imin_imax(imin, imax, Shl_dsgtr) ! find subshell designators

   Z = 0
   do while (Z .NE. Z_needed)
      !print*, 'Z:', Z, Z_needed
      call read_throu_EPDL(FN2, Z, A, Yi, Yo, AW, Date, Iflag, C, I, S, X1, EndCheck)
      call How_many_lines(FN2, counter, EndCheck) ! to know how many lines are there
      if (Z .GE. 100) print*, 'The element is not found in the EPDL97 database...'
      if (Z .GE. 100) exit
   enddo

   if (Z .EQ. Z_needed) then ! this is the element we were looking for
     run_i = 0.0e0
     do while (Z .LT. Z_needed + 1) ! do throughout the whole data for this element: all shells
        call read_throu_EPDL(FN2, Z, A, Yi, Yo, AW, Date, Iflag, C, I, S, X1, EndCheck)
        !write(*,'(a,i,i,i,i)') 'test 1.5', Z, C, I, S
        if ((C .EQ. C_needed) .AND. (I .EQ. I_needed) .AND. (S .EQ. S_needed)) then ! read it:
           SELECT CASE (S_needed) ! which data are these:
              CASE (0) ! total cross-section
                 if (present(Photoabs_sigma_tot)) then
                    call Find_value_while_reading(FN2, Iflag, counter, EndCheck, E_photon, READ2)
                    Photoabs_sigma_tot = READ2*1d-8 ! [A^2]
                 endif
              CASE (91) ! cross-section by shell
                 call Find_value_while_reading(FN2, Iflag, counter, EndCheck, E_photon, READ2)
                 run_i = run_i + 1
                 !write(*,'(a,f7.2,i,i,i)') 'run_i', X1, Shl_dsgtr, imin, imax
                 if ((X1 .EQ. Shl_dsgtr) .OR. ( (X1 .GE. imin) .AND. (X1 .LE. imax) ) )then
                     !print*, 'summed', Photoabs_sigma, READ2
                     Photoabs_sigma = Photoabs_sigma + READ2*1d-8 ! [A^2]
                 endif
              CASE DEFAULT
                 call How_many_lines(FN2, counter, EndCheck) ! to know how many lines are there
                 print*, 'Who knows what is going on...'
           END SELECT
        else ! just skip all these lines:
           call How_many_lines(FN2, counter, EndCheck) ! to know how many lines are there
        endif ! ((C .EQ. C_needed) .AND. (I .EQ. I_needed) .AND. (S .EQ. S_needed))
     enddo !while (Z .LT. Z_needed + 1) ! do throughout the whole data for this element: all shells
   endif ! (Z .EQ. Z_needed)
!    flush(FN2)
    rewind(FN2)

9999   format(I3,I3,I2,I2,E11.4,I6,I1)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
end subroutine READ_EPDL_TYPE_FILE_real


subroutine Check_EADP_EPDL97(path_sep, File_name, File_name2, INFO)
    character(1), intent(in) :: path_sep
    character(100), intent(inout) :: File_name
    character(100), intent(inout) :: File_name2
    integer INFO
    character(100) :: Folder_name
    logical file_exist
    Folder_name = 'INPUT_EADL'  ! here we keep databases
    File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'eadl.all'
    File_name2 = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'epdl97.all'
    inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
    if (.not. file_exist) then
       write(*,'(a,a,a)') 'File ',trim(adjustl(File_name)),' is not found.'
       INFO = 1
       !print*, 'The program cannot continue.'
       !goto 9999
    else
       INFO = 0
    endif
    inquire(file=trim(adjustl(File_name2)),exist=file_exist) ! check if input file is there
    if (.not. file_exist) then
       write(*,'(a,a,a)') 'File ',trim(adjustl(File_name2)),' is not found.'
       !print*, 'The program cannot continue.'
       !goto 9999
       INFO = 2
    endif
end subroutine Check_EADP_EPDL97




! Reading first 2 lines with the descriptions in EPDL database: 
subroutine read_throu_EPDL(FN, Z, A, Yi, Yo, AW, Date1, Iflag, C, I, S, X1, EndCheck)
   integer, intent (inout) :: FN, Z, A, Yi, Yo, Date1, Iflag, C, I, S, EndCheck
   real(8), intent (inout) :: AW, X1
   real(8) skipread
   integer :: INFO	! info whether file read well
   EndCheck = 0
   ! First header line:
   !read(FN,'(i3, i3,i2, i2, E11.4, 2x, i8, i1)') Z, A, Yi, Yo, AW, Date1, Iflag  	! 1st line
   read(FN,'(i3, i3, 1X, i2, 1X, i2, 1X, E11.4, 1x, i6, i1)',IOSTAT=INFO) Z, A, Yi, Yo, AW, Date1, Iflag
   ! Second header line:
   !read(FN,'(i2,i3,i3,E11.4,E11.4)') C, I, S, skipread, X1  			! 2st line
   read(FN,'(i2, i3, i3, 14X, E11.4)',IOSTAT=INFO) C, I, S, X1
   read(FN,9997) EndCheck  			! 3st line, cheking...
   
9009   format(i3,i3,i2,i2,e11.4,i6)	! 1st line
9999   format(I3,I3,I2,I2,E11.4,I6,I1)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
end subroutine read_throu_EPDL


! Reading first 2 lines with the descriptions in EADL database:
subroutine read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
   integer, intent (inout) :: FN, Z, A, Yi, Yo, Date, C, I, S, EndCheck
   real(8), intent (inout) :: AW, X1
   EndCheck = 0
   read(FN,9999) Z, A, Yi, Yo, AW, Date  	! 1st line
   read(FN,9998) C, I, S, X1  			! 2st line
   read(FN,9997) EndCheck  			! 3st line, cheking...

9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
end subroutine read_throu 


! ! Reading first 2 lines with the descriptions in EPDL database: 
! subroutine read_throu_EPDL(FN, Z, A, Yi, Yo, AW, Date1, Iflag, C, I, S, X1, EndCheck)
!    integer, intent (inout) :: FN, Z, A, Yi, Yo, Date1, Iflag, C, I, S, EndCheck
!    real(8), intent (inout) :: AW, X1
!    real(8) skipread
!    EndCheck = 0
!    read(FN,'(i3, i3,i2, i2, E11.4, 2x, i8, i1)') Z, A, Yi, Yo, AW, Date1, Iflag  	! 1st line
!    !read(FN,'(i2,i3,i3,E11.4,E11.4)') C, I, S, skipread, X1  			! 2st line
!    read(FN,'(i2,i3,i3,13X,E11.4)') C, I, S, X1  			! 2st line
!    read(FN,9997) EndCheck  			! 3st line, cheking...
! 9009   format(i3,i3,i2,i2,e11.4,i6)	! 1st line
! 9999   format(I3,I3,I2,I2,E11.4,I6,I1)	! 1st line
! 9998   format(I2,I3,I3,E11.4)		! 2st line
! 9997   format(71X,I1)			! last line
! end subroutine read_throu_EPDL
! 
! ! Reading first 2 lines with the descriptions in EADL database:
! subroutine read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
!    integer, intent (inout) :: FN, Z, A, Yi, Yo, Date, C, I, S, EndCheck
!    real(8), intent (inout) :: AW, X1
!    EndCheck = 0
!    read(FN,9999) Z, A, Yi, Yo, AW, Date  	! 1st line
!    read(FN,9998) C, I, S, X1  			! 2st line
!    read(FN,9997) EndCheck  			! 3st line, cheking...
! 9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
! 9998   format(I2,I3,I3,E11.4)		! 2st line
! 9997   format(71X,I1)			! last line
! end subroutine read_throu 


! Counts how many lines contain the data in a block of EADL, EPDL-databases
subroutine How_many_lines(FN,counter, EndCheck)
   integer, intent(in) :: FN
   integer, intent(inout) :: counter, EndCheck
   counter = 0
   do while (EndCheck .NE. 1) ! until end of data for this sub-shell/input
      counter = counter + 1 ! how many lines are there
      read(FN,9997) EndCheck  			! 3st line, cheking...
      !print*, 'How_many_lines ', counter, EndCheck
      if (counter .GE. 5000) EndCheck = -1	! indicates some ERROR
      if (counter .GE. 5000) exit	! if something went wrong...
   enddo
9997   format(71X,I1)			! last line
end subroutine How_many_lines


! EADL database, that's how to find the cross-section for given photon energy for each subshell:
subroutine Find_value_while_reading(FN, Iflag, counter, EndCheck, E_photon, OUT_value)
   integer, intent(in) :: FN, Iflag
   integer, intent(inout) :: counter, EndCheck
   real(8), intent(in) :: E_photon
   real(8), intent(out) :: OUT_value
   real(8) READ1, READ2, E1, E2, Sigma1, Sigma2
   integer :: Reason
   LOGICAL :: Firsttime
   Firsttime = .true.
   counter = 0
   do while (EndCheck .NE. 1) ! until end of data for this sub-shell/input
      counter = counter + 1 ! how many lines are there
      read(FN,9997) EndCheck  			! 3st line, cheking...
      if (EndCheck .NE. 1) then
         backspace(FN)
!          read(FN,9991) READ1, READ2
         read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
         if (Reason /= 0) then	! possible old format is used, try it:
            backspace(FN)	! to try the same line
            read(FN,9991) READ1, READ2
         endif 
      endif
      if ((READ1 .GT. E_photon) .AND. (Firsttime)) then
         Firsttime = .false.
         backspace(FN)
         backspace(FN)
         !read(FN,9991) READ1, READ2
         read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
         if (Reason /= 0) then	! possible old format is used, try it:
            backspace(FN)	! to try the same line
            read(FN,9991) READ1, READ2
         endif 
         if (counter .EQ. 1) then ! too small energy!
            OUT_value = 0.0e0
         else ! photon energy is above the ionization potential
            E1 = READ1
            Sigma1 = READ2
            !read(FN,9991) READ1, READ2
            read(FN,*, IOSTAT=Reason) READ1, READ2  ! Shell designator; No. of electrons in the shell
            if (Reason /= 0) then	! possible old format is used, try it:
               backspace(FN)	! to try the same line
               read(FN,9991) READ1, READ2
            endif 
            E2 = READ1
            Sigma2 = READ2
            call Interpolate_EPDL(Iflag, E1, E2, Sigma1, Sigma2, E_photon, OUT_value)
         endif
         !OUT_value = READ2
      endif
      if (counter .GE. 5000) EndCheck = -1	! indicates some ERROR
      if (counter .GE. 5000) exit	! if something went wrong...
   enddo
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
9997   format(71X,I1)			! last line
end subroutine Find_value_while_reading


! Interpolation of photoabsorbtion cross-section according to EADL database:
subroutine Interpolate_EPDL(Iflag, E1, E2, Sigma1, Sigma2, E_needed, OUT_value)
   integer, intent(in) :: Iflag
   real(8), intent(in) :: E1, E2, Sigma1, Sigma2, E_needed
   real(8), intent(out) :: OUT_value
   real(8) E2log, E1log, E_needed_log, Sigma1log, Sigma2log
   select case(Iflag) ! what interpolation to use:
      case(0,2) ! linear x and y
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1)
      case(3)	! logarithmic x, linear y
         E2log = log(E2)
         E1log = log(E1)
         E_needed_log = log(E_needed)
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2log - E1log)*(E_needed_log - E1log)
      case(4)	! linear x, logarithmic y
         Sigma1log = log(Sigma1)
         Sigma2log = log(Sigma2)
         OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2 - E1)*(E_needed - E1)
         OUT_value = exp(OUT_value)
      case(5)	! logarithmic x and y
         E2log = log(E2)
         E1log = log(E1)
         E_needed_log = log(E_needed)
         Sigma1log = log(Sigma1)
         Sigma2log = log(Sigma2)
         OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2log - E1log)*(E_needed_log - E1log)
         OUT_value = exp(OUT_value)
      case default ! linear x and y
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1) 
   end select
end subroutine Interpolate_EPDL


! According to EADL-data format, the atomic shells are numerated as:
subroutine define_PQN(READ1, Shell_name, PQN)
   integer, intent(in) :: READ1 ! shell designator
   integer, INTENT(inout), optional :: PQN !Principal quantum number
   character(11), intent(inout), optional :: Shell_name ! names
   SELECT CASE(INT(READ1))
   CASE ( : 1)
      if (present(PQN)) PQN = 1 ! K-shell
      if (present(Shell_name)) Shell_name = 'K-shell'
   CASE (2)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L-shell'
   CASE (3)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L1-shell'
   CASE (4)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L23-shell'
   CASE (5)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L2-shell'
   CASE (6)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L3-shell'
   CASE (7)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M-shell'
   CASE (8)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M1-shell'
   CASE (9)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M23-shell'
   CASE (10)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M2-shell'
   CASE (11)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M3-shell'
   CASE (12)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M45-shell'
   CASE (13)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M4-shell'
   CASE (14)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M5-shell'
   CASE (15)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N-shell'
   CASE (16)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N1-shell'
   CASE (17)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N23-shell'
   CASE (18)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N2-shell'
   CASE (19)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N3-shell'
   CASE (20)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N45-shell'
   CASE (21)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N4-shell'
   CASE (22)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N5-shell'
   CASE (23)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N67-shell'
   CASE (24)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N6-shell'
   CASE (25)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N7-shell'
   CASE (26)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O-shell'
   CASE (27)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O1-shell'
   CASE (28)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O23-shell'
   CASE (29)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O2-shell'
   CASE (30)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O3-shell'
   CASE (31)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O45-shell'
   CASE (32)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O4-shell'
   CASE (33)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O5-shell'
   CASE (34)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O67-shell'
   CASE (35)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O6-shell'
   CASE (36)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O7-shell'
   CASE (37)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O89-shell'
   CASE (38)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O8-shell'
   CASE (39)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O9-shell'
   CASE (40)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P-shell'
   CASE (41)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P1-shell'
   CASE (42)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P23-shell'
   CASE (43)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P2-shell'
   CASE (44)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P3-shell'
   CASE (45)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P45-shell'
   CASE (46)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P4-shell'
   CASE (47)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P5-shell'
   CASE (48)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P67-shell'
   CASE (49)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P6-shell'
   CASE (50)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P7-shell'
   CASE (51)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P89-shell'
   CASE (52)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P8-shell'
   CASE (53)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P9-shell'
   CASE (54)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P1011-shell'
   CASE (55)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P10-shell'
   CASE (56)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P11-shell'
   CASE (57)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q-shell'
   CASE (58)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q1-shell'
   CASE (59)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q23-shell'
   CASE (60)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q2-shell'
   CASE (61)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q3-shell'
   CASE (62)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q...-shell'
   CASE (63:)
      if (present(PQN)) PQN = 8 ! R-shell
      if (present(Shell_name)) Shell_name = 'Valence'
   CASE DEFAULT
      if (present(PQN)) PQN = 1 ! K-shell, if something else happened...
      if (present(Shell_name)) Shell_name = 'K-shell'
   END SELECT
end subroutine define_PQN


subroutine Find_element_name(Z, Name, Full_Name, M)
   integer, intent (in) :: Z	! atomic number
   real(8), intent(out) :: M	! mass [proton mass]
   Character(3), intent(out) :: Name ! short name of element
   Character(30), intent(out) :: Full_Name ! full name of element

   select case (Z)
   case(1)
 	Name = 'H'
 	Full_Name = 'Hydrogen' 	
 	M = 1.0082e0
   case(2) 	
	Name = 'He' 	
	Full_Name = 'Helium' 	
	M = 4.002602e0
   case(3)	
	Name = 'Li' 	
	Full_Name = 'Lithium' 	
	M = 6.942e0
   case(4)	
	Name = 'Be' 	
	Full_Name = 'Beryllium' 
	M = 9.012182e0
   case(5)
 	Name = 'B' 	
	Full_Name = 'Boron' 
	M = 10.812e0
   case(6)	
 	Name = 'C' 	
	Full_Name = 'Carbon' 
	M = 12.0112e0
   case(7)
	Name = 'N' 	
	Full_Name = 'Nitrogen'  	
 	M = 14.0072e0
   case(8)
	Name = 'O' 	
	Full_Name = 'Oxygen'  	
 	M = 15.9992e0
   case(9)
	Name = 'F' 	
	Full_Name = 'Fluorine'  	
 	M = 18.9984032e0
   case(10)
	Name = 'Ne' 	
	Full_Name = 'Neon'  	
 	M = 20.1797e0
   case(11)
	Name = 'Na' 	
	Full_Name = 'Sodium'  	
 	M = 22.98976928e0
   case(12)
	Name = 'Mg' 	
	Full_Name = 'Magnesium'  	
 	M = 24.3059e0
   case(13)
	Name = 'Al' 	
	Full_Name = 'Aluminium' 	
 	M = 26.9815386e0
   case(14)
	Name = 'Si' 	
	Full_Name = 'Silicon'  	
 	M = 28.0854e0
   case(15)
	Name = 'P' 	
	Full_Name = 'Phosphorus'  	
 	M = 30.973762e0
   case(16)
	Name = 'S' 	
	Full_Name = 'Sulfur'  	
 	M = 32.062e0
   case(17)
	Name = 'Cl' 	
	Full_Name = 'Chlorine' 	
 	M = 35.452e0
   case(18)
	Name = 'Ar' 	
	Full_Name = 'Argon'
	M = 39.948e0
   case(19)
	Name = 'K' 	
	Full_Name = 'Potassium'  	
 	M = 39.0983e0
   case(20)
	Name = 'Ca' 	
	Full_Name = 'Calcium'  	
 	M = 40.078e0
   case(21)
	Name = 'Sc' 	
	Full_Name = 'Scandium'  	
 	M = 44.955912e0
   case(22)
	Name = 'Ti' 	
	Full_Name = 'Titanium'  	
 	M = 47.867e0
   case(23)
	Name = 'V'	
	Full_Name = 'Vanadium'  	
 	M = 50.9415e0
   case(24)
	Name = 'Cr' 	
	Full_Name = 'Chromium'  	
 	M = 51.9961e0
   case(25)
	Name = 'Mn' 	
	Full_Name = 'Manganese' 	
 	M = 54.938045e0
   case(26)
	Name = 'Fe' 	
	Full_Name = 'Iron'  	
 	M = 55.845e0
   case(27)
	Name = 'Co' 	
	Full_Name = 'Cobalt' 	
 	M = 58.933195e0
   case(28)
	Name = 'Ni' 	
	Full_Name = 'Nickel'  	
 	M = 58.6934e0
   case(29)
	Name = 'Cu' 	
	Full_Name = 'Copper'  	
 	M = 63.546e0 
   case(30)
	Name = 'Zn' 	
	Full_Name = 'Zinc'  	
 	M = 65.38e0
   case(31)
	Name = 'Ga' 	
	Full_Name = 'Gallium'  	
 	M = 69.723e0
   case(32)
	Name = 'Ge' 	
	Full_Name = 'Germanium'  	
 	M = 72.630e0
   case(33)
	Name = 'As' 	
	Full_Name = 'Arsenic' 	
 	M = 74.92160e0
   case(34)
	Name = 'Se' 	
	Full_Name = 'Selenium' 	
 	M = 78.96e0
   case(35)
	Name = 'Br' 	
	Full_Name = 'Bromine'  	
 	M = 79.9049e0
   case(36)
	Name = 'Kr' 	
	Full_Name = 'Krypton'  	
 	M = 83.798e0
   case(37)
	Name = 'Rb' 	
	Full_Name = 'Rubidium'  	
 	M = 85.4678e0
   case(38)
	Name = 'Sr' 	
	Full_Name = 'Strontium'  	
 	M = 87.62e0
   case(39)
	Name = 'Y' 	
	Full_Name = 'Yttrium'  	
 	M = 88.90585e0
   case(40)
	Name = 'Zr' 	
	Full_Name = 'Zirconium'  	
 	M = 91.224e0
   case(41)
	Name = 'Nb' 	
	Full_Name = 'Niobium'  	
 	M = 92.90638e0
   case(42)
	Name = 'Mo' 	
	Full_Name = 'Molybdenum'  	
 	M = 95.96e0
   case(43)
	Name = 'Tc' 	
	Full_Name = 'Technetium' 
 	M = 98e0
   case(44)
	Name = 'Ru' 	
	Full_Name = 'Ruthenium'  	
 	M = 101.07e0
   case(45)
	Name = 'Rh' 	
	Full_Name = 'Rhodium'  	
 	M = 102.90550e0
   case(46)
	Name = 'Pd' 	
	Full_Name = 'Palladium'  	
 	M = 106.42e0
   case(47)
	Name = 'Ag' 	
	Full_Name = 'Silver'  	
 	M = 107.8682e0
   case(48)
	Name = 'Cd' 	
	Full_Name = 'Cadmium'  	
 	M = 112.411e0
   case(49)
	Name = 'In' 	
	Full_Name = 'Indium'  	
 	M = 114.818e0
   case(50)
	Name = 'Sn' 	
	Full_Name = 'Tin'
 	M = 118.710e0
   case(51)
	Name = 'Sb' 	
	Full_Name = 'Antimony'  	
 	M = 121.760e0
   case(52)
	Name = 'Te' 	
	Full_Name = 'Tellurium'  	
 	M = 127.60e0
   case(53)
	Name = 'I'	
	Full_Name = 'Iodine'  	
 	M = 126.90447e0
   case(54)
	Name = 'Xe' 	
	Full_Name = 'Xenon'  	
 	M = 131.293e0
   case(55)
	Name = 'Cs' 	
	Full_Name = 'Caesium'  	
 	M = 132.9054519e0
   case(56)
	Name = 'Ba' 	
	Full_Name = 'Barium'  	
 	M = 137.327e0
   case(57)
	Name = 'La' 	
	Full_Name = 'Lanthanum'  	
 	M = 138.90547e0
   case(58)
	Name = 'Ce' 	
	Full_Name = 'Cerium'  	
 	M = 140.116e0
   case(59)
	Name = 'Pr' 	
	Full_Name = 'Praseodymium'  	
 	M = 140.90765e0
   case(60)
	Name = 'Nd' 	
	Full_Name = 'Neodymium'  	
 	M = 144.242e0
   case(61)
	Name = 'Pm' 	
	Full_Name = 'Promethium' 
 	M = 145e0
   case(62)
	Name = 'Sm' 	
	Full_Name = 'Samarium'  	
 	M = 150.36e0
   case(63)
	Name = 'Eu' 	
	Full_Name = 'Europium' 	
 	M = 151.964e0
   case(64)
	Name = 'Gd' 	
	Full_Name = 'Gadolinium'  	
 	M = 157.25e0
   case(65)
	Name = 'Tb' 	
	Full_Name = 'Terbium'  	
 	M = 158.92535e0
   case(66)
	Name = 'Dy' 	
	Full_Name = 'Dysprosium'  	
 	M = 162.500e0
   case(67)
	Name = 'Ho' 	
	Full_Name = 'Holmium'  	
 	M = 164.93032e0
   case(68)
	Name = 'Er' 	
	Full_Name = 'Erbium'  	
 	M = 167.259e0
   case(69)
	Name = 'Tm' 	
	Full_Name = 'Thulium'  	
 	M = 168.93421e0
   case(70)
	Name = 'Yb' 	
	Full_Name = 'Ytterbium'  	
 	M = 173.054e0
   case(71)
	Name = 'Lu' 	
	Full_Name = 'Lutetium'  	
 	M = 174.9668e0
   case(72)
	Name = 'Hf' 	
	Full_Name = 'Hafnium'  	
 	M = 178.49e0
   case(73)
	Name = 'Ta' 	
	Full_Name = 'Tantalum'  	
 	M = 180.94788e0
   case(74)
	Name = 'W' 	
	Full_Name = 'Tungsten'  	
 	M = 183.84e0 
   case(75)
	Name = 'Re' 	
	Full_Name = 'Rhenium'  	
 	M = 186.207e0 
   case(76)
	Name = 'Os' 	
	Full_Name = 'Osmium' 	
 	M = 190.23e0
   case(77)
	Name = 'Ir' 	
	Full_Name = 'Iridium'  	
 	M = 192.217e0
   case(78)
	Name = 'Pt' 	
	Full_Name = 'Platinum'  	
 	M = 195.084e0
   case(79)
	Name = 'Au' 	
	Full_Name = 'Gold'  	
 	M = 196.966569e0
   case(80)
	Name = 'Hg' 	
	Full_Name = 'Mercury'  	
 	M = 200.592e0 
   case(81)
	Name = 'Tl' 	
	Full_Name = 'Thallium' 	
 	M = 204.389e0
   case(82)
	Name = 'Pb' 	
	Full_Name = 'Lead' 
	M = 207.2e0
   case(83)
	Name = 'Bi' 	
	Full_Name = 'Bismuth'  	
 	M = 208.98040e0
   case(84)
	Name = 'Po' 	
	Full_Name = 'Polonium' 
 	M = 209.0e0
   case(85)
	Name = 'At' 	
	Full_Name = 'Astatine' 
 	M = 210.0e0
   case(86)
	Name = 'Rn' 	
	Full_Name = 'Radon' 	
 	M = 222.0e0
   case(87)
	Name = 'Fr' 	
	Full_Name = 'Francium' 	
 	M = 223.0e0
   case(88)
	Name = 'Ra' 	
	Full_Name = 'Radium' 
 	M = 226.0e0
   case(89)
	Name = 'Ac' 	
	Full_Name = 'Actinium' 
 	M = 227.0e0
   case(90)
	Name = 'Th' 	
	Full_Name = 'Thorium'  	
 	M = 232.03806e0
   case(91)
	Name = 'Pa' 	
	Full_Name = 'Protactinium'
 	M = 231.035880e0
   case(92)
	Name = 'U' 	
	Full_Name = 'Uranium'  	
 	M = 238.02891e0
   case(93)
	Name = 'Np' 	
	Full_Name = 'Neptunium'
 	M = 237.0e0
   case(94)
	Name = 'Pu' 	
	Full_Name = 'Plutonium'
 	M = 244.0e0
   case(95)
	Name = 'Am' 	
	Full_Name = 'Americium' 
	M = 243.0e0
   case(96)
	Name = 'Cm' 	
	Full_Name = 'Curium' 
 	M = 247.0e0
   case(97)
	Name = 'Bk' 	
	Full_Name = 'Berkelium' 
 	M = 247.0e0
   case(98)
	Name = 'Cf' 	
	Full_Name = 'Californium' 
 	M = 251.0e0
   case(99)
	Name = 'Es' 	
	Full_Name = 'Einsteinium' 
 	M = 252.0e0
   case(100)
	Name = 'Fm' 	
	Full_Name = 'Fermium'
 	M = 257.0e0
   case(101)
	Name = 'Md' 	
	Full_Name = 'Mendelevium' 
 	M = 258.0e0
   case(102)
	Name = 'No' 	
	Full_Name = 'Nobelium'
 	M = 259.0e0
   case(103)
	Name = 'Lr' 	
	Full_Name = 'Lawrencium'
 	M = 262.0e0
   case(104)
	Name = 'Rf' 	
	Full_Name = 'Rutherfordium' 
 	M = 267.0e0
   case(105)
	Name = 'Db' 	
	Full_Name = 'Dubnium' 
 	M = 268.0e0
   case(106)
	Name = 'Sg' 	
	Full_Name = 'Seaborgium' 
 	M = 269.0e0
   case(107)
	Name = 'Bh' 	
	Full_Name = 'Bohrium'
 	M = 270.0e0
   case(108)
	Name = 'Hs' 	
	Full_Name = 'Hassium' 
 	M = 269.0e0
   case(109)
	Name = 'Mt' 	
	Full_Name = 'Meitnerium' 
 	M = 278.0e0
   case(110)
	Name = 'Ds' 	
	Full_Name = 'Darmstadtium' 	
 	M = 281.0e0
   case(111)
	Name = 'Rg' 	
	Full_Name = 'Roentgenium'
 	M = 281.0e0
   case(112)
	Name = 'Cn' 	
	Full_Name = 'Copernicium'
 	M = 285.0e0
   case(113)
	Name = 'Uut' 	
	Full_Name = 'Ununtrium'
 	M = 286.0e0
   case(114)
	Name = 'Fl' 	
	Full_Name = 'Flerovium' 
 	M = 289.0e0
   case(115)
	Name = 'Uup' 	
	Full_Name = 'Ununpentium' 
 	M = 288.0e0
   case(116)
	Name = 'Lv' 	
	Full_Name = 'Livermorium'
 	M = 293.0e0 
   case(117)
	Name = 'Uus' 	
	Full_Name = 'Ununseptium'	
 	M = 294.0e0
   case(118)
	Name = 'Uuo' 	
	Full_Name = 'Ununoctium'
 	M = 294.0e0
   case(119:)
	Name = 'UFO'
	Full_Name = 'Unknown'
	print*, 'Mass of this element is not in our database.'
	print*, 'Please, provide the value in the units of the proton mass:'
	read*, M
end select
end subroutine Find_element_name


END MODULE Dealing_with_EADL
