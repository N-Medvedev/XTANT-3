PROGRAM Entropy
! Compilation:
!
! for DEBUG:
! ifort.exe -c /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_entropy.f90 /link /stack:9999999999
!
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs *.obj -o XTANT_entropy.exe /link /stack:9999999999

! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_entropy.f90 -o XTANT_entropy.exe /link /stack:9999999999

!
!
! for RELEASE:
! ifort.exe -c /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays XTANT_entropy.f90 /link /stack:9999999999
!
! ifort.exe /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays *.obj -o XTANT_entropy.exe /link /stack:9999999999
!
! To execute:
! XTANT_entropy.exe
!<===========================================
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

real(8), dimension(:), allocatable :: Ei
real(8), dimension(:), allocatable :: distr, flnf
real(8), dimension(:), allocatable :: distr_eq
real(8) :: S, S_eq, tim, g_kb
integer :: FN_distr, FN_out, FN_out1
integer :: INFO, Reason, i, Te, Nsiz
logical :: read_well, equilibr
character(100) :: File_distr, File_out, File_fig, Gnu_file, Gnu_script
character(1) :: path_sep, temp_ch

parameter (g_kb = 11604.0d0)  ! Boltzmann constant [K/eV]

call Path_separator(path_sep)  ! Objects_and_types

! Set defaults:
INFO = 0
FN_out = 1000
FN_distr = 9999

File_distr = 'OUTPUT_electron_distribution.dat'
File_out = 'OUT_entropy.dat'


!---------------------------------------
print*, '******************************************************************************'
print*, 'For analysis of entropy, call XTANT_entropy.exe'
print*, 'File OUTPUT_electron_distribution.dat must be present'
print*, '******************************************************************************'

!-----------------------------------
! Get the time grid from the file:
open (unit=FN_distr, file=trim(adjustl(File_distr)), status = 'old', readonly)
! Get the size of the array for the distribution and allocate arrays:
call get_size_of_distribution(FN_distr, Ei, distr, distr_eq, equilibr)  ! below
allocate(flnf(size(Ei)))

! Output file:
open (unit=FN_out, file=trim(adjustl(File_out)))
write(FN_out,'(a)') '# Time  Entropy  Equilibrium_entropy'
write(FN_out,'(a)') '# [fs]  [eV/K]   [eV/K]'

i = 0 ! block counter
read_well = .true. ! to start with
do while (read_well)
   i = i + 1
   ! Read the file with distribution block by block:
   read(FN_distr,*,IOSTAT=Reason) temp_ch, tim  ! read first line with timeprint
   if (Reason /= 0) then ! wrong format or end of file reached
       read_well = .false.
       exit
   else   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   end if

   ! Read the distributions at this timestep:
   call read_distributions(FN_distr, Ei, distr, distr_eq, read_well, equilibr)  ! below
   if (.not. read_well) exit

   ! Demarcation between blocks:
   read(FN_distr,*,IOSTAT=Reason)   ! skip empty line
   if (Reason /= 0) then ! wrong format or end of file reached
       read_well = .false.
       exit
   else   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   end if
   read(FN_distr,*,IOSTAT=Reason)   ! skip second empty line
   if (Reason /= 0) then ! wrong format or end of file reached
       read_well = .false.
       exit
   else   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   end if

   !-------------
   ! Now, process the data and print them out:
   S = 0.0d0   ! to start with
   flnf(:) = 0.0d0 ! to start with
   where(distr(:) > 0.0d0) flnf(:) = distr(:)*log(distr(:)/2.0d0)
   S = S + SUM(flnf) ! first part
   flnf(:) = 0.0d0   ! to restart
   where(distr(:) < 2.0d0) flnf(:) = (2.0d0 - distr(:))*log((2.0d0-distr(:))/2.0d0)
   S = S + SUM(flnf) ! second part
   !S = -g_kb * S  ! normalization factor
   S = -g_kb_EV * S  ! normalization factor [eV/K]

   S_eq = 0.0d0   ! to start with
   flnf(:) = 0.0d0 ! to start with
   where(distr_eq(:) > 0.0d0) flnf(:) = distr_eq(:)*log(distr_eq(:)/2.0d0)
   S_eq = S_eq + SUM(flnf) ! first part
   flnf(:) = 0.0d0   ! to restart
   where(distr_eq(:) < 2.0d0) flnf(:) = (2.0d0 - distr_eq(:))*log((2.0d0-distr_eq(:))/2.0d0)
   S_eq = S_eq + SUM(flnf) ! second part
   !S_eq = -g_kb * S_eq  ! normalization factor
   S_eq = -g_kb_EV * S_eq  ! normalization factor [eV/K]

   write(FN_out,'(es,es,es)') tim, S, S_eq
enddo
close(FN_distr)
close(FN_out)


!-------------------------
! Make gnuplot script:
FN_out1 = 9998
Gnu_script = 'OUT_entropy'
Gnu_file = 'OUT_entropy.png'
if (path_sep .EQ. '\') then	! if it is Windows
   open (unit=FN_out1, file=trim(adjustl(Gnu_script))//'.cmd')
   write(FN_out1, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   write(FN_out1, '(a)') 'set terminal png font arial'
   write(FN_out1, '(a)') 'set output "'//trim(adjustl(Gnu_file))//'"'
   write(FN_out1, '(a)') 'set xlabel "Time (fs)"'
   write(FN_out1, '(a)') 'set ylabel "Electron entropy (eV/K)"'
   write(FN_out1, '(a,a,a)') 'plot[][] "', trim(adjustl(File_out)), '" u 1:2 w l lw 3 title "Nonequilibrium" ,\'
   write(FN_out1, '(a,a,a)') '"', trim(adjustl(File_out)), '" u 1:3 w l lw 3 title "Equilibrium"'
else ! it is linux
   open (unit=FN_out1, file=trim(adjustl(Gnu_script))//'.sh')
   write(FN_out1, '(a)') '#!/bin/bash'
   write(FN_out1, '(a)') ''
   write(FN_out1, '(a)') 'NAME='//trim(adjustl(Gnu_file))
   write(FN_out1, '(a)') ' echo "'
   write(FN_out1, '(a)') 'set terminal png font arial'
   write(FN_out1, '(a)') 'set output \"$NAME\"'
   write(FN_out1, '(a)') 'set xlabel \"Time (fs) \" '
   write(FN_out1, '(a)') 'set ylabel \"Electron entropy \" '
   write(FN_out1, '(a,a,a)') 'plot[][] \"', trim(adjustl(File_out)), '\" u 1:2 w l lw 3 title \"Nonequilibrium\" ,\'
   write(FN_out1, '(a,a,a)') '\"', trim(adjustl(File_out)), '\" u 1:3 w l lw 3 title \"Equilibrium\"'
   write(FN_out1, '(a)') 'reset'
   write(FN_out1, '(a)') '" | gnuplot '
   CALL system('chmod +x '//trim(adjustl(Gnu_script))//'.sh') ! make the output-script executable
endif
close(FN_out1)

if (path_sep .EQ. '\') then	! if it is Windows
   CALL system(trim(adjustl(Gnu_script))//'.cmd')
else ! linux:
   CALL system(trim(adjustl(Gnu_script))//'.sh')
endif

2012 continue   ! to exit the program
STOP
!---------------------
 contains


subroutine read_distributions(FN_distr, Ei, distr, distr_eq, read_well, equilibr)  ! below
   integer, intent(in) :: FN_distr ! file number of the file with electronic distribution
   real(8), dimension(:), allocatable, intent(inout) :: Ei
   real(8), dimension(:), allocatable, intent(inout) :: distr
   real(8), dimension(:), allocatable, intent(inout) :: distr_eq
   logical, intent(inout) :: read_well
   logical, intent(in) :: equilibr  ! flag
   !-------------
   integer :: i, Reason

   Ei = 0.0d0  ! to start with
   distr = 0.0d0  ! to start with
   distr_eq = 0.0d0  ! to start with
   read_well = .true. ! to start with
   do i = 1, size(Ei)
      ! Read the file with distribution block by block:
      !read(FN_distr,*,IOSTAT=Reason) Ei(i), distr(i), distr_eq(i)
      if (equilibr) then
         read(FN_distr,*,IOSTAT=Reason) Ei(i), distr(i)
         distr_eq(i) = distr(i)
      else
         read(FN_distr,'(e,e,e)',IOSTAT=Reason) Ei(i), distr(i), distr_eq(i)
      endif

      if (Reason < 0) then ! end of file
         read_well = .false.
         exit
      endif
   enddo
end subroutine read_distributions



subroutine get_size_of_distribution(FN_distr, Ei, distr, distr_eq, equilibr)  ! below
   integer, intent(in) :: FN_distr ! file number of the file with electronic distribution
   real(8), dimension(:), allocatable, intent(inout) :: Ei
   real(8), dimension(:), allocatable, intent(inout) :: distr
   real(8), dimension(:), allocatable, intent(inout) :: distr_eq
   logical, intent(inout) :: equilibr  ! flag
   !-------------
   real(8) :: temp
   integer :: i, Reason, Ncol
   logical :: read_well

   ! Find if there is nonequilibrium distribution or not:
   call Count_columns_in_file(FN_distr, Ncol, skip_lines=1)   ! below
   if (Ncol == 2) then
      equilibr = .true. ! only equilibrium distribution is given
   else
      equilibr = .false.   ! both, non- and equilibrium distributions are there
   endif

   ! Skip first line:
   read(FN_distr,*,IOSTAT=Reason)

   i = 0
   read_well = .true. ! to start with
   do while (read_well)
      i = i + 1
      ! Read the file with distribution block by block:
      read(FN_distr,*,IOSTAT=Reason) temp
      if (Reason /= 0) then ! wrong format: probably new block
         exit
      endif
   enddo
   i = i - 1

   allocate(Ei(i))
   allocate(distr(i))
   allocate(distr_eq(i))
   !print*, 'i:', i, size(distr_eq)

   rewind(FN_distr) ! to restart reading into the arrays from the start
   !pause 'get_size_of_distribution'
end subroutine get_size_of_distribution



subroutine Count_columns_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of columns in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    real(8) temp
    character(1000) temp_ch
    integer i, Reason
    integer :: temp_i
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=605)
       enddo
       605 continue
    endif

    read(File_num,'(a)', IOSTAT=Reason) temp_ch ! count columns in this line
    N = number_of_columns(trim(adjustl(temp_ch))) ! see below

    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
end subroutine Count_columns_in_file



pure function number_of_columns(line)
   integer :: number_of_columns
   character(*), intent(in) :: line
   integer i, n
   logical :: same_space
   same_space = .false.
   i = 0
   n = len(line)
   number_of_columns = 0
   do while(i < n) ! scan through all the line
      i = i + 1
      selectcase (line(i:I))
      case (' ', '	') ! space or tab can be a separator between the columns
         if (.not.same_space) number_of_columns = number_of_columns + 1
         same_space = .true. ! in case columns are separated by more than one space or tab
      case default ! column data themselves, not a space inbetween
         same_space = .false.
      endselect
   enddo
   number_of_columns = number_of_columns + 1	! number of columns is by 1 more than number of spaces inbetween
end function number_of_columns



! Find out which OS it is:
subroutine Path_separator(path_sep)
   CHARACTER(len=1), intent(out) :: path_sep
   CHARACTER(len = 100) :: path
   CALL get_environment_variable("PATH",path)
   if (path(1:1) .EQ. '/') then        !unix based OS
       path_sep = '/'
   else if (path(3:3) .EQ. '\') then   !Windows OS
       path_sep = '\'
   else
       path_sep = ''
       print*, 'Path separator is not defined'    !Unknown OS
   endif 
end subroutine Path_separator

END PROGRAM Entropy
