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

PROGRAM XTANT_average_DOS
! Compilation:
!
! for DEBUG:
! ifort.exe -c /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_average_DOS.f90 /link /stack:9999999999
!
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_average_DOS.obj -o XTANT_average_DOS.exe /link /stack:9999999999
!
!
! for RELEASE:
! ifort.exe -c /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays XTANT_average_DOS.f90 /link /stack:9999999999
!
! ifort.exe /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays *.obj -o XTANT_average_DOS.exe /link /stack:9999999999
!
! To execute:
! XTANT_average_DOS.exe filename
! here "filename" is the file to be averaged.
! These files must be in the directories with identical names, except for version nmber at the end of them.
! The file must contain a few columns, the first of which must be time.
! It DOES NOT work for files with block-data, such as distribution function, powder-diffraction spectra, etc.
!<===========================================

USE IFLPORT, only : system, chdir

type Instant_data   ! from all files
   real(8), dimension(:), allocatable :: NRG            ! energy [eV]
   real(8), dimension(:,:), allocatable :: DOS          ! DOS
end type Instant_data

integer :: N_data  ! number of files
integer :: N_peaks  ! number of peaks

character(400) :: error_message, File_DOS, File_out
character(200) :: Gnu_script, Gnu_file, command, chtemp, File_name
character(32), dimension(10) :: char_var
character(10) :: temp_ch, call_slash, sh_cmd
character(1) :: path_sep

character(200), dimension(:), allocatable :: Folders_with_data
type(Instant_data), dimension(:), allocatable :: Read_data_DOS           ! All data and parameters at this timestep

real(8), dimension(:,:), allocatable :: ave_DOS         ! DOS
real(8), dimension(:), allocatable :: ave_NRG_grid      ! grid for averaged data

real(8) :: temp

integer :: FN_in, FN_out
integer :: INFO, Reason, i, j, k, INFO_flag
logical :: read_well, file_exist

call Path_separator(path_sep)  ! Objects_and_types

! Set defaults:
INFO = 0
FN_in = 9999
FN_out = 1000

File_DOS = 'OUTPUT_DOS.dat'   ! defaul name of XTANT out file with DOS
File_out = 'OUT_DOS.dat'   ! name of outpout file with averaged DOS

!---------------------------------------
print*, '******************************************************************************'
print*, 'Averaging DOS'
print*, '******************************************************************************'


! Get all the output folders names:
call collect_all_output(Folders_with_data)	! below
if (allocated(Folders_with_data)) then
    N_data = size(Folders_with_data)
else
    N_data = 0
endif

if (N_data < 1) then
    print*, 'No folders with data found. Terminating.'
    goto 2012   ! no data, nothing to do
endif

! Allocate arrays with input data:
allocate(Read_data_DOS(N_data))


! Now, read the data from each directory:
INFO_flag = 0   ! to start with
do i = 1, N_data   ! for all output data files
    File_name = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_DOS))

    inquire(file=trim(adjustl(File_name)),exist=file_exist)
    if (.not.file_exist) then
        print*, 'File not found: '//trim(adjustl(File_name))
        print*, 'Cannot continue, terminating.'
        goto 2012   ! no data, nothing to do
    endif

    open (unit=FN_in, file=trim(adjustl(File_name)), status = 'old', readonly) ! diffraction peaks
    print*, 'File #', i, trim(adjustl(File_name))

    ! Read from this file:
    call read_DOS_data(FN_in, Read_data_DOS, i, INFO)   ! below
    close(FN_in)    ! don't need this file anymore

    select case (INFO)
    case (2)
        print*, 'Error while reading File: '//trim(adjustl(File_name))
        print*, 'Cannot continue, terminating.'
        goto 2012   ! no data, nothing to do
    case (1)
        print*, 'Different number of columns in File: '//trim(adjustl(File_name))
        print*, 'Cannot continue, terminating.'
        goto 2012   ! no data, nothing to do
    case (-1)
        if (INFO_flag == 0) then ! print only the first time
           print*, 'Different number of lines in File: '//trim(adjustl(File_name))
           print*, 'Using different time grids...'
        endif
        INFO_flag = INFO    ! save it for later
        !print*, 'Cannot continue, terminating.'
        !goto 2012   ! no data, nothing to do
    case (0)
        !print*, 'Read well File: '//trim(adjustl(File_name))
    end select

    !pause 'i = 1, N_data'
enddo ! i = 1, N_data
if (INFO_flag == -1) INFO = INFO_flag   ! make sure the INFO accounts for different grids

! Now, average the data:
select case (INFO)
    case (-1)   ! different time grids used
        ! Create time grid appropriate for all used data:
        call create_grid(ave_DOS, ave_NRG_grid, Read_data_DOS)  ! below

        ! Now, average the data, taking into account different time grids:
        call get_average_data(ave_DOS, ave_NRG_grid, Read_data_DOS)  ! below

    case default ! identical time grids, simple averaging:
        ! Allocate the averaged data array:
        allocate(ave_DOS(size(Read_data_DOS(1)%DOS,1),size(Read_data_DOS(1)%DOS,2)), source = 0.0d0)
        do i = 1, size(ave_DOS,1)    ! for all peaks
            do j = 1, size(ave_DOS,2)    ! for all timesteps
                temp = 0.0d0    ! to start with
                do k = 1, size(Read_data_DOS)  ! collect all data:
                    temp = temp + Read_data_DOS(k)%DOS(i,j)
                    !print*, i, j, k, Read_data_DOS(k)%DOS(i,j)
                enddo
                ave_DOS(i,j) = temp/dble(N_data)
            enddo ! j
            !pause 'test 0'
        enddo ! i
end select


! Now, print out the data:
open (unit=FN_out, file=trim(adjustl(File_out))) ! averaged diffraction peaks
! Write the data:
select case (INFO)
    case (-1)   ! different time grids used
      do i = 1, size(ave_DOS,2)    ! for all timesteps
         write(FN_out, '(es24.16,$)') ave_NRG_grid(i), ave_DOS(:,i)    ! Time; Data for all peaks
         write(FN_out, '(a)') '' ! end line
      enddo ! i
    case default ! identical time grids
      do i = 1, size(ave_DOS,2)    ! for all timesteps
         write(FN_out, '(es24.16,$)') Read_data_DOS(1)%NRG(i), ave_DOS(:,i)    ! Time; Data for all peaks
         write(FN_out, '(a)') '' ! end line
      enddo ! i
end select
close(FN_out)    ! don't need this file anymore


2012 continue   ! to exit the program
STOP
!---------------------
 contains



subroutine get_average_data(ave_DOS, ave_NRG_grid, Read_data_DOS)
   real(8), dimension(:,:), allocatable, intent(inout) :: ave_DOS          ! diffraction peak intensity
   real(8), dimension(:), allocatable, intent(inout) :: ave_NRG_grid       ! time grid for averaged data
   type(Instant_data), dimension(:), intent(in) :: Read_data_DOS           ! All data and parameters at this timestep
   !----------------------
   real(8) :: temp, t_cur, temp_raw, t_first(size(ave_DOS,1))
   integer :: i, j, k, N_data, i_t

   N_data = size(Read_data_DOS) ! number of datasets

   do i = 1, size(ave_DOS,2)    ! for all NRG steps
      t_cur = ave_NRG_grid(i)  ! current NRG point on the grid
      do j = 1, size(ave_DOS,1)    ! for all peaks
         temp = 0.0d0    ! to start with
         do k = 1, N_data  ! collect all data:
             ! temp = temp + Read_data_DOS(k)%DOS(i,j)

             ! Find the point on the grid of raw data:
             call Find_in_monotonous_1D_array(Read_data_DOS(k)%NRG, t_cur, i_t) ! below

             !print*, i_t, t_cur, Read_data_DOS(k)%NRG(i_t)

             ! Interpolate raw data for the average-time-grid:
             if ( i_t < size(Read_data_DOS(k)%NRG) ) then
                temp_raw = Read_data_DOS(k)%DOS(j,i_t) + (Read_data_DOS(k)%DOS(j,i_t+1) - Read_data_DOS(k)%DOS(j,i_t)) / &
                            ( Read_data_DOS(k)%NRG(i_t+1) - Read_data_DOS(k)%NRG(i_t) ) * (t_cur-Read_data_DOS(k)%NRG(i_t))
             else ! last grid point
                temp_raw = Read_data_DOS(k)%DOS(j,i_t-1) + (Read_data_DOS(k)%DOS(j,i_t) - Read_data_DOS(k)%DOS(j,i_t-1)) / &
                            ( Read_data_DOS(k)%NRG(i_t) - Read_data_DOS(k)%NRG(i_t-1) ) * (t_cur-Read_data_DOS(k)%NRG(i_t-1))
             endif
             !print*, i, j, k, temp_raw

             temp = temp + temp_raw
         enddo
         ave_DOS(j,i) = temp/dble(N_data)
         !pause 'get_average_data'
     enddo ! j
   enddo ! i

end subroutine get_average_data




subroutine Find_in_monotonous_1D_array(Array, Value0, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2

   N = size(Array)
   i_1 = 1
   val_1 = Array(i_1)
   i_2 = N
   val_2 = Array(i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(i_cur)
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(i_cur), Array(i_1), Array(i_2)
        pause 'STOPPED WORKING...'
   else
       if (Value0 .LT. Array(1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LE. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = Array(i_1)
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                else
                   i_2 = i_cur
                   !val_2 = Array(i_2)
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                endif
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_1D_array', coun
                    write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(i_cur), Array(i_1), Array(i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_1D_array




subroutine create_grid(ave_DOS, ave_NRG_grid, Read_data_DOS)
   real(8), dimension(:,:), allocatable, intent(inout) :: ave_DOS          ! diffraction peak intensity
   real(8), dimension(:), allocatable, intent(inout) :: ave_NRG_grid       ! time grid for averaged data
   type(Instant_data), dimension(:), intent(in) :: Read_data_DOS           ! All data and parameters at this timestep
   !----------------------
   integer :: N_1, N_2, i
   real(8) :: start_time, end_time, cur_time, dt, dt_cur

   start_time = Read_data_DOS(1)%NRG(1)    ! to start with
   do i = 1, size(Read_data_DOS) ! find largest start-time
      cur_time = Read_data_DOS(i)%NRG(1)
      if (cur_time > start_time) start_time = cur_time
      !print*, 'start:', i, cur_time, start_time
   enddo

   end_time = Read_data_DOS(1)%NRG(size(Read_data_DOS(1)%NRG))    ! to start with
   do i = 1, size(Read_data_DOS) ! find smallest end-time
      cur_time = Read_data_DOS(i)%NRG(size(Read_data_DOS(i)%NRG))
      if (cur_time < end_time) end_time = cur_time
      !print*, 'end:', i, cur_time, end_time
   enddo

   dt = Read_data_DOS(1)%NRG(2) - Read_data_DOS(1)%NRG(1) ! to start with
   do i = 1, size(Read_data_DOS) ! find smallest end-time
      dt_cur = Read_data_DOS(i)%NRG(2) - Read_data_DOS(i)%NRG(1)
      if (dt_cur < dt) dt = dt_cur
      !print*, 'dr:', i, dt_cur, dt
   enddo

   ! Now, set the grid for average data:
   N_1 = size(Read_data_DOS(1)%DOS,1) ! number of peaks
   N_2 = int( (end_time - start_time)/dt )  ! number of timesteps

   allocate(ave_NRG_grid(N_2), source = 0.0d0)
   allocate(ave_DOS(N_1, N_2), source = 0.0d0)

   ! Set the time grid:
   ave_NRG_grid(1) = start_time
   do i = 2, N_2
      ave_NRG_grid(i) = ave_NRG_grid(i-1) + dt
   enddo

   !print*, start_time, end_time, dt, N_1, N_2
   !pause 'create_grid'
end subroutine create_grid


subroutine read_DOS_data(FN_in, Read_data_DOS, i, INFO)
    integer, intent(in) :: FN_in ! file number (must be already openned)
    type(Instant_data), dimension(:), intent(inout) :: Read_data_DOS           ! All data and parameters at this timestep
    integer, intent(in) :: i    ! data file index
    integer, intent(inout) :: INFO  ! info flag:
    ! 0     =   all good
    ! -1    =   error: number of lines (time grid is not identical in the files)
    ! 1     =   error: number of columns (different number of diffraction peaks)
    ! 2     =   error while reading file
    !-----------------------
    integer :: j, i_siz_1, i_col_1, N_col, N_lines, Reason
    character(200) :: read_string
    character(10) :: temp_ch

    ! Get the nuimber of diffraction peaks in the file:
    call Count_columns_in_file(FN_in, N_col) ! below
    N_col = N_col - 1   ! exclude the first column (with timestep)

    ! Get the nuimber of diffraction peaks in the file:
    call Count_lines_in_file(FN_in, N_lines, skip_lines=0) ! below
    !print*, 'N_col:', N_col
    !print*, 'N_lines:', N_lines

    ! Allocate the arrays:
    allocate(Read_data_DOS(i)%NRG(N_lines))
    allocate(Read_data_DOS(i)%DOS(N_col,N_lines))
    ! Compare it to the first data file:
    i_siz_1 = size(Read_data_DOS(1)%NRG)
    i_col_1 = size(Read_data_DOS(1)%DOS,1)
    if (N_lines /= i_siz_1) then
        INFO = -1
        !return ! different time grids are allowed!
    endif
    if (N_col /= i_col_1) then
        INFO = 1
        return
    endif

    ! Read the data:
    do j = 1, N_lines
        read(FN_in, *, IOSTAT=Reason) Read_data_DOS(i)%NRG(j), Read_data_DOS(i)%DOS(:,j)
        if (Reason .LT. 0) then
            print*, 'Error in read #1: ', Reason, 'Line #:', j
            INFO = 2
            return
        endif
        !print*, j, Read_data_DOS(i)%NRG(j), Read_data_DOS(i)%DOS(:,j)
    enddo
end subroutine read_DOS_data





subroutine collect_all_output(Folders_with_data)
   character(200), dimension(:), allocatable, intent(inout) :: Folders_with_data
   character(1) :: path_sep
   character(500) :: File_name, command, read_line, temp_file
   integer :: FN, open_status, leng, Reason, count_lines, i
   FN = 1300
   ! Find out which OS it is:
   call Path_separator(path_sep)  ! Objects_and_types

   ! Create a temporary file to store the list of all data files from the folder:
   temp_file = 'List_of_folders.txt'
   File_name = trim(adjustl(temp_file))

   ! Get the names of all data files in the folder using system commands:
   if (path_sep .EQ. '\') then	! if it is Windows
      command = 'dir OUTPUT_* /b >'//trim(adjustl(File_name))
   else
      command = "ls -t | grep 'OUTPUT_' >"//trim(adjustl(File_name))
   endif
   !call system(trim(adjustl(command))) ! execute the command
   i = system(trim(adjustl(command))) ! execute the command

   ! Read file names:
   open(UNIT=FN, file=trim(adjustl(File_name)), iostat=open_status, action='read')
   TEMP:if ( open_status /= 0 ) then
      print *, 'Could not open ',trim(adjustl(File_name)),' to get the list of folders with all data'
   else TEMP
      ! Read all lines in the file one by one:
      call Count_lines_in_file(FN, count_lines)

      if (allocated(Folders_with_data)) deallocate(Folders_with_data)
      allocate(Folders_with_data(count_lines))
      do i =1,count_lines
         read(FN,'(a)',IOSTAT=Reason) Folders_with_data(i)  ! all folders with the outputs of XTANT
         if (Reason < 0) exit
      enddo
      close(FN, status='delete') ! temp file is not needed anymore, erase it
   endif TEMP
end subroutine collect_all_output




subroutine order_of_time(tim, text, gnu_text, x_tics)
   real(8), intent(in) :: tim ! time to find its order
   character(*), intent(out) :: text ! fs, ps, ns, mks, ms, s
   character(*), intent(out), optional :: gnu_text ! culomn to set in gnuplot
   real(8), intent(out), optional :: x_tics ! tics for gnuplot
   integer :: time_ord
   time_ord = find_order_of_number_real(tim) ! module "Little_subroutines"
   if (present(x_tics)) then
      x_tics = 10.0d0**(time_ord) ! set tics for gnuplot
      if (tim/dble(x_tics) > 0.5) then
         x_tics = 10.0d0**(time_ord-1) ! set tics for gnuplot
      else if (tim/dble(x_tics) > 0.2) then
         x_tics = 0.5d0*10.0d0**(time_ord-1) ! set tics for gnuplot
      else
         x_tics = 10.0d0**(time_ord-2) ! set tics for gnuplot
      endif
   endif

   if (time_ord > 1e15) then ! s
      text = '(s)'
      if (present(gnu_text)) gnu_text = '($1/1e15)'
   else if (time_ord > 1e12) then ! ms
      text = '(ms)'
      if (present(gnu_text)) gnu_text = '($1/1e12)'
   else if (time_ord > 1e9) then ! mks
      text = '(mks)'
      if (present(gnu_text)) gnu_text = '($1/1e9)'
   else if (time_ord > 1e6) then ! ns
      text = '(ns)'
      if (present(gnu_text)) gnu_text = '($1/1e6)'
   else if (time_ord > 1e3) then ! ps
      text = '(ps)'
      if (present(gnu_text)) gnu_text = '($1/1e3)'
   else ! fs
      text = '(fs)'
      if (present(gnu_text)) gnu_text = '($1)'
   endif
end subroutine order_of_time



pure function find_order_of_number_real(num)
   integer find_order_of_number_real
   real(8), intent(in) :: num
   character(64) :: temp
   write(temp,'(i8)') CEILING(num) ! make it a string
   find_order_of_number_real = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
end function find_order_of_number_real

pure function find_order_of_number_int(num)
   integer find_order_of_number_int
   integer, intent(in) :: num
   character(64) :: temp
   write(temp,'(i8)') num ! make it a string
   find_order_of_number_int = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
end function find_order_of_number_int



subroutine write_gnuplot_script_header_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, path_sep, setkey, set_x_log, set_y_log, fontsize)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   integer, intent(in), optional :: setkey, fontsize
   logical, intent(in), optional :: set_x_log, set_y_log
   !---------------------------------
   integer :: font_size, set_key
   logical :: x_log, y_log


   if (present(fontsize)) then   ! user-set font size
      font_size = fontsize
   else  ! default font size
      font_size = 14
   endif

   if (present(setkey)) then
      set_key = setkey
   else  ! default
      set_key = 1
   endif

   if (present(set_x_log)) then
      x_log = set_x_log
   else  ! default
      x_log = .false.
   endif

   if (present(set_y_log)) then
      y_log = set_y_log
   else  ! default
      y_log = .false.
   endif

   if (path_sep .EQ. '\') then	! if it is Windows
      call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, set_key, x_log, y_log, font_size)
   else ! it is linux
      call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, set_key, x_log, y_log, font_size)
   endif
end subroutine write_gnuplot_script_header_new


subroutine write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, set_x_log, set_y_log, font_size)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey, font_size
   logical, intent(in), optional :: set_x_log, set_y_log
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp, temp2

   if (present(font_size)) then
      write(temp2,'(i0)') font_size
   else
      write(temp2,'(i0)') 14
   endif

   select case (ind)
   case(1:) ! any file format
      write(FN, '(a)') '#!/bin/bash'
      write(FN, '(a)') ''
      write(FN, '(a)') 'NAME='//trim(adjustl(Out_file))
   end select
   write(FN, '(a,f3.1)') 'LW=', LW
   write(FN, '(a)') 'LABL="'//trim(adjustl(labl))//'"'
   write(temp, '(f12.2)') x_tics
   write(FN, '(a)') 'TICSIZ='//trim(adjustl(temp))
   write(FN, '(a)') 'echo " '
   select case (ind)
      case (1)  ! eps
         write(FN, '(a)') 'set terminal postscript enhanced \"Helvetica\" 16 color '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (2)  ! jpeg
         write(FN, '(a)') 'set terminal jpeg font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (4)  ! png
         !write(FN, '(a)') 'set terminal png font \"arial,14\" '
         write(FN, '(a)') 'set terminal pngcairo dashed font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (6)  ! animated gif
         write(FN, '(a)') 'set terminal gif animate delay 10 font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (0)
         write(FN, '(a)') 'set terminal x11 persist'
         write(FN, '(a)') 'unset label'
   endselect
!    write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//' \"        font \"Helvetica,20\" '
!    write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//' \"      font \"Helvetica,20\" '
   write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//'\" font \"arial,18\" '
   write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//'\" font \"arial,18\" '

   !write(FN, '(a)') 'set label \"$LABL\" at 150,-8 font \"Helvetica,22\" '
   if (present(setkey)) then
      select case(setkey)
      case (1)
         write(FN, '(a)') 'set key right bottom '
      case (2)
         write(FN, '(a)') 'set key left top '
      case (3)
         write(FN, '(a)') 'set key left bottom '
      case (4)
         write(FN, '(a)') 'unset key '
      case default
         write(FN, '(a)') 'set key right top '
      endselect
   else
      write(FN, '(a)') 'set key right top '
   endif

   if (present(set_x_log)) then
      if (set_x_log) then
         write(FN, '(a)') "set logscale x"
         write(FN, '(a)') 'set format x \"10^{\%L}\"'
      endif
   endif

   if (present(set_y_log)) then
      if (set_y_log) then
         write(FN, '(a)') "set logscale y"
         write(FN, '(a)') 'set format y \"10^{\%L}\"'
      endif
   endif

   write(FN, '(a)') 'set xtics \"$TICSIZ\" '
end subroutine write_gnuplot_script_header_linux_new



subroutine write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, set_x_log, set_y_log, font_size)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey, font_size
   logical, intent(in), optional :: set_x_log, set_y_log
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp, temp2


   if (present(font_size)) then
      write(temp2,'(i0)') font_size
   else
      write(temp2,'(i0)') 14
   endif


   select case (ind)
   case(1:)	! eps
      write(FN, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   end select
   write(FN, '(a,f3.1)') 'LW=', LW

    select case (ind)
      case (1)  ! eps
         write(FN, '(a)') 'set terminal postscript enhanced "Helvetica,'//trim(adjustl(temp2))//'" color '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (2)  ! gpeg
         write(FN, '(a)') 'set terminal jpeg large font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif large font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (4)  ! png
         !write(FN, '(a)') 'set terminal png font "arial,14" '
         write(FN, '(a)') 'set terminal pngcairo dashed font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (6)  ! animated gif
         write(FN, '(a)') 'set terminal gif animate delay 10 font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (0)
         write(FN, '(a)') 'set terminal x11 persist'
         write(FN, '(a)') 'unset label'
   endselect
   write(FN, '(a)') 'set xlabel "'//trim(adjustl(xlabl))//'" font "arial,18"'
   write(FN, '(a)') 'set ylabel "'//trim(adjustl(ylabl))//'" font "arial,18"'

   !write(FN, '(a)') 'set label \"$LABL\" at 150,-8 font \"Helvetica,22\" '
   if (present(setkey)) then
      select case(setkey)
      case (1)
         write(FN, '(a)') 'set key right bottom '
      case (2)
         write(FN, '(a)') 'set key left top '
      case (3)
         write(FN, '(a)') 'set key left bottom '
      case (4)
         write(FN, '(a)') 'unset key '
      case default
         write(FN, '(a)') 'set key right top '
      endselect
   else
      write(FN, '(a)') 'set key right top '
   endif

   if (present(set_x_log)) then
      if (set_x_log) then
         write(FN, '(a)') "set logscale x"
         write(FN, '(a)') 'set format x "10^{%L}"'
      endif
   endif

   if (present(set_y_log)) then
      if (set_y_log) then
         write(FN, '(a)') "set logscale y"
         write(FN, '(a)') 'set format y "10^{%L}"'
      endif
   endif

   write(temp, '(f12.2)') x_tics
   write(FN, '(a,a)') 'set xtics ', trim(adjustl(temp))
end subroutine write_gnuplot_script_header_windows_new


subroutine write_gnuplot_script_ending(FN, File_name, ind, path_sep)
   integer, intent(in) :: FN, ind
   character(*), intent(in) :: File_name
   character(*), intent(in) :: path_sep
   character(100) :: command
   integer :: iret

   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      select case (ind)
      case (:1)
         write(FN, '(a)') 'reset'
         write(FN, '(a)') '" | gnuplot '
         !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
         command = 'chmod +x '//trim(adjustl(File_name))
         iret = system(command)

         !call system(trim(adjustl(File_name))) ! execute the prepared script
      case (2:)
         write(FN, '(a)') 'reset'
         write(FN, '(a)') '" | gnuplot '
         !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
         command = 'chmod +x '//trim(adjustl(File_name))
         iret = system(command)
         !call system(trim(adjustl(File_name))) ! execute the prepared script
      endselect
   endif
end subroutine write_gnuplot_script_ending



subroutine create_file_header(FN, text)
   integer, intent(in) :: FN		! file number
   character(*), intent(in) :: text	! what to write in this file
   write(FN,'(a)') trim(adjustl(text))
end subroutine create_file_header



subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i, Reason
    real(8) :: temp

    if (present(skip_lines)) then ! in case you want to skip some comment lines and count only lines with data
       do i=1,skip_lines
          read(File_num, *, end=604)
       enddo
       604 continue

       ! assume real nmumbers in the file:
       i = 0
       do
           read(File_num, *, IOSTAT=Reason, end=605) temp
           if (Reason /= 0) exit ! found something that's not real number, probably end of file
           i = i + 1
       enddo
       605 continue
       rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
       N = i
    else ! nothing to skip, use universal routine:
       i = 0
       do
           read(File_num,*, end=603)
           i = i + 1
       enddo
       603 continue
       rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
       N = i
    endif
end subroutine Count_lines_in_file



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

END PROGRAM XTANT_average_DOS
