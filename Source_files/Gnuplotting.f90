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
module Gnuplotting

USE IFLPORT, only : system, chdir   ! library, allowing to operate with directories in intel fortran
use Dealing_with_files, only: Count_lines_in_file


implicit none 
PRIVATE

public :: write_gnu_printout, write_gnuplot_script_header_new, collect_gnuplots


 contains



!===================================================
! Gnuplotting all the scripts:
subroutine collect_gnuplots(path_sep, out_path)
   character(*), intent(in) :: path_sep, out_path    ! folder with the cmd-files
   !------------------------
   character(200) :: File_name, command, Gnuplot_all_file
   integer :: FN, N_f, i, n_slash
   integer :: open_status, iret, idir, leng
   character(200), dimension(:), allocatable :: All_files
   character(300) :: output_path
   character(5) ::  call_slash, sh_cmd

   ! In which directory to collect all gnu scripts:
   output_path = out_path

   ! Create a temporary file:
   Gnuplot_all_file = 'OUTPUT_Gnuplot_all'

   ! Find the extension of the gnuplot scripts:
   call cmd_vs_sh(trim(adjustl(path_sep)), call_slash, sh_cmd)  ! below
   ! Include the extension of the script:
   Gnuplot_all_file = trim(adjustl(Gnuplot_all_file))//trim(adjustl(sh_cmd))

   ! Include the path to the directory:
   File_name = trim(adjustl(output_path))//trim(adjustl(path_sep))//trim(adjustl(Gnuplot_all_file))

   ! Save the names of all gnuplot scripts into this file:
   if (trim(adjustl(path_sep)) == '\') then  ! if it is Windows
      command = 'dir '//trim(adjustl(output_path))//'\*'//trim(adjustl(sh_cmd))//' /b >'//trim(adjustl(File_name))
   else ! linux:
      command = "ls -t "//trim(adjustl(output_path))//" | grep '"//trim(adjustl(sh_cmd))//"' >"//trim(adjustl(File_name))
   endif

   iret = system(trim(adjustl(command)))   ! execute the command to save file names in the temp file
   !call system(trim(adjustl(command))) ! execute the command to save file names in the temp file

   ! Open the files with gnuplot script names:
   open(NEWUNIT=FN, file=trim(adjustl(File_name)), iostat=open_status)
   if ( open_status /= 0 ) then
      print *, 'Could not open ',trim(adjustl(File_name)),' for gnuplotting.', ' Unit = ', FN
   endif

   ! Find out how many there are:
   call Count_lines_in_file(FN, N_f) ! module "Dealing_with_files"

   ! Allocate array with them:
   allocate(All_files(N_f)) ! array with all relevant file names

   ! Read file names:
   do i = 1,N_f
      read(FN,*) All_files(i)
   enddo

   ! Rewind file to overwrite including the calls:
   rewind(FN)
   ! Make the script executable:
   if (trim(adjustl(path_sep)) == '\') then  ! if it is Windows
      write(FN,'(a)') '@echo off'
   else
      write(FN,'(a)') '#!/bin/bash'
   endif
   do i = 1,N_f
      if (trim(adjustl(All_files(i))) /= trim(adjustl(Gnuplot_all_file))) then ! exclude the file itself
         if (trim(adjustl(path_sep)) == '\') then  ! if it is Windows
            write(FN,'(a)') trim(adjustl(call_slash))//' '//trim(adjustl(All_files(i)))
         else
            leng = LEN(trim(adjustl(All_files(i))))
            if (trim(adjustl(All_files(i)(leng-2:))) == '.sh') then  ! to exclude other files possible containing "sh"
               write(FN,'(a)') trim(adjustl(call_slash))//trim(adjustl(All_files(i)))
            endif
         endif
      endif
   enddo
   close (FN)
   if (trim(adjustl(path_sep)) /= '\') then	! if it is Linux
      iret = system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
   !pause 'Execute all'

   !--------------
   ! Execute all the gnuplot scripts:
   idir = chdir(trim(adjustl(output_path))) ! go into the directory with output files
   !call chdir(trim(adjustl(output_path))) ! go into the directory with output files

   if (trim(adjustl(path_sep)) == '\') then	! if it is Windows
      iret = system( '@echo off' )   ! create the folder
      iret = system(trim(adjustl(call_slash))//' '//trim(adjustl(Gnuplot_all_file)))   ! create the folder
      !call system( '@echo off' )
      !call system(trim(adjustl(call_slash))//' '//trim(adjustl(Gnuplot_all_file)))
   else ! linux:
      iret = system( '#!/bin/bash' )
      iret = system(trim(adjustl(call_slash))//trim(adjustl(Gnuplot_all_file)))
      !call system( '#!/bin/bash' )
      !call system(trim(adjustl(call_slash))//trim(adjustl(Gnuplot_all_file)))
   endif

   ! Count how many times the system has to go out of the directory to get back into the original directory:
   ! Defined by the number of slashes in the path given:
   n_slash = count( (/ (trim(adjustl(output_path(i:i))), i=1,len_trim(output_path)) /) == trim(adjustl(path_sep)) )
   do i = 1, n_slash+1  ! go up in the directory tree as many times as needed
      idir = chdir("../")    ! exit the directory with output files
      !call chdir("../")    ! exit the directory with output files
   enddo
end subroutine collect_gnuplots


pure subroutine cmd_vs_sh(path_sep, call_slash, sh_cmd)
   character(*), intent(in) :: path_sep
   character(*), intent(out) :: call_slash, sh_cmd
   if (path_sep .EQ. '\') then	! if it is Windows
      call_slash = 'call '
      sh_cmd = '.cmd'
   else ! It is linux
      call_slash = './'
      sh_cmd = '.sh'
   endif
end subroutine cmd_vs_sh

 
subroutine write_gnu_printout(FN, first_line, last_line, file_name, x_start, x_end, y_start, y_end, col_x, col_y, lw, title, additional_info, linux_s)
   integer, intent(in) :: FN ! file nuumber to write to
   logical, intent(in) :: first_line, last_line ! is it the first line to plot? last line? for multiple curves on the same plot
   character(*) file_name ! file to plot from
   real(8), intent(in), optional :: x_start, x_end, y_start, y_end ! starting and ending points to put on axes
   character(*), intent(in), optional :: col_x, col_y ! which columns from the file to plot
   integer, intent(in), optional :: lw ! line width
   character(*), intent(in), optional :: title ! title to print on the plot
   character(*), intent(in), optional :: additional_info ! any additional gnuplot operators to pass
   logical, intent(in), optional :: linux_s ! is it linux system?
   character(500) :: format_var
   character(50) :: temp
   
   format_var = ''
   if (first_line) then 
      format_var = 'p [' 
      if (present(x_start)) then
         write(temp,'(es24.16)') x_start
         format_var = trim(adjustl(format_var))//trim(adjustl(temp))//':'
      endif
      if (present(x_end)) then
         write(temp,'(es24.16)') x_end
         format_var = trim(adjustl(format_var))//trim(adjustl(temp))
      endif
      format_var = trim(adjustl(format_var))//']['
      if (present(y_start)) then
         write(temp,'(es24.16)') y_start
         format_var = trim(adjustl(format_var))//trim(adjustl(temp))//':'
      endif
      if (present(y_end)) then
         write(temp,'(es24.16)') y_end
         format_var = trim(adjustl(format_var))//trim(adjustl(temp))
      endif
      format_var = trim(adjustl(format_var))//']'
   endif
   ! File name to get data from:
   if (present(linux_s)) then
      format_var = trim(adjustl(format_var))//' \"'//trim(adjustl(file_name))//'\" '
   else ! windows
      format_var = trim(adjustl(format_var))//' "'//trim(adjustl(file_name))//'" '
   endif
   ! Column to use for x
   format_var = trim(adjustl(format_var))//' u '//trim(adjustl(col_x))
   ! Column to use for y
   format_var = trim(adjustl(format_var))//':'//trim(adjustl(col_y))
   
   ! Any additional settings
   if (present(additional_info)) format_var = trim(adjustl(format_var))//trim(adjustl(additional_info))
   
   if (present(lw)) then ! line width
      write(temp,'(i2)') lw
   else ! default value:
      write(temp,'(i1)') 3
   endif
   format_var = trim(adjustl(format_var))//' w l lw '//trim(adjustl(temp))
   if (present(title)) then 
      if (present(linux_s)) then
         format_var = trim(adjustl(format_var))//' title \"'//trim(adjustl(title))//'\"'
      else ! windows
         format_var = trim(adjustl(format_var))//' title "'//trim(adjustl(title))//'"'
      endif
   endif
   if (.not.last_line) format_var = trim(adjustl(format_var))//' ,\ '
   
   write(FN, '(a)') trim(adjustl(format_var))
end subroutine write_gnu_printout
 
 

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


subroutine write_gnuplot_script_ending_new(FN, File_name, path_sep)
   integer, intent(in) :: FN
   character(*), intent(in) :: File_name
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   !--------------
   integer :: iret

   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      write(FN, '(a)') 'reset'
      write(FN, '(a)') '" | gnuplot '
      !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
      iret = system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
end subroutine write_gnuplot_script_ending_new



subroutine get_extension_index(text_ext, ind)
   character(*), intent(in) :: text_ext   ! extension
   integer, intent(out) :: ind   ! internal index
   !--------------------
   select case (trim(adjustl(text_ext)))
   case ('eps', 'EPS')  ! eps
      ind = 1
   case ('JPG', 'jpg', 'JPEG', 'jpeg')  ! jpeg
      ind = 2
   case ('GIF', 'gif')  ! gif
      ind = 3
   case ('PNG', 'png')  ! png
      ind = 4
   case ('PDF', 'pdf')  ! pdf
      ind = 5
   case ('animated_gif')  ! animated gif
      ind = 6
   case default ! default
      print*, 'Using default gnuplot format for plots: jpeg'
      ind = 2 ! exclude
   endselect
end subroutine get_extension_index

end module Gnuplotting
