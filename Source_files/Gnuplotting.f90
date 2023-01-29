! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2021 Nikita Medvedev
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

implicit none 

 contains
 
 
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
 
 

subroutine write_gnuplot_script_header_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, path_sep, setkey)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   integer, intent(in), optional :: setkey
      
   if (present(setkey)) then
      if (path_sep .EQ. '\') then	! if it is Windows
         call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey)
      else ! it is linux
         call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey)
      endif
   else
      if (path_sep .EQ. '\') then	! if it is Windows
         call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file)
      else ! it is linux
         call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file)
      endif
   endif
end subroutine write_gnuplot_script_header_new


subroutine write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp
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
      case (2)  ! gpeg
         write(FN, '(a)') 'set terminal jpeg font \"arial,14\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif font \"arial,14\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (4)  ! png
         write(FN, '(a)') 'set terminal png font \"arial,14\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font \"arial,14\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (6)  ! animated gif
         write(FN, '(a)') 'set terminal gif animate delay 10 font \"arial,14\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (0)
         write(FN, '(a)') 'set terminal x11 persist'
         write(FN, '(a)') 'unset label'
   endselect
!    write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//' \"        font \"Helvetica,20\" '
!    write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//' \"      font \"Helvetica,20\" '
   write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//' \" '
   write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//' \" '
   
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
   write(FN, '(a)') 'set xtics \"$TICSIZ\" '
end subroutine write_gnuplot_script_header_linux_new



subroutine write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp
   select case (ind)
   case(1:)	! eps
      write(FN, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   end select
   write(FN, '(a,f3.1)') 'LW=', LW

    select case (ind)
      case (1)  ! eps
         write(FN, '(a)') 'set terminal postscript enhanced "Helvetica,14" color '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (2)  ! gpeg
         write(FN, '(a)') 'set terminal jpeg large font "arial,14" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif large font "arial,14" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (4)  ! png
         write(FN, '(a)') 'set terminal png font "arial,14" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font "arial,14" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (6)  ! animated gif
         write(FN, '(a)') 'set terminal gif animate delay 10 font "arial,14" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (0)
         write(FN, '(a)') 'set terminal x11 persist'
         write(FN, '(a)') 'unset label'
   endselect
   write(FN, '(a)') 'set xlabel "'//trim(adjustl(xlabl))//' " '
   write(FN, '(a)') 'set ylabel "'//trim(adjustl(ylabl))//' " '
   
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
   !write(FN, '(a,f6.2)') 'set xtics ', x_tics
   write(temp, '(f12.2)') x_tics
   write(FN, '(a,a)') 'set xtics ', trim(adjustl(temp))
end subroutine write_gnuplot_script_header_windows_new


subroutine write_gnuplot_script_ending_new(FN, File_name, path_sep)
   integer, intent(in) :: FN
   character(*), intent(in) :: File_name
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   
   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      write(FN, '(a)') 'reset'
      write(FN, '(a)') '" | gnuplot '
      call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
end subroutine write_gnuplot_script_ending_new

 
 

end module Gnuplotting
