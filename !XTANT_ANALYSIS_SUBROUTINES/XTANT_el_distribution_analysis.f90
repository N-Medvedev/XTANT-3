PROGRAM XTANT_el_distribution_analysis
! Compilation:
!
! for DEBUG:
! ifort.exe -c /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_el_distribution_analysis.f90 /link /stack:9999999999
!
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs *.obj -o XTANT_el_distribution_analysis.exe /link /stack:9999999999

! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_el_distribution_analysis.f90 -o XTANT_el_distribution_analysis.exe /link /stack:9999999999

! for RELEASE:
! ifort.exe -c /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays XTANT_el_distribution_analysis.f90 /link /stack:9999999999
!
! ifort.exe /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays *.obj -o XTANT_el_distribution_analysis.exe /link /stack:9999999999
!
! To execute:
! XTANT_el_distribution_analysis.exe
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
real(8), dimension(:), allocatable :: distr
real(8), dimension(:), allocatable :: distr_conv
real(8), dimension(:), allocatable :: distr_aver
real(8), dimension(:), allocatable :: weights_read
real(8) :: tim, tim_Start, Conv_dE, temp_r
integer :: FN_distr, FN_out_conv, FN_out_average, FN_weights, weights_COL
integer :: INFO, Reason, i, Nsiz, tim_counter
logical :: file_exist, read_well, print_conv, file_opened
character(100) :: File_distr, File_out_conv, File_out_average, File_weights, File_in
character(100) :: File_out_conv_gnu, File_out_average_gnu
character(1) :: path_sep, temp_ch

call Path_separator(path_sep)  ! below

! File names:
File_in           = 'EL_DISTR.txt'  ! file with parameters to use for post-processing (see XTANT Manual)
File_distr        = 'OUTPUT_electron_distribution_on_grid.dat'
File_out_conv     = 'OUT_el_distr_vs_time.dat'
File_out_average  = 'OUT_el_distr_average.dat'
File_out_conv_gnu     = 'OUT_el_distr_vs_time'
File_out_average_gnu  = 'OUT_el_distr_average'


! File numbers:
FN_weights = 9998
FN_distr = 9999
FN_out_conv = 1000
FN_out_average = 2000

! Set defaults:
INFO = 0

!---------------------------------------
print*, '******************************************************************************'
print*, 'For analysis of entropy, call XTANT_el_distribution_analysis.exe'
print*, 'Files required by the program:'
print*, 'EL_DISTR.txt -- file with parameters (see XTANT manual)'
print*, 'OUTPUT_electron_distribution_on_grid.dat -- file with fe(E)'
print*, '******************************************************************************'

!-----------------------------------
! Read input:
inquire(file=trim(adjustl(File_in)),exist=file_exist) ! check if input file is there
if (file_exist) then ! try to read it, if there is a grid provided
   call read_EL_DISTR_txt(File_in, Conv_dE, File_weights, FN_weights, weights_COL) ! below
   ! Knowing the parameters, define flag to mark if convolution is required:
   print_conv = (Conv_dE > 0.0d0)   ! only if convolution is used
   ! Open output file(s):
   call creat_output_files(FN_out_conv, File_out_conv, FN_out_average, File_out_average, print_conv)  ! below
   print*, 'Output files are created, starting the analysis'
else
   print*, 'File ', trim(adjustl(File_in)), 'not found, program terminates.'
   print*, 'For execution of the program, this file must be present,'
   print*, 'containing the following 2 lines:'
   print*, 'Line #1: Energy width in [eV]. Set a negative number to exclude convolution in the energy space, if not needed.'
   print*, 'Line #2: File name, providing the data on the widths in time to average the distribution with;'
   print*, 'and the index of column to use in this file.'
   print*, 'Set a non-existing file name to exclude averaging in time, if not needed.'
   print*, '******************************************************************************'
   goto 2012   ! exit, nothing else to do
endif

! Get the energy grid from the file:
open (unit=FN_distr, file=trim(adjustl(File_distr)), status = 'old', readonly)
! Get the size of the array for the distribution and allocate arrays:
call get_size_of_distribution(FN_distr, Ei, distr)  ! below
print*, 'Distribution size is known, proceeding to reading files'
! Knowing the size, define convolved distribution array:
allocate(distr_conv(size(distr)), source = 0.0d0)
allocate(distr_aver(size(distr)), source = 0.0d0)

! Read weights from the file:
if (weights_COL > 0) then
   call read_weigths_file(FN_weights, weights_read, weights_COL)  ! below
   print*, 'Weights are read from the file: ', trim(adjustl(File_weights))
endif

! Output files:
if (print_conv) open (unit=FN_out_conv, file=trim(adjustl(File_out_conv)))
open (unit=FN_out_average, file=trim(adjustl(File_out_average)))

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
   if (print_conv) write(FN_out_conv,*) temp_ch, tim  ! copy the same markerline (timestemp)
   if (i == 1) tim_Start = tim ! save starting time for gnuplotting below

   print*, 'Reading and analysing distribution, block #', i

   ! Read the distributions at this timestep:
   call read_distribution(FN_distr, Ei, distr, read_well)  ! below
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
   ! If the convolution is required:
   if (print_conv) then
      ! Convolve with Gaussian:
      call convolve_with_Gaussian(Ei, distr, Conv_dE, distr_conv) ! below
      ! Printout the convolved function:
      call print_convolved(FN_out_conv, Ei, distr_conv)   ! below

      ! Average convolved function over time:
      call average_distr(distr_aver, distr_conv, weights_read, weights_COL, i)   ! below
   else
      ! Average original function over time:
      call average_distr(distr_aver, distr, weights_read, weights_COL, i)   ! below
   endif
enddo
! Average distribution over time:
distr_aver(:) = distr_aver(:) / dble(i)
! And let's normalize it to the peak height:
temp_r = maxval(distr_aver)
distr_aver(:) = distr_aver(:) / temp_r

! Clean up:
close(FN_distr)
if (print_conv) close(FN_out_conv)

! Printout the avereaged distribution:
call print_averaged(FN_out_average, Ei, distr_aver)   ! below
close(FN_out_average)

print*, 'Analysis is done, starting gnuplotting'

!-------------------------
! Make gnuplot script:
call gnuplot_figures(path_sep, File_out_conv, File_out_average, File_out_conv_gnu, File_out_average_gnu, print_conv, tim_Start)   ! below
print*, 'Everything is done, check the output files.'


2012 continue   ! to exit the program
STOP
!---------------------
 contains


subroutine gnuplot_figures(path_sep, File_out_conv, File_out_average, File_out_conv_gnu, File_out_average_gnu, print_conv, tim_Start)
   character(1), intent(in) :: path_sep
   character(*), intent(in) :: File_out_conv, File_out_average, File_out_conv_gnu, File_out_average_gnu
   logical, intent(in) :: print_conv
   real(8), intent(in) :: tim_Start
   !----------------
   character(200) :: Gnu_script, Plot_file, ch_temp, ch_temp2, ch_temp3, ch_temp4
   integer :: FN_gnu_script


  !-------------------
  ! 1) Time evolution of convolved distribution:
  if (print_conv) then

   FN_gnu_script = 9996
   if (path_sep .EQ. '\') then	! if it is Windows
      Gnu_script = trim(adjustl(File_out_conv_gnu))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font "arial,14" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_conv_gnu))//'.gif'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron density (1/box)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      write(FN_gnu_script, '(a)') 'set format y "%2.0tx10^{%L}"'
   else
      Gnu_script = trim(adjustl(File_out_conv_gnu))//'.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_conv_gnu))//'.gif'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font \"arial,14\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron density (1/box)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      write(FN_gnu_script, '(a)') 'set format y "%2.0tx10^{%L}"'
   endif

   ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
   write(ch_temp,'(f)') 100.0d0      ! Scell(NSC)%E_top
   write(ch_temp2,'(f)') abs(tim_Start)
   if (tim_Start > 0.0d0) then
      ch_temp2 = '+'//trim(adjustl(ch_temp2))
   else
      ch_temp2 = '-'//trim(adjustl(ch_temp2))
   endif
   write(ch_temp3,'(f)') 1.0


   ! minimal energy grid:
   write(ch_temp4,'(f)') -25.0d0  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)
   if (path_sep .EQ. '\') then	! if it is Windows
      write(FN_gnu_script, '(a)') 'stats "'//trim(adjustl(File_out_conv))//'" nooutput'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') 'do for [i=1:int(STATS_blocks)] {'
      write(FN_gnu_script, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] "'// &
         trim(adjustl(File_out_conv))// &
         '" index (i-1) u 1:2 w l lw 3 title sprintf("%i fs",(i-1'// &
         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
   else  ! Linux
      write(FN_gnu_script, '(a)') 'stats \"'//trim(adjustl(File_out_conv))//'\" nooutput'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') 'do for [i=1:int(STATS_blocks)] {'
      write(FN_gnu_script, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] \"'// &
         trim(adjustl(File_out_conv))// &
         '\" index (i-1) u 1:2 w l lw 3 title sprintf(\"%i fs\",(i-1'// &
         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
   endif
   write(FN_gnu_script, '(a)') '}'

   if (path_sep .EQ. '\') then	! if it is Windows
      ! nothing to do
   else
      write(FN_gnu_script, '(a)') 'reset'
      write(FN_gnu_script, '(a)') '" | gnuplot '
   endif

   ! Done, clean up:
   close (FN_gnu_script)
  endif



   !--------------------
   ! 2) Average distribution:
   FN_gnu_script = 9995
   if (path_sep .EQ. '\') then   ! if it is Windows
      Gnu_script = trim(adjustl(File_out_average_gnu))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font "arial,14" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_average_gnu))//'.png'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron density (1/box)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      write(FN_gnu_script, '(a)') '#set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'p [][] "'//trim(adjustl(File_out_average))//'" u 1:2 w l lw LW title "Average"'
   else
      Gnu_script = trim(adjustl(File_out_average_gnu))//'.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_average_gnu))//'.png'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font \"arial,14\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron density (1/box)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      write(FN_gnu_script, '(a)') '#set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'p [][] \"'//trim(adjustl(File_out_average))//'\" u 1:3 w l lw \"$LW\" title \"Average\"'
      write(FN_gnu_script, '(a)') 'reset'
      write(FN_gnu_script, '(a)') '" | gnuplot '
   endif

   ! Done, clean up:
   close (FN_gnu_script)
end subroutine gnuplot_figures



subroutine read_weigths_file(FN_weights, weights_read, weights_COL)
   integer, intent(in) :: FN_weights      ! file number
   real(8), dimension(:), allocatable, intent(inout) :: weights_read  ! weights
   integer, intent(inout) :: weights_COL  ! column index
   !-------------------------
   logical :: file_opened
   integer :: i, Nsiz
   real(8) :: temp_r(weights_COL)
   character(1000) :: temp_ch

   call Count_lines_in_file(FN_weights, Nsiz)  ! below
   allocate(weights_read(Nsiz))

   do i = 1, Nsiz
      read(FN_weights,*,IOSTAT=Reason) temp_ch
      if (trim(adjustl(temp_ch(1:1))) == '#') then
         ! skip comment line
      else
         backspace(FN_weights)   ! to read the line again
         read(FN_weights,*,IOSTAT=Reason) temp_r
         weights_read(i) = temp_r(weights_COL)  ! number of this column is the weight
         if (Reason /= 0) then ! wrong format
            print*, 'Could not read line', i, 'in file provided file with weights'
            print*, 'Averaging without weights will be performed'
            weights_COL = -1  ! cannot be used, ignor weights then
            exit
         endif
      endif
   enddo

   inquire(unit = FN_weights, opened=file_opened)
   if (file_opened) close(FN_weights)
end subroutine read_weigths_file



subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    character(100) temp_ch
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=604) temp_ch
       enddo
       604 continue
    endif
    i = 0
    do
        read(File_num,*, end=603) temp_ch
        if (trim(adjustl(temp_ch(1:1))) /= '#') i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines_in_file



subroutine average_distr(distr_aver, distr_in, weights_read, weights_COL, i)
   real(8), dimension(:), intent(in) :: distr_in   ! input distribution function
   real(8), dimension(:), intent(inout) :: distr_aver ! averaged distribution function
   real(8), dimension(:), intent(in) :: weights_read  ! weights read from file
   integer, intent(in) :: weights_COL, i  ! column index; time index
   !------------------------------
   if (weights_COL > 0) then ! weights are provided
      distr_aver(:) = distr_aver(:) + distr_in(:) * weights_read(i)
   else ! no weights, just average:
      distr_aver(:) = distr_aver(:) + distr_in(:)
   endif
end subroutine average_distr



subroutine print_convolved(FN_out_conv, Ei, distr_conv)
   integer, intent(in) :: FN_out_conv  ! file to print into (must be opened)
   real(8), dimension(:), intent(in) :: Ei, distr_conv ! convolved function
   !-----------------------
   integer :: i, Nsiz
   Nsiz = size(distr_conv)

   do i = 1, Nsiz
      write(FN_out_conv,*) Ei(i), distr_conv(i)
   enddo

   ! Demarcation between blocks:
   write(FN_out_conv,*)   ! first empty line
   write(FN_out_conv,*)   ! second empty line
end subroutine print_convolved


subroutine print_averaged(FN_out_average, Ei, distr_aver)
   integer, intent(in) :: FN_out_average  ! file to print into (must be opened)
   real(8), dimension(:), intent(in) :: Ei, distr_aver ! convolved function
   !-----------------------
   integer :: i, Nsiz
   Nsiz = size(distr_aver)

   do i = 1, Nsiz
      write(FN_out_average,*) Ei(i), distr_aver(i)
   enddo
end subroutine print_averaged



subroutine convolve_with_Gaussian(Ei, distr, Conv_dE, distr_conv)
   real(8), dimension(:), intent(in) :: Ei, distr ! energy grid, distribution on this grid
   real(8), intent(in) :: Conv_dE   ! Gaussian width to convolve with [eV]
   real(8), dimension(:), intent(out) :: distr_conv ! convolved function
   !--------------------
   integer :: i, j, Nsiz
   real(8) :: Gaus_fact
   Nsiz = size(Ei)
   distr_conv = 0.0d0   ! to start with
   do i = 1, Nsiz ! for all points on the grid
      CONVD:do j = 1, Nsiz ! convolution with all the points
         Gaus_fact = Gaussian(Ei(j), Conv_dE, Ei(i)) ! below
         if ((Ei(j) > Ei(i)) .and. (Gaus_fact < 1.0d-10)) exit CONVD
         distr_conv(i) = distr_conv(i) + distr(j) * Gaus_fact
      enddo CONVD
   enddo
end subroutine convolve_with_Gaussian




pure function Gaussian(mu, sigma, x, normalized_max) result (Gaus)
   real(8) :: Gaus
   real(8), intent(in) :: mu, sigma, x ! position, width, variable
   real(8), intent(in), optional :: normalized_max
   real(8), parameter :: Pi = 3.1415926535897932384626433832795d0
   if (present(normalized_max)) then
      Gaus = exp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))*normalized_max ! it will be equal to "normalized_max" at the maximum
   else ! normalized as integral(Gaus) = 1
      Gaus = 1.0d0/(sqrt(2.0d0*Pi)*sigma)*exp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma)) ! it will be normalized to integral=1
   endif
end function Gaussian



subroutine read_EL_DISTR_txt(File_in, Conv_dE, File_weights, FN_weights, weights_COL)
   character(*), intent(in) :: File_in ! file name of the file with parameters for post-processing
   real(8), intent(out) :: Conv_dE  ! gaussian width for convlution in energy space
   character(*), intent(out) :: File_weights ! file name with weights for averaging of distribution over time
   integer, intent(in) :: FN_weights   ! file number with weights, to check if it exists and opened it
   integer, intent(out) :: weights_COL  ! column index in this file to use
   !-----------------------
   integer :: FN_in  ! file number of the file with parameters for post-processing
   real(8) :: temp_r
   integer :: temp_i, Reason
   character(100) :: temp_ch

   FN_in = 9997
   open (unit=FN_in, file=trim(adjustl(File_in)), status = 'old', readonly)

   ! Defaults to start with:
   Conv_dE = -1.0d0
   File_weights = '0.txt'
   weights_COL = -1
   ! Read the file:
   read(FN_in,*,IOSTAT=Reason) temp_r   ! energy convolution parameter
   if (Reason == 0) then   ! read well
      Conv_dE = temp_r
   else
      print*, 'Could not interprete line #1 in file ', trim(adjustl(File_in))
      print*, 'Convolution in energy space will not be performed'
   endif

   read(FN_in,*,IOSTAT=Reason) temp_ch, temp_i   ! file name and column index
   if (Reason == 0) then   ! read well
      File_weights = temp_ch
      weights_COL = temp_i
      open (unit=FN_weights, file=trim(adjustl(File_weights)), status = 'old', readonly)
   else
      print*, 'Could not interprete line #2 in file ', trim(adjustl(File_in))
      print*, 'Averaeging in time will be performed without weights'
   endif

   close(FN_in)
end subroutine read_EL_DISTR_txt


subroutine creat_output_files(FN_out_conv, File_out_conv, FN_out_average, File_out_average, print_conv)  ! below
   integer, intent(in) :: FN_out_conv, FN_out_average ! file numbers
   character(*), intent(in) :: File_out_conv, File_out_average ! file names
   logical :: print_conv   ! flag for convolved output
   !-----------------
   if (print_conv) then
      open (unit=FN_out_conv, file=trim(adjustl(File_out_conv)))
   endif
   open (unit=FN_out_average, file=trim(adjustl(File_out_average)))
end subroutine creat_output_files


subroutine get_size_of_distribution(FN_distr, Ei, distr)  ! below
   integer, intent(in) :: FN_distr ! file number of the file with electronic distribution
   real(8), dimension(:), allocatable, intent(inout) :: Ei
   real(8), dimension(:), allocatable, intent(inout) :: distr
   !-------------
   real(8) :: temp
   integer :: i, Reason
   logical :: read_well

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

   rewind(FN_distr) ! to restart reading into the arrays from the start
end subroutine get_size_of_distribution



subroutine read_distribution(FN_distr, Ei, distr, read_well)
   integer, intent(in) :: FN_distr ! file number of the file with electronic distribution
   real(8), dimension(:), allocatable, intent(inout) :: Ei
   real(8), dimension(:), allocatable, intent(inout) :: distr
   logical, intent(inout) :: read_well
   !-------------
   integer :: i, Reason

   Ei = 0.0d0  ! to start with
   distr = 0.0d0  ! to start with
   read_well = .true. ! to start with
   do i = 1, size(Ei)
      ! Read the file with distribution block by block:
      read(FN_distr,*,IOSTAT=Reason) Ei(i), distr(i)
      if (Reason < 0) then ! end of file
         read_well = .false.
         exit
      endif
   enddo
end subroutine read_distribution





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

END PROGRAM XTANT_el_distribution_analysis
