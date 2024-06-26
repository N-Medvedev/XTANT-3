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


real(8), dimension(:), allocatable :: Ei, Ei_fe, Ei_fe_aver, Ei_gridded
real(8), dimension(:), allocatable :: fe, fe_gridded  ! distribution function on molecular orbitals
real(8), dimension(:), allocatable :: fe_aver, fe_gridded_aver ! averaged fe with given wieghts
real(8), dimension(:), allocatable :: distr, distr_norm  ! distribution function of the equidistant grid
real(8), dimension(:), allocatable :: distr_conv, distr_norm_conv
real(8), dimension(:), allocatable :: distr_aver, distr_norm_aver
real(8), dimension(:), allocatable :: weights_read
real(8) :: tim, tim_Start, Conv_dE, temp_r, temp_r_norm, exclude_limit, Grid_dE
integer :: FN_distr, FN_out_conv, FN_out_average, FN_weights, weights_COL
integer :: INFO, Reason, i, Nsiz, tim_counter
integer :: FN_fe, FN_out_fe_aver, FN_out_fe_gridded, FN_out_fe_gridded_aver
logical :: file_exist, read_well, file_opened, norm_present
logical :: print_conv, print_fe, print_on_grid
character(100) :: File_distr, File_out_conv, File_out_average, File_weights, File_in
character(100) :: File_fe, File_out_fe_average, File_out_fe_gridded, File_out_fe_gridded_average
character(100) :: File_out_conv_gnu, File_out_average_gnu
character(1) :: path_sep, temp_ch

call Path_separator(path_sep)  ! below

! File names:
File_in           = 'EL_DISTR.txt'  ! file with parameters to use for post-processing (see XTANT Manual)
File_fe           = 'OUTPUT_electron_distribution.dat'
File_out_fe_average  = 'OUT_fe_average'
File_out_fe_gridded  = 'OUT_fe_gridded'
File_out_fe_gridded_average  = 'OUT_fe_gridded_average'
File_distr        = 'OUTPUT_electron_distribution_on_grid.dat'
File_out_conv     = 'OUT_el_distr_vs_time.dat'
File_out_average  = 'OUT_el_distr_average.dat'
File_out_conv_gnu     = 'OUT_el_distr_vs_time'
File_out_average_gnu  = 'OUT_el_distr_average'


! File numbers:
FN_weights = 9998
FN_distr = 9999
FN_fe = 9994
FN_out_conv = 1000
FN_out_average = 2000
FN_out_fe_aver = 3000
FN_out_fe_gridded = 3001
FN_out_fe_gridded_aver = 3002

! Set defaults:
INFO = 0

!---------------------------------------
print*, '******************************************************************************'
print*, 'For analysis of distribution, call XTANT_el_distribution_analysis.exe'
print*, 'Files required by the program:'
print*, 'EL_DISTR.txt -- file with parameters (see XTANT manual)'
print*, 'OUTPUT_electron_distribution.dat -- file with fe(E)'
print*, 'OUTPUT_electron_distribution_on_grid.dat -- file with fe(E) on grid'
print*, '******************************************************************************'


! Check which distributions are present:
inquire(file=trim(adjustl(File_distr)),exist=print_on_grid) ! check if input file is there
inquire(file=trim(adjustl(File_fe)),exist=print_fe) ! check if input file is there


!-----------------------------------
! Read input:
inquire(file=trim(adjustl(File_in)),exist=file_exist) ! check if input file is there
if (file_exist) then ! try to read it, if there is a grid provided
   call read_EL_DISTR_txt(File_in, Conv_dE, Grid_dE, File_weights, FN_weights, weights_COL, exclude_limit) ! below
   ! Knowing the parameters, define flag to mark if convolution is required:
   print_conv = (Conv_dE > 0.0d0)   ! only if convolution is used
   ! Open output file(s):
   call creat_output_files(FN_out_conv, File_out_conv, FN_out_average, File_out_average, &
                           FN_out_fe_aver, trim(adjustl(File_out_fe_average))//'.dat', &
                           FN_out_fe_gridded, trim(adjustl(File_out_fe_gridded))//'.dat', &
                           FN_out_fe_gridded_aver, trim(adjustl(File_out_fe_gridded_average))//'.dat', &
                           print_conv, print_on_grid, print_fe)  ! below
   print*, 'Output files are created, starting the analysis'
else
   print*, 'File ', trim(adjustl(File_in)), 'not found, program terminates.'
   print*, 'For execution of the program, this file must be present,'
   print*, 'containing the following lines:'
   print*, 'Line #1: Energy width in [eV]. Set a negative number to exclude convolution in the energy space, if not needed.'
   print*, 'Line #2: File name, providing the data on the widths in time to average the distribution with;'
   print*, 'and the index of column to use in this file.'
   print*, 'Set a non-existing file name to exclude averaging in time, if not needed.'
   print*, 'Line #3: Cut off for plotting distribution: limit, below which data should be excluded from plotting.'
   print*, '******************************************************************************'
   goto 2012   ! exit, nothing else to do
endif

! Get the energy grid from the file:
if (print_on_grid) then
   open (unit=FN_distr, file=trim(adjustl(File_distr)), status = 'old', readonly)
endif
if (print_fe) then
   open (unit=FN_fe, file=trim(adjustl(File_fe)), status = 'old', readonly)
endif

! Get the size of the array for the distribution and allocate arrays:
call get_size_of_distribution(FN_distr, FN_fe, Ei, distr, Ei_fe, fe, Ei_gridded, Grid_dE, norm_present, print_fe, print_on_grid)  ! below
print*, 'Distribution size is known, proceeding to reading files'
if (print_fe) then
   print*, '1) Distribution on energy levels is provided to average'
else
   print*, '1) No distribution on energy levels provided'
endif
if (print_on_grid) then
   print*, '2) Distribution on grid is provided to average'
   if (.not.norm_present) then
      print*, 'Normalized distribution is NOT provided, only spectrum'
   else
      print*, 'Normalized distribution is provided, including it in the analysis'
   endif
else
   print*, '2) No distribution on grid provided'
endif

! Knowing the size, define convolved distribution array:
if (print_on_grid) then ! distribution on grid
   allocate(distr_conv(size(distr)), source = 0.0d0)
   allocate(distr_aver(size(distr)), source = 0.0d0)
   if (norm_present) then
      allocate(distr_norm(size(distr)), source = 0.0d0)
      allocate(distr_norm_conv(size(distr)), source = 0.0d0)
      allocate(distr_norm_aver(size(distr)), source = 0.0d0)
   endif
endif
if (print_fe) then   ! distribution on MO
   allocate(fe_aver(size(fe)), source = 0.0d0)
   allocate(Ei_fe_aver(size(fe)), source = 0.0d0)
   allocate(fe_gridded_aver(size(Ei_gridded)), source = 0.0d0)
   allocate(fe_gridded(size(Ei_gridded)), source = 0.0d0)
endif


! Read weights from the file:
if (weights_COL > 0) then
   call read_weigths_file(FN_weights, weights_read, weights_COL)  ! below
   print*, 'Weights are read from the file: ', trim(adjustl(File_weights))
endif

! Output files:
if (print_on_grid) then
   if (print_conv) open (unit=FN_out_conv, file=trim(adjustl(File_out_conv)))
   open (unit=FN_out_average, file=trim(adjustl(File_out_average)))
endif
if (print_fe) then
   open (unit=FN_out_fe_aver, file=trim(adjustl(File_out_fe_average))//'.dat')
   open (unit=FN_out_fe_gridded, file=trim(adjustl(File_out_fe_gridded))//'.dat')
   open (unit=FN_out_fe_gridded_aver, file=trim(adjustl(File_out_fe_gridded_average))//'.dat')
endif

i = 0 ! block counter
read_well = .true. ! to start with
do while (read_well)
   i = i + 1
   print*, 'Reading and analysing distribution, block #', i

   ! Read the file with distribution on grid block by block:
   if (print_on_grid) then
      read(FN_distr,*,IOSTAT=Reason) temp_ch, tim  ! read first line with timeprint
      if (Reason /= 0) then ! wrong format or end of file reached
         read_well = .false.
         exit
      else   ! normal reading
         read_well = .true.  ! it read well, nothing to report
      end if
      if (print_conv) write(FN_out_conv,*) temp_ch, tim  ! copy the same markerline (timestemp)
      if (i == 1) tim_Start = tim ! save starting time for gnuplotting below
   endif
   ! Distribution on MO:
   if (print_fe) then
      read(FN_fe,*,IOSTAT=Reason) temp_ch, tim  ! read first line with timeprint
      if (Reason /= 0) then ! wrong format or end of file reached
         read_well = .false.
         exit
      else   ! normal reading
         read_well = .true.  ! it read well, nothing to report
      end if
      if (i == 1) tim_Start = tim ! save starting time for gnuplotting below
   endif

   ! Read the distributions on grid at this timestep:
   if (print_on_grid) then
      call read_distribution(FN_distr, Ei, distr, distr_norm, norm_present, read_well)  ! below
      if (.not. read_well) exit
   endif
   ! Read the distributions on MO at this timestep:
   if (print_fe) then
      call read_distribution(FN_fe, Ei_fe, fe, distr_norm, .false., read_well)  ! below
      if (.not. read_well) exit
   endif

   ! Demarcation between blocks:
   if (print_on_grid) then
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
   endif
   if (print_fe) then
      read(FN_fe,*,IOSTAT=Reason)   ! skip empty line
      if (Reason /= 0) then ! wrong format or end of file reached
         read_well = .false.
         exit
      else   ! normal reading
         read_well = .true.  ! it read well, nothing to report
      end if
      read(FN_fe,*,IOSTAT=Reason)   ! skip second empty line
      if (Reason /= 0) then ! wrong format or end of file reached
         read_well = .false.
         exit
      else   ! normal reading
         read_well = .true.  ! it read well, nothing to report
      end if
   endif

   !-------------
   ! Now, process the data and print them out:
   ! If the convolution is required:
   if (print_on_grid) then
      if (print_conv) then
         ! Convolve electronic spetrum with Gaussian:
         call convolve_with_Gaussian(Ei, distr, Conv_dE, distr_conv) ! below
         !if (i == 1) temp_r = maxval(distr_conv)
         !distr_conv(:) = distr_conv(:) / (temp_r)

         ! Convolve normalized electronic distribution with Gaussian:
         if (norm_present) then
            call convolve_with_Gaussian(Ei, distr_norm, Conv_dE, distr_norm_conv) ! below
            ! And let's normalize it to the peak height:
            if (i == 1) then
               temp_r_norm = maxval(distr_norm_conv) * 0.5d0
               if (temp_r_norm <= 1.0d-14) temp_r_norm = 1.0d0   ! no normalization, if no data
            endif
            distr_norm_conv(:) = distr_norm_conv(:) / (temp_r_norm)
         endif

         ! Printout the convolved function:
         call print_convolved(FN_out_conv, Ei, distr_conv, distr_norm_conv, norm_present)   ! below

         ! Average convolved spectrum function over time:
         call average_distr(distr_aver, distr_conv, weights_read, weights_COL, i)   ! below
         ! Average convolved distribution function over time:
         if (norm_present) then
            call average_distr(distr_norm_aver, distr_norm_conv, weights_read, weights_COL, i)   ! below
         endif
      else
         ! Average original Spectrum over time:
         call average_distr(distr_aver, distr, weights_read, weights_COL, i)   ! below
         ! Average original Distribution function over time:
         if (norm_present) then
            call average_distr(distr_norm_aver, distr_norm, weights_read, weights_COL, i)   ! below
         endif
      endif
   endif ! (print_on_grid)

   if (print_fe) then
      ! Average original Spectrum over time:
      call average_distr(fe_aver, fe, weights_read, weights_COL, i)   ! below
      Ei_fe_aver(:) = Ei_fe_aver(:) + Ei_fe(:)*weights_read(i)
      call make_gridded_fe(Ei_gridded, fe_gridded_aver, fe_gridded, Ei_fe, fe, weights_read, weights_COL, i)
      ! Printout the convolved function:
      call print_convolved(FN_out_fe_gridded, Ei_gridded, fe_gridded, distr_norm_conv, .false.)   ! below
   endif ! (print_fe)

enddo
! Average distribution over time:
if (print_on_grid) then ! distribution on grid
   distr_aver(:) = distr_aver(:) / max(1.0d0, dble(i))
   ! And let's normalize it to the number of particles:
   temp_r = SUM(distr_aver)
   if (temp_r > 1.0d-14) distr_aver(:) = distr_aver(:) / temp_r

   if (norm_present) then
      distr_norm_aver(:) = distr_norm_aver(:) / max(1.0d0, dble(i))
      temp_r = maxval(distr_norm_aver)
      if (temp_r > 1.0d-14) distr_norm_aver(:) = distr_norm_aver(:) / (temp_r*0.5d0)
   endif
endif
if (print_fe) then   ! distribution on MO
   temp_r = SUM(weights_read)
   fe_aver(:) = fe_aver(:) / max(1.0d0, temp_r)
   Ei_fe_aver(:) = Ei_fe_aver(:) / max(1.0d0, temp_r)
   fe_gridded_aver(:) = fe_gridded_aver(:) / max(1.0d0, temp_r)
endif

! Clean up:
if (print_on_grid) then ! distribution on grid
   close(FN_distr)
endif
if (print_fe) then ! distribution on grid
   close(FN_fe)
endif

! Printout the avereaged distribution:
if (print_on_grid) then ! distribution on grid
   call print_averaged(FN_out_average, Ei, distr_aver, distr_norm_aver, norm_present, exclude_limit)   ! below
   close(FN_out_average)
   if (print_conv) close(FN_out_conv)
endif
if (print_fe) then ! distribution on grid
   call print_averaged(FN_out_fe_aver, Ei_fe_aver, fe_aver, distr_norm_aver, .false., exclude_limit)   ! below
   close(FN_out_fe_aver)
   close(FN_out_fe_gridded)
   call print_averaged(FN_out_fe_gridded_aver, Ei_gridded, fe_gridded_aver, distr_norm_aver, .false., exclude_limit)   ! below
   close(FN_out_fe_gridded_aver)
endif
print*, 'Analysis is done, starting gnuplotting'

!-------------------------
! Make gnuplot script:
if (print_on_grid) then ! distribution on grid
   call gnuplot_figures(path_sep, File_out_conv, File_out_average, File_out_conv_gnu, File_out_average_gnu, &
                     print_conv, tim_Start, norm_present)   ! below
endif
if (print_fe) then ! distribution on grid
   call gnuplot_figures_MO(path_sep, File_out_fe_average, exclude_limit)  ! below
   call gnuplot_figures_gridded(path_sep, File_out_fe_gridded, File_out_fe_gridded_average, exclude_limit)  ! below
endif
print*, 'Scripts are prepared, starting to gnuplot them...'

! execute gnuplot files:
call execute_gnuplots(path_sep, File_out_conv_gnu, File_out_average_gnu, norm_present, File_out_fe_average, &
                     File_out_fe_gridded, File_out_fe_gridded_average, print_on_grid, print_fe) ! below
print*, 'Everything is done, check the output files.'


2012 continue   ! to exit the program
STOP
!---------------------
 contains


subroutine execute_gnuplots(path_sep, File_out_conv_gnu, File_out_average_gnu, norm_present, File_out_fe_average, &
                  File_out_fe_gridded, File_out_fe_gridded_average, print_on_grid, print_fe)
   character(*), intent(in) :: path_sep, File_out_conv_gnu, File_out_average_gnu, File_out_fe_average, &
                           File_out_fe_gridded, File_out_fe_gridded_average
   logical, intent(in) :: norm_present, print_on_grid, print_fe
   !---------------------------
   character(100) :: command
   integer :: iret

   if (print_fe) then
      if (path_sep .EQ. '\') then	! if it is Windows
         !call system("OUTPUT_Gnuplot_all.cmd")
         command = trim(adjustl(File_out_fe_average))//'.cmd'
         iret = system(command)
         command = trim(adjustl(File_out_fe_gridded_average))//'.cmd'
         iret = system(command)
      else ! linux:
         !call system("./OUTPUT_Gnuplot_all.sh")
         command = "./"//trim(adjustl(File_out_fe_average))//'.sh'
         iret = system(command)
         command = "./"//trim(adjustl(File_out_fe_gridded_average))//'.sh'
         iret = system(command)
      endif
   endif

   if (print_on_grid) then
      if (path_sep .EQ. '\') then	! if it is Windows
         !call system("OUTPUT_Gnuplot_all.cmd")
         command = trim(adjustl(File_out_average_gnu))//'.cmd'
         iret = system(command)
         if (norm_present) then
            command = trim(adjustl(File_out_average_gnu))//'_norm.cmd'
            iret = system(command)
         endif
      else ! linux:
         !call system("./OUTPUT_Gnuplot_all.sh")
         command = "./"//trim(adjustl(File_out_average_gnu))//'.sh'
         iret = system(command)
         if (norm_present) then
            command = "./"//trim(adjustl(File_out_average_gnu))//'_norm.sh'
            iret = system(command)
         endif
      endif
   endif

   ! Time-dependent:
   if (print_fe) then
      if (path_sep .EQ. '\') then	! if it is Windows
         command = trim(adjustl(File_out_fe_gridded))//'.cmd'
         iret = system(command)
      else ! linux:
         command = "./"//trim(adjustl(File_out_fe_gridded))//'.sh'
         iret = system(command)
      endif
   endif

   if (print_on_grid) then
      if (path_sep .EQ. '\') then	! if it is Windows
         command = trim(adjustl(File_out_conv_gnu))//'.cmd'
         iret = system(command)
         if (norm_present) then
            command = trim(adjustl(File_out_conv_gnu))//'_norm.cmd'
            iret = system(command)
         endif
      else ! linux:
         command = "./"//trim(adjustl(File_out_conv_gnu))//'.sh'
         iret = system(command)
         if (norm_present) then
            command = "./"//trim(adjustl(File_out_conv_gnu))//'_norm.sh'
            iret = system(command)
         endif
      endif
   endif

end subroutine execute_gnuplots





subroutine gnuplot_figures(path_sep, File_out_conv, File_out_average, File_out_conv_gnu, File_out_average_gnu, &
            print_conv, tim_Start, norm_present)
   character(1), intent(in) :: path_sep
   character(*), intent(in) :: File_out_conv, File_out_average, File_out_conv_gnu, File_out_average_gnu
   logical, intent(in) :: print_conv, norm_present
   real(8), intent(in) :: tim_Start
   !----------------
   character(200) :: Gnu_script, Plot_file, ch_temp, ch_temp2, ch_temp3, ch_temp4
   integer :: FN_gnu_script


  !-------------------
  ! 1) Time evolution of convolved distribution and spectrum:
if (print_conv) then

   ! 1.a) Electronic spectrum:
   FN_gnu_script = 9996
   if (path_sep .EQ. '\') then	! if it is Windows
      Gnu_script = trim(adjustl(File_out_conv_gnu))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font "arial,16" '
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
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font \"arial,16\" '
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


   ! 1.b) Electronic distribution:
  if (norm_present) then
   FN_gnu_script = 9996
   if (path_sep .EQ. '\') then	! if it is Windows
      Gnu_script = trim(adjustl(File_out_conv_gnu))//'_norm.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font "arial,16" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_conv_gnu))//'_norm.gif'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron distribution'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      !write(FN_gnu_script, '(a)') 'set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{%L}"'
   else
      Gnu_script = trim(adjustl(File_out_conv_gnu))//'_norm.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_conv_gnu))//'_norm.gif'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font \"arial,16\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron distribution'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      !write(FN_gnu_script, '(a)') 'set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{\%L}"'
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
         '" index (i-1) u 1:3 w l lw 3 title sprintf("%i fs",(i-1'// &
         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
   else  ! Linux
      write(FN_gnu_script, '(a)') 'stats \"'//trim(adjustl(File_out_conv))//'\" nooutput'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') 'do for [i=1:int(STATS_blocks)] {'
      write(FN_gnu_script, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] \"'// &
         trim(adjustl(File_out_conv))// &
         '\" index (i-1) u 1:3 w l lw 3 title sprintf(\"%i fs\",(i-1'// &
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
endif



   !--------------------
   ! 2) Average distribution:

   ! 2.a) Electron spectrum:
   FN_gnu_script = 9995
   if (path_sep .EQ. '\') then   ! if it is Windows
      Gnu_script = trim(adjustl(File_out_average_gnu))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font "arial,16" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_average_gnu))//'.png'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron density (1/eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      write(FN_gnu_script, '(a)') '#set logscale y'
      !write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') '#set format y "10^{%L}"'
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
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font \"arial,16\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron density (1/eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      write(FN_gnu_script, '(a)') '#set logscale y'
      !write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') '#set format y "10^{\%L}"'
      write(FN_gnu_script, '(a)') 'p [][] \"'//trim(adjustl(File_out_average))//'\" u 1:2 w l lw \"$LW\" title \"Average\"'
      write(FN_gnu_script, '(a)') 'reset'
      write(FN_gnu_script, '(a)') '" | gnuplot '
   endif
   ! Done, clean up:
   close (FN_gnu_script)


   ! 2.b) Electron distribution:
  if (norm_present) then
   FN_gnu_script = 9995
   if (path_sep .EQ. '\') then   ! if it is Windows
      Gnu_script = trim(adjustl(File_out_average_gnu))//'_norm.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font "arial,16" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_average_gnu))//'_norm.png'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron distribution'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      write(FN_gnu_script, '(a)') '#set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'p [][] "'//trim(adjustl(File_out_average))//'" u 1:3 w l lw LW title "Average"'
   else
      Gnu_script = trim(adjustl(File_out_average_gnu))//'_norm.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_average_gnu))//'_norm.png'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font \"arial,16\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron distribution'//'\" font \"arial,18\" '
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
  endif
end subroutine gnuplot_figures




subroutine gnuplot_figures_MO(path_sep, File_out_fe_average, exclude_limit)
   character(1), intent(in) :: path_sep
   character(*), intent(in) :: File_out_fe_average
   real(8), intent(in) :: exclude_limit
   !----------------
   character(200) :: Gnu_script
   integer :: FN_gnu_script

   !--------------------
   ! Average distribution:
   FN_gnu_script = 9995
   if (path_sep .EQ. '\') then   ! if it is Windows
      Gnu_script = trim(adjustl(File_out_fe_average))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font "arial,16" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_fe_average))//'.png'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron distribution'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{%L}"'
      write(FN_gnu_script, '(a)') 'p [][] "'//trim(adjustl(File_out_fe_average))//'.dat'//'" u 1:2 w l lw LW title "Average"'
   else
      Gnu_script = trim(adjustl(File_out_fe_average))//'.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_fe_average))//'.png'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font \"arial,16\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron distribution'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "\%2.0tx10^{\%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{\%L}"'
      write(FN_gnu_script, '(a)') 'p [][] \"'//trim(adjustl(File_out_fe_average))//'.dat'//'\" u 1:2 w l lw \"$LW\" title \"Average\"'
      write(FN_gnu_script, '(a)') 'reset'
      write(FN_gnu_script, '(a)') '" | gnuplot '
   endif
   ! Done, clean up:
   close (FN_gnu_script)
end subroutine gnuplot_figures_MO


subroutine gnuplot_figures_gridded(path_sep, File_out_fe_gridded, File_out_fe_gridded_average, exclude_limit)  ! below
   character(1), intent(in) :: path_sep
   character(*), intent(in) :: File_out_fe_gridded, File_out_fe_gridded_average
   real(8), intent(in) :: exclude_limit
   !----------------
   character(200) :: Gnu_script, ch_temp, ch_temp2, ch_temp3, ch_temp4
   integer :: FN_gnu_script


   FN_gnu_script = 9996
   if (path_sep .EQ. '\') then	! if it is Windows
      Gnu_script = trim(adjustl(File_out_fe_gridded))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font "arial,16" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_fe_gridded))//'.gif'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron distribution'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      !write(FN_gnu_script, '(a)') 'set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{%L}"'
   else
      Gnu_script = trim(adjustl(File_out_fe_gridded))//'.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_fe_gridded))//'.gif'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal gif animate delay 10 font \"arial,16\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron distribution'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      !write(FN_gnu_script, '(a)') 'set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{\%L}"'
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
      write(FN_gnu_script, '(a)') 'stats "'//trim(adjustl(File_out_fe_gridded))//'.dat'//'" nooutput'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') 'do for [i=1:int(STATS_blocks)] {'
      write(FN_gnu_script, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] "'// &
         trim(adjustl(File_out_fe_gridded))//'.dat'// &
         '" index (i-1) u 1:2 w l lw 3 title sprintf("%i fs",(i-1'// &
         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
   else  ! Linux
      write(FN_gnu_script, '(a)') 'stats \"'//trim(adjustl(File_out_fe_gridded))//'.dat'//'\" nooutput'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') 'do for [i=1:int(STATS_blocks)] {'
      write(FN_gnu_script, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] \"'// &
         trim(adjustl(File_out_fe_gridded))//'.dat'// &
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


   !--------------------
   ! Average distribution:
   FN_gnu_script = 9995
   if (path_sep .EQ. '\') then   ! if it is Windows
      Gnu_script = trim(adjustl(File_out_fe_gridded_average))//'.cmd'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      write(FN_gnu_script, '(a,f3.1)') 'LW=', 3.0
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font "arial,16" '
      write(FN_gnu_script, '(a)') 'set output "'//trim(adjustl(File_out_fe_gridded_average))//'.png'//'"'
      write(FN_gnu_script, '(a)') 'set xlabel "'//'Energy (eV)'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set ylabel "'//'Electron distribution'//'" font "arial,18"'
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics 10'
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "%2.0tx10^{%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{%L}"'
      write(FN_gnu_script, '(a)') 'p [][] "'//trim(adjustl(File_out_fe_gridded_average))//'.dat'//'" u 1:2 w l lw LW title "Average"'
   else
      Gnu_script = trim(adjustl(File_out_fe_gridded_average))//'.sh'
      open (unit=FN_gnu_script, file=trim(adjustl(Gnu_script )))
      write(FN_gnu_script, '(a)') '#!/bin/bash'
      write(FN_gnu_script, '(a)') ''
      write(FN_gnu_script, '(a)') 'NAME='//trim(adjustl(File_out_fe_gridded_average))//'.png'
      write(FN_gnu_script, '(a)') 'LABL="Distribution"'
      write(FN_gnu_script, '(a)') 'TICSIZ=10.00'
      write(FN_gnu_script, '(a)') 'echo "'
      write(FN_gnu_script, '(a)') 'set terminal pngcairo font \"arial,16\" '
      write(FN_gnu_script, '(a)') 'set output \"$NAME\"'
      write(FN_gnu_script, '(a)') 'set xlabel \"'//'Energy (eV)'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set ylabel \"'//'Electron distribution'//'\" font \"arial,18\" '
      write(FN_gnu_script, '(a)') 'set key right top '
      write(FN_gnu_script, '(a)') 'set xtics \"$TICSIZ\" '
      write(FN_gnu_script, '(a)') 'set logscale y'
      write(FN_gnu_script, '(a)') '#set format y "\%2.0tx10^{\%L}"'
      write(FN_gnu_script, '(a)') 'set format y "10^{\%L}"'
      write(FN_gnu_script, '(a)') 'p [][] \"'//trim(adjustl(File_out_fe_gridded_average))//'.dat'//'\" u 1:2 w l lw \"$LW\" title \"Average\"'
      write(FN_gnu_script, '(a)') 'reset'
      write(FN_gnu_script, '(a)') '" | gnuplot '
   endif
   ! Done, clean up:
   close (FN_gnu_script)

end subroutine gnuplot_figures_gridded



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


subroutine make_gridded_fe(Ei_gridded, fe_gridded_aver, fe_gridded, Ei_fe, fe, weights_read, weights_COL, i_weights)
   real(8), dimension(:), intent(inout) :: fe_gridded_aver, fe_gridded
   real(8), dimension(:), intent(in) :: Ei_gridded, Ei_fe, fe
   real(8), dimension(:), intent(in) :: weights_read  ! weights read from file
   integer, intent(in) :: weights_COL, i_weights  ! column index; time index
   !--------------------------
   integer :: i, Nsiz, count_i, j, NgrdSiz

   NgrdSiz = size(Ei_gridded)
   Nsiz = size(Ei_fe)
   fe_gridded = 0.0d0   ! to start with
   ! All points in distribution
   i = 1 ! to start with
   GRD:do j = 1, NgrdSiz
      count_i = 0 ! to start with
      do while (Ei_fe(i) < Ei_gridded(j))
         i = i + 1
         if (i > Nsiz) exit GRD
         count_i = count_i + 1   ! count how many points are within this grid step
         fe_gridded(j) = fe_gridded(j) + fe(i)
      enddo
      fe_gridded(j) = fe_gridded(j) / dble(max(count_i,1))
   enddo GRD

   !------------------------------
   ! Average gridded fe:
   if (weights_COL > 0) then ! weights are provided
      fe_gridded_aver(:) = fe_gridded_aver(:) + fe_gridded(:) * weights_read(i_weights)
   else ! no weights, just average:
      fe_gridded_aver(:) = fe_gridded_aver(:) + fe_gridded(:)
   endif
end subroutine make_gridded_fe



subroutine print_convolved(FN_out_conv, Ei, distr_conv, distr_norm_conv, norm_present)
   integer, intent(in) :: FN_out_conv  ! file to print into (must be opened)
   real(8), dimension(:), intent(in) :: Ei, distr_conv, distr_norm_conv ! convolved function
   logical, intent(in) :: norm_present
   !-----------------------
   integer :: i, Nsiz
   Nsiz = size(distr_conv)

   if (norm_present) then ! new format
      do i = 1, Nsiz
         write(FN_out_conv,*) Ei(i), distr_conv(i), distr_norm_conv(i)
      enddo
   else  ! legacy format
      do i = 1, Nsiz
         write(FN_out_conv,*) Ei(i), distr_conv(i)
      enddo
   endif
   ! Demarcation between blocks:
   write(FN_out_conv,*)   ! first empty line
   write(FN_out_conv,*)   ! second empty line
end subroutine print_convolved


subroutine print_averaged(FN_out_average, Ei, distr_aver, distr_norm_aver, norm_present, exclude_limit)
   integer, intent(in) :: FN_out_average  ! file to print into (must be opened)
   real(8), dimension(:), intent(in) :: Ei, distr_aver, distr_norm_aver ! convolved function
   logical, intent(in) :: norm_present
   real(8), intent(in) :: exclude_limit   ! limit, below which the data should be excluded
   !-----------------------
   integer :: i, Nsiz
   Nsiz = size(distr_aver)

   if (norm_present) then
      do i = 1, Nsiz
         if (distr_aver(i) > exclude_limit) then
            write(FN_out_average,*) Ei(i), distr_aver(i), distr_norm_aver(i)
         endif
      enddo
   else ! legacy format
      do i = 1, Nsiz
         if (distr_aver(i) > exclude_limit) then
            write(FN_out_average,*) Ei(i), distr_aver(i)
         endif
      enddo
   endif
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



subroutine read_EL_DISTR_txt(File_in, Conv_dE, Grid_dE, File_weights, FN_weights, weights_COL, exclude_limit)
   character(*), intent(in) :: File_in ! file name of the file with parameters for post-processing
   real(8), intent(out) :: Conv_dE, Grid_dE ! gaussian width for convlution in energy space; grid to set fe on
   character(*), intent(out) :: File_weights ! file name with weights for averaging of distribution over time
   integer, intent(in) :: FN_weights   ! file number with weights, to check if it exists and opened it
   integer, intent(out) :: weights_COL  ! column index in this file to use
   real(8), intent(out) :: exclude_limit  ! limit, to exclude data points below it (if provided)
   !-----------------------
   integer :: FN_in  ! file number of the file with parameters for post-processing
   real(8) :: temp_r, temp_r2
   integer :: temp_i, Reason
   character(100) :: temp_ch

   FN_in = 9997
   open (unit=FN_in, file=trim(adjustl(File_in)), status = 'old', readonly)

   ! Defaults to start with:
   Conv_dE = -1.0d0
   Grid_dE = 0.1d0   ! default grid step
   File_weights = '0.txt'
   weights_COL = -1
   ! Read the file:
   read(FN_in,*,IOSTAT=Reason) temp_r, temp_r2  ! energy convolution parameter
   if (Reason == 0) then   ! read well
      Conv_dE = temp_r
      Grid_dE = temp_r2
   else
      print*, 'Could not interprete line #1 in file ', trim(adjustl(File_in))
      print*, 'Convolution in energy space will not be performed'
   endif

   read(FN_in,*,IOSTAT=Reason) temp_ch, temp_i   ! file name and column index
   if (Reason == 0) then   ! read well
      File_weights = temp_ch
      inquire(file=trim(adjustl(File_weights)),exist=file_exist) ! check if input file is there
      if (file_exist) then ! file available, use it:
         weights_COL = temp_i
         open (unit=FN_weights, file=trim(adjustl(File_weights)), status = 'old', readonly)
      else  ! no file available => no averaging weights
         weights_COL = -1
      endif
   else
      print*, 'Could not interprete line #2 in file ', trim(adjustl(File_in))
      print*, 'Averaeging in time will be performed without weights'
   endif

   read(FN_in,*,IOSTAT=Reason) temp_r
   if (Reason == 0) then   ! read well
      exclude_limit = temp_r
   else  ! didn't read well, or not provided
      exclude_limit = -1.0d10 ! don't exclude anything
   endif

   close(FN_in)
end subroutine read_EL_DISTR_txt


subroutine creat_output_files(FN_out_conv, File_out_conv, FN_out_average, File_out_average, &
                              FN_out_fe_aver, File_out_fe_average, &
                              FN_out_fe_gridded, File_out_fe_gridded, &
                              FN_out_fe_gridded_aver, File_out_fe_gridded_average, &
                              print_conv, print_on_grid, print_fe)  ! below
   integer, intent(in) :: FN_out_conv, FN_out_average, FN_out_fe_aver, FN_out_fe_gridded, FN_out_fe_gridded_aver ! file numbers
   character(*), intent(in) :: File_out_conv, File_out_average, File_out_fe_average, &
                               File_out_fe_gridded, File_out_fe_gridded_average ! file names
   logical :: print_conv, print_on_grid, print_fe   ! flag for convolved output
   !-----------------
   if (print_on_grid) then
      if (print_conv) then
         open (unit=FN_out_conv, file=trim(adjustl(File_out_conv)))
      endif
      open (unit=FN_out_average, file=trim(adjustl(File_out_average)))
   endif
   if (print_fe) then
      open (unit=FN_out_fe_aver, file=trim(adjustl(File_out_fe_average)))
      open (unit=FN_out_fe_gridded, file=trim(adjustl(File_out_fe_gridded)))
      open (unit=FN_out_fe_gridded_aver, file=trim(adjustl(File_out_fe_gridded_average)))
   endif
end subroutine creat_output_files


subroutine get_size_of_distribution(FN_distr, FN_fe, Ei, distr, Ei_fe, fe, Ei_gridded, Grid_dE, norm_present, print_fe, print_on_grid)  ! below
   integer, intent(in) :: FN_distr, FN_fe ! file number of the file with electronic distribution
   real(8), dimension(:), allocatable, intent(inout) :: Ei, Ei_fe
   real(8), dimension(:), allocatable, intent(inout) :: distr, fe, Ei_gridded
   real(8), intent(in) :: Grid_dE   ! user-defined grid step size
   logical, intent(inout) :: norm_present, print_fe, print_on_grid
   ! in the old format, there was no third column with normalized distribution
   !-------------
   real(8) :: temp, Ei_Start, Ei_End
   integer :: i, Reason, Ncol, Nsiz
   logical :: read_well

   ! Find out if it is legacy format or new:
   if (print_on_grid) then
      call Count_columns_in_file(FN_distr, Ncol, skip_lines = 1)  ! below
      if (Ncol > 2) then
         norm_present = .true.
      else
         norm_present = .false.
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

      rewind(FN_distr) ! to restart reading into the arrays from the start
   endif ! print_on_grid

   ! Distribution on the electronic energy levels (molecular orbitals)
   if (print_fe) then
      ! Skip first line:
      read(FN_fe,*,IOSTAT=Reason)

      i = 0
      read_well = .true. ! to start with
      do while (read_well)
         i = i + 1
         ! Read the file with distribution block by block:
         read(FN_fe,*,IOSTAT=Reason) temp
         if (i == 1) Ei_Start = temp   ! save it for using below to set gridded distribution
         if (Reason /= 0) then ! wrong format: probably new block
            exit
         else
            Ei_End = temp  ! save it for using below to set gridded distribution
         endif
      enddo
      i = i - 1

      allocate(Ei_fe(i))
      allocate(fe(i))

      Nsiz = INT( (min(80.0e0,Ei_End) - Ei_Start)/max(Grid_dE,0.01e0) )
      allocate(Ei_gridded(Nsiz))
      ! Set the grid for gridded fe:
      Ei_gridded(1) = dble(FLOOR(Ei_Start))
      do i = 2, Nsiz
         Ei_gridded(i) = Ei_gridded(i-1) + max(Grid_dE,0.01e0)
      enddo

      rewind(FN_fe) ! to restart reading into the arrays from the start
   endif ! print_on_grid

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




subroutine read_distribution(FN_distr, Ei, distr, distr_norm, norm_present, read_well)
   integer, intent(in) :: FN_distr ! file number of the file with electronic distribution
   real(8), dimension(:), intent(inout) :: Ei
   real(8), dimension(:), intent(inout) :: distr, distr_norm
   logical, intent(in) :: norm_present ! new or legacy format
   logical, intent(inout) :: read_well
   !-------------
   integer :: i, Reason

   Ei = 0.0d0  ! to start with
   distr = 0.0d0  ! to start with
   read_well = .true. ! to start with
   if (norm_present) then
      distr_norm = 0.0d0  ! to start with
      do i = 1, size(Ei)
         ! Read the file with distribution block by block:
         read(FN_distr,*,IOSTAT=Reason) Ei(i), distr(i), distr_norm(i)
         if (Reason < 0) then ! end of file
            read_well = .false.
            exit
         endif
      enddo
   else  ! legacy format, no normalized Distribution
      do i = 1, size(Ei)
         ! Read the file with distribution block by block:
         read(FN_distr,*,IOSTAT=Reason) Ei(i), distr(i)
         if (Reason < 0) then ! end of file
            read_well = .false.
            exit
         endif
      enddo
   endif

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
