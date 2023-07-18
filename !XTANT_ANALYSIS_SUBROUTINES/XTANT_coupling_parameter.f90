PROGRAM Coupling_parameter
! Compilation:
!
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_coupling_parameter.f90 -o XTANT_coupling_parameter.exe /link /stack:9999999999
!
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays XTANT_coupling_parameter.f90 -o XTANT_coupling_parameter.exe /link /stack:9999999999
!
! To execute:
! XTANT_coupling_parameter.exe time Te_max
! where "time" is an optional argument meaning the time from which to start averaging the coupling parameter (used to exclude early time where atoms are not equilibrated yet)
!<===========================================
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
USE IFLPORT
use Universal_constants

character(200) :: File_electron, File_temperatures, File_coupling_partial, File_pressure, File_energies, File_Ce_partial
character(200) :: File_name1, File_name2, File_name3, File_name4, File_name5, File_name6, char_var
character(200) :: File_name_out, File_name_out2, File_name_out3, File_name_out4, File_name_out5, File_name_out6
character(200) :: File_numpar, File_matter, Matter_name, File_INPUT
character(200), dimension(:), allocatable :: Folders_with_data
character(200) :: gnu_script, gnu_fig
character(10) :: call_slash, sh_cmd
character(1) :: path_sep
real(8) :: starting_time, ending_time, scaling_G, Up_lim_Te
real(8), dimension(:,:,:), allocatable :: G_mean_part, Ce_mean_part
real(8), dimension(:,:), allocatable :: G_mean, mu_mean, Ce_mean
real(8), dimension(:,:), allocatable :: P_mean, E_mean, Grun_mean
real(8), dimension(:), allocatable :: Te_grid, T_ave, G_ave, G_err, Ce_err
real(8), dimension(:,:), allocatable :: G_part_ave, Ce_part_ave
real(8), dimension(:), allocatable :: mu_ave, Ce_ave, E_ave, P_ave, Grun_ave
integer :: FN1, FN2, FN3, FN4, FN5, FN6, FN7, FN8
integer :: FN_out, FN_out2, FN_out3, FN_out4, FN_out5, FN_out6  ! file number
integer :: Reason, i, j, siz, Tsiz, N_arg, i_arg, INFO
logical :: read_well, file_exist, file_exist2


!---------------------------------------
print*, '******************************************************************************'
print*, 'For analysis of coulping, call it with the options:'
print*, 'XTANT_coupling_parameter.exe t0 Te_max'
print*, 'where you set numbers: t0, Te_max'
print*, 't0 is the starting time [fs] (default t=t0-2*FWHM fs)'
print*, 'Te_max is the maximal electron temerature [K] (default Te_max=25000 K)'
print*, '******************************************************************************'

INFO = 0 ! to start with


call Path_separator(path_sep)  ! Objects_and_types

! Set defaults:
FN_out6 = 9987
FN_out5 = 9989
FN_out4 = 9993
FN_out3 = 9994
FN_out2 = 9996
FN_out  = 9997
FN1 = 9999
FN2 = 9998
FN3 = 9995
FN4 = 9992
FN5 = 9991
FN6 = 9990
FN7 = 9988
FN8 = 9986
File_electron = 'OUTPUT_electron_properties.dat'
File_temperatures = 'OUTPUT_temperatures.dat'
File_coupling_partial = 'OUTPUT_coupling.dat'
File_pressure = 'OUTPUT_pressure_and_stress.dat'
File_energies = 'OUTPUT_energies.dat'
File_Ce_partial = 'OUTPUT_electron_Ce.dat'
File_numpar = 'NUMERICAL_PARAMETERS'
File_matter = 'INPUT_MATERIAL'
File_INPUT = 'INPUT'

File_name_out =  'OUT_average_coupling.dat'
File_name_out2 = 'OUT_average_parameters.dat'
File_name_out3 = 'OUT_average_partial_couplings.dat'
File_name_out4 = 'OUT_average_pressure.dat'
File_name_out5 = 'OUT_average_partial_Ce.dat'

! Get the starting time, if user defined it:
! if (IARGC() >= 1) then 	! there was at least 1 argument passes by the user:
!    call getarg(1, char_var, status=Reason)	! read only the first one
!    if (Reason <= 0) then
!       starting_time = -1.0d20	! start from -infinity, i.e. include all data
!    else
!       read(char_var,*) starting_time
!    endif
! endif

! Default values (may be overwritten by the user):
starting_time = -1.0d20	! start from -infinity, i.e. include all data
Up_lim_Te = 25e3  ! uper limit for the grid on Te

! Count how many arguments the user provided:
N_arg = COMMAND_ARGUMENT_COUNT() ! Fortran intrinsic function
ALLARG:do i_arg = 1, N_arg ! read all the arguments passed
   ! Read the argument provided:
   call GET_COMMAND_ARGUMENT(i_arg, char_var)  ! intrinsic

   select case (i_arg)
   case (1) ! variable #1 must be starting_time
      read(char_var,*) starting_time
   case (2) ! uper limit tfor Te
      read(char_var,*) Up_lim_Te
   end select
enddo ALLARG


! print*, 'starting_time', starting_time

! Get all the output folders names:
call collect_all_output(Folders_with_data)	! below

! Read the factor for G rescaling set by the user:
call read_NUMPAR_file(Folders_with_data, path_sep, FN7, File_numpar, FN8, File_matter, File_INPUT, scaling_G, &
                        starting_time, Matter_name, INFO)  ! below
! If something went wrong, cannot execute the analysis:
if (INFO /= 0) goto 2012   ! go to the exit

File_name_out6 = 'OUT_XTANT3_'//trim(adjustl(Matter_name))//'_partial_Ce_G.dat'


! Create the electron temeprature grid:
call create_grid(Te_grid, G_mean, size(Folders_with_data), Up_lim_Te)
allocate(mu_mean( size(G_mean,1) , size(G_mean,2) ) )
allocate(Ce_mean( size(G_mean,1) , size(G_mean,2) ) )
allocate(P_mean( size(G_mean,1) , size(G_mean,2) ) )
allocate(E_mean( size(G_mean,1) , size(G_mean,2) ) )
allocate(Grun_mean( size(G_mean,1) , size(G_mean,2) ) )

! Get all the data and set them on the grid:
siz = size(Folders_with_data)
do i = 1, siz	! for all output data files
   File_name1 = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_electron))
   File_name2 = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_temperatures))
   File_name3 = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_coupling_partial))
   File_name4 = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_pressure))
   File_name5 = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_energies))
   File_name6 = trim(adjustl(Folders_with_data(i)))//path_sep//trim(adjustl(File_Ce_partial))
   open (unit=FN1, file=trim(adjustl(File_name1)), status = 'old', readonly) ! electron properties
   open (unit=FN2, file=trim(adjustl(File_name2)), status = 'old', readonly) ! temperatures
   open (unit=FN4, file=trim(adjustl(File_name4)), status = 'old', readonly) ! pressure
   open (unit=FN5, file=trim(adjustl(File_name5)), status = 'old', readonly) ! energies
   inquire(file=trim(adjustl(File_name3)),exist=file_exist)
   inquire(file=trim(adjustl(File_name6)),exist=file_exist2)
   if (file_exist) then ! also partial coupling data ! coupling
      open (unit=FN3, file=trim(adjustl(File_name3)), status = 'old', readonly)

      if (file_exist2) then ! is there particla Ce
         open (unit=FN6, file=trim(adjustl(File_name6)), status = 'old', readonly)
         call  read_and_set_on_grid(FN1, FN2, FN3, FN4, FN5, FN6, starting_time, Te_grid(:), G_mean(i,:), mu_mean(i,:), Ce_mean(i,:), &
               E_mean(i,:), P_mean(i,:), Grun_mean(i,:), G_mean_part=G_mean_part, &
               Ce_mean_part=Ce_mean_part, i_fold=i, siz=siz)     ! below
         close(FN6)
      else ! no partial Ce
         call  read_and_set_on_grid(FN1, FN2, FN3, FN4, FN5, FN6, starting_time, Te_grid(:), G_mean(i,:), mu_mean(i,:), Ce_mean(i,:), &
               E_mean(i,:), P_mean(i,:), Grun_mean(i,:), G_mean_part=G_mean_part, i_fold=i, siz=siz)     ! below
      endif
      close(FN3)
   else ! no partial coupling
       if (file_exist2) then ! is there particla Ce
         open (unit=FN6, file=trim(adjustl(File_name6)), status = 'old', readonly)
         call  read_and_set_on_grid(FN1, FN2, FN3, FN4, FN5, FN6, starting_time, Te_grid(:), G_mean(i,:), mu_mean(i,:), Ce_mean(i,:), &
               E_mean(i,:), P_mean(i,:), Grun_mean(i,:), Ce_mean_part=Ce_mean_part, i_fold=i, siz=siz)     ! below
         close(FN6)
      else ! no partial Ce
         call  read_and_set_on_grid(FN1, FN2, FN3, FN4, FN5, FN6, starting_time, Te_grid(:), G_mean(i,:), mu_mean(i,:), Ce_mean(i,:), &
               E_mean(i,:), P_mean(i,:), Grun_mean(i,:))     ! below
      endif
   endif
   
   close(FN1)
   close(FN2)
   close(FN4)
   close(FN5)
enddo

! Get the average results to print out:
Tsiz = size(Te_grid)
allocate(T_ave(Tsiz))
allocate(G_ave(Tsiz), mu_ave(Tsiz), Ce_ave(Tsiz))
allocate(E_ave(Tsiz), P_ave(Tsiz), Grun_ave(Tsiz))
allocate(G_err(Tsiz))
allocate(Ce_err(Tsiz))
if (allocated(G_mean_part)) then
   allocate(G_part_ave( Tsiz, size(G_mean_part,3) ) )
endif
if (allocated(Ce_mean_part)) then
   allocate(Ce_part_ave( Tsiz, size(Ce_mean_part,3) ) )
endif

do i = 1, Tsiz
   T_ave(i) = Te_grid(i)
   G_ave(i) = SUM(G_mean(:,i))/dble(siz)
   mu_ave(i) = SUM(mu_mean(:,i))/dble(siz)
   Ce_ave(i) = SUM(Ce_mean(:,i))/dble(siz)
   E_ave(i) = SUM(E_mean(:,i))/dble(siz)
   P_ave(i) = SUM(P_mean(:,i))/dble(siz)
   Grun_ave(i) = SUM(Grun_mean(:,i))/dble(siz)

   if (size(G_mean,1) > 1) then
      G_err(i) = sqrt(SUM( (G_mean(:,i) - G_ave(i))*(G_mean(:,i) - G_ave(i)) ) ) / (dble(size(G_mean,1))-1.0d0)
      Ce_err(i) = sqrt(SUM( (Ce_mean(:,i) - Ce_ave(i))*(Ce_mean(:,i) - Ce_ave(i)) ) ) / (dble(size(Ce_mean,1))-1.0d0)
   else
      G_err(i) = sqrt(SUM( (G_mean(:,i) - G_ave(i))*(G_mean(:,i) - G_ave(i)) ) )
      Ce_err(i) = sqrt(SUM( (Ce_mean(:,i) - Ce_ave(i))*(Ce_mean(:,i) - Ce_ave(i)) ) )
   endif
   if (allocated(G_mean_part)) then
      do j = 1, size(G_part_ave,2)
         G_part_ave(i,j) = SUM(G_mean_part(:,i,j))/dble(siz)
      enddo
   endif
   if (allocated(Ce_mean_part)) then
      do j = 1, size(Ce_part_ave,2)
         Ce_part_ave(i,j) = SUM(Ce_mean_part(:,i,j))/dble(siz)
      enddo
   endif
enddo
! Check if all the datapoints are present and printout:
open (unit=FN_out, file=trim(adjustl(File_name_out)))
open (unit=FN_out2, file=trim(adjustl(File_name_out2)))
if (allocated(G_mean_part)) then
   open (unit=FN_out3, file=trim(adjustl(File_name_out3)))
endif
if (allocated(Ce_mean_part)) then
   open (unit=FN_out5, file=trim(adjustl(File_name_out5)))
endif
open (unit=FN_out4, file=trim(adjustl(File_name_out4)))
open (unit=FN_out6, file=trim(adjustl(File_name_out6)))

write(FN_out,'(a)') 'Te    G    sigma(G)'
write(FN_out,'(a)') 'K    W/(m^3K)    -'

write(FN_out2,'(a)') 'Te    mu  Ce  sigma(Ce)'
write(FN_out2,'(a)') 'K    eV    J/(m^3K) -'

write(FN_out3,'(a)') 'Te    G_total    G_atoms-  G_shell_pairs-'
write(FN_out3,'(a)') 'K    W/(m^3K) '

write(FN_out4,'(a)') 'Te    E  P Gruneisen'
write(FN_out4,'(a)') 'K    eV/atom    GPa Pa/(J/atom)'

write(FN_out6,'(a)') '# Te    Ce_total Ce_s  Ce_p  Ce_d   G_total G_s-s G_s-p G_s-d G_p-p G_p-d G_d-d'
write(FN_out6,'(a)') '# K    J/(m^3K)  -- W/(m^3K) --'

!oooooooooooooooooooooooooooooooooo
! OUT files:
do i =1,Tsiz

   write(FN_out,'(es,es,es)') T_ave(i), G_ave(i), G_err(i)

   write(FN_out2,'(es,es,es,es)') T_ave(i), mu_ave(i), Ce_ave(i), Ce_err(i)

   if (allocated(G_mean_part)) then
      write(FN_out3,'(es)', advance='no') T_ave(i)
       do j = 1, size(G_part_ave,2)
          write(FN_out3,'(es)', advance='no') G_part_ave(i,j)
       enddo
      write(FN_out3,'(a)') ''
   endif

   write(FN_out4,'(es,es,es,es)') T_ave(i), E_ave(i), P_ave(i), Grun_ave(i)*(1d9/g_e)

   if (allocated(Ce_mean_part)) then
      write(FN_out5,'(es)', advance='no') T_ave(i)
       do j = 1, size(Ce_part_ave,2)
          write(FN_out5,'(es)', advance='no') Ce_part_ave(i,j)
       enddo
      write(FN_out5,'(a)') ''
   endif

   ! The output for partial Ce and G (rescaled!) in one file:
   if (allocated(Ce_mean_part)) then
      write(FN_out6,'(es)', advance='no') T_ave(i)
      do j = 1, size(Ce_part_ave,2)
         write(FN_out6,'(es)', advance='no') Ce_part_ave(i,j)
      enddo
      ! Depending on the basis size, we have different number of columns:
      ! UNFINISHED: only a single atom and a few common basis sets are accounted for!!!
      if (size(G_part_ave,2) == 3) then ! s: total, atom-atom, s-s
         write(FN_out6,'(es)', advance='no') G_part_ave(i,1) * scaling_G   ! total
      elseif (size(G_part_ave,2) == 6) then ! sp3: total, atom-atom, s-s, s-p, p-s, pp
         write(FN_out6,'(es)', advance='no') G_part_ave(i,1) * scaling_G  ! total
         write(FN_out6,'(es)', advance='no') G_part_ave(i,3) * scaling_G  ! s-s
         write(FN_out6,'(es)', advance='no') (G_part_ave(i,4)+G_part_ave(i,5)) * scaling_G   ! s-p + p-s
         write(FN_out6,'(es)', advance='no') G_part_ave(i,6) * scaling_G  ! p-p
      elseif (size(G_part_ave,2) == 11 ) then ! sp3d5: total, atom-atom, s-s, s-p, s-d, p-s, p-p, p-d, d-s, d-p, d-d
                                                      !  1        2       3    4    5    6    7    8    9    10   11
         write(FN_out6,'(es)', advance='no') G_part_ave(i,1) * scaling_G  ! total
         write(FN_out6,'(es)', advance='no') G_part_ave(i,3) * scaling_G  ! s-s
         write(FN_out6,'(es)', advance='no') (G_part_ave(i,4)+G_part_ave(i,6)) * scaling_G  ! s-p + p-s
         write(FN_out6,'(es)', advance='no') (G_part_ave(i,5)+G_part_ave(i,9)) * scaling_G  ! s-d + p-d
         write(FN_out6,'(es)', advance='no') G_part_ave(i,7) * scaling_G  ! p-p
         write(FN_out6,'(es)', advance='no') (G_part_ave(i,8)+G_part_ave(i,10)) * scaling_G ! p-d + d-p
         write(FN_out6,'(es)', advance='no') G_part_ave(i,11) * scaling_G ! d-d
      endif
      write(FN_out6,'(a)') ''
   endif

enddo


! Gnuplot of Ce and G:
if (path_sep .EQ. '\') then	! if it is Windows
   call_slash = 'call '
   sh_cmd = '.cmd'
else ! It is linux
   call_slash = './'
   sh_cmd = '.sh'
endif



call gnu_plot(File_name_out6, sh_cmd, 'OUT_XTANT3_'//trim(adjustl(Matter_name))//'_partial_Ce.gif', 'OUT_XTANT3_'//trim(adjustl(Matter_name))//'_partial_G.gif', G_part_ave, Ce_part_ave, Up_lim_Te)  ! below


close (FN_out)
close (FN_out2)
close (FN_out4)
if (allocated(G_mean_part)) close (FN_out3)
if (allocated(Ce_mean_part)) close (FN_out5)

2012  print*, 'XTANT: analysis with Coupling_parameter is executed'


!---------------------
 contains


subroutine gnu_plot(File_name_out, sh_cmd, File_Ce_gif, File_G_gif, G_part_ave, Ce_part_ave, Up_lim_Te)
   character(*), intent(in) :: File_name_out, sh_cmd, File_Ce_gif, File_G_gif
   real(8), dimension(:,:), intent(in) :: G_part_ave, Ce_part_ave
   real(8), intent(in) :: Up_lim_Te
   !-------------------
   character(300) :: File_name, File_name2, command, char_Te
   integer :: FN, iret
   real(8) :: x_tics

   if (Up_lim_Te < 10000.0e0) then
      x_tics = 1000.0d0  ! [K]
   elseif (Up_lim_Te < 50000.0e0) then
      x_tics = 5000.0d0  ! [K]
   else
      x_tics = 10000.0d0  ! [K]
   endif

   write(char_Te, '(i0)') CEILING(Up_lim_Te/1000.0e0)*1000

   ! Script for Ce:
   File_name  = 'OUT_gnuplot_Ce'//trim(adjustl(sh_cmd))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   call write_gnuplot_script_header_new(FN, 3, 3.0d0, x_tics,  'Ce', 'Electron temperature (K)', 'Electron heat capacity (J/(m^3K))',  trim(adjustl(File_Ce_gif)), path_sep, 2)

   if (path_sep .EQ. '\') then	! if it is Windows
      if (size(Ce_part_ave,2) < 3) then ! s
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), ' "u 1:2 w l lw LW title "Total" '
      else ! more than one line
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), ' "u 1:2 w l lw LW title "Total" ,\'
      endif

!       print*, 'Ce:', size(Ce_part_ave,2)

      if (size(Ce_part_ave,2) == 3) then ! sp3
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:3 w l lw LW title "s" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:4 w l lw LW title "p"'
      elseif (size(Ce_part_ave,2) == 4) then ! sp3d5
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:3 w l lw LW title "s" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:4 w l lw LW title "p" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:5 w l lw LW title "d" '
      endif
   else
      if (size(Ce_part_ave,2) < 3) then ! s
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] \"' , trim(adjustl(File_name_out)), '\"u 1:4 w l lw \"$LW\" title \"Total\" '
      else
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] \"' , trim(adjustl(File_name_out)), '\"u 1:4 w l lw \"$LW\" title \"Total\" ,\'
      endif

      if (size(Ce_part_ave,2) == 3) then ! sp3
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:3 w l lw \"$LW\" title "s" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:4 w l lw \"$LW\" title "p"'
      elseif (size(Ce_part_ave,2) == 4) then ! sp3d5
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:3 w l lw \"$LW\" title "s" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:4 w l lw \"$LW\" title "p" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:5 w l lw \"$LW\" title "d" '
      endif
   endif
   call write_gnuplot_script_ending_new(FN, File_name, path_sep)
   close(FN)


   ! Coupling parameter:
   File_name2  = 'OUT_gnuplot_G'//trim(adjustl(sh_cmd))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name2)), action="write", status="replace")

   call write_gnuplot_script_header_new(FN, 3, 3.0d0, x_tics,  'G', 'Electron temperature (K)', 'Electron-phonon coupling (W/(m^3K))',  trim(adjustl(File_G_gif)), path_sep, 0)

   if (path_sep .EQ. '\') then	! if it is Windows
      if (size(G_part_ave,2) < 4) then ! s
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), ' "u 1:3 w l lw LW title "Total" '
      elseif (size(G_part_ave,2) == 6) then ! sp3
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), ' "u 1:5 w l lw LW title "Total" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:6 w l lw LW title "s-s" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:7 w l lw LW title "s-p" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:8 w l lw LW title "p-p"'
      elseif (size(G_part_ave,2) > 6) then ! sp3d5
         write(FN, '(a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), ' "u 1:6 w l lw LW title "Total" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:7 w l lw LW title "s-s" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:8 w l lw LW title "s-p" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:9 w l lw LW title "s-d" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:10 w l lw LW title "p-p" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:11 w l lw LW title "p-d" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(File_name_out)), ' "u 1:12 w l lw LW title "d-d" '
      endif
   else
      if (size(G_part_ave,2) < 4) then ! s
         write(FN, '(a,es25.16,a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), '\"u 1:3 w l lw \"$LW\" title \"Total\" '
      elseif (size(G_part_ave,2) == 6) then ! sp3
         write(FN, '(a,es25.16,a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), '\"u 1:5 w l lw \"$LW\" title \"Total\" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:6 w l lw \"$LW\" title "s-s" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:7 w l lw \"$LW\" title "s-p" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:8 w l lw \"$LW\" title "p-p"'
      elseif (size(G_part_ave,2) > 6) then ! sp3d5
         write(FN, '(a,es25.16,a,a,a)') 'p [0.0:'//trim(adjustl(char_Te))//'][] "' , trim(adjustl(File_name_out)), '\"u 1:6 w l lw \"$LW\" title \"Total\" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:7 w l lw \"$LW\" title "s-s" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:8 w l lw \"$LW\" title "s-p" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:9 w l lw \"$LW\" title "s-d" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:10 w l lw \"$LW\" title "p-p" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:11 w l lw \"$LW\" title "p-d" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(File_name_out)), '\"u 1:12 w l lw \"$LW\" title "d-d" '
      endif
   endif
   call write_gnuplot_script_ending_new(FN, File_name2, path_sep)
   close(FN)

   ! Execute gnuplot scripts to make figures:
   command = "."//trim(adjustl(path_sep))//trim(adjustl(File_name))
   iret = system(command)
   command = "."//trim(adjustl(path_sep))//trim(adjustl(File_name2))
   iret = system(command)
end subroutine gnu_plot



subroutine read_NUMPAR_file(Folders_with_data, path_sep, FN7, File_numpar, FN8, File_matter, File_INPUT, &
                            scaling_G, starting_time, Matter_name, INFO)
   character(*), dimension(:), intent(in) :: Folders_with_data
   character(1), intent(in) :: path_sep
   integer, intent(in) :: FN7, FN8
   character(*), intent(inout) :: File_numpar, File_matter, File_INPUT
   real(8), intent(inout) :: scaling_G, starting_time
   character(*), intent(out) :: Matter_name
   integer, intent(inout) :: INFO   ! error flag
   !------------------
   character(300) :: File_name, File_name2, string
   character(10) :: temp, temp2
   logical :: file_exists, short_named
   integer :: open_status, v, i, Reason
   real(8) :: r_temp, t0, t_FWHM

   scaling_G = 4.0d0 ! to start with, no known scaling
   short_named = .false.   ! to start with
   INFO = 0 ! to start with, no error

   ! Find filename for NUMPAR:
   file_exists = .false.
   v = 0
   do while (.not.file_exists)
      if (v > 0) then
         write(temp2, '(i6)') v
         write(temp, '(a)') '_'//trim(adjustl(temp2))//'.txt'
      else
         write(temp, '(a)') '.txt'
      endif
      File_name = trim(adjustl(Folders_with_data(1)))//path_sep//trim(adjustl(File_numpar))//trim(adjustl(temp))
      inquire(file=trim(adjustl(File_name)),exist=file_exists)

      ! Check if the long-named file exists:
      if (file_exists) then   ! two files with input
         short_named = .false.
         File_name2 = trim(adjustl(Folders_with_data(1)))//path_sep//trim(adjustl(File_matter))//trim(adjustl(temp))
      else  ! check one file with short name
         short_named = .true.
         File_name2 = trim(adjustl(Folders_with_data(1)))//path_sep//trim(adjustl(File_INPUT))//trim(adjustl(temp))
         inquire(file=trim(adjustl(File_name2)),exist=file_exists)
      endif
      !print*, v, trim(adjustl(File_name)), file_exists

      ! If not:
      v = v + 1   ! to check the next number
      if (v > 100) then
         print*, 'Could not find file ', trim(adjustl(File_name)), ', cannot use scaling factor!'
         exit
      endif
   enddo
   ! Once correct file found, save the name:
   if (.not.short_named) then
      File_numpar = trim(adjustl(File_name))
      ! Open file:
      open(UNIT=FN7, file=trim(adjustl(File_numpar)), iostat=open_status, action='read')
   endif

   File_matter = trim(adjustl(File_name2))
   inquire(file=trim(adjustl(File_matter)),exist=file_exists)
   if (.not.file_exists) then
      print*, 'File not found: ', trim(adjustl(File_matter))
      INFO = -1   ! file not found
      goto 2015
   endif
   open(UNIT=FN8, file=trim(adjustl(File_matter)), iostat=open_status, action='read')
   !print*, 'Matter:', trim(adjustl(File_matter)), file_exists, open_status

   ! Read material name:
   read(FN8,*) Matter_name
   do i = 2, 12   ! read some more to get the  pulse start
      read(FN8,*) string
      if (i == 10) then ! line for pulse duration:
         read(string, *) t_FWHM
      endif
      if (i == 12) then ! line for pulse center:
         read(string, *) t0
      endif
   enddo
   if (starting_time < -1.0d10) then ! user did not provide the starting time, use automatic estimation
      starting_time = t0 - 2.0d0*t_FWHM
   endif
   !print*, scaling_G, trim(adjustl(Matter_name)), starting_time
   print*, 'Scaling factor used: ', scaling_G
   print*, 'Staring time used: ', starting_time

   ! Read the file with numpar to find the sacling factor:
   if (.not.short_named) then
      do i = 1, 16
         read(FN7,*) ! skip first 16 lines
      enddo
      ! read line with scaling factor:
      read(FN7,*) r_temp, scaling_G
      close(FN7)
      close(FN8)
   else ! short-named inupt file:
      RDFL:do
         read(FN8,'(a)',IOSTAT=Reason) string
         if (Reason < 0) then ! end of file reached ...
            print*, 'Could not find Numerical parameters in File:', trim(adjustl(File_matter)), ', cannot use scaling factor!'
            close(FN8)
            goto 2015
         else
            select case (trim(adjustl(string)))
            case ('NUMPAR', 'NumPar', 'numpar', &
               'NUMERICAL_PARAMETERS', 'Numerical_Parameters', 'numerical_parameters', &
               'NUMERICS', 'Numerics', 'numerics')
               do i = 1, 16
                  read(FN8,*) ! skip first 16 lines
               enddo
               exit RDFL
            endselect
         endif
      enddo RDFL
      ! read line with scaling factor:
      read(FN8,*) r_temp, scaling_G
      close(FN8)
   endif

   ! Use the default value to renormalize it:
2015   scaling_G = 4.0d0/scaling_G ! to multiply G with this factor
end subroutine read_NUMPAR_file



subroutine read_and_set_on_grid(FN1, FN2, FN3, FN4, FN5, FN6, starting_time, Te_grid, G_mean, mu_mean, Ce_mean, &
            E_mean, P_mean, Grun_mean, G_mean_part, Ce_mean_part, i_fold, siz)
   ! file numbers (must be already opened):
   integer, intent(in) :: FN1 ! electron properties
   integer, intent(in) :: FN2 ! temperatures
   integer, intent(in) :: FN3 ! coupling
   integer, intent(in) :: FN4 ! pressure
   integer, intent(in) :: FN5 ! energies
   integer, intent(in) :: FN6 ! electron Ce
   real(8), intent(in) :: starting_time	! to cut off the times that needs to be excluded (during thermalization)
   real(8), dimension(:), intent(inout) :: Te_grid, G_mean, mu_mean, Ce_mean
   real(8), dimension(:), intent(inout) :: E_mean, P_mean, Grun_mean
   real(8), dimension(:,:,:), allocatable, intent(inout), optional :: G_mean_part, Ce_mean_part
   integer, intent(in), optional :: i_fold, siz
   !000000000000000000000000000
   integer :: i, Nsiz, i_num, j, j_low, j_high, j_l, j_h, Ncol, k
   real(8) :: Te(2), G(6), T_temp, E_file_read(9), P_file_read(2)
   real(8) :: G_temp, G_low, G_high
   real(8) :: mu_temp, mu_low, mu_high
   real(8) :: Ce_temp, Ce_low, Ce_high
   real(8) :: E_temp, E_low, E_high
   real(8) :: P_temp, P_low, P_high
   real(8) :: Grun_temp, Grun_low, Grun_high
   real(8), dimension(:), allocatable :: G_part_temp, G_part_low, G_part_high
   real(8), dimension(:), allocatable :: Ce_part_temp, Ce_part_low, Ce_part_high
   real(8), dimension(:), allocatable :: Te_read, G_read, tim, mu_read, Ce_read
   real(8), dimension(:), allocatable :: E_read, P_read, Grun_read
   real(8), dimension(:,:), allocatable :: G_part_read, Ce_part_read
   ! Read all the data:
   call Count_lines_in_file(FN1, Nsiz, skip_lines=2)
   allocate(Te_read(Nsiz), G_read(Nsiz), tim(Nsiz), mu_read(Nsiz), Ce_read(Nsiz))
   allocate(E_read(Nsiz), P_read(Nsiz), Grun_read(Nsiz))
   if (present(G_mean_part)) then
      call Count_columns_in_file(FN3, Ncol, skip_lines=1)
      allocate(G_part_read(Nsiz+1, Ncol))
      allocate(G_part_temp(Ncol))
      allocate(G_part_low(Ncol))
      allocate(G_part_high(Ncol))
      if (.not.allocated(G_mean_part)) allocate( G_mean_part(siz, size(G_mean), Ncol-1) )
      read(FN3,*,IOSTAT=Reason)
   endif
   if (present(Ce_mean_part)) then
      call Count_columns_in_file(FN6, Ncol, skip_lines=1)
      allocate(Ce_part_read(Nsiz+1, Ncol))
      allocate(Ce_part_temp(Ncol))
      allocate(Ce_part_low(Ncol))
      allocate(Ce_part_high(Ncol))
      if (.not.allocated(Ce_mean_part)) allocate( Ce_mean_part(siz, size(Ce_mean), Ncol-1) )
      read(FN3,*,IOSTAT=Reason)
   endif
   
   ! Skip first lines with comments:
   read(FN1,*,IOSTAT=Reason)
   read(FN1,*,IOSTAT=Reason)
   read(FN2,*,IOSTAT=Reason)
   read(FN2,*,IOSTAT=Reason)
   read(FN4,*,IOSTAT=Reason)
   read(FN4,*,IOSTAT=Reason)
   read(FN5,*,IOSTAT=Reason)
   read(FN5,*,IOSTAT=Reason)
   
   do i = 1, Nsiz
      read(FN1,*,IOSTAT=Reason) G(:)	! coupling parameter
      IF (Reason .LT. 0) exit
      read(FN2,*,IOSTAT=Reason) Te(:)	! electron temeprature
      IF (Reason .LT. 0) exit
      if (present(G_mean_part)) then
         read(FN3,*,IOSTAT=Reason) G_part_read(i,:)
         IF (Reason .LT. 0) exit
!         print*, i, G_part_read(i,2:4)
      endif
      read(FN4,*,IOSTAT=Reason) P_file_read(:)	! pressures
      read(FN5,*,IOSTAT=Reason) E_file_read(:)	! energies
      if (present(Ce_mean_part)) then
         read(FN6,*,IOSTAT=Reason) Ce_part_read(i,:)
         IF (Reason .LT. 0) exit
!         print*, i, G_part_read(i,2:4)
      endif

      ! Now extract the variables to be used:
      tim(i) = Te(1)	! time
      ! By construction of the XTANT output file, the columns are:
      Te_read(i) = Te(2)
      G_read(i) = G(6)
      mu_read(i) = G(3)
      Ce_read(i) = G(5)
      E_read(i) = E_file_read(8)
      P_read(i) = P_file_read(2)
      if (i /= 1) then
        if ( abs(E_read(i) - E_read(i-1)) > 1.0d-6 ) then
           Grun_read(i) = (P_read(i) - P_read(i-1))/(E_read(i) - E_read(i-1))
         else
            Grun_read(i) = 0.0d0
         endif
      else
         Grun_read(i) = 0.0d0
      endif

!       print*, i, G_read(i), G_part_read(i,1:4)
   enddo
   
   ! Sort the data by electronic temperature:
   if (present(G_mean_part)) then
      if (present(Ce_mean_part)) then
         call sort_array(Te_read, G_read, mu_read, Ce_read, E_read, P_read, Grun_read, Ar2d=G_part_read, Ar2d_2=Ce_part_read)   ! below
      else
         call sort_array(Te_read, G_read, mu_read, Ce_read, E_read, P_read, Grun_read, Ar2d=G_part_read)   ! below
      endif
   else
      if (present(Ce_mean_part)) then
         call sort_array(Te_read, G_read, mu_read, Ce_read, E_read, P_read, Grun_read, Ar2d_2=Ce_part_read)   ! below
      else
         call sort_array(Te_read, G_read, mu_read, Ce_read, E_read, P_read, Grun_read)   ! below
      endif
   endif
   
   ! Average the data onto the grid:
   j = 1	! to start with
   do i = 1, size(Te_grid)
      G_temp = 0.0d0	! to start with
      mu_temp = 0.0d0	! to start with
      Ce_temp = 0.0d0	! to start with
      E_temp = 0.0d0	! to start with
      P_temp = 0.0d0	! to start with
      Grun_temp = 0.0d0	! to start with
      if (present(G_mean_part)) G_part_temp = 0.0d0 ! to start with
      if (present(Ce_mean_part)) Ce_part_temp = 0.0d0 ! to start with
      i_num = 0
      AV:do
         j = j + 1
         if (j > size(Te_read)) exit AV
         if (Te_read(j) <= Te_grid(i)) then
            i_num = i_num + 1
            ! Add it up only if it is within the allowed time:
            if (tim(j) > starting_time) then
               G_temp = G_temp + G_read(j)
               mu_temp = mu_temp + mu_read(j)
               Ce_temp = Ce_temp + Ce_read(j)
               E_temp = E_temp + E_read(j)
               P_temp = P_temp + P_read(j)
               Grun_temp = Grun_temp + Grun_read(j)
               if (present(G_mean_part)) then
                  G_part_temp(:) = G_part_temp(:) + G_part_read(j,:)
               endif
               if (present(Ce_mean_part)) then
                  Ce_part_temp(:) = Ce_part_temp(:) + Ce_part_read(j,:)
               endif
!                print*, j, G_read(j), G_part_read(j,1:4)
!                pause
            else
!                print*, 'The value is excluded:', j, tim(j), starting_time
            endif
         else
            j = j - 1
            exit AV
         endif
      enddo AV
      if (i_num > 0) then
         G_mean(i) = G_temp/dble(i_num)
         mu_mean(i) = mu_temp/dble(i_num)
         Ce_mean(i) = Ce_temp/dble(i_num)
         E_mean(i) = E_temp/dble(i_num)
         P_mean(i) = P_temp/dble(i_num)
         Grun_mean(i) = Grun_temp/dble(i_num)
         if (present(G_mean_part)) then
            do k = 1, size(G_mean_part,3)
               G_mean_part(i_fold,i,k) = G_part_temp(k+1)/dble(i_num) 
            enddo
         endif
         if (present(Ce_mean_part)) then
            do k = 1, size(Ce_mean_part,3)
               Ce_mean_part(i_fold,i,k) = Ce_part_temp(k+1)/dble(i_num)
            enddo
         endif
      else
         G_low = 0.0d0
         mu_low = 0.0d0
         Ce_low = 0.0d0
         E_low = 0.0d0
         P_low = 0.0d0
         Grun_low = 0.0d0
         if (present(G_mean_part)) G_part_low = 0.0d0 ! to start with
         if (present(Ce_mean_part)) Ce_part_low = 0.0d0 ! to start with
         j_l = 1
         do j_low = min(j,size(Te_read)), 1, -1
            G_low = G_read(j_low)
            mu_low = mu_read(j_low)
            Ce_low = Ce_read(j_low)
            E_low = E_read(j_low)
            P_low = P_read(j_low)
            Grun_low = Grun_read(j_low)
            if (present(G_mean_part)) G_part_low(:) = G_part_read(j_low,:)
            if (present(Ce_mean_part)) Ce_part_low(:) = Ce_part_read(j_low,:)
            if (G_low > 0) then
               j_l = j_low
               exit
            endif
         enddo
         G_high = 0.0d0
         mu_high = 0.0d0
         Ce_high = 0.0d0
         E_high = 0.0d0
         P_high = 0.0d0
         Grun_high = 0.0d0
         if (present(G_mean_part)) G_part_high = 0.0d0 ! to start with
         if (present(Ce_mean_part)) Ce_part_high = 0.0d0 ! to start with
         j_h = size(Te_read)
         do j_high = j+1, size(Te_read)
            G_high = G_read(j_high)
            mu_high = mu_read(j_high)
            Ce_high = Ce_read(j_high)
            E_high = E_read(j_high)
            P_high = P_read(j_high)
            Grun_high = Grun_read(j_high)
            if (present(G_mean_part)) G_part_high(:) = G_part_read(j_high,:)
            if (present(Ce_mean_part)) Ce_part_high(:) = Ce_part_read(j_high,:)
            !if (G_low > 0) then
            if (G_high > 0) then
               j_h = j_high
               exit
            endif
         enddo
         ! Add it up only if it is within the allowed time:
         if (tim(j_l) > starting_time) then   
            if (j_h == j_l) then
               G_mean(i) = G_read(j_h)
               mu_mean(i) = mu_read(j_h)
               Ce_mean(i) = Ce_read(j_h)
               E_mean(i) = E_read(j_h)
               P_mean(i) = P_read(j_h)
               Grun_mean(i) = Grun_read(j_h)
               if (present(G_mean_part)) then
                  do k = 1, size(G_mean_part,3)
                     G_mean_part(i_fold,i,k) = G_part_read(j_h,k+1)
                  enddo
               endif
               if (present(Ce_mean_part)) then
                  do k = 1, size(Ce_mean_part,3)
                     Ce_mean_part(i_fold,i,k) = Ce_part_read(j_h,k+1)
                  enddo
               endif
            else
               T_temp = (Te_grid(i) - Te_read(j_l)) / (Te_read(j_h) - Te_read(j_l))
               !G_mean(i) = G_low + (G_high - G_low)/(Te_read(j_h) - Te_read(j_l))*(Te_grid(i) - Te_read(j_l))
               G_mean(i) = G_low + (G_high - G_low) * T_temp 
               mu_mean(i) = mu_low + (mu_high - mu_low) * T_temp 
               Ce_mean(i) = Ce_low + (Ce_high - Ce_low) * T_temp
               E_mean(i) = E_low + (E_high - E_low) * T_temp
               P_mean(i) = P_low + (P_high - P_low) * T_temp
               Grun_mean(i) = Grun_low + (Grun_high - Grun_low) * T_temp
               if (present(G_mean_part)) then
                  do k = 1, size(G_mean_part,3)
                     G_mean_part(i_fold,i,k) = G_part_low(k+1) + (G_part_high(k+1) - G_part_low(k+1)) * T_temp
                  enddo
               endif
               if (present(Ce_mean_part)) then
                  do k = 1, size(Ce_mean_part,3)
                     Ce_mean_part(i_fold,i,k) = Ce_part_low(k+1) + (Ce_part_high(k+1) - Ce_part_low(k+1)) * T_temp
                  enddo
               endif
            endif
         else
!             print*, 'The value is excluded2:', j_l, tim(j_l), starting_time
         endif
      endif
   enddo
   
   deallocate(Te_read, G_read, mu_read, Ce_read, E_read, P_read, Grun_read, tim)
   if (present(G_mean_part)) then
      deallocate(G_part_low, G_part_temp, G_part_high, G_part_read)
   endif
   if (present(Ce_mean_part)) then
      deallocate(Ce_part_low, Ce_part_temp, Ce_part_high, Ce_part_read)
   endif
end subroutine read_and_set_on_grid


subroutine sort_array(Ev, Ar2, Ar3, Ar4, Ar5, Ar6, Ar7, Ar2d, Ar2d_2)	! bubble method
   real(8), dimension(:), intent(inout) :: Ev	! array to sort
   real(8), dimension(:), intent(inout) :: Ar2
   real(8), dimension(:), intent(inout), optional :: Ar3, Ar4, Ar5, Ar6, Ar7
   real(8), dimension(:,:), intent(inout), optional :: Ar2d
   real(8), dimension(:,:), intent(inout), optional :: Ar2d_2
   real(8) :: temp_c, temp_c2, temp_c3, temp_c4, temp_c5, temp_c6, temp_c7
   integer N,i,j,k
   logical :: swapped
   N = size(Ev)
   do j = N-1, 1, -1
      do i = 1, j
         IF (Ev(i) > Ev(i+1)) THEN
            temp_c = Ev(i)
            Ev(i) = Ev(i+1)
            Ev(i+1) = temp_c
            ! And the associated arrays too:
            temp_c2 = Ar2(i)
            Ar2(i) = Ar2(i+1)
            Ar2(i+1) = temp_c2
            if (present(Ar3)) then
               temp_c3 = Ar3(i)
               Ar3(i) = Ar3(i+1)
               Ar3(i+1) = temp_c3
            endif
            if (present(Ar4)) then
               temp_c4 = Ar4(i)
               Ar4(i) = Ar4(i+1)
               Ar4(i+1) = temp_c4
            endif
            if (present(Ar5)) then
               temp_c5 = Ar5(i)
               Ar5(i) = Ar5(i+1)
               Ar5(i+1) = temp_c5
            endif
            if (present(Ar6)) then
               temp_c6 = Ar6(i)
               Ar6(i) = Ar6(i+1)
               Ar6(i+1) = temp_c6
            endif
            if (present(Ar7)) then
               temp_c7 = Ar7(i)
               Ar7(i) = Ar7(i+1)
               Ar7(i+1) = temp_c7
            endif
            if (present(Ar2d)) then
               do k = 1, size(Ar2d,2)
                  temp_c5 = Ar2d(i,k)
                  Ar2d(i,k) = Ar2d(i+1,k)
                  Ar2d(i+1,k) = temp_c5
               enddo
            endif
            if (present(Ar2d_2)) then
               do k = 1, size(Ar2d_2,2)
                  temp_c5 = Ar2d_2(i,k)
                  Ar2d_2(i,k) = Ar2d_2(i+1,k)
                  Ar2d_2(i+1,k) = temp_c5
               enddo
            endif

            swapped = .TRUE.
         END IF
      enddo
      IF (.NOT. swapped) EXIT
   enddo
end subroutine sort_array



subroutine create_grid(Te_grid, G_mean, N_data, Up_lim_in)
   real(8), dimension(:), allocatable, intent(inout) :: Te_grid
   real(8), dimension(:,:), allocatable, intent(inout) :: G_mean
   integer, intent(in) :: N_data
   real(8), intent(in), optional :: Up_lim_in
   real(8) :: Up_lim, step
   integer :: i, Nsiz
   
   if (allocated(Te_grid)) deallocate(Te_grid)
   if (allocated(G_mean)) deallocate(G_mean)
   
   if (present(Up_lim_in)) then  ! user-defined uer limit of Te
      Up_lim = Up_lim_in
   else  ! default limit of Te
      Up_lim = 25000.0d0
   endif
   step = 100.0d0
   Nsiz = int(Up_lim/step)
   allocate(Te_grid(Nsiz), source = 0.0d0)
   allocate(G_mean(N_data, Nsiz), source = 0.0d0)
   ! fill the grid with the values:
   do i = 1, Nsiz
      Te_grid(i) = dble(i)*step
   enddo
end subroutine create_grid
 
 

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
         read(FN,'(a)',IOSTAT=Reason) Folders_with_data(i)	! all folders with the outputs of XTANT
         if (Reason < 0) exit
      enddo
      close(FN, status='delete') ! temp file is not needed anymore, erase it
   endif TEMP
end subroutine collect_all_output

 

pure function linear_interpolation(x, x1, x0, y1, y0) result ( y )
   real(8), intent(in) :: x, x1, x0, y1, y0
   real(8) :: y
   y = y0 + (y1 - y0)/(x1 - x0)*(x - x0)
end function linear_interpolation


! Reads additional data from the command line passed along:
subroutine get_add_data(File_name, hw, read_well)
   character(*), intent(inout) :: File_name
   real(8), intent(inout) :: hw
   logical, intent(inout) :: read_well  ! did we read ok?
   !---------------------------------------
   integer, dimension(:), allocatable :: stop_markers ! how many different optionas are passed?
   integer :: i, count_dash
   character(500) :: string, char1
   read_well = .true. ! to start with
   
   string = '' ! to start with empty line
   ! Read all arguments passed to the code:
   do i = 1, iargc() ! for as many arguments as we passed
      call getarg(i, char1)
      string = trim(adjustl(string))//' '//trim(adjustl(char1)) ! collect them all into one line
   enddo
   
   read(string,*) hw	! start with reading the energy required [eV]
   
2017   if (.not.read_well) then
      write(*,'(a)') '***************************************************************************'
      write(*,'(a)') 'Could not interpret the passed string: '//trim(adjustl(string))
   endif
end subroutine get_add_data



subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then ! in case you want to skip some comment lines and count only lines with data
       do i=1,skip_lines
          read(File_num,*, end=604) 
       enddo
       604 continue
    endif
    i = 0
    do
        read(File_num,*, end=603)
        i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines_in_file



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


subroutine Count_columns_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of columns in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    real(8) temp
    character(10000) temp_ch
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
   case(1:)	! eps
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
         write(FN, '(a)') 'set terminal jpeg large font arial '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif large font arial'
         write(FN, '(a)') 'set output \"$NAME\"'
      case (4)  ! png
         write(FN, '(a)') 'set terminal png font arial '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color '
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
         write(FN, '(a)') 'set terminal postscript enhanced "Helvetica" 16 color '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (2)  ! gpeg
         write(FN, '(a)') 'set terminal jpeg large font "arial,16" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif large font "arial,16" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (4)  ! png
         write(FN, '(a)') 'set terminal png font "arial,16" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font "arial,16" '
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
   write(FN, '(a,a)') 'set xtics ', trim(adjustl(temp))//' font "arial,14"'
   write(FN, '(a)') 'set ytics font "arial,14"'
end subroutine write_gnuplot_script_header_windows_new


subroutine write_gnuplot_script_ending_new(FN, File_name, path_sep)
   integer, intent(in) :: FN
   character(*), intent(in) :: File_name
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   integer :: i

   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      write(FN, '(a)') 'reset'
      write(FN, '(a)') '" | gnuplot '
      !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
      i = system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
end subroutine write_gnuplot_script_ending_new



END PROGRAM Coupling_parameter
