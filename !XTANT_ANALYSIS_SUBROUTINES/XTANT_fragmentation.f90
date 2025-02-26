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

include 'Periodic_table.f90'     ! include periodic Periodic_Table


PROGRAM Fragmentation
! Compilation:
!
! for DEBUG:
! ifort.exe -c /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_fragmentation.f90 /link /stack:9999999999
!
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs XTANT_fragmentation.obj -o XTANT_fragmentation.exe /link /stack:9999999999
!
!
! for RELEASE:
! ifort.exe -c /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays XTANT_fragmentation.f90 /link /stack:9999999999
!
! ifort.exe /F9999999999 /O3 /Qipo /fpp /Qopenmp /heap-arrays *.obj -o XTANT_fragmentation.exe /link /stack:9999999999
!
! To execute:
! XTANT_fragmentation.exe r dt
! here r is the cut-off radius in [A]
! dt is time-step for printout analyzed fragments [fs]
! Default values are: r=5 A, dt=100 fs
!<===========================================


use Universal_constants
use periodic_table
USE IFLPORT, only : system, chdir

type Instant_data
   real(8) :: Tim   ! time instant [fs]
   integer, dimension(:), allocatable :: KOA ! kind of atom, index
   character(3), dimension(:), allocatable :: name    ! element name
   real(8), dimension(:,:), allocatable :: R    ! absolute coordinates for all atoms [A]
   real(8), dimension(:,:), allocatable :: S    ! relative coordinates for all atoms
   real(8), dimension(3,3) :: Supce
   real(8), dimension(:), allocatable :: M    ! Atomic mass [a.m.u.]
   real(8), dimension(:), allocatable :: q    ! Atomic charge [electron charge]
end type Instant_data

!==============================================
type Displacement_analysis
   character(100) :: mask_name   ! name of the mask
   real(8) :: MSD_power          ! power of the mean displacement for this particular analysis
   real(8) :: mean_disp          ! [A^MSD_power] mean displacement of all atoms
   real(8), dimension(:), allocatable :: mean_disp_sort     ! [A^MSD_power] mean displacement of atoms by sort
   logical, dimension(3) :: print_r       ! along which exis the user requirested the analysis
   real(8), dimension(3) :: mean_disp_r   ! [A^MSD_power] mean displacements along X, Y, Z axes
   real(8), dimension(:,:), allocatable :: mean_disp_r_sort   ! [A^MSD_power] mean displacements along X, Y, Z axes by sort
   logical, dimension(:), allocatable :: Atomic_mask  ! atomic mask to be used
   ! Definition of spatial-section mask:
   ! Currently, only 2 sections are supported, connected by a logical operator 'and' or 'or':
   logical :: logical_and, logical_or  ! which one is used
   integer, dimension(2) :: axis_ind   ! index of the axis for the section
   real(8), dimension(2,3) :: r_start, r_end   ! section starting and ending points
end type Displacement_analysis

character(400) :: error_message, File_XYZ, File_supce, File_atomic_masks, read_line
character(200) :: File_fragments, Gnu_script, Gnu_file, command, File_m_over_z, chtemp, File_name
character(100), dimension(:), allocatable :: file_sect_displ, file_sect_displ_short
character(100) :: ChemFormula
character(32), dimension(10) :: char_var
character(10) :: temp_ch, chtemp11, sh_cmd
character(1) :: path_sep

type(Instant_data), dimension(:), allocatable :: Step           ! All data and parameters at this timestep
type(Displacement_analysis), dimension(:), allocatable :: Displ ! [A] mean displacements for atoms masked

real(8) :: cut_r, time_dt, cur_t, dt, time_print, cur_mass, cur_q, time
real(8) :: MSD	! [A^MSD_power] mean displacements average over all atoms
real(8), dimension(:), allocatable :: MSDP ! [A] mean displacement of atoms for each sort of atoms in a compound
integer :: MSD_power ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>)

character(3), dimension(:), allocatable :: at_short_names ! name of the element
real(8), dimension(:), allocatable :: at_masses ! mass of each element [a.m.u.]

integer, dimension(:), allocatable :: indices, FN_MSD
real(8), dimension(:), allocatable :: fragment_masses, fragment_m_over_z
real(8), dimension(:), allocatable :: mass_grid
real(8), dimension(:,:), allocatable :: output_fragment_array, output_m_over_z_array
integer :: FN1, FN2, FN3, FN_out, FN_out1, FN_out_mz, FN_out2  ! file number
integer :: INFO, Reason, i, ii, j, Tsiz, Nat, existing_elem, at_num, cur_j, iret, KOA, Nsiz
logical :: read_well

call Path_separator(path_sep)  ! Objects_and_types

! Set defaults:
INFO = 0
FN_out = 1000
FN_out_mz = 1001
!FN_MSD = 1002
FN1 = 9999
FN2 = 9998


File_supce = 'OUTPUT_supercell.dat'	! defaul name of XTANT out file with supercell sizes
File_XYZ = 'OUTPUT_atomic_coordinates.xyz'
File_fragments = 'OUT_fragments_spectrum.dat'   ! default name, for start
File_m_over_z = 'OUT_m_over_z_spectrum.dat'     ! default name, for start
File_atomic_masks = 'Atomic_masks.txt'          ! default name of the file with atomic masks for diplacement analysis
!File_displacements = 'OUT_displacements'        ! part of the name with atomic displacements

! Default values:
cut_r = 5.0d0     ! [A] default cut off radius for separation of fragments
time_dt = 100.0  ! calculate autocorrelators every 200 fs
time_print = 0.0d0  ! to start with

!---------------------------------------
print*, '******************************************************************************'
print*, 'For analysis of fragmentation, call it with the options:'
print*, 'XTANT_fragments.exe r dt'
print*, 'where you set numbers: r, dt'
print*, 'r is the cut off radius in [A] (default r=5 A)'
print*, 'dt is the printout timesteps [fs] (default dt=100 fs)'
print*, '******************************************************************************'
print*, 'For displacement analysis of section of the sample, a file'
print*, 'Atomic_masks.txt must be provided, in the same format as input files with masks.'

! Get r and timestep, if user defined it:
do i = 1, iargc()
   call getarg(i, char_var(i), status=Reason)	! read only the first one
   select case (i)
   case (1)
       read(char_var(i),*) cut_r   ! cutoff radius [A]
   case (2)
       read(char_var(i),*) time_print  ! timestep for calculation of fragments [fs]
   case default
       write(*,'(a)') char_var(i)
   endselect
end do
print*, cut_r, time_print

!-----------------------------------
! Get the time grid from the file:
open (unit=FN1, file=trim(adjustl(File_supce)), status = 'old', readonly)
open (unit=FN2, file=trim(adjustl(File_XYZ)), status = 'old', readonly)
call read_time_step(FN1, FN2, Step, time_print)    ! below
Tsiz = size(Step)   ! time grid size
close(FN1)
close(FN2)

print*, 'Output timesteps to be analyzed: ', Tsiz

!-----------------------------------
! Construct chemical formula of the compound:
ChemFormula = ''  ! to start with
KOA = 0  ! to start with
Nat = size(Step(1)%R(1,:)) ! number of atoms in the simulation box
do i = 1, Nat
   ! Find if this element repeats, or a new one:
   existing_elem = INDEX(trim(adjustl(ChemFormula)), trim(adjustl(Step(1)%name(i))))   ! intrinsic function
   if (existing_elem == 0) then ! new element
      KOA = KOA + 1  ! new element index
      ! Add it to the list (construct full chemical formula):
      ChemFormula = trim(adjustl(ChemFormula))//trim(adjustl(Step(1)%name(i)))
   endif
   ! Kind of atom:
   Step(1)%KOA(i) = KOA
   !print*, i, existing_elem, trim(adjustl(ChemFormula)), ' : ', trim(adjustl(Step(1)%name(i)))
enddo
print*, 'The following elements are in the compound: ', trim(adjustl(ChemFormula))

!-----------------------------------
! Find the masses of each element:
call Decompose_compound('', ChemFormula, path_sep, INFO, error_message, at_num, at_short_names=at_short_names, at_masses=at_masses) ! molude 'Periodic_Table'
if (INFO .NE. 0) then
   print*, trim(adjustl(error_message))
   goto 2012
endif
print*, 'Interpreted as:'
do i = 1, at_num
   write(*,'(i3,a,a,a,f)') i, ') element: ', trim(adjustl(at_short_names(i))), ', mass: ', at_masses(i)
enddo

! Knowing the masses, set the mass grid for fragments:
cur_j = 0
do i = 1, Nat  ! save mass of each atom:
   INNDloop:do j = 1, at_num
      if ( trim(adjustl(Step(1)%name(i))) == trim(adjustl(at_short_names(j))) ) then
         cur_j = j
         exit INNDloop
      endif
   enddo INNDloop
   Step(1)%M(i) = at_masses(cur_j)
enddo
allocate(mass_grid( CEILING( SUM( Step(1)%M(:) ) ) ) , source=0.0d0)
allocate(fragment_masses( CEILING( SUM( Step(1)%M(:) ) ) ) , source=0.0d0)
allocate(fragment_m_over_z( CEILING( SUM( Step(1)%M(:) ) ) ) , source=0.0d0)
! Create mass grid:
do i = 1, size(mass_grid)
   mass_grid(i) = dble(i)
enddo


!-----------------------------------
call interprete_displacement_command(File_atomic_masks, Displ, MSD_power, Reason)  ! below


!-----------------------------------
! Output file:
open (unit=FN_out, file=trim(adjustl(File_fragments)))
open (unit=FN_out_mz, file=trim(adjustl(File_m_over_z)))

Nsiz = size(Displ)   ! how many masks
if (.not. allocated(FN_MSD)) then  ! allocate and open files
   allocate(FN_MSD(Nsiz))
endif
if (.not. allocated(file_sect_displ)) then  ! allocate and open files
   allocate(file_sect_displ(Nsiz))
   !allocate(file_sect_displ_short(Nsiz))
endif
do i = 1, Nsiz ! for all masks
   file_sect_displ(i) = 'OUT_displacements_'//trim(adjustl(Displ(i)%mask_name))//'.dat'
   open(NEWUNIT=FN_MSD(i), FILE = trim(adjustl(file_sect_displ(i))))

   N_at = at_num    ! number of kinds of atoms
   chtemp = ''  ! to start with
   do j = 1, N_at
      chtemp = trim(adjustl(chtemp))//'   '//trim(adjustl(at_short_names(j)))//':total   X  Y  Z'
   enddo
   call create_file_header(FN_MSD(i), '#Time  Total X  Y  Z  '//trim(adjustl(chtemp)) )
   if (INT(Displ(i)%MSD_power) > 1) then
      write(chtemp11,'(i2)') INT(Displ(i)%MSD_power)
      call create_file_header(FN_MSD(i), '#[fs]   [A^'//trim(adjustl(chtemp11))//'](:)')
   else
      call create_file_header(FN_MSD(i), '#[fs]   [A](:)')
   endif
enddo

! Set the output array:
allocate(output_fragment_array(Tsiz, size(mass_grid)))
allocate(output_m_over_z_array(Tsiz, size(mass_grid)))

!-----------------------------------
! Analyse fragments:
cur_t = Step(1)%Tim    ! start from here
do i = 1, Tsiz ! time steps
   fragment_masses = 0.0d0 ! reset for each timestep
   fragment_m_over_z = 0.0d0  ! reset for each timestep

   ! Get the number of fragments and their indices:
   call get_fragments_indices(Step, i, cut_r, indices)  ! below
   print*, 'time:', i, 'fragments:', maxval(indices)

   ! Sort atoms to fragments, and save the fragments parametes:
   do j = 1, maxval(indices)  ! all fragments
      cur_mass = SUM(Step(1)%M(:), MASK=( indices(:) == j ) )  ! mass of the fragment
      fragment_masses(ANINT(cur_mass)) = fragment_masses(ANINT(cur_mass)) + 1

      if (allocated(Step(i)%q)) then
         cur_q = SUM(Step(i)%q(:), MASK=( indices(:) == j ) )  ! charge of the fragment
         call sort_fragment_m_over_z(cur_mass, cur_q, fragment_m_over_z) ! below
      endif
   enddo

   ! Save for output:
   output_fragment_array(i, :) = fragment_masses(:)
   if (allocated(Step(i)%q)) output_m_over_z_array(i, :) = fragment_m_over_z(:)


   ! Calculate the atomic displacements for each mask:
   call get_mean_square_displacement(Step, i, at_num, Displ, MSD, MSDP, MSD_power) ! below
   Nsiz = size(Displ)   ! how many masks
   time = Step(i)%Tim
   do ii = 1, Nsiz    ! for all masks
      write(FN_MSD(ii), '(es25.16,$)') time, Displ(ii)%mean_disp, Displ(ii)%mean_disp_r(:)
      ! Now for kinds of atoms:
      N_at = at_num    ! number of kinds of atoms
      do j = 1, N_at
         write(FN_MSD(ii), '(es25.16,$)') Displ(ii)%mean_disp_sort(j), Displ(ii)%mean_disp_r_sort(j,:)
      enddo
      write(FN_MSD(ii),'(a)') ! make a new line
   enddo ! ii
enddo


!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! Printout:
!write(FN_out, '(f)', advance='no') 0.0e0  ! corner
if (Tsiz > 1) then
   write(FN_out, '(f)', advance='no') Step(1)%Tim - (Step(2)%Tim - Step(1)%Tim)  ! step "-1"
   if (allocated(Step(1)%q)) write(FN_out_mz, '(f)', advance='no') Step(1)%Tim - (Step(2)%Tim - Step(1)%Tim)  ! step "-1"
else
   write(FN_out, '(f)', advance='no') Step(1)%Tim  ! step "-1"
   if (allocated(Step(1)%q)) write(FN_out_mz, '(f)', advance='no') Step(1)%Tim   ! step "-1"
endif

do i = 1, Tsiz ! time steps
   write(FN_out, '(f)', advance='no') Step(i)%Tim  ! time grid as columns
   if (allocated(Step(1)%q)) write(FN_out_mz, '(f)', advance='no') Step(i)%Tim  ! time grid as columns
enddo
write(FN_out, '(a)') ''
if (allocated(Step(1)%q)) write(FN_out_mz, '(a)') ''

do j = 1, size(fragment_masses) ! Save output file
   write(FN_out, '(f10.2)', advance='no') mass_grid(j)
   if (allocated(Step(1)%q)) write(FN_out_mz, '(f10.2)', advance='no') (mass_grid(j)-1.0d0)
   do i = 1, Tsiz ! time steps
      write(FN_out, '(f10.2)', advance='no') output_fragment_array(i,j)
      if (allocated(Step(1)%q)) write(FN_out_mz, '(f10.2)', advance='no') output_m_over_z_array(i,j)
   enddo
   write(FN_out, '(a)') ''
   if (allocated(Step(1)%q)) write(FN_out_mz, '(a)') ''
enddo

close(FN_out)
if (allocated(Step(1)%q)) close(FN_out_mz)




!-------------------------
! Make gnuplot script for mass spectrum:
FN_out1 = 9997
Gnu_script = 'OUT_fragments'
Gnu_file = 'OUT_fragments.png'
if (path_sep .EQ. '\') then	! if it is Windows
   open (unit=FN_out1, file=trim(adjustl(Gnu_script))//'.cmd')
   write(FN_out1, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   write(FN_out1, '(a)') 'LABL="Fragments mass spectrum"'
   write(FN_out1, '(a)') 'set terminal png font arial'
   write(FN_out1, '(a)') 'set output "OUT_fragments.png"'
   write(FN_out1, '(a)') 'set xlabel "Time (fs)"'
   write(FN_out1, '(a)') 'set ylabel "Mass spectrum (a.m.u.)"'
   write(FN_out1, '(a)') 'set autoscale xfix'
   write(FN_out1, '(a)') 'set autoscale cbfix'
   write(FN_out1, '(a)') 'set palette defined (0 "white",\'
   write(FN_out1, '(a)') '0.01 "blue",\'
   write(FN_out1, '(a)') '0.4 "purple",\'
   write(FN_out1, '(a)') '0.6 "red",\'
   write(FN_out1, '(a)') '0.8 "yellow",\'
   write(FN_out1, '(a)') '1.0 "light-green",\'
   write(FN_out1, '(a)') '1.2 "green")'
   write(FN_out1, '(a)') "plot[][0:50] 'OUT_fragments_spectrum.dat' matrix nonuniform with image"//' title"Fragments mass spectrum"'
else ! it is linux
   open (unit=FN_out1, file=trim(adjustl(Gnu_script))//'.sh')
   write(FN_out1, '(a)') '#!/bin/bash'
   write(FN_out1, '(a)') ''
   write(FN_out1, '(a)') 'NAME='//trim(adjustl(Gnu_file))
   write(FN_out1, '(a)') 'LABL="Fragments mass spectrum"'
   write(FN_out1, '(a)') ' echo "'
   write(FN_out1, '(a)') 'set terminal png font arial'
   write(FN_out1, '(a)') 'set output \"$NAME\"'
   write(FN_out1, '(a)') 'set xlabel \"Time (fs) \" '
   write(FN_out1, '(a)') 'set ylabel \"Mass spectrum (a.m.u.) \" '
   write(FN_out1, '(a)') 'set autoscale xfix'
   write(FN_out1, '(a)') 'set autoscale cbfix'
   write(FN_out1, '(a)') 'set palette defined (0 "white",\'
   write(FN_out1, '(a)') '0.01 "blue",\'
   write(FN_out1, '(a)') '0.4 "purple",\'
   write(FN_out1, '(a)') '0.6 "red",\'
   write(FN_out1, '(a)') '0.8 "yellow",\'
   write(FN_out1, '(a)') '1.0 "light-green",\'
   write(FN_out1, '(a)') '1.2 "green")'
   write(FN_out1, '(a)') "plot[][0:50] 'OUT_fragments_spectrum.dat' matrix nonuniform with image title"//'\"Fragments mass spectrum\"'
   write(FN_out1, '(a)') 'reset'
   write(FN_out1, '(a)') '" | gnuplot '
   !call system('chmod +x '//trim(adjustl(Gnu_script))//'.sh') ! make the output-script executable
   !command = 'chmod +x '//trim(adjustl(File_name))
   command = 'chmod +x '//trim(adjustl(Gnu_script))//'.sh'
   iret = system(command)
endif
close(FN_out1)

if (path_sep .EQ. '\') then	! if it is Windows
   !call system(trim(adjustl(Gnu_script))//'.cmd')
!    command = "OUTPUT_Gnuplot_all.cmd"
   command = 'call '//trim(adjustl(Gnu_script))//'.cmd'
   iret = system(command)
else ! linux:
   !call system(trim(adjustl(Gnu_script))//'.sh')
!    command = "./OUTPUT_Gnuplot_all.sh"
   command = './'//trim(adjustl(Gnu_script))//'.sh'
   iret = system(command)
endif



!-------------------------
! Make gnuplot script for m/z spectrum:
FN_out2 = 9996
Gnu_script = 'OUT_m_over_z'
Gnu_file = 'OUT_m_over_z.png'
if (path_sep .EQ. '\') then	! if it is Windows
   open (unit=FN_out2, file=trim(adjustl(Gnu_script))//'.cmd')
   write(FN_out2, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   write(FN_out2, '(a)') 'LABL="Fragments mass spectrum m/z"'
   write(FN_out2, '(a)') 'set terminal png font arial'
   write(FN_out2, '(a)') 'set output "OUT_m_over_z.png"'
   write(FN_out2, '(a)') 'set xlabel "Time (fs)"'
   write(FN_out2, '(a)') 'set ylabel "Mass spectrum (m/z)"'
   write(FN_out2, '(a)') 'set autoscale xfix'
   write(FN_out2, '(a)') 'set autoscale cbfix'
   write(FN_out2, '(a)') 'set palette defined (0 "white",\'
   write(FN_out2, '(a)') '0.01 "blue",\'
   write(FN_out2, '(a)') '0.4 "purple",\'
   write(FN_out2, '(a)') '0.6 "red",\'
   write(FN_out2, '(a)') '0.8 "yellow",\'
   write(FN_out2, '(a)') '1.0 "light-green",\'
   write(FN_out2, '(a)') '1.2 "green")'
   write(FN_out2, '(a)') "plot[][0:50] 'OUT_m_over_z_spectrum.dat' matrix nonuniform with image"//' title"Fragments mass spectrum m/z"'
else ! it is linux
   open (unit=FN_out2, file=trim(adjustl(Gnu_script))//'.sh')
   write(FN_out2, '(a)') '#!/bin/bash'
   write(FN_out2, '(a)') ''
   write(FN_out2, '(a)') 'NAME='//trim(adjustl(Gnu_file))
   write(FN_out2, '(a)') 'LABL="Fragments mass spectrum m/z"'
   write(FN_out2, '(a)') ' echo "'
   write(FN_out2, '(a)') 'set terminal png font arial'
   write(FN_out2, '(a)') 'set output \"$NAME\"'
   write(FN_out2, '(a)') 'set xlabel \"Time (fs) \" '
   write(FN_out2, '(a)') 'set ylabel \"Mass spectrum (m/z) \" '
   write(FN_out2, '(a)') 'set autoscale xfix'
   write(FN_out2, '(a)') 'set autoscale cbfix'
   write(FN_out2, '(a)') 'set palette defined (0 "white",\'
   write(FN_out2, '(a)') '0.01 "blue",\'
   write(FN_out2, '(a)') '0.4 "purple",\'
   write(FN_out2, '(a)') '0.6 "red",\'
   write(FN_out2, '(a)') '0.8 "yellow",\'
   write(FN_out2, '(a)') '1.0 "light-green",\'
   write(FN_out2, '(a)') '1.2 "green")'
   write(FN_out2, '(a)') "plot[][0:50] 'OUT_m_over_z_spectrum.dat' matrix nonuniform with image title"//'\"Fragments mass spectrum m/z\"'
   write(FN_out2, '(a)') 'reset'
   write(FN_out2, '(a)') '" | gnuplot '
   !call system('chmod +x '//trim(adjustl(Gnu_script))//'.sh') ! make the output-script executable
   !command = 'chmod +x '//trim(adjustl(File_name))
   command = 'chmod +x '//trim(adjustl(Gnu_script))//'.sh'
   iret = system(command)
endif
close(FN_out2)

if (path_sep .EQ. '\') then	! if it is Windows
   !call system(trim(adjustl(Gnu_script))//'.cmd')
!    command = "OUTPUT_Gnuplot_all.cmd"
   sh_cmd = '.cmd'
   command = trim(adjustl(Gnu_script))//'.cmd'
   iret = system(command)
else ! linux:
   !call system(trim(adjustl(Gnu_script))//'.sh')
!    command = "./OUTPUT_Gnuplot_all.sh"
   sh_cmd = '.sh'
   command = trim(adjustl(Gnu_script))//'.sh'
   iret = system(command)
endif



!-------------------------
! Make gnuplot script for displacements:
! Atomic masks for sectional displacements:
if (allocated(Displ)) then
   Nsiz = size(Displ)   ! how many masks
   do j = 1, Nsiz    ! for all masks
      File_name = 'OUT_displacements_'//trim(adjustl(Displ(j)%mask_name))//'_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_displacements(File_name, file_sect_displ(j), Step(1)%Tim, Step(size(Step))%Tim, 'OUT_mean_displacement_'// &
            trim(adjustl(Displ(j)%mask_name))//'.jpg', &
            Displ(j)%MSD_power, path_sep) ! below
      iret = system(File_name)

      ! Partial by elements, if there is more than one:
      File_name = 'OUT_displacements_'//trim(adjustl(Displ(j)%mask_name))//'_partial_Gnuplot'//trim(adjustl(sh_cmd))
      if (at_num > 1) then
         call gnu_displacements_partial(File_name, file_sect_displ(j), Step(1)%Tim, Step(size(Step))%Tim, 'OUT_mean_displacement_'// &
            trim(adjustl(Displ(j)%mask_name))//'_partial.jpg', &
            Displ(j)%MSD_power, at_num, at_short_names, path_sep) ! below
         iret = system(File_name)
      endif
   enddo ! j
endif



2012 continue   ! to exit the program
STOP
!---------------------
 contains




subroutine gnu_displacements_partial(File_name, file_MSD, t0, t_last, eps_name, MSD_power, at_num, at_short_names, path_sep)
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_MSD	! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   real(8), intent(in) :: MSD_power ! power of MSD
   !type(Solid), intent(in) :: matter     ! material parameters
   integer, intent(in) :: at_num    ! number of kinds of atoms
   character(*), dimension(:), intent(in) :: at_short_names
   character(*), intent(in) :: path_sep
   !------------------------
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp, MSD_text

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   if (MSD_power > 1) then
      ! Power of MSD:
      write(MSD_text,'(i2)') int(MSD_power)

      call write_gnuplot_script_header_new(FN, 2, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), path_sep, 2)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, 2, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), path_sep, 2)	! module "Gnuplotting"
   endif

   i_start = 6
   do i = 1, at_num
      chtemp = trim(adjustl(at_short_names(i)))
      if (path_sep .EQ. '\') then	! if it is Windows
         if (i == 1) then
            write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:', i_start , &
                     ' w l lw LW title "'//trim(adjustl(chtemp))//'" ,\'
         else
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start + (i-1)*4,' w l lw LW title "'//trim(adjustl(chtemp))//'" ,\'
         endif
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+1+(i-1)*4 ,' w l lw LW title " ', trim(adjustl(chtemp))//':X'  ,' " ,\'
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+2+(i-1)*4,' w l lw LW title " ', trim(adjustl(chtemp))//':Y'  ,' " ,\'
         if (i /= at_num) then
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+3+(i-1)*4,' w l lw LW title " ', trim(adjustl(chtemp))//':Z'  ,' " ,\'
         else
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+3+(i-1)*4,' w l lw LW title " ', trim(adjustl(chtemp))//':Z'  ,' " '
         endif

      else
         if (i == 1) then
            write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:', i_start, &
                                          ' w l lw \"$LW\" title \"'//trim(adjustl(chtemp))//'\" ,\'
         else
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+(i-1)*4, ' w l lw \"$LW\" title title \"'//trim(adjustl(chtemp))//'\" ,\'
         endif
         write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+1+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':X' ,'\" ,\'
         write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+2+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':Y' ,'\" ,\'

         if (i /= at_num) then
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+3+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':Z' ,'\" ,\'
         else
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+3+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':Z' ,'\"'
         endif
      endif
   enddo
   call write_gnuplot_script_ending(FN, File_name, 1, path_sep)
   close(FN)
end subroutine gnu_displacements_partial


subroutine gnu_displacements(File_name, file_MSD, t0, t_last, eps_name, MSD_power, path_sep)
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_MSD, path_sep
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   real(8), intent(in) :: MSD_power ! power of MSD
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp, MSD_text

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   if (MSD_power > 1) then
      ! Power of MSD:
      write(MSD_text,'(i2)') int(MSD_power)

      call write_gnuplot_script_header_new(FN, 2, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), path_sep, 2)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, 2, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), path_sep, 2)	! module "Gnuplotting"
   endif

   i_start = 2
   if (path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:', i_start ,' w l lw LW title "Average" ,\'
      write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+1 ,' w l lw LW title " ', 'X'  ,' " ,\'
      write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+2 ,' w l lw LW title " ', 'Y'  ,' " ,\'
      write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+3 ,' w l lw LW title " ', 'Z'  ,' " '
   else
      write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:', i_start, &
                                          ' w l lw \"$LW\" title \"Average\" ,\'
      write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+1, ' w l lw \"$LW\" title \" ', 'X' ,'\" ,\'
      write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+2, ' w l lw \"$LW\" title \" ', 'Y' ,'\" ,\'
      write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+3, ' w l lw \"$LW\" title \" ', 'Z' ,'\"'
   endif

   call write_gnuplot_script_ending(FN, File_name, 1, path_sep)
   close(FN)
end subroutine gnu_displacements


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


subroutine get_mean_square_displacement(MDdata, i_step, at_num, Displ, MSD, MSDP, MSD_power)	! currently, it calculates mean displacement, without sqaring it
   !type(Super_cell), dimension(:), intent(inout), target :: Scell	! super-cell with all the atoms inside
   !type(Solid), intent(in) :: matter     ! material parameters
   type(Instant_data), dimension(:), intent(inout), target :: MDdata   ! atomic coordinates and data
   integer, intent(in) :: i_step   ! step number
   integer, intent(in) :: at_num ! how many kinds of atoms
   type(Displacement_analysis), dimension(:), allocatable, intent(inout) :: Displ ! [A] mean displacements for atoms masked
   real(8), intent(out) :: MSD	! [A^MSD_power] mean displacements average over all atoms
   real(8), dimension(:), allocatable, intent(out) :: MSDP ! [A] mean displacement of atoms for each sort of atoms in a compound
   integer, intent(in) :: MSD_power ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>)
   !-------------------------
   integer :: N, iat, ik, i,  j, k, Nat, Nsiz, i_masks
   integer, pointer :: KOA
   real(8) :: zb(3), x, y, z, a_r, r1, x0, y0, z0
   real(8), dimension(:), pointer :: S, S0

   if (.not.allocated(MSDP)) allocate(MSDP(at_num))
   N = size(MDdata(1)%name)   ! number of atoms (assume constant)

   ! Check if user defined any atomic masks:
   if (allocated(Displ)) then
      call update_atomic_masks_displ(MDdata, i_step, Displ, at_num) ! below
      ! Restart counting for this step:
      Nsiz = size(Displ)
      do i_masks = 1, Nsiz ! for all requested masks
         Displ(i_masks)%mean_disp = 0.0d0
         Displ(i_masks)%mean_disp_sort(:) = 0.0d0
         Displ(i_masks)%mean_disp_r(:) = 0.0d0
         Displ(i_masks)%mean_disp_r_sort = 0.0d0
      enddo
   endif

   ! Get equilibrium relative coordinates from given absolute coordinates inside of the current supercell:
   call get_coords_in_new_supce(MDdata, i_step)	! below (S_eq arae updated, R_eq do not change)

   MSD = 0.0d0	! to start with
   MSDP = 0.0d0 ! to start with
   do iat = 1, N	! for all atoms
      KOA => MDdata(1)%KOA(iat)
      S => MDdata(i_step)%S(:,iat)
      S0 => MDdata(1)%S(:,iat)
      a_r = 1.0d31	! just to start from
      ! For the case of periodical boundaries:
      do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
         ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell:
         ! periodic boundary conditions
         zb(1) = dble(i)
         do j = -1,1
            zb(2) = dble(j)
            do k = -1,1
               zb(3) = dble(k)
               x0 = 0.0d0
               y0 = 0.0d0
               z0 = 0.0d0
               do ik = 1,3
                  x0 = x0 + (S(ik) - S0(ik) + zb(ik))*MDdata(i_step)%supce(ik,1) ! correct
                  y0 = y0 + (S(ik) - S0(ik) + zb(ik))*MDdata(i_step)%supce(ik,2)
                  z0 = z0 + (S(ik) - S0(ik) + zb(ik))*MDdata(i_step)%supce(ik,3)
               enddo ! ik
               r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
               if (r1 .LT. a_r) then
                  x = x0
                  y = y0
                  z = z0
                  a_r = r1
               endif !  (r1 .LT. a_r)
            enddo ! k
         enddo ! j
      enddo ! i
      MSD = MSD + a_r**MSD_power ! mean displacement^N
      MSDP(KOA) = MSDP(KOA) + a_r**MSD_power    ! mean displacement^N

      ! Section of atoms according to masks, if any:
      if (allocated(Displ)) then
         Nsiz = size(Displ)
         do i_masks = 1, Nsiz ! for all requested masks
            if (Displ(i_masks)%Atomic_mask(iat)) then ! this atom is included in the mask
               r1 = a_r**Displ(i_masks)%MSD_power  ! [A^N] displacement
               Displ(i_masks)%mean_disp = Displ(i_masks)%mean_disp + r1
               Displ(i_masks)%mean_disp_sort(KOA) = Displ(i_masks)%mean_disp_sort(KOA) + r1
               ! Along axes:
               r1 = x**Displ(i_masks)%MSD_power
               Displ(i_masks)%mean_disp_r(1) = Displ(i_masks)%mean_disp_r(1) + r1
               Displ(i_masks)%mean_disp_r_sort(KOA,1) = Displ(i_masks)%mean_disp_r_sort(KOA,1) + r1
               r1 = y**Displ(i_masks)%MSD_power
               Displ(i_masks)%mean_disp_r(2) = Displ(i_masks)%mean_disp_r(2) + r1
               Displ(i_masks)%mean_disp_r_sort(KOA,2) = Displ(i_masks)%mean_disp_r_sort(KOA,2) + r1
               r1 = z**Displ(i_masks)%MSD_power
               Displ(i_masks)%mean_disp_r(3) = Displ(i_masks)%mean_disp_r(3) + r1
               Displ(i_masks)%mean_disp_r_sort(KOA,3) = Displ(i_masks)%mean_disp_r_sort(KOA,3) + r1
            endif
         enddo ! i_masks
      endif ! (allocated(Displ))
   enddo ! iat


   MSD = MSD/dble(N)	! averaged over all atoms
   ! Section of atoms according to masks, if any:
   if (allocated(Displ)) then
      Nsiz = size(Displ)
      do i_masks = 1, Nsiz ! for all requested masks
         Nat = COUNT(MASK = Displ(i_masks)%Atomic_mask)
         if (Nat > 0) then
            Displ(i_masks)%mean_disp = Displ(i_masks)%mean_disp / Nat
            Displ(i_masks)%mean_disp_r(:) = Displ(i_masks)%mean_disp_r(:) / Nat
         else
            Displ(i_masks)%mean_disp = 0.0d0
            Displ(i_masks)%mean_disp_r(:) = 0.0d0
         endif
      enddo
   endif


   ! For all elements:
   do i = 1, at_num
      ! how many atoms of this kind are in the supercell:
      Nat = COUNT(MASK = (MDdata(1)%KOA(:) == i))
      if (Nat > 0) then
         MSDP(i) = MSDP(i) / dble(Nat)
      else
         MSDP(i) = 0.0d0
      endif

      ! Section of atoms according to masks, if any:
      if (allocated(Displ)) then
         Nsiz = size(Displ)
         do i_masks = 1, Nsiz ! for all requested masks
            Nat = COUNT(MASK = (Displ(i_masks)%Atomic_mask(:) .and. (MDdata(1)%KOA(:) == i) ))
            if (Nat > 0) then
               Displ(i_masks)%mean_disp_sort(i) = Displ(i_masks)%mean_disp_sort(i) / Nat
               Displ(i_masks)%mean_disp_r_sort(i,:) = Displ(i_masks)%mean_disp_r_sort(i,:) / Nat
            else
               Displ(i_masks)%mean_disp_sort(i) = 0.0d0
               Displ(i_masks)%mean_disp_r_sort(i,:) = 0.0d0
            endif
         enddo
      endif
   enddo

!    do i_masks = 1, Nsiz ! for all requested masks
!       print*, trim(adjustl(Displ(i_masks)%mask_name)), i_masks, MSD, Displ(i_masks)%mean_disp
!       print*, MSDP(:)
!       print*, Displ(i_masks)%mean_disp_sort(:)
!       print*, Displ(i_masks)%mean_disp_r(:)
!       print*, 'K', Displ(i_masks)%mean_disp_r_sort
!    enddo

   nullify(S,S0,KOA)
end subroutine get_mean_square_displacement



subroutine get_coords_in_new_supce(MDdata, i_step) !  (S_eq are updated, R_eq do not change)
   !type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   !integer, intent(in) :: NSC ! number of super-cell
   type(Instant_data), dimension(:), intent(inout) :: MDdata   ! atomic coordinates and data
   integer, intent(in) :: i_step   ! step number
   !---------------
   real(8) :: sx, sy, sz
   real(8), dimension(3,3) :: supce_inv
   integer j, N, ik
   N = size(MDdata(i_step)%name)
   call Invers_3x3(MDdata(i_step)%supce, supce_inv, 'get_coords_in_new_supce')	! module "Algebra_tools"
   do j = 1,N	! all atoms:
      sx = 0.0d0
      sy = 0.0d0
      sz = 0.0d0
      do ik = 1,3
         sx = sx + MDdata(1)%R(ik,j)*supce_inv(ik,1) ! correct
         sy = sy + MDdata(1)%R(ik,j)*supce_inv(ik,2)
         sz = sz + MDdata(1)%R(ik,j)*supce_inv(ik,3)
      enddo ! ik
      MDdata(1)%S(1,j) = sx
      MDdata(1)%S(2,j) = sy
      MDdata(1)%S(3,j) = sz
   enddo ! j
end subroutine get_coords_in_new_supce



subroutine update_atomic_masks_displ(MDdata, i_step, Displ, at_num)
   !type(Super_cell), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(Instant_data), dimension(:), intent(in) :: MDdata
   integer, intent(in) :: i_step   ! step number
   type(Displacement_analysis), dimension(:), intent(inout) :: Displ ! [A] mean displacements for atoms masked
   !type(Solid), intent(in) :: matter     ! material parameters
   integer, intent(in) :: at_num ! number of elements in the compound
   !-----------------
   integer :: N_at, Nsiz, i, iat, N_KAO
   logical :: mask_1, mask_2

   N_at = size(MDdata(1)%KOA)	! number of atoms
   N_KAO = at_num

   Nsiz = size(Displ)
   do i = 1, Nsiz ! for all requested masks
      ! Make sure the arrays are allocated:
      if (.not.allocated(Displ(i)%mean_disp_sort)) allocate(Displ(i)%mean_disp_sort(at_num))
      if (.not.allocated(Displ(i)%mean_disp_r_sort)) allocate(Displ(i)%mean_disp_r_sort(at_num,3))
      if (.not.allocated(Displ(i)%Atomic_mask)) allocate(Displ(i)%Atomic_mask(N_at))

      ! Create or update the masks:
      ! What type of mask is it:
      select case( trim(adjustl(Displ(i)%mask_name(1:7))) )
      case default ! all atoms, no selection
         Displ(i)%Atomic_mask = .true. ! all atoms included

      case ('Section', 'section', 'SECTION') ! spatial section of atoms
         Displ(i)%Atomic_mask = .false. ! to start with
         do iat = 1, N_at  ! for all atoms
            mask_1 = .false.  ! to start with
            mask_2 = .false.  ! to start with

            ! Mask #1:
            if ( (MDdata(i_step)%R(1,iat) > Displ(i)%r_start(1, 1) ) .and. &
                 (MDdata(i_step)%R(1,iat) < Displ(i)%r_end(1, 1) )  .and. & ! X
                 (MDdata(i_step)%R(2,iat) > Displ(i)%r_start(1, 2) ) .and. &
                 (MDdata(i_step)%R(2,iat) < Displ(i)%r_end(1, 2) )  .and. & ! Y
                 (MDdata(i_step)%R(3,iat) > Displ(i)%r_start(1, 3) ) .and. &
                 (MDdata(i_step)%R(3,iat) < Displ(i)%r_end(1, 3) ) ) then ! Z
               mask_1 = .true.
            endif

            ! Mask #2, if present:
            if (Displ(i)%logical_and .or. Displ(i)%logical_or) then
               if (  (MDdata(i_step)%R(1,iat) > Displ(i)%r_start(2, 1) ) .and. &
                     (MDdata(i_step)%R(1,iat) < Displ(i)%r_end(2, 1) )  .and. & ! X
                     (MDdata(i_step)%R(2,iat) > Displ(i)%r_start(2, 2) ) .and. &
                     (MDdata(i_step)%R(2,iat) < Displ(i)%r_end(2, 2) )  .and. & ! Y
                     (MDdata(i_step)%R(3,iat) > Displ(i)%r_start(2, 3) ) .and. &
                     (MDdata(i_step)%R(3,iat) < Displ(i)%r_end(2, 3) ) ) then ! Z

                  Displ(i)%Atomic_mask(iat) = Displ(i)%Atomic_mask(iat)
               endif
            endif

            ! Combine masks:
            if (Displ(i)%logical_and) then  ! both
               Displ(i)%Atomic_mask(iat) = (mask_1 .and. mask_2)
            elseif (Displ(i)%logical_or) then  ! either
               Displ(i)%Atomic_mask(iat) = (mask_1 .or. mask_2)
            else  ! only one mask:
               Displ(i)%Atomic_mask(iat) = mask_1
            endif

         enddo ! iat = 1, N_at
      end select

   enddo
end subroutine update_atomic_masks_displ




subroutine interprete_displacement_command(read_line, Displ, MSD_power, Reason)
   character(*), intent(in) :: read_line
   !type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Displacement_analysis), dimension(:), allocatable, intent(inout) :: Displ ! [A] mean displacements for atoms masked
   !type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   integer, intent(out) :: MSD_power
   integer, intent(inout) :: Reason
   !-------------------------------
   character(200) :: Folder_name, File_name, sting_temp
   character(10) :: ch_temp
   logical :: file_exists
   integer :: FN1, ch_int, i

   ! Try to interprete it as a number (power of mean displacement)
   ! to make back-compatible with the legacy format:
   read(read_line,*,IOSTAT=Reason) MSD_power
   if (Reason == 0) then
      print*, 'Atomic displacement analysis is set in legacy format'
      return ! it was a number, nothing more to do
   endif

   !---------------
   ! If it was not a command, check if it was a file
   !Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep   ! directory with input files
   Folder_name = ''  ! Directory with output files!
   read(read_line,*,IOSTAT=Reason) sting_temp   ! read first block
   File_name = trim(adjustl(Folder_name))//trim(adjustl(sting_temp))  ! file name with masks
   inquire(file=trim(adjustl(File_name)),exist=file_exists)
   if (file_exists) then
      FN1 = 4001
      open(UNIT=FN1, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      read(FN1,*,IOSTAT=Reason) ch_temp, ch_int
      if (Reason == 0) then   ! if read number of masks well
         ! Allocate displacements data objects:
         allocate(Displ(ch_int))
         ! Read parameters of the masks:
         do i = 1, ch_int  ! for all masks:
            read(FN1,'(a)',IOSTAT=Reason) sting_temp
            if (Reason == 0) then
               call read_displacement_command(trim(adjustl(sting_temp)), Displ, Reason, i, FN1) ! below

            else
               print*, 'Could not read atomic masks from file: ', trim(adjustl(File_name))
               print*, 'Using default mean displacement only'
               Reason = -1
               exit
            endif
         enddo ! i
      endif ! (Reason == 0)
      close (FN1)
      if (Reason == 0) then
         print*, 'Atomic masks for displacements are read from file: ', trim(adjustl(File_name))
         return ! read all, nothing more to do
      endif
   else ! check if it is a command:
      ! if it is not a number, check if it is a command:
      call read_displacement_command(trim(adjustl(read_line)), Displ, Reason, 1) ! below
      if (Reason == 0) then
         print*, 'Atomic masks for displacements are set by command'
         return ! it was a number, nothing more to do
      else
         print*, 'Could not interprete command for atomic masks'
         print*, 'Using default mean displacement only'
      endif
   endif

   !---------------
   ! 4) if it was none of the above, just use default: mean displacement N=1:
   MSD_power = 1
   Reason = 0
   print*, 'Default atomic displacement analysis is used'
end subroutine interprete_displacement_command


subroutine read_displacement_command(read_line, Displ, Reason, mask_num, FN1)
   character(*), intent(in) :: read_line
   !type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Displacement_analysis), dimension(:), allocatable, intent(inout) :: Displ ! [A] mean displacements for atoms masked
   integer, intent(inout) :: Reason
   integer, intent(in) :: mask_num  ! number of mask
   integer, intent(in), optional :: FN1 ! file number with masks definition
   !-------------------------------
   integer :: Nsiz, N_char, N_len, i, j
   character(100) :: ch_comm, ch_val
   character(200) :: ch_temp

   Reason = -1 ! to start with

   !---------------
   ! Make sure the array is allocated
   if (.not.allocated(Displ)) then
      allocate(Displ(1))
   endif

   !---------------
   ! Set defaults:
   Displ(mask_num)%MSD_power = 1.0d0   ! default: linear mean displacement
   write(ch_temp, '(i6)') mask_num
   Displ(mask_num)%mask_name = 'mask_'//trim(adjustl(ch_temp)) ! default namae
   Displ(mask_num)%print_r = .true.   ! all axis-resolved displacement

   !---------------
   ! Interprete the input line:
   N_len = LEN(read_line)  ! length of the string
   N_char = 1  ! to start with
   do while (N_char < N_len)  ! read the entire line
      read(read_line(N_char:N_len),*,IOSTAT=Reason) ch_temp ! read this block
      N_char = N_char + LEN(trim(ch_temp)) + 1  ! mark begining of the next block to read

      ! Split the string into command and its value:
      call split_command_separator(trim(adjustl(ch_temp)), ':', ch_comm, ch_val)  ! below
      ! If the command is incorrect, nothing more to do
      if (LEN(trim(adjustl(ch_comm))) == 0) cycle
      if (LEN(trim(adjustl(ch_val))) == 0) cycle

      select case(trim(adjustl(ch_comm)))
      !==========
      case ('name', 'Name', 'NAME')   ! mask name
         ! Make sure the name does not repeat:
         call number_atomic_mask(ch_val, Displ, mask_num) ! below
         Displ(mask_num)%mask_name = trim(adjustl(ch_val))
         Reason = 0

         ! Identify the mask format, and read extra parameters if any:
         select case (trim(adjustl(Displ(mask_num)%mask_name(1:3)))) ! define section
         case ('sec', 'Sec', 'SEC') ! Spatial section
            ! Read one more line with the definition of the section:
            if (present(FN1)) then ! read next line from this file
               call read_and_define_section(FN1, Displ, mask_num) ! below
            endif

         case ('all', 'All', 'ALL') ! all atoms
            ! Nothing to do, all atoms are included
         end select

      !==========
      case ('power', 'Power', 'POWER')   ! power of mean displacement
         read(ch_val,*,IOSTAT=Reason) Displ(mask_num)%MSD_power
         Reason = 0

      !==========
      case ('axis', 'Axis', 'AXIS', 'axes', 'Axes', 'AXES')   ! axis-resolved data
      ! (specification unused, all axis are printed out)
         do i = 1, LEN(trim(adjustl(ch_val)))   ! interprete all letters
            select case(trim(adjustl(ch_val(i:i))))
            case ('X', 'x')
               Displ(mask_num)%print_r(1) = .true.
            case ('Y', 'y')
               Displ(mask_num)%print_r(2) = .true.
            case ('Z', 'z')
               Displ(mask_num)%print_r(3) = .true.
            endselect
         enddo
         Reason = 0
      end select

      if (LEN(trim(adjustl(ch_temp))) < 1) exit ! nothing more to read
   enddo

   !pause 'read_displacement_command'
end subroutine read_displacement_command


subroutine number_atomic_mask(mask_name, Displ, mask_num)
   character(*), intent(inout) :: mask_name
   !type(Super_cell), intent(in) :: Scell  ! supercell with all the atoms as one object
   type(Displacement_analysis), dimension(:), intent(in) :: Displ ! [A] mean displacements for atoms masked
   integer, intent(in) :: mask_num
   !---------------
   integer :: j, rep_num, Reason
   character(100) :: string_part1, string_part2, string_temp, ch_mask_num

   do j = 1, mask_num
      string_temp = trim(adjustl(Displ(j)%mask_name))  ! to work with
      ! If the name repeats, add a number to it at the end:
      if ( trim(adjustl(mask_name)) == trim(adjustl(string_temp)) ) then
         ! Find if there is already a number assigned:
         call split_command_separator(trim(adjustl(string_temp)), '_', string_part1, string_part2, back=.true.)  ! below

         if (LEN(trim(adjustl(string_part2))) == 0) then ! no numnber in the name
            mask_name = trim(adjustl(mask_name))//'_1'   ! make it the first
         else  ! there is a number, add to it
            read(string_part2, *, iostat=Reason) rep_num
            if (Reason == 0) then ! read well
               write(ch_mask_num, '(i5)') rep_num+1   ! next number
               mask_name = trim(adjustl(string_part1))//'_'//trim(adjustl(ch_mask_num))   ! add this number to the name
            else  ! not a number, just underscore in the name
               mask_name = trim(adjustl(string_temp))//'_1' ! add the first number
            endif
         endif
      endif
   enddo
end subroutine number_atomic_mask



subroutine read_and_define_section(FN1, Displ, mask_num)
   integer, intent(in) :: FN1 ! file number
   !type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Displacement_analysis), dimension(:), intent(inout) :: Displ ! [A] mean displacements for atoms masked
   integer, intent(in) :: mask_num  ! mask number
   !---------------------------
   character(200) :: sting_temp, sting_temp2
   character(200), dimension(3) :: string_part
   character(100) :: ch_comm, ch_val, block1, block2
   integer :: Reason, i

   read(FN1,'(a)',IOSTAT=Reason) sting_temp
   if (Reason /= 0) then   ! couldn't read, use default
      ! No mask parameters, cannot do section, use all instead
      Displ(mask_num)%mask_name = 'all'
      return
   else  ! read something, attempt to interprete it:
      ! Define default parameters of the section:
      Displ(mask_num)%logical_and = .false.  ! to start with
      Displ(mask_num)%logical_or = .false.   ! to start with
      Displ(mask_num)%r_start = -1.0d10   ! start at -infinity
      Displ(mask_num)%r_end = 1.0d10      ! end at +infinity

      !---------------
      ! Check if there is a separator:
      string_part = ''  ! to start with

      call split_command_separator(trim(adjustl(sting_temp)), ';', string_part(1), string_part(2))  ! below
      if (LEN(trim(adjustl(string_part(2)))) /= 0) then
         sting_temp2 = string_part(2)
         call split_command_separator(trim(adjustl(sting_temp2)), ';', string_part(2), string_part(3))  ! below
      endif

      ! Read and interprete all 3 parts:
      do i = 1, 3
         block1 = '' ! to start with
         block2 = '' ! to start with
         ch_comm = '' ! to start with
         ch_val = '' ! to start with

         if (LEN(trim(adjustl(string_part(i)))) > 0) then
            sting_temp = trim(adjustl(string_part(i)))
         else  ! no text to interprete here
            if (i > 1) cycle  ! nothing else to do, skip it
         endif

         !---------------
         ! Check if there is an 'and':
         call split_command_separator(trim(adjustl(sting_temp)), 'and', ch_comm, ch_val)  ! below
         ! Check other possible ways:
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'AND', ch_comm, ch_val)  ! below
         endif
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'And', ch_comm, ch_val)  ! below
         endif
         ! If there was 'and', save two blocks:
         if (LEN(trim(adjustl(ch_comm))) /= 0) then
            Displ(mask_num)%logical_and = .true.
            block1 = ch_comm
            block2 = ch_val
         endif

         !---------------
         ! Check if there is an 'or':
         call split_command_separator(trim(adjustl(sting_temp)), 'or', ch_comm, ch_val)  ! below
         ! Check other possible ways:
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'OR', ch_comm, ch_val)  ! below
         endif
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'Or', ch_comm, ch_val)  ! below
         endif
         ! If there was 'or', save two blocks:
         if (LEN(trim(adjustl(ch_comm))) /= 0) then
            Displ(mask_num)%logical_or = .true.
            block1 = ch_comm
            block2 = ch_val
         endif

         !---------------
         ! interprete the line:
         !if (Displ(mask_num)%logical_and .or. Displ(mask_num)%logical_or) then ! interprete 2 blocks
         if ((LEN(trim(adjustl(block1))) > 0) .and. (LEN(trim(adjustl(block2))) > 0) ) then  ! interprete 2 blocks
            ! block 1:
            call identify_section_axis(Displ, block1, mask_num, 1) ! below

            ! block 2:
            call identify_section_axis(Displ, block2, mask_num, 2) ! below

         else ! interprete single line
            ! whole line:
            call identify_section_axis(Displ, sting_temp, mask_num, 1) ! below
         endif
      enddo ! i = 1,3
      !print*, 'read_and_define_section: ', Displ(mask_num)%logical_and, Displ(mask_num)%logical_or, &
      !trim(adjustl(block1))//' ', trim(adjustl(block2))
!       print*, 'sta', mask_num, Displ(mask_num)%r_start(1, :)
!       print*, 'end', mask_num, Displ(mask_num)%r_end(1, :)
!       print*, 'sta', mask_num, Displ(mask_num)%r_start(2, :)
!       print*, 'end', mask_num, Displ(mask_num)%r_end(2, :)

   endif
!    pause 'PAUSE read_and_define_section'
end subroutine read_and_define_section


subroutine identify_section_axis(Displ, read_line, mask_num, axis_ind)
   !type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Displacement_analysis), dimension(:), intent(inout) :: Displ ! [A] mean displacements for atoms masked
   character(*), intent(in) :: read_line  ! line with two command separated by a given symbol
   integer, intent(in) :: mask_num, axis_ind
   !-----------------
   character(100) :: command1, command2, ch_num
   integer :: Reason, len1, len2, axis_num, ind_larger, ind_smaller
   real(8) :: r_temp

   ! section along X:
   call split_command_separator(trim(adjustl(read_line)), 'x', command1, command2)  ! below
   if (LEN(trim(adjustl(command1))) == 0) then   ! another possible way of writing it
      call split_command_separator(trim(adjustl(read_line)), 'X', command1, command2)  ! below
   endif
   command1 = trim(adjustl(command1))
   command2 = trim(adjustl(command2))
   len1 = LEN(trim(adjustl(command1)))
   len2 = LEN(trim(adjustl(command2)))

   if ( (len1 /= 0) .or. (len2 /=0) ) then   ! section along X:
      axis_num = 1
      Displ(mask_num)%axis_ind(axis_ind) = axis_num  ! X

      ! Block #1:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command1)), '>')
      ind_smaller = INDEX(trim(adjustl(command1)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command1, Displ, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command1, Displ, mask_num, axis_ind, axis_num, 2)   ! below

      ! Block #2:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command2)), '>')
      ind_smaller = INDEX(trim(adjustl(command2)), '<')

      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command2, Displ, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command2, Displ, mask_num, axis_ind, axis_num, 2)   ! below
   endif

   ! section along Y:
   call split_command_separator(trim(adjustl(read_line)), 'y', command1, command2)  ! below
   if (LEN(trim(adjustl(command1))) == 0) then   ! another possible way of writing it
      call split_command_separator(trim(adjustl(read_line)), 'Y', command1, command2)  ! below
   endif
   command1 = trim(adjustl(command1))
   command2 = trim(adjustl(command2))
   len1 = LEN(trim(adjustl(command1)))
   len2 = LEN(trim(adjustl(command2)))
   if ( (len1 /= 0) .or. (len2 /=0) ) then ! section along Y
      axis_num = 2
      Displ(mask_num)%axis_ind(axis_ind) = axis_num  ! Y

      ! Block #1:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command1)), '>')
      ind_smaller = INDEX(trim(adjustl(command1)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command1, Displ, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command1, Displ, mask_num, axis_ind, axis_num, 2)   ! below

      ! Block #2:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command2)), '>')
      ind_smaller = INDEX(trim(adjustl(command2)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command2, Displ, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command2, Displ, mask_num, axis_ind, axis_num, 2)   ! below
   endif

   ! section along Z:
   call split_command_separator(trim(adjustl(read_line)), 'z', command1, command2)  ! below
   if (LEN(trim(adjustl(command1))) == 0) then   ! another possible way of writing it
      call split_command_separator(trim(adjustl(read_line)), 'Z', command1, command2)  ! below
   endif
   command1 = trim(adjustl(command1))
   command2 = trim(adjustl(command2))
   len1 = LEN(trim(adjustl(command1)))
   len2 = LEN(trim(adjustl(command2)))
   if ( (len1 /= 0) .or. (len2 /=0) ) then   ! section along Z
      axis_num = 3
      Displ(mask_num)%axis_ind(axis_ind) = axis_num  ! Z

      ! Block #1:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command1)), '>')
      ind_smaller = INDEX(trim(adjustl(command1)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command1, Displ, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command1, Displ, mask_num, axis_ind, axis_num, 2)   ! below

      ! Block #2:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command2)), '>')
      ind_smaller = INDEX(trim(adjustl(command2)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command2, Displ, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command2, Displ, mask_num, axis_ind, axis_num, 2)   ! below
   endif
end subroutine identify_section_axis



subroutine assign_atomic_section(ind_given, command, Displ, mask_num, ind_sec, ind_axis, start_or_end)
   !type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Displacement_analysis), dimension(:), intent(inout) :: Displ ! [A] mean displacements for atoms masked
   integer, intent(in) :: ind_given, mask_num, ind_sec, ind_axis, start_or_end
   character(*), intent(in) :: command
   !-----------------
   integer :: Reason, len1
   real(8) :: r_temp
   logical :: left_boundary

   left_boundary = .false. ! to start with

   if (ind_given > 0) then ! there is a '>' or '<'
      Reason = -1 ! to start with
      len1 = LEN(trim(adjustl(command)))

      if (ind_given == 1) then  ! it is on the left of the number
         read(command(2:len1), *, IOSTAT=Reason) r_temp
         if (start_or_end == 1) then   ! '>' it is the lower boundary
            left_boundary = .true.
         endif
      elseif (ind_given == len1) then ! it is on the right of the number
         read(command(1:len1-1), *, IOSTAT=Reason) r_temp
         if (start_or_end == 2) then   ! '<' it is the lower boundary
            left_boundary = .true.
         endif
      endif ! (ind_given == len1)
      !print*, 'r_temp', r_temp

      if (Reason == 0) then   ! read it well, use it
         if (left_boundary) then ! '>'
            Displ(mask_num)%r_start(ind_sec, ind_axis) = max ( r_temp, Displ(mask_num)%r_start(ind_sec, ind_axis) )
         else  ! '<'
            Displ(mask_num)%r_end(ind_sec, ind_axis) = min ( r_temp, Displ(mask_num)%r_end(ind_sec, ind_axis) )
         endif
      endif ! (Reason == 0)
      !print*, Displ(mask_num)%r_start(ind_sec, ind_axis), Displ(mask_num)%r_end(ind_sec, ind_axis)
   endif ! (ind_given > 0)
end subroutine assign_atomic_section


subroutine split_command_separator(read_line, separator, command1, command2, back)
   character(*), intent(in) :: read_line  ! line with two command separated by a given symbol
   character(*), intent(in) :: separator  ! separator symbol
   character(*), intent(out) :: command1  ! command #1, before separator
   character(*), intent(out) :: command2  ! command #2, after separator
   logical, intent(in), optional :: back  ! to search from the backend
   !-----------------------
   integer :: N_sep, sep_len
   logical :: back_search

   if (present(back)) then ! read what user set
      back_search = back
   else  ! by default, search from the start
      back_search = .false.
   endif

   ! Find the position of the separator:
   N_sep = 0   ! to start with
   N_sep = INDEX(read_line, separator, back=back_search) ! intrinsic
   sep_len = LEN(trim(adjustl(separator)))

   ! If separator is there, read commands:
   if (N_sep > 0) then
      command1 = read_line(1:N_sep-1)
      command2 = read_line(N_sep+sep_len:)
   else ! no separator in the string
      command1 = ''  ! undefined
      command2 = ''  ! undefined
   endif
end subroutine split_command_separator




subroutine sort_fragment_m_over_z(cur_mass, cur_q, fragment_m_over_z)
   real(8), intent(in) :: cur_mass, cur_q  ! average mass and charge of given fragment
   real(8), dimension(:), intent(inout) :: fragment_m_over_z
   !----------------------
   integer :: Nsiz, i1, i2
   real(8) :: floor_q, ceiling_q

   Nsiz = size(fragment_m_over_z)

   if (cur_q < 0.0d0) then ! neutral or negative, save here:
      fragment_m_over_z(1) = fragment_m_over_z(1) + 1
   else  ! positive ion:
      floor_q = dble(FLOOR(cur_q))
      ceiling_q = dble(ceiling(cur_q))
      i2 = dble(ANINT(cur_mass)) / ceiling_q
      if (floor_q <= 0.5d0) then ! neutrals present:
         i1 = 1
      else  ! only charges:
         i1 = dble(ANINT(cur_mass)) / floor_q
      endif
      fragment_m_over_z(i1) = fragment_m_over_z(i1) + ( ceiling_q - cur_q )
      fragment_m_over_z(i2) = fragment_m_over_z(i2) + ( cur_q - floor_q )
   endif
end subroutine sort_fragment_m_over_z


subroutine get_fragments_indices(Step, i_step, cut_r, indices)
   type(Instant_data), dimension(:), allocatable, intent(in) :: Step           ! All data and parameters at this timestep
   integer, intent(in) :: i_step ! time step
   real(8), intent(in) :: cut_r  ! [A] cut-off radius
   integer, dimension(:), allocatable, intent(inout) :: indices ! working array of indices
   real(8) a_r, dm
   integer i, j, coun, ind_i, ind_same, Na

   dm = cut_r ! set cut-off radius

   Na = size(Step(i_step)%name) ! corresponding to number of atoms
   if (allocated(indices)) deallocate(indices)
   allocate(indices(Na), source=0)

   coun = 0 ! to start with
   do i = 1, Na ! Check each atom - to which fragment it belongs
      if (indices(i)==0) then ! this atom is not allocated within any fragment yet
         coun = coun + 1   ! next fragment
         indices(i) = coun ! atom is within this fragment
      endif
      ind_i = indices(i)
      do j = 1, Na ! check if there are other atoms within the same fragment
         if (i /= j) then
            call shortest_distance(Step(i_step)%S, Step(i_step)%Supce, i, j, a_r)
            if (a_r < dm) then ! this atom is within the same fragment as the atom 'i'
               if (indices(j) == 0) then ! this atom has not been allocated to a fragment yet:
                  indices(j) = indices(i) ! mark this atoms with the same index
               else if (indices(j) == indices(i)) then ! this atom is already within this fragment
                  ! do nothing, it's already been done
               else ! it turns out, it's a part of the same fragment, who could've guessed?!
                  ! change indices of the second fragment to the first fragment, since it's actually the same:
                  ind_same = indices(j)
                  ind_i = MIN(ind_i,ind_same)
                  where(indices == ind_same) indices = ind_i ! renumber these atoms to this fragment
               endif
            endif
         endif
      enddo
   enddo

   !Check that there are no empty fragments left:
   do i = 1, maxval(indices(:))
      j = count(indices(:) == i)	! that's how many atoms are in this fragment
      do while (j < 1)
         where(indices >= i) indices = indices - 1	! shift them by one
         if (i >= maxval(indices(:))) exit	! all indices are shifted down to here, nothing to do more
         j = count(indices(:) == i)	! that's how many atoms are in this fragment
!       print*, 'i,j', i, j, maxval(indices(:))
      enddo
   enddo

end subroutine get_fragments_indices



subroutine shortest_distance(S, Supce, i1, j1, a_r)
   real(8), dimension(:,:), intent(in) :: S  ! relative coordinates
   real(8), dimension(:,:), intent(in) :: Supce  ! supercell
   integer, intent(in) :: i1, j1 ! atomic numbers
   real(8), intent(out) ::  a_r	! [A] shortest distance between the two atoms within supercell with periodic boundaries
   real(8) x, y, z, zb(3), x0, y0, z0, r1
   integer i, j, k, ik

   x = 0.0d0
   y = 0.0d0
   z = 0.0d0
   if (i1 == j1) then ! it's the same atom:
      a_r = 0.0d0
   else
      ! For the case of periodic boundaries:
      do ik = 1,3
         x = x + (S(ik,i1) - S(ik,j1))*Supce(ik,1)
         y = y + (S(ik,i1) - S(ik,j1))*Supce(ik,2)
         z = z + (S(ik,i1) - S(ik,j1))*Supce(ik,3)
      enddo ! ik
      a_r = DSQRT(x*x + y*y + z*z)

      do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
         ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell:
         ! periodic boundary conditions
         zb(1) = dble(i)
         do j =-1,1
            zb(2) = dble(j)
            do k = -1,1
               zb(3) = dble(k)
               x0 = 0.0d0
               y0 = 0.0d0
               z0 = 0.0d0
               do ik = 1,3
                  x0 = x0 + (S(ik,i1) - S(ik,j1) + zb(ik))*Supce(ik,1)
                  y0 = y0 + (S(ik,i1) - S(ik,j1) + zb(ik))*Supce(ik,2)
                  z0 = z0 + (S(ik,i1) - S(ik,j1) + zb(ik))*Supce(ik,3)
               enddo ! ik
               r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
               if (r1 <= a_r) then
                  x = x0
                  y = y0
                  z = z0
                  a_r = r1
               endif
            enddo ! k
         enddo ! j
      enddo ! i
   endif ! i1 = j1
end subroutine shortest_distance




subroutine read_time_step(FN, FN2, Step, time_print)
   integer, intent(in) :: FN, FN2    ! file number with data for time grid and supercell parameters
   type(Instant_data), dimension(:), allocatable, intent(inout) :: Step           ! All data and parameters at this timestep
   real(8), intent(in) :: time_print
   integer :: Reason, Nsiz, Nat, i, j, counter
   real(8) eps, temp, cur_t, Vol
   character(500) :: read_line
   logical :: there_is_q

   there_is_q = .false.
   eps = 1.0d-6
   call Count_lines_in_file(FN, Nsiz)  ! below
   Nsiz = Nsiz - 2  ! skip first two lines with comments
   read(FN,*,IOSTAT=Reason) ! skip first two lines with comments
   read(FN,*,IOSTAT=Reason)
   cur_t = -1.0d20   ! to start with
   counter = 0
   do i = 1, Nsiz
      read(FN,*,IOSTAT=Reason) temp ! time step
      if (temp > cur_t-1.0d-6) then
         cur_t = temp + time_print
         counter = counter + 1   ! count printable steps
      endif
      !print*, 'temp=', i, temp
   enddo
   rewind(FN)

   ! Knowing how many steps, allocate the data:
   allocate (Step(counter))

   ! Read temperatures:
   read(FN,*,IOSTAT=Reason) ! skip first two lines with comments
   read(FN,*,IOSTAT=Reason)
   cur_t = -1.0d20   ! to start with
   counter = 0
   do i = 1, Nsiz
      read(FN,*,IOSTAT=Reason) temp ! time step
      if (temp > cur_t-1.0d-6) then
         cur_t = temp + time_print
         counter = counter + 1   ! count printable steps
         backspace(FN) ! to reread the same line, but with supercell parameters:
         read(FN,*,IOSTAT=Reason) Step(counter)%Tim, Vol, Step(counter)%Supce(:,:)
         if (Reason .LT. 0) exit
         !print*, 'Supcell: ', Step(counter)%Tim, Vol, Step(counter)%Supce(:,:)

         ! Read the atomic coordinates (in XYZ format):
         read(FN2,*,IOSTAT=Reason) Nat
         read(FN2,*,IOSTAT=Reason) ! skip comment line
         ! Allocate the coordinates array:
         if (.not.allocated(Step(counter)%name)) then
            allocate(Step(counter)%name(Nat))
            if (.not.allocated(Step(counter)%KOA)) allocate(Step(counter)%KOA(Nat))
            if (.not.allocated(Step(counter)%R)) allocate(Step(counter)%R(3,Nat))
            if (.not.allocated(Step(counter)%S)) allocate(Step(counter)%S(3,Nat))
            if (.not.allocated(Step(counter)%M)) allocate(Step(counter)%M(Nat))
            if (.not.allocated(Step(counter)%q)) allocate(Step(counter)%q(Nat), source = -1.0d10)
            ! Check if there is charge state provided:
            there_is_q = .true.  ! mark the first step

            read(FN2,'(a)',IOSTAT=Reason) read_line
            read(read_line,*,IOSTAT=Reason) Step(counter)%name(1), Step(counter)%R(:,1), Step(counter)%q(1)
            if (Reason /= 0) then   ! there is no
               there_is_q = .false. ! there is no q in the file
               if (allocated(Step(counter)%q)) deallocate(Step(counter)%q)
            endif
            backspace(FN2) ! to reread the line below
         endif

         do j = 1, Nat
            read(FN2,'(a)',IOSTAT=Reason) read_line
            if (there_is_q) then ! there is q
               read(read_line,*,IOSTAT=Reason) Step(counter)%name(j), Step(counter)%R(:,j), Step(counter)%q(j)
               if (Reason /= 0) then
                  there_is_q = .false. ! there is no q in the file
                  read(read_line,*,IOSTAT=Reason) Step(counter)%name(j), Step(counter)%R(:,j)
               endif
!                ! Test:
!                print*, j, Step(counter)%q(j)
            else  ! there is no q
               read(read_line,*,IOSTAT=Reason) Step(counter)%name(j), Step(counter)%R(:,j)
!                ! Test:
!                print*, j, 'no q'
            endif
         enddo

         ! Get the relative coordinates too:
         call coords_abs_to_rel(Step(counter)%R, Step(counter)%Supce, Step(counter)%S) ! below

      else ! skip-read the XYZ block:
         read(FN2,*,IOSTAT=Reason) Nat
         read(FN2,*,IOSTAT=Reason) ! skip comment line
         do j = 1, Nat
            read(FN2,*,IOSTAT=Reason) ! skip line
         enddo
      endif
   enddo
end subroutine read_time_step

 

subroutine coords_abs_to_rel(R, Supce, S)
   real(8), dimension(:,:), intent(in) :: R ! absolute coordinates
   real(8), dimension(3,3), intent(in) :: Supce
   real(8), dimension(:,:), intent(out) :: S ! relative coordinates
   ! -----------------------
   real(8) v(3), dsupce(3,3)
   integer i, ik, N
   N = size(R,2)
   !Relative velocities:
   call Invers_3x3(Supce, dsupce, 'coords_abs_to_rel') ! below
   do i = 1, N
      v = 0.0d0
      do ik = 1,3
         v(:) = v(:) + R(ik,i)*dsupce(ik,:)
      enddo ! ik
      S(:,i) = v(:)
   enddo
end subroutine coords_abs_to_rel


! This subroutine calculates the inverse of a 3x3 matrix:
subroutine Invers_3x3(A, InvA, ref_sub)
   REAL(8), DIMENSION(3,3), INTENT(in) :: A ! input matrix
   REAL(8), DIMENSION(3,3), INTENT(out) :: InvA ! output matrix, inverse of A(3x3)
   character(*), intent(in) :: ref_sub  ! from which subroutine it is called
   real(8) detA ! determinant of A
   call Det_3x3(A,detA) ! find determinant of A, below
   InvA(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)    ! A
   InvA(1,2) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3)) ! B
   InvA(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)    ! C
   InvA(2,1) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3)) ! D
   InvA(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)    ! E
   InvA(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1)) ! F
   InvA(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)    ! G
   InvA(3,2) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2)) ! H
   InvA(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)    ! K
   if (detA .LE. 0.0e0) then
      print*, 'Det A = 0 in Invers_3x3, called from '//trim(adjustl(ref_sub))
      InvA = 1.0d29
   else
      InvA = InvA/detA
   endif
   if (ABS(InvA(1,1)) .GT. 1.0e30) print*, 'InvA(1,1) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(1,1)
   if (ABS(InvA(1,2)) .GT. 1.0e30) print*, 'InvA(1,2) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(1,2)
   if (ABS(InvA(1,3)) .GT. 1.0e30) print*, 'InvA(1,3) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(1,2)
   if (ABS(InvA(2,1)) .GT. 1.0e30) print*, 'InvA(2,1) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(2,1)
   if (ABS(InvA(2,2)) .GT. 1.0e30) print*, 'InvA(2,2) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(2,2)
   if (ABS(InvA(2,3)) .GT. 1.0e30) print*, 'InvA(2,3) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(2,3)
   if (ABS(InvA(3,1)) .GT. 1.0e30) print*, 'InvA(3,1) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(3,1)
   if (ABS(InvA(3,2)) .GT. 1.0e30) print*, 'InvA(3,2) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(3,2)
   if (ABS(InvA(3,3)) .GT. 1.0e30) print*, 'InvA(3,3) = infinity in Invers_3x3, called from '//trim(adjustl(ref_sub)), A(3,3)
end subroutine Invers_3x3


subroutine Det_3x3(A,detA)
   REAL(8), DIMENSION(3,3), INTENT(in) :: A ! input matrix
   REAL(8), INTENT(out) :: detA ! determinant of A, output
   detA = A(1,1)*( A(2,2)*A(3,3) - A(2,3)*A(3,2) ) - &
          A(1,2)*( A(2,1)*A(3,3) - A(2,3)*A(3,1) ) + &
          A(1,3)*( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
   if (detA .LE. 0.0d0) print*, 'Det A = 0 in Det_3x3'
end subroutine Det_3x3 ! checked!




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

END PROGRAM Fragmentation
