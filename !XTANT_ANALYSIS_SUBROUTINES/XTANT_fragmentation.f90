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
   character(3), dimension(:), allocatable :: name    ! element name
   real(8), dimension(:,:), allocatable :: R    ! absolute coordinates for all atoms [A]
   real(8), dimension(:,:), allocatable :: S    ! relative coordinates for all atoms
   real(8), dimension(3,3) :: Supce
   real(8), dimension(:), allocatable :: M    ! Atomic mass [a.m.u.]
   real(8), dimension(:), allocatable :: q    ! Atomic charge [electron charge]
end type Instant_data

character(400) :: error_message, File_XYZ, File_supce, read_line
character(200) :: File_fragments, Gnu_script, Gnu_file, command, File_m_over_z
character(100) :: ChemFormula
character(32), dimension(10) :: char_var
character(10) :: temp_ch
character(1) :: path_sep

type(Instant_data), dimension(:), allocatable :: Step           ! All data and parameters at this timestep

real(8) :: cut_r, time_dt, cur_t, dt, time_print, cur_mass, cur_q

character(3), dimension(:), allocatable :: at_short_names ! name of the element
real(8), dimension(:), allocatable :: at_masses ! mass of each element [a.m.u.]

integer, dimension(:), allocatable :: indices
real(8), dimension(:), allocatable :: fragment_masses, fragment_m_over_z
real(8), dimension(:), allocatable :: mass_grid
real(8), dimension(:,:), allocatable :: output_fragment_array, output_m_over_z_array
integer :: FN1, FN2, FN_out, FN_out1, FN_out_mz, FN_out2 ! file number
integer :: INFO, Reason, i, j, Tsiz, Nat, existing_elem, at_num, cur_j, iret
logical :: read_well

call Path_separator(path_sep)  ! Objects_and_types

! Set defaults:
INFO = 0
FN_out = 1000
FN_out_mz = 1001
FN1 = 9999
FN2 = 9998


File_supce = 'OUTPUT_supercell.dat'	! defaul name of XTANT out file with supercell sizes
File_XYZ = 'OUTPUT_atomic_coordinates.xyz'
File_fragments = 'OUT_fragments_spectrum.dat'   ! default name, for start
File_m_over_z = 'OUT_m_over_z_spectrum.dat'     ! default name, for start

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
Nat = size(Step(1)%R(1,:)) ! number of atoms in the simulation box
do i = 1, Nat
   ! Find if this element repeats, or a new one:
   existing_elem = INDEX(trim(adjustl(ChemFormula)), trim(adjustl(Step(1)%name(i))))   ! intrinsic function
   if (existing_elem == 0) then ! new element
      ! Add it to the list (construct full chemical formula):
      ChemFormula = trim(adjustl(ChemFormula))//trim(adjustl(Step(1)%name(i)))
   endif
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
! Output file:
open (unit=FN_out, file=trim(adjustl(File_fragments)))
open (unit=FN_out_mz, file=trim(adjustl(File_m_over_z)))

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
enddo

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
   command = trim(adjustl(Gnu_script))//'.cmd'
   iret = system(command)
else ! linux:
   !call system(trim(adjustl(Gnu_script))//'.sh')
!    command = "./OUTPUT_Gnuplot_all.sh"
   command = trim(adjustl(Gnu_script))//'.sh'
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
   command = trim(adjustl(Gnu_script))//'.cmd'
   iret = system(command)
else ! linux:
   !call system(trim(adjustl(Gnu_script))//'.sh')
!    command = "./OUTPUT_Gnuplot_all.sh"
   command = trim(adjustl(Gnu_script))//'.sh'
   iret = system(command)
endif



2012 continue   ! to exit the program
STOP
!---------------------
 contains


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
   call Det_3x3(A,detA) ! find determinant of A, see above
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
