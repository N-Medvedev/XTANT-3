PROGRAM Atomic_data_analysis
! Compilation:
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs XTANT_atomic_data_analysis.f90 -o XTANT_atomic_data_analysis.exe /link /stack:9999999999 
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays XTANT_atomic_data_analysis.f90 -o XTANT_atomic_data_analysis.exe /link /stack:9999999999 
!<===========================================

use Universal_constants

implicit none

character(200) :: File_name_coords, File_name_NRGs, File_name_T, File_name_supce, File_parameters, File_name_out, command, File_name
character(1) :: path_sep
character(10) :: temp_ch
integer :: FN, FN_out, FN_coord, FN_supce, FN_NRG, FN_T, FN_Param
integer :: NOMP, Nat, Reason, open_status
real(8), dimension(:,:,:), allocatable :: R, V, SR, SV
logical :: read_well, file_opened, do_PCF, do_autocorr, do_Cv

!<===========================================
! get the file names:
FN_coord = 9999
FN_supce = 9998
FN_out = 9997
FN_NRG = 9996
FN_T = 9995
FN_Param = 9994
File_name_coords = 'OUTPUT_coordinates_and_velosities.dat'	! defaul name of XTANT out file with coordinates and velosities
File_name_supce = 'OUTPUT_supercell.dat'	! defaul name of XTANT out file with supercell sizes
File_name_NRGs = 'OUTPUT_energies.dat'
File_name_T = 'OUTPUT_temperatures.dat'

call Path_separator(path_sep)
! Get the name of file with atomic parameters:
File_name = 'Temp.txt'
! Get the names of all data files in the folder using system commands:
if (path_sep .EQ. '\') then	! if it is Windows
   command = 'dir *.txt /b > '//trim(adjustl(File_name))  ! 'Temp.txt'
else
   command = "ls -t | grep '.txt' > "//trim(adjustl(File_name)) ! "Temp.txt"
endif
call system(trim(adjustl(command))) ! execute the command to save file names in the temp file
FN = 1100
open(UNIT=FN, file=trim(adjustl(File_name)), iostat=open_status, action='read')
File_parameters = ''
do
   read(FN,'(A)',IOSTAT=Reason) temp_ch
   if (Reason /= 0) exit
   if (trim(adjustl(temp_ch(1:1))) == '!') then	! we found the required name
       backspace(FN)
       read(FN,'(A)',IOSTAT=Reason) File_parameters
      exit
   endif
enddo
close(FN, status='delete') ! temp file is not needed anymore, erase it

open (unit=FN_coord, file=trim(adjustl(File_name_coords)), status = 'old', readonly)
open (unit=FN_supce, file=trim(adjustl(File_name_supce)), status = 'old', readonly)
open (unit=FN_NRG, file=trim(adjustl(File_name_NRGs)), status = 'old', readonly)
open (unit=FN_T, file=trim(adjustl(File_name_T)), status = 'old', readonly)
open (unit=FN_Param, file=trim(adjustl(File_parameters)), status = 'old', readonly)

!<===========================================
! get the parameters:
call get_add_data(do_PCF, do_autocorr, do_Cv, read_well)	! below

#ifdef OMP_inside
   NOMP = omp_get_max_threads()	! number of processors available by default
   call omp_set_dynamic(0)	! set static number of threads
   call omp_set_num_threads(NOMP)	! define the number of threads for openmp = number of cores
#else ! if you set to use OpenMP in compiling: 'make OMP=no'
   NOMP = 1	! unparallelized by default
#endif

!<===========================================
! Perform analyses:
! 1) Pair correlation function, if needed:
if (do_PCF) then
   File_name_out =  'OUTPUT_PCF.dat'
   open (unit=FN_out, file=trim(adjustl(File_name_out))) 
   call Get_pair_correlation_function(FN_coord, FN_supce, FN_out, do_PCF)	! below
   close(FN_out)
endif

! 2) Autocorrelation function (phonon spectrum):
if (do_autocorr) then
   File_name_out =  'OUTPUT_Autocorrelation.dat'
   open (unit=FN_out, file=trim(adjustl(File_name_out))) 
   call Get_autocorrelation(FN_coord, FN_out, do_autocorr)	! below
   close(FN_out)
endif

! 3) Atomic heat capacity:
if (do_Cv) then
   File_name_out =  'OUTPUT_atomic_Cv.dat'
   open (unit=FN_out, file=trim(adjustl(File_name_out))) 
   call Get_heat_capacity(FN_NRG, FN_T, FN_Param, FN_out, do_Cv)	! below
   close(FN_out)
endif

!<===========================================
! Program ends:
inquire(file=trim(adjustl(File_name_coords)),opened=file_opened)
if (file_opened) close(FN_coord)
inquire(file=trim(adjustl(File_name_supce)),opened=file_opened)
if (file_opened) close(FN_supce)
inquire(file=trim(adjustl(File_name_NRGs)),opened=file_opened)
if (file_opened) close(FN_NRG)
inquire(file=trim(adjustl(File_name_T)),opened=file_opened)
if (file_opened) close(FN_T)
inquire(file=trim(adjustl(File_parameters)),opened=file_opened)
if (file_opened) close(FN_Param)
inquire(file=trim(adjustl(File_name_out)),opened=file_opened)
if (file_opened) close(FN_out)


!<===========================================
 contains

! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Main subroutines
 
! 1) Pair correlation function, if needed:
subroutine Get_pair_correlation_function(FN_coord, FN_supce, FN_out, do_PCF)	! below
   integer, intent(in) :: FN_coord, FN_supce, FN_out	! file numbers: input and output files
   logical, intent(in) :: do_PCF
   !--------------------------------------
   real(8), dimension(:,:), allocatable :: PCF	! pair correlation function
   real(8), dimension(:,:), allocatable :: Coord(:,:), V(:,:), SR(:,:), SV(:,:)
   real(8), dimension(3,3) :: supce
   real(8) r, dr, a_r, tim, Vol
   integer n, i, m, k, j, Reason, count_lines, count_lines2
   
   m = 500		! grid points
   dr = 0.05d0	! [nm] grid step
   
   if (do_PCF) then	! if user needs it only
      print*, 'Calculating Pair Correlation Function...'
   
      if (.not. allocated(PCF)) then
         allocate(PCF(2,m))
         r = 0.0d0	! starting point
         do i = 1,m
            r = r + dr 
            PCF(1,i) = r
         enddo
      else
         m = size(PCF,2)
      endif
      PCF(2,:) = 0.0d0	! to start with
      
      ! how many atoms we had:
      call count_Nat(FN_coord, n)
      print*, 'Nat=', n
      allocate(Coord(n,3))
      allocate(V(n,3))
      allocate(SR(n,3))
      allocate(SV(n,3))
      
      count_lines = 0
      count_lines2 = 0
      ! skip comments lines:
      read(FN_supce,*,IOSTAT=Reason)  
      call read_file(Reason, count_lines, read_well)
      if (.not.read_well) goto 2012
      read(FN_supce,*,IOSTAT=Reason)  
      call read_file(Reason, count_lines, read_well)
      if (.not.read_well) goto 2012
      
      do 	! time points
         ! Read supercell parameters:
         read(FN_supce,*,IOSTAT=Reason)  tim, Vol, supce(:,:)
         call read_file(Reason, count_lines, read_well)
         if (.not.read_well) goto 2012
!          print*, 'SUPCE',  tim, Vol
         
         ! Read atomic coordinates at this time instant:
         if (count_lines2 > 0) then	! skip 2 empty lines between 2 time points
             read(FN_coord,*,IOSTAT=Reason)
             if (.not.read_well) goto 2012
             read(FN_coord,*,IOSTAT=Reason)
             if (.not.read_well) goto 2012
         endif
         do i = 1, n	! read all atoms
            read(FN_coord,*,IOSTAT=Reason) Coord(i,:), V(i,:), SR(i,:), SV(i,:)
            call read_file(Reason, count_lines2, read_well)
            if (.not.read_well) goto 2012
         enddo
         
!$omp PARALLEL private(i,j,a_r,k)
!$omp do schedule(dynamic)
         do i = 1, n	! trace all atoms
            do j = 1, n	! trace all neighbours for PCF
               if (i .NE. j) then
                  call shortest_distance(Supce, SR, i, j, a_r)
                  if (a_r .GE. PCF(1,m)) then
                     k = m
                  else
                     call Find_in_monotonous_2D_array(PCF, a_r, 1, k)	! module "Little_subroutines"
                  endif
                  PCF(2,k) = PCF(2,k) + 1.0d0
               endif
            enddo ! j
         enddo ! i
!$omp end do
!$omp do schedule(dynamic)
         do k = 1,m	! all points of the PCF
            PCF(2,k) = PCF(2,k)/(PCF(1,k)*PCF(1,k)) 
         enddo
!$omp end do
!$omp end parallel
         PCF(2,:) = PCF(2,:)/(4.0d0*g_Pi*dr)*Vol/dble(n*n) ! normalizing per volume
         ! Write the PCF out:
         do i = 1, m
            write(FN_out, '(f25.16,es25.16)') PCF(1,i), PCF(2,i)
         enddo
         write(FN_out, '(a)') ''
         write(FN_out, '(a)') ''
         
      enddo ! tim
      
2012  deallocate(PCF, Coord, V, SR, SV)
   else
      print*, 'Pair Correlation Function is not required'
   endif	! if (do_PCF)
end subroutine Get_pair_correlation_function



! 2) Autocorrelation function (phonon spectrum):
subroutine Get_autocorrelation(FN_coord, FN2,  do_autocorr)	! below
   integer, intent(in) :: FN_coord, FN2	! file numbers: input and output files
   logical, intent(in) ::  do_autocorr
   integer :: N, count_lines2, Reason, i
   logical :: read_well
   real(8), dimension(:,:), allocatable :: Coord(:,:), V(:,:), SR(:,:), SV(:,:)
   
   if ( do_autocorr) then	! if user needs it only
      print*, 'Calculating Autocorrelation Function...'
      
   else
      print*, 'Autocorrelation Function is not required'
   endif
end subroutine Get_autocorrelation



! 3) Atomic heat capacity:
subroutine Get_heat_capacity(FN_NRG, FN_T, FN_Param, FN_out, do_Cv)	! below
   integer, intent(in) :: FN_NRG, FN_T, FN_Param, FN_out	! file numbers: input and output files
   logical, intent(in) :: do_Cv
   logical :: read_well
   integer :: N, count_lines, count_lines2, Reason, i, j, k, N_per
   real(8) :: temp1, temp2, temp3, dt, Tmean, Ekinmean, E2mean, Cv, dble_N, g_kb_J2, T2, Mass
   real(8), dimension(:), allocatable :: Tim, Ta, Ekin
      
   if (do_Cv) then	! if user needs it only
      print*, 'Calculating Atomic Heat Capacity...'
      
      ! how many time grid points we had:
      call Count_lines_in_file(FN_T, N, skip_lines=2)
      print*, 'Time points=', N
      allocate(Tim(N))
      allocate(Ta(N))
      allocate(Ekin(N))
      dble_N = dble(N)
      g_kb_J2 = g_kb_J*g_kb_J
      
      call get_atomic_mass(FN_Param, Mass)	! below
!        print*, 'Average mass=', Mass
!        pause
      
      ! Skip 2 lines with titles:
      read(FN_NRG,*)
      read(FN_NRG,*)
      read(FN_T,*)
      read(FN_T,*)
      count_lines = 0
      count_lines2 = 0
      do i = 1, N	! time points
         ! Read kinetic energy
         read(FN_NRG,*,IOSTAT=Reason) Tim(i), temp1, temp2, temp3, Ekin(i)
         call read_file(Reason, count_lines, read_well)
         ! Read temperature
         read(FN_T,*,IOSTAT=Reason) temp1, temp2, Ta(i)
         if (.not.read_well) then
!             print*, 'Could not read file with energies, line ', count_lines
            exit
         endif
         call read_file(Reason, count_lines2, read_well)
         if (.not.read_well) then
!             print*, 'Could not read file with temperatures, line ', count_lines2
            exit
         endif
!          print*, Tim(i), Ta(i), Ekin(i)
      enddo ! tim
      dt = Tim(2) - Tim(1)	! [fs] time step
      N_per = 100.0d0/dt	! number of points for 100 fs period
      if (N_per +1 > N) N_per = N - 1
      
      do i = N_per+1, N	! all time points starting from the minimum period needed to estimate average values
         ! Get mean parameters:
         Tmean = SUM(Ta(i-N_per:i))/dble(N_per)	! [K]
         T2 = Tmean*Tmean	! [K^2]
         Ekinmean = SUM(Ekin(i-N_per:i))/dble(N_per)	! [eV]
         E2mean = SUM(Ekin(i-N_per:i)*Ekin(i-N_per:i))/dble(N_per)	! [eV^2]
         ! Heat capacity:
!          Cv = 1.5d0*dble_N*dble_N*g_kb_J*g_kb_J2*T2/(dble_N*g_kb_J2*T2 - 2.0d0*g_e*g_e/3.0d0*(E2mean - Ekinmean*Ekinmean) )/Mass
         Cv = 1.5d0*g_kb_J / (1.0d0 - (E2mean - Ekinmean*Ekinmean)*g_kb*g_kb/(1.50d0*T2) ) / Mass
         write(*,'(i,f,f,f,f)') i, Tim(i), Cv, Cv*Mass/(1.5d0*g_kb_J)
         write(FN_out,'(f,f)') Tim(i), Cv	! [fs], [J/(K*g)]
      enddo ! i = N_per, N
      
2013  deallocate(Tim, Ta, Ekin)

   else
      print*, 'Atomic Heat Capacity is not required'
   endif ! if (do_autocorr)
end subroutine Get_heat_capacity








! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Supporting subroutines

subroutine get_atomic_mass(FN_Param, Mass)
   integer, intent(in) :: FN_Param
   real(8), intent(out) :: Mass	! [g]
   real(8) :: rho, at_rho
   integer :: Reason, count_lines
   logical :: read_well
   character(34) :: temp
   count_lines = 0
   do
      read(FN_Param,'(A)',IOSTAT=Reason) temp
!       print*, 'temp :: ', trim(adjustl(temp(1:24)))
      call read_file(Reason, count_lines, read_well)
      if (.not.read_well) goto 2015
      if (trim(adjustl(temp(1:24))) == 'Density of the material') exit
   enddo
   backspace(FN_Param)
   read(FN_Param,'(56X,A)',IOSTAT=Reason) temp
   read(temp,*) rho
   read(FN_Param,'(45X,A)',IOSTAT=Reason) temp
   read(temp,*) at_rho
!    print*, 'rho=', rho
!    print*, 'at_rho=', at_rho
   Mass = rho/at_rho	! [g]
2015 continue
end subroutine get_atomic_mass



subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then
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

subroutine Find_in_monotonous_2D_array(Array, Value0, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2

   N = size(Array,2)
   i_1 = 1
   val_1 = Array(Indx,i_1)
   i_2 = N
   val_2 = Array(Indx,i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(Indx,i_cur)

   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_2D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
        pause 'STOPPED WORKING...'
   else
       if (Value0 .LT. Array(Indx,1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(Indx,N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(Indx,i_cur)) .AND. (Value0 .LE. Array(Indx,i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                endif
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_2D_array', coun
                    write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_2D_array




subroutine shortest_distance(Supce, S, i1, j1, a_r, x1, y1, z1)
   real(8), dimension(:,:), intent(in) :: Supce	! supercell size
   real(8), dimension(:,:), intent(in) :: S		! relative coordinates of atoms in the supercell
   integer, intent(in) :: i1, j1 ! atomic numbers
   real(8), intent(out) ::  a_r	! [A] shortest distance between the two atoms within supercell with periodic boundaries
   real(8), intent(out), optional :: x1, y1, z1		! [A] projections of the shortest distance
   real(8) x, y, z, zb(3), r, x0, y0, z0, r1
   integer i, j, k, ik

   x = 0.0d0
   y = 0.0d0
   z = 0.0d0
  if (i1 .EQ. j1) then ! it's the same atom:
     a_r = 0.0d0
     if (present(x1) .and. present(y1) .and. present(z1)) then
        x1 = x
        y1 = y
        z1 = z
     endif
  else
   ! For the case of periodic boundaries:
   do ik = 1,3
      x = x + (S(i1,ik) - S(j1,ik))*supce(ik,1)
      y = y + (S(i1,ik) - S(j1,ik))*supce(ik,2)
      z = z + (S(i1,ik) - S(j1,ik))*supce(ik,3)
   enddo ! ik
   a_r = SQRT(x*x + y*y + z*z)
   if (present(x1) .and. present(y1) .and. present(z1)) then
      x1 = x
      y1 = y
      z1 = z
   endif
   do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
      ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell: 
      ! periodic boundary conditions.
      zb(1) = dble(i)
      do j =-1,1
         zb(2) = dble(j)
         do k = -1,1
            zb(3) = dble(k)
            x0 = 0.0d0
            y0 = 0.0d0
            z0 = 0.0d0
            do ik = 1,3
               x0 = x0 + (S(i1,ik) - S(j1,ik) + zb(ik))*supce(ik,1)
               y0 = y0 + (S(i1,ik) - S(j1,ik) + zb(ik))*supce(ik,2)
               z0 = z0 + (S(i1,ik) - S(j1,ik) + zb(ik))*supce(ik,3)
            enddo ! ik
            r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
            if (r1 < a_r) then
               x = x0
               y = y0
               z = z0
               a_r = r1
               if (present(x1) .and. present(y1) .and. present(z1)) then
                  x1 = x
                  y1 = y
                  z1 = z
               endif
            endif
         enddo ! k
      enddo ! j
   enddo ! i
  endif ! i1 = j1
end subroutine shortest_distance



subroutine count_Nat(FN, n)
   integer, intent(in) :: FN	! file to read from
   integer, intent(out) :: n	! number of atoms
   integer :: i, Reason
   character(6) temp
   i = 0
   do
      i = i + 1
      read(FN,'(a)',IOSTAT=Reason) temp
      if ((Reason /= 0) .or. (len_trim(adjustl(temp)) == 0)) exit	! we reached a break point
   enddo
   rewind(FN)
   n = i-1	! number of atoms
end subroutine count_Nat



 subroutine read_file(Reason, i, read_well, Error_descript)
   integer, intent(in) :: Reason    ! description of error
   integer, intent(inout) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   character(*), intent(out), optional :: Error_descript ! what is wrong
   character(16) chi
   i = i + 1    ! it's next line
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       if (present(Error_descript)) then
          write(chi, '(i12)') i
          write(Error_descript,'(a,a,a)') 'Problem reading input file in line ', trim(adjustl(chi)), ' ,wrong type of variable'
       endif
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       if (present(Error_descript)) then
          write(chi, '(i12)') i
          write(Error_descript,'(a,a,a)') 'Problem reading input file in line ', trim(adjustl(chi)), ' ,unexpected END of file'
       endif
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file
 
 
 
! Reads additional data from the command line passed along:
subroutine get_add_data(do_PCF, do_autocorr, do_Cv, read_well)
   logical, intent(inout) :: do_PCF, do_autocorr, do_Cv
   logical, intent(inout) :: read_well  ! did we read ok?
   !---------------------------------------
   integer, dimension(:), allocatable :: stop_markers ! how many different optionas are passed?
   integer :: i, count_dash
   character(500) :: string, char1
   ! Set defaults:
   do_PCF = .false.
   do_autocorr = .false.
   do_Cv = .false.
   read_well = .true. ! to start with
   
   string = '' ! to start with empty line
   ! Read all arguments passed to the code:
   do i = 1, iargc() ! for as many arguments as we passed
      call getarg(i, char1)
      string = trim(adjustl(string))//' '//trim(adjustl(char1)) ! collect them all into one line
   enddo
   print*, 'Interpreting line: ', trim(adjustl(string))
   
   ! Find out how many parameters user provided:
   call parse_add_param_line(string, stop_markers)
   
   ! Interpret them, for each command - read it, interpret and execute:
   do i = 1, size(stop_markers)
      if (i < size(stop_markers)) then
!          print*, 'READ:', trim(adjustl(string(stop_markers(i):stop_markers(i+1)-1)))
         call react_to_command_passed(trim(adjustl(string(stop_markers(i):stop_markers(i+1)-1))), do_PCF, do_autocorr, do_Cv) ! below
      else
         call react_to_command_passed(trim(adjustl(string(stop_markers(i):LEN(trim(adjustl(string)))+1))), do_PCF, do_autocorr, do_Cv) ! below
!          print*, 'READ:', trim(adjustl(string(stop_markers(i):LEN(trim(adjustl(string)))+1)))
      endif
   enddo
  
2017   if (.not.read_well) then
      write(*,'(a)') '***************************************************************************'
      write(*,'(a)') 'Could not interpret the passed string: '//trim(adjustl(string))
   endif
end subroutine get_add_data


subroutine react_to_command_passed(string, do_PCF, do_autocorr, do_Cv)
   character(*), intent(in) :: string	! line with the command to analyze
   logical, intent(inout) :: do_PCF, do_autocorr, do_Cv	! do the following parameters or not
   ! Read if user changed the defaults:
   select case (trim(adjustl(string(2:))))
   case ('PCF', 'pcf')
      do_PCF = .true.
   case ('Autocorrelation', 'Autocorr', 'AUTOCORR', 'AUTO', 'Auto', 'auto', 'autocorr')
      do_autocorr = .true.
   case ('Cv', 'CV', 'cv', 'Heat_capacity')
      do_Cv = .true.
   case ('help', 'Help', 'HELP')
      print*, 'The line '//trim(adjustl(string))//' could not be interpreted.'
      print*, 'The following parameters are supported:'
      print*, '* /PCF - which stands for pair correlation function'
      print*, '* /AUTO - which stands for autocorrelation function'
      print*, '* /Cv - which stands for atomic heat capacity'
   end select
end subroutine react_to_command_passed



pure subroutine parse_add_param_line(string, stop_markers)
   character(*), intent(in) :: string ! the line we read from the 
   integer, dimension(:), allocatable, intent(inout) :: stop_markers ! how many different optionas are passed?
   integer :: leng, i, count_dash
   character(LEN(trim(adjustl(string)))) :: string_cur
   character(1) :: separator
   separator = '/'
   string_cur = trim(adjustl(string))
   leng = LEN(trim(adjustl(string_cur)))
   count_dash = 0
   do i = 1,leng
      if (trim(adjustl(string_cur(i:i))) == separator) then
         count_dash = count_dash + 1
      endif
   enddo
   if (allocated(stop_markers)) deallocate(stop_markers)
   allocate(stop_markers(count_dash))
   count_dash = 0
   do i = 1,leng
      if (trim(adjustl(string_cur(i:i))) == separator) then
         count_dash = count_dash + 1  ! next marker
         stop_markers(count_dash) = i ! save where the separation is
      endif
   enddo
end subroutine parse_add_param_line


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


END PROGRAM Atomic_data_analysis
