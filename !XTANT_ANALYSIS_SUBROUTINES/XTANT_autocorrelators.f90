PROGRAM Autocorrelators
! Compilation:
!
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs XTANT_autocorrelators.f90 -o XTANT_autocorrelators.exe /link /stack:9999999999 
!
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays XTANT_autocorrelators.f90 -o XTANT_autocorrelators.exe /link /stack:9999999999 
!
! To execute:
! XTANT_autocorrelators.exe  alpha  time_step
! "alpha" is the exponential factor suppressing correlation at long times exp(-alph a*t)
! "time_step" sets how often to print our autocorrelators and phonon spectra
! Default values are: alpha=1e3 fs, dt=200 fs
!<===========================================

use Universal_constants

type Instant_data
   real(8) :: Tim   ! time instant [fs]
   real(8), dimension(:,:), allocatable :: V    ! (Vx,Vy,Vz) velosity vectors for all atoms
end type Instant_data

character(200) :: File_temperatures, File_velosites
character(200) :: File_name1, File_name2, File_name_out, File_name_out1
character(32), dimension(10) :: char_var
character(10) :: temp_ch
character(1) :: path_sep

type(Instant_data), dimension(:), allocatable :: Step           ! All data and parameters at this timestep
real(8), dimension(:), allocatable :: psi   ! VAF
real(8) :: alpha, time_dt, cur_t, dt, df
! complex(8), dimension(:), allocatable :: Vibration_spectrum
real(8), dimension(:), allocatable :: Vibration_spectrum
integer :: FN1, FN2, FN_out, FN_out1	! file number
integer :: Reason, i, j, Tsiz, time_print, Nfreq
logical :: read_well

call Path_separator(path_sep)  ! Objects_and_types

! Set defaults:
FN_out = 1000
FN_out1 = 1001
FN1 = 9999
FN2 = 9998
File_temperatures = 'OUTPUT_temperatures.dat'	! default name, for start
File_velosites = 'OUTPUT_coordinates_and_velocities.dat'    ! default name

File_name_out1 = 'OUT_VAF'
File_name_out =  'OUT_vibrational_spectrum'

! Default values:
alpha = 1.0d20	   ! start from infinity, i.e. no cut off
time_dt = 200.0  ! calculate autocorrelators every 200 fs
time_print = 0  ! to start with 

!---------------------------------------
print*, '**********************************************************************************************'
print*, 'Set two numbers: alpha and number of time-points'
print*, 'for analysis of autocorrelation functions'
print*, 'Alpha is the exponential factor suppressing correlation at long times exp(-alpha*t)'
print*, 'Time-points sets how often to print out autocorrelators and phonon spectra'
print*, 'Default values are: alpha=1e3 fs, TP sets number to printout every 200 fs'
print*, '**********************************************************************************************'


! Get alpha and time prints, if user defined it:
! print*, 'Arguments passed:', iargc()
do i = 1, iargc()
   call getarg(i, char_var(i), status=Reason)	! read only the first one
   if (i == 1) then
       read(char_var(i),*) alpha   ! alpha coefficient
   elseif (i == 2) then
       read(char_var(i),*) time_print  ! number of grid points, at which to calculate autocorrelators
   else
      write(*,'(a)') char_var(i)
   endif
end do

alpha = 1.0d0/(alpha*alpha) ! convert into 1/tau

! Get the time grid from the file:
open (unit=FN1, file=trim(adjustl(File_temperatures)), status = 'old', readonly)
call read_time_grid(FN1, Step, time_print, time_dt)    ! below
close(FN1)
Tsiz = size(Step)   ! time grid size
! print*, 'Tsiz', Tsiz

! Get the velosities from the file:
open (unit=FN2, file=trim(adjustl(File_velosites)), status = 'old', readonly)
call read_velosities(FN2, Step)    ! below
close(FN2)

dt = Step(2)%Tim - Step(1)%Tim
Nfreq = 5000     ! points in frequency space
df = 1.0d-4 ! omega [1/fs]

! Do calculations for the chosen times:
 cur_t = Step(1)%Tim    ! start from here
do i = 1, Tsiz
   if (Step(i)%Tim > cur_t) then    ! calcualte
      cur_t = cur_t + time_dt   ! next step will be at this time
      print*, 'Step:', Step(i)%Tim
      
      ! Get velosity autocorrelation function (VAF) :
      call get_VAF(Step, psi, alpha, i) ! below
      
      ! Get the vibrational spectrum from VAF:
      !call make_FFT(psi, Vibration_spectrum) ! Fast Fourier Transform
      
      if (.not.allocated(Vibration_spectrum)) allocate(Vibration_spectrum(Nfreq))	! the same dimension   
      call make_Fourier(Step, psi, Vibration_spectrum, df ,dt, i) ! below
      
      ! Printout the results:
!       print*, i, size(Step)
      write(temp_ch,'(i6)') FLOOR(Step(i)%Tim)
      open (unit=FN_out1, file=trim(adjustl(File_name_out1))//'_'//trim(adjustl(temp_ch))//'.dat')
      open (unit=FN_out, file=trim(adjustl(File_name_out))//'_'//trim(adjustl(temp_ch))//'.dat')
      do j =1,size(Step)
         write(FN_out1,'(es,es)') Step(j)%Tim, psi(j)
      enddo
      do j =1,size(Vibration_spectrum)
         write(FN_out,'(es,es)') (df*dble(j)*1.0d3/(2.0d0*g_Pi)), dble(Vibration_spectrum(j))    ! [THz]
      enddo
      close (FN_out)
      close (FN_out1)
      if (allocated(psi)) deallocate(psi)
      if (allocated(Vibration_spectrum)) deallocate(Vibration_spectrum)
      
   endif
enddo


2012 continue   ! to exit the program
STOP
!---------------------
 contains
 


subroutine make_Fourier(Step, Vector, FFT_vector, df ,dt, i_t)
   type(Instant_data), dimension(:), intent(in) :: Step           ! All data and parameters at this timestep
   real(8), dimension(:), intent(in) :: Vector	! input
   real(8), dimension(:), intent(out) :: FFT_vector	! output: discrete Fourier-transformed Vector
   real(8), intent(in) :: df, dt
   integer, intent(in) :: i_t    ! starting time index
   real(8) :: f, tmp, time, term
   integer :: i, r, Nsiz
   Nsiz = size(Vector)
   
   f = 0.0d0
   do r = 1, size(FFT_vector)
      f = f + df
      term = 0.0d0
      ! transform into frequency space
      do i = 1, Nsiz    ! time steps
!          time = (i-1) * dt
         time = Step(i)%Tim - Step(i_t)%Tim
         tmp  = f * time  ! frequency in [1/fs] and time in [fs]
         term = term + Vector(i) * cos(tmp)
      enddo
      FFT_vector(r) = term * dt
   enddo
!    FFT_vector(:) = FFT_vector(:) !* 1.0d5    ! units conversion
end subroutine make_Fourier

 
 
subroutine get_VAF(Step, psi, alpha, i_t)
   type(Instant_data), dimension(:), intent(in) :: Step           ! All data and parameters at this timestep
   real(8), dimension(:), allocatable :: psi    ! VAF
   integer, intent(in) :: i_t    ! starting time index
   real(8), intent(in) :: alpha ! coefficient
   integer :: i, Nsiz, j, Nat
   real(8) :: term, time
   Nsiz = size(Step)
   if (allocated(psi)) deallocate(psi)
   allocate(psi(Nsiz), source=0.0d0)
   Nat = size(Step(1)%V,1)
   
   do j = i_t, Nsiz   ! all time steps
      term = 0.0d0
      ! sum_i^natom v(tau) . v(tau+t)>
      do i = 1, Nat ! all atoms
         term = term + SUM( Step(j)%V(i,:)*Step(i_t)%V(i,:) )
      enddo ! j = 1, natom
      time = Step(j)%Tim - Step(i_t)%Tim
      psi(j) = term / dble(Nat) * exp(-alpha*time*time)
!       print*, Step(j)%Tim, term,  exp(-alpha*time*time)
   enddo
!    pause 'get_VAF'
end subroutine get_VAF


 
subroutine read_velosities(FN, Step)
   integer, intent(in) :: FN    ! file number with data for time grid
   type(Instant_data), dimension(:), allocatable, intent(inout) :: Step           ! All data and parameters at this timestep
   character(500) :: read_line
   real(8), dimension(12) :: read_line_r
   integer :: i, j, Nsiz
   i = 0
   do
      i = i + 1
      read_line = ''
      read(FN,'(a)',IOSTAT=Reason) read_line
      if ((Reason /= 0) .or. (len_trim(read_line) < 1)) then
!          print*, 'Nat=', i-1
         exit
      endif
   enddo
   Nsiz = i - 1 ! number of atoms
   rewind(FN)   ! start from the start
   
   do j = 1, size(Step) ! for all time steps
      ! Allocate velosities:
      allocate(Step(j)%V(Nsiz,3))
      ! Read out velosities:
      if (j > 1) then ! skip two empty lines
         read(FN,'(a)',IOSTAT=Reason)
         read(FN,'(a)',IOSTAT=Reason)
      endif
      read_line = ''
      do i = 1, Nsiz    ! for all atoms
         read(FN,'(a)',IOSTAT=Reason) read_line
         if (Reason == 0) then
            read(read_line,*) read_line_r
            Step(j)%V(i,1:3) = read_line_r(4:6)
!             print*, 'V',  j, Step(j)%V(i,1:3)
         endif
      enddo
   enddo
end subroutine read_velosities


 

subroutine read_time_grid(FN, Step, time_print, time_dt)
   integer, intent(in) :: FN    ! file number with data for time grid
   type(Instant_data), dimension(:), allocatable, intent(inout) :: Step           ! All data and parameters at this timestep
   integer, intent(inout) :: time_print
   real(8), intent(inout) :: time_dt
   integer :: Reason, Nsiz
   real(8) eps
   eps = 1.0d-6
   call Count_lines_in_file(FN, Nsiz)  ! below
   Nsiz = Nsiz - 2  ! skip first two lines with comments
   allocate (Step(Nsiz))
   ! Read temperatures:
   read(FN,*,IOSTAT=Reason) ! skip first two lines with comments
   read(FN,*,IOSTAT=Reason)
   do i = 1, Nsiz
      read(FN,*,IOSTAT=Reason) Step(i)%Tim
      if (Reason .LT. 0) exit
   enddo
   
   if (time_print > 0) then  ! redefine printout times
      time_dt = (Step(Nsiz)%Tim - Step(1)%Tim)/dble(time_print)
   endif
   
   if (time_dt < Step(2)%Tim - Step(1)%Tim) time_dt = Step(2)%Tim - Step(1)%Tim
end subroutine read_time_grid

 


subroutine sort_array(Ev, Ar2)	! bubble method
   real(8), dimension(:), intent(inout) :: Ev	! array to sort
   real(8), dimension(:), intent(inout) :: Ar2
   real(8) :: temp_c, temp_c2
   integer N,i,j
   logical :: swapped
   N = size(Ev)
   do j = N-1, 1, -1
      do i = 1, j
         IF (Ev(i) > Ev(i+1)) THEN
            temp_c = Ev(i)
            Ev(i) = Ev(i+1)
            Ev(i+1) = temp_c
            ! And the associated array too:
            temp_c2 = Ar2(i)
            Ar2(i) = Ar2(i+1)
            Ar2(i+1) = temp_c2
            swapped = .TRUE.
         END IF
      enddo
      IF (.NOT. swapped) EXIT
   enddo
end subroutine sort_array
 

pure function linear_interpolation(x, x1, x0, y1, y0) result ( y )
   real(8), intent(in) :: x, x1, x0, y1, y0
   real(8) :: y
   y = y0 + (y1 - y0)/(x1 - x0)*(x - x0)
end function linear_interpolation


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





subroutine make_DFT(Vector, FFT_vector) ! Discrete Fourier Transform
   real(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: FFT_vector	! output: discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, i, j
   real(8) :: two_pi_N, fac
   N = size(Vector)
!    if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant

   if (.not.allocated(FFT_vector)) allocate(FFT_vector(N))	! the same dimension
   FFT_vector = 0.0d0
   two_pi_N = 2.0d0*g_Pi/dble(N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, N-1
!     tmp  = 2.d0 * PI * (lightspeed * f * 1.0d2) * (time * 1.0d-15) ! frequency in [cm^-1] and time in [fs]
      fac = two_pi_N*dble(i)
      do j = 0, N-1
         FFT_vector(i+1) = FFT_vector(i+1) + Vector(j+1) * exp(-g_CI*fac*dble(j))
      enddo
   enddo
end subroutine make_DFT




subroutine make_FFT(Vector, FFT_vector) ! Fast Fourier Transform
   real(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: FFT_vector	! output: discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, half_N, i, j, k
   real(8), dimension(:), allocatable :: odd_vec, even_vec
   complex(8), dimension(:,:), allocatable :: exp_vec
   complex(8) :: exp_prefac,  part_1, part_2
   real(8) :: two_pi_N, fac
   N = size(Vector)
   if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant
   half_N = N/2
      
   if (.not.allocated(FFT_vector)) allocate(FFT_vector(N))	! the same dimension
   ! Sub-arrays used in FFT:
   allocate(odd_vec(half_N))
   allocate(even_vec(half_N))
   allocate(exp_vec(half_N,half_N))
   
   two_pi_N = 2.0d0*g_Pi/dble(half_N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, half_N-1
      even_vec(i+1) =  Vector(2*i+1)	! X0
      odd_vec(i+1) = Vector(2*(i+1))	! X1
!       even_vec(i+1) =  Vector(2*i+1)	! X0
!       odd_vec(i+1) = Vector(2*(i))	! X1
      fac = two_pi_N*dble(i)
      FORALL(j=0:half_N-1)  exp_vec(i+1,j+1) = exp(-g_CI*fac*dble(j))
   enddo
   
   ! Construct re FFT:
   FFT_vector = (0.0d0,0.0d0)
   do i = 0, half_N-1
      part_1 = (0.0d0,0.0d0)
      part_2 = (0.0d0,0.0d0)
      do j = 0, half_N-1
          !j = k - (half_N-1)
          part_1 = part_1 + even_vec(j+1)*exp_vec(i+1,j+1)
          part_2 = part_2 + odd_vec(j+1)*exp_vec(i+1,j+1)
      enddo
      fac = 2.0d0*g_Pi/dble(N)*dble(i)
      exp_prefac =  exp(-g_CI*fac)
      ! Since Fourier transform is periodic function in time, we can shift all the outcome by half of a period:
      FFT_vector(i+1) = part_1 + exp_prefac*part_2
      FFT_vector(i+1+half_N) = part_1 - exp_prefac*part_2
   enddo
end subroutine make_FFT


! This function determines whether given integer N is odd (returns .false.) or even (returns .true.)
pure function is_even(N) result(evn)
   logical :: evn
   integer, intent(in) :: N
   if (MOD(N,2) == 0) then
      evn = .true.
   else
      evn = .false.
   endif
end function is_even



subroutine inverse_FFT (Vector, IFFT_vector) ! Fast Fourier Transform
   complex(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: IFFT_vector	! output: inverse discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, half_N, i, j, k
   complex(8), dimension(:), allocatable :: odd_vec, even_vec
   complex(8), dimension(:,:), allocatable :: exp_vec
   complex(8) :: exp_prefac,  part_1, part_2
   real(8) :: two_pi_N, fac
   N = size(Vector)
   if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant
   half_N = N/2
      
   if (.not.allocated(IFFT_vector)) allocate(IFFT_vector(N))	! the same dimension
   ! Sub-arrays used in FFT:
   allocate(odd_vec(half_N))
   allocate(even_vec(half_N))
   allocate(exp_vec(half_N,half_N))
   
   two_pi_N = 2.0d0*g_Pi/dble(half_N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, half_N-1
      even_vec(i+1) =  Vector(2*i+1)	! X0
      odd_vec(i+1) = Vector(2*(i+1))	! X1
      fac = two_pi_N*dble(i)
      FORALL(j=0:half_N-1)  exp_vec(i+1,j+1) = exp(g_CI*fac*dble(j))
   enddo
   
   ! Construct re FFT:
   IFFT_vector = (0.0d0,0.0d0)
   do i = 0, half_N-1
      part_1 = (0.0d0,0.0d0)
      part_2 = (0.0d0,0.0d0)
      do j = 0, half_N-1
          !j = k - (half_N-1)
          part_1 = part_1 + even_vec(j+1)*exp_vec(i+1,j+1)
          part_2 = part_2 + odd_vec(j+1)*exp_vec(i+1,j+1)
      enddo
      fac = 2.0d0*g_Pi/dble(N)*dble(i)
      exp_prefac =  exp(g_CI*fac)
      ! Since Fourier transform is periodic function in time, we can shift all the outcome by half of a period:
      IFFT_vector(i+1) = part_1 + exp_prefac*part_2
      IFFT_vector(i+1+half_N) = part_1 - exp_prefac*part_2
   enddo
   IFFT_vector = IFFT_vector/dble(N)
end subroutine inverse_FFT



END PROGRAM Autocorrelators
