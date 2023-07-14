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
module Little_subroutines
use Universal_constants
use Objects
use Dealing_with_files, only : Count_columns_in_file, Count_lines_in_file
implicit none
PRIVATE

! this interface finds by itself which of the two subroutine to use depending on the parameters passed:
interface extend_array_size ! extend array size
    module procedure extend_array_size_int
    module procedure extend_array_size_real
end interface extend_array_size

! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array ! search cheaking one by one
    module procedure Find_in_1D_array
    module procedure Find_in_2D_array
end interface Find_in_array
! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array_monoton ! search with bisection method
    module procedure Find_in_monotonous_1D_array
    module procedure Find_in_monotonous_2D_array
end interface Find_in_array_monoton

! this is a function that returns the order of a passed number:
interface find_order_of_number
   module procedure find_order_of_number_real
   module procedure find_order_of_number_int
end interface find_order_of_number   

! this is a function that deallocates array after checking if it is allocated:
interface deallocate_array
   module procedure deallocate_array_int1d
   module procedure deallocate_array_int2d
   module procedure deallocate_array_int3d
   module procedure deallocate_array_r1
   module procedure deallocate_array_r2
   module procedure deallocate_array_r3
   module procedure deallocate_array_ch1
   module procedure deallocate_array_c1
   module procedure deallocate_array_c2
end interface deallocate_array   


! private :: ! hides items not listed on public statement 
public :: Find_in_array, Find_in_array_monoton, extend_array_size, deallocate_array, Find_monotonous_LE, Fermi_interpolation, &
linear_interpolation, Find_in_monotonous_1D_array, Gaussian, print_time_step, fast_pow, count_3d, print_progress, &
interpolate_data_on_grid, number_of_types_of_orbitals, name_of_orbitals, order_of_time, set_starting_time, convolution, &
sample_gaussian, Fermi_function, d_Fermi_function, print_time, parse_yes_no, parse_time




 contains
 
 
subroutine deallocate_array_int1d(X)
   integer, dimension(:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_int1d

subroutine deallocate_array_int2d(X)
   integer, dimension(:,:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_int2d

subroutine deallocate_array_int3d(X)
   integer, dimension(:,:,:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_int3d

 subroutine deallocate_array_r1(X)
   real(8), dimension(:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_r1

subroutine deallocate_array_r2(X)
   real(8), dimension(:,:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_r2

subroutine deallocate_array_r3(X)
   real(8), dimension(:,:,:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_r3

subroutine deallocate_array_ch1(X)
   character(*), dimension(:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_ch1

 subroutine deallocate_array_c1(X)
   complex, dimension(:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_c1

subroutine deallocate_array_c2(X)
   complex, dimension(:,:), allocatable, intent(inout) :: X
   if (allocated(X)) deallocate(X)
end subroutine deallocate_array_c2

 
 
 
function number_of_types_of_orbitals(norb, cartesian) result(n_types)
   integer :: n_types
   integer, intent(in) :: norb  ! number of orbitals per atom
   logical, intent(in), optional :: cartesian
   logical :: do_cartesian
   if (present(cartesian)) then ! user defines if cartesian basis is used
      do_cartesian = cartesian
   else ! spherical basis by default
      do_cartesian = .false.
   endif

   if (.not.do_cartesian) then ! Spherical
      select case (norb)   ! identify basis set
      case (1) ! s
         n_types = 1    ! s
      case (4) ! sp3
         n_types = 2    ! s and p
      case (5) ! sp3s*
         n_types = 3    ! s, p and s*
      case (9) ! sp3d5
         n_types = 3    ! s, p and d
      case (10) ! sp3d5s*
        n_types = 4    ! s, p, d, and s*
      endselect
   else ! Cartesian
      select case (norb)   ! identify basis set
      case (1) ! s
         n_types = 1    ! s
      case (4) ! sp3
         n_types = 2    ! s and p
      case (5) ! sp3s*
         n_types = 3    ! s, p and s*
      case (10) ! sp3d6
         n_types = 3    ! s, p and d
      case (11) ! sp3d6s*
        n_types = 4    ! s, p, d, and s*
      endselect
   endif
end function number_of_types_of_orbitals


function name_of_orbitals(norb, i_orb) result(orb_name)
   character(2) :: orb_name
   integer :: norb  ! number of orbitals per atom
   integer :: i_orb ! index of orbital
   select case (norb)   ! identify basis set
   case (1) ! s
      orb_name = 's'
   case (4) ! sp3
      select case (i_orb)
      case (2)
         orb_name = 'p'
      case default
         orb_name = 's'
      endselect
   case (5) ! sp3s*
      select case (i_orb)
      case (3)
         orb_name = 's*'
      case (2)
         orb_name = 'p'
      case default
         orb_name = 's'
      endselect
   case (9) ! sp3d5
      select case (i_orb)
      case (3)
         orb_name = 'd'
      case (2)
         orb_name = 'p'
      case default
         orb_name = 's'
      endselect
   endselect
end function name_of_orbitals
 
 

! Set the right starting time, depending on whether we use the pulse or not:
subroutine set_starting_time(laser, tim, t_start, t_NA, t_Te_Ee)
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   real(8), intent(in) :: t_start   ! user-provided value for the starting time [fs]
   real(8), intent(inout) :: tim    ! defined starting time step [fs]
   real(8), intent(inout), optional :: t_Te_Ee ! time when we switch from Te=const, to Ee=const [fs] (<0 if unused)
   real(8), intent(inout), optional :: t_NA ! time when we switch on nonadiabatic terms [fs]
   tim = t_start !-50.0d0 + dble(CEILING(min(minval(laser(:)%t0-laser(:)%t*2.35d0), 0.0d0)))  ! [fs]
   if (present(t_NA)) t_NA = tim + t_NA + 1d-6           ! [fs]
   if (present(t_Te_Ee)) t_Te_Ee = tim + t_Te_Ee         ! [fs]
   tim = min(tim,t_start)    ! check what user set
end subroutine set_starting_time

 
pure function Fermi_function(rcut, d, r) result(F)
   real(8) F
   real(8), intent(in) :: rcut, d, r
   if (r >= rcut + 10.0d0*d) then
      F = 0.0d0
   else
      F  = 1.0d0/(1.0d0 + exp( (r - rcut)/d))
   endif
end function Fermi_function


pure function d_Fermi_function(rcut, d, r) result(F)
   real(8) F
   real(8), intent(in) :: rcut, d, r   ! [A]
   real(8) :: exprd, exprd1
   if (r >= rcut + 10.0d0*d) then
      F = 0.0d0
   else
      exprd = exp( (r - rcut)/d )
      exprd1 = 1.0d0 + exprd
      F = -exprd/(d*exprd1*exprd1)
   endif
end function d_Fermi_function

 

pure function find_symbol_in_string(string, symb) result(N)
   integer :: N    ! position of the symbol in the string (0 if there is no such symbol)
   character(*), intent(in) :: string   ! given string of symbols
   character(1), intent(in) :: symb ! given symbol to be found in the string
   integer :: i, ilen
   if (verify(symb,trim(adjustl(string))) == 0) then !there is such a symbol
      ilen = LEN(string)
      do i = 1, ilen
         if (symb == string(i:i)) then  ! found it
            N = i   ! save it for output
            exit    ! no need to scan more
         endif
      enddo
   else ! there is no such symbol
      N = 0
   endif
end function find_symbol_in_string


 
pure function count_3d(Nx,Ny,Nz,Cx,Cy,Cz) ! function counts linear index number of cell on 3-d mash:
   integer, intent(in) :: Nx,Ny,Nz,Cx,Cy,Cz
   integer :: count_3d
   count_3d = 1 + (2*Nx+1)*(2*Ny+1)*(Nz+Cz) + (2*Nx+1)*(Ny+Cy) + (Nx+Cx)
end function count_3d
 
 
function fast_pow(x,n)
! this function calculates integer power of a variable x for powers up to 13 in a fast manner:
   real(8) :: fast_pow
   real(8), intent(in) :: x
   integer, intent(in) :: n
   real(8) x2, x3, x4, x5, x6, x12
   select case(n)
   case(0)
      fast_pow = 1.0d0
   case (1)
      fast_pow = x
   case (2)
      fast_pow = x*x
   case (3)
      x2 = x*x
      fast_pow = x2*x
   case (4)
      x2 = x*x
      fast_pow = x2*x2
   case (5)
      x2 = x*x
      x3 = x2*x
      fast_pow = x3*x2
   case (6)
      x2 = x*x
      x3 = x2*x
      fast_pow = x3*x3
   case (7)
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      fast_pow = x4*x3
   case (8)
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      fast_pow = x4*x4
   case (9)
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      x5 = x2*x3
      fast_pow = x5*x4
   case (10)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      x5 = x2*x3
      fast_pow = x5*x5
   case (11)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      x5 = x2*x3
      x6 = x3*x3
      fast_pow = x6*x5
   case (12)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      !x5 = x2*x3
      x6 = x3*x3
      fast_pow = x6*x6
   case (13)
      x2 = x*x
      x3 = x2*x
      !x4 = x2*x2
      !x5 = x2*x3
      x6 = x3*x3
      fast_pow = x6*x6*x
   case (24)
      x2 = x*x
      x3 = x2*x
      x6 = x3*x3
      x12 = x6*x6
      fast_pow = x12*x12
   case default
      print*, 'This power is not supported in function "fast_pow"'
      pause 'What now?'
   end select
end function fast_pow


subroutine convolution(FN, Gaus_conv)
   integer, intent(in) :: FN ! file number
   real(8), intent(in) :: Gaus_conv      ! [fs] gaussian FWHM of the probe pulse
   !=====================================
   integer i, j, FN2, N, M, Reason, i_start
   real(8), dimension(:,:), allocatable :: Spectr
   real(8), dimension(:,:), allocatable :: Conv_Spectr
   real(8) temp, temp2(3), Gaus, t_g, dt
   character(500) :: File_name, File_name2, First_line, Second_line
   !character(100), dimension(:), allocatable :: First_line
   logical file_exists, file_opened, file_named

   t_g = Gaus_conv ! [fs] FWHM of gaussian to convolve with

   !open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(UNIT=FN, opened=file_opened, named=file_named, name=File_name)
   
   exists:if (file_opened .and. file_named) then
      ! Input file:
      rewind(FN) ! start reading file from the first line
      call Count_columns_in_file(FN, M, 2)   ! module "Dealing_with_files"
      call Count_lines_in_file(FN, N)  ! module "Dealing_with_files"
      N = N - 2
      allocate(Spectr(N,M))
      allocate(Conv_Spectr(N,M))

      Spectr = 0.0d0
      Conv_Spectr = 0.0d0
      read(FN,'(a)',IOSTAT=Reason) First_line
      read(FN,'(a)',IOSTAT=Reason) Second_line
      if (Second_line(1:1) /= '#') then ! it is a functional line, not a comment
         backspace(FN) ! get back and read this line then
      endif

      do i = 1, N
         read(FN,*,IOSTAT=Reason) Spectr(i,:)
         IF (Reason .GT. 0)  THEN ! ... something wrong ...
            write(*,'(a,a,a,i5,a)') 'Problem reading file ', trim(adjustl(File_name)), ' in line ', i, ', wrong type of variable'
            goto 1010
         ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
            write(*,'(a,a,a,i5,a)') 'Problem reading file ' ,trim(adjustl(File_name)), ' in line ', i, ', unexpected END of file'
            goto 1010
         ELSE   ! normal reading
         END IF
      enddo
      Conv_Spectr(:,1) = Spectr(:,1)
      dt = Spectr(2,1) - Spectr(1,1) ! [fs], timestep, assuming equal steps everywhere

      ! Analysis:
      File_name2 = trim(adjustl(File_name(1:len(trim(adjustl(File_name)))-4)))//'_CONVOLVED.dat'

      open (newunit=FN2, file=trim(adjustl(File_name2))) ! output file
      write(FN2,'(a)')  First_line  !e.g. '#Time[fs]	Reflectivity	Transmittion	Absorbtion  n   k   Re(e)   Im(e)'
      i_start = min(int(20/dt),N)
      do i = 1,N ! all time steps
         if (Spectr(i,1) <= Spectr(1,1)+t_g*3.0+dt*0.5d0) then ! incomplete data at the start
            do j = 1,min(CEILING((Spectr(1,1)+t_g*3.0-Spectr(i,1))/dt),N)  ! convolve with initial value
               call Gaussian(mu=Spectr(i,1), sigma=t_g, x=(Spectr(i,1)-t_g*3.0)+(Spectr(j,1)-Spectr(1,1)), Gaus=Gaus)  ! Gaussian weight
               Conv_Spectr(i,2:) = Conv_Spectr(i,2:) + Spectr(i_start,2:)*Gaus   ! Convolution of all arrays
            enddo
         endif

         do j = 1,N  ! convolve for complete data
            ! do convolution only within 3*sigma interval:
            if ((Spectr(j,1) >= Spectr(i,1)-t_g*3.0-dt*0.5d0) .and. (Spectr(j,1) <= Spectr(i,1)+t_g*3.0+dt*0.5d0)) then
               call Gaussian(mu=Spectr(i,1), sigma=t_g, x=Spectr(j,1), Gaus=Gaus)  ! Gaussian weight
               Conv_Spectr(i,2:) = Conv_Spectr(i,2:) + Spectr(j,2:)*Gaus   ! Convolution of all arrays
            endif
         enddo
         if (Spectr(i,1) >= Spectr(N,1)-t_g*3.0+dt*0.5d0) then ! incomplete data at the end
            do j = 1,min(N,CEILING((Spectr(i,1)-(Spectr(N,1)-t_g*3.0))/dt))  ! convolve with final value
               call Gaussian(mu=Spectr(i,1), sigma=t_g, x=Spectr(N,1)+(Spectr(j,1)-Spectr(1,1))+dt, Gaus=Gaus)  ! Gaussian weight
               Conv_Spectr(i,2:) = Conv_Spectr(i,2:) + Spectr(N-1,2:)*Gaus   ! Convolution of all arrays
            enddo
         endif
         WRITE(FN2,'(e25.16,$)') Conv_Spectr(i,:)
         write(FN2,'(a)') ' '
      enddo
      close(FN2) ! close the new file
1010 continue
   else exists
      print*, 'In subroutine convolution'
      print*, 'file ', trim(adjustl(File_name)), ' not found.'
   endif exists
! print*, 'convolution is done!'
end subroutine convolution



subroutine order_of_time(tim, text, gnu_text, x_tics)
   real(8), intent(in) :: tim ! time to find its order
   character(*), intent(out) :: text ! fs, ps, ns, mks, ms, s
   character(*), intent(out), optional :: gnu_text ! culomn to set in gnuplot
   real(8), intent(out), optional :: x_tics ! tics for gnuplot
   integer :: time_ord
   time_ord = find_order_of_number(tim) ! module "Little_subroutines"
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


subroutine parse_time(chtest, sec_in, c0_in, c1_in)
   character(*), intent(out) :: chtest ! split it into miuns, hours, days...
   real(8), intent(inout), optional :: sec_in   ! time interval in [sec]
   integer, dimension(8), intent(inout), optional :: c1_in, c0_in ! time stamp
   !-------------------------
   character(100) temp
   integer, dimension(8) :: c1 ! time stamp
   real(8) days, hours, mins, sec

   if (present(sec_in)) then  ! data provided in total number of seconds
      sec = sec_in
   elseif (present(c1_in) .and. present(c0_in)) then   ! data provided in fortran time-stamp format
      sec = get_seconds_from_timestamp(c0_in, c1_in)   ! below
   elseif (present(c0_in)) then
      call date_and_time(values=c1) ! current time
      sec = get_seconds_from_timestamp(c0_in, c1)   ! below
   else
      sec = 0.0d0 ! no data provided, nothing to printout
   endif

   days = 0.0d0
   hours = 0.0d0
   mins = 0.0d0
   if (sec .GE. 60.0d0) then
      mins = FLOOR(sec/60.0d0)
      sec = sec - mins*60.0d0
      if (mins .GT. 60.0d0) then
         hours = FLOOR(mins/60.0d0)
         mins = mins - hours*60.0d0
         if (hours .GT. 24.0d0) then
            days = FLOOR(hours/24.0d0)
            hours = hours - days*24.0d0
         endif
      endif
   endif
   chtest = ''
   temp = ''
   if (days .GT. 1.0d0) then
      write(temp, '(i9)') int(days)
      write(chtest, '(a,a)') trim(adjustl(temp)), ' days'
   else if (days .GT. 0.0d0) then
      write(temp, '(i9)') int(days)
      write(chtest, '(a,a)') trim(adjustl(temp)), ' day'
   endif
   if (hours .GT. 1.0d0) then
      write(temp, '(i9)') int(hours)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' hours'
   else if (hours .GT. 0.0d0) then
      write(temp, '(i9)') int(hours)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' hour'
   endif
   if (mins .GT. 1.0d0) then
      write(temp, '(i9)') int(mins)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' mins'
   else if (mins .GT. 0.0d0) then
      write(temp, '(i9)') int(mins)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' min'
   endif
   write(temp, '(f7.3)') sec
   write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' sec'
end subroutine parse_time


pure function get_seconds_from_timestamp(c0, c1) result(Sec)
   real(8) Sec
   integer, dimension(8), intent(in) :: c1, c0 ! timestamp, current and starting
   Sec = dble(24.0d0*3600.0d0*(c1(3)-c0(3)) + &    ! days
               3600.0d0*(c1(5)-c0(5)) + &          ! hours
               60.0d0*(c1(6)-c0(6)) + &            ! minutes
               (c1(7)-c0(7)) + &                   ! seconds
               (c1(8)-c0(8))*0.001d0)              ! milliseconds
end function get_seconds_from_timestamp



subroutine print_progress(string,ndone,ntotal)
    implicit none
    character*(*) string
    character(255) prog,oldprog
    integer ndone,ntotal,i
    if (100.0*ndone/ntotal .GE. 100.0d0) then
        write(0,'(a,$)') '                                                                                   ',char(13)
    else
        write(prog,'(a25,1x,''['')') string
        do i=1,40
            prog(27+i:27+i)=' '
        enddo
        write(prog(44:51),'(f7.2,''%'')') 100.0*ndone/ntotal
        do i=1,40
            if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
                if (prog(27+i:27+i).eq.' ') prog(27+i:27+i)='-'
            endif
        enddo
        prog(67:67)=']'
        write(0,'(a,a,$)') prog(1:77),char(13)
        return
    endif
end subroutine print_progress


pure subroutine Gaussian(mu, sigma, x, Gaus, normalized_max) ! at the time x according to Gaussian shape
   real(8), intent(in) :: mu, sigma, x
   real(8), intent(out) :: Gaus
   real(8), parameter :: Pi = 3.1415926535897932384626433832795d0
   real(8), intent(in), optional :: normalized_max
   if (present(normalized_max)) then
      Gaus = exp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))*normalized_max ! it will be equal to "normalized_max" at the maximum
   else ! normalized as integral(Gaus) = 1
      Gaus = 1.0d0/(sqrt(2.0d0*Pi)*sigma)*exp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma)) ! it will be normalized to integral=1
   endif
end subroutine Gaussian


function sample_gaussian(x0, sigma, if_FWHM) result(x_out)
   real(8) x_out
   real(8), intent(in) :: x0, sigma ! gaussian mean and sigma (by default, it's not FWHM)
   logical, intent(in), optional :: if_FWHM   ! if given sigma is FWHM (instead of standard deviation)
   !---------------------
   real(8), parameter :: Pi = 3.1415926535897932384626433832795d0
   real(8) :: rvt1, rvt2, t0, s0

   s0 = sigma
   if (present(if_FWHM)) then ! check if it is FWHM instead of sigma
      if (if_FWHM) s0 = sigma/( 2.0d0*sqrt( 2.0d0*log(2.0d0) ) )
   endif

   t0 = -10.0d0 * s0 ! to start with
   do while (abs(t0) > 3.0e0*s0) ! accept only within 3*sigma interval
      call random_number(rvt1)   ! Random value 1 for generated time t0
      call random_number(rvt2)   ! Random value 2 for generated time t0
      t0 = s0 * COS(2.0d0*Pi*rvt2)*SQRT(-2.0d0*LOG(rvt1))    ! Gaussian distribution for t0 around 0
      x_out = x0 + t0   ! around the given mean value x0
   enddo
end function sample_gaussian



subroutine interpolate_data_on_grid(given_array_x, given_array_y, given_grid, array_to_fill)
   real(8), dimension(:), intent(in) :: given_array_x, given_array_y, given_grid
   real(8), dimension(:), intent(out) :: array_to_fill
   !----------------------------------------------
   integer :: i, N, i_closest
   real(8) :: eps
   eps = 1.0d-3
   N = min(size(given_grid), size(array_to_fill))
   do i = 1, N
      call Find_monotonous_LE(given_array_x,  given_grid(i), i_closest)	! see below
      if (given_grid(i) >= given_array_x(i_closest)) then
         if (i_closest > 1) then
            do while (ABS(given_array_x(i_closest) - given_array_x(i_closest-1)) < eps) ! exclude the same value
               i_closest = i_closest + 1
               if (i_closest == N) exit
            enddo
         endif
         call linear_interpolation(given_array_x, given_array_y, given_grid(i), array_to_fill(i), i_closest)	! see below
      else ! no photoabsorption below the edge:
         array_to_fill(i) = 0.0d0
      endif
   enddo
end subroutine interpolate_data_on_grid



pure subroutine linear_interpolation(xarray, yarray, x, y, i, x0, y0, replac)
   real(8), dimension(:), intent(in) :: xarray, yarray	! x-array, y-array
   real(8), intent(in) :: x   ! input
   real(8), intent(out) :: y  ! output
   integer, intent(in) :: i   ! index for x-array
   real(8), intent(in), optional :: x0, y0 ! assume initial value for i = 1
   logical, intent(in), optional :: replac ! replace x0 by the given one independantly on which i is it?
   REDO: if (.not.present(replac)) then
    if (i .GT. 1) then
      if (x - xarray(i-1) .GE. 0.0d0) then
         y = yarray(i-1) + (yarray(i) - yarray(i-1))/(xarray(i) - xarray(i-1))*(x - xarray(i-1))
      else
         if (present(y0) .and. present(x0)) then
            y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
         else
            if (present(x0)) then
               y = (yarray(i) - 0)/(xarray(i) - x0)*(x - x0)
            else
               y = (yarray(i) - 0)/(xarray(i) - 0)*(x - 0)
            endif
         endif
      endif
    else
      if (present(y0) .and. present(x0)) then
         y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
      else
         if (present(x0)) then
            y = (yarray(i) - 0)/(xarray(i) - x0)*(x - x0)
         else
            y = (yarray(i) - 0)/(xarray(i) - 0)*(x - 0)
         endif
      endif
    endif
   else REDO ! in this case use just the given values:
    if ((replac) .ANd. present(x0) .AND. present(y0)) then
       y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
    else
      if (i .GT. 1) then
         if (x - xarray(i-1) .GE. 0.0d0) then
            y = yarray(i-1) + (yarray(i) - yarray(i-1))/(xarray(i) - xarray(i-1))*(x - xarray(i-1))
         else
            if (present(y0) .and. present(x0)) then
               y = y0 + (yarray(i) - y0)/(xarray(i) - x0)*(x - x0)
            else
               if (present(x0)) then
                  y = (yarray(i) - 0)/(xarray(i) - x0)*(x - x0)
               else
                  y = (yarray(i) - 0)/(xarray(i) - 0)*(x - 0)
               endif
            endif
         endif
      endif
    endif
   endif REDO
end subroutine linear_interpolation



subroutine Fermi_interpolation(xarray, yarray, x, y, i)
   real(8), dimension(:), intent(in) :: xarray, yarray	! x-array, y-array
   real(8), intent(in) :: x	! input
   real(8), intent(out) :: y	! output
   integer, intent(in) :: i	! index for x-array
   real(8) :: mu, T ! parameters of the fermi function to be found and used
   real(8) :: temp
   if (i > 1) then
      if ((yarray(i) == yarray(i-1)) .or. (xarray(i) == xarray(i-1))) then ! no need for doing anything
         y = yarray(i)
      else
         if ((yarray(i-1) < tiny(x)) .or. (yarray(i) < tiny(x)) .or. (yarray(i-1) == 2.0d0) .or. (yarray(i) == 2.0d0)) then ! just in case it's 0
            call linear_interpolation(xarray, yarray, x, y, i)
         else
            temp = log(2.0d0/yarray(i)-1.0d0)
            T = (temp - log(2.0d0/yarray(i-1)-1.0d0))/(xarray(i)-xarray(i-1))
            mu = xarray(i) - T*temp
            y = 2.0d0/(1.0d0 + exp((x - mu)/T))
         endif
      endif
   else ! use default value
      y = yarray(1)
   endif   
end subroutine Fermi_interpolation



subroutine Find_in_1D_array(Array, Value, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(i) .LT. Value-1.0d-10)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_1D_array

subroutine Find_in_2D_array(Array, Value, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(Indx,i) .LT. Value)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_2D_array




subroutine Find_monotonous_LE(Array, Value0, Number)    ! monotoneausly increasing array
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2

   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f25.16,f25.16,f25.16,f25.16)') Value0, Array(i_cur), Array(i_1), Array(i_2)
        pause 'STOPPED WORKING...'
   else
       N = size(Array)
       if (Value0 .LE. Array(1)) then ! it's the first value, no need to search
           i_cur = 1
       else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
           i_cur = N
       else
           i_1 = 1
           i_2 = N
           i_cur = FLOOR((i_1+i_2)/2.0)
           temp_val = Array(i_cur)
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LT. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LT. Value0) then
                   i_1 = i_cur
                else
                   i_2 = i_cur
                endif
                i_cur = FLOOR((i_1+i_2)/2.0)
                temp_val = Array(i_cur)
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
   Number = i_cur
end subroutine Find_monotonous_LE





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



subroutine extend_array_size_int(array1)
   integer, dimension(:), allocatable, intent(inout) :: array1
   integer N
   integer, dimension(:), allocatable :: array2
   N = size(array1)
   allocate(array2(N))
   array2 = array1
   deallocate(array1)
   allocate(array1(2*N))
   array1(1:N) = array2(1:N)
   deallocate(array2)
end subroutine extend_array_size_int

subroutine extend_array_size_real(array1)
   real(8), dimension(:), allocatable, intent(inout) :: array1
   integer N
   real(8), dimension(:), allocatable :: array2
   N = size(array1)
   allocate(array2(N))
   array2 = array1
   deallocate(array1)
   allocate(array1(2*N))
   array1(1:N) = array2(1:N)
   deallocate(array2)
end subroutine extend_array_size_real



subroutine print_time_step(text, num, msec)
   CHARACTER(len=*) :: text	! text to print out
   real(8), intent(in), optional :: num	! to print out this number
   logical, intent(in), optional :: msec ! print msec or not?
   character(len=100) :: var 
   integer c1(8) ! time stamps
   call date_and_time(values=c1) ! standard FORTRAN time and date
   if (present(num) .and. present(msec)) then
      write(var,'(f8.2)') num
      write(*, 1002) trim(adjustl(text))//' '//trim(adjustl(var))//' at ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
   elseif (present(msec)) then
      write(*, 1002) trim(adjustl(text))//' at ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
   elseif (present(num)) then
      write(var,'(f8.2)') num
      write(*, 1001) trim(adjustl(text))//' '//trim(adjustl(var))//' at ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
   else
      write(*, 1001) trim(adjustl(text))//' at ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
   endif
   ! Formats for printing:
   1001 format (a, i2.2, ':', i2.2, ':', i2.2, ' on ', i2.2, '/', i2.2, '/', i4.4)
   1002 format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, ' on ', i2.2, '/', i2.2, '/', i4.4)
end subroutine print_time_step


subroutine trim_zeros(in_line)
   character(len=*), intent(inout) :: in_line
   integer :: i, j, leng
   i = 0
   leng = LEN(trim(in_line))
   j = leng
   do while (in_line(j:j) == '0') ! trim extra zeros from the end of the line
      in_line(j:j) = ' ' ! delete last characters
      i = i + 1
      j = leng-i
   enddo
   if (in_line(j:j) == '.') in_line(j+1:j+1) = '00'
end subroutine trim_zeros


subroutine print_time(text, ind, iter, ctim)
   CHARACTER(len=*) :: text	! text to print out
   integer, intent(in), optional :: ind, iter	! ind = to print out miliseconds or not; iter = integer number to print out
   integer, intent(in), optional ::  ctim(8)	! given time to print it out
   integer c1(8) ! time stamps
   character(len=100) :: var 
   call date_and_time(values=c1) ! standard FORTRAN time and date
   if (present(ind)) then
      if (present(ctim)) then
         if (present(iter)) then
            write(var,'(i12)') iter
            write(*, 1002) trim(adjustl(text))//' '//trim(adjustl(var))//' ', ctim(5), ctim(6), ctim(7), ctim(8), ctim(3), ctim(2), ctim(1)
         else
            write(*, 1002) trim(adjustl(text))//' ', ctim(5), ctim(6), ctim(7), ctim(8), ctim(3), ctim(2), ctim(1)
         endif
      else
         if (present(iter)) then
            write(var,'(i12)') iter
            write(*, 1002) trim(adjustl(text))//' '//trim(adjustl(var))//' ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
         else
            write(*, 1002) trim(adjustl(text))//' ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
         endif
      endif
   else
      if (present(ctim)) then
         if (present(iter)) then
            write(var,'(i12)') iter
            write(*, 1001) trim(adjustl(text))//' '//trim(adjustl(var))//' ', ctim(5), ctim(6), ctim(7), ctim(8), ctim(3), ctim(2), ctim(1)
         else
            write(*, 1001) trim(adjustl(text))//' ', ctim(5), ctim(6), ctim(7), ctim(3), ctim(2), ctim(1)
         endif
      else
         if (present(iter)) then
            write(var,'(i12)') iter
            write(*, 1001) trim(adjustl(text))//' '//trim(adjustl(var))//' ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
         else
            write(*, 1001) trim(adjustl(text))//' ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
         endif
      endif
   endif
   ! Formats for printing:
   1001 format (a, i2.2, ':', i2.2, ':', i2.2, ' on ', i2.2, '/', i2.2, '/', i4.4)
   1002 format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, ' on ', i2.2, '/', i2.2, '/', i4.4)
end subroutine print_time


subroutine parse_yes_no(string, yesno)
   character(*), intent(in) :: string ! figure out whether it is 'yes' or 'no'
   logical, intent(out) :: yesno ! yes = true, no = false
   select case (trim(adjustl(string))) 
   case ('y', 'Y', 'yes', 'YES', 'Yeah', 'Yes', 'yEs', 'yeS', 'YEs', 'yES', 'YeS', '1', 'true', 'Da', 'Ja', 'da', 'aha', 'ja', 'Aha')
      yesno = .true.
   case default
      yesno = .false.
   end select
end subroutine parse_yes_no


end module Little_subroutines
