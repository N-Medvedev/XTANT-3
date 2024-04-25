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
! This module includes some tools for performing vector algebra operations:
MODULE Algebra_tools
use Universal_constants

USE OMP_LIB, only : OMP_GET_THREAD_NUM

implicit none
PRIVATE

! Interface to automatically chose from the bubble array-sorting subroutines
interface sort_array
   module procedure sort_array_r ! for real arrays
   module procedure sort_array_c ! for complex arrays
end interface sort_array

! this interface finds by itself which of the two subroutine to use depending on the array passed:
interface mkl_matrix_mult
   module procedure mkl_matrix_mult_r	! mulitply two real matrices
   module procedure mkl_matrix_mult_c	! mulitply two complex matrices
end interface mkl_matrix_mult

! this interface finds by itself which of the two subroutine to use depending on the array passed:
interface sym_diagonalize
    module procedure r_diagonalize ! diagonilize symmetric real matrix
    module procedure c_diagonalize ! diagonilize symmetric complex matrix
end interface sym_diagonalize

! this interface finds by itself which of the subroutine to use depending on the array passed:
interface nonsym_diagonalize
    module procedure nonsym_diagonalize_r ! diagonilize non-symmetric general real matrix
    module procedure nonsym_diagonalize_c ! diagonilize non-symmetric general complex matrix
    module procedure nonsym_diagonalize_cc ! diagonilize non-symmetric general complex matrix
end interface nonsym_diagonalize

! this interface finds by itself which of the subroutine to use depending on the array passed:
interface check_Ha
    module procedure check_Ha_r
    module procedure check_Ha_c
    module procedure check_Ha_cc
end interface check_Ha

interface Two_Matr_mult
   module procedure Two_Matr_mult_r
   module procedure Two_Matr_mult_c
end interface Two_Matr_mult


real(8), parameter :: m_34 = 3.0d0/4.0d0
real(8), parameter :: m_16 = 1.0d0/6.0d0
real(8), parameter :: m_one_over_sqrtPi = 1.0d0/sqrt(g_Pi)


!private  ! hides items not listed on public statement
public :: sym_diagonalize, nonsym_diagonalize, check_Ha, Kronecker_delta, sort_array, Invers_3x3, Transpose_M
public :: double_factorial, Heavyside_tau, get_factorial, Two_Vect_Matr, Det_3x3, Cross_Prod, Matrix_Vec_Prod
public :: mkl_matrix_mult, Reciproc, check_hermiticity, Laguerre_up_to_6, d_Laguerre_up_to_6, check_symmetry
public :: d_detH_d_h_a_b, Two_Matr_mult, get_eigenvalues_from_eigenvectors, fit_parabola_to_3points
public :: make_cubic_splines, cubic_function, d_cubic_function, c8_diagonalize, mkl_matrix_mult_c8, numerical_delta
!=======================================

 contains


pure function numerical_delta(x, eta) result(delta)
   real(8), intent(in) :: x, eta
   real(8) delta  ! numerical approximation for Dirac's delta function

   if (abs(eta) < 1.0d-10) then  ! near zero width, detla -> infinity
      delta = 0.95d30 ! approximate infinity with a distinct number
   else
      delta = m_one_over_sqrtPi * exp(-(x/eta)**2)/eta
   endif
end function numerical_delta


subroutine make_cubic_splines(x_array, y_array, a, b, c, d)
   ! The subroutine for creating an array of cubic splines on arbitrary (non necesserely equidistant) grid
   ! (the grid, however, must be orderred and cannot contain degeneracies [i.e., x(i)/=x(i+1), for any i])
   ! https://en.wikipedia.org/wiki/Spline_(mathematics)
   real(8), dimension(:), intent(in) :: x_array, y_array   ! array of data to construct splines for
   real(8), dimension(:), allocatable, intent(inout) :: a, b, c, d  ! spline coefficients for the given array
   !-----------------------------
   integer :: Nsiz, i
   real(8), dimension(size(x_array)) :: h, h_inv, f_prime, f_prime2

   Nsiz = size(x_array)
   if (Nsiz > 1) then ! it makes sense to cobic-spline the array:
      ! Define the difference:
      do i = 1, Nsiz-1
         h(i) = x_array(i+1) - x_array(i)
      enddo
      h(Nsiz) = h(Nsiz - 1)
      h_inv = 1.0d0/h

      ! Define the first and second derivatives:
      do i = 1, Nsiz
         f_prime(i) = get_first_derivative(x_array, y_array, i)   ! below
         f_prime2(i) = get_second_derivative(x_array, y_array, i) ! below
      enddo

      ! Ensure that size of arrays is correct:
      if (allocated(a)) then
         if (size(a) /= Nsiz) deallocate(a)
      endif
      if (allocated(b)) then
         if (size(b) /= Nsiz) deallocate(b)
      endif
      if (allocated(c)) then
         if (size(c) /= Nsiz) deallocate(c)
      endif
      if (allocated(d)) then
         if (size(d) /= Nsiz) deallocate(d)
      endif
      allocate(a(Nsiz), source = y_array) ! a-coefficients are defined
      allocate(b(Nsiz), source = 0.0d0)
      allocate(c(Nsiz), source = 0.0d0)
      allocate(d(Nsiz), source = 0.0d0)

      ! Define b-coefficients:
      b(Nsiz) = f_prime(Nsiz)
      b(Nsiz-1) = f_prime(Nsiz-1)
      do i = Nsiz-2, 1, -1
         b(i) = f_prime(i+1) + m_34 * f_prime2(i+1)*h(i)
      enddo

      ! Define c-coefficients:
      do i = Nsiz-2, 1, -1
         c(i) = h_inv(i) * ( f_prime(i+1) - b(i) - 0.5d0 *h(i)* f_prime2(i+1) )
      enddo

      ! Define d-coefficients:
      do i = Nsiz-2, 1, -1
         d(i) = m_16*h_inv(i) * ( f_prime2(i+1) - 2.0d0*c(i) )
      enddo

   else ! does not make sence to cubic-spline it
      if (allocated(a)) deallocate(a)
      if (allocated(b)) deallocate(b)
      if (allocated(c)) deallocate(c)
      if (allocated(d)) deallocate(d)
   endif
end subroutine make_cubic_splines


pure function get_first_derivative(x_array, y_array, i) result(df)
   ! The subroutine for calculation of the finite-difference derivative on arbitrary (non necesserely equidistant) grid
   ! (the grid, however, must be orderred and cannot contain degeneracies [i.e., x(i)/=x(i+1), for any i])
   real(8) df
   real(8), dimension(:), intent(in) :: x_array, y_array
   integer, intent(in) :: i   ! index of the point, where to take derivative
   !-----------------------
   integer :: Nsiz, j

   ! Find the index for the derivative (including exceptional cases):
   if (i <= 1) then ! too small
      j = 1
   else ! (i <= 1)
      ! Define the size:
      Nsiz = size(x_array)
      if (size(y_array) < Nsiz) Nsiz = size(y_array)

      if (i >= Nsiz) then  ! too large
         j = Nsiz - 1
      else
         j = i
      endif ! (i >= Nsiz)
   endif ! (i <= 1)

   ! Get the first derivative:
   df = (y_array(j+1) - y_array(j))/(x_array(j+1)-x_array(j))
end function get_first_derivative


pure function get_second_derivative(x_array, y_array, i) result(df)
   ! The subroutine for calculation of the finite-difference second derivative on arbitrary (non necesserely equidistant) grid
   ! (the grid, however, must be orderred and cannot contain degeneracies [i.e., x(i)/=x(i+1), for any i])
   real(8) df
   real(8), dimension(:), intent(in) :: x_array, y_array
   integer, intent(in) :: i   ! index of the point, where to take derivative
   !-----------------------
   integer :: Nsiz, j
   real(8) :: xnn1, xn1n2, xnn2

   ! Find the index for the derivative (including exceptional cases):
   if (i <= 2) then ! too small
      j = 2
   else ! (i <= 2)
      ! Define the size:
      Nsiz = size(x_array)
      if (size(y_array) < Nsiz) Nsiz = size(y_array)

      if (i >= Nsiz-1) then  ! too large
         j = Nsiz - 1
      else
         j = i
      endif ! (i >= Nsiz)
   endif ! (i <= 1)

   ! Get the second derivative:
   xnn1 = x_array(j+1) - x_array(j)
   xn1n2 = x_array(j) - x_array(j-1)
   xnn2 = x_array(j+1) - x_array(j-1)
   df = 2.0d0*(y_array(j+1)*xn1n2 - y_array(j)*xnn2 + y_array(j-1)*xnn1)/(xnn2*xn1n2*xnn1)
end function get_second_derivative



pure function cubic_function(x, a, b, c, d) result(f)
   real(8) f
   real(8), intent(in) :: x, a, b, c, d
   !---------------------
   real(8) :: x2
   x2 = x*x
   f = a + b*x + c*x2 + d*x*x2
end function cubic_function


pure function d_cubic_function(x, b, c, d) result(f)
   real(8) f
   real(8), intent(in) :: x, b, c, d
   !---------------------
   real(8) :: x2
   x2 = x*x
   f = b + 2.0d0*c*x + 3.0d0*d*x2
end function d_cubic_function



pure subroutine Laguerre_up_to_6(d, L, ind_max)
   ! https://en.wikipedia.org/wiki/Laguerre_polynomials#The_first_few_polynomials
   real(8), intent(in) :: d   ! variable
   real(8), dimension(6), intent(out) :: L   ! Laguerre polinomial with indices from 0 to 5
   integer, intent(in), optional :: ind_max  ! up to which polinomial we need it?
   integer :: ind
   real(8) :: d2, d3, d4, d5
   ! set how many polynomials we need:
   if (present(ind_max)) then
      ind = ind_max
   else
      ind = 6  ! all 6 by default
   endif

   L(:) = 0.0d0   ! to start with

   ! to be reused below:
   d2 = d*d
   d3 = d2*d
   d4 = d2*d2
   d5 = d4*d
   ! Polynomials:
   L(1) = 1.0d0
   L(2) = 1.0d0 - d
   L(3) = 0.5d0*(d2 - 4.0d0*d + 2.0d0)
   L(4) = 1.0d0/6.0d0*(-d3 + 9.0d0*d2 - 18.0d0*d + 6.0d0)
   L(5) = 1.0d0/24.0d0*(d4 - 16.0d0*d3 + 72.0d0*d2 - 96.0d0*d + 24.0d0)
   L(6) = 1.0d0/120.0d0*(-d5 + 25.0d0*d4 - 200.0d0*d3 + 600.0d0*d2 - 600.0d0*d + 120.0d0)
   ! if (ind > 6) L(7) = 1.0d0/720.0d0*(d5*d - 36.0d0*d5 + 450.0d0*d4 - 2400.0d0*d3 + 5400.0d0*d2 - 4320.0d0*d + 720.0d0)
end subroutine Laguerre_up_to_6


pure subroutine d_Laguerre_up_to_6(d, L, ind_max)
   ! derivatives of the first 6 Laguerre polynomials
   real(8), intent(in) :: d   ! variable
   real(8), dimension(6), intent(out) :: L   ! derivatives of Laguerre polinomial with indices from 0 to 5
   integer, intent(in), optional :: ind_max  ! up to which polinomial we need it?
   integer :: ind
   real(8) :: d2, d3, d4, d5

   ! set how many polynomials we need:
   if (present(ind_max)) then
      ind = ind_max
   else
      ind = 6  ! all 6 by default
   endif

   L(:) = 0.0d0   ! to start with

   ! to be reused below:
   d2 = d*d
   d3 = d2*d
   d4 = d2*d2
   !d5 = d4*d
   ! Polynomials' derivatives:
   L(1) = 0.0d0
   L(2) = -1.0d0
   L(3) = 0.5d0*(2.0d0*d - 4.0d0)
   L(4) = 1.0d0/6.0d0*(-3.0d0*d2 + 18.0d0*d - 18.0d0)
   L(5) = 1.0d0/24.0d0*(4.0d0*d3 - 48.0d0*d2 + 144.0d0*d - 96.0d0)
   L(6) = 1.0d0/120.0d0*(-5.0d0*d4 + 100.0d0*d3 - 600.0d0*d2 + 1200.0d0*d - 600.0d0)
   !if (ind > 6) L(7) = 1.0d0/720.0d0*(6.0d0*d5 - 180.0d0*d4 + 1800.0d0*d3 - 7200.0d0*d2 + 10800.0d0*d - 4320.0d0)
end subroutine d_Laguerre_up_to_6



pure function Moivre_cos(m, gamm) result(cos_mg)
   real(8) :: cos_mg
   integer, intent(in) :: m   ! coefficient
   real(8), intent(in) :: gamm ! variable
   !-------------------
   integer :: k, m2, k2
   real(8) :: summ, bcoef, cos_gam, sin_gam

   if (m == 0) then
      cos_mg = 1.0d0
   else
      m2 = FLOOR(abs(m)/2.0d0)
      cos_gam = cos(gamm)
      sin_gam = sin(gamm)
      summ = 0.0d0 ! to start with
      do k = 0, m2
         k2 = 2*k
         bcoef = binomial_coef(m, k2)   ! below
         summ = summ + (-1)**k * bcoef * cos_gam**(abs(m)-k2) * sin_gam**k2
      enddo
      cos_mg = summ
   endif
end function Moivre_cos



pure function Moivre_sin(m, gamm) result(sin_mg)
   real(8) :: sin_mg
   integer, intent(in) :: m   ! coefficient
   real(8), intent(in) :: gamm ! variable
   !-------------------
   integer :: k, m2, k2
   real(8) :: summ, bcoef, cos_gam, sin_gam

   if (m == 0) then
      sin_mg = 0.0d0
   else
      m2 = FLOOR((abs(m)-1)/2.0d0)
      cos_gam = cos(gamm)
      sin_gam = sin(gamm)
      summ = 0.0d0 ! to start with
      do k = 0, m2
         k2 = 2*k
         bcoef = binomial_coef(m, k2+1)   ! below
         summ = summ + (-1)**k * bcoef * cos_gam**(abs(m)-k2-1) * sin_gam**(k2+1)
      enddo
      sin_mg = summ
   endif
end function Moivre_sin



pure function binomial_coef(n, k) result(bc)
   real(8) :: bc  ! binomial coefficient
   integer, intent(in) :: k, n
   integer :: nk, i
   real(8) :: nk_f, k_f, n_f

   if ((n <= k) .or. (k <= 0))then
      bc = 1
   else
      nk = n - k
      k_f = 1.0d0    ! to start with
      nk_f = 1.0d0   ! to start with
      n_f = 1.0d0    ! to start with
      do i = 2, n
         n_f = n_f*dble(i) ! factorial
         if (i == k) then ! save this factorial
            k_f = n_f
         endif
         if (i == nk) then ! save this factorial
            nk_f = n_f
         endif
      enddo
      ! Collect the terms of the binomial coefficient:
      bc = n_f/(k_f * nk_f)
   endif
end function binomial_coef



pure function get_factorial(n, step) result(fact)   ! factorial of an integer number "n"
   real(8) :: fact
   integer, intent(in) :: n
   integer, intent(in), optional :: step    ! factorial step
   integer :: i, i_start, i_by
   
   fact = 1.0d0 ! to start from
   if (n > 1) then
      i_start = 2
      if (present (step)) then
         i_by = step
      else
         i_by = 1
      endif
      do i = i_start, n, i_by
         fact = fact*dble(i)
      enddo
   endif
end function get_factorial
 

pure function double_factorial(n)
   real(8) :: double_factorial
   integer, intent(in) :: n
   integer :: i, i_start
   double_factorial = 1.0d0
   if (n > 1) then
      if (mod( n, 2 ) == 0) then ! odd number
         i_start = 1
      else ! even number
         i_start = 2
      endif
      do i = i_start, n, 2
         double_factorial = double_factorial*dble(i)
      enddo
   endif
end function double_factorial
 
 

subroutine fit_parabola_to_3points(x1, y1, x2, y2, x3, y3, A, B, C)  ! find coefficients of parabola by given three points
   real(8), intent(in) :: x1, x2, x3    ! given three points x coordinates
   real(8), intent(in) :: y1, y2, y3    ! function values for the given three points
   real(8), intent(out) :: A, B, C      ! parabola coefficients
   real(8) :: denom, x1x2, x2x3, x1x3, x3y1y2, x2y1y3, x1y2y3
   x1x2 = x1 - x2
   x2x3 = x2 - x3
   x1x3 = x1 - x3
   x3y1y2 = x3*(y1 - y2)
   x2y1y3 = x2*(y1 - y3)
   x1y2y3 = x1*(y2 - y3)
   denom = x1x2*x1x3*x2x3
   A = (-x3y1y2 + x2y1y3 - x1y2y3) / denom
   B = (x3*x3y1y2 - x2*x2y1y3 + x1*x1y2y3) / denom
   C = (x2*x3*x2x3*y1 - x3*x1*x1x3*y2 + x1*x2*x1x2*y3) / denom
end subroutine fit_parabola_to_3points



subroutine fit_line_to_2points(x1, y1, x2, y2, A, B)  ! find coefficients of a line by given two points
   real(8), intent(in) :: x1, x2    ! given 2 points x coordinates
   real(8), intent(in) :: y1, y2    ! function values for the given 2 points
   real(8), intent(out) :: A, B      ! line coefficients
   real(8) :: x1x2
   x1x2 = x2 - x1
   A = (x2*y1 - x1*y2) / x1x2
   B = (y2 - y1) / x1x2
end subroutine fit_line_to_2points


 
pure function Kronecker_delta(i, j) result(delta)
   real(8) :: delta
   integer, intent(in) :: i, j
   if (i .EQ. j) then
      delta = 1.0d0
   else
      delta = 0.0d0
   endif
end function Kronecker_delta


pure function Heavyside_tau(m) result(tau)
   real(8) :: tau
   integer, intent(in) :: m
   if (m >= 0) then
      tau = 1.0d0
   else
      tau = 0.0d0
   endif
end function Heavyside_tau





subroutine nonsym_diagonalize_r(M, Ev, Error_descript, print_Ei, check_M)
   real(8), dimension(:,:), intent(inout) :: M	! matrix
   real(8), dimension(:), intent(out) :: Ev	! eigenvalues
   character(*), intent(inout) :: Error_descript	! error description
   logical, intent(in), optional :: print_Ei ! print out eigenvalues or not
   logical, intent(in), optional :: check_M  ! check the diagonalized matrix
   integer :: LWORK, N, INFO
   real(8), dimension(:), allocatable :: WORK
   real(8), dimension(:,:), allocatable :: M_save	! matrix
   real(8), dimension(:,:), allocatable :: VL	! matrix of left eigenvectors
   real(8), dimension(:,:), allocatable :: VR	! matrix of right eigenvectors
   real(8), dimension(:), allocatable :: Im_Ev	! imaginary parts of eigenvalues
   integer :: i_countin, FN, i, j

   N = size(M,1)
   allocate(M_save(N,N))
   allocate(VR(N,N))
   allocate(Im_Ev(N))
   Im_Ev = 0.0d0
   LWORK = 34*N
   allocate(WORK(LWORK))

   !$OMP WORKSHARE
   M_save = M ! save matrix before diagonalization just in case
   !$OMP END WORKSHARE

   call DGEEV('N','V', N, M, N, Ev, Im_Ev, VL, 1, VR, N, WORK, LWORK, INFO) ! library MKL (or LAPACK)

   !$OMP WORKSHARE
   M = VR ! save eigenvectors in the former matrix (Hamiltonian) assuming real values
   !$OMP END WORKSHARE

   if (INFO < 0) then ! if LAPACK diagonalization procidure failed:
      write(Error_descript,'(a,i5)') 'Module Algebra_tools: nonsym_diagonalize_r failed! Illigal value line ', ABS(INFO)
      print*, trim(adjustl(Error_descript))
   else if (INFO > 0) then
      write(Error_descript,'(a,i5)') 'Module Algebra_tools: nonsym_diagonalize_r failed! Did not converge, line ', ABS(INFO)
      print*, trim(adjustl(Error_descript))
   endif

   if (present(check_M)) then
      i = 1
      do while (i <= N)
         if (ABS(Im_Ev(i)) /= 0.0d0) then 
            !write(*,'(a,i4,f12.5,a,es12.5)') 'Eigenvalue #', i, Ev(i), ' has imaginary part ', Im_Ev(i)
            M(:,i)   = VR(:,i) !, VR(:,j+1)) ! v(j) = VR(:,j) + i*VR(:,j+1)
            M(:,i+1) = VR(:,i) !,-VR(:,j+1)) ! v(j+1) = VR(:,j) - i*VR(:,j+1)
            i = i + 2
         else
            M(:,i) = VR(:,i) ! only real part
            i = i + 1
         endif
      enddo
      call check_Ha(M_save, M, Ev)   ! to make sure diagonalization went well
   endif

   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif

   deallocate(Im_Ev)
   deallocate(M_save)
   deallocate(WORK)
   deallocate(VR)
end subroutine nonsym_diagonalize_r


subroutine nonsym_diagonalize_c(M, VR_c, Ev, Error_descript, print_Ei, check_M)
   real(8), dimension(:,:), intent(inout) :: M	! matrix
   complex, dimension(:,:), intent(out) :: VR_c	! matrix of complex right eigenvectors
   complex, dimension(:), intent(out) :: Ev	! eigenvalues
   character(*), intent(inout) :: Error_descript	! error description
   logical, intent(in), optional :: print_Ei ! print out eigenvalues or not
   logical, intent(in), optional :: check_M  ! check the diagonalized matrix
   integer :: LWORK, N, INFO
   real(8), dimension(:), allocatable :: WORK
   real(8), dimension(:,:), allocatable :: M_save	! matrix
   real(8), dimension(:,:), allocatable :: VL	! matrix of left eigenvectors
   real(8), dimension(:,:), allocatable :: VR	! matrix of right eigenvectors
   real(8), dimension(:), allocatable :: Re_Ev	! real parts of eigenvalues
   real(8), dimension(:), allocatable :: Im_Ev	! imaginary parts of eigenvalues
   integer :: i_countin, FN, i, j

   N = size(M,1)
   allocate(M_save(N,N))
   allocate(VR(N,N))
   allocate(Re_Ev(N))
   allocate(Im_Ev(N))
   Im_Ev = 0.0d0
   LWORK = 34*N
   allocate(WORK(LWORK))

   M_save = M ! save matrix before diagonalization just in case

   call DGEEV('N','V', N, M, N, Re_Ev, Im_Ev, VL, 1, VR, N, WORK, LWORK, INFO)

   if (INFO < 0) then ! if LAPACK diagonalization procidure failed:
      write(Error_descript,'(a,i5)') 'Module Algebra_tools: nonsym_diagonalize_r failed! Illigal value line ', ABS(INFO)
      print*, trim(adjustl(Error_descript))
   else if (INFO > 0) then
      write(Error_descript,'(a,i5)') 'Module Algebra_tools: nonsym_diagonalize_r failed! Did not converge, line ', ABS(INFO)
      print*, trim(adjustl(Error_descript))
   endif

   Ev = dcmplx(Re_Ev, Im_Ev) ! construct complex eigenvalues
   i = 1
   do while (i <= N)
      if (ABS(Im_Ev(i)) /= 0.0d0) then
         VR_c(:,i)   = dcmplx(VR(:,i), VR(:,i+1)) ! v(j) = VR(:,j) + i*VR(:,j+1)
         VR_c(:,i+1) = dcmplx(VR(:,i),-VR(:,i+1)) ! v(j+1) = VR(:,j) - i*VR(:,j+1)
         i = i + 2
      else
         VR_c(:,i) = dcmplx(VR(:,i),0.0d0) ! only real part
         i = i + 1
      endif
   enddo
   call sort_eigenvalues(Ev, VR_c) ! sort them accending

   if (present(check_M)) then
!       call check_Ha(M_save, VR_c, Ev)   ! make sure diagonalization and sorting went well
   endif

   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif

   deallocate(M_save)
   deallocate(VR)
   deallocate(Re_Ev)
   deallocate(Im_Ev)
   deallocate(WORK)
end subroutine nonsym_diagonalize_c



subroutine nonsym_diagonalize_cc(M, Ev, Error_descript, print_Ei, check_M)
   complex, dimension(:,:), intent(inout) :: M	! matrix
   complex, dimension(:), intent(out) :: Ev	! eigenvalues
   character(*), intent(inout) :: Error_descript	! error description
   logical, intent(in), optional :: print_Ei ! print out eigenvalues or not
   logical, intent(in), optional :: check_M  ! check the diagonalized matrix
   integer :: LWORK, N, INFO
   complex, dimension(:), allocatable :: WORK
   real(8), dimension(:), allocatable :: RWORK
   complex, dimension(:,:), allocatable :: M_save	! matrix
   complex, dimension(:,:), allocatable :: M_0	! matrix
   complex, dimension(:,:), allocatable :: VR_c	! matrix of complex right eigenvectors
   complex, dimension(:,:), allocatable :: VL_c	! matrix of complex right eigenvectors
   integer :: i_countin, FN, i, j

   N = size(M,1)
   allocate(M_save(N,N))
   allocate(M_0(N,N))
   allocate(VL_c(1,1))
   allocate(VR_c(N,N))
   LWORK = 2*N
   allocate(WORK(LWORK))
   allocate(RWORK(2*N))

   M_save = M ! save matrix before diagonalization just in case
   M_0 = M_save

   !call DGEEV('N','V', N, M_0, N, Re_Ev, Im_Ev, VL, 1, VR, N, WORK, LWORK, INFO) ! LAPACK
   call CGEEV('N','V', N, M, N, Ev, VL_c, 1, VR_c, N, WORK, LWORK, RWORK, INFO) ! LAPACK
   M = VR_c ! eigenvectors

   if (INFO < 0) then ! if LAPACK diagonalization procidure failed:
      write(Error_descript,'(a,i5)') 'Module Algebra_tools: nonsym_diagonalize_r failed! Illigal value line ', ABS(INFO)
      print*, trim(adjustl(Error_descript))
   else if (INFO > 0) then
      write(Error_descript,'(a,i5)') 'Module Algebra_tools: nonsym_diagonalize_r failed! Did not converge, line ', ABS(INFO)
      print*, trim(adjustl(Error_descript))
   endif

   call sort_eigenvalues(Ev, M) ! sort them accending, below

   if (present(check_M)) then
      call check_Ha(M_save, M, Ev) ! make sure diagonalization and sorting went well, below
   endif

   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif

   deallocate(M_save)
   deallocate(M_0)
   deallocate(VL_c)
   deallocate(VR_c)
   deallocate(WORK)
   deallocate(RWORK)
end subroutine nonsym_diagonalize_cc


subroutine sort_eigenvalues(Ev, Evec)
   complex, dimension(:), intent(inout) :: Ev	! eigenvalues
   complex, dimension(:,:), intent(inout), optional :: Evec	! eigenvectors
   real(8) :: temp
   complex :: temp_c
   complex, dimension(:), allocatable :: temp_c2
   integer N,i,j
   LOGICAL :: swapped
   N = size(Ev)
   allocate(temp_c2(N))
   do j = N-1, 1, -1
      swapped = .false. ! nothing swapped at the start
      do i = 1, j
         IF (real(Ev(i)) > real(Ev(i+1))) THEN
            ! Eigenvalue:
            temp_c = Ev(i)
            Ev(i) = Ev(i+1)
            Ev(i+1) = temp_c
            ! Eigenvector:
            if (present(Evec)) then
               temp_c2(:) = Evec(:,i)
               Evec(:,i) = Evec(:,i+1)
               Evec(:,i+1) = temp_c2(:)
            endif
            swapped = .true.  ! at least one pair of elements needed swapping
         END IF
      enddo
      IF (.NOT. swapped) EXIT
   enddo
end subroutine sort_eigenvalues


subroutine sort_array_r(array_in)  ! bubble sorting algorithm for real 1d array
   real(8), dimension(:), intent(inout) :: array_in
   real(8) :: temp
   integer N,i,j
   logical :: swapped
   N = size(array_in)
   do j = N-1, 1, -1
      swapped = .false. ! nothing swapped at the start
      do i = 1, j
         if (array_in(i) > array_in(i+1)) then ! swap elements
            temp = array_in(i)
            array_in(i) = array_in(i+1)
            array_in(i+1) = temp
            swapped = .true.  ! at least one pair of elements needed swapping
         end if
      enddo
      if (.not. swapped) exit
   enddo
end subroutine sort_array_r

subroutine sort_array_c(array_in)  ! bubble sorting algorithm for complex 1d array (sorted by real part)
   complex, dimension(:), intent(inout) :: array_in
   complex :: temp
   integer N,i,j
   logical :: swapped
   N = size(array_in)
   do j = N-1, 1, -1
      swapped = .false. ! nothing swapped at the start
      do i = 1, j
         if (real(array_in(i)) > real(array_in(i+1))) then ! swap elements
            temp = array_in(i)
            array_in(i) = array_in(i+1)
            array_in(i+1) = temp
            swapped = .true.  ! at least one pair of elements needed swapping
         end if
      enddo
      if (.not. swapped) exit
   enddo
end subroutine sort_array_c



recursive subroutine QsortC(A)
  real, intent(inout), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real, intent(inout), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine Partition



subroutine r_diagonalize(M, Ev, Error_descript, print_Ei, check_M, use_DSYEV)
   real(8), dimension(:,:), intent(inout) :: M	! matrix
   real(8), dimension(:), intent(out) :: Ev	! eigenvalues
   character(*), intent(inout) :: Error_descript	! error description
   logical, intent(in), optional :: print_Ei ! print out eigenvalues or not
   logical, intent(in), optional :: check_M  ! check the diagonalized matrix
   logical, intent(in), optional :: use_DSYEV	! do not use divide-and-conquare algorithm
   integer :: LIWORK, LWORK, N, INFO
   real(8), dimension(:), allocatable :: LAPWORK
   INTEGER, dimension(:), allocatable :: IWORK
   !real(8), dimension(size(M,1),size(M,2)) :: M_save	! matrix
   real(8), dimension(:,:), allocatable :: M_save	! matrix
   integer :: i_countin, FN, i, j

   N = size(M,1)
   allocate(M_save(N,N))
   LIWORK = 30 + 5*N
   LWORK = 10 + 6*N + 2*N*N
   
   allocate(LAPWORK(LWORK))
   allocate(IWORK(LIWORK))

   !$OMP WORKSHARE
   M_save = M ! save matrix before diagonalization just in case
   !$OMP END WORKSHARE

   if (.not.present(use_DSYEV)) then
      call dsyevd('V','U', N, M, N, Ev, LAPWORK, LWORK, IWORK, LIWORK, INFO) ! LAPACK
!       call dsyevd_2stage('V','U', N, M, N, Ev, LAPWORK, LWORK, IWORK, LIWORK, INFO) ! LAPACK
   else
      call DSYEV('V', 'U', N, M, N, Ev, LAPWORK, LWORK, INFO)
   endif
   
   if (INFO .NE. 0) then ! if divide-n-conquare LAPACK diagonalization procidure failed, try regular one:
      !$OMP WORKSHARE
      M = M_save
      !$OMP END WORKSHARE
      call DSYEV('V', 'U', N, M, N, Ev, LAPWORK, LWORK, INFO)
      if (INFO .NE. 0) then ! if LAPACK diagonalization procidure failed:
         Error_descript = 'Module Algebra_tools: real matrix diagonalization failed!'
         print*, trim(adjustl(Error_descript))
      endif
      if (present(check_M)) then
         call check_Ha(M_save, M, Ev)   ! to make sure diagonalization went well
      endif
   endif
   
   if (present(check_M)) then
      call check_Ha(M_save, M, Ev)   ! to make sure diagonalization went well
   endif

   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif
end subroutine r_diagonalize


subroutine check_Ha_r(Mat, Eigenvec, Eigenval) ! real matrix real eigenstates
   real(8), dimension(:,:), intent(in) :: Mat	! matrix
   real(8), dimension(:,:), intent(in) :: Eigenvec	! eigenvectors of this matrix
   real(8), dimension(:), intent(in) :: Eigenval	! eigenvalues of this matrix
   real(8), dimension(size(Mat,1), size(Mat,2)) :: M_temp
   real(8) Ev, Evec(size(Eigenval))
   integer i, j, k, M

   M = size(Eigenval)
   Ev = 0.0d0
   do i = 1,M
      do k = 1,M
!          Evec(k) = SUM(Mat(k,:)*Eigenvec(:,i)) ! also correct, for Hermitian matrices
          Evec(k) = SUM(Mat(:,k)*Eigenvec(:,i))
      enddo
      Ev = SUM(Evec(:)*Eigenvec(:,i))
      if (ABS(Ev - Eigenval(i))/min(ABS(Ev),ABS(Eigenval(i))) .GT. 1.0d-10) then ! diagonalization went wrong:
         print*, 'Ev:', i, Ev, Eigenval(i)
      endif
   enddo
end subroutine check_Ha_r



subroutine get_eigenvalues_from_eigenvectors(Mat, Eigenvec, Eigenval) ! real matrix real eigenstates
   real(8), dimension(:,:), intent(in) :: Mat   ! matrix
   real(8), dimension(:,:), intent(in) :: Eigenvec ! eigenvectors of this matrix
   real(8), dimension(:), intent(out) :: Eigenval  ! eigenvalues of this matrix
   real(8) Evec(size(Eigenval))
   integer i, k, Nsiz

   Nsiz = size(Eigenval)   ! size of eigenvalues array

   ! 1) Calculate all the eigenvalues from the eigenvectors of the matrix:
   do i = 1,Nsiz
      do k = 1,Nsiz
          Evec(k) = SUM(Mat(:,k)*Eigenvec(:,i))
      enddo
      Eigenval(i) = SUM(Evec(:)*Eigenvec(:,i))
   enddo

   ! 2) Sort the eigenvalues array in increasing order:
   call sort_array(Eigenval)  ! above
end subroutine get_eigenvalues_from_eigenvectors



subroutine c_diagonalize(M, Ev, Error_descript, print_Ei, check_M)
   complex, dimension(:,:), intent(inout) :: M	! matrix
   real(8), dimension(:), intent(out) :: Ev	! eigenvalues
   character(*), intent(inout) :: Error_descript	! error description
   logical, intent(in), optional :: print_Ei ! print out eigenvalues or not
   logical, intent(in), optional :: check_M ! chech diagonalization
   !-------------------------
   integer :: LIWORK, LWORK, N, INFO, LRWORK
!    complex(16), dimension(:), allocatable :: LAPWORK
   complex, dimension(:), allocatable :: LAPWORK
   integer, dimension(:), allocatable :: IWORK
   real(8), dimension(:), allocatable ::  RWORK
   !real, dimension(:), allocatable ::  RWORK
   complex, dimension(:,:), allocatable :: M_save, M_work   ! matrix
!    complex, dimension(size(M,1),size(M,2)) :: M_save	! matrix
!    complex, dimension(size(M,1),size(M,2)) :: M_work	! matrix in 16 bit format
!    complex(16), dimension(size(M,1),size(M,2)) :: M_work	! matrix in 16 bit format :: WRONG!
   
   integer :: i_countin, FN, i, j
   N = size(M,1)
   allocate(M_save(N,N))
   allocate(M_work(N,N))
   LIWORK = 3 + 5*N  ! 3 + 5*N
   LRWORK = 1 + 5*N + 2*N*N !  1 + 5*N + 2*N**2
   LWORK = 2*N + N*N ! 2*N + N**2
   allocate(LAPWORK(MAX(1,LWORK))) ! (MAX(1,LWORK))
   allocate(IWORK(LIWORK)) ! (MAX(1,LIWORK))
   allocate(RWORK(LRWORK)) ! LRWORK

   !!$OMP WORKSHARE
   M_save = M ! to make sure it's fine
   M_work = M ! convert it into complex*16
   !!$OMP END WORKSHARE
   
   ! Test it before the start:
   CH1:do i = 1, size(M,1) 
      do j = 1, size(M,2) 
         if (isnan(dble(M_work(i,j))) .or. isnan(AIMAG(M_work(i,j)))) then
            Error_descript = 'Module Algebra_tools: complex matrix diagonalization: M is NaN before'
            print*, trim(adjustl(Error_descript))
            exit CH1
         endif
      enddo ! j
   enddo CH1 ! i
   
   !call ZHEEVD('V','U', N, M, N, Ev, LAPWORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
   call ZHEEVD('V','U', N, M_work, N, Ev, LAPWORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
   !print*, OMP_GET_THREAD_NUM(), 'After ZHEEVD'
   
   if (INFO .NE. 0) then ! if divide-n-conquare LAPACK diagonalization procidure failed, try regular one:
      print*, 'ZHEEVD did not work, INFO=', INFO
      !!$OMP WORKSHARE
      M = M_save ! to make sure it's fine
      M_work = M_save
      !!$OMP END WORKSHARE
      !call ZHEEV('V', 'U', N, M, N, Ev, LAPWORK, LWORK, RWORK,INFO)
      call ZHEEV('V', 'U', N, M_work, N, Ev, LAPWORK, LWORK, RWORK,INFO)
      if (INFO .NE. 0) then ! if LAPACK diagonalization procidure failed:
         Error_descript = 'Module Algebra_tools: complex matrix diagonalization failed!'
         print*, trim(adjustl(Error_descript))
      endif
   endif
   deallocate(LAPWORK, IWORK, RWORK)
   
   !!$OMP WORKSHARE
   M = M_work ! save processed matrix back into the output matrix
   !!$OMP END WORKSHARE
   
   if (present(check_M)) then
      if (check_M) then
         call check_Ha(M_save, M_work, Ev)   ! to make sure diagonalization went well
!          pause 'check_M'
      endif
   endif
   deallocate(M_save, M_work) ! clean up
   
   ! test the resulting diagonalized Hamiltonian:
   CH2:do i = 1, size(M,1)
      if (isnan(Ev(i))) then
         Error_descript = 'Module Algebra_tools: complex matrix diagonalization: Ev is NaN after'
         print*, trim(adjustl(Error_descript)), INFO   
         exit CH2
      endif
      do j = 1, size(M,2) 
         if (isnan(dble(M(i,j))) .or. isnan(AIMAG(M(i,j)))) then
            Error_descript = 'Module Algebra_tools: complex matrix diagonalization: M is NaN after'
            print*, trim(adjustl(Error_descript))
            exit CH2
         endif
      enddo ! j
   enddo CH2 ! i
   
   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif

   !print*, OMP_GET_THREAD_NUM(), 'c_diagonalize done'
end subroutine c_diagonalize




subroutine c8_diagonalize(M, Ev, Error_descript, print_Ei, check_M) ! double precision complex
   complex(8), dimension(:,:), intent(inout) :: M  ! matrix
   real(8), dimension(:), intent(inout) :: Ev      ! eigenvalues
   character(*), intent(inout) :: Error_descript   ! error description
   logical, intent(in), optional :: print_Ei       ! print out eigenvalues or not
   logical, intent(in), optional :: check_M        ! chech diagonalization
   !-------------------------
   integer :: LIWORK, LWORK, N, INFO, LRWORK
   complex(8), dimension(:), allocatable :: LAPWORK
   integer, dimension(:), allocatable :: IWORK
   real(8), dimension(:), allocatable ::  RWORK
   !complex(8), dimension(size(Ev),size(Ev)) :: M_save ! This way of defining size breaks down!
   complex(8), dimension(:,:), allocatable :: M_save, M_work
   integer :: i_countin, FN, i, j

   N = size(M,1)
   allocate(M_save(N,N))
   allocate(M_work(N,N))

   LIWORK = 3 + 5*N  ! 3 + 5*N
   LRWORK = 1 + 5*N + 2*N*N !  1 + 5*N + 2*N**2
   LWORK = 2*N + N*N ! 2*N + N**2
   allocate(LAPWORK(MAX(1,LWORK))) ! (MAX(1,LWORK))
   allocate(IWORK(LIWORK)) ! (MAX(1,LIWORK))
   allocate(RWORK(LRWORK)) ! LRWORK

   !$OMP WORKSHARE
   M_save = M ! to make sure it's fine
   M_work = M ! convert it into complex*16
   !$OMP END WORKSHARE

   ! Test it before the start:
   CH1:do i = 1, N
      do j = 1, N
         if (isnan(dble(M_work(i,j))) .or. isnan(AIMAG(M_work(i,j)))) then
            Error_descript = 'Module Algebra_tools: complex matrix diagonalization: M is NaN before'
            print*, trim(adjustl(Error_descript))
            exit CH1
         endif
      enddo ! j
   enddo CH1 ! i

   call ZHEEVD('V','U', N, M_work, N, Ev, LAPWORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

   if (INFO .NE. 0) then ! if divide-n-conquare LAPACK diagonalization procidure failed, try regular one:
      print*, 'ZHEEVD did not work, INFO=', INFO
      !$OMP WORKSHARE
      M = M_save ! to make sure it's fine
      M_work = M_save
      !$OMP END WORKSHARE
      !call ZHEEV('V', 'U', N, M, N, Ev, LAPWORK, LWORK, RWORK,INFO)
      call ZHEEV('V', 'U', N, M_work, N, Ev, LAPWORK, LWORK, RWORK,INFO)
      if (INFO .NE. 0) then ! if LAPACK diagonalization procidure failed:
         Error_descript = 'Module Algebra_tools: complex matrix diagonalization failed!'
         print*, trim(adjustl(Error_descript))
      endif
   endif

   deallocate(LAPWORK)
   deallocate(IWORK)
   deallocate(RWORK)

   !!$OMP WORKSHARE
   M = M_work ! save processed matrix back into the output matrix
   !!$OMP END WORKSHARE

   print*, OMP_GET_THREAD_NUM(), 'Before check_M'

   if (present(check_M)) then
      if (check_M) then
         call check_Ha_c8(M_save, M_work, Ev)   ! to make sure diagonalization went well
!          pause 'check_M'
      endif
   endif

   print*, OMP_GET_THREAD_NUM(), 'After check_M'

   ! Clean up temporary arrays:
   deallocate(M_save, M_work)

   ! test the resulting diagonalized Hamiltonian:
   CH2:do i = 1, N
      if (isnan(Ev(i))) then
         Error_descript = 'Module Algebra_tools: complex matrix diagonalization: Ev is NaN after'
         print*, trim(adjustl(Error_descript)), INFO
         exit CH2
      endif
      do j = 1, N
         if (isnan(dble(M(i,j))) .or. isnan(AIMAG(M(i,j)))) then
            Error_descript = 'Module Algebra_tools: complex matrix diagonalization: M is NaN after'
            print*, trim(adjustl(Error_descript))
            exit CH2
         endif
      enddo ! j
   enddo CH2 ! i

   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif
end subroutine c8_diagonalize



subroutine check_Ha_c(Mat, Eigenvec, Eigenval) ! real matrix, complex eigenstates
   complex, dimension(:,:), intent(in) :: Mat	! matrix
   complex, dimension(:,:), intent(in) :: Eigenvec	! eigenvectors of this matrix
   real(8), dimension(:), intent(in) :: Eigenval	! eigenvalues of this matrix
   real(8), dimension(size(Mat,1), size(Mat,2)) :: M_temp
   complex Evec(size(Eigenval))
   real(8) :: Ev
   integer i, j, k, M

   M = size(Eigenval)

   Ev = 0.0d0
   do i = 1,M
      do k = 1,M
          Evec(k) = SUM(Mat(k,:)*Eigenvec(:,i)) ! correct
!           Evec(k) = SUM(Mat(:,k)*Eigenvec(:,i))
      enddo
!       Ev = dble(SUM(DCONJG(Eigenvec(:,i))*Evec(:)))
      Ev = dble(SUM(CONJG(Eigenvec(:,i))*Evec(:)))
      if (ABS(Ev - Eigenval(i)) .GT. min(ABS(Ev),ABS(Eigenval(i)))*1.0d-10) then ! diagonalization went wrong:
         print*, 'Ev_c:', i, Ev, Eigenval(i) !, SUM(CONJG(Eigenvec(:,i))*Eigenvec(:,i))
      endif
   enddo
end subroutine check_Ha_c



subroutine check_Ha_c8(Mat, Eigenvec, Eigenval) ! real matrix, complex eigenstates
   complex(8), dimension(:,:), intent(in) :: Mat	! matrix
   complex(8), dimension(:,:), intent(in) :: Eigenvec	! eigenvectors of this matrix
   real(8), dimension(:), intent(in) :: Eigenval	! eigenvalues of this matrix
   real(8), dimension(size(Mat,1), size(Mat,2)) :: M_temp
   complex(8) :: Evec(size(Eigenval))
   real(8) :: Ev
   integer :: i, j, k, M

   M = size(Eigenval)

   Ev = 0.0d0
   do i = 1,M
      do k = 1,M
          Evec(k) = SUM(Mat(k,:)*Eigenvec(:,i)) ! correct
!           Evec(k) = SUM(Mat(:,k)*Eigenvec(:,i))
      enddo
!       Ev = dble(SUM(DCONJG(Eigenvec(:,i))*Evec(:)))
      Ev = dble(SUM(CONJG(Eigenvec(:,i))*Evec(:)))
      if (ABS(Ev - Eigenval(i))/min(ABS(Ev),ABS(Eigenval(i))) .GT. 1.0d-10) then ! diagonalization went wrong:
         print*, 'Ev_c:', i, Ev, Eigenval(i) !, SUM(CONJG(Eigenvec(:,i))*Eigenvec(:,i))
      endif
   enddo
end subroutine check_Ha_c8



subroutine check_Ha_cc(Mat, Eigenvec, Eigenval) ! complex matrix, complex eigenstates
   complex, dimension(:,:), intent(in) :: Mat	! matrix
   complex, dimension(:,:), intent(in) :: Eigenvec	! eigenvectors of this matrix
   complex, dimension(:), intent(in) :: Eigenval	! eigenvalues of this matrix
   complex, dimension(size(Mat,1), size(Mat,2)) :: M_temp
   complex Ev, Evec(size(Eigenval))
   integer i, j, k, M

   M = size(Eigenval)
   Ev = 0.0d0
   do i = 1,M
      do k = 1,M
         Evec(k) = SUM(Mat(k,:)*Eigenvec(:,i)) ! correct
         !Evec(k) = SUM(Mat(:,k)*Eigenvec(:,i))
      enddo
      Ev = SUM(CONJG(Eigenvec(:,i))*Evec(:))
      if (ABS(Ev - Eigenval(i))/min(ABS(Ev),ABS(Eigenval(i))) .GT. 1.0d-10) then ! diagonalization went wrong:
         print*, 'Ev_cc:', i, Ev, Eigenval(i)
      endif
   enddo
end subroutine check_Ha_cc


subroutine check_symmetry(Ha)
   real(8), dimension(:,:), intent(in) :: Ha	! matrix to check
   integer i, j
   !$OMP PARALLEL DO PRIVATE(i,j), SHARED(Ha)
   do i = 1, size(Ha,1)
      do j = i, size(Ha,2)
         if ((i /= j) .and. (Ha(i,j) /= 0.0d0)) then
            if ( ABS((Ha(i,j)-Ha(j,i))/Ha(i,j)) >= 1.0d-8) then
               write(*,'(a, i4, i4, es25.16, es25.16)') 'Nonsymmetric element found:', i, j, Ha(i,j), Ha(j,i)
            else
!              write(*,'(a, i4, i4, es25.16, es25.16)') 'Symmetric element:', i, j, Ha(i,j), Ha(j,i)
            endif
         endif
      enddo
   enddo
   !$OMP END PARALLEL DO
end subroutine check_symmetry


subroutine check_hermiticity(Ha)
   complex, dimension(:,:), intent(in) :: Ha	! matrix to check
   integer i, j
   real(8) :: epsylon
   complex temp
   epsylon = 1.0d-8
   !$OMP PARALLEL DO PRIVATE(i,j,temp), SHARED(Ha)
   do i = 1, size(Ha,1)
      do j = i, size(Ha,2)
!          if (i /= j) then
            temp = Ha(i,j) - CONJG(Ha(j,i))
            R:if (ABS(dble(Ha(i,j))) > epsylon) then   ! real part
               if ( ABS(dble(temp)/dble(Ha(i,j))) >= epsylon) then
                  write(*,'(a, i4, i4, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') 'Nonhermitian element Re:', i, j, Ha(i,j), Ha(j,i), temp
               endif
            endif R
            IR:if (ABS(aimag(Ha(i,j))) > epsylon) then   ! imaginary part
               if ( ABS(aimag(temp)/aimag(Ha(i,j))) >= epsylon) then
                  write(*,'(a, i4, i4, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') 'Nonhermitian element Im:', i, j, Ha(i,j), Ha(j,i), temp
               endif
            endif IR
!          endif !  (i /= j)
      enddo ! j = i, size(Ha,2)
   enddo ! i = 1, size(Ha,1)
   !$OMP END PARALLEL DO
end subroutine check_hermiticity


! This subroutine multiplies two real matrices:
subroutine mkl_matrix_mult_r(TRANSA, TRANSB, A, B, ResultM)
   !DGEMM  performs one of the matrix-matrix operations
   !C := alpha*op( A )*op( B ) + beta*C,
   !where  op( X ) is one of
   !op( X ) = X   or   op( X ) = X**T,
   real(8), dimension(:,:), intent(in) :: A, B
   real(8), dimension(:,:), intent(inout) ::ResultM
   character(1), intent(in) :: TRANSA, TRANSB
   real(8) :: ALPHA, BETA
   integer :: M, N, K, ka, kb
   integer :: LDA, LDB, LDC
   !DDDDDDDDDDDDDDDDDDDDDDDDD

   ALPHA = 1.0d0
   BETA = 0.0d0

   M = size(A,1)
   N = size(A,2)
   K = size(B,1)
   LDA = max(1,K,M)
   LDB = max(1,K,M)
   LDC = max(1,N)
#ifdef CUBLAS	
! if CUBLAS library for GPU computing is used
!    print*, 'CUBLAS'
   CALL cublas_dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, ResultM, LDC)
#else	! if standard MKL or BLAS library is used
   CALL dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, ResultM, LDC)
#endif
end subroutine mkl_matrix_mult_r


! This subroutine multiplies two complex matrices:
subroutine mkl_matrix_mult_c(TRANSA, TRANSB, A, B, ResultM)
   !DGEMM  performs one of the matrix-matrix operations
   !C := alpha*op( A )*op( B ) + beta*C,
   !where  op( X ) is one of
   !op( X ) = X   or   op( X ) = X**T,
   complex, dimension(:,:), intent(in) :: A, B
   complex, dimension(:,:), intent(inout) :: ResultM
   character(1), intent(in) :: TRANSA, TRANSB
   complex :: ALPHA, BETA
   integer :: M, N, K, ka, kb
   integer :: LDA, LDB, LDC
   !DDDDDDDDDDDDDDDDDDDDDDDDD
   !print*, 'mkl_matrix_mult_c start'

   ALPHA = dcmplx(1.0d0,0.0d0)
   BETA = dcmplx(0.0d0,0.0d0)

   M = size(A,1)
   N = size(B,2)
   K = size(A,2)
   !LDA = max(1,K,M)
   select case (TRANSA)
   case ('N', 'n')
      LDA = max(1,M)
   case default
      LDA = max(1,K)
   endselect
   !LDB = max(1,K,N)
   select case (TRANSB)
   case ('N', 'n')
      LDB = max(1,K)
   case default
      LDB = max(1,N)
   endselect
   LDC = max(1,M)

   !print*, 'calling zgemm '
   !CALL cgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, ResultM, LDC)  ! LAPACK, single precision
   CALL zgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, ResultM, LDC) ! LAPACK, double precision

   !print*, 'mkl_matrix_mult_c done'
end subroutine mkl_matrix_mult_c



subroutine mkl_matrix_mult_c8(TRANSA, TRANSB, A, B, ResultM)
   !DGEMM  performs one of the matrix-matrix operations
   !C := alpha*op( A )*op( B ) + beta*C,
   !where  op( X ) is one of
   !op( X ) = X   or   op( X ) = X**T,
   complex(8), dimension(:,:), intent(in) :: A, B
   complex(8), dimension(:,:), intent(inout) :: ResultM
   character(1), intent(in) :: TRANSA, TRANSB
   complex :: ALPHA, BETA
   integer :: M, N, K, ka, kb
   integer :: LDA, LDB, LDC
   !DDDDDDDDDDDDDDDDDDDDDDDDD

   ALPHA = cmplx(1.0d0,0.0d0)
   BETA = cmplx(0.0d0,0.0d0)

   M = size(A,1)
   N = size(B,2)
   K = size(A,2)
   LDA = max(1,K,M)
   LDB = max(1,K,N)
   LDC = max(1,M)
   CALL zgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, ResultM, LDC)
end subroutine mkl_matrix_mult_c8


! Create a matrix of reciprocal vectors:
subroutine Reciproc(M1, R1)
   real(8), dimension(3,3), intent(in) :: M1    ! matrix of vectors
   real(8), dimension(3,3), intent(out) :: R1   ! matrix of reciprocal-space vectors
   real(8), dimension(3) :: vec
   real(8) :: Vol, Pi2_Vol

   ! Get volume:
   call Det_3x3(M1,Vol) ! below

   if (isnan(Vol)) then
      print*, 'Error #1 in Reciproc, undefined volume: ', Vol
      R1(:,:) = 0.0d0
   elseif (Vol == 0.0d0) then
      print*, 'Error #2 in Reciproc, zero volume: ', Vol
      R1(:,:) = 0.0d0
   else
      Pi2_Vol = 2.0d0 * g_Pi / Vol

      call Cross_Prod(M1(2,:), M1(3,:), vec) ! below
      R1(1,:) = Pi2_Vol * vec(:)

      call Cross_Prod(M1(3,:), M1(1,:), vec) ! below
      R1(2,:) = Pi2_Vol * vec(:)

      call Cross_Prod(M1(1,:), M1(2,:), vec) ! below
      R1(3,:) = Pi2_Vol * vec(:)
   endif
end subroutine Reciproc


subroutine Reciproc_old(M1, R1)
   real(8), dimension(3,3), intent(in) :: M1 ! matrix of vectors
   real(8), dimension(3,3), intent(out) :: R1 ! matrix of reciprocal-space vectors
   real(8), dimension(3) :: vec
   real(8) :: a
   call Cross_Prod(M1(2,:),M1(3,:),vec)
   call Dot_prod(M1(1,:),vec,a)
   R1(1,:) = 2.0d0*g_Pi*vec(:)/a
   call Cross_Prod(M1(3,:),M1(1,:),vec)
   call Dot_prod(M1(2,:),vec,a)
   R1(2,:) = 2.0d0*g_Pi*vec(:)/a
   call Cross_Prod(M1(1,:),M1(2,:),vec)
   call Dot_prod(M1(3,:),vec,a)
   R1(3,:) = 2.0d0*g_Pi*vec(:)/a
   if (isnan(R1(1,1))) print*, 'R1(1,1) is nan in "Reciproc"'
   if (isnan(R1(1,2))) print*, 'R1(1,2) is nan in "Reciproc"'
   if (isnan(R1(1,3))) print*, 'R1(1,3) is nan in "Reciproc"'
   if (isnan(R1(2,1))) print*, 'R1(2,1) is nan in "Reciproc"'
   if (isnan(R1(2,2))) print*, 'R1(2,2) is nan in "Reciproc"'
   if (isnan(R1(2,3))) print*, 'R1(2,3) is nan in "Reciproc"'
   if (isnan(R1(3,1))) print*, 'R1(3,1) is nan in "Reciproc"'
   if (isnan(R1(3,2))) print*, 'R1(3,2) is nan in "Reciproc"'
   if (isnan(R1(3,3))) print*, 'R1(3,3) is nan in "Reciproc"'
end subroutine Reciproc_old


! This subroutine calculates the cross-product of two 3-d vectors:
subroutine Cross_Prod(x, y, z)
   REAL(8), DIMENSION(:), INTENT(in) :: x, y ! input vectors of size 3
   REAL(8), DIMENSION(:), INTENT(out) :: z   ! cross-product, vector(3)
   if ((size(x) /= 3) .or. (size(y) /= 3) .or. (size(z) /= 3)) then
      print*, 'Subroutine Cross_Prod works only for vectors of length 3'
   else
      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
   endif
end subroutine Cross_Prod ! checked!

subroutine Dot_prod(x, y, z)
   REAL(8), DIMENSION(:), INTENT(in) :: x, y ! input vectors of size 3
   REAL(8), INTENT(out) :: z   ! dot-product, scalar
   z = SUM(x(:)*y(:))
   if (isnan(z)) print*, 'z is nan in "Dot_prod"'
end subroutine Dot_prod


! This subroutine calculates the determinant of a 3x3 matrix:
subroutine Det_3x3(A,detA)
   REAL(8), DIMENSION(3,3), INTENT(in) :: A ! input matrix
   REAL(8), INTENT(out) :: detA ! determinant of A, output
   detA = A(1,1)*( A(2,2)*A(3,3) - A(2,3)*A(3,2) ) - &
          A(1,2)*( A(2,1)*A(3,3) - A(2,3)*A(3,1) ) + &
          A(1,3)*( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
   if (detA .LE. 0.0d0) print*, 'Det A = 0 in Det_3x3'
end subroutine Det_3x3 ! checked!


! This subroutine calculates the transpose of a NxN matrix:
subroutine Transpose_M(M,TransM)
   REAL(8), DIMENSION(:,:), INTENT(in) :: M ! input matrix
   REAL(8), DIMENSION(:,:), INTENT(out) :: TransM ! output matrix, inverse of A(3x3)
!    integer i,j
!    do i = 1,size(M,1)
!       do j = 1,size(M,2)
!          TransM(i,j) = M(j,i)
!          print*, i, j, TransM(i,j)
!       enddo ! j
!    enddo ! i

   ! Use intrinsic fortran function:
   !$OMP WORKSHARE
   TransM = TRANSPOSE(M)
   !$OMP END WORKSHARE
end subroutine Transpose_M ! checked!


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
end subroutine Invers_3x3 ! checked!


! This subroutine calculates a product of NxN matrix with a N-dimensional vector:
subroutine Matrix_Vec_Prod(M,V,OutV)
   REAL(8), DIMENSION(:,:), INTENT(in) :: M ! input matrix (NxN)
   REAL(8), DIMENSION(:), INTENT(in) :: V ! input vector (size N)
   REAL(8), DIMENSION(:), INTENT(out) :: OutV ! output product, vector(N)
   integer i, j, n, l
   l = size(M,1)
   n = size(V)
   if (n .NE. l) print*,'(Martix x Vector) failed - different dimensions!'
   !$OMP WORKSHARE
   OutV = 0.0d0
   !$OMP END WORKSHARE
   !$omp PARALLEL DO private(i,j), SHARED(M,V)
   do i = 1,n
      do j = 1,n
         OutV(i) = OutV(i) + M(j,i)*V(j)
      enddo ! j 
   enddo ! i
   !$omp END PARALLEL DO
end subroutine Matrix_Vec_Prod ! checked



! This subroutine calculates a derivative of a determinant of matrix H(3x3) by element of a matrix h_a_b:
subroutine d_detH_d_h_a_b(H,a,b,dHdh)
   REAL(8), DIMENSION(3,3), INTENT(in) :: H ! input matrix (3x3)
   integer, INTENT(in) :: a ! number of raw of the element
   integer, INTENT(in) :: b ! number of column of the element
   REAL(8), INTENT(out) :: dHdh ! output, the derivative of the determinant
   dHdh = 0.0d0
   select case(a)
   case (1)
      select case (b)
      case (1)
         dHdh = H(2,2)*H(3,3) - H(2,3)*H(3,2)
      case (2)
         dHdh = -H(3,3)*H(2,1) + H(2,3)*H(3,1)
      case (3)
         dHdh = H(2,1)*H(3,2) - H(2,2)*H(3,1)
      end select
   case (2)
      select case (b)
      case (1)
         dHdh = -H(3,3)*H(1,2) + H(1,3)*H(3,2)
      case (2)
         dHdh = H(1,1)*H(3,3) - H(3,1)*H(1,3)
      case (3)
         dHdh = -H(1,1)*H(3,2) + H(3,1)*H(1,2)
      end select
   case (3)
      select case (b)
      case (1)
         dHdh = H(1,2)*H(2,3) - H(1,3)*H(2,2)
      case (2)
         dHdh = -H(1,1)*H(2,3) + H(1,3)*H(2,1)
      case (3)
         dHdh = H(1,1)*H(2,2) - H(1,2)*H(2,1)
      end select
   end select
end subroutine d_detH_d_h_a_b ! checked!



! This subroutine calculates a derivative of a determinant of matrix H(3x3) by element of a matrix h_a_b:
subroutine d_detH_d_h_a_b_OLD(H,a,b,dHdh)
   REAL(8), DIMENSION(3,3), INTENT(in) :: H ! input matrix (3x3)
   integer, INTENT(in) :: a ! number of raw of the element
   integer, INTENT(in) :: b ! number of column of the element
   REAL(8), INTENT(out) :: dHdh ! output, the derivative of the determinant
   dHdh = 0.0d0
   if ((a .EQ. 1) .AND. (b .EQ. 1)) dHdh = H(2,2)*H(3,3) - H(2,3)*H(3,2)
   if ((a .EQ. 1) .AND. (b .EQ. 2)) dHdh = -H(3,3)*H(2,1) + H(2,3)*H(3,1)
   if ((a .EQ. 1) .AND. (b .EQ. 3)) dHdh = H(2,1)*H(3,2) - H(2,2)*H(3,1)
   if ((a .EQ. 2) .AND. (b .EQ. 1)) dHdh = -H(3,3)*H(1,2) + H(1,3)*H(3,2)
   if ((a .EQ. 2) .AND. (b .EQ. 2)) dHdh = H(1,1)*H(3,3) - H(3,1)*H(1,3)
   if ((a .EQ. 2) .AND. (b .EQ. 3)) dHdh = -H(1,1)*H(3,2) + H(3,1)*H(1,2)
   if ((a .EQ. 3) .AND. (b .EQ. 1)) dHdh = H(1,2)*H(2,3) - H(1,3)*H(2,2)
   if ((a .EQ. 3) .AND. (b .EQ. 2)) dHdh = -H(1,1)*H(2,3) + H(1,3)*H(2,1)
   if ((a .EQ. 3) .AND. (b .EQ. 3)) dHdh = H(1,1)*H(2,2) - H(1,2)*H(2,1)
end subroutine d_detH_d_h_a_b_OLD ! checked!


! Construct diadic matrix out of two vectors: v1*v2^T = M
subroutine Two_Vect_Matr(V1,V2,M)
   REAL(8), DIMENSION(:), INTENT(in) :: V1 ! input vector 1 (column)
   REAL(8), DIMENSION(:), INTENT(in) :: V2 ! input vector 2 (raw)
   REAL(8), DIMENSION(:,:), INTENT(out) :: M ! output matrix
   integer i,j,n,l
   n = size(v1)
   l = size(v2)
!$omp PARALLEL private(i,j)
!$omp do
   do i = 1,n
      do j = 1,l
         M(i,j) = V1(i)*V2(j)
      enddo ! j
   enddo ! i
!$omp end do
!$omp end parallel
end subroutine Two_Vect_Matr ! checked


! This subroutine multiplies two sqare matricis:
subroutine Two_Matr_mult_r(M1,M2,Mout)
   REAL(8), DIMENSION(:,:), INTENT(in) :: M1 ! input matrix 1
   REAL(8), DIMENSION(:,:), INTENT(in) :: M2 ! input matrix 2
   REAL(8), DIMENSION(:,:), INTENT(out) :: Mout ! output matrix
   integer i,j,k,n,l,n2,l2
   n = size(M1,1)
   n2 = size(M2,1)
   l = size(M1,2)
   l2 = size(M2,2)
   if ((n .NE. n2)) print*, 'Two_Matr_mult Failed! (1) Ranges of matricis are different!'
   if ((n2 .NE. l)) print*, 'Two_Matr_mult Failed! (2) Ranges of matricis are different!'
   if ((l .NE. l2)) print*, 'Two_Matr_mult Failed! (3) Ranges of matricis are different!'
   Mout = 0.0d0
   do i = 1,l2
      do j = 1,l
         !Mout(i,j) = Mout(i,j) + SUM(M1(i,:)*M2(:,j)) ! incorrect
         Mout(i,j) = Mout(i,j) + SUM(M1(:,i)*M2(j,:)) ! correct
      enddo ! j
   enddo ! i
end subroutine Two_Matr_mult_r ! checked!


subroutine Two_Matr_mult_c(M1,M2,Mout)
   complex, DIMENSION(:,:), INTENT(in) :: M1 ! input matrix 1
   complex, DIMENSION(:,:), INTENT(in) :: M2 ! input matrix 2
   complex, DIMENSION(:,:), INTENT(out) :: Mout ! output matrix
   integer i,j,k,n,l,n2,l2
   n = size(M1,1)
   n2 = size(M2,1)
   l = size(M1,2)
   l2 = size(M2,2)
   if ((n .NE. n2)) print*, 'Two_Matr_mult_c Failed! (1) Ranges of matricis are different!'
   if ((n2 .NE. l)) print*, 'Two_Matr_mult_c Failed! (2) Ranges of matricis are different!'
   if ((l .NE. l2)) print*, 'Two_Matr_mult_c Failed! (3) Ranges of matricis are different!'
   Mout = 0.0d0
   do i = 1,l2
      do j = 1,l
         Mout(i,j) = Mout(i,j) + SUM(M1(i,:)*M2(:,j)) ! tested against MATMUL, correct
      enddo ! j
   enddo ! i
end subroutine Two_Matr_mult_c



END MODULE Algebra_tools
