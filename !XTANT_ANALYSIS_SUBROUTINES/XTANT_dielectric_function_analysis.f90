PROGRAM Dielectric_function_analysis
! Compilation:
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs XTANT_dielectric_function_analysis.f90 -o XTANT_dielectric_function_analysis.exe /link /stack:9999999999 
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays XTANT_dielectric_function_analysis.f90 -o XTANT_dielectric_function_analysis.exe /link /stack:9999999999 
!<===========================================
! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2022 Nikita Medvedev
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
! REFERENCES:
! [1] P. Yeh, "Optical waves in layered media"


use Universal_constants

character(200) :: File_name, File_name_2, File_name_out, char_var, File_name_out2
real(8) :: ReEps	 ! real part of the dielectric constant
real(8) :: ImEps	! imaginary part of the dielectric constant
real(8) :: n, k ! optical coefficients

real(8) :: n1, k1, n3, k3  ! optical coefficients of the first layer (vacuum) and substrate

! real(8) :: R, T, A ! reflectivity, transmittion, absorbtion
! real(8) :: R_all, T_all, A_all ! reflectivity, transmittion, absorbtion
real(8) :: R_s, T_s, A_s, R_p, T_p, A_p, R_s_all, T_s_all, A_s_all, R_p_all, T_p_all, A_p_all

real(8) :: ll, dd, teta0    ! wavelength, layer thikness, angle of incidence
real(8), dimension(9) :: EPS1, EPS0, EPS_cur
real(8) :: tim	! time [fs]
real(8) :: hw, hw0, hw1	! [eV]
integer :: FN, FN2, FN3, FN4	! file number
integer :: Reason, i, j
logical :: read_well

! Set defaults:
FN = 9999
FN2 = 9998
FN3 = 9997
FN4 = 9996

print*, '*********************************************************************'
print*, 'This program constructs optical coefficients from dielectric function calculated with XTANT'
print*, 'The following files must be present:'
print*, '1) OUTPUT_dielectric_function.dat - output file from XTANT calculations'
print*, '2) OPTICAL_PARAMETERS.txt - file with parameters of the target material'
print*, 'This file must contain the following lines (in this order):' 
print*, 'Layer thickness in [nm]'
print*, 'Probe incident angle in [degrees to normal]'
print*, 'Optical n and k of material above the layer (1.0 0.0 for vacuum or air)'
print*, 'Optical n and k of material below the layer (1.0 0.0 for vacuum or air)'
print*, '*********************************************************************'
print*, 'To execute the program, call:'
print*, 'XTANT_dielectric_function_analysis.exe  hw'
print*, 'Where hw is the probe photon energy in [eV]'
print*, '*********************************************************************'

File_name = 'OUTPUT_dielectric_function.dat'	! defaul name, for start
open (unit=FN, file=trim(adjustl(File_name)), status = 'old', readonly) 

File_name_2 = 'OPTICAL_PARAMETERS.txt'  ! file with parametrs of the probe pulse and material
open (unit=FN4, file=trim(adjustl(File_name_2)), status = 'old', readonly) 

!<===========================================
! get the files and parameters:
call get_add_data(File_name, hw, read_well)	! below
write(char_var,*) hw

! Set layer parameteres:
call set_layer_parameters(FN4, hw, ll, dd, teta0, n1, k1, n3, k3)   ! below


File_name_out =  'OUTPUT_dielectric_function_'//trim(adjustl( char_var ))//'.dat'
File_name_out2 =  'OUTPUT_dielectric_function_'//trim(adjustl( char_var ))//'_RTA.dat'
open (unit=FN2, file=trim(adjustl(File_name_out))) 
open (unit=FN3, file=trim(adjustl(File_name_out2))) 

print*, 'hw=', hw

write(FN2,'(a)') 'time  hw  Re_eps  Im_eps  LF  R   T   A   n   k'
write(FN3,'(a)') 'time   R_s    T_s   A_s  R_p T_p    A_p   R_s_all  T_s_all A_s_all    R_p_all   T_p_all  A_p_all '


EPS1 = 0.0d0
EPS0 = 0.0d0
EPS_cur = 0.0d0
i = 0
do
   i = i + 1
   !EPS1(:) = hw , Re_eps , Im_eps , LF , R , T , A  ,n , k
   read(FN,*,IOSTAT=Reason) EPS1(:)

   IF (Reason .LT. 0)  exit	! ... end of file reached ...
   
   if ( (EPS1(1) > hw) .and. (EPS0(1) <= hw) ) then ! here is the value we need
      forall (j=2:9) EPS_cur(j) = linear_interpolation(hw, EPS1(1), EPS0(1), EPS1(j), EPS0(j))
      write(FN2,'(f,es,es,es,es,es,es,es,es,es)') dble(i), EPS_cur(:)
      
!       call get_RTA_from_CDF(ReEps, ImEps, ll, dd, teta0, 1, R, T, A)    ! below
!       call get_RTA_from_CDF(ReEps, ImEps, ll, dd, teta0, 2, R_all, T_all, A_all)    ! below
      
      ! Get optical coefficients:
!       call get_R_T_A_total(hw, dd, dcmplx(1.0d0,0.0d0), teta0, dcmplx(EPS_cur(8), EPS_cur(9)), dcmplx(n3,k3), &
!                                        R_s, T_s, A_s, R_p, T_p, A_p, R_s_all, T_s_all, A_s_all, R_p_all, T_p_all, A_p_all)   ! below
      call get_R_T_A_total(hw, dd, dcmplx(n1,k1), teta0, dcmplx(EPS_cur(8), EPS_cur(9)), dcmplx(n3,k3), &
                                       R_s, T_s, A_s, R_p, T_p, A_p, R_s_all, T_s_all, A_s_all, R_p_all, T_p_all, A_p_all)   ! below
      
      write(FN3,'(f,es,es,es,es,es,es,es,es,es,es,es,es)') dble(i),  R_s, T_s, A_s, R_p, T_p, A_p, R_s_all, T_s_all, A_s_all, R_p_all, T_p_all, A_p_all
      !    write(*,*) dble(i), EPS_cur(:), 'a', EPS1(:)
!       pause
   endif
   
   ! for the next step
   EPS0 = EPS1
enddo

close(FN)
close(FN2)
close(FN3)

 !---------------------
 contains


subroutine set_layer_parameters(FN4, hw, ll, dd, teta0, n1, k1, n3, k3)
   integer, intent(in) :: FN4   ! file with parameters
   real(8), intent(in) :: hw    ! [eV]
   real(8), intent(out) :: ll, dd, teta0    ! [nm], [nm], angle of the probe pulse [radians]
   real(8), intent(out) :: n1, k1, n3, k3
   !------------------------------
   integer :: Reason
   ! Wavelength:
   ll = 2.0d0*g_Pi*g_h*g_cvel / (hw*g_e) * 1.0d9    ! [nm]
   ! Read other parameters from the file:
   read(FN4,*,IOSTAT=Reason) dd
   read(FN4,*,IOSTAT=Reason) teta0
   read(FN4,*,IOSTAT=Reason) n1, k1
   read(FN4,*,IOSTAT=Reason) n3, k3
   
!    print*, 'Enter the layer thickness in [nm]: '
!    read(*,*) dd
!    print*, 'Enter the angle of incidence of the prob pulse with respect to normal in [degree]: '
!    read(*,*) teta0
   ! Make it with respect to normal:
   teta0 = 90.0d0 - teta0
   ! Convert to radians:
   teta0 = teta0*g_Pi/(180.0d0) ! -> [radians]
   print*, 'The data you provided are: '
   print*, 'Thickness: ', dd, ' [nm]'
   print*, 'Angle: ', teta0, ' [rad]'
   print*, 'Probe wavelength: ', ll, ' [nm]'
end subroutine set_layer_parameters 
 
 

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW:


subroutine get_R_T_A_total(hw, dd, n1, teta1, n2, n3, R_s, T_s, A_s, R_p, T_p, A_p, R_s_all, T_s_all, A_s_all, R_p_all, T_p_all, A_p_all)
   real(8), intent(in) :: hw     ! [eV] photon energy
   real(8), intent(in) :: dd     ! [nm] layer thickness
   complex(8), intent(in) :: n1    ! real part of refractive index in media "1"
   real(8), intent(in) :: teta1  ! angle of light propagation with respect to normal in media 1
   complex(8), intent(in) :: n2    ! real part of refractive index in media "2"
   complex(8), intent(in) :: n3    ! real part of refractive index in media "3"
   real(8), intent(out) :: R_s, T_s, A_s, R_p, T_p, A_p     ! first ray reflection, transmission and absorption
   real(8), intent(out) :: R_s_all, T_s_all, A_s_all, R_p_all, T_p_all, A_p_all  ! all-rays reflection, transmission and absorption
   real(8) :: ll, w
   complex(8) ::  r12_s, t12_s, r12_p, t12_p, r23_s, t23_s, r23_p, t23_p, phi
   complex(8) :: k1x_out, k2x_out, k3x_out, R_amp, T_amp, A_amp
   
   ! Get frequency:
   w = frequency_from_enenrgy(hw)    ! below
!    print*, 'w=', w

   ! Get first layer coefficients:
   call get_r12_t12_complex(w, n1, teta1, n2, r12_s, t12_s, r12_p, t12_p, k1x_out=k1x_out, k2x_out=k2x_out)     ! below
   
   ! Get second layer coefficients:
   call get_r12_t12_complex(w, n2, teta1, n3, r23_s, t23_s, r23_p, t23_p, n0_in=n1, k2x_out=k3x_out)     ! below
   
   ! Get first ray coefficients:
   ! s polarization:
   R_s = dble(r12_s) * dble(r12_s)      ! Eq. (4.3-1) [1]
   !T_s = dble(t12_s) * dble(t12_s) * dble(k1x_out) / dble(k3x_out)  ! generalized Eq.(4.3-2) [1]
   T_s = dble(t12_s) * dble(t12_s) * dble(k3x_out) / dble(k1x_out)  ! generalized Eq.(4.3-2) [1]
   A_s = 1.0d0 - R_s - T_s      ! Eq.(4.3-3) [1]
   ! p polarization:
   R_p = dble(r12_p) * dble(r12_p)
   !T_p = dble(t12_p) * dble(t12_p) * dble(k1x_out) / dble(k3x_out)
   T_p = dble(t12_p) * dble(t12_p) * dble(k3x_out) / dble(k1x_out)
   A_p = 1.0d0 - R_p - T_p
   
   ! Get all rays coefficients:
   ! Get phi:
   phi = get_phi_complex(k2x_out, dd)  ! below
!    print*, 'phi=', phi, k2x_out, dd
!    print*, 'r23_p=', r23_p

   ! Get all rays reflection, transmission and absorption:
   ! s polarization:
   call total_R(r12_s, r23_s, t12_s, t23_s, phi, R_amp, T_amp)   ! below
   R_s_all= dble(R_amp) * dble(R_amp)
!    T_s_all = dble(T_amp) * dble(T_amp) * dble(k1x_out) / dble(k3x_out)
   T_s_all = dble(T_amp) * dble(T_amp) * dble(k3x_out) / dble(k1x_out)
   A_s_all= 1.0d0 - R_s_all - T_s_all
   ! p polarization:
   call total_R(r12_p, r23_p, t12_p, t23_p, phi, R_amp, T_amp)   ! below
   R_p_all = dble(R_amp) * dble(R_amp)
!    T_p_all = dble(T_amp) * dble(T_amp) * dble(k1x_out) / dble(k3x_out)
   T_p_all = dble(T_amp) * dble(T_amp) * dble(k3x_out) / dble(k1x_out)
   A_p_all = 1.0d0 - R_p_all - T_p_all
   
!    print*, 'R_p=', R_p, R_p_all
end subroutine get_R_T_A_total



pure subroutine total_R(r12, r23, t12, t23, Phi, R_all, T_all)
   complex(8), intent(in) :: r12, r23, t12,t23, Phi
   complex(8), intent(out) :: R_all, T_all
   complex(8) :: exp_arg1, exp_arg2, denom, R, T
   exp_arg1 = exp(-g_CI * Phi)
   exp_arg2 = exp(-2.0d0 * g_CI * Phi)
   denom = 1.0d0 / (1.0d0 + r12*r23*exp_arg2)
   R_all =  (r12 + r23*exp_arg2) * denom     ! Eq.(4.1-19) [1]
   T_all =  (t12 * t23 * exp_arg1) * denom      ! Eq.(4.1-18) [1]
end subroutine total_R



pure subroutine get_r12_t12_complex(w, n1, teta1, n2, r12_s, t12_s, r12_p, t12_p, n0_in, k1x_out, k2x_out)  ! for complex coefficients
   real(8), intent(in) :: w     ! [s] frequency of incident light
!    real(8), intent(in) :: dd    ! [nm] layer thickness
   complex(8), intent(in) :: n1    ! complex refractive index in media "1"
   real(8), intent(in) :: teta1  ! angle of light propagation with respect to normal in media 1
   complex(8), intent(in) :: n2    ! complex refractive index in media "2"
   complex(8), intent(out) :: r12_s, t12_s, r12_p, t12_p    ! first ray reflection and transmission amplitudes for s and p polarization
   complex(8), intent(in), optional :: n0_in    ! complex refraction index of the very first media (vacuum), if we are calculating interface inside of material
   complex(8), intent(out), optional :: k1x_out, k2x_out
   complex(8) :: k1z, k2z
   complex(8) :: k1x, k2x
   complex(8) :: k1xk2x, n1k2, n2k1
   ! Get z xomponents:
   if (present(n0_in)) then ! use the provided "air" values:
      k1z = get_kiz_complex(w, dble(n0_in), aimag(n0_in), teta1) ! below
   else     ! use the first media as the first one:
      k1z = get_kiz_complex(w, dble(n1), aimag(n1), teta1) ! below
   endif
   k2z = k1z    ! Eq.(3.4-2)
   ! Get x companents:
   k1x = get_kix_complex(w, dble(n1), aimag(n1), k1z) ! below
   k2x = get_kix_complex(w, dble(n2), aimag(n2), k2z) ! below
   
   if (present(k1x_out)) k1x_out = k1x
   if (present(k2x_out)) k2x_out = k2x
   
   ! Get s polarization coefficients:
   k1xk2x = 1.0d0 / (k1x + k2x)
   r12_s = (k1x - k2x) * k1xk2x     ! Eq.(5.1-13)
   t12_s = 2.0d0 * k1x * k1xk2x     ! Eq.(5.1-14)
   
   ! Get p polarization coefficients:
!    n1k2 = n1*n1*k2x
!    n2k1 = n2*n2*k1x
   n1k2 = dconjg(n1)*dconjg(n1)*k2x
   n2k1 = dconjg(n2)*dconjg(n2)*k1x
   k1xk2x = 1.0d0 / (n1k2 + n2k1)
   r12_p = (n1k2 - n2k1) * k1xk2x     ! Eq.(5.1-13)
   t12_p = 2.0d0 * n1k2 * k1xk2x     ! Eq.(5.1-14)
end subroutine get_r12_t12_complex



pure function get_phi_complex(kix, dd) result(Phi)
   complex(8) Phi
   complex(8), intent(in) :: kix    ! [1/m]
   real(8), intent(in) :: dd            ! [nm]
   Phi = kix * dd    ! Eq.(4.1-20) [1]
   Phi = Phi * 1.0d-9   ! converting dd from [nm] -> [m]
end function get_phi_complex


pure function get_kix_complex(w, ni, ki, kiz) result(kix) ! for complex index of refraction
   complex(8) kix    ! [1/m]
   real(8), intent(in) :: w     ! [s] frequency of incident light
   real(8), intent(in) :: ni    ! real part of refractive index in media "i"
   real(8), intent(in) :: ki    ! imaginary part of refractive index in media "i"
   complex(8), intent(in) :: kiz    ! z component
   real(8) wc
   complex(8) :: n
   wc = w/g_cvel
   n = dcmplx(ni,-ki)   ! definition of complex refractive index, Eq.(3.4-1) [1]
   kix = sqrt( wc*wc * (n*n) - kiz*kiz )  ! Eq.(3.4-5), P. 78 [1]
end function get_kix_complex


pure function get_kiz_complex(w, ni, ki, teta1) result(kiz)
   complex(8) kiz   ! [1/m]
   real(8), intent(in) :: w     ! [s] frequency of incident light
   real(8), intent(in) :: ni    ! real part of refractive index in media "i"
   real(8), intent(in) :: ki    ! imaginary part of refractive index in media "i"
   real(8), intent(in) :: teta1  ! angle of light propagation with respect to normal
   real(8) wc
   wc = w/g_cvel
   !kiz = dcmplx(ni * wc * sin(teta1), -ki) ! line before Eq.(3.4-5) [1]
   kiz = dcmplx(ni, -ki) * wc * sin(teta1) ! line before Eq.(3.4-2) [1]
end function get_kiz_complex



pure function Snell_sin_teta2(teta1, n1, n2) result(teta2)
    real(8) teta2   ! [rad] angle in media 2
    real(8), intent(in) :: teta1    ! [rad] incident angle in media 1
    real(8), intent(in) :: n1, n2   ! real refractive indices in media 1 and media 2
    real(8) :: sin_teta1, sin_teta2
    real(8) :: eps
    eps = 1.0d-10
    if (abs(n2) > eps) then
       sin_teta2 = sin_teta1 * n1 / n2
    endif
    if (abs(sin_teta2) <= 1.0d0) then   ! transmission into the media is possible
       teta2 = ASIN(sin_teta2)
    else
       teta2 = -1.0d10    ! absolute reflection
    endif
end function Snell_sin_teta2



subroutine get_n_k(ReEps, ImEps, n, k)
   real(8), intent(in) :: ReEps, ImEps	! real and imagenary parts of dialectric function 
   real(8), intent(out) :: n, k	! n and k optical coefficients
   n = sqrt((ReEps + sqrt(ReEps*ReEps+ImEps*ImEps))*0.5d0)
   k =sqrt((-ReEps + sqrt(ReEps*ReEps+ImEps*ImEps))*0.5d0)
end subroutine get_n_k

subroutine get_Re_and_Im(ReEps, ImEps, n, k)
   real(8), intent(out) :: ReEps, ImEps	! real and imagenary parts of dialectric function 
   real(8), intent(in) :: n, k	! n and k optical coefficients
   ReEps = n*n - k*k
   ImEps = 2.0d0*n*k
end subroutine get_Re_and_Im


pure function wavelength_from_frequency(w) result (ll)
   real(8) ll   ! [nm]
   real(8), intent(in) :: w ! [1/sec] frequency
   ll = 2.0d0*g_Pi*g_h*g_cvel / (hw*g_e) * 1.0d9    ! [nm]
end function wavelength_from_frequency


pure function frequency_from_enenrgy(hw) result(w)
   real(8) w    ! [1/sec] frequency of light
   real(8), intent(in) :: hw    ! [eV] photon energy
   w = hw * g_e / g_h    ! [1/sec] frequency
end function frequency_from_enenrgy





pure subroutine get_r12_t12(w, n1, teta1, n2, teta2, r12, t12)  ! for real coefficients
   real(8), intent(in) :: w     ! [s] frequency of incident light
   real(8), intent(in) :: n1    ! real part of refractive index in media "1"
   real(8), intent(in) :: teta1  ! angle of light propagation with respect to normal in media 1
   real(8), intent(in) :: n2    ! real part of refractive index in media "2"
   real(8), intent(out) :: teta2    ! [rad] angle of light propagation in the second media (according to Snell's law)
   real(8), intent(out) :: r12, t12 ! first ray reflection and transmission
   real(8) :: k1x, k2x, k1xk2x
   ! get the angle of light propagation in the second media:
   teta2 = Snell_sin_teta2(teta1, n1, n2)    ! below
   if (teta2 < -1.d-9) then ! total reflection
      r12 = 1.0d0
      t12 = 0.0d0
   else ! partial transmission
      ! Get kix coefficients accroding to Eq.(4.1-4) [1]
      k1x = get_kix(w, n1, teta1)    ! below
      k2x = get_kix(w, n2, teta2)    ! below
      ! Get r12 and t12 coefficients (first ray reflection and transmission):
      k1xk2x = 1.0d0/(k1x + k2x)
      r12 = (k1x - k2x) * k1xk2x     ! Eq.(4.1-14)
      t12 = 2.0d0 * k1x * k1xk2x    ! Eq.(4.1-16)
   endif
end subroutine get_r12_t12

pure function get_kix(w, ni, teta) result(kix)  ! for real index of refraction
   real(8) kix
   real(8), intent(in) :: w     ! [s] frequency of incident light
   real(8), intent(in) :: ni    ! real part of refractive index in media "i"
   real(8), intent(in) :: teta  ! angle of light propagation with respect to normal
   kix = w/g_cvel * ni * cos(teta)  ! Eq.(4.1-4), P. 86 [1]
end function get_kix

pure function get_optical_Phi(n2, ll, dd, teta2) result(Phi)
   real(8) Phi
   real(8), intent(in) :: n2, ll, dd, teta2 ! index of refraction, wavelength [nm], layer thickness [nm], angle of light propagation [rad]
   Phi = 2.0d0 * g_Pi / ll * dd * n2 * cos(teta2)
end function get_optical_Phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OLD:
! General optical coefficients subroutines:
subroutine get_RTA_from_CDF(ReEps, ImEps, ll, dd, teta0, drude_ray, R, T, A)
   real(8), intent(in) :: ReEps, ImEps	! Re(e) and Im(e) of dielectric function
   real(8), intent(in) :: ll  ! [nm] probe-pulse wavelength
   real(8), intent(in) :: dd  ! [nm] experimental layer thickness
   real(8), intent(in) :: teta0 ! angle of the probe pulse [radians]
   integer, intent(in) :: drude_ray	! over how many rays to sum up: (1)=1 or (>1)=all
   real(8), intent(out) :: R, T, A	! reflection, transmittion and absorption coefficients
   complex :: cxn1, cxn2, ckax, ckaz
   complex cs(5), cr12, ct12, cr23, ct23, cxki, cd, cref, ctran
   real(8) :: A2, B2, n, k, alpha, temp
   real(8) :: sai, xkaz, xkax, nm1, np1, exp2ak, dak, dan

   call get_n_k(ReEps, ImEps, n, k)	! transfer Re(e) and Im(e) into n and k
   alpha = 2.0d0*g_Pi*dd/ll
   nm1 = n - 1.0d0
   np1 = n + 1.0d0
   dak = 2.0d0*alpha*k
   if (dak .GT. 300) then
      exp2ak = 1.0d200
   else
      exp2ak = exp(dak)
   endif

   dan = 2.0d0*alpha*n
!    A2 = ((nm1*(1.0d0-exp2ak*cos(dan)-k*exp2ak*sin(dan)))**2+(k*(1.0d0-exp2ak)*cos(dan)+nm1*sin(dan)*exp2ak)**2)/(nm1**2+k**2)
!    B2 = (((np1*np1+k*k)**2-exp2ak*(((n*n+k*k-1.0d0)**2-4.0d0*k*k)*cos(dan)+4.0d0*k*((n*n-1.0d0)+k*k)*sin(dan)))**2 + (exp2ak*(cos(dan)*(-4.0d0*k*(n*n-1.0d0+k*k))+((n*n+k*k-1.0d0)**2-4.0d0*k*k)*sin(dan)))**2 )/(np1*np1+k*k)**4

   cxn1 = cmplx(n,-k)
   sai   = 1.0d0	! air
   cxn2 = cmplx(sai,0.0d0) ! air
   ! Eq (4.4.3) from "multilayer book", p.97:
   temp = 2.0d0*g_Pi*sai/ll
   xkaz = temp*sin(teta0)
   xkax = temp*cos(teta0)
   ckax = cmplx(xkax,0.0d0)
   ckaz = cmplx(xkaz,0.0d0)
   cxki  = 2.0d0*g_Pi*cxn1/ll
   cs(3) = sqrt(cxki*cxki-ckaz*ckaz)

   cd    = cmplx(dd,0.0d0)
   cs(2) = cphase(cs(3),cd) ! function below
   ckaz  = cmplx(xkaz,0.0d0)
   cxki  = 2.0d0*g_Pi*cxn2/ll
   cs(5) = sqrt(cxki*cxki-ckaz*ckaz)

   cr12 = crefl(ckax,cs(3))
   ct12 = ctrans(ckax,cs(3))

   cr23 = crefl(cs(3),cs(5))
   ct23 = ctrans(cs(3),cs(5))
   if (drude_ray .GT. 1) then  ! total summed up:
      cref = (cr12+cr23*cs(2)*cs(2))/(1.0d0+cr12*cr23*cs(2)*cs(2))
      ctran = (ct12*ct23*cs(2))/(1.0d0+cr12*cr23*cs(2)*cs(2))
      R = cabs(cref)*cabs(cref)
      T = cabs(ctran)*cabs(ctran)*cs(5)/ckax
      A = 1.0d0 - R - T
   else ! first ray propagation:
      R = cabs(cr12)*cabs(cr12)
      T = cabs(ct12*ct23*cs(2))
      T = T*T*cs(5)/ckax
      A = 1.0d0 - R - T
   endif
end subroutine get_RTA_from_CDF


!c creflection coefficient between two layers as function of two wave vectors:
function crefl(cxki,cxkj)
   complex :: crefl, cxki, cxkj
   crefl = (cxki-cxkj)/(cxki+cxkj)
end function crefl

!c transmission coefficient between two layers as function of two wave vectors:
function ctrans(cxki,cxkj)
   complex :: ctrans, cxki, cxkj
   ctrans = 2.0d0*cxki/(cxki+cxkj)
end function ctrans

!c phase factor within the layer i:
function cphase(cxkx,cd)
   complex :: cphase, cxkx, cd
   complex :: cc, carg
   cc = (0.0d0,1.0d0)
   carg = -cc*cxkx*cd
   cphase=exp(carg)
end function cphase

END PROGRAM Dielectric_function_analysis
