! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2021 Nikita Medvedev
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
! This module contains subroutines to deal with Gaussian functions in Spherical coordinates,
! transformations are from:
! https://theochem.github.io/horton/2.0.1/tech_ref_gaussian_basis.html

module BS_Spherical_Gaussians

use Universal_constants
use Objects, only: Basis_set
!use Algebra_tools

implicit none

real(8) :: m_sqrt2, m_sqrt3, m_sqrt5, m_sqrt6, m_sqrt7, m_sqrt10, m_sqrt14, m_sqrt21, m_sqrt30, m_sqrt35, m_sqrt70, m_sqrt105
real(8) :: m_sqrt_1_Pi, m_sqrt_3_4Pi, m_sqrt_15_Pi


parameter(m_sqrt2 = sqrt(2.0d0) )
parameter(m_sqrt3 = sqrt(3.0d0) )
parameter(m_sqrt5 = sqrt(5.0d0) )
parameter(m_sqrt6 = sqrt(6.0d0) )
parameter(m_sqrt7 = sqrt(7.0d0) )
parameter(m_sqrt10 = sqrt(10.0d0) )
parameter(m_sqrt14 = sqrt(14.0d0) )
parameter(m_sqrt21 = sqrt(21.0d0) )
parameter(m_sqrt30 = sqrt(30.0d0) )
parameter(m_sqrt35 = sqrt(35.0d0) )
parameter(m_sqrt70 = sqrt(70.0d0) )
parameter(m_sqrt105 = sqrt(105.0d0) )
parameter(m_sqrt_1_Pi = sqrt(1.0d0/g_Pi) )
parameter(m_sqrt_3_4Pi = sqrt(3.0d0/(4.0d0*g_Pi)) )
parameter(m_sqrt_15_Pi = sqrt(15.0d0/g_Pi) )

 
 contains
 
!=======================================================
! Real spherical harmonics: s, p, d:
! https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics

pure function Y00() result(Y_out) ! l=0, Y 0 0
   real(8) Y_out
   Y_out = 0.5d0 * m_sqrt_1_Pi
end function Y00

pure function Y1(x, r) result(Y_out) ! l=1
! identical expressions for: 
! Y 1 -1, py
! Y 1 0, pz
! Y 1 1, px
   real(8) Y_out
   real(8), intent(in) :: x, r
   Y_out = m_sqrt_3_4Pi * x/r
end function Y1

pure function Y2_xy(x, y, r) result(Y_out) ! l=2
! identical expressions for: 
! Y 2 -2, dxy
! Y 2 -1, dyz
! Y 2 1, dxz
   real(8) Y_out
   real(8), intent(in) :: x, y, r
   Y_out = 0.5d0 * m_sqrt_15_Pi * x*y/(r*r)
end function Y2_xy

pure function Y2_z2(x, y, z, r) result(Y_out) ! l=2
! Y 2 0, dz^2
   real(8) Y_out
   real(8), intent(in) :: x, y, z, r
   Y_out = 0.25d0 * m_sqrt5 * m_sqrt_1_Pi * (-x*x - y*y + 2.0d0*z*z)/(r*r)
end function Y2_z2


pure function Y2_x2_y2(x, y, r) result(Y_out) ! l=2
! Y 2 2, dx^2-y^2
   real(8) Y_out
   real(8), intent(in) :: x, y, r
   Y_out = 0.25d0 * m_sqrt_15_Pi * (x*x - y*y)/(r*r)
end function Y2_x2_y2

!======================================================= 
pure function find_TM_size(Cart_Basis) result(siz) ! find transformation matrix size
   integer :: siz
   type(Basis_set), dimension(:), intent(in) :: Cart_Basis
   integer :: tot_siz, i
   tot_siz = 0
   do i = 1, size(Cart_Basis)
      tot_siz = tot_siz + find_AM_size(Cart_Basis(i)%AM(:))   ! below
   enddo
   siz = tot_siz
end function find_TM_size
 

pure function find_AM_size(AM) result(siz) ! find size of the angular momentum matrix
   integer :: siz
   integer, dimension(:), intent(in) :: AM  ! angular momenta
   integer :: i, cur, T_AM
   T_AM = SUM(AM(:))    ! total angular momentum
   selectcase(AM(i))
   case (0)  ! s
      cur = 1
   case (1)  ! p
      cur = 3
   case (2)  ! d
      cur = 5
   case (3)  ! f
      cur = 7
   case (4)  ! g
      cur = 9
   case default
      cur = 0
   end select
   siz = cur
end function find_AM_size


pure function find_AM_size_cartesian(AM) result(siz) ! find size of the angular momentum matrix
   integer :: siz
   integer, dimension(:), intent(in) :: AM  ! angular momenta
   integer :: i, cur, T_AM
   T_AM = SUM(AM(:))    ! total angular momentum
   selectcase(AM(i))
   case (0)  ! s
      cur = 1
   case (1)  ! p
      cur = 3
   case (2)  ! d
      cur = 6
   case (3)  ! f
      cur = 10
   case (4)  ! g
      cur = 15
   case default
      cur = 0
   end select
   siz = cur
end function find_AM_size_cartesian


 
subroutine transformatrix_Sph_C(AM, C_transform, i_start, j_start)
   integer, dimension(:), intent(in) :: AM  ! angular momenta
   integer, intent(in) :: i_start, j_start  ! indices where to start placing the new elements
   real(8), dimension(:,:), intent(inout) :: C_transform    ! transformation matrix from cartesian to spherical gaussians
   integer :: i, j, siz_i, siz_j, T_AM
   
   ! Find the sizes:
   siz_i = find_AM_size(AM)     ! above
   siz_j = find_AM_size_cartesian(AM)   ! above
   T_AM = SUM(AM(:))    ! total angular momentum
   
   ! Define the transformation matrix from crtesian to spherical gaussians
   ! https://theochem.github.io/horton/2.0.1/tech_ref_gaussian_basis.html
   selectcase(AM(i))    ! set only non-zero elements
   case (0)  ! s
      C_transform(i_start+1,j_start+1) = 1.0d0
   case (1)  ! p
      C_transform(i_start+1,j_start+3) = 1.0d0
      C_transform(i_start+2,j_start+1) = 1.0d0
      C_transform(i_start+3,j_start+2) = 1.0d0
   case (2)  ! d
      C_transform(i_start+1,j_start+1) = -0.5d0
      C_transform(i_start+1,j_start+4) = -0.5d0
      C_transform(i_start+1,j_start+6) = 1.0d0
      C_transform(i_start+2,j_start+3) = 1.0d0
      C_transform(i_start+3,j_start+5) = 1.0d0
      C_transform(i_start+4,j_start+1) = 0.5d0*m_sqrt3
      C_transform(i_start+4,j_start+4) = -0.5d0*m_sqrt3
      C_transform(i_start+5,j_start+2) = 1.0d0
   case (3)  ! f
      C_transform(i_start+1,j_start+3) = -3.0d0/10.0d0*m_sqrt5
      C_transform(i_start+1,j_start+8) = C_transform(i_start+1,j_start+3)
      C_transform(i_start+1,j_start+10) = 1.0d0
      C_transform(i_start+2,j_start+1) = -0.25d0*m_sqrt6
      C_transform(i_start+2,j_start+4) = -1.0d0/20.0d0*m_sqrt30
      C_transform(i_start+2,j_start+6) = 0.2d0*m_sqrt30
      C_transform(i_start+3,j_start+2) = -1.0d0/20.0d0*m_sqrt30
      C_transform(i_start+3,j_start+7) = -0.25d0*m_sqrt6
      C_transform(i_start+3,j_start+9) = 0.5d0*m_sqrt30
      C_transform(i_start+4,j_start+3) = 0.5d0*m_sqrt3
      C_transform(i_start+4,j_start+8) = -0.5d0*m_sqrt3
      C_transform(i_start+5,j_start+5) = 1.0d0
      C_transform(i_start+6,j_start+1) = 0.25d0*m_sqrt10
      C_transform(i_start+6,j_start+4) = -0.75d0*m_sqrt2
      C_transform(i_start+7,j_start+2) = 0.75d0*m_sqrt2
      C_transform(i_start+7,j_start+7) = -0.25d0*m_sqrt10
   case (4)  ! g
      ! Not ready yet!
   end select
end subroutine transformatrix_Sph_C


 
end module BS_Spherical_Gaussians
