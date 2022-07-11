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
! This module contains subroutines to deal with ZBL potential
! https://en.wikipedia.org/wiki/Stopping_power_(particle_radiation)#Repulsive_interatomic_potentials

MODULE ZBL_potential
use Universal_constants

implicit none

! Modular parameters:
real(8), parameter :: m_a_u = 0.8854d0
real(8), parameter :: m_phi1 = 0.1818d0
real(8), parameter :: m_phi2 = 0.5099d0
real(8), parameter :: m_phi3 = 0.2802d0
real(8), parameter :: m_phi4 = 0.02817d0
real(8), parameter :: m_exp1 = -3.2d0
real(8), parameter :: m_exp2 = -0.9423d0
real(8), parameter :: m_exp3 = -0.4028d0
real(8), parameter :: m_exp4 = -0.2016d0
real(8), parameter :: m_k = 1.0d0/(4.0d0 * g_Pi * g_e0)


 contains


pure function ZBL_pot(Z1, Z2, r) result(V_ZBL)
   real(8) V_ZBL    ! ZBL potential [eV]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) phi
   phi = ZBL_phi(Z1, Z2, r) ! see below
   V_ZBL = m_k * Z1 * Z2 * g_e/(1.0d-10*r) * phi    ! [eV] ZBL potential
end function ZBL_pot


pure function d_ZBL_pot(Z1, Z2, r) result(d_V_ZBL)
   real(8) d_V_ZBL  ! derivative of ZBL potential [eV/A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) phi, d_phi, V_ZBL_part
   
   V_ZBL_part = m_k * Z1 * Z2 * g_e/(1.0d-10*r) ! part without phi
   phi = ZBL_phi(Z1, Z2, r) ! see below
   d_phi = d_ZBL_phi(Z1, Z2, r) ! see below
   d_V_ZBL = V_ZBL_part* (d_phi - phi/r)  ! [eV/A] derivative of ZBL potential
end function d_ZBL_pot


pure function d_ZBL_phi(Z1, Z2, r) result(phi)
   real(8) phi
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) :: a, x
   a = ZBL_a(Z1, Z2)    ! see below
   x = r / a    ! normalized distance
   phi = (m_phi1 * m_exp1 * exp(m_exp1 * x) + &
          m_phi2 * m_exp2 * exp(m_exp2 * x) + &
          m_phi3 * m_exp3 * exp(m_exp3 * x) + &
          m_phi4 * m_exp4 * exp(m_exp4 * x)) / a
end function d_ZBL_phi


pure function ZBL_phi(Z1, Z2, r) result(phi)
   real(8) phi
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) :: a, x
   a = ZBL_a(Z1, Z2)    ! see below
   x = r / a    ! normalized distance
   phi = m_phi1 * exp(m_exp1 * x) + m_phi2 * exp(m_exp2 * x) + &
         m_phi3 * exp(m_exp3 * x) + m_phi4 * exp(m_exp4 * x)
end function ZBL_phi


pure function ZBL_a(Z1, Z2) result (a)
   real(8) a    ! [A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   a = m_a_u * g_a0/(Z1**0.23d0 + Z2**0.23d0)
end function ZBL_a




END MODULE ZBL_potential
