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
! This module contains Koster-Slater integrals for all shells:

MODULE TB_Koster_Slater

use Algebra_tools, only: Kronecker_delta

implicit none

! Modular parameters:
real(8) :: m_sqrt3, m_sqrt3_half

parameter (m_sqrt3 = sqrt(3.0d0))
parameter (m_sqrt3_half = m_sqrt3*0.5d0)

 contains

!---------------------------------------------------------------------
! Derivatives of r_{i,j} and d_{i,j} by relative coordinates s_{k}

! Derivative of the relative distance between the two atoms by a coordinate of a third atom:
! d r_{ij} / d s_{k,a}
pure function drij_dska(i, j, k, xij, yij, zij, rij, supce, alpha, with_delta)
   real(8) :: drij_dska	! derivative : d r_{ij} / d s_{k,a}
   integer, intent (in) :: i, j, k	! atomic indices
   real(8), intent (in) :: xij, yij, zij, rij	! interatomic distance and its projections
   real(8), dimension(3,3), intent(in) :: supce	! matrix of supercell vectors
   integer, intent(in) :: alpha	! along which coordinate we are calculating the derivative: x=1, y=2, z=3
   logical, intent(in), optional :: with_delta	! calculate the part independent of delta, meaning dependent of k-atom
   real(8) :: dik, djk, delt
   delt = 1.0d0
   if (present(with_delta)) then
      if (with_delta) then
         dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
         djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
         delt = (dik - djk)
      endif
   endif
   drij_dska = delt/rij*( xij*supce(alpha,1) + yij*supce(alpha,2) + zij*supce(alpha,3) )
!    drij_dska = delt/rij*( xij*supce(1,alpha) + yij*supce(2,alpha) + zij*supce(3,alpha) )
end function drij_dska



! Derivative of the relative distance between the two atoms by a coordinate of a third atom:
! d d_{ij} / d s_{k,a} , if the drij_dsk is precalculated:
pure function ddija_dskb_kd(i, j, k, xij, yij, zij, rij, supce, alpha, beta, drij_dsk) result(ddija_dskb)
   real(8) :: ddija_dskb	! derivative : d d_{ij,a} / d s_{k,b}
   integer, intent (in) :: i, j, k	! atomic indices
   real(8), intent (in) :: xij, yij, zij, rij	! interatomic distance and its projections
   real(8), dimension(3,3), intent(in) :: supce	! matrix of supercell vectors
   integer, intent(in) :: alpha, beta	! which cosine and along which coordinate we are calculating the derivative: x=1, y=2, z=3
   real(8), intent(in) :: drij_dsk	! directional cosine is known
   real(8) :: dik, djk, delt, rija
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   delt = (dik - djk)
   select case (alpha)	! which cosine: 1=x, 2=y, 3=z
   case (1)	! x
      rija = xij
   case (2)	! y
      rija = yij
   case (3)	! z
      rija = zij
   end select
   ddija_dskb = (delt*supce(beta,alpha) - rija/rij*drij_dsk)/rij
end function ddija_dskb_kd


! Derivative of the relative distance between the two atoms by a coordinate of a third atom:
! d d_{ij} / d s_{k,a}  --  CHECKED, CORRECT
pure function ddija_dskb(i, j, k, xij, yij, zij, rij, supce, alpha, beta, with_delta)
   real(8) :: ddija_dskb	! derivative : d d_{ij,a} / d s_{k,b}
   integer, intent (in) :: i, j, k	! atomic indices
   real(8), intent (in) :: xij, yij, zij, rij	! interatomic distance and its projections
   real(8), dimension(3,3), intent(in) :: supce	! matrix of supercell vectors
   integer, intent(in) :: alpha, beta	! which cosine and along which coordinate we are calculating the derivative: x=1, y=2, z=3
   logical, intent(in), optional :: with_delta	! calculate the part independent of delta, meaning dependent of k-atom
   real(8) :: dik, djk, delt, rija, drij_dsk, rij2
   delt = 1.0d0
   if (present(with_delta)) then
      if (with_delta) then
         dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
         djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
         delt = (dik - djk)
         drij_dsk = drij_dska(i, j, k, xij, yij, zij, rij, supce, beta, with_delta=.true.)
      else
         drij_dsk = drij_dska(i, j, k, xij, yij, zij, rij, supce, beta)
      endif
   else
      drij_dsk = drij_dska(i, j, k, xij, yij, zij, rij, supce, beta)
   endif
   select case (alpha)	! which cosine: 1=x, 2=y, 3=z
   case (1)	! x
      rija = xij
   case (2)	! y
      rija = yij
   case (3)	! z
      rija = zij
   end select
   rij2 = rij*rij
!    ddija_dskb = (delt*supce(alpha,beta) - rija/rij*drij_dsk)/rij
   ddija_dskb = (delt*supce(beta,alpha) - rija/rij*drij_dsk)/rij
end function ddija_dskb


!---------------------------------------------------------------------
! Derivatives of r_{i,j} and d_{i,j} by supercell coordinates h_{a,b}

pure function  drij_dhab(rija, sijb, rij)
   real(8) :: drij_dhab
   real(8), intent(in) :: rija, sijb, rij
   drij_dhab = rija*sijb/rij
end function  drij_dhab


pure function  dda_dhgd(alpha, gama, rda, rijg, sijd, rij)
   real(8) :: dda_dhgd
   integer, intent(in) :: alpha, gama
   real(8), intent(in) :: rda, rijg, sijd, rij
   real(8) :: r2, delt
   r2 = rij*rij
   delt = Kronecker_delta(alpha,gama) 	! module "Algebra_tools"
   dda_dhgd = sijd/(r2*rij)*(r2*delt - rda*rijg)
end function  dda_dhgd


!---------------------------------------------------------------------
! Derivatives of r_{i,j} and d_{i,j} by absolute coordinates r_{k}

! Derivative of the relative distance between the two atoms by a coordinate of a third atom:
! d R_{ij} / d x_k,    or d R_{ij} / d y_k,    or d R_{ij} / d z_k : 
pure function drij_drka(i, j, k, rija, rij)
   real(8) :: drij_drka
   integer, intent (in) :: i, j, k
   real(8), intent (in) :: rija, rij
   real(8) :: dik, djk
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   drij_drka = (dik - djk)*rija/rij
end function drij_drka


! Second derivative of the relative distance between the two atoms by coordinates of a third and fourth atom:
! d2 R_{ij} / (d r_{k,a} d r_{l,b} )
pure function d2rij_drka_drlb(i, j, k, l, alpha, beta, rija, rijb, rij) result(d2r)
   real(8) :: d2r
   integer, intent (in) :: i, j, k, l, alpha, beta
   real(8), intent (in) :: rija, rijb, rij
   real(8) :: dik, djk, dil, djl, dab
   ! Delta symbols:
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   dil = Kronecker_delta(i,l) 	! module "Algebra_tools"
   djl = Kronecker_delta(j,l) 	! module "Algebra_tools"
   dab= Kronecker_delta(alpha, beta)	! module "Algebra_tools"
   d2r = (dik - djk)*(dil - djl)/rij*(dab - rija*rijb/(rij*rij))
end function d2rij_drka_drlb


! Derivative of the direction cosine by a coordinate of a third atom:
! d d_{ij,a} / d r_{k,b}
pure function ddija_drkb(i, j, k, alpha, beta, rija, rijb, rij)
   real(8) :: ddija_drkb
   integer, intent (in) :: i, j, k, alpha, beta
   real(8), intent (in) :: rija, rij, rijb
   real(8) :: dik, djk, dab
   ! Delta-symbols:
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   dab= Kronecker_delta(alpha, beta)	! module "Algebra_tools"
   ddija_drkb = (dik - djk)/rij*(dab - rija*rijb/(rij*rij))
end function ddija_drkb


! Second derivatives of the direction cosine by a coordinate of a third atom:
! d^2 d_{ij,a} / (d r_{k,b} d r_{l,g})1
pure function d2dija_drkb_drlg(i, j, k, l, alpha, beta, gamma, rija, rijb, rijg, rij) result(d2d)
   real(8) :: d2d
   integer, intent (in) :: i, j, k, l, alpha, beta, gamma
   real(8), intent (in) :: rija, rijb, rijg, rij
   real(8) :: dik, djk, dil, djl, dab, dag, dbg, rij2, rij3
   ! Delta-symbols:
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   dil = Kronecker_delta(i,l) 	! module "Algebra_tools"
   djl = Kronecker_delta(j,l) 	! module "Algebra_tools"
   dab= Kronecker_delta(alpha, beta)	! module "Algebra_tools"
   dag= Kronecker_delta(alpha, gamma)	! module "Algebra_tools"
   dbg= Kronecker_delta(beta, gamma)	! module "Algebra_tools"
   rij2 = rij*rij
   rij3 = rij*rij2
   d2d = -(dik - djk)*(dil - djl)/rij3 * ( dab*rijg + dag*rijb + dbg*rija - 3.0d0*rija*rijb*rijg/rij2 )
end function d2dija_drkb_drlg


! Second derivatives of the directional cosine by a coordinate of a third atom:
! d^2 d_{ij,a} / d r^2_{k,b}
pure function d2dija_drkb2(i, j, k, alpha, beta, rija, rijb, rij)
   real(8) :: d2dija_drkb2
   integer, intent (in) :: i, j, k, alpha, beta
   real(8), intent (in) :: rija, rijb, rij
   real(8) :: dik, djk, dab, rij2, rij3
   ! Delta-symbols:
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   dab= Kronecker_delta(alpha, beta)	! module "Algebra_tools"
   rij2 = rij*rij
   rij3 = rij*rij2
   d2dija_drkb2 = -(dik - djk)/rij3 * ( 2.0d0*rijb*(dab - rija*rijb/rij2) + rija*(1.0d0 - rijb*rijb/rij2) )
end function d2dija_drkb2


!---------------------------------------------------------------------
! All Koster-Slater hopping integrals as individual functions:
! 1) 
pure function t_s_s(Vss_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: Vss_sigma
   Ecc = Vss_sigma
end function t_s_s


! 2) 
pure function t_s_px(da, Vsp_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, Vsp_sigma
   Ecc = da*Vsp_sigma
end function t_s_px


! 3) 
pure function t_pa_pa(da, Vpp_sigma, Vpp_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, Vpp_sigma, Vpp_pi
   real(8) :: da2
   da2 = da*da
   Ecc = (Vpp_sigma - Vpp_pi)*da2 + Vpp_pi
end function t_pa_pa


! 4)
pure function t_pa_pb(da, db, Vpp_sigma, Vpp_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, db, Vpp_sigma, Vpp_pi
   real(8) :: dadb
   dadb = da*db
   Ecc = (Vpp_sigma - Vpp_pi)*dadb
end function t_pa_pb


! 5)
pure function t_s_dab(da, db, Vsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, db, Vsd_sigma
   Ecc = m_sqrt3*da*db*Vsd_sigma
end function t_s_dab


! 6)
pure function t_s_dx2_y2(da, db, Vsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, db, Vsd_sigma
   Ecc = m_sqrt3_half*(da*da - db*db)*Vsd_sigma
end function t_s_dx2_y2


! 7)
pure function t_s_dz2_r2(l, m, n, Vsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vsd_sigma
   Ecc = ( n*n - 0.5d0*(l*l + m*m) )*Vsd_sigma
end function t_s_dz2_r2


! 8)
pure function t_px_dxy(l, m, Vpd_sigma, Vpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, Vpd_sigma, Vpd_pi
   real(8) :: l2
   l2 = l*l
   Ecc = (m_sqrt3*l2*Vpd_sigma + (1.0d0 - 2.0d0*l2)*Vpd_pi)*m
end function t_px_dxy


! 9)
pure function t_px_dyz(l, m, n, Vpd_sigma, Vpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vpd_sigma, Vpd_pi
   Ecc = (m_sqrt3*Vpd_sigma - 2.0d0*Vpd_pi)*l*m*n
end function t_px_dyz


! 10)
pure function t_px_dx2_y2(l, m, Vpd_sigma, Vpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, Vpd_sigma, Vpd_pi
   real(8) :: l2_m2
   l2_m2 = l*l - m*m
   Ecc = ( (m_sqrt3_half*Vpd_sigma - Vpd_pi)*l2_m2 + Vpd_pi )*l
end function t_px_dx2_y2


! 11)
pure function t_py_dx2_y2(l, m, Vpd_sigma, Vpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, Vpd_sigma, Vpd_pi
   real(8) :: l2_m2
   l2_m2 = l*l - m*m
   Ecc = ( (m_sqrt3_half*Vpd_sigma - Vpd_pi)*l2_m2 - Vpd_pi )*m
end function t_py_dx2_y2


! 12)
pure function t_pz_dx2_y2(l, m, n, Vpd_sigma, Vpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vpd_sigma, Vpd_pi
   real(8) :: l2_m2
   l2_m2 = l*l - m*m
   Ecc = (m_sqrt3_half*Vpd_sigma - Vpd_pi)*l2_m2*n
end function t_pz_dx2_y2


! 13)
pure function t_pa_d3z2_r2(da, l, m, n, Vpd_sigma, Vpd_pi) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: da, l, m, n, Vpd_sigma, Vpd_pi
   real(8) :: n2
   n2 = n*n
   Ecc = ( (n2 - 0.5d0*(l*l + m*m))*Vpd_sigma - m_sqrt3*n2*Vpd_pi )*da
end function t_pa_d3z2_r2


! 14)
pure function t_pz_d3z2_r2(l, m, n, Vpd_sigma, Vpd_pi) result (Ecc) ! only pz
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vpd_sigma, Vpd_pi
   real(8) :: l2_m2, n2
   l2_m2 = l*l + m*m
   n2 = n*n
   Ecc = ( (n2 - 0.5d0*l2_m2)*Vpd_sigma + m_sqrt3*l2_m2*Vpd_pi )*n
end function t_pz_d3z2_r2


! 15)
pure function t_dab_dab(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: l2, m2, l2m2
   l2 = l*l
   m2 = m*m
   l2m2 = l2*m2
   Ecc = 3.0d0*l2m2*Vdd_sigma + (l2 + m2 - 4.0d0*l2m2)*Vdd_pi + (n*n + l2m2)*Vdd_delta
end function t_dab_dab


! 16)
pure function t_dab_dbg(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: ln, m2
   ln = l*n
   m2 = m*m
   Ecc = ( 3.0d0*m2*Vdd_sigma + (1.0d0 - 4.0d0*m2)*Vdd_pi + (m2 - 1.0d0)*Vdd_delta )*ln
end function t_dab_dbg


! 17)
pure function t_dxy_dx2_y2(l, m, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: lm, l2m2
   lm = l*m
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*Vdd_sigma - 2.0d0*Vdd_pi + 0.5d0*Vdd_delta )*lm*l2m2
end function t_dxy_dx2_y2


! 18)
pure function t_dyz_dx2_y2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: mn, l2m2
   mn = m*n
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*l2m2*Vdd_sigma - (1.0d0+2.0d0*l2m2)*Vdd_pi + (1.0d0+0.5d0*l2m2)*Vdd_delta )*mn
end function t_dyz_dx2_y2


! 19)
pure function t_dxz_dx2_y2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: ln, l2m2
   ln = l*n
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*l2m2*Vdd_sigma + (1.0d0-2.0d0*l2m2)*Vdd_pi - (1.0d0-0.5d0*l2m2)*Vdd_delta )*ln
end function t_dxz_dx2_y2


! 20)
pure function t_dxy_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: lm, l2m2, n2
   lm = l*m
   l2m2 = l*l + m*m
   n2 = n*n
   Ecc = ( m_sqrt3*( (n2 - 0.5d0*l2m2)*Vdd_sigma - 2.0d0*n2*Vdd_pi) + m_sqrt3_half*(1.0d0+n2)*Vdd_delta )*lm
end function t_dxy_d3z2_r2


! 21)
pure function t_dyz_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: mn, l2m2, n2
   mn = m*n
   l2m2 = l*l + m*m
   n2 = n*n
   Ecc = ( (n2 - 0.5d0*l2m2)*Vdd_sigma + (l2m2 - n2)*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*m_sqrt3*mn
end function t_dyz_d3z2_r2


! 22)
pure function t_dxz_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: ln, l2m2, n2
   ln = l*n
   l2m2 = l*l + m*m
   n2 = n*n
   Ecc = ( (n2 - 0.5d0*l2m2)*Vdd_sigma + (l2m2 - n2)*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*m_sqrt3*ln
end function t_dxz_d3z2_r2


! 23)
pure function t_dx2_y2_dx2_y2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: l2m2, l2_m2, l2_m2_2, n2
   l2m2 = l*l + m*m
   l2_m2 = l*l - m*m
   l2_m2_2 = l2_m2*l2_m2
   n2 = n*n
   Ecc = 0.75d0*l2_m2_2*Vdd_sigma + (l2m2 - l2_m2_2)*Vdd_pi + (n2+0.25d0*l2_m2_2)*Vdd_delta
end function t_dx2_y2_dx2_y2


! 24)
pure function t_dx2_y2_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: l2m2, l2_m2, n2
   l2m2 = l*l + m*m
   l2_m2 = l*l - m*m
   n2 = n*n
   Ecc = ( 0.5d0*(n2 - 0.5d0*l2m2)*Vdd_sigma - n2*Vdd_pi + 0.25d0*(1.0d0+n2)*Vdd_delta )*m_sqrt3*l2_m2
end function t_dx2_y2_d3z2_r2


! 25)
pure function t_d3z2_r2_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: l2m2, n2, temp
   l2m2 = l*l + m*m
   n2 = n*n
   temp = n2 - 0.5d0*l2m2
   Ecc = temp*temp*Vdd_sigma + 3.0d0*n2*l2m2*Vdd_pi + 0.75d0*l2m2*l2m2*Vdd_delta
end function t_d3z2_r2_d3z2_r2


!------------------------------------------------------------------------
! Derivatives of the Koster-Slater angular functions.
! This functions do NOT include the derivatives of dr/ds, or dr/dh, or dr/dr_k,
! only the derivatives of the other parts:
! 1) 
pure function d_t_s_s(dVss_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: dVss_sigma
   Ecc = dVss_sigma
end function d_t_s_s


! 2) 
pure function d_t_s_px(da, dda, Vsp_sigma, dVsp_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, Vsp_sigma, dVsp_sigma
   Ecc = dda*Vsp_sigma + da*dVsp_sigma
end function d_t_s_px


! 3) 
pure function d_t_pa_pa(da, dda, Vpp_sigma, dVpp_sigma, Vpp_pi, dVpp_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, Vpp_sigma, dVpp_sigma, Vpp_pi, dVpp_pi
   real(8) :: da2
   da2 = da*da
   Ecc = (dVpp_sigma - dVpp_pi)*da2 + (Vpp_sigma - Vpp_pi)*2.0d0*da*dda + dVpp_pi
end function d_t_pa_pa


! 4)
pure function d_t_pa_pb(da, dda, db, ddb, Vpp_sigma, dVpp_sigma, Vpp_pi, dVpp_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, db, ddb, Vpp_sigma, dVpp_sigma, Vpp_pi, dVpp_pi
   Ecc = (dVpp_sigma - dVpp_pi)*(da*db) + (Vpp_sigma - Vpp_pi)*(dda*db+da*ddb)
end function d_t_pa_pb


! 5)
pure function d_t_s_dab(da, dda, db, ddb, Vsd_sigma, dVsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, db, ddb, Vsd_sigma, dVsd_sigma
   Ecc = m_sqrt3*( (dda*db + da*ddb)*Vsd_sigma + da*db*dVsd_sigma)
end function d_t_s_dab


! 6)
pure function d_t_s_dx2_y2(da, dda, db, ddb, Vsd_sigma, dVsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, db, ddb, Vsd_sigma, dVsd_sigma
   Ecc = 2.0d0*(dda*da - ddb*db)*Vsd_sigma + (da*da - db*db)*dVsd_sigma
   Ecc = m_sqrt3_half*Ecc
end function d_t_s_dx2_y2


! 7)
pure function d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsd_sigma, dVsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vsd_sigma, dVsd_sigma
   Ecc = ( 2.0d0*dn*n - (dl*l + dm*m) )*Vsd_sigma + ( n*n - 0.5d0*(l*l + m*m) )*dVsd_sigma
end function d_t_s_dz2_r2


! 8)
pure function d_t_px_dxy(l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2, temp, dll
   l2 = l*l
   dll = 2.0d0*dl*l
   temp = (1.0d0 - 2.0d0*l2)
   Ecc = (m_sqrt3*l2*Vpd_sigma + temp*Vpd_pi)*dm
   Ecc = Ecc  + (m_sqrt3*(dll*Vpd_sigma + l2*dVpd_sigma) - 2.0d0*dll*Vpd_pi + temp*dVpd_pi)*m
end function d_t_px_dxy


! 9)
pure function d_t_px_dyz(l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: mn
   mn = m*n
   Ecc = (m_sqrt3*dVpd_sigma - 2.0d0*dVpd_pi)*l*mn
   Ecc = Ecc + (m_sqrt3*Vpd_sigma - 2.0d0*Vpd_pi)*(dl*mn + l*dm*n + l*m*dn)
end function d_t_px_dyz


! 10)
pure function d_t_px_dx2_y2(l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2, temp
   l2_m2 = l*l - m*m
   dl2_m2 = 2.0d0*( dl*l - dm*m )
   temp = (m_sqrt3_half*Vpd_sigma - Vpd_pi)
   Ecc = ( temp*l2_m2 + Vpd_pi )*dl
   Ecc = Ecc + ( (m_sqrt3_half*dVpd_sigma - dVpd_pi)*l2_m2 + temp*dl2_m2 + dVpd_pi )*l
end function d_t_px_dx2_y2


! 11)
pure function d_t_py_dx2_y2(l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2, temp
   l2_m2 = l*l - m*m
   dl2_m2 = 2.0d0*( dl*l - dm*m )
   temp = (m_sqrt3_half*Vpd_sigma - Vpd_pi)
   Ecc = ( temp*l2_m2 - Vpd_pi )*dm
   Ecc = Ecc + ( (m_sqrt3_half*dVpd_sigma - dVpd_pi)*l2_m2 + temp*dl2_m2 - dVpd_pi )*m
end function d_t_py_dx2_y2


! 12)
pure function d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2
   l2_m2 = l*l - m*m
   dl2_m2 = 2.0d0*(dl*l - dm*m)
   Ecc = (m_sqrt3_half*dVpd_sigma - dVpd_pi)*l2_m2*n
   Ecc = Ecc + (m_sqrt3_half*Vpd_sigma - Vpd_pi)*(dl2_m2*n + l2_m2*dn)
end function d_t_pz_dx2_y2


! 13)
pure function d_t_pa_d3z2_r2(da, dda, l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: n2, temp, dn2
   n2 = n*n
   dn2 = 2.0d0*dn*n
   temp = (n2 - 0.5d0*(l*l + m*m))
   Ecc = ( temp*Vpd_sigma - m_sqrt3*n2*Vpd_pi )*dda
   Ecc = Ecc + ( (dn2 - (dl*l + dm*m))*Vpd_sigma + temp*dVpd_sigma - m_sqrt3*(dn2*Vpd_pi + n2*dVpd_pi) )*da
end function d_t_pa_d3z2_r2


! 14)
pure function d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc) ! only pz
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2, n2, temp
   l2_m2 = l*l + m*m
   dl2_m2 = 2.0d0*(dl*l + dm*m)
   n2 = n*n
   temp = (n2 - 0.5d0*l2_m2)
   Ecc = ( temp*Vpd_sigma + m_sqrt3*l2_m2*Vpd_pi )*dn
   Ecc = Ecc  + ( (2.0d0*dn*n - 0.5d0*dl2_m2)*Vpd_sigma + temp*dVpd_sigma + m_sqrt3*(dl2_m2*Vpd_pi + l2_m2*dVpd_pi) )*n
end function d_t_pz_d3z2_r2


! 15)
pure function d_t_dab_dab(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2, m2, l2m2, dl2m2
   l2 = l*l
   m2 = m*m
   l2m2 = l2*m2
   dl2m2 = 2.0d0*(dl*l*m2 + l2*m*dm)
   Ecc = 3.0d0*(dl2m2*Vdd_sigma + l2m2*dVdd_sigma)
   Ecc = Ecc + 2.0d0*(l*dl + m*dm - 2.0d0*dl2m2)*Vdd_pi + (l2 + m2 - 4.0d0*l2m2)*dVdd_pi
   Ecc = Ecc + (2.0d0*n*dn + dl2m2)*Vdd_delta +  (n*n + l2m2)*dVdd_delta
end function d_t_dab_dab


! 16)
pure function d_t_dab_dbg(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: ln, m2, dm2, temp, temp1
   ln = l*n
   m2 = m*m
   dm2 = 2.0d0*m*dm
   temp = (1.0d0 - 4.0d0*m2)
   temp1 = (m2 - 1.0d0)
   Ecc = 3.0d0*(dm2*Vdd_sigma + m2*dVdd_sigma)
   Ecc = Ecc + temp*dVdd_pi - 4.0d0*dm2*Vdd_pi 
   Ecc = Ecc +  temp1*dVdd_delta + dm2*Vdd_delta
   Ecc = Ecc*ln
   Ecc = Ecc + ( 3.0d0*m2*Vdd_sigma + temp*Vdd_pi + temp1*Vdd_delta )*(dl*n+l*dn)
end function d_t_dab_dbg


! 17)
pure function d_t_dxy_dx2_y2(l, dl, m, dm, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: lm, l2m2
   lm = l*m
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*dVdd_sigma - 2.0d0*dVdd_pi + 0.5d0*dVdd_delta )*lm*l2m2
   Ecc = Ecc + ( 1.50d0*Vdd_sigma - 2.0d0*Vdd_pi + 0.5d0*Vdd_delta )*((dl*m+l*dm)*l2m2 + 2.0d0*lm*(dl*l-dm*m))
end function d_t_dxy_dx2_y2


! 18)
pure function d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: mn, l2m2, dl2m2, temp, temp1
   mn = m*n
   l2m2 = l*l - m*m
   dl2m2 = 2.0d0*(dl*l - dm*m)
   temp = (1.0d0+2.0d0*l2m2)
   temp1 = (1.0d0+0.5d0*l2m2)
   Ecc = ( 1.50d0*(dl2m2*Vdd_sigma+l2m2*dVdd_sigma) - (2.0d0*dl2m2*Vdd_pi + temp*dVdd_pi) + temp1*dVdd_delta + 0.5d0*dl2m2*Vdd_delta )*mn
   Ecc = Ecc + ( 1.50d0*l2m2*Vdd_sigma - temp*Vdd_pi + temp1*Vdd_delta )*(dm*n+m*dn)
end function d_t_dyz_dx2_y2


! 19)
pure function d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: ln, l2m2, dl2m2, temp, temp1
   ln = l*n
   l2m2 = l*l - n*n
   dl2m2 = 2.0d0*(dl*l - dm*m)
   temp = (1.0d0-2.0d0*l2m2)
   temp1 = (1.0d0-0.5d0*l2m2)
   Ecc = ( 1.50d0*(dl2m2*Vdd_sigma+l2m2*dVdd_sigma) + (-2.0d0*dl2m2*Vdd_pi + temp*dVdd_pi) - (temp1*dVdd_delta - 0.5d0*dl2m2*Vdd_delta) )*ln
   Ecc = Ecc + ( 1.50d0*l2m2*Vdd_sigma + temp*Vdd_pi + temp1*Vdd_delta )*(dl*n+l*dn)
end function d_t_dxz_dx2_y2


! 20)
pure function d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: lm, l2m2, n2, dl2m2, temp, temp1, dn2
   lm = l*m
   l2m2 = l*l + m*m
   dl2m2 = 2.0d0*(dl*l + dm*m)
   n2 = n*n
   dn2 = 2.0d0*n*dn
   temp = (n2 - 0.5d0*l2m2)
   temp1 = (1.0d0+n2)
   Ecc = ( m_sqrt3*( ((dn2 - 0.5d0*dl2m2)*Vdd_sigma + temp*dVdd_sigma) - 2.0d0*(dn2*Vdd_pi+n2*dVdd_pi) ) + m_sqrt3_half*(dn2*Vdd_delta + temp1*dVdd_delta) )*lm
   Ecc = Ecc + ( m_sqrt3*(temp*Vdd_sigma - 2.0d0*n2*Vdd_pi) + m_sqrt3_half*temp1*Vdd_delta )*(dl*m + l*dm)
end function d_t_dxy_d3z2_r2


! 21)
pure function d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: mn, l2m2, n2, temp, temp1, dn2, dl2m2
   mn = m*n
   l2m2 = l*l + m*m
   dl2m2 = 2.0d0*(l*dl + m*dm)
   n2 = n*n
   dn2 = 2.0d0*n*dn
   temp = (n2 - 0.5d0*l2m2)
   temp1 = (l2m2 - n2)
   Ecc = ( temp*Vdd_sigma + temp1*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*(dm*n+m*dn)
   Ecc = Ecc + ( (dn2 - 0.5d0*dl2m2)*Vdd_sigma + temp*dVdd_sigma + (dl2m2 - dn2)*Vdd_pi + temp1*dVdd_pi - 0.5d0*(dl2m2*Vdd_delta + l2m2*dVdd_delta) )*mn
   Ecc = Ecc*m_sqrt3
end function d_t_dyz_d3z2_r2


! 22)
pure function d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: ln, l2m2, n2, temp, temp1, dn2, dl2m2
   ln = l*n
   l2m2 = l*l + m*m
   dl2m2 = 2.0d0*(l*dl + m*dm)
   n2 = n*n
   dn2 = 2.0d0*n*dn
   temp = (n2 - 0.5d0*l2m2)
   temp1 = (l2m2 - n2)
   Ecc = ( temp*Vdd_sigma + temp1*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*(dl*n+l*dn)
   Ecc = Ecc + ( (dn2 - 0.5d0*dl2m2)*Vdd_sigma + temp*dVdd_sigma + (dl2m2 - dn2)*Vdd_pi + temp1*dVdd_pi - 0.5d0*(dl2m2*Vdd_delta + l2m2*dVdd_delta) )*ln
   Ecc = Ecc*m_sqrt3
end function d_t_dxz_d3z2_r2


! 23)
pure function d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2m2, l2_m2, l2_m2_2, n2, dl2m2, dl2_m2_2, l2, m2, ldl, mdm, dn2
   l2 = l*l
   m2 = m*m
   l2m2 = l2 + m2
   ldl = l*dl
   mdm = m*dm
   dl2m2 = 2.0d0*(ldl + mdm)
   l2_m2 = l2 - m2
   l2_m2_2 = l2_m2*l2_m2
   dl2_m2_2 = 4.0d0*l2_m2*(ldl - mdm)
   n2 = n*n
   dn2 = 2.0d0*n*dn
   Ecc = 0.75d0*(dl2_m2_2*Vdd_sigma + l2_m2_2*dVdd_sigma)
   Ecc = Ecc + (dl2m2 - dl2_m2_2)*Vdd_pi + (l2m2 - l2_m2_2)*dVdd_pi  
   Ecc = Ecc + (dn2+0.25d0*dl2_m2_2)*Vdd_delta + (n2+0.25d0*l2_m2_2)*dVdd_delta
end function d_t_dx2_y2_dx2_y2


! 24)
pure function d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2m2, l2_m2, n2, dn2, l2, m2, dl2, dm2, dl2m2, dl2_m2, temp, temp1
   l2 = l*l
   dl2 = 2.0d0*l*dl
   m2 = m*m
   dm2 = 2.0d0*m*dm
   l2m2 = l2 + m2
   dl2m2 = dl2 + dm2
   l2_m2 = l2 - m2
   dl2_m2 = dl2 - dm2
   n2 = n*n
   dn2 = 2.0d0*n*dn
   temp = (n2 - 0.5d0*l2m2)
   temp1 = (1.0d0+n2)
   Ecc = ( 0.5d0*temp*Vdd_sigma - n2*Vdd_pi + 0.25d0*temp1*Vdd_delta )*dl2_m2
   Ecc = Ecc + ( 0.5d0*((dn2 - 0.5d0*dl2m2)*Vdd_sigma+temp*dVdd_sigma) - (dn2*Vdd_pi + n2*dVdd_pi) + 0.25d0*(dn2*Vdd_delta + temp1*dVdd_delta) )*l2_m2
   Ecc = Ecc*m_sqrt3
end function d_t_dx2_y2_d3z2_r2


! 25)
pure function d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2m2, dl2m2, n2, temp, dtemp, dn2
   l2m2 = l*l + m*m
   dl2m2 = 2.0d0*( dl*l + dm*m )
   n2 = n*n
   dn2 = 2.0d0*n*dn
   temp = n2 - 0.5d0*l2m2
   dtemp = dn2 - 0.5d0*dl2m2
   Ecc = (2.0d0*dtemp*Vdd_sigma + temp*dVdd_sigma)*temp
   Ecc = Ecc + 3.0d0*((dn2*l2m2 + n2*dl2m2)*Vdd_pi + n2*l2m2*dVdd_pi)
   Ecc = Ecc + 0.75d0*(2.0d0*dl2m2*Vdd_delta + l2m2*dVdd_delta)*l2m2
end function d_t_d3z2_r2_d3z2_r2




!------------------------------------------------------------------------
! Hetero-nuclear case:


! Constructing Koster Slater part for s basis set:
pure subroutine KS_s_hetero(Vsr12, ts)
   real(8), intent(in) :: Vsr12 ! coefficient for this pair of atoms
   real(8), dimension(1,1), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)

   ! functions above:
   ts(1,1) = t_s_s(Vsr12)
end subroutine KS_s_hetero

pure subroutine d_KS_s_hetero(dVsr12, ts)
   real(8), intent(in) :: dVsr12 ! coefficient for this pair of atoms
   real(8), dimension(1,1), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12)
end subroutine d_KS_s_hetero



! Constructing Koster Slater part for ss* basis set:
pure subroutine KS_ss_hetero(Vsr12, Vsr21, ts)
   real(8), dimension(3), intent(in) :: Vsr12, Vsr21   ! two sets of coefficients for this pair of atoms
   real(8), dimension(2,2), intent(out) :: ts   ! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s* s sigma)
!    V(3) = (s* s* sigma)
   ! functions above:
   ts(1,1) = t_s_s(Vsr12(1))
   ts(1,2) = t_s_s(Vsr12(2))
   ts(2,1) = t_s_s(Vsr21(2))
   ts(2,2) = t_s_s(Vsr12(3))
end subroutine KS_ss_hetero

pure subroutine d_KS_ss_hetero(dVsr12, dVsr21, ts)
   real(8), dimension(3), intent(in) :: dVsr12, dVsr21   ! two sets of coefficients for this pair of atoms
   real(8), dimension(2,2), intent(out) :: ts   ! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s* s sigma)
!    V(3) = (s* s* sigma)
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))
   ts(1,2) = d_t_s_s(dVsr12(2))
   ts(2,1) = d_t_s_s(dVsr21(2))
   ts(2,2) = d_t_s_s(dVsr12(3))
end subroutine d_KS_ss_hetero


! Constructing Koster Slater part for sp3 basis set:
pure subroutine KS_sp3_hetero(Vsr12, Vsr21, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(4), intent(in) :: Vsr12, Vsr21    ! two sets of coefficients for this pair of atoms
   real(8), dimension(4,4), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)

   ! functions above:
   ts(1,1) = t_s_s(Vsr12(1))
   ts(1,2) = -t_s_px(l, Vsr12(2))
   ts(1,3) = -t_s_px(m, Vsr12(2))
   ts(1,4) = -t_s_px(n, Vsr12(2))
   
   ts(2,1) = t_s_px(l, Vsr21(2))
   ts(2,2) = t_pa_pa(l, Vsr12(3), Vsr12(4))
   ts(2,3) = t_pa_pb(l, m, Vsr12(3), Vsr12(4))
   ts(2,4) = t_pa_pb(l, n, Vsr12(3), Vsr12(4))
   
   ts(3,1) = t_s_px(m, Vsr21(2))
   ts(3,2) = t_pa_pb(l, m, Vsr21(3), Vsr21(4))
   ts(3,3) = t_pa_pa(m, Vsr12(3), Vsr12(4))
   ts(3,4) = t_pa_pb(m, n, Vsr12(3), Vsr12(4))
   
   ts(4,1) = t_s_px(n, Vsr21(2))
   ts(4,2) = t_pa_pb(l, n, Vsr21(3), Vsr21(4))
   ts(4,3) = t_pa_pb(m, n, Vsr21(3), Vsr21(4))
   ts(4,4) = t_pa_pa(n, Vsr12(3), Vsr12(4))
end subroutine KS_sp3_hetero


! Constructing Koster Slater part for sp3 basis set:
pure subroutine d_KS_sp3_hetero(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives by r_{ij}
   real(8), dimension(4), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21	! coefficients and derivatives by r_{ij} for this pair of atoms
   real(8), dimension(4,4), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))
   ts(1,2) = d_t_s_px(l, dl, Vsr12(2), dVsr12(2))
   ts(1,3) = d_t_s_px(m, dm, Vsr12(2), dVsr12(2))
   ts(1,4) = d_t_s_px(n, dn, Vsr12(2), dVsr12(2))
   
   ts(2,1) = -d_t_s_px(l, dl, Vsr21(2), dVsr21(2))
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   
   ts(3,1) = -d_t_s_px(m, dm, Vsr21(2), dVsr21(2))
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   
   ts(4,1) = -d_t_s_px(n, dn, Vsr21(2), dVsr21(2))
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
end subroutine d_KS_sp3_hetero



! Constructing Koster Slater part for sp3s* basis set:
pure subroutine KS_sp3s_hetero(Vsr12, Vsr21, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(7), intent(in) :: Vsr12, Vsr21		! coefficients for this pair of atoms
   real(8), dimension(5,5), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)
!    V(5) = (s* s sigma)
!    V(6) = (s* p sigma)
!    V(7) = (s* s* sigma)

   ! functions above:
   ts(1,1) = t_s_s(Vsr12(1))
   ts(1,2) = -t_s_px(l, Vsr12(2))
   ts(1,3) = -t_s_px(m, Vsr12(2))
   ts(1,4) = -t_s_px(n, Vsr12(2))
   ts(1,5) = t_s_s(Vsr12(5))
   
   ts(2,1) = t_s_px(l, Vsr21(2))
   ts(2,2) = t_pa_pa(l, Vsr12(3), Vsr12(4))
   ts(2,3) = t_pa_pb(l, m, Vsr12(3), Vsr12(4))
   ts(2,4) = t_pa_pb(l, n, Vsr12(3), Vsr12(4))
   ts(2,5) = -t_s_px(l, Vsr12(6))
   
   ts(3,1) = t_s_px(m, Vsr21(2))
   ts(3,2) = t_pa_pb(l, m, Vsr21(3), Vsr21(4))
   ts(3,3) = t_pa_pa(m, Vsr12(3), Vsr12(4))
   ts(3,4) = t_pa_pb(m, n, Vsr12(3), Vsr12(4))
   ts(3,5) = -t_s_px(m, Vsr12(6))
   
   ts(4,1) = t_s_px(n, Vsr21(2))
   ts(4,2) = t_pa_pb(l, n, Vsr21(3), Vsr21(4))
   ts(4,3) = t_pa_pb(m, n, Vsr21(3), Vsr21(4))
   ts(4,4) = t_pa_pa(n, Vsr12(3), Vsr12(4))
   ts(4,5) = -t_s_px(n, Vsr12(6))
   
   ts(5,1) = t_s_s(Vsr21(5))
   ts(5,2) = t_s_px(l, Vsr21(6))
   ts(5,3) = t_s_px(m, Vsr21(6))
   ts(5,4) = t_s_px(n, Vsr21(6))
   ts(5,5) = t_s_s(Vsr12(7))
end subroutine KS_sp3s_hetero


pure subroutine d_KS_sp3s_hetero(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and their derivatives
   real(8), dimension(7), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21   ! coefficients for this pair of atoms and derivatives
   real(8), dimension(5,5), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)
!    V(5) = (s* s sigma)
!    V(6) = (s* p sigma)
!    V(7) = (s* s* sigma)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))
   ts(1,2) = d_t_s_px(l, dl, Vsr12(2), dVsr12(2))
   ts(1,3) = d_t_s_px(m, dm, Vsr12(2), dVsr12(2))
   ts(1,4) = d_t_s_px(n, dn, Vsr12(2), dVsr12(2))
   ts(1,5) = d_t_s_s(dVsr12(5))
   
   ts(2,1) = -d_t_s_px(l, dl, Vsr21(2), dVsr21(2))
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,5) = d_t_s_px(l, dl, Vsr12(6), dVsr12(6))
   
   ts(3,1) = -d_t_s_px(m, dm, Vsr21(2), dVsr21(2))
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(3,5) = d_t_s_px(m, dm, Vsr12(6), dVsr12(6))
   
   ts(4,1) = -d_t_s_px(n, dn, Vsr21(2), dVsr21(2))
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(4,5) = d_t_s_px(n, dn, Vsr12(6), dVsr12(6))
   
   ts(5,1) = d_t_s_s(dVsr21(5))
   ts(5,2) = -d_t_s_px(l, dl, Vsr21(6), dVsr21(6))
   ts(5,3) = -d_t_s_px(m, dm, Vsr21(6), dVsr21(6))
   ts(5,4) = -d_t_s_px(n, dn, Vsr21(6), dVsr21(6))
   ts(5,5) = d_t_s_s(dVsr12(7))
end subroutine d_KS_sp3s_hetero


! Constructing Koster Slater part for sp3d5 basis set:
pure subroutine KS_sp3d5_hetero(Vsr12, Vsr21, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(10), intent(in) :: Vsr12, Vsr21		! coefficients for this pair of atoms
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = t_s_s(Vsr12(1))		! s s
   ts(1,2) = -t_s_px(l, Vsr12(2))	! s px
   ts(1,3) = -t_s_px(m, Vsr12(2))	! s py
   ts(1,4) = -t_s_px(n, Vsr12(2))	! s pz
   ts(1,5) = t_s_dab(l, m, Vsr12(3))	! s dxy
   ts(1,6) = t_s_dab(l, n, Vsr12(3))	! s dxz
   ts(1,7) = t_s_dab(m, n, Vsr12(3))	! s dyz
   ts(1,8) = t_s_dx2_y2(l, m, Vsr12(3))	! s dx2-y2
   ts(1,9) = t_s_dz2_r2(l, m, n, Vsr12(3))	! s d3z2-r2
   
   
   ts(2,1) = t_s_px(l, Vsr21(2))		! px s
   ts(2,2) = t_pa_pa(l, Vsr12(4), Vsr12(5))		! px px
   ts(2,3) = t_pa_pb(l, m, Vsr12(4), Vsr12(5))	! px py
   ts(2,4) = t_pa_pb(l, n, Vsr12(4), Vsr12(5))	! px pz
   ts(2,5) = -t_px_dxy(l, m, Vsr12(6), Vsr12(7))	! px dxy
   ts(2,6) = -t_px_dxy(l, n, Vsr12(6), Vsr12(7))	! px dxz
   ts(2,7) = -t_px_dyz(l, m, n, Vsr12(6), Vsr12(7))	! px dyz
   ts(2,8) = -t_px_dx2_y2(l, m, Vsr12(6), Vsr12(7))	! px dx2-y2
   ts(2,9) = -t_pa_d3z2_r2(l, l, m, n, Vsr12(6), Vsr12(7))	! px d3z2-r2
   
   
   ts(3,1) = t_s_px(m, Vsr21(2))	! py s
   ts(3,2) = t_pa_pb(l, m, Vsr21(4), Vsr21(5))	! py px
   ts(3,3) = t_pa_pa(m, Vsr12(4), Vsr12(5))		! py py
   ts(3,4) = t_pa_pb(m, n, Vsr12(4), Vsr12(5))	! py pz
   ts(3,5) = -t_px_dxy(m, l, Vsr12(6), Vsr12(7))	! py dxy
   ts(3,6) = -t_px_dyz(m, l, n, Vsr12(6), Vsr12(7))	! py dxz
   ts(3,7) = -t_px_dxy(m, n,  Vsr12(6), Vsr12(7))	! py dyz
   ts(3,8) = -t_py_dx2_y2(l, m, Vsr12(6), Vsr12(7))	! py dx2-y2
   ts(3,9) = -t_pa_d3z2_r2(m, l, m, n, Vsr12(6), Vsr12(7))	! py d3z2-r2
   
   
   ts(4,1) = t_s_px(n, Vsr21(2))	! pz s
   ts(4,2) = t_pa_pb(l, n, Vsr21(4), Vsr21(5))	! pz px 
   ts(4,3) = t_pa_pb(m, n, Vsr21(4), Vsr21(5))	! pz py 
   ts(4,4) = t_pa_pa(n, Vsr12(4), Vsr12(5))		! pz pz
   ts(4,5) = -t_px_dyz(n, l, m, Vsr12(6), Vsr12(7))	! pz dxy
   ts(4,6) = -t_px_dxy(n, l, Vsr12(6), Vsr12(7))	! pz dxz
   ts(4,7) = -t_px_dxy(n, m, Vsr12(6), Vsr12(7))	! pz dyz
   ts(4,8) = -t_pz_dx2_y2(l, m, n, Vsr12(6), Vsr12(7))	! pz dx2-y2
   ts(4,9) = -t_pz_d3z2_r2(l, m, n, Vsr12(6), Vsr12(7))	! pz d3z2-r2

   
   ts(5,1) = t_s_dab(l, m, Vsr21(3))	! dxy s 
   ts(5,2) = t_px_dxy(l, m, Vsr21(6), Vsr21(7))	! dxy px 
   ts(5,3) = t_px_dxy(m, l, Vsr21(6), Vsr21(7))	! dxy py 
   ts(5,4) = t_px_dyz(n, l, m, Vsr21(6), Vsr21(7))	! dxy pz 
   ts(5,5) = t_dab_dab(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dxy
   ts(5,6) = t_dab_dbg(m, l, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dxz
   ts(5,7) = t_dab_dbg(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dyz
   ts(5,8) = t_dxy_dx2_y2(l, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dx2-y2
   ts(5,9) = t_dxy_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = t_s_dab(l, n, Vsr21(3))	! dxz s 
   ts(6,2) = t_px_dxy(l, n, Vsr21(6), Vsr21(7))	! dxz px 
   ts(6,3) = t_px_dyz(m, l, n, Vsr21(6), Vsr21(7))	! dxz py 
   ts(6,4) = t_px_dxy(n, l, Vsr21(6), Vsr21(7))	! dxz pz 
   ts(6,5) = t_dab_dbg(m, l, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dxz dxy 
   ts(6,6) = t_dab_dab(l, n, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dxz
   ts(6,7) = t_dab_dbg(l, n, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dyz
   ts(6,8) = t_dxz_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dx2-y2
   ts(6,9) = t_dxz_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dxz d3z2-r2

   
   ts(7,1) = t_s_dab(m, n, Vsr21(3))	! dyz s 
   ts(7,2) = t_px_dyz(l, m, n, Vsr21(6), Vsr21(7))	! dyz px 
   ts(7,3) = t_px_dxy(m, n,  Vsr21(6), Vsr21(7))	! dyz py 
   ts(7,4) = t_px_dxy(n, m, Vsr21(6), Vsr21(7))	! dyz pz 
   ts(7,5) = t_dab_dbg(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dyz dxy 
   ts(7,6) = t_dab_dbg(l, n, m, Vsr21(8), Vsr21(9), Vsr21(10))	! dyz dxz 
   ts(7,7) = t_dab_dab(m, n, l, Vsr12(8), Vsr12(9), Vsr12(10))	! dyz dyz
   ts(7,8) = t_dyz_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dyz dx2-y2
   ts(7,9) = t_dyz_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dyz d3z2-r2

   
   ts(8,1) = t_s_dx2_y2(l, m, Vsr21(3))	! dx2-y2 s 
   ts(8,2) = t_px_dx2_y2(l, m, Vsr21(6), Vsr21(7))	! dx2-y2 px 
   ts(8,3) = t_py_dx2_y2(l, m, Vsr21(6), Vsr21(7))	! dx2-y2 py
   ts(8,4) = t_pz_dx2_y2(l, m, n, Vsr21(6), Vsr21(7))	! dx2-y2 pz
   ts(8,5) = t_dxy_dx2_y2(l, m, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dxy 
   ts(8,6) = t_dxz_dx2_y2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dxz 
   ts(8,7) = t_dyz_dx2_y2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dyz
   ts(8,8) = t_dx2_y2_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = t_dx2_y2_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = t_s_dz2_r2(l, m, n, Vsr21(3))	! d3z2-r2 s 
   ts(9,2) = t_pa_d3z2_r2(l, l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 px 
   ts(9,3) = t_pa_d3z2_r2(m, l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 py 
   ts(9,4) = t_pz_d3z2_r2(l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 pz 
   ts(9,5) = t_dxy_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dxy 
   ts(9,6) = t_dxz_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dxz 
   ts(9,7) = t_dyz_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dyz 
   ts(9,8) = t_dx2_y2_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dx2-y2
   ts(9,9) = t_d3z2_r2_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! d3z2-r2 d3z2-r2
end subroutine KS_sp3d5_hetero



pure subroutine d_KS_sp3d5_hetero(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives
   real(8), dimension(10), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21	! coefficients for this pair of atoms and derivatives
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
 
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))	! s s
   ts(1,2) = d_t_s_px(l, dl, Vsr12(2), dVsr12(2))	! s px
   ts(1,3) = d_t_s_px(m, dm, Vsr12(2), dVsr12(2))	! s py
   ts(1,4) = d_t_s_px(n, dn, Vsr12(2), dVsr12(2))	! s pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr12(3), dVsr12(3))	! s dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr12(3), dVsr12(3))	! s dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr12(3), dVsr12(3))	! s d3z2-r2
   
   
   ts(2,1) = -d_t_s_px(l, dl, Vsr21(2), dVsr21(2))	! px s 
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! px px
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px py
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px pz
   ts(2,5) = d_t_px_dxy(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxy
   ts(2,6) = d_t_px_dxy(l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxz
   ts(2,7) = d_t_px_dyz(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dyz
   ts(2,8) = d_t_px_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dx2-y2
   ts(2,9) = d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px d3z2-r2
   
   
   ts(3,1) = -d_t_s_px(m, dm, Vsr21(2), dVsr21(2))	! py s 
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! py px 
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! py py
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! py pz
   ts(3,5) = d_t_px_dxy(m, dm, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxy
   ts(3,6) = d_t_px_dyz(m, dm, l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxz
   ts(3,7) = d_t_px_dxy(m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dyz
   ts(3,8) = d_t_py_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dx2-y2
   ts(3,9) = d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py d3z2-r2
   
   
   ts(4,1) = -d_t_s_px(n, dn, Vsr21(2), dVsr21(2))	! pz s 
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz px 
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz py 
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! pz pz
   ts(4,5) = d_t_px_dyz(n, dn, l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxy
   ts(4,6) = d_t_px_dxy(n, dn, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxz
   ts(4,7) = d_t_px_dxy(n, dn, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dyz
   ts(4,8) = d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dx2-y2
   ts(4,9) = d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz d3z2-r2

   
   ts(5,1) = d_t_s_dab(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dxy s 
   ts(5,2) = -d_t_px_dxy(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy px 
   ts(5,3) = -d_t_px_dxy(m, dm, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy py 
   ts(5,4) = -d_t_px_dyz(n, dn, l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy pz 
   ts(5,5) = d_t_dab_dab(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxy
   ts(5,6) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxz
   ts(5,7) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dyz
   ts(5,8) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dx2-y2
   ts(5,9) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = d_t_s_dab(l, dl, n, dn, Vsr21(3), dVsr21(3))	! dxz s 
   ts(6,2) = -d_t_px_dxy(l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz px 
   ts(6,3) = -d_t_px_dyz(m, dm, l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz py 
   ts(6,4) = -d_t_px_dxy(n, dn, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz pz 
   ts(6,5) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dxz dxy 
   ts(6,6) = d_t_dab_dab(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dxz
   ts(6,7) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dyz
   ts(6,8) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dx2-y2
   ts(6,9) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxz d3z2-r2

   
   ts(7,1) = d_t_s_dab(m, dm, n, dn, Vsr21(3), dVsr21(3))	! dyz s 
   ts(7,2) = -d_t_px_dyz(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz px 
   ts(7,3) = -d_t_px_dxy(m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz py 
   ts(7,4) = -d_t_px_dxy(n, dn, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz pz 
   ts(7,5) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxy 
   ts(7,6) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxz 
   ts(7,7) = d_t_dab_dab(m, dm, n, dn, l, dl, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dyz
   ts(7,8) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dx2-y2
   ts(7,9) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dyz d3z2-r2

   
   ts(8,1) = d_t_s_dx2_y2(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dx2-y2 s 
   ts(8,2) = -d_t_px_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 px 
   ts(8,3) = -d_t_py_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 py 
   ts(8,4) = -d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 pz 
   ts(8,5) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxy 
   ts(8,6) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxz 
   ts(8,7) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dyz 
   ts(8,8) = d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr21(3), dVsr21(3))	! d3z2-r2 s 
   ts(9,2) = -d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 px 
   ts(9,3) = -d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 py 
   ts(9,4) = -d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 pz 
   ts(9,5) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxy 
   ts(9,6) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxz 
   ts(9,7) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dyz 
   ts(9,8) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dx2-y2 
   ts(9,9) = d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! d3z2-r2 d3z2-r2
end subroutine d_KS_sp3d5_hetero




! Constructing Koster-Slater part for sp3d5s* basis set:
pure subroutine KS_sp3d5ss_hetero(Vsr12, Vsr21, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(14), intent(in) :: Vsr12, Vsr21		! coefficients for this pair of atoms
   real(8), dimension(10,10), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
!    V(11) = (s* s sigma)
!    V(12) = (s* p sigma)
!    V(13) = (s* d sigma)
!    V(14) = (s* s* sigma)

   ! functions above:
   ts(1,1) = t_s_s(Vsr12(1))		! s s
   ts(1,2) = -t_s_px(l, Vsr12(2))		! s px
   ts(1,3) = -t_s_px(m, Vsr12(2))	! s py
   ts(1,4) = -t_s_px(n, Vsr12(2))	! s pz
   ts(1,5) = t_s_dab(l, m, Vsr12(3))	! s dxy
   ts(1,6) = t_s_dab(l, n, Vsr12(3))	! s dxz
   ts(1,7) = t_s_dab(m, n, Vsr12(3))	! s dyz
   ts(1,8) = t_s_dx2_y2(l, m, Vsr12(3))	! s dx2-y2
   ts(1,9) = t_s_dz2_r2(l, m, n, Vsr12(3))	! s d3z2-r2
   ts(1,10) = t_s_s(Vsr12(11))		! s s*


   ts(2,1) = t_s_px(l, Vsr21(2))		! px s
   ts(2,2) = t_pa_pa(l, Vsr12(4), Vsr12(5))		! px px
   ts(2,3) = t_pa_pb(l, m, Vsr12(4), Vsr12(5))	! px py
   ts(2,4) = t_pa_pb(l, n, Vsr12(4), Vsr12(5))	! px pz
   ts(2,5) = -t_px_dxy(l, m, Vsr12(6), Vsr12(7))	! px dxy
   ts(2,6) = -t_px_dxy(l, n, Vsr12(6), Vsr12(7))	! px dxz
   ts(2,7) = -t_px_dyz(l, m, n, Vsr12(6), Vsr12(7))	! px dyz
   ts(2,8) = -t_px_dx2_y2(l, m, Vsr12(6), Vsr12(7))	! px dx2-y2
   ts(2,9) = -t_pa_d3z2_r2(l, l, m, n, Vsr12(6), Vsr12(7))	! px d3z2-r2
   ts(2,10) = -t_s_px(l, Vsr12(12))   ! px s*


   ts(3,1) = t_s_px(m, Vsr21(2))	! py s
   ts(3,2) = t_pa_pb(l, m, Vsr21(4), Vsr21(5))	! py px
   ts(3,3) = t_pa_pa(m, Vsr12(4), Vsr12(5))		! py py
   ts(3,4) = t_pa_pb(m, n, Vsr12(4), Vsr12(5))	! py pz
   ts(3,5) = -t_px_dxy(m, l, Vsr12(6), Vsr12(7))	! py dxy
   ts(3,6) = -t_px_dyz(m, l, n, Vsr12(6), Vsr12(7))	! py dxz
   ts(3,7) = -t_px_dxy(m, n,  Vsr12(6), Vsr12(7))	! py dyz
   ts(3,8) = -t_py_dx2_y2(l, m, Vsr12(6), Vsr12(7))	! py dx2-y2
   ts(3,9) = -t_pa_d3z2_r2(m, l, m, n, Vsr12(6), Vsr12(7))	! py d3z2-r2
   ts(3,10) = -t_s_px(m, Vsr12(12))   ! py s*


   ts(4,1) = t_s_px(n, Vsr21(2))	! pz s
   ts(4,2) = t_pa_pb(l, n, Vsr21(4), Vsr21(5))	! pz px
   ts(4,3) = t_pa_pb(m, n, Vsr21(4), Vsr21(5))	! pz py
   ts(4,4) = t_pa_pa(n, Vsr12(4), Vsr12(5))		! pz pz
   ts(4,5) = -t_px_dyz(n, l, m, Vsr12(6), Vsr12(7))	! pz dxy
   ts(4,6) = -t_px_dxy(n, l, Vsr12(6), Vsr12(7))	! pz dxz
   ts(4,7) = -t_px_dxy(n, m, Vsr12(6), Vsr12(7))	! pz dyz
   ts(4,8) = -t_pz_dx2_y2(l, m, n, Vsr12(6), Vsr12(7))	! pz dx2-y2
   ts(4,9) = -t_pz_d3z2_r2(l, m, n, Vsr12(6), Vsr12(7))	! pz d3z2-r2
   ts(4,10) = -t_s_px(m, Vsr12(12))   ! pz s*


   ts(5,1) = t_s_dab(l, m, Vsr21(3))	! dxy s
   ts(5,2) = t_px_dxy(l, m, Vsr21(6), Vsr21(7))	! dxy px
   ts(5,3) = t_px_dxy(m, l, Vsr21(6), Vsr21(7))	! dxy py
   ts(5,4) = t_px_dyz(n, l, m, Vsr21(6), Vsr21(7))	! dxy pz
   ts(5,5) = t_dab_dab(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dxy
   ts(5,6) = t_dab_dbg(m, l, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dxz
   ts(5,7) = t_dab_dbg(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dyz
   ts(5,8) = t_dxy_dx2_y2(l, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dx2-y2
   ts(5,9) = t_dxy_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dxy d3z2-r2
   ts(5,10) = t_s_dab(l, m, Vsr12(13))	! dxy s*


   ts(6,1) = t_s_dab(l, n, Vsr21(3))	! dxz s
   ts(6,2) = t_px_dxy(l, n, Vsr21(6), Vsr21(7))	! dxz px
   ts(6,3) = t_px_dyz(m, l, n, Vsr21(6), Vsr21(7))	! dxz py
   ts(6,4) = t_px_dxy(n, l, Vsr21(6), Vsr21(7))	! dxz pz
   ts(6,5) = t_dab_dbg(m, l, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dxz dxy
   ts(6,6) = t_dab_dab(l, n, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dxz
   ts(6,7) = t_dab_dbg(l, n, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dyz
   ts(6,8) = t_dxz_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dx2-y2
   ts(6,9) = t_dxz_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dxz d3z2-r2
   ts(6,10) = t_s_dab(l, n, Vsr12(13))	! dxz s*


   ts(7,1) = t_s_dab(m, n, Vsr21(3))	! dyz s
   ts(7,2) = t_px_dyz(l, m, n, Vsr21(6), Vsr21(7))	! dyz px
   ts(7,3) = t_px_dxy(m, n,  Vsr21(6), Vsr21(7))	! dyz py
   ts(7,4) = t_px_dxy(n, m, Vsr21(6), Vsr21(7))	! dyz pz
   ts(7,5) = t_dab_dbg(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dyz dxy
   ts(7,6) = t_dab_dbg(l, n, m, Vsr21(8), Vsr21(9), Vsr21(10))	! dyz dxz
   ts(7,7) = t_dab_dab(m, n, l, Vsr12(8), Vsr12(9), Vsr12(10))	! dyz dyz
   ts(7,8) = t_dyz_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dyz dx2-y2
   ts(7,9) = t_dyz_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dyz d3z2-r2
   ts(7,10) = t_s_dab(m, n, Vsr12(13))	! dyz s*


   ts(8,1) = t_s_dx2_y2(l, m, Vsr21(3))	! dx2-y2 s
   ts(8,2) = t_px_dx2_y2(l, m, Vsr21(6), Vsr21(7))	! dx2-y2 px
   ts(8,3) = t_py_dx2_y2(l, m, Vsr21(6), Vsr21(7))	! dx2-y2 py
   ts(8,4) = t_pz_dx2_y2(l, m, n, Vsr21(6), Vsr21(7))	! dx2-y2 pz
   ts(8,5) = t_dxy_dx2_y2(l, m, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dxy
   ts(8,6) = t_dxz_dx2_y2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dxz
   ts(8,7) = t_dyz_dx2_y2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dyz
   ts(8,8) = t_dx2_y2_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = t_dx2_y2_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dx2-y2 d3z2-r2
   ts(8,10) = t_s_dx2_y2(l, m, Vsr12(13))	! dx2-y2 s*


   ts(9,1) = t_s_dz2_r2(l, m, n, Vsr21(3))	! d3z2-r2 s
   ts(9,2) = t_pa_d3z2_r2(l, l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 px
   ts(9,3) = t_pa_d3z2_r2(m, l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 py
   ts(9,4) = t_pz_d3z2_r2(l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 pz
   ts(9,5) = t_dxy_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dxy
   ts(9,6) = t_dxz_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dxz
   ts(9,7) = t_dyz_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dyz
   ts(9,8) = t_dx2_y2_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dx2-y2
   ts(9,9) = t_d3z2_r2_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! d3z2-r2 d3z2-r2
   ts(9,10) = t_s_dz2_r2(l, m, n, Vsr12(13))	! d3z2-r2 s*


   ts(10,1) = t_s_s(Vsr21(11))      ! s* s
   ts(10,2) = t_s_px(l, Vsr21(12))  ! s* px
   ts(10,3) = t_s_px(m, Vsr21(12))  ! s* py
   ts(10,4) = t_s_px(n, Vsr21(12))  ! s* pz
   ts(10,5) = t_s_dab(l, m, Vsr21(13))  ! s* dxy
   ts(10,6) = t_s_dab(l, n, Vsr21(13))  ! s* dxz
   ts(10,7) = t_s_dab(m, n, Vsr21(13))  ! s* dyz
   ts(10,8) = t_s_dx2_y2(l, m, Vsr21(13))   ! s* dx2-y2
   ts(10,9) = t_s_dz2_r2(l, m, n, Vsr21(13))    ! s* d3z2-r2
   ts(10,10) = t_s_s(Vsr12(14)) ! s* s*
end subroutine KS_sp3d5ss_hetero


pure subroutine d_KS_sp3d5ss_hetero(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives
   real(8), dimension(14), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21	! coefficients for this pair of atoms and derivatives
   real(8), dimension(10,10), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
!    V(11) = (s* s sigma)
!    V(12) = (s* p sigma)
!    V(13) = (s* d sigma)
!    V(14) = (s* s* sigma)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))	! s s
   ts(1,2) = d_t_s_px(l, dl, Vsr12(2), dVsr12(2))	! s px
   ts(1,3) = d_t_s_px(m, dm, Vsr12(2), dVsr12(2))	! s py
   ts(1,4) = d_t_s_px(n, dn, Vsr12(2), dVsr12(2))	! s pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr12(3), dVsr12(3))	! s dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr12(3), dVsr12(3))	! s dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr12(3), dVsr12(3))	! s d3z2-r2
   ts(1,10) = d_t_s_s(dVsr12(11))	! s s*


   ts(2,1) = -d_t_s_px(l, dl, Vsr21(2), dVsr21(2))	! px s
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! px px
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px py
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px pz
   ts(2,5) = d_t_px_dxy(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxy
   ts(2,6) = d_t_px_dxy(l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxz
   ts(2,7) = d_t_px_dyz(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dyz
   ts(2,8) = d_t_px_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dx2-y2
   ts(2,9) = d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px d3z2-r2
   ts(2,10) = -d_t_s_px(l, dl, Vsr12(12), dVsr12(12))	! px s


   ts(3,1) = -d_t_s_px(m, dm, Vsr21(2), dVsr21(2))	! py s
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! py px
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! py py
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! py pz
   ts(3,5) = d_t_px_dxy(m, dm, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxy
   ts(3,6) = d_t_px_dyz(m, dm, l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxz
   ts(3,7) = d_t_px_dxy(m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dyz
   ts(3,8) = d_t_py_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dx2-y2
   ts(3,9) = d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py d3z2-r2
   ts(3,10) = -d_t_s_px(m, dm, Vsr12(12), dVsr12(12))	! py s*


   ts(4,1) = -d_t_s_px(n, dn, Vsr21(2), dVsr21(2))	! pz s
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz px
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz py
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! pz pz
   ts(4,5) = d_t_px_dyz(n, dn, l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxy
   ts(4,6) = d_t_px_dxy(n, dn, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxz
   ts(4,7) = d_t_px_dxy(n, dn, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dyz
   ts(4,8) = d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dx2-y2
   ts(4,9) = d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz d3z2-r2
   ts(4,10) = -d_t_s_px(n, dn, Vsr12(12), dVsr12(12))	! pz s*


   ts(5,1) = d_t_s_dab(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dxy s
   ts(5,2) = -d_t_px_dxy(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy px
   ts(5,3) = -d_t_px_dxy(m, dm, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy py
   ts(5,4) = -d_t_px_dyz(n, dn, l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy pz
   ts(5,5) = d_t_dab_dab(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxy
   ts(5,6) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxz
   ts(5,7) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dyz
   ts(5,8) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dx2-y2
   ts(5,9) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxy d3z2-r2
   ts(5,10) = d_t_s_dab(l, dl, m, dm, Vsr12(13), dVsr12(13))	! dxy s*


   ts(6,1) = d_t_s_dab(l, dl, n, dn, Vsr21(3), dVsr21(3))	! dxz s
   ts(6,2) = -d_t_px_dxy(l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz px
   ts(6,3) = -d_t_px_dyz(m, dm, l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz py
   ts(6,4) = -d_t_px_dxy(n, dn, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz pz
   ts(6,5) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dxz dxy
   ts(6,6) = d_t_dab_dab(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dxz
   ts(6,7) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dyz
   ts(6,8) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dx2-y2
   ts(6,9) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxz d3z2-r2
   ts(6,10) = d_t_s_dab(l, dl, n, dn, Vsr12(13), dVsr12(13))	! dxz s*


   ts(7,1) = d_t_s_dab(m, dm, n, dn, Vsr21(3), dVsr21(3))	! dyz s
   ts(7,2) = -d_t_px_dyz(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz px
   ts(7,3) = -d_t_px_dxy(m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz py
   ts(7,4) = -d_t_px_dxy(n, dn, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz pz
   ts(7,5) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxy
   ts(7,6) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxz
   ts(7,7) = d_t_dab_dab(m, dm, n, dn, l, dl, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dyz
   ts(7,8) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dx2-y2
   ts(7,9) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dyz d3z2-r2
   ts(7,10) = d_t_s_dab(m, dm, n, dn, Vsr12(13), dVsr12(13))	! dyz s*


   ts(8,1) = d_t_s_dx2_y2(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dx2-y2 s
   ts(8,2) = -d_t_px_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 px
   ts(8,3) = -d_t_py_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 py
   ts(8,4) = -d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 pz
   ts(8,5) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxy
   ts(8,6) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxz
   ts(8,7) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dyz
   ts(8,8) = d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dx2-y2 d3z2-r2
   ts(8,10) = d_t_s_dx2_y2(l, dl, m, dm, Vsr12(13), dVsr12(13))	! dx2-y2 s*


   ts(9,1) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr21(3), dVsr21(3))	! d3z2-r2 s
   ts(9,2) = -d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 px
   ts(9,3) = -d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 py
   ts(9,4) = -d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 pz
   ts(9,5) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxy
   ts(9,6) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxz
   ts(9,7) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dyz
   ts(9,8) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dx2-y2
   ts(9,9) = d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! d3z2-r2 d3z2-r2
   ts(9,10) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr12(13), dVsr12(13))	! d3z2-r2 s*

   ts(1,1) = d_t_s_s(dVsr21(1))	! s* s
   ts(1,2) = d_t_s_px(l, dl, Vsr21(2), dVsr21(2))	! s* px
   ts(1,3) = d_t_s_px(m, dm, Vsr21(2), dVsr21(2))	! s* py
   ts(1,4) = d_t_s_px(n, dn, Vsr21(2), dVsr21(2))	! s* pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr21(3), dVsr21(3))	! s* dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr21(3), dVsr21(3))	! s* dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr21(3), dVsr21(3))	! s* dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr21(3), dVsr21(3))	! s* dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr21(3), dVsr21(3))	! s* d3z2-r2
   ts(1,10) = d_t_s_s(dVsr12(14))	! s* s*
end subroutine d_KS_sp3d5ss_hetero




pure subroutine d_KS_sp3d5_hetero_OLD(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives
   real(8), dimension(10), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21	! coefficients for this pair of atoms and derivatives
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))	! s s
   ts(1,2) = -d_t_s_px(l, dl, Vsr12(2), dVsr12(2))	! s px
   ts(1,3) = -d_t_s_px(m, dm, Vsr12(2), dVsr12(2))	! s py
   ts(1,4) = -d_t_s_px(n, dn, Vsr12(2), dVsr12(2))	! s pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr12(3), dVsr12(3))	! s dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr12(3), dVsr12(3))	! s dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr12(3), dVsr12(3))	! s d3z2-r2
   
   
   ts(2,1) = d_t_s_px(l, dl, Vsr21(2), dVsr21(2))	! px s 
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! px px
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px py
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px pz
   ts(2,5) = -d_t_px_dxy(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxy
   ts(2,6) = -d_t_px_dxy(l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxz
   ts(2,7) = -d_t_px_dyz(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dyz
   ts(2,8) = -d_t_px_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dx2-y2
   ts(2,9) = -d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px d3z2-r2
   
   
   ts(3,1) = d_t_s_px(m, dm, Vsr21(2), dVsr21(2))	! py s 
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! py px 
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! py py
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! py pz
   ts(3,5) = -d_t_px_dxy(m, dm, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxy
   ts(3,6) = -d_t_px_dyz(m, dm, l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxz
   ts(3,7) = -d_t_px_dxy(m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dyz
   ts(3,8) = -d_t_py_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dx2-y2
   ts(3,9) = -d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py d3z2-r2
   
   
   ts(4,1) = d_t_s_px(n, dn, Vsr21(2), dVsr21(2))	! pz s 
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz px 
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz py 
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! pz pz
   ts(4,5) = -d_t_px_dyz(n, dn, l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxy
   ts(4,6) = -d_t_px_dxy(n, dn, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxz
   ts(4,7) = -d_t_px_dxy(n, dn, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dyz
   ts(4,8) = -d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dx2-y2
   ts(4,9) = -d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz d3z2-r2

   
   ts(5,1) = d_t_s_dab(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dxy s 
   ts(5,2) = d_t_px_dxy(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy px 
   ts(5,3) = d_t_px_dxy(m, dm, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy py 
   ts(5,4) = d_t_px_dyz(n, dn, l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy pz 
   ts(5,5) = d_t_dab_dab(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxy
   ts(5,6) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxz
   ts(5,7) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dyz
   ts(5,8) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dx2-y2
   ts(5,9) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = d_t_s_dab(l, dl, n, dn, Vsr21(3), dVsr21(3))	! dxz s 
   ts(6,2) = d_t_px_dxy(l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz px 
   ts(6,3) = d_t_px_dyz(m, dm, l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz py 
   ts(6,4) = d_t_px_dxy(n, dn, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz pz 
   ts(6,5) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dxz dxy 
   ts(6,6) = d_t_dab_dab(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dxz
   ts(6,7) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dyz
   ts(6,8) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dx2-y2
   ts(6,9) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxz d3z2-r2

   
   ts(7,1) = d_t_s_dab(m, dm, n, dn, Vsr21(3), dVsr21(3))	! dyz s 
   ts(7,2) = d_t_px_dyz(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz px 
   ts(7,3) = d_t_px_dxy(m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz py 
   ts(7,4) = d_t_px_dxy(n, dn, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz pz 
   ts(7,5) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxy 
   ts(7,6) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxz 
   ts(7,7) = d_t_dab_dab(m, dm, n, dn, l, dl, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dyz
   ts(7,8) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dx2-y2
   ts(7,9) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dyz d3z2-r2

   
   ts(8,1) = d_t_s_dx2_y2(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dx2-y2 s 
   ts(8,2) = d_t_px_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 px 
   ts(8,3) = d_t_py_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 py 
   ts(8,4) = d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 pz 
   ts(8,5) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxy 
   ts(8,6) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxz 
   ts(8,7) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dyz 
   ts(8,8) = d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr21(3), dVsr21(3))	! d3z2-r2 s 
   ts(9,2) = d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 px 
   ts(9,3) = d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 py 
   ts(9,4) = d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 pz 
   ts(9,5) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxy 
   ts(9,6) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxz 
   ts(9,7) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dyz 
   ts(9,8) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dx2-y2 
   ts(9,9) = d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! d3z2-r2 d3z2-r2
end subroutine d_KS_sp3d5_hetero_OLD



!------------------------------------------------------------------------
! Homo-nuclear case:
! Constructing Koster Slater part for s basis set (e.g. H and He atoms):
pure subroutine KS_s(Vsr, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(1), intent(in) :: Vsr		! coefficients for this pair of atoms
   real(8), dimension(1,1), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
   !    V(1) = (s s sigma)
   ! functions above:
   ts(1,1) = t_s_s(Vsr(1))
end subroutine KS_s

! Constructing derivatives of the Koster Slater part for s basis set:
pure subroutine d_KS_s(dVsr, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(1), intent(in) :: dVsr		! coefficients for this pair of atoms
   real(8), dimension(1,1), intent(out) :: ts	! constructed hopping integrals
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr(1))
end subroutine d_KS_s



! Constructing Koster Slater part for sp3 basis set: (TESTED, CORRECT)
pure subroutine KS_sp3(Vsr, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(4), intent(in) :: Vsr		! coefficients for this pair of atoms
   real(8), dimension(4,4), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)

   ! functions above:
   ts(1,1) = t_s_s(Vsr(1))
   ts(1,2) = -t_s_px(l, Vsr(2))
   ts(1,3) = -t_s_px(m, Vsr(2))
   ts(1,4) = -t_s_px(n, Vsr(2))
   
   ts(2,1) = -ts(1,2)
   ts(2,2) = t_pa_pa(l, Vsr(3), Vsr(4))
   ts(2,3) = t_pa_pb(l, m, Vsr(3), Vsr(4))
   ts(2,4) = t_pa_pb(l, n, Vsr(3), Vsr(4))
   
   ts(3,1) = -ts(1,3)
   ts(3,2) = ts(2,3)
   ts(3,3) = t_pa_pa(m, Vsr(3), Vsr(4))
   ts(3,4) = t_pa_pb(m, n, Vsr(3), Vsr(4))
   
   ts(4,1) = -ts(1,4)
   ts(4,2) = ts(2,4)
   ts(4,3) = ts(3,4)
   ts(4,4) = t_pa_pa(n, Vsr(3), Vsr(4))
end subroutine KS_sp3


! Constructing Koster Slater part for sp3 basis set: (TESTED, CORRECT)
pure subroutine d_KS_sp3(Vsr, dVsr, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives by r_{ij}
   real(8), dimension(4), intent(in) :: Vsr, dVsr	! coefficients and derivatives by r_{ij} for this pair of atoms
   real(8), dimension(4,4), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr(1))
   ts(1,2) = d_t_s_px(l, dl, Vsr(2), dVsr(2))
   ts(1,3) = d_t_s_px(m, dm, Vsr(2), dVsr(2))
   ts(1,4) = d_t_s_px(n, dn, Vsr(2), dVsr(2))
   
   ts(2,1) = -ts(1,2)
   ts(2,2) = d_t_pa_pa(l, dl, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   
   ts(3,1) = -ts(1,3)
   ts(3,2) = ts(2,3)
   ts(3,3) = d_t_pa_pa(m, dm, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   
   ts(4,1) = -ts(1,4)
   ts(4,2) = ts(2,4)
   ts(4,3) = ts(3,4)
   ts(4,4) = d_t_pa_pa(n, dn, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
end subroutine d_KS_sp3

!------------------------------------------------------------------------
! Constructing Koster Slater part for sp3s* basis set:
pure subroutine KS_sp3s(Vsr, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(7), intent(in) :: Vsr		! coefficients for this pair of atoms
   real(8), dimension(5,5), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)
!    V(5) = (s* s sigma)
!    V(6) = (s* p sigma)
!    V(7) = (s* s* sigma)

   ! functions above:
   ts(1,1) = t_s_s(Vsr(1))
   ts(1,2) = -t_s_px(l, Vsr(2))
   ts(1,3) = -t_s_px(m, Vsr(2))
   ts(1,4) = -t_s_px(n, Vsr(2))
   ts(1,5) = t_s_s(Vsr(5))
   
   ts(2,1) = -ts(1,2)
   ts(2,2) = t_pa_pa(l, Vsr(3), Vsr(4))
   ts(2,3) = t_pa_pb(l, m, Vsr(3), Vsr(4))
   ts(2,4) = t_pa_pb(l, n, Vsr(3), Vsr(4))
   ts(2,5) = -t_s_px(l, Vsr(6))
   
   ts(3,1) = -ts(1,3)
   ts(3,2) = ts(2,3)
   ts(3,3) = t_pa_pa(m, Vsr(3), Vsr(4))
   ts(3,4) = t_pa_pb(m, n, Vsr(3), Vsr(4))
   ts(3,5) = -t_s_px(m, Vsr(6))
   
   ts(4,1) = -ts(1,4)
   ts(4,2) = ts(2,4)
   ts(4,3) = ts(3,4)
   ts(4,4) = t_pa_pa(n, Vsr(3), Vsr(4))
   ts(4,5) = -t_s_px(n, Vsr(6))
   
   ts(5,1) = ts(1,5)
   ts(5,2) = -ts(2,5)
   ts(5,3) = -ts(3,5)
   ts(5,4) = -ts(4,5)
   ts(5,5) = t_s_s(Vsr(7))
end subroutine KS_sp3s


pure subroutine d_KS_sp3s(Vsr, dVsr, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and their derivatives
   real(8), dimension(7), intent(in) :: Vsr, dVsr	! coefficients for this pair of atoms and derivatives
   real(8), dimension(5,5), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)
!    V(5) = (s* s sigma)
!    V(6) = (s* p sigma)
!    V(7) = (s* s* sigma)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr(1))
   ts(1,2) = d_t_s_px(l, dl, Vsr(2), dVsr(2))
   ts(1,3) = d_t_s_px(m, dm, Vsr(2), dVsr(2))
   ts(1,4) = d_t_s_px(n, dn, Vsr(2), dVsr(2))
   ts(1,5) = d_t_s_s(dVsr(5))
   
   ts(2,1) = -ts(1,2)
   ts(2,2) = d_t_pa_pa(l, dl, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(2,5) = d_t_s_px(l, dl, Vsr(6), dVsr(6))
   
   ts(3,1) = -ts(1,3)
   ts(3,2) = ts(2,3)
   ts(3,3) = d_t_pa_pa(m, dm, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(3,5) = d_t_s_px(m, dm, Vsr(6), dVsr(6))
   
   ts(4,1) = -ts(1,4)
   ts(4,2) = ts(2,4)
   ts(4,3) = ts(3,4)
   ts(4,4) = d_t_pa_pa(n, dn, Vsr(3), dVsr(3), Vsr(4), dVsr(4))
   ts(4,5) = d_t_s_px(n, dn, Vsr(6), dVsr(6))
   
   ts(5,1) = ts(1,5)
   ts(5,2) = -ts(2,5)
   ts(5,3) = -ts(3,5)
   ts(5,4) = -ts(4,5)
   ts(5,5) = d_t_s_s(dVsr(7))
end subroutine d_KS_sp3s

!------------------------------------------------------------------------
! Constructing Koster Slater part for sp3d5 basis set:
pure subroutine KS_sp3d5(Vsr, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(10), intent(in) :: Vsr		! coefficients for this pair of atoms
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = t_s_s(Vsr(1))		! s s
   ts(1,2) = -t_s_px(l, Vsr(2))		! s px
   ts(1,3) = -t_s_px(m, Vsr(2))	! s py
   ts(1,4) = -t_s_px(n, Vsr(2))	! s pz
   ts(1,5) = t_s_dab(l, m, Vsr(3))	! s dxy
   ts(1,6) = t_s_dab(l, n, Vsr(3))	! s dxz
   ts(1,7) = t_s_dab(m, n, Vsr(3))	! s dyz
   ts(1,8) = t_s_dx2_y2(l, m, Vsr(3))	! s dx2-y2
   ts(1,9) = t_s_dz2_r2(l, m, n, Vsr(3))	! s d3z2-r2
   
   
   ts(2,1) = -ts(1,2)	! s px
   ts(2,2) = t_pa_pa(l, Vsr(4), Vsr(5))		! px px
   ts(2,3) = t_pa_pb(l, m, Vsr(4), Vsr(5))	! px py
   ts(2,4) = t_pa_pb(l, n, Vsr(4), Vsr(5))	! px pz
   ts(2,5) = -t_px_dxy(l, m, Vsr(6), Vsr(7))	! px dxy
   ts(2,6) = -t_px_dxy(l, n, Vsr(6), Vsr(7))	! px dxz
   ts(2,7) = -t_px_dyz(l, m, n, Vsr(6), Vsr(7))	! px dyz
   ts(2,8) = -t_px_dx2_y2(l, m, Vsr(6), Vsr(7))	! px dx2-y2
   ts(2,9) = -t_pa_d3z2_r2(l, l, m, n, Vsr(6), Vsr(7))	! px d3z2-r2
   
   
   ts(3,1) = -ts(1,3)	! s py
   ts(3,2) = ts(2,3)	! px py
   ts(3,3) = t_pa_pa(m, Vsr(4), Vsr(5))		! py py
   ts(3,4) = t_pa_pb(m, n, Vsr(4), Vsr(5))	! py pz
   ts(3,5) = -t_px_dxy(m, l, Vsr(6), Vsr(7))	! py dxy
   ts(3,6) = -t_px_dyz(m, l, n, Vsr(6), Vsr(7))	! py dxz
   ts(3,7) = -t_px_dxy(m, n,  Vsr(6), Vsr(7))	! py dyz
   ts(3,8) = -t_py_dx2_y2(l, m, Vsr(6), Vsr(7))	! py dx2-y2
   ts(3,9) = -t_pa_d3z2_r2(m, l, m, n, Vsr(6), Vsr(7))	! py d3z2-r2
   
   
   ts(4,1) = -ts(1,4)	! s pz
   ts(4,2) = ts(2,4)	! px pz
   ts(4,3) = ts(3,4)	! py pz
   ts(4,4) = t_pa_pa(n, Vsr(4), Vsr(5))		! pz pz
   ts(4,5) = -t_px_dyz(n, l, m, Vsr(6), Vsr(7))	! pz dxy
   ts(4,6) = -t_px_dxy(n, l, Vsr(6), Vsr(7))	! pz dxz
   ts(4,7) = -t_px_dxy(n, m, Vsr(6), Vsr(7))	! pz dyz
   ts(4,8) = -t_pz_dx2_y2(l, m, n, Vsr(6), Vsr(7))	! pz dx2-y2
   ts(4,9) = -t_pz_d3z2_r2(l, m, n, Vsr(6), Vsr(7))	! pz d3z2-r2

   
   ts(5,1) = ts(1,5)	! s dxy
   ts(5,2) = -ts(2,5)	! px dxy
   ts(5,3) = -ts(3,5)	! py dxy
   ts(5,4) = -ts(4,5)	! pz dxy
   ts(5,5) = t_dab_dab(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dxy dxy
   ts(5,6) = t_dab_dbg(m, l, n, Vsr(8), Vsr(9), Vsr(10))	! dxy dxz
   ts(5,7) = t_dab_dbg(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dxy dyz
   ts(5,8) = t_dxy_dx2_y2(l, m, Vsr(8), Vsr(9), Vsr(10))	! dxy dx2-y2
   ts(5,9) = t_dxy_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = ts(1,6)	! s dxz
   ts(6,2) = -ts(2,6)	! px dxz
   ts(6,3) = -ts(3,6)	! py dxz
   ts(6,4) = -ts(4,6)	! pz dxz
   ts(6,5) = ts(5,6)	! dxy dxz
   ts(6,6) = t_dab_dab(l, n, m, Vsr(8), Vsr(9), Vsr(10))	! dxz dxz
   ts(6,7) = t_dab_dbg(l, n, m, Vsr(8), Vsr(9), Vsr(10))	! dxz dyz
   ts(6,8) = t_dxz_dx2_y2(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dxz dx2-y2
   ts(6,9) = t_dxz_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dxz d3z2-r2

   
   ts(7,1) = ts(1,7)	! s dyz
   ts(7,2) = -ts(2,7)	! px dyz
   ts(7,3) = -ts(3,7)	! py dyz
   ts(7,4) = -ts(4,7)	! pz dyz
   ts(7,5) = ts(5,7)	! dxy dyz
   ts(7,6) = ts(6,7)	! dyz dxz
   ts(7,7) = t_dab_dab(m, n, l, Vsr(8), Vsr(9), Vsr(10))	! dyz dyz
   ts(7,8) = t_dyz_dx2_y2(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dyz dx2-y2
   ts(7,9) = t_dyz_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dyz d3z2-r2

   
   ts(8,1) = ts(1,8)	! s dx2-y2
   ts(8,2) = -ts(2,8)	! px dx2-y2
   ts(8,3) = -ts(3,8)	! py dx2-y2
   ts(8,4) = -ts(4,8)	! pz dx2-y2
   ts(8,5) = ts(5,8)	! dxy dx2-y2
   ts(8,6) = ts(6,8)	! dyz dx2-y2
   ts(8,7) = ts(7,8)	! dyz dx2-y2
   ts(8,8) = t_dx2_y2_dx2_y2(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dx2-y2 dx2-y2
   ts(8,9) = t_dx2_y2_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = ts(1,9)	! s d3z2-r2
   ts(9,2) = -ts(2,9)	! px d3z2-r2
   ts(9,3) = -ts(3,9)	! py d3z2-r2
   ts(9,4) = -ts(4,9)	! pz d3z2-r2
   ts(9,5) = ts(5,9)	! dxy d3z2-r2
   ts(9,6) = ts(6,9)	! dyz d3z2-r2
   ts(9,7) = ts(7,9)	! dyz d3z2-r2
   ts(9,8) = ts(8,9)	! dx2-y2 d3z2-r2
   ts(9,9) = t_d3z2_r2_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! d3z2-r2 d3z2-r2
end subroutine KS_sp3d5



pure subroutine d_KS_sp3d5(Vsr, dVsr, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives
   real(8), dimension(10), intent(in) :: Vsr, dVsr	! coefficients for this pair of atoms and derivatives
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr(1))	! s s
   ts(1,2) = -d_t_s_px(l, dl, Vsr(2), dVsr(2))	! s px
   ts(1,3) = -d_t_s_px(m, dm, Vsr(2), dVsr(2))	! s py
   ts(1,4) = -d_t_s_px(n, dn, Vsr(2), dVsr(2))	! s pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr(3), dVsr(3))	! s dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr(3), dVsr(3))	! s dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr(3), dVsr(3))	! s dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr(3), dVsr(3))	! s dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr(3), dVsr(3))	! s d3z2-r2
   
   
   ts(2,1) = -ts(1,2)	! s px
   ts(2,2) = d_t_pa_pa(l, dl, Vsr(4), dVsr(4), Vsr(5), dVsr(5))		! px px
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr(4), dVsr(4), Vsr(5), dVsr(5))	! px py
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr(4), dVsr(4), Vsr(5), dVsr(5))	! px pz
   ts(2,5) = -d_t_px_dxy(l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dxy
   ts(2,6) = -d_t_px_dxy(l, dl, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dxz
   ts(2,7) = -d_t_px_dyz(l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dyz
   ts(2,8) = -d_t_px_dx2_y2(l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dx2-y2
   ts(2,9) = -d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px d3z2-r2
   
   
   ts(3,1) = -ts(1,3)	! s py
   ts(3,2) = ts(2,3)	! px py
   ts(3,3) = d_t_pa_pa(m, dm, Vsr(4), dVsr(4), Vsr(5), dVsr(5))		! py py
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr(4), dVsr(4), Vsr(5), dVsr(5))	! py pz
   ts(3,5) = -d_t_px_dxy(m, dm, l, dl, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dxy
   ts(3,6) = -d_t_px_dyz(m, dm, l, dl, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dxz
   ts(3,7) = -d_t_px_dxy(m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dyz
   ts(3,8) = -d_t_py_dx2_y2(l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dx2-y2
   ts(3,9) = -d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py d3z2-r2
   
   
   ts(4,1) = -ts(1,4)	! s pz
   ts(4,2) = ts(2,4)	! px pz
   ts(4,3) = ts(3,4)	! py pz
   ts(4,4) = d_t_pa_pa(n, dn, Vsr(4), dVsr(4), Vsr(5), dVsr(5))		! pz pz
   ts(4,5) = -d_t_px_dyz(n, dn, l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dxy
   ts(4,6) = -d_t_px_dxy(n, dn, l, dl, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dxz
   ts(4,7) = -d_t_px_dxy(n, dn, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dyz
   ts(4,8) = -d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dx2-y2
   ts(4,9) = -d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz d3z2-r2

   
   ts(5,1) = ts(1,5)	! s dxy
   ts(5,2) = -ts(2,5)	! px dxy
   ts(5,3) = -ts(3,5)	! py dxy
   ts(5,4) = -ts(4,5)	! pz dxy
   ts(5,5) = d_t_dab_dab(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dxy
   ts(5,6) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dxz
   ts(5,7) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dyz
   ts(5,8) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dx2-y2
   ts(5,9) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = ts(1,6)	! s dxz
   ts(6,2) = -ts(2,6)	! px dxz
   ts(6,3) = -ts(3,6)	! py dxz
   ts(6,4) = -ts(4,6)	! pz dxz
   ts(6,5) = ts(5,6)	! dxy dxz
   ts(6,6) = d_t_dab_dab(l, dl, n, dn, m, dm, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxz dxz
   ts(6,7) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxz dyz
   ts(6,8) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxz dx2-y2
   ts(6,9) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dxz d3z2-r2

   
   ts(7,1) = ts(1,7)	! s dyz
   ts(7,2) = -ts(2,7)	! px dyz
   ts(7,3) = -ts(3,7)	! py dyz
   ts(7,4) = -ts(4,7)	! pz dyz
   ts(7,5) = ts(5,7)	! dxy dyz
   ts(7,6) = ts(6,7)	! dyz dxz
   ts(7,7) = d_t_dab_dab(m, dm, n, dn, l, dl, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dyz dyz
   ts(7,8) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dyz dx2-y2
   ts(7,9) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dyz d3z2-r2

   
   ts(8,1) = ts(1,8)	! s dx2-y2
   ts(8,2) = -ts(2,8)	! px dx2-y2
   ts(8,3) = -ts(3,8)	! py dx2-y2
   ts(8,4) = -ts(4,8)	! pz dx2-y2
   ts(8,5) = ts(5,8)	! dxy dx2-y2
   ts(8,6) = ts(6,8)	! dyz dx2-y2
   ts(8,7) = ts(7,8)	! dyz dx2-y2
   ts(8,8) = d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dx2-y2 dx2-y2
   ts(8,9) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = ts(1,9)	! s d3z2-r2
   ts(9,2) = -ts(2,9)	! px d3z2-r2
   ts(9,3) = -ts(3,9)	! py d3z2-r2
   ts(9,4) = -ts(4,9)	! pz d3z2-r2
   ts(9,5) = ts(5,9)	! dxy d3z2-r2
   ts(9,6) = ts(6,9)	! dyz d3z2-r2
   ts(9,7) = ts(7,9)	! dyz d3z2-r2
   ts(9,8) = ts(8,9)	! dx2-y2 d3z2-r2
   ts(9,9) = d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! d3z2-r2 d3z2-r2
end subroutine d_KS_sp3d5



!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! Different notations for the same subroutines:

pure subroutine KS_sp3d5_TEST(Vsr, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(10), intent(in) :: Vsr		! coefficients for this pair of atoms
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = t_s_s(Vsr(1))		! s s
   ts(1,2) = t_s_px(l, Vsr(2))		! s px
   ts(1,3) = t_s_px(m, Vsr(2))	! s py
   ts(1,4) = t_s_px(n, Vsr(2))	! s pz
   ts(1,5) = t_s_dab(l, m, Vsr(3))	! s dxy
   ts(1,6) = t_s_dab(l, n, Vsr(3))	! s dxz
   ts(1,7) = t_s_dab(m, n, Vsr(3))	! s dyz
   ts(1,8) = t_s_dx2_y2(l, m, Vsr(3))	! s dx2-y2
   ts(1,9) = t_s_dz2_r2(l, m, n, Vsr(3))	! s d3z2-r2
   
   
   ts(2,1) = -ts(1,2)	! s px
   ts(2,2) = t_pa_pa(l, Vsr(4), Vsr(5))		! px px
   ts(2,3) = t_pa_pb(l, m, Vsr(4), Vsr(5))	! px py
   ts(2,4) = t_pa_pb(l, n, Vsr(4), Vsr(5))	! px pz
   ts(2,5) = t_px_dxy(l, m, Vsr(6), Vsr(7))	! px dxy
   ts(2,6) = t_px_dxy(l, n, Vsr(6), Vsr(7))	! px dxz
   ts(2,7) = t_px_dyz(l, m, n, Vsr(6), Vsr(7))	! px dyz
   ts(2,8) = t_px_dx2_y2(l, m, Vsr(6), Vsr(7))	! px dx2-y2
   ts(2,9) = t_pa_d3z2_r2(l, l, m, n, Vsr(6), Vsr(7))	! px d3z2-r2
   
   
   ts(3,1) = -ts(1,3)	! s py
   ts(3,2) = ts(2,3)	! px py
   ts(3,3) = t_pa_pa(m, Vsr(4), Vsr(5))		! py py
   ts(3,4) = t_pa_pb(m, n, Vsr(4), Vsr(5))	! py pz
   ts(3,5) = t_px_dxy(m, l, Vsr(6), Vsr(7))	! py dxy
   ts(3,6) = t_px_dyz(m, l, n, Vsr(6), Vsr(7))	! py dxz
   ts(3,7) = t_px_dxy(m, n,  Vsr(6), Vsr(7))	! py dyz
   ts(3,8) = t_py_dx2_y2(l, m, Vsr(6), Vsr(7))	! py dx2-y2
   ts(3,9) = t_pa_d3z2_r2(m, l, m, n, Vsr(6), Vsr(7))	! py d3z2-r2
   
   
   ts(4,1) = -ts(1,4)	! s pz
   ts(4,2) = ts(2,4)	! px pz
   ts(4,3) = ts(3,4)	! py pz
   ts(4,4) = t_pa_pa(n, Vsr(4), Vsr(5))		! pz pz
   ts(4,5) = t_px_dyz(n, l, m, Vsr(6), Vsr(7))	! pz dxy
   ts(4,6) = t_px_dxy(n, l, Vsr(6), Vsr(7))	! pz dxz
   ts(4,7) = t_px_dxy(n, m, Vsr(6), Vsr(7))	! pz dyz
   ts(4,8) = t_pz_dx2_y2(l, m, n, Vsr(6), Vsr(7))	! pz dx2-y2
   ts(4,9) = t_pz_d3z2_r2(l, m, n, Vsr(6), Vsr(7))	! pz d3z2-r2

   
   ts(5,1) = ts(1,5)	! s dxy
   ts(5,2) = -ts(2,5)	! px dxy
   ts(5,3) = -ts(3,5)	! py dxy
   ts(5,4) = -ts(4,5)	! pz dxy
   ts(5,5) = t_dab_dab(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dxy dxy
   ts(5,6) = t_dab_dbg(m, l, n, Vsr(8), Vsr(9), Vsr(10))	! dxy dxz
   ts(5,7) = t_dab_dbg(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dxy dyz
   ts(5,8) = t_dxy_dx2_y2(l, m, Vsr(8), Vsr(9), Vsr(10))	! dxy dx2-y2
   ts(5,9) = t_dxy_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = ts(1,6)	! s dxz
   ts(6,2) = -ts(2,6)	! px dxz
   ts(6,3) = -ts(3,6)	! py dxz
   ts(6,4) = -ts(4,6)	! pz dxz
   ts(6,5) = ts(5,6)	! dxy dxz
   ts(6,6) = t_dab_dab(l, n, m, Vsr(8), Vsr(9), Vsr(10))	! dxz dxz
   ts(6,7) = t_dab_dbg(l, n, m, Vsr(8), Vsr(9), Vsr(10))	! dxz dyz
   ts(6,8) = t_dxz_dx2_y2(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dxz dx2-y2
   ts(6,9) = t_dxz_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dxz d3z2-r2

   
   ts(7,1) = ts(1,7)	! s dyz
   ts(7,2) = -ts(2,7)	! px dyz
   ts(7,3) = -ts(3,7)	! py dyz
   ts(7,4) = -ts(4,7)	! pz dyz
   ts(7,5) = ts(5,7)	! dxy dyz
   ts(7,6) = ts(6,7)	! dyz dxz
   ts(7,7) = t_dab_dab(m, n, l, Vsr(8), Vsr(9), Vsr(10))	! dyz dyz
   ts(7,8) = t_dyz_dx2_y2(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dyz dx2-y2
   ts(7,9) = t_dyz_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dyz d3z2-r2

   
   ts(8,1) = ts(1,8)	! s dx2-y2
   ts(8,2) = -ts(2,8)	! px dx2-y2
   ts(8,3) = -ts(3,8)	! py dx2-y2
   ts(8,4) = -ts(4,8)	! pz dx2-y2
   ts(8,5) = ts(5,8)	! dxy dx2-y2
   ts(8,6) = ts(6,8)	! dyz dx2-y2
   ts(8,7) = ts(7,8)	! dyz dx2-y2
   ts(8,8) = t_dx2_y2_dx2_y2(l, m, n, Vsr(8), Vsr(9), Vsr(10))	! dx2-y2 dx2-y2
   ts(8,9) = t_dx2_y2_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = ts(1,9)	! s d3z2-r2
   ts(9,2) = -ts(2,9)	! px d3z2-r2
   ts(9,3) = -ts(3,9)	! py d3z2-r2
   ts(9,4) = -ts(4,9)	! pz d3z2-r2
   ts(9,5) = ts(5,9)	! dxy d3z2-r2
   ts(9,6) = ts(6,9)	! dyz d3z2-r2
   ts(9,7) = ts(7,9)	! dyz d3z2-r2
   ts(9,8) = ts(8,9)	! dx2-y2 d3z2-r2
   ts(9,9) = t_d3z2_r2_d3z2_r2(l, m, n, Vsr(8), Vsr(9), Vsr(10)) 	! d3z2-r2 d3z2-r2
end subroutine KS_sp3d5_TEST

pure subroutine d_KS_sp3d5_TEST(Vsr, dVsr, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives
   real(8), dimension(10), intent(in) :: Vsr, dVsr	! coefficients for this pair of atoms and derivatives
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr(1))	! s s
   ts(1,2) = d_t_s_px(l, dl, Vsr(2), dVsr(2))	! s px
   ts(1,3) = d_t_s_px(m, dm, Vsr(2), dVsr(2))	! s py
   ts(1,4) = d_t_s_px(n, dn, Vsr(2), dVsr(2))	! s pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr(3), dVsr(3))	! s dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr(3), dVsr(3))	! s dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr(3), dVsr(3))	! s dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr(3), dVsr(3))	! s dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr(3), dVsr(3))	! s d3z2-r2
   
   
   ts(2,1) = -ts(1,2)	! s px
   ts(2,2) = d_t_pa_pa(l, dl, Vsr(4), dVsr(4), Vsr(5), dVsr(5))		! px px
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr(4), dVsr(4), Vsr(5), dVsr(5))	! px py
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr(4), dVsr(4), Vsr(5), dVsr(5))	! px pz
   ts(2,5) = d_t_px_dxy(l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dxy
   ts(2,6) = d_t_px_dxy(l, dl, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dxz
   ts(2,7) = d_t_px_dyz(l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dyz
   ts(2,8) = d_t_px_dx2_y2(l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px dx2-y2
   ts(2,9) = d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! px d3z2-r2
   
   
   ts(3,1) = -ts(1,3)	! s py
   ts(3,2) = ts(2,3)	! px py
   ts(3,3) = d_t_pa_pa(m, dm, Vsr(4), dVsr(4), Vsr(5), dVsr(5))		! py py
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr(4), dVsr(4), Vsr(5), dVsr(5))	! py pz
   ts(3,5) = d_t_px_dxy(m, dm, l, dl, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dxy
   ts(3,6) = d_t_px_dyz(m, dm, l, dl, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dxz
   ts(3,7) = d_t_px_dxy(m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dyz
   ts(3,8) = d_t_py_dx2_y2(l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py dx2-y2
   ts(3,9) = d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! py d3z2-r2
   
   
   ts(4,1) = -ts(1,4)	! s pz
   ts(4,2) = ts(2,4)	! px pz
   ts(4,3) = ts(3,4)	! py pz
   ts(4,4) = d_t_pa_pa(n, dn, Vsr(4), dVsr(4), Vsr(5), dVsr(5))		! pz pz
   ts(4,5) = d_t_px_dyz(n, dn, l, dl, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dxy
   ts(4,6) = d_t_px_dxy(n, dn, l, dl, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dxz
   ts(4,7) = d_t_px_dxy(n, dn, m, dm, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dyz
   ts(4,8) = d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz dx2-y2
   ts(4,9) = d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr(6), dVsr(6), Vsr(7), dVsr(7))	! pz d3z2-r2

   
   ts(5,1) = ts(1,5)	! s dxy
   ts(5,2) = -ts(2,5)	! px dxy
   ts(5,3) = -ts(3,5)	! py dxy
   ts(5,4) = -ts(4,5)	! pz dxy
   ts(5,5) = d_t_dab_dab(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dxy
   ts(5,6) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dxz
   ts(5,7) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dyz
   ts(5,8) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxy dx2-y2
   ts(5,9) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = ts(1,6)	! s dxz
   ts(6,2) = -ts(2,6)	! px dxz
   ts(6,3) = -ts(3,6)	! py dxz
   ts(6,4) = -ts(4,6)	! pz dxz
   ts(6,5) = ts(5,6)	! dxy dxz
   ts(6,6) = d_t_dab_dab(l, dl, n, dn, m, dm, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxz dxz
   ts(6,7) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxz dyz
   ts(6,8) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dxz dx2-y2
   ts(6,9) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dxz d3z2-r2

   
   ts(7,1) = ts(1,7)	! s dyz
   ts(7,2) = -ts(2,7)	! px dyz
   ts(7,3) = -ts(3,7)	! py dyz
   ts(7,4) = -ts(4,7)	! pz dyz
   ts(7,5) = ts(5,7)	! dxy dyz
   ts(7,6) = ts(6,7)	! dyz dxz
   ts(7,7) = d_t_dab_dab(m, dm, n, dn, l, dl, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dyz dyz
   ts(7,8) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dyz dx2-y2
   ts(7,9) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dyz d3z2-r2

   
   ts(8,1) = ts(1,8)	! s dx2-y2
   ts(8,2) = -ts(2,8)	! px dx2-y2
   ts(8,3) = -ts(3,8)	! py dx2-y2
   ts(8,4) = -ts(4,8)	! pz dx2-y2
   ts(8,5) = ts(5,8)	! dxy dx2-y2
   ts(8,6) = ts(6,8)	! dyz dx2-y2
   ts(8,7) = ts(7,8)	! dyz dx2-y2
   ts(8,8) = d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10))	! dx2-y2 dx2-y2
   ts(8,9) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = ts(1,9)	! s d3z2-r2
   ts(9,2) = -ts(2,9)	! px d3z2-r2
   ts(9,3) = -ts(3,9)	! py d3z2-r2
   ts(9,4) = -ts(4,9)	! pz d3z2-r2
   ts(9,5) = ts(5,9)	! dxy d3z2-r2
   ts(9,6) = ts(6,9)	! dyz d3z2-r2
   ts(9,7) = ts(7,9)	! dyz d3z2-r2
   ts(9,8) = ts(8,9)	! dx2-y2 d3z2-r2
   ts(9,9) = d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vsr(8), dVsr(8), Vsr(9), dVsr(9), Vsr(10), dVsr(10)) 	! d3z2-r2 d3z2-r2
end subroutine d_KS_sp3d5_TEST


pure subroutine KS_sp3_hetero_TEST(Vsr12, Vsr21, l, m, n, ts)
   real(8), intent(in) :: l, m, n		! direction cosines
   real(8), dimension(4), intent(in) :: Vsr12, Vsr21    ! two sets of coefficients for this pair of atoms
   real(8), dimension(4,4), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)

   ! functions above:
   ts(1,1) = t_s_s(Vsr12(1))
   ts(1,2) = t_s_px(l, Vsr12(2))
   ts(1,3) = t_s_px(m, Vsr12(2))
   ts(1,4) = t_s_px(n, Vsr12(2))
   
   ts(2,1) = -t_s_px(l, Vsr21(2))
   ts(2,2) = t_pa_pa(l, Vsr12(3), Vsr12(4))
   ts(2,3) = t_pa_pb(l, m, Vsr12(3), Vsr12(4))
   ts(2,4) = t_pa_pb(l, n, Vsr12(3), Vsr12(4))
   
   ts(3,1) = -t_s_px(m, Vsr21(2))
   ts(3,2) = t_pa_pb(l, m, Vsr21(3), Vsr21(4))
   ts(3,3) = t_pa_pa(m, Vsr12(3), Vsr12(4))
   ts(3,4) = t_pa_pb(m, n, Vsr12(3), Vsr12(4))
   
   ts(4,1) = -t_s_px(n, Vsr21(2))
   ts(4,2) = t_pa_pb(l, n, Vsr21(3), Vsr21(4))
   ts(4,3) = t_pa_pb(m, n, Vsr21(3), Vsr21(4))
   ts(4,4) = t_pa_pa(n, Vsr12(3), Vsr12(4))
end subroutine KS_sp3_hetero_TEST


! Constructing Koster Slater part for sp3 basis set:
pure subroutine d_KS_sp3_hetero_TEST(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives by r_{ij}
   real(8), dimension(4), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21	! coefficients and derivatives by r_{ij} for this pair of atoms
   real(8), dimension(4,4), intent(out) :: ts	! constructed hopping integrals
   !    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (p p sigma)
!    V(4) = (p p pi)

   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))
   ts(1,2) = d_t_s_px(l, dl, Vsr12(2), dVsr12(2))
   ts(1,3) = d_t_s_px(m, dm, Vsr12(2), dVsr12(2))
   ts(1,4) = d_t_s_px(n, dn, Vsr12(2), dVsr12(2))
   
   ts(2,1) = -d_t_s_px(l, dl, Vsr21(2), dVsr21(2))
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   
   ts(3,1) = -d_t_s_px(m, dm, Vsr21(2), dVsr21(2))
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
   
   ts(4,1) = -d_t_s_px(n, dn, Vsr21(2), dVsr21(2))
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(3), dVsr21(3), Vsr21(4), dVsr21(4))
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(3), dVsr12(3), Vsr12(4), dVsr12(4))
end subroutine d_KS_sp3_hetero_TEST


pure subroutine d_KS_sp3d5_hetero_TEST(Vsr12, Vsr21, dVsr12, dVsr21, l, m, n, dl, dm, dn, ts)
   real(8), intent(in) :: l, m, n	, dl, dm, dn	! direction cosines and derivatives
   real(8), dimension(10), intent(in) :: Vsr12, Vsr21, dVsr12, dVsr21	! coefficients for this pair of atoms and derivatives
   real(8), dimension(9,9), intent(out) :: ts	! constructed hopping integrals
!    Reminder:
!    V(1) = (s s sigma)
!    V(2) = (s p sigma)
!    V(3) = (s d sigma)
!    V(4) = (p p sigma)
!    V(5) = (p p pi)
!    V(6) = (p d sigma)
!    V(7) = (p d pi)
!    V(8) = (d d sigma)
!    V(9) = (d d pi)
!    V(10) = (d d delta)
   
   ! functions above:
   ts(1,1) = d_t_s_s(dVsr12(1))	! s s
   ts(1,2) = d_t_s_px(l, dl, Vsr12(2), dVsr12(2))	! s px
   ts(1,3) = d_t_s_px(m, dm, Vsr12(2), dVsr12(2))	! s py
   ts(1,4) = d_t_s_px(n, dn, Vsr12(2), dVsr12(2))	! s pz
   ts(1,5) = d_t_s_dab(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dxy
   ts(1,6) = d_t_s_dab(l, dl, n, dn, Vsr12(3), dVsr12(3))	! s dxz
   ts(1,7) = d_t_s_dab(m, dm, n, dn, Vsr12(3), dVsr12(3))	! s dyz
   ts(1,8) = d_t_s_dx2_y2(l, dl, m, dm, Vsr12(3), dVsr12(3))	! s dx2-y2
   ts(1,9) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr12(3), dVsr12(3))	! s d3z2-r2
   
   
   ts(2,1) = -d_t_s_px(l, dl, Vsr21(2), dVsr21(2))	! px s 
   ts(2,2) = d_t_pa_pa(l, dl, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! px px
   ts(2,3) = d_t_pa_pb(l, dl, m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px py
   ts(2,4) = d_t_pa_pb(l, dl, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! px pz
   ts(2,5) = d_t_px_dxy(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxy
   ts(2,6) = d_t_px_dxy(l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dxz
   ts(2,7) = d_t_px_dyz(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dyz
   ts(2,8) = d_t_px_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px dx2-y2
   ts(2,9) = d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! px d3z2-r2
   
   
   ts(3,1) = -d_t_s_px(m, dm, Vsr21(2), dVsr21(2))	! py s 
   ts(3,2) = d_t_pa_pb(l, dl, m, dm, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! py px 
   ts(3,3) = d_t_pa_pa(m, dm, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! py py
   ts(3,4) = d_t_pa_pb(m, dm, n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))	! py pz
   ts(3,5) = d_t_px_dxy(m, dm, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxy
   ts(3,6) = d_t_px_dyz(m, dm, l, dl, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dxz
   ts(3,7) = d_t_px_dxy(m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dyz
   ts(3,8) = d_t_py_dx2_y2(l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py dx2-y2
   ts(3,9) = d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! py d3z2-r2
   
   
   ts(4,1) = -d_t_s_px(n, dn, Vsr21(2), dVsr21(2))	! pz s 
   ts(4,2) = d_t_pa_pb(l, dl, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz px 
   ts(4,3) = d_t_pa_pb(m, dm, n, dn, Vsr21(4), dVsr21(4), Vsr21(5), dVsr21(5))	! pz py 
   ts(4,4) = d_t_pa_pa(n, dn, Vsr12(4), dVsr12(4), Vsr12(5), dVsr12(5))		! pz pz
   ts(4,5) = d_t_px_dyz(n, dn, l, dl, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxy
   ts(4,6) = d_t_px_dxy(n, dn, l, dl, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dxz
   ts(4,7) = d_t_px_dxy(n, dn, m, dm, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dyz
   ts(4,8) = d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz dx2-y2
   ts(4,9) = d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(6), dVsr12(6), Vsr12(7), dVsr12(7))	! pz d3z2-r2

   
   ts(5,1) = d_t_s_dab(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dxy s 
   ts(5,2) = -d_t_px_dxy(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy px 
   ts(5,3) = -d_t_px_dxy(m, dm, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy py 
   ts(5,4) = -d_t_px_dyz(n, dn, l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxy pz 
   ts(5,5) = d_t_dab_dab(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxy
   ts(5,6) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dxz
   ts(5,7) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dyz
   ts(5,8) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxy dx2-y2
   ts(5,9) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxy d3z2-r2
   
   
   ts(6,1) = d_t_s_dab(l, dl, n, dn, Vsr21(3), dVsr21(3))	! dxz s 
   ts(6,2) = -d_t_px_dxy(l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz px 
   ts(6,3) = -d_t_px_dyz(m, dm, l, dl, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz py 
   ts(6,4) = -d_t_px_dxy(n, dn, l, dl, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dxz pz 
   ts(6,5) = d_t_dab_dbg(m, dm, l, dl, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dxz dxy 
   ts(6,6) = d_t_dab_dab(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dxz
   ts(6,7) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dyz
   ts(6,8) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dxz dx2-y2
   ts(6,9) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dxz d3z2-r2

   
   ts(7,1) = d_t_s_dab(m, dm, n, dn, Vsr21(3), dVsr21(3))	! dyz s 
   ts(7,2) = -d_t_px_dyz(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz px 
   ts(7,3) = -d_t_px_dxy(m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz py 
   ts(7,4) = -d_t_px_dxy(n, dn, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dyz pz 
   ts(7,5) = d_t_dab_dbg(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxy 
   ts(7,6) = d_t_dab_dbg(l, dl, n, dn, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dyz dxz 
   ts(7,7) = d_t_dab_dab(m, dm, n, dn, l, dl, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dyz
   ts(7,8) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dyz dx2-y2
   ts(7,9) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dyz d3z2-r2

   
   ts(8,1) = d_t_s_dx2_y2(l, dl, m, dm, Vsr21(3), dVsr21(3))	! dx2-y2 s 
   ts(8,2) = -d_t_px_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 px 
   ts(8,3) = -d_t_py_dx2_y2(l, dl, m, dm, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 py 
   ts(8,4) = -d_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! dx2-y2 pz 
   ts(8,5) = d_t_dxy_dx2_y2(l, dl, m, dm, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxy 
   ts(8,6) = d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dxz 
   ts(8,7) = d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10))	! dx2-y2 dyz 
   ts(8,8) = d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! dx2-y2 d3z2-r2

   
   ts(9,1) = d_t_s_dz2_r2(l, dl, m, dm, n, dn, Vsr21(3), dVsr21(3))	! d3z2-r2 s 
   ts(9,2) = -d_t_pa_d3z2_r2(l, dl, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 px 
   ts(9,3) = -d_t_pa_d3z2_r2(m, dm, l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 py 
   ts(9,4) = -d_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(6), dVsr21(6), Vsr21(7), dVsr21(7))	! d3z2-r2 pz 
   ts(9,5) = d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxy 
   ts(9,6) = d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dxz 
   ts(9,7) = d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dyz 
   ts(9,8) = d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vsr21(8), dVsr21(8), Vsr21(9), dVsr21(9), Vsr21(10), dVsr21(10)) 	! d3z2-r2 dx2-y2 
   ts(9,9) = d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vsr12(8), dVsr12(8), Vsr12(9), dVsr12(9), Vsr12(10), dVsr12(10)) 	! d3z2-r2 d3z2-r2
end subroutine d_KS_sp3d5_hetero_TEST

END MODULE TB_Koster_Slater
