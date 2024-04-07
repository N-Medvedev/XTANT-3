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
! This module contains Koster-Slater integrals for all shells:
! References used:
! [1] P. Koskinen, V. Mäkinen, Computational Materials Science 47, 237-253 (2009)
! [2] A.V. Podolskiy and P. Vogl, Phys. Rev. B. 69, 233101 (2004)
! =============================================================

MODULE TB_Koster_Slater

use Algebra_tools, only: Kronecker_delta, Heavyside_tau, get_factorial

implicit none
PRIVATE

! Modular parameters:
real(8) :: m_sqrt_inv_2, m_sqrt3, m_sqrt3_half

parameter (m_sqrt_inv_2 = 1.0d0 / sqrt(2.0d0))
parameter (m_sqrt3 = sqrt(3.0d0))
parameter (m_sqrt3_half = m_sqrt3*0.5d0)


public :: KS_s, KS_ss_hetero, KS_sp3_hetero, KS_sp3d5_hetero, drij_dska, ddija_dskb_kd, &
         d_KS_s, d_KS_sp3_hetero, d_KS_sp3d5_hetero, dda_dhgd, drij_dhab, drij_drka, &
         ddija_drkb, d2dija_drkb2, m_sqrt3, t_s_dab, t_s_dx2_y2, t_s_dz2_r2, t_s_s, t_pa_pa, &
         t_dab_dab, t_dx2_y2_dx2_y2, t_d3z2_r2_d3z2_r2, KS_sp3_hetero_TEST, KS_sp3d5_hetero_TEST, &
         KS_sp3d5, d_KS_sp3d5, d_KS_sp3d5_TEST


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
   dab = Kronecker_delta(alpha, beta)	! module "Algebra_tools"
   ddija_drkb = (dik - djk)/rij*(dab - rija*rijb/(rij*rij))
end function ddija_drkb


! Second derivatives of the direction cosine by a coordinate of a third  nad fourth atom:
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
   real(8) :: dija, dijb
   ! Delta-symbols:
   dik = Kronecker_delta(i,k) 	! module "Algebra_tools"
   djk = Kronecker_delta(j,k) 	! module "Algebra_tools"
   dab = Kronecker_delta(alpha, beta)	! module "Algebra_tools"

   dija = rija / rij
   dijb = rijb / rij
   rij2 = rij*rij
   !rij3 = rij*rij2

   d2dija_drkb2 = -(dik - djk)/rij2 * ( 2.0d0*(dab - dija*dijb)*dijb + dija*(1.0d0 - dijb*dijb) )
   !d2dija_drkb2 = -(dik - djk)/rij3 * ( 2.0d0*rijb*(dab - rija*rijb/rij2) + rija*(1.0d0 - rijb*rijb/rij2) )   ! OLD
end function d2dija_drkb2


!---------------------------------------------------------------------
! Koster-Slater rotational functions for arbitrary (n,l,m), following
! [A.V. Podolskiy and P. Vogl, Phys. Rev. B. 69, 233101 (2004)]

! Computes rotational matrix elements according to Eq.(7):
pure function coefs_dlmm(N, l, m1, m2) result(dlmm)
   real(8) :: dlmm
   real(8), intent(in) :: N ! directional cosine relative to z-axis
   integer, intent(in) :: l, m1, m2 ! orbital and magnetic quantum numbers
   !---------------------
   integer :: t
   real(8) :: prefN, prefR, summ
   real(8) :: fact_1, fact_2, fact_3, fact_4

   ! Prefactors with factorials:
   if (abs(N) < 1.0d0) then
      prefN = (0.5d0 + 0.5d0*N)**l * ((1.0d0 - N)/(1.0d0 + N))**(0.5d0*(m1 - m2))

      ! Get factorials:
      fact_1 = get_factorial(l + m2)  ! module "Algebra_tools"
      fact_2 = get_factorial(l - m2)  ! module "Algebra_tools"
      fact_3 = get_factorial(l + m1)  ! module "Algebra_tools"
      fact_4 = get_factorial(l - m1)  ! module "Algebra_tools"
      prefR = sqrt(fact_1 * fact_2 * fact_3 * fact_4)

      summ = 0.0d0   ! to start with
      do t = 0, (2*l + 1)
         if ((0 <= l + m2 - t) .and. (0 <= l - m1 - t) .and. (0 <= t + m1 - m2)) then
            fact_1 = get_factorial(l + m2 - t)  ! module "Algebra_tools"
            fact_2 = get_factorial(l - m1 - t)  ! module "Algebra_tools"
            fact_3 = get_factorial(t)  ! module "Algebra_tools"
            fact_4 = get_factorial(t + m1 - m2)  ! module "Algebra_tools"
            summ = summ + (-1.0d0)**t*((1.0d0 - N)/(1.0d0 + N))**t / (fact_1*fact_2*fact_3*fact_4)
         endif
      enddo
   else ! along Z, no contribution
      prefN = 1.0d0
      prefR = 1.0d0
      summ = 0.0d0
   endif

   ! Collect the terms:
   dlmm = prefR * prefN * summ
end function coefs_dlmm


! Coefficients of S and T from Eqs.(24,25):
pure subroutine coefs_S_T()
   ! UNFINISHED
end subroutine coefs_S_T


! Both Am abd Bm coefficients from Eq.(21):
pure subroutine coefs_Am_Bm(m, gam, Am, Bm)
   integer, intent(in) :: m
   real(8), intent(in) :: gam
   real(8), intent(out) :: Am, Bm  ! Eq.(21)
   !----------------------
   real(8) :: abs_m, tau_p, tau_m, prefac, sin_m, cos_m
   if (m == 0) then
      Am = m_sqrt_inv_2
      Bm = 0.0d0  ! added for convenience of Eq.(23-25)
   else
      abs_m = abs(m)
      tau_p = Heavyside_tau(m)   ! module "Algebra_tools"
      tau_m = Heavyside_tau(-m)  ! module "Algebra_tools"
      prefac = (-1.0d0)**abs_m
      sin_m = sin(abs_m * gam)
      cos_m = cos(abs_m * gam)
      Am = prefac * (tau_p * cos_m - tau_m * sin_m)
      Bm = prefac * (tau_p * sin_m + tau_m * cos_m)
   endif
end subroutine coefs_Am_Bm


! Am and Bm as individual functions:
pure function coef_Am(m, gam) result(Am)
   real(8) Am  ! Eq.(21)
   integer, intent(in) :: m
   real(8), intent(in) :: gam
   !---------------
   real(8) :: abs_m, tau_p, tau_m
   if (m == 0) then
      Am = m_sqrt_inv_2
   else
      abs_m = abs(m)
      tau_p = Heavyside_tau(m)  ! module "Algebra_tools"
      tau_m = Heavyside_tau(-m)  ! module "Algebra_tools"
      Am = (-1.0d0)**abs_m * (tau_p * cos(abs_m * gam) - tau_m * sin(abs_m * gam))
   endif
end function coef_Am


pure function coef_Bm(m, gam) result(Bm)
   real(8) Bm  ! Eq.(21), with added case Bm(m=0)
   integer, intent(in) :: m
   real(8), intent(in) :: gam
   !---------------
   real(8) :: abs_m, tau_p, tau_m
   if (m == 0) then
      Bm = 0.0d0  ! added for convenience of Eq.(23-25)
   else
      abs_m = abs(m)
      tau_p = Heavyside_tau(m)  ! module "Algebra_tools"
      tau_m = Heavyside_tau(-m)  ! module "Algebra_tools"
      Bm = (-1.0d0)**abs_m * (tau_p * sin(abs_m * gam) + tau_m * cos(abs_m * gam))
   endif
end function coef_Bm


!---------------------------------------------------------------------
! Individual Koster-Slater hopping integrals as individual functions:
! Following the original work [Slater and Koster, PRB 94, 1498 (1954)]

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
pure function t_dab_dab(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: l2, m2, l2m2
   l2 = l*l
   m2 = m*m
   l2m2 = l2*m2
   Ecc = 3.0d0*l2m2*Vdd_sigma + (l2 + m2 - 4.0d0*l2m2)*Vdd_pi + (n*n + l2m2)*Vdd_delta
end function t_dab_dab


! 16)
pure function t_dab_dbg(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: ln, m2
   ln = l*n
   m2 = m*m
   Ecc = ( 3.0d0*m2*Vdd_sigma + (1.0d0 - 4.0d0*m2)*Vdd_pi + (m2 - 1.0d0)*Vdd_delta )*ln
end function t_dab_dbg


! 17)
pure function t_dxy_dx2_y2(l, m, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: lm, l2m2
   lm = l*m
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*Vdd_sigma - 2.0d0*Vdd_pi + 0.5d0*Vdd_delta )*lm*l2m2
end function t_dxy_dx2_y2


! 18)
pure function t_dyz_dx2_y2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: mn, l2m2
   mn = m*n
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*l2m2*Vdd_sigma - (1.0d0+2.0d0*l2m2)*Vdd_pi + (1.0d0+0.5d0*l2m2)*Vdd_delta )*mn
end function t_dyz_dx2_y2


! 19)
pure function t_dxz_dx2_y2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: ln, l2m2
   ln = l*n
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*l2m2*Vdd_sigma + (1.0d0-2.0d0*l2m2)*Vdd_pi - (1.0d0-0.5d0*l2m2)*Vdd_delta )*ln
end function t_dxz_dx2_y2


! 20)
pure function t_dxy_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: lm, l2m2, n2
   lm = l*m
   l2m2 = l*l + m*m
   n2 = n*n
   Ecc = ( m_sqrt3*( (n2 - 0.5d0*l2m2)*Vdd_sigma - 2.0d0*n2*Vdd_pi) + m_sqrt3_half*(1.0d0+n2)*Vdd_delta )*lm
end function t_dxy_d3z2_r2


! 21)
pure function t_dyz_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: mn, l2m2, n2
   mn = m*n
   l2m2 = l*l + m*m
   n2 = n*n
   Ecc = ( (n2 - 0.5d0*l2m2)*Vdd_sigma + (l2m2 - n2)*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*m_sqrt3*mn
end function t_dyz_d3z2_r2


! 22)
pure function t_dxz_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: ln, l2m2, n2
   ln = l*n
   l2m2 = l*l + m*m
   n2 = n*n
   Ecc = ( (n2 - 0.5d0*l2m2)*Vdd_sigma + (l2m2 - n2)*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*m_sqrt3*ln
end function t_dxz_d3z2_r2


! 23)
pure function t_dx2_y2_dx2_y2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
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
pure function t_dx2_y2_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta
   real(8) :: l2m2, l2_m2, n2
   l2m2 = l*l + m*m
   l2_m2 = l*l - m*m
   n2 = n*n
   Ecc = ( 0.5d0*(n2 - 0.5d0*l2m2)*Vdd_sigma - n2*Vdd_pi + 0.25d0*(1.0d0+n2)*Vdd_delta )*m_sqrt3*l2_m2
end function t_dx2_y2_d3z2_r2


! 25)
pure function t_d3z2_r2_d3z2_r2(l, m, n, Vdd_sigma, Vdd_pi, Vdd_delta) result (Ecc)
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
pure function d_t_dab_dab(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dab_dbg(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dxy_dx2_y2(l, dl, m, dm, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: lm, l2m2
   lm = l*m
   l2m2 = l*l - m*m
   Ecc = ( 1.50d0*dVdd_sigma - 2.0d0*dVdd_pi + 0.5d0*dVdd_delta )*lm*l2m2
   Ecc = Ecc + ( 1.50d0*Vdd_sigma - 2.0d0*Vdd_pi + 0.5d0*Vdd_delta )*((dl*m+l*dm)*l2m2 + 2.0d0*lm*(dl*l-dm*m))
end function d_t_dxy_dx2_y2


! 18)
pure function d_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
pure function d_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
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
! Second derivatives of the Koster-Slater angular functions.
! This functions do NOT include the derivatives of dr/ds, or dr/dh, or dr/dr_k,
! meaning, they MUST be already included in the coefficients passed to the functions:
! dV = dV/dr_ij * dr_ij/dr_k,g


! Preporatories:
! first derivative dV/dr_k,g = dV/dr_ij * dr_ij/dr_k,g:
pure function d_V(dVss_sigma, rij, rijg, i, j, k) result(dV)
   real(8) :: dV
   real(8), intent(in) :: dVss_sigma   ! dV/dr_ij (without dr_ij/dr_k,g)
   real(8), intent(in) :: rij, rijg
   integer, intent(in) :: i, j, k   ! atom indices
   !---------------------
   real(8) :: da, delt, dik, djk

   ! Delta-symbols:
   dik = Kronecker_delta(i,k) ! module "Algebra_tools"
   djk = Kronecker_delta(j,k) ! module "Algebra_tools"
   delt = dik - djk
   ! d rij / d r_ij,g:
   da = delt * rijg / rij
   ! dV / dr_{k,g} = dV/dr_ij * dt_ij/dr_{k,g}:
   dV = dVss_sigma*da
end function d_V


! second derivative d^2V/dr^2_k,g = d^2V/dr^2_ij * (dr_ij/dr_k,g)^2 + dV/dr_ij * d^2r_ij/dr^2_k,g
pure function d2_V(dVss_sigma, d2Vss_sigma, rij, rijg, i, j, k) result(d2V)
   real(8) :: d2V
   real(8), intent(in) :: dVss_sigma, d2Vss_sigma, rij, rijg
   integer, intent(in) :: i, j, k   ! atom indices
   !---------------------
   real(8) :: da, dda, delt, dik, djk

   ! Delta-symbols:
   dik = Kronecker_delta(i,k) ! module "Algebra_tools"
   djk = Kronecker_delta(j,k) ! module "Algebra_tools"
   delt = dik - djk
   ! d rij / d r_ij,g:
   da = delt * rijg / rij
   ! d^2r_ij/dr^2_k,g:
   dda = delt * (1.0d0 - da**2)
   ! d^2V/dr^2_k,g:
   d2V = d2Vss_sigma*da**2 + dVss_sigma*dda
end function d2_V


! 1)
pure function d2_t_s_s(d2Vss_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: d2Vss_sigma  ! d^2 V / d r^2_k,g
   Ecc = d2Vss_sigma
end function d2_t_s_s


! 2)
pure function d2_t_s_px(da, dda, d2da, Vsp_sigma, dVsp_sigma, d2Vss_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, d2da, Vsp_sigma, dVsp_sigma, d2Vss_sigma
   Ecc = d2da*Vsp_sigma + 2.0d0*da*dVsp_sigma + da*d2Vss_sigma
end function d2_t_s_px


! 3)
pure function d2_t_pa_pa(da, dda, d2da, Vpp_sigma, dVpp_sigma, d2Vpp_sigma, Vpp_pi, dVpp_pi, d2Vpp_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, d2da, Vpp_sigma, dVpp_sigma, d2Vpp_sigma, Vpp_pi, dVpp_pi, d2Vpp_pi
   real(8) :: da2, dda2
   da2 = da**2
   dda2 = dda**2
   Ecc = 2.0d0*(Vpp_sigma - Vpp_pi)*(dda2 + da*d2da) + &
         4.0d0*(dVpp_sigma - dVpp_pi)*da*dda + &
         da2*d2Vpp_sigma + (1.0d0 - da2)*d2Vpp_pi
end function d2_t_pa_pa


! 4)
pure function d2_t_pa_pb(da, dda, d2da, db, ddb, d2db, Vpp_sigma, dVpp_sigma, d2Vpp_sigma, Vpp_pi, dVpp_pi, d2Vpp_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, d2da, db, ddb, d2db, Vpp_sigma, dVpp_sigma, d2Vpp_sigma, Vpp_pi, dVpp_pi, d2Vpp_pi
   Ecc = (d2da*db + d2db*da + 2.0d0*dda*ddb)*(Vpp_sigma - Vpp_pi) + &
         2.0d0*(dda*db+da*ddb)*(dVpp_sigma - dVpp_pi) + &
         da*db*(d2Vpp_sigma - d2Vpp_pi)
end function d2_t_pa_pb


! 5)
pure function d2_t_s_dab(da, dda, d2da, db, ddb, d2db, Vsd_sigma, dVsd_sigma, d2Vsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, d2da, db, ddb, d2db, Vsd_sigma, dVsd_sigma, d2Vsd_sigma

   Ecc = (d2da*db + da*d2db + 2.0d0*dda*ddb)*Vsd_sigma + &
         2.0d0*(dda*db + da*ddb)*dVsd_sigma + &
         da*db*d2Vsd_sigma
   Ecc = Ecc * m_sqrt3
end function d2_t_s_dab


! 6)
pure function d2_t_s_dx2_y2(da, dda, d2da, db, ddb, d2db, Vsd_sigma, dVsd_sigma, d2Vsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, d2da, db, ddb, d2db, Vsd_sigma, dVsd_sigma, d2Vsd_sigma
   !---------------------
   Ecc = 2.0d0*(d2da*da + dda**2 - ddb**2 - d2db*db)*Vsd_sigma + &
         4.0d0*(da*dda - db*ddb)*dVsd_sigma + &
         (da*da - db*db)*d2Vsd_sigma
   Ecc = m_sqrt3_half*Ecc
end function d2_t_s_dx2_y2


! 7)
pure function d2_t_s_dz2_r2(l, dl, d2l, m, dm, d2m, n, dn, d2n, Vsd_sigma, dVsd_sigma, d2Vsd_sigma) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, d2l, m, dm, d2m, n, dn, d2n, Vsd_sigma, dVsd_sigma, d2Vsd_sigma
   Ecc = ( 2.0d0*(d2n*n + dn**2) - (d2l*l + dl**2 + d2m*m + dm**2) )*Vsd_sigma + &
         2.0d0*( 2.0d0*(dn*n) - (dl*l + dm*m) )*dVsd_sigma + &
         ( n*n - 0.5d0*(l*l + m*m) )*d2Vsd_sigma
end function d2_t_s_dz2_r2


! 8)
pure function d2_t_px_dxy(l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2, temp, dll
   ! NOT READY
   Ecc = 0.0d0
!    l2 = l*l
!    dll = 2.0d0*dl*l
!    temp = (1.0d0 - 2.0d0*l2)
!    Ecc = (m_sqrt3*l2*Vpd_sigma + temp*Vpd_pi)*dm
!    Ecc = Ecc  + (m_sqrt3*(dll*Vpd_sigma + l2*dVpd_sigma) - 2.0d0*dll*Vpd_pi + temp*dVpd_pi)*m
end function d2_t_px_dxy


! 9)
pure function d2_t_px_dyz(l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: mn
   ! NOT READY
   Ecc = 0.0d0
!    mn = m*n
!    Ecc = (m_sqrt3*dVpd_sigma - 2.0d0*dVpd_pi)*l*mn
!    Ecc = Ecc + (m_sqrt3*Vpd_sigma - 2.0d0*Vpd_pi)*(dl*mn + l*dm*n + l*m*dn)
end function d2_t_px_dyz


! 10)
pure function d2_t_px_dx2_y2(l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2, temp
   ! NOT READY
   Ecc = 0.0d0
!    l2_m2 = l*l - m*m
!    dl2_m2 = 2.0d0*( dl*l - dm*m )
!    temp = (m_sqrt3_half*Vpd_sigma - Vpd_pi)
!    Ecc = ( temp*l2_m2 + Vpd_pi )*dl
!    Ecc = Ecc + ( (m_sqrt3_half*dVpd_sigma - dVpd_pi)*l2_m2 + temp*dl2_m2 + dVpd_pi )*l
end function d2_t_px_dx2_y2


! 11)
pure function d2_t_py_dx2_y2(l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2, temp
   ! NOT READY
   Ecc = 0.0d0
!    l2_m2 = l*l - m*m
!    dl2_m2 = 2.0d0*( dl*l - dm*m )
!    temp = (m_sqrt3_half*Vpd_sigma - Vpd_pi)
!    Ecc = ( temp*l2_m2 - Vpd_pi )*dm
!    Ecc = Ecc + ( (m_sqrt3_half*dVpd_sigma - dVpd_pi)*l2_m2 + temp*dl2_m2 - dVpd_pi )*m
end function d2_t_py_dx2_y2


! 12)
pure function d2_t_pz_dx2_y2(l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2
   ! NOT READY
   Ecc = 0.0d0
!    l2_m2 = l*l - m*m
!    dl2_m2 = 2.0d0*(dl*l - dm*m)
!    Ecc = (m_sqrt3_half*dVpd_sigma - dVpd_pi)*l2_m2*n
!    Ecc = Ecc + (m_sqrt3_half*Vpd_sigma - Vpd_pi)*(dl2_m2*n + l2_m2*dn)
end function d2_t_pz_dx2_y2


! 13)
pure function d2_t_pa_d3z2_r2(da, dda, l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc) ! only px or py
   real(8) :: Ecc
   real(8), intent(in) :: da, dda, l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: n2, temp, dn2
   ! NOT READY
   Ecc = 0.0d0
!    n2 = n*n
!    dn2 = 2.0d0*dn*n
!    temp = (n2 - 0.5d0*(l*l + m*m))
!    Ecc = ( temp*Vpd_sigma - m_sqrt3*n2*Vpd_pi )*dda
!    Ecc = Ecc + ( (dn2 - (dl*l + dm*m))*Vpd_sigma + temp*dVpd_sigma - m_sqrt3*(dn2*Vpd_pi + n2*dVpd_pi) )*da
end function d2_t_pa_d3z2_r2


! 14)
pure function d2_t_pz_d3z2_r2(l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi) result (Ecc) ! only pz
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vpd_sigma, dVpd_sigma, Vpd_pi, dVpd_pi
   real(8) :: l2_m2, dl2_m2, n2, temp
   ! NOT READY
   Ecc = 0.0d0
!    l2_m2 = l*l + m*m
!    dl2_m2 = 2.0d0*(dl*l + dm*m)
!    n2 = n*n
!    temp = (n2 - 0.5d0*l2_m2)
!    Ecc = ( temp*Vpd_sigma + m_sqrt3*l2_m2*Vpd_pi )*dn
!    Ecc = Ecc  + ( (2.0d0*dn*n - 0.5d0*dl2_m2)*Vpd_sigma + temp*dVpd_sigma + m_sqrt3*(dl2_m2*Vpd_pi + l2_m2*dVpd_pi) )*n
end function d2_t_pz_d3z2_r2


! 15)
pure function d2_t_dab_dab(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2, m2, l2m2, dl2m2
   ! NOT READY
   Ecc = 0.0d0
!    l2 = l*l
!    m2 = m*m
!    l2m2 = l2*m2
!    dl2m2 = 2.0d0*(dl*l*m2 + l2*m*dm)
!    Ecc = 3.0d0*(dl2m2*Vdd_sigma + l2m2*dVdd_sigma)
!    Ecc = Ecc + 2.0d0*(l*dl + m*dm - 2.0d0*dl2m2)*Vdd_pi + (l2 + m2 - 4.0d0*l2m2)*dVdd_pi
!    Ecc = Ecc + (2.0d0*n*dn + dl2m2)*Vdd_delta +  (n*n + l2m2)*dVdd_delta
end function d2_t_dab_dab


! 16)
pure function d2_t_dab_dbg(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: ln, m2, dm2, temp, temp1
   ! NOT READY
   Ecc = 0.0d0
!    ln = l*n
!    m2 = m*m
!    dm2 = 2.0d0*m*dm
!    temp = (1.0d0 - 4.0d0*m2)
!    temp1 = (m2 - 1.0d0)
!    Ecc = 3.0d0*(dm2*Vdd_sigma + m2*dVdd_sigma)
!    Ecc = Ecc + temp*dVdd_pi - 4.0d0*dm2*Vdd_pi
!    Ecc = Ecc +  temp1*dVdd_delta + dm2*Vdd_delta
!    Ecc = Ecc*ln
!    Ecc = Ecc + ( 3.0d0*m2*Vdd_sigma + temp*Vdd_pi + temp1*Vdd_delta )*(dl*n+l*dn)
end function d2_t_dab_dbg


! 17)
pure function d2_t_dxy_dx2_y2(l, dl, m, dm, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: lm, l2m2
   ! NOT READY
   Ecc = 0.0d0
!    lm = l*m
!    l2m2 = l*l - m*m
!    Ecc = ( 1.50d0*dVdd_sigma - 2.0d0*dVdd_pi + 0.5d0*dVdd_delta )*lm*l2m2
!    Ecc = Ecc + ( 1.50d0*Vdd_sigma - 2.0d0*Vdd_pi + 0.5d0*Vdd_delta )*((dl*m+l*dm)*l2m2 + 2.0d0*lm*(dl*l-dm*m))
end function d2_t_dxy_dx2_y2


! 18)
pure function d2_t_dyz_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: mn, l2m2, dl2m2, temp, temp1
   ! NOT READY
   Ecc = 0.0d0
!    mn = m*n
!    l2m2 = l*l - m*m
!    dl2m2 = 2.0d0*(dl*l - dm*m)
!    temp = (1.0d0+2.0d0*l2m2)
!    temp1 = (1.0d0+0.5d0*l2m2)
!    Ecc = ( 1.50d0*(dl2m2*Vdd_sigma+l2m2*dVdd_sigma) - (2.0d0*dl2m2*Vdd_pi + temp*dVdd_pi) + temp1*dVdd_delta + 0.5d0*dl2m2*Vdd_delta )*mn
!    Ecc = Ecc + ( 1.50d0*l2m2*Vdd_sigma - temp*Vdd_pi + temp1*Vdd_delta )*(dm*n+m*dn)
end function d2_t_dyz_dx2_y2


! 19)
pure function d2_t_dxz_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: ln, l2m2, dl2m2, temp, temp1
   ! NOT READY
   Ecc = 0.0d0
!    ln = l*n
!    l2m2 = l*l - n*n
!    dl2m2 = 2.0d0*(dl*l - dm*m)
!    temp = (1.0d0-2.0d0*l2m2)
!    temp1 = (1.0d0-0.5d0*l2m2)
!    Ecc = ( 1.50d0*(dl2m2*Vdd_sigma+l2m2*dVdd_sigma) + (-2.0d0*dl2m2*Vdd_pi + temp*dVdd_pi) - (temp1*dVdd_delta - 0.5d0*dl2m2*Vdd_delta) )*ln
!    Ecc = Ecc + ( 1.50d0*l2m2*Vdd_sigma + temp*Vdd_pi + temp1*Vdd_delta )*(dl*n+l*dn)
end function d2_t_dxz_dx2_y2


! 20)
pure function d2_t_dxy_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: lm, l2m2, n2, dl2m2, temp, temp1, dn2
   ! NOT READY
   Ecc = 0.0d0
!    lm = l*m
!    l2m2 = l*l + m*m
!    dl2m2 = 2.0d0*(dl*l + dm*m)
!    n2 = n*n
!    dn2 = 2.0d0*n*dn
!    temp = (n2 - 0.5d0*l2m2)
!    temp1 = (1.0d0+n2)
!    Ecc = ( m_sqrt3*( ((dn2 - 0.5d0*dl2m2)*Vdd_sigma + temp*dVdd_sigma) - 2.0d0*(dn2*Vdd_pi+n2*dVdd_pi) ) + m_sqrt3_half*(dn2*Vdd_delta + temp1*dVdd_delta) )*lm
!    Ecc = Ecc + ( m_sqrt3*(temp*Vdd_sigma - 2.0d0*n2*Vdd_pi) + m_sqrt3_half*temp1*Vdd_delta )*(dl*m + l*dm)
end function d2_t_dxy_d3z2_r2


! 21)
pure function d2_t_dyz_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: mn, l2m2, n2, temp, temp1, dn2, dl2m2
   ! NOT READY
   Ecc = 0.0d0
!    mn = m*n
!    l2m2 = l*l + m*m
!    dl2m2 = 2.0d0*(l*dl + m*dm)
!    n2 = n*n
!    dn2 = 2.0d0*n*dn
!    temp = (n2 - 0.5d0*l2m2)
!    temp1 = (l2m2 - n2)
!    Ecc = ( temp*Vdd_sigma + temp1*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*(dm*n+m*dn)
!    Ecc = Ecc + ( (dn2 - 0.5d0*dl2m2)*Vdd_sigma + temp*dVdd_sigma + (dl2m2 - dn2)*Vdd_pi + temp1*dVdd_pi - 0.5d0*(dl2m2*Vdd_delta + l2m2*dVdd_delta) )*mn
!    Ecc = Ecc*m_sqrt3
end function d2_t_dyz_d3z2_r2


! 22)
pure function d2_t_dxz_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: ln, l2m2, n2, temp, temp1, dn2, dl2m2
   ! NOT READY
   Ecc = 0.0d0
!    ln = l*n
!    l2m2 = l*l + m*m
!    dl2m2 = 2.0d0*(l*dl + m*dm)
!    n2 = n*n
!    dn2 = 2.0d0*n*dn
!    temp = (n2 - 0.5d0*l2m2)
!    temp1 = (l2m2 - n2)
!    Ecc = ( temp*Vdd_sigma + temp1*Vdd_pi - 0.5d0*l2m2*Vdd_delta )*(dl*n+l*dn)
!    Ecc = Ecc + ( (dn2 - 0.5d0*dl2m2)*Vdd_sigma + temp*dVdd_sigma + (dl2m2 - dn2)*Vdd_pi + temp1*dVdd_pi - 0.5d0*(dl2m2*Vdd_delta + l2m2*dVdd_delta) )*ln
!    Ecc = Ecc*m_sqrt3
end function d2_t_dxz_d3z2_r2


! 23)
pure function d2_t_dx2_y2_dx2_y2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2m2, l2_m2, l2_m2_2, n2, dl2m2, dl2_m2_2, l2, m2, ldl, mdm, dn2
   ! NOT READY
   Ecc = 0.0d0
!    l2 = l*l
!    m2 = m*m
!    l2m2 = l2 + m2
!    ldl = l*dl
!    mdm = m*dm
!    dl2m2 = 2.0d0*(ldl + mdm)
!    l2_m2 = l2 - m2
!    l2_m2_2 = l2_m2*l2_m2
!    dl2_m2_2 = 4.0d0*l2_m2*(ldl - mdm)
!    n2 = n*n
!    dn2 = 2.0d0*n*dn
!    Ecc = 0.75d0*(dl2_m2_2*Vdd_sigma + l2_m2_2*dVdd_sigma)
!    Ecc = Ecc + (dl2m2 - dl2_m2_2)*Vdd_pi + (l2m2 - l2_m2_2)*dVdd_pi
!    Ecc = Ecc + (dn2+0.25d0*dl2_m2_2)*Vdd_delta + (n2+0.25d0*l2_m2_2)*dVdd_delta
end function d2_t_dx2_y2_dx2_y2


! 24)
pure function d2_t_dx2_y2_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2m2, l2_m2, n2, dn2, l2, m2, dl2, dm2, dl2m2, dl2_m2, temp, temp1
   ! NOT READY
   Ecc = 0.0d0
!    l2 = l*l
!    dl2 = 2.0d0*l*dl
!    m2 = m*m
!    dm2 = 2.0d0*m*dm
!    l2m2 = l2 + m2
!    dl2m2 = dl2 + dm2
!    l2_m2 = l2 - m2
!    dl2_m2 = dl2 - dm2
!    n2 = n*n
!    dn2 = 2.0d0*n*dn
!    temp = (n2 - 0.5d0*l2m2)
!    temp1 = (1.0d0+n2)
!    Ecc = ( 0.5d0*temp*Vdd_sigma - n2*Vdd_pi + 0.25d0*temp1*Vdd_delta )*dl2_m2
!    Ecc = Ecc + ( 0.5d0*((dn2 - 0.5d0*dl2m2)*Vdd_sigma+temp*dVdd_sigma) - (dn2*Vdd_pi + n2*dVdd_pi) + 0.25d0*(dn2*Vdd_delta + temp1*dVdd_delta) )*l2_m2
!    Ecc = Ecc*m_sqrt3
end function d2_t_dx2_y2_d3z2_r2


! 25)
pure function d2_t_d3z2_r2_d3z2_r2(l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta) result (Ecc)
   real(8) :: Ecc
   real(8), intent(in) :: l, dl, m, dm, n, dn, Vdd_sigma, dVdd_sigma, Vdd_pi, dVdd_pi, Vdd_delta, dVdd_delta
   real(8) :: l2m2, dl2m2, n2, temp, dtemp, dn2
   ! NOT READY
   Ecc = 0.0d0
!    l2m2 = l*l + m*m
!    dl2m2 = 2.0d0*( dl*l + dm*m )
!    n2 = n*n
!    dn2 = 2.0d0*n*dn
!    temp = n2 - 0.5d0*l2m2
!    dtemp = dn2 - 0.5d0*dl2m2
!    Ecc = (2.0d0*dtemp*Vdd_sigma + temp*dVdd_sigma)*temp
!    Ecc = Ecc + 3.0d0*((dn2*l2m2 + n2*dl2m2)*Vdd_pi + n2*l2m2*dVdd_pi)
!    Ecc = Ecc + 0.75d0*(2.0d0*dl2m2*Vdd_delta + l2m2*dVdd_delta)*l2m2
end function d2_t_d3z2_r2_d3z2_r2



!------------------------------------------------------------------------
! Arrenge the SK-terms by parts (easier to read):
pure function KS_single_term(l1, l2, dx, dy, dz, V) result(t_kk)
   real(8) :: t_kk
   integer, intent(in) :: l1, l2 ! orbital indices
   real(8), intent(in) :: dx, dy, dz, V   ! direction cosines, and the radial function
   !---------------------
   t_kk = 0.0d0
end function KS_single_term



pure function KS_Cmnj_orbital(l1, l2, j, l, n, m) result(Cmnj) ! for up to sp3d5 (no f-orbitals currently)
   real(8) :: Cmnj   ! orbital coefficient for a single term; Table IV in Ref. [1]
   integer, intent(in) :: l1, l2, j ! orbital indices; j=sigma, pi, or delta
   real(8), intent(in) :: l, n, m   ! direction cosines
   ! Reminder:
   ! l=1: s
   ! l=2: p_x
   ! l=3: p_y
   ! l=4: p_z
   ! l=5: d_xy
   ! l=6: d_yz
   ! l=7: d_xz
   ! l=8: d_x^2-y^2
   ! l=9: d_3z^2-r^2
   !-----------
   ! j=0: sigma
   ! j=1: pi
   ! j=2: delta
   !===============
   real(8) :: alpha, beta

   Cmnj = 0.0d0   ! for all but those defined below
   if ((j < 0) .or. (j > 2)) return ! those do not exist, something must be wrong

   select case (l1)
   case (1) ! l=0: s

      if (j /= 0) return ! s-orbitals are non-zero only for j=0, nothing to do here
      ! If j=0, get the value:
      select case (l2) ! s s
      !---------------
      case (1) ! l=0: s s
         Cmnj = 1.0d0

      !---------------
      case (2) ! l=1: s p_x
         Cmnj = l

      !---------------
      case (3) ! l=1: s p_y
         Cmnj = m

      !---------------
      case (4) ! l=1: s p_z
         Cmnj = n

      !---------------
      case (5) ! l=2: s d_xy
         Cmnj = m_sqrt3 * l*m    ! d1,2 = dx, dy, or dz

      !---------------
      case (6) ! l=2: s d_yz
         Cmnj = m_sqrt3 * m*n    ! d1,2 = dx, dy, or dz

      !---------------
      case (7) ! l=2: s d_xz
         Cmnj = m_sqrt3 * l*n    ! d1,2 = dx, dy, or dz

      !---------------
      case (8) ! l=3: s d_x^2-y^2
         beta = l**2 - m**2
         Cmnj = m_sqrt3_half * beta

      !---------------
      case (9) ! l=4: s d_3z^2-r^2
         alpha = l**2 + m**2
         Cmnj = n**2 - 0.5d0*alpha
      endselect ! select case (l2)

   !===============
   case (2) ! l=1: p_x

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_x s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -l
         endselect

      !---------------
      case (2) ! l=1: p_x p_x
         select case (j)
         case (0) ! sigma
            Cmnj = l**2
         case (1) ! pi
            Cmnj = 1.0d0 - l**2
         endselect

      !---------------
      case (3) ! l=1: p_x p_y
         select case (j)
         case (0) ! sigma
            Cmnj = l*m
         case (1) ! pi
            Cmnj = -l*m
         endselect

      !---------------
      case (4) ! l=1: p_x p_z
         select case (j)
         case (0) ! sigma
            Cmnj = l*n
         case (1) ! pi
            Cmnj = -l*n
         endselect

      !---------------
      case (5) ! l=2: p_x d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * l**2 * m
         case (1) ! pi
            Cmnj = m*(1.0d0 - 2.0d0*l**2)
         endselect

      !---------------
      case (6) ! l=2: p_x d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * l * m * n
         case (1) ! pi
            Cmnj = -2.0d0*l*m*n
         endselect

      !---------------
      case (7) ! l=2: p_x d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * l**2 * n
         case (1) ! pi
            Cmnj = n*(1.0d0 - 2.0d0*l**2)
         endselect

      !---------------
      case (8) ! l=3: p_x d_x^2-y^2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * l*beta
         case (1) ! pi
            Cmnj = l*(1.0d0 - beta)
         endselect

      !---------------
      case (9) ! l=4: p_x d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = l*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*l*n**2
         endselect
      endselect ! select case (l2)

   !===============
   case (3) ! l=1: p_y

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_y s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -m
         endselect

      !---------------
      case (2) ! l=1: p_y p_x
         select case (j)
         case (0) ! sigma
            Cmnj = l*m
         case (1) ! pi
            Cmnj = -l*m
         endselect

      !---------------
      case (3) ! l=1: p_y p_y
         select case (j)
         case (0) ! sigma
            Cmnj = m**2
         case (1) ! pi
            Cmnj = 1.0d0 - m**2
         endselect

      !---------------
      case (4) ! l=1: p_y p_z
         select case (j)
         case (0) ! sigma
            Cmnj = m*n
         case (1) ! pi
            Cmnj = -m*n
         endselect

      !---------------
      case (5) ! l=2: p_y d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * m**2 * l
         case (1) ! pi
            Cmnj = l*(1.0d0 - 2.0d0*m**2)
         endselect

      !---------------
      case (6) ! l=2: p_y d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * m**2 * n
         case (1) ! pi
            Cmnj = n*(1.0d0 - 2.0d0*m**2)
         endselect

      !---------------
      case (7) ! l=2: p_y d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * l * m * n
         case (1) ! pi
            Cmnj = -2.0d0*l*m*n
         endselect

      !---------------
      case (8) ! l=3: p_y d_x^2-y^2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * m*beta
         case (1) ! pi
            Cmnj = -m*(1.0d0 + beta)
         endselect

      !---------------
      case (9) ! l=4: p_y d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = m*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*m*n**2
         endselect
      endselect ! select case (l2)

   !===============
   case (4) ! l=1: p_z

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_z s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -n
         endselect

      !---------------
      case (2) ! l=1: p_z p_x
         select case (j)
         case (0) ! sigma
            Cmnj = l*n
         case (1) ! pi
            Cmnj = -l*n
         endselect

      !---------------
      case (3) ! l=1: p_z p_y
         select case (j)
         case (0) ! sigma
            Cmnj = m*n
         case (1) ! pi
            Cmnj = -m*n
         endselect

      !---------------
      case (4) ! l=1: p_z p_z
          select case (j)
         case (0) ! sigma
            Cmnj = n**2
         case (1) ! pi
            Cmnj = 1.0d0 - n**2
         endselect

      !---------------
      case (5) ! l=2: p_z d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l * m * n
         case (1) ! pi
            Cmnj = 2.0d0*l*m*n
         endselect

      !---------------
      case (6) ! l=2: p_z d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * n**2 * m
         case (1) ! pi
            Cmnj = m*(1.0d0 - 2.0d0*n**2)
         endselect

      !---------------
      case (7) ! l=2: p_z d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * n**2 * l
         case (1) ! pi
            Cmnj = l*(1.0d0 - 2.0d0*n**2)
         endselect

      !---------------
      case (8) ! l=3: p_z d_x^2-y^2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * n*beta
         case (1) ! pi
            Cmnj = -n*beta
         endselect

      !---------------
      case (9) ! l=4: p_z d_3z^2-r^2
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = n*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*n*alpha
         endselect
      endselect ! select case (l2)

   !===============
   case (5) ! l=2: d_xy

      select case (l2)
      case (1) ! l=0: d_xy s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * l*m   ! d1,2 = dx, dy, or dz
         endselect

      !---------------
      case (2) ! l=1: d_xy p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l**2 * m
         case (1) ! pi
            Cmnj = -m*(1.0d0 - 2.0d0*l**2)
         endselect

      !---------------
      case (3) ! l=1: d_xy p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l * m**2
         case (1) ! pi
            Cmnj = -l*(1.0d0 - 2.0d0*m**2)
         endselect

      !---------------
      case (4) ! l=1: d_xy p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l * m * n
         case (1) ! pi
            Cmnj = 2.0d0*l*m*n
         endselect

      !---------------
      case (5) ! l=2: d_xy d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * l**2 * m**2
         case (1) ! pi
            alpha = l**2 + m**2
            Cmnj = alpha - 4.0d0*l**2 * m**2
         case (2) ! delta
            Cmnj = n**2 + l**2 * m**2
         endselect

      !---------------
      case (6) ! l=2: d_xy d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * l*n * m**2
         case (1) ! pi
            Cmnj = l*n*(1.0d0-4.0d0*m**2)
         case (2) ! delta
            Cmnj = l*n*(m**2 - 1.0d0)
         endselect

      !---------------
      case (7) ! l=2: d_xy d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * m*n * l**2
         case (1) ! pi
            Cmnj = m*n*(1.0d0-4.0d0*l**2)
         case (2) ! delta
            Cmnj = m*n*(l**2 - 1.0d0)
         endselect

      !---------------
      case (8) ! l=3: d_xy d_x^2-y^2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * l*m * beta
         case (1) ! pi
            Cmnj = -2.0d0 * l*m * beta
         case (2) ! delta
            Cmnj = 0.5d0 * l*m * beta
         endselect

      !---------------
      case (9) ! l=4: d_xy d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = m_sqrt3 * l*m * (n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = -2.0d0*m_sqrt3 * l*m * n**2
         case (2) ! delta
            Cmnj = m_sqrt3_half * l*m * (1.0d0 + n**2)
         endselect
      endselect ! select case (l2)


   !===============
   case (6) ! l=2: d_yz

      select case (l2)
      case (1) ! l=0: d_yz s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * m*n
         endselect

      !---------------
      case (2) ! l=1: d_yz p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l * m * n
         case (1) ! pi
            Cmnj = 2.0d0*l*m*n
         endselect

      !---------------
      case (3) ! l=1: d_yz p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * m**2 * n
         case (1) ! pi
            Cmnj = -n*(1.0d0 - 2.0d0*m**2)
         endselect

      !---------------
      case (4) ! l=1: d_yz p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * n**2 * m
         case (1) ! pi
            Cmnj = -m*(1.0d0 - 2.0d0*n**2)
         endselect

      !---------------
      case (5) ! l=2: d_yz d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * l*n * m**2
         case (1) ! pi
            Cmnj = l*n*(1.0d0-4.0d0*m**2)
         case (2) ! delta
            Cmnj = l*n*(m**2 - 1.0d0)
         endselect

      !---------------
      case (6) ! l=2: d_yz d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * m**2 * n**2
         case (1) ! pi
            alpha = n**2 + m**2
            Cmnj = alpha - 4.0d0*n**2 * m**2
         case (2) ! delta
            Cmnj = l**2 + n**2 * m**2
         endselect

      !---------------
      case (7) ! l=2: d_yz d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * m*l * n**2
         case (1) ! pi
            Cmnj = m*l*(1.0d0-4.0d0*n**2)
         case (2) ! delta
            Cmnj = m*l*(n**2 - 1.0d0)
         endselect

      !---------------
      case (8) ! l=3: d_yz d_x^2-y^2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * m*n * beta
         case (1) ! pi
            Cmnj = -m*n * (1.0d0 + 2.0d0*beta)
         case (2) ! delta
            Cmnj = m*n *(1.0d0 + 0.5d0*beta)
         endselect

      !---------------
      case (9) ! l=4: d_yz d_3z^2-r^2
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 *m*n * (n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3 *m*n * (alpha - n**2)
         case (2) ! delta
            Cmnj = -m_sqrt3_half * m*n * alpha
         endselect
      endselect ! select case (l2)


   !===============
   case (7) ! l=2: d_xz

      select case (l2)
      case (1) ! l=0: d_xz s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * l*n
         endselect

      !---------------
      case (2) ! l=1: d_xz p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l**2 * n
         case (1) ! pi
            Cmnj = -n*(1.0d0 - 2.0d0*l**2)
         endselect

      !---------------
      case (3) ! l=1: d_xz p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * l * m * n
         case (1) ! pi
            Cmnj = 2.0d0*l*m*n
         endselect

      !---------------
      case (4) ! l=1: d_xz p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * n**2 * l
         case (1) ! pi
            Cmnj = -l*(1.0d0 - 2.0d0*n**2)
         endselect

      !---------------
      case (5) ! l=2: d_xz d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * m*n * l**2
         case (1) ! pi
            Cmnj = m*n*(1.0d0-4.0d0*l**2)
         case (2) ! delta
            Cmnj = m*n*(l**2 - 1.0d0)
         endselect

      !---------------
      case (6) ! l=2: d_xz d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * m*l * n**2
         case (1) ! pi
            Cmnj = m*l*(1.0d0-4.0d0*n**2)
         case (2) ! delta
            Cmnj = m*l*(n**2 - 1.0d0)
         endselect

      !---------------
      case (7) ! l=2: d_xz d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * l**2 * n**2
         case (1) ! pi
            alpha = n**2 + l**2
            Cmnj = alpha - 4.0d0*n**2 * l**2
         case (2) ! delta
            Cmnj = m**2 + n**2 * l**2
         endselect

      !---------------
      case (8) ! l=3: d_xz d_x^2-y^2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * l*n * beta
         case (1) ! pi
            Cmnj = l*n * (1.0d0 - 2.0d0*beta)
         case (2) ! delta
            Cmnj = -l*n *(1.0d0 - 0.5d0*beta)
         endselect

      !---------------
      case (9) ! l=4: d_xz d_3z^2-r^2
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 *l*n * (n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3 *l*n * (alpha - n**2)
         case (2) ! delta
            Cmnj = -m_sqrt3_half * l*n * alpha
         endselect
      endselect ! select case (l2)

   !===============
   case (8) ! l=3: d_x^2-y^2

      select case (l2)
      case (1) ! l=0: d_x^2-y^2 s
         select case (j)
         case (0) ! sigma
            beta = l**2 - m**2
            Cmnj = m_sqrt3_half * beta
         endselect

      !---------------
      case (2) ! l=1: d_x^2-y^2 p_x
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * l*beta
         case (1) ! pi
            Cmnj = -l*(1.0d0 - beta)
         endselect

      !---------------
      case (3) ! l=1: d_x^2-y^2 p_y
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * m*beta
         case (1) ! pi
            Cmnj = m*(1.0d0 + beta)
         endselect

      !---------------
      case (4) ! l=1: d_x^2-y^2 p_z
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * n*beta
         case (1) ! pi
            Cmnj = n*beta
         endselect

      !---------------
      case (5) ! l=2: d_x^2-y^2 d_xy
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * l*m * beta
         case (1) ! pi
            Cmnj = -2.0d0 * l*m * beta
         case (2) ! delta
            Cmnj = 0.5d0 * l*m * beta
         endselect

      !---------------
      case (6) ! l=2: d_x^2-y^2 d_yz
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * m*n * beta
         case (1) ! pi
            Cmnj = -m*n * (1.0d0 + 2.0d0*beta)
         case (2) ! delta
            Cmnj = m*n *(1.0d0 + 0.5d0*beta)
         endselect

      !---------------
      case (7) ! l=2: d_x^2-y^2 d_xz
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * l*n * beta
         case (1) ! pi
            Cmnj = l*n * (1.0d0 - 2.0d0*beta)
         case (2) ! delta
            Cmnj = -l*n *(1.0d0 - 0.5d0*beta)
         endselect

      !---------------
      case (8) ! l=3: d_x^2-y^2 d_x^2-y^2
         alpha = l**2 + m**2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = 0.75d0*beta**2
         case (1) ! pi
            Cmnj = alpha - beta**2
         case (2) ! delta
            Cmnj = n**2 + 0.25d0*beta**2
         endselect

      !---------------
      case (9) ! l=4: d_x^2-y^2 d_3z^2-r^2
         alpha = l**2 + m**2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * beta*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*n**2*beta
         case (2) ! delta
            Cmnj = 0.25d0*m_sqrt3*(1.0d0 + n**2)*beta
         endselect
      endselect ! select case (l2)

   !===============
   case (9) ! l=4: d_3z^2-r^2

      select case (l2)
      case (1) ! l=0: d_3z^2-r^2 s
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = n**2 - 0.5d0*alpha
         endselect

      !---------------
      case (2) ! l=1: d_3z^2-r^2 p_x
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = -l*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*l*n**2
         endselect

      !---------------
      case (3) ! l=1: d_3z^2-r^2 p_y
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = -m*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*m*n**2
         endselect

      !---------------
      case (4) ! l=1: d_3z^2-r^2 p_z
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = n*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*n*alpha
         endselect

      !---------------
      case (5) ! l=2: d_3z^2-r^2 d_xy
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            Cmnj = m_sqrt3 * l*m * (n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = -2.0d0*m_sqrt3 * l*m * n**2
         case (2) ! delta
            Cmnj = m_sqrt3_half * l*m * (1.0d0 + n**2)
         endselect

      !---------------
      case (6) ! l=2: d_3z^2-r^2 d_yz
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 *m*n * (n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3 *m*n * (alpha - n**2)
         case (2) ! delta
            Cmnj = -m_sqrt3_half * m*n * alpha
         endselect

      !---------------
      case (7) ! l=2: d_3z^2-r^2 d_xz
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 *l*n * (n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = m_sqrt3 *l*n * (alpha - n**2)
         case (2) ! delta
            Cmnj = -m_sqrt3_half * l*n * alpha
         endselect

      !---------------
      case (8) ! l=3: d_3z^2-r^2 d_x^2-y^2
         alpha = l**2 + m**2
         beta = l**2 - m**2
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * beta*(n**2 - 0.5d0*alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*n**2*beta
         case (2) ! delta
            Cmnj = 0.25d0*m_sqrt3*(1.0d0 + n**2)*beta
         endselect

      !---------------
      case (9) ! l=4: d_3z^2-r^2 d_3z^2-r^2
         alpha = l**2 + m**2
         select case (j)
         case (0) ! sigma
            Cmnj = (n**2 - 0.5d0*alpha)**2
         case (1) ! pi
            Cmnj = 3.0d0*m**2 * alpha
         case (2) ! delta
            Cmnj = 0.75d0*alpha**2
         endselect
      endselect ! select case (l2)
   endselect ! select case (l1)

end function KS_Cmnj_orbital



! First derivative of KS-orbital angular part:
pure function d_KS_Cmnj_orbital(l1, l2, j, l, n, m, dl, dm, dn) result(Cmnj) ! for up to sp3d5 (no f-orbitals currently)
   real(8) :: Cmnj   ! orbital coefficient for a single term; Table IV in Ref. [1]
   integer, intent(in) :: l1, l2, j ! orbital indices; j=sigma, pi, or delta
   real(8), intent(in) :: l, n, m   ! direction cosines
   real(8), intent(in) :: dl, dm, dn   ! derivatives of the direction cosines
   ! Reminder:
   ! l=1: s
   ! l=2: p_x
   ! l=3: p_y
   ! l=4: p_z
   ! l=5: d_xy
   ! l=6: d_yz
   ! l=7: d_xz
   ! l=8: d_x^2-y^2
   ! l=9: d_3z^2-r^2
   !-----------
   ! j=0: sigma
   ! j=1: pi
   ! j=2: delta
   !===============
   real(8) :: alpha, beta, d_alpha, d_beta

   Cmnj = 0.0d0   ! for all but those defined below
   if ((j < 0) .or. (j > 2)) return ! those do not exist, something must be wrong

   select case (l1)
   case (1) ! l=0: s

      if (j /= 0) return ! s-orbitals are non-zero only for j=0, nothing to do here
      ! If j=0, get the value:
      select case (l2) ! s s
      !---------------
      case (1) ! l=0: s s
         Cmnj = 0.0d0

      !---------------
      case (2) ! l=1: s p_x
         Cmnj = dl

      !---------------
      case (3) ! l=1: s p_y
         Cmnj = dm

      !---------------
      case (4) ! l=1: s p_z
         Cmnj = dn

      !---------------
      case (5) ! l=2: s d_xy
         Cmnj = m_sqrt3 * (dl*m + l*dm)    ! d1,2 = dx, dy, or dz

      !---------------
      case (6) ! l=2: s d_yz
         Cmnj = m_sqrt3 * (dm*n + m*dn)    ! d1,2 = dx, dy, or dz

      !---------------
      case (7) ! l=2: s d_xz
         Cmnj = m_sqrt3 * (dl*n + l*dn)    ! d1,2 = dx, dy, or dz

      !---------------
      case (8) ! l=3: s d_x^2-y^2
         d_beta = 2.0d0*(l*dl - m*dm)
         Cmnj = m_sqrt3_half * d_beta

      !---------------
      case (9) ! l=4: s d_3z^2-r^2)
         Cmnj = 2.0d0*n*dn - (l*dl + m*dm)
      endselect ! select case (l2)

   !===============
   case (2) ! l=1: p_x

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_x s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -dl
         endselect

      !---------------
      case (2) ! l=1: p_x p_x
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*l*dl
         case (1) ! pi
            Cmnj = -2.0d0*l*dl
         endselect

      !---------------
      case (3) ! l=1: p_x p_y
         select case (j)
         case (0) ! sigma
            Cmnj = (dl*m + l*dm)
         case (1) ! pi
            Cmnj = -(dl*m + l*dm)
         endselect

      !---------------
      case (4) ! l=1: p_x p_z
         select case (j)
         case (0) ! sigma
            Cmnj = (dl*n + l*dn)
         case (1) ! pi
            Cmnj = -(dl*n + l*dn)
         endselect

      !---------------
      case (5) ! l=2: p_x d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*l*dl * m + l**2 * dm)
         case (1) ! pi
            Cmnj = dm*(1.0d0 - 2.0d0*l**2) + m*(-4.0d0*l*dl)
         endselect

      !---------------
      case (6) ! l=2: p_x d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (dl*m*n + l*dm*n + l*m*dn)
         case (1) ! pi
            Cmnj = -2.0d0 * (dl*m*n + l*dm*n + l*m*dn)
         endselect

      !---------------
      case (7) ! l=2: p_x d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*l*dl*n + l**2*dn)
         case (1) ! pi
            Cmnj = dn*(1.0d0 - 2.0d0*l**2) + n*(-4.0d0*l*dl)
         endselect

      !---------------
      case (8) ! l=3: p_x d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * (dl*beta + l*d_beta)
         case (1) ! pi
            Cmnj = dl*(1.0d0 - beta) + l*(-d_beta)
         endselect

      !---------------
      case (9) ! l=4: p_x d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = dl*(n**2 - 0.5d0*alpha) + l*(2.0d0*n*dn - 0.5d0*d_alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*(dl*n**2 + 2.0d0*l*n*dn)
         endselect
      endselect ! select case (l2)

   !===============
   case (3) ! l=1: p_y

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here


      select case (l2) ! p_y s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -dm
         endselect

      !---------------
      case (2) ! l=1: p_y p_x
         select case (j)
         case (0) ! sigma
            Cmnj = (dl*m + l*dm)
         case (1) ! pi
            Cmnj = -(dl*m + l*dm)
         endselect

      !---------------
      case (3) ! l=1: p_y p_y
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*m*dm
         case (1) ! pi
            Cmnj = -2.0d0*m*dm
         endselect

      !---------------
      case (4) ! l=1: p_y p_z
         select case (j)
         case (0) ! sigma
            Cmnj = (dm*n + m*dn)
         case (1) ! pi
            Cmnj = -(dm*n + m*dn)
         endselect

      !---------------
      case (5) ! l=2: p_y d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*m*dm * l + m**2 * dl)
         case (1) ! pi
            Cmnj = dl*(1.0d0 - 2.0d0*m**2) + l*(-4.0d0*m*dm)
         endselect

      !---------------
      case (6) ! l=2: p_y d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*m*dm * n + m**2 * dn)
         case (1) ! pi
            Cmnj = dn*(1.0d0 - 2.0d0*m**2) + n*(-4.0d0*m*dm)
         endselect

      !---------------
      case (7) ! l=2: p_y d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (dl*m*n + l*dm*n + l*m*dn)
         case (1) ! pi
            Cmnj = -2.0d0*(dl*m*n + l*dm*n + l*m*dn)
         endselect

      !---------------
      case (8) ! l=3: p_y d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * (dm*beta + m*d_beta)
         case (1) ! pi
            Cmnj = -(dm*(1.0d0 + beta) + m*d_beta )
         endselect

      !---------------
      case (9) ! l=4: p_y d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = dm*(n**2 - 0.5d0*alpha) + m*(2.0d0*n*dn - 0.5d0*d_alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*(dm*n**2 + 2.0d0*m*n*dn)
         endselect
      endselect ! select case (l2)

   !===============
   case (4) ! l=1: p_z

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_z s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -dn
         endselect

      !---------------
      case (2) ! l=1: p_z p_x
         select case (j)
         case (0) ! sigma
            Cmnj = dl*n + l*dn
         case (1) ! pi
            Cmnj = -(dl*n + l*dn)
         endselect

      !---------------
      case (3) ! l=1: p_z p_y
         select case (j)
         case (0) ! sigma
            Cmnj = dm*n + m*dn
         case (1) ! pi
            Cmnj = -(dm*n + m*dn)
         endselect

      !---------------
      case (4) ! l=1: p_z p_z
          select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*n*dn
         case (1) ! pi
            Cmnj = -2.0d0*n*dn
         endselect

      !---------------
      case (5) ! l=2: p_z d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (dl*m*n + l*dm*n + l*m*dn)
         case (1) ! pi
            Cmnj = 2.0d0*(dl*m*n + l*dm*n + l*m*dn)
         endselect

      !---------------
      case (6) ! l=2: p_z d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*n*dn * m + n**2 * dm)
         case (1) ! pi
            Cmnj = dm*(1.0d0 - 2.0d0*n**2) + m*(-4.0d0*n*dn)
         endselect

      !---------------
      case (7) ! l=2: p_z d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*n*dn * l + n**2 * dl)
         case (1) ! pi
            Cmnj = dl*(1.0d0 - 2.0d0*n**2) + l*(-4.0d0*n*dn)
         endselect

      !---------------
      case (8) ! l=3: p_z d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * (dn*beta + n*d_beta)
         case (1) ! pi
            Cmnj = -(dn*beta + n*d_beta)
         endselect

      !---------------
      case (9) ! l=4: p_z d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = dn*(n**2 - 0.5d0*alpha) + n*(2.0d0*n*dn - 0.5d0*d_alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*(dn*alpha + n*d_alpha)
         endselect
      endselect ! select case (l2)

   !===============
   case (5) ! l=2: d_xy

      select case (l2)
      case (1) ! l=0: d_xy s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 *(dl*m + l*dm)
         endselect

      !---------------
      case (2) ! l=1: d_xy p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*l*dl * m + l**2 * dm)
         case (1) ! pi
            Cmnj = -( dm*(1.0d0 - 2.0d0*l**2) + m*(-4.0d0*l*dl))
         endselect

      !---------------
      case (3) ! l=1: d_xy p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (dl*m**2 + 2.0d0*l*m*dm)
         case (1) ! pi
            Cmnj = -( dl*(1.0d0 - 2.0d0*m**2) + l*(-4.0d0*m*dm) )
         endselect

      !---------------
      case (4) ! l=1: d_xy p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (dl*m*n + l*dm*n + l*m*dn)
         case (1) ! pi
            Cmnj = 2.0d0*(dl*m*n + l*dm*n + l*m*dn)
         endselect

      !---------------
      case (5) ! l=2: d_xy d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 6.0d0 * (l*dl*m**2 + l**2*m*dm)
         case (1) ! pi
            !alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = d_alpha - 8.0d0*(l*dl * m*dm)
         case (2) ! delta
            Cmnj = 2.0d0*(n*dn + l*dl*m**2 + l**2*m*dm)
         endselect

      !---------------
      case (6) ! l=2: d_xy d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * (dl*n*m**2 + l*dn*m**2 + 2.0d0*l*n*m*dm)
         case (1) ! pi
            Cmnj = (dl*n+ l*dn)*(1.0d0-4.0d0*m**2) + l*n*(-8.0d0*m*dm)
         case (2) ! delta
            Cmnj = (dl*n + l*dn)*(m**2 - 1.0d0) + 2.0d0*l*n*m*dm
         endselect

      !---------------
      case (7) ! l=2: d_xy d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (dm*n + m*dn)*l**2 + 2.0d0*m*n*l*dl )
         case (1) ! pi
            Cmnj = (dm*n + m*dn)*(1.0d0-4.0d0*l**2) + m*n*(-8.0d0*l*dl)
         case (2) ! delta
            Cmnj = (dm*n + m*dn)*(l**2 - 1.0d0) + 2.0d0*m*n*l*dl
         endselect

      !---------------
      case (8) ! l=3: d_xy d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ((dl*m + l*dm)*beta + l*m*d_beta)
         case (1) ! pi
            Cmnj = -2.0d0 * ((dl*m + l*dm)*beta + l*m*d_beta)
         case (2) ! delta
            Cmnj = 0.5d0 * ((dl*m + l*dm)*beta + l*m*d_beta)
         endselect

      !---------------
      case (9) ! l=4: d_xy d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = m_sqrt3 * ( (dl*m + l*dm)*(n**2 - 0.5d0*alpha) + l*m*(2.0d0*n*dn - 0.5d0*d_alpha))
         case (1) ! pi
            Cmnj = -2.0d0*m_sqrt3 * ( (dl*m + l*dm)*n**2 + 2.0d0*l*m*n*dn)
         case (2) ! delta
            Cmnj = m_sqrt3_half * ( (dl*m + l*dm)*(1.0d0 + n**2) + l*m*(2.0d0*n*dn) )
         endselect
      endselect ! select case (l2)


   !===============
   case (6) ! l=2: d_yz

      select case (l2)
      case (1) ! l=0: d_yz s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (dm*n + m*dn)
         endselect

      !---------------
      case (2) ! l=1: d_yz p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * ( dl*m*n + l*dm*n + l*m*dn )
         case (1) ! pi
            Cmnj = 2.0d0*( dl*m*n + l*dm*n + l*m*dn )
         endselect

      !---------------
      case (3) ! l=1: d_yz p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * ( 2.0d0*m*dm*n + m**2*dn )
         case (1) ! pi
            Cmnj = -( dn*(1.0d0 - 2.0d0*m**2) + n*(-4.0d0*m*dm) )
         endselect

      !---------------
      case (4) ! l=1: d_yz p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*n*dn*m + n**2*dm)
         case (1) ! pi
            Cmnj = -( dm*(1.0d0 - 2.0d0*n**2) + m*(-4.0d0*n*dn) )
         endselect

      !---------------
      case (5) ! l=2: d_yz d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (dl*n + l*dn)*m**2 + 2.0d0*l*n*m*dm)
         case (1) ! pi
            Cmnj = (dl*n + l*dn)*(1.0d0-4.0d0*m**2) + l*n*(-8.0d0*m*dm)
         case (2) ! delta
            Cmnj = (dl*n + l*dn)*(m**2 - 1.0d0) + l*n*(2.0d0*m*dm)
         endselect

      !---------------
      case (6) ! l=2: d_yz d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 6.0d0 * (m*dm * n*dn)
         case (1) ! pi
            !alpha = n**2 + m**2
            d_alpha = 2.0d0*(n*dn + m*dm)
            Cmnj = d_alpha - 8.0d0*(n*dn*m**2 + n**2*m*dm)
         case (2) ! delta
            Cmnj = 2.0d0*(l*dl + n*dn*m**2 + n**2*m*dm)
         endselect

      !---------------
      case (7) ! l=2: d_yz d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (dm*l + m*dl)*n**2 + 2.0d0*m*l*n*dn )
         case (1) ! pi
            Cmnj = (dm*l + m*dl)*(1.0d0-4.0d0*n**2) + m*l*(-8.0d0*n*dn)
         case (2) ! delta
            Cmnj = (dm*l + m*dl)*(n**2 - 1.0d0) + m*l*(2.0d0*n*dn)
         endselect

      !---------------
      case (8) ! l=3: d_yz d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ( (dm*n + m*dn)*beta + m*n*d_beta )
         case (1) ! pi
            Cmnj = -( (dm*n + m*dn)*(1.0d0 + 2.0d0*beta) + m*n*2.0d0*d_beta )
         case (2) ! delta
            Cmnj = (dm*n + m*dn)*(1.0d0 + 0.5d0*beta) + m*n*0.5d0*d_beta
         endselect

      !---------------
      case (9) ! l=4: d_yz d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (dm*n + m*dn)*(n**2 - 0.5d0*alpha) + m*n*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (dm*n + m*dn)*(alpha - n**2) + m*n*(d_alpha - 2.0d0*n*dn))
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ((dm*n + m*dn)*alpha + m*n*d_alpha)
         endselect
      endselect ! select case (l2)


   !===============
   case (7) ! l=2: d_xz

      select case (l2)
      case (1) ! l=0: d_xz s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (dl*n + l*dn)
         endselect

      !---------------
      case (2) ! l=1: d_xz p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*l*dl*n + l**2*dn)
         case (1) ! pi
            Cmnj = -( dn*(1.0d0 - 2.0d0*l**2) + n*(-4.0d0*l*dl) )
         endselect

      !---------------
      case (3) ! l=1: d_xz p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (dl*m*n + l*dm*n + l*m*dn)
         case (1) ! pi
            Cmnj = 2.0d0*(dl*m*n + l*dm*n + l*m*dn)
         endselect

      !---------------
      case (4) ! l=1: d_xz p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*n*dn*l + n**2*dl)
         case (1) ! pi
            Cmnj = -( dl*(1.0d0 - 2.0d0*n**2) + l*(-4.0d0*n*dn) )
         endselect

      !---------------
      case (5) ! l=2: d_xz d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ((dm*n + m*dn)*l**2 + 2.0d0*m*n*l*dl)
         case (1) ! pi
            Cmnj = (dm*n + m*dn)*(1.0d0-4.0d0*l**2) + m*n*(-8.0d0*l*dl)
         case (2) ! delta
            Cmnj = (dm*n + m*dn)*(l**2 - 1.0d0) + 2.0d0*m*n*l*dl
         endselect

      !---------------
      case (6) ! l=2: d_xz d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (dm*l + m*dl)*n**2 + 2.0d0*m*l*n*dn )
         case (1) ! pi
            Cmnj = (dm*l + m*dl)*(1.0d0-4.0d0*n**2) + m*l*(-8.0d0*n*dn)
         case (2) ! delta
            Cmnj = (dm*l + m*dl)*(n**2 - 1.0d0) + 2.0d0*m*l*n*dn
         endselect

      !---------------
      case (7) ! l=2: d_xz d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 6.0d0*(l*dl*n**2 + l**2*n*dn)
         case (1) ! pi
            !alpha = n**2 + l**2
            d_alpha = 2.0d0*(n*dn + l*dl)
            Cmnj = d_alpha - 8.0d0*(n*dn*l**2 + n**2*l*dl)
         case (2) ! delta
            Cmnj = 2.0d0*(m*dm + n*dn*l**2 + n**2*l*dl)
         endselect

      !---------------
      case (8) ! l=3: d_xz d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ((dl*n + l*dn)*beta + l*n*d_beta)
         case (1) ! pi
            Cmnj = (dl*n + l*dn)*(1.0d0 - 2.0d0*beta) + l*n*(-2.0d0*d_beta)
         case (2) ! delta
            Cmnj = -( (dl*n + l*dn)*(1.0d0 - 0.5d0*beta) + l*n*(-0.5d0*d_beta) )
         endselect

      !---------------
      case (9) ! l=4: d_xz d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ((dl*n + l*dn)*(n**2 - 0.5d0*alpha) + l*n*(2.0d0*n*dn - 0.5d0*d_alpha))
         case (1) ! pi
            Cmnj = m_sqrt3 * ((dl*n + l*dn)*(alpha - n**2) + l*n*(d_alpha - 2.0d0*n*dn))
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ( (dl*n + l*dn)*alpha + l*n*d_alpha )
         endselect
      endselect ! select case (l2)

   !===============
   case (8) ! l=3: d_x^2-y^2

     select case (l2)
      case (1) ! l=0: d_x^2-y^2 s
         select case (j)
         case (0) ! sigma
            !beta = d1**2 - d2**2
            d_beta = 2.0d0*(l*dl - m*dm)
            Cmnj = m_sqrt3_half * d_beta
         endselect

      !---------------
      case (2) ! l=1: d_x^2-y^2 p_x
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * ( dl*beta + l*d_beta )
         case (1) ! pi
            Cmnj = -( dl*(1.0d0 - beta) + l*(-d_beta) )
         endselect

      !---------------
      case (3) ! l=1: d_x^2-y^2 p_y
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * ( dm*beta + m*d_beta )
         case (1) ! pi
            Cmnj = dm*(1.0d0 + beta) + m*d_beta
         endselect

      !---------------
      case (4) ! l=1: d_x^2-y^2 p_z
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * (dn*beta + n*d_beta)
         case (1) ! pi
            Cmnj = dn*beta + n*d_beta
         endselect

      !---------------
      case (5) ! l=2: d_x^2-y^2 d_xy
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ((dl*m + l*dm)*beta + l*m*d_beta)
         case (1) ! pi
            Cmnj = -2.0d0 * ((dl*m + l*dm)*beta + l*m*d_beta)
         case (2) ! delta
            Cmnj = 0.5d0 * ((dl*m + l*dm)*beta + l*m*d_beta)
         endselect

      !---------------
      case (6) ! l=2: d_x^2-y^2 d_yz
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ((dm*n + m*dn)*beta + m*n*d_beta)
         case (1) ! pi
            Cmnj = -( (dm*n + m*dn)*(1.0d0 + 2.0d0*beta) + m*n*2.0d0*d_beta )
         case (2) ! delta
            Cmnj = (dm*n + m*dn)*(1.0d0 + 0.5d0*beta) + m*n*0.5d0*d_beta
         endselect

      !---------------
      case (7) ! l=2: d_x^2-y^2 d_xz
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ((dl*n + l*dn)*beta + l*n*d_beta)
         case (1) ! pi
            Cmnj = (dl*n + l*dn)*(1.0d0 - 2.0d0*beta) + l*n*(-2.0d0*d_beta)
         case (2) ! delta
            Cmnj = -(dl*n + l*dn)*(1.0d0 - 0.5d0*beta) + l*n*0.5d0*d_beta
         endselect

      !---------------
      case (8) ! l=3: d_x^2-y^2 d_x^2-y^2
         !alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0*beta*d_beta
         case (1) ! pi
            Cmnj = d_alpha - 2.0d0*beta*d_beta
         case (2) ! delta
            Cmnj = 2.0d0*n*dn + 0.5d0*beta*d_beta
         endselect

      !---------------
      case (9) ! l=4: d_x^2-y^2 d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * ( d_beta*(n**2 - 0.5d0*alpha) + beta*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = -m_sqrt3*( 2.0d0*n*dn*beta + n**2*d_beta )
         case (2) ! delta
            Cmnj = 0.25d0*m_sqrt3*( 2.0d0*n*dn*beta + (1.0d0 + n**2)*d_beta )
         endselect
      endselect ! select case (l2)

   !===============
   case (9) ! l=4: d_3z^2-r^2

      select case (l2)
      case (1) ! l=0: d_3z^2-r^2 s
         select case (j)
         case (0) ! sigma
            !alpha = d1**2 + d2**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = 2.0d0*n*dn - 0.5d0*d_alpha
         endselect

      !---------------
      case (2) ! l=1: d_3z^2-r^2 p_x
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = -( dl*(n**2 - 0.5d0*alpha) + l*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3*( dl*n**2 + 2.0d0*l*n*dn )
         endselect

      !---------------
      case (3) ! l=1: d_3z^2-r^2 p_y
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = -( dm*(n**2 - 0.5d0*alpha) + m*(2.0d0*n*dn - 0.5d0*d_alpha))
         case (1) ! pi
            Cmnj = m_sqrt3*( dm*n**2 + 2.0d0*m*n*dn )
         endselect

      !---------------
      case (4) ! l=1: d_3z^2-r^2 p_z
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = dn*(n**2 - 0.5d0*alpha) + n*(2.0d0*n*dn - 0.5d0*d_alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*(dn*alpha + n*d_alpha)
         endselect

      !---------------
      case (5) ! l=2: d_3z^2-r^2 d_xy
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            alpha = 2.0d0*(l*dl + m*dm)
            Cmnj = m_sqrt3 * ( (dl*m + l*dm)*(n**2 - 0.5d0*alpha) + l*m*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = -2.0d0*m_sqrt3 * ( (dl*m + l*dm)*n**2 + 2.0d0*l*m*n*dn )
         case (2) ! delta
            Cmnj = m_sqrt3_half * ( (dl*m + l*dm)*(1.0d0 + n**2) + l*m*(2.0d0*n*dn) )
         endselect

      !---------------
      case (6) ! l=2: d_3z^2-r^2 d_yz
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (dm*n + m*dn)*(n**2 - 0.5d0*alpha) + m*n*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (dm*n + m*dn)*(alpha - n**2) + m*n*(d_alpha - 2.0d0*n*dn) )
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ((dm*n + m*dn)*alpha + m*n*d_alpha)
         endselect

      !---------------
      case (7) ! l=2: d_3z^2-r^2 d_xz
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (dl*n + l*dn)*(n**2 - 0.5d0*alpha) + l*n*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (dl*n + l*dn)*(alpha - n**2) + l*n*(d_alpha - 2.0d0*n*dn) )
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ( (dl*n + l*dn)*alpha + l*n*d_alpha )
         endselect

      !---------------
      case (8) ! l=3: d_3z^2-r^2 d_x^2-y^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * ( d_beta*(n**2 - 0.5d0*alpha) + beta*(2.0d0*n*dn - 0.5d0*d_alpha) )
         case (1) ! pi
            Cmnj = -m_sqrt3*( 2.0d0*n*dn*beta + n**2*d_beta )
         case (2) ! delta
            Cmnj = 0.25d0*m_sqrt3*( 2.0d0*n*dn*beta + (1.0d0 + n**2)*d_beta )
         endselect

      !---------------
      case (9) ! l=4: d_3z^2-r^2 d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*(n**2 - 0.5d0*alpha)*(2.0d0*n*dn - 0.5d0*d_alpha)
         case (1) ! pi
            Cmnj = 3.0d0*( 2.0d0*m*dm*alpha + m**2*d_alpha )
         case (2) ! delta
            Cmnj = 1.5d0*alpha*d_alpha
         endselect
      endselect ! select case (l2)
   endselect ! select case (l1)

end function d_KS_Cmnj_orbital




! Second derivative of KS-orbital angular part:
pure function d2_KS_Cmnj_orbital(l1, l2, j, l, n, m, dl, dm, dn, d2l, d2m, d2n) result(Cmnj) ! for up to sp3d5 (no f-orbitals currently)
   real(8) :: Cmnj   ! orbital coefficient for a single term; Table IV in Ref. [1]
   integer, intent(in) :: l1, l2, j ! orbital indices; j=sigma, pi, or delta
   real(8), intent(in) :: l, n, m   ! direction cosines
   real(8), intent(in) :: dl, dm, dn   ! derivatives of the direction cosines
   real(8), intent(in) :: d2l, d2m, d2n   ! second derivatives of the direction cosines
   ! Reminder:
   ! l=1: s
   ! l=2: p_x
   ! l=3: p_y
   ! l=4: p_z
   ! l=5: d_xy
   ! l=6: d_yz
   ! l=7: d_xz
   ! l=8: d_x^2-y^2
   ! l=9: d_3z^2-r^2
   !-----------
   ! j=0: sigma
   ! j=1: pi
   ! j=2: delta
   !===============
   real(8) :: alpha, beta, d_alpha, d_beta, d2_alpha, d2_beta

   Cmnj = 0.0d0   ! for all but those defined below
   if ((j < 0) .or. (j > 2)) return ! those do not exist, something must be wrong

   select case (l1)
   case (1) ! l=0: s

      if (j /= 0) return ! s-orbitals are non-zero only for j=0, nothing to do here
      ! If j=0, get the value:
      select case (l2) ! s s
      !---------------
      case (1) ! l=0: s s
         Cmnj = 0.0d0

      !---------------
      case (2) ! l=1: s p_x
         Cmnj = d2l

      !---------------
      case (3) ! l=1: s p_y
         Cmnj = d2m

      !---------------
      case (4) ! l=1: s p_z
         Cmnj = d2n

      !---------------
      case (5) ! l=2: s d_xy
         Cmnj = m_sqrt3 * (d2l*m + 2.0d0*dl*dm + l*d2m)

      !---------------
      case (6) ! l=2: s d_yz
         Cmnj = m_sqrt3 * (d2m*n + 2.0d0*dm*dn + m*d2n)

      !---------------
      case (7) ! l=2: s d_xz
         Cmnj = m_sqrt3 * (d2l*n + 2.0d0*dl*dn + l*d2n)

      !---------------
      case (8) ! l=3: s d_x^2-y^2
         !d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(l*d2l + dl**2 - dm**2 - m*d2m)
         Cmnj = m_sqrt3_half * d2_beta

      !---------------
      case (9) ! l=4: s d_3z^2-r^2)
         Cmnj = 2.0d0*(dn**2 + n*d2n) - (dl**2 + l*d2l + dm**2 + m*d2m)

      endselect ! select case (l2)

   !===============
   case (2) ! l=1: p_x

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_x s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -d2l
         endselect

      !---------------
      case (2) ! l=1: p_x p_x
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*(dl**2 + l*d2l)
         case (1) ! pi
            Cmnj = -2.0d0*(dl**2 + l*d2l)
         endselect

      !---------------
      case (3) ! l=1: p_x p_y
         select case (j)
         case (0) ! sigma
            Cmnj = (d2l*m + 2.0d0*dl*dm + l*d2m)
         case (1) ! pi
            Cmnj = -(d2l*m + 2.0d0*dl*dm + l*d2m)
         endselect

      !---------------
      case (4) ! l=1: p_x p_z
         select case (j)
         case (0) ! sigma
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)
         case (1) ! pi
            Cmnj = -(d2l*n + 2.0d0*dl*dn + l*d2n)
         endselect

      !---------------
      case (5) ! l=2: p_x d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( 2.0d0*(dl*dl*m + l*d2l*m + l*dl*dm) + (l**2*d2m + 2.0d0*l*dl*dm) )
         case (1) ! pi
            Cmnj = d2m*(1.0d0 - 2.0d0*l**2) + dm*(-4.0d0*l*dl) + dm*(-4.0d0*l*dl) + m*(-4.0d0*(dl**2+l*d2l))
         endselect

      !---------------
      case (6) ! l=2: p_x d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*dl*dm*n + 2.0d0*dl*m*dn + 2.0d0*l*dm*dn)
         case (1) ! pi
            Cmnj = -2.0d0 * (d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*dl*dm*n + 2.0d0*dl*m*dn + 2.0d0*l*dm*dn)
         endselect

      !---------------
      case (7) ! l=2: p_x d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*(dl**2*n + l*d2l*n + l*dl*dn) + &
                              2.0d0*l*dl*dn + l**2*d2n)
         case (1) ! pi
            Cmnj = d2n*(1.0d0 - 2.0d0*l**2) + 2.0d0*dn*(-4.0d0*l*dl) + &
                   n*(-4.0d0*(dl**2 + l*d2l))
         endselect

      !---------------
      case (8) ! l=3: p_x d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * (d2l*beta + 2.0d0*dl*d_beta + l*d2_beta)
         case (1) ! pi
            Cmnj = d2l*(1.0d0 - beta) + 2.0d0*l*(-d_beta) + l*(-d2_beta)
         endselect

      !---------------
      case (9) ! l=4: p_x d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
            Cmnj = d2l*(n**2 - 0.5d0*alpha) + 2.0d0*dl*(2.0d0*n*dn - 0.5d0*d_alpha) + &
                   l*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*(d2l*n**2 + 4.0d0*l*n*dn + 2.0d0*l*n*d2n)
         endselect
      endselect ! select case (l2)

   !===============
   case (3) ! l=1: p_y

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here


      select case (l2) ! p_y s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -d2m
         endselect

      !---------------
      case (2) ! l=1: p_y p_x
         select case (j)
         case (0) ! sigma
            Cmnj = (d2l*m + 2.0d0*dl*dm + l*d2m)
         case (1) ! pi
            Cmnj = -(d2l*m + 2.0d0*dl*dm + l*d2m)
         endselect

      !---------------
      case (3) ! l=1: p_y p_y
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*(dm**2 + m*d2m)
         case (1) ! pi
            Cmnj = -2.0d0*(dm**2 + m*d2m)
         endselect

      !---------------
      case (4) ! l=1: p_y p_z
         select case (j)
         case (0) ! sigma
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)
         case (1) ! pi
            Cmnj = -(d2m*n + 2.0d0*dm*dn + m*d2n)
         endselect

      !---------------
      case (5) ! l=2: p_y d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*(dm**2+m*d2m)*l + 4.0d0*m*dm*dl + m**2*d2l)
         case (1) ! pi
            Cmnj = d2l*(1.0d0 - 2.0d0*m**2) + 2.0d0*dl*(-4.0d0*m*dm) + l*(-4.0d0*(dm**2 + m*d2m))
         endselect

      !---------------
      case (6) ! l=2: p_y d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*(dm**2 + m*d2m)*n + 4.0d0*m*dm*dn + m**2*d2n)
         case (1) ! pi
            Cmnj = d2n*(1.0d0 - 2.0d0*m**2) + 2.0d0*dn*(-4.0d0*m*dm) + n*(-4.0d0*(dm**2 + m*d2m))
         endselect

      !---------------
      case (7) ! l=2: p_y d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn) )
         case (1) ! pi
            Cmnj = -2.0d0*( d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn) )
         endselect

      !---------------
      case (8) ! l=3: p_y d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(l*d2l + dl**2 - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * (d2m*beta + 2.0d0*dm*d_beta + m*d2_beta)
         case (1) ! pi
            Cmnj = -(d2m*(1.0d0 + beta) + 2.0d0*dm*d_beta + m*d2_beta)
         endselect

      !---------------
      case (9) ! l=4: p_y d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(l*d2l + dl**2 + m*d2m + dm**2)
            Cmnj = d2m*(n**2 - 0.5d0*alpha) + 2.0d0*dm*(2.0d0*n*dn - 0.5d0*d_alpha) + m*(2.0d0*(dn**2+n*d2n) - 0.5d0*d2_alpha)
         case (1) ! pi
            Cmnj = -m_sqrt3*( d2m*n**2 + 2.0d0*dm*n*dn + 2.0d0*m*(dn**2 + n*d2n) )
         endselect
      endselect ! select case (l2)

   !===============
   case (4) ! l=1: p_z

      if (j > 1) return ! p-orbitals are non-zero only for j=0 or 1; nothing to do here

      select case (l2) ! p_z s
      case (1) ! l=0: s
         select case (j)
         case (0) ! sigma
            Cmnj = -d2n
         endselect

      !---------------
      case (2) ! l=1: p_z p_x
         select case (j)
         case (0) ! sigma
            Cmnj = d2l*n + 2.0d0*dl*dn + l*d2n
         case (1) ! pi
            Cmnj = -(d2l*n + 2.0d0*dl*dn + l*d2n)
         endselect

      !---------------
      case (3) ! l=1: p_z p_y
         select case (j)
         case (0) ! sigma
            Cmnj = d2m*n + 2.0d0*dm*dn + m*d2n
         case (1) ! pi
            Cmnj = -(d2m*n + 2.0d0*dm*dn + m*d2n)
         endselect

      !---------------
      case (4) ! l=1: p_z p_z
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*(dn**2 + n*d2n)
         case (1) ! pi
            Cmnj = -2.0d0*(dn**2 + n*d2n)
         endselect

      !---------------
      case (5) ! l=2: p_z d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         case (1) ! pi
            Cmnj = 2.0d0*(d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         endselect

      !---------------
      case (6) ! l=2: p_z d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*((dn**2 + n*d2n)*m + n*dn*dm) + 2.0d0*n*dn*dm + n**2*d2m)
         case (1) ! pi
            Cmnj = d2m*(1.0d0 - 2.0d0*n**2) + 2.0d0*dm*(-4.0d0*n*dn) + m*(-4.0d0*(dn**2 + n*d2n))
         endselect

      !---------------
      case (7) ! l=2: p_z d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (2.0d0*((dn**2+n*d2n)*l + n*dn*dl) + 2.0d0*n*dn*dl + n**2*d2l)
         case (1) ! pi
            Cmnj = d2l*(1.0d0 - 2.0d0*n**2) + 2.0d0*dl*(-4.0d0*n*dn) + l*(-4.0d0*(dn**2 + n*d2n))
         endselect

      !---------------
      case (8) ! l=3: p_z d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * (d2n*beta + 2.0d0*dn*d_beta + n*d2_beta)
         case (1) ! pi
            Cmnj = -(d2n*beta + 2.0d0*dn*d_beta + n*d2_beta)
         endselect

      !---------------
      case (9) ! l=4: p_z d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = d2n*(n**2 - 0.5d0*alpha) + 2.0d0*dn*(2.0d0*n*dn - 0.5d0*d_alpha) + n*(2.0d0*(dn**2+n*d2n) - 0.5d0*d2_alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*(d2n*alpha + 2.0d0*dn*d_alpha + n*d2_alpha)
         endselect
      endselect ! select case (l2)

   !===============
   case (5) ! l=2: d_xy

      select case (l2)
      case (1) ! l=0: d_xy s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 *(d2l*m + 2.0d0*dl*dm + l*d2m)
         endselect

      !---------------
      case (2) ! l=1: d_xy p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*((dl**2 + l*d2l)*m + l*dl*dm) + 2.0d0*l*dl*dm + l**2*d2m)
         case (1) ! pi
            Cmnj = -( d2m*(1.0d0 - 2.0d0*l**2) + 2.0d0*dm*(-4.0d0*l*dl) + m*(-4.0d0*(dl**2 + l*d2l)))
         endselect

      !---------------
      case (3) ! l=1: d_xy p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (d2l*m**2 + 4.0d0*dl*m*dm + 2.0d0*l*(dm**2 + m*d2m))
         case (1) ! pi
            Cmnj = -( d2l*(1.0d0 - 2.0d0*m**2) + 2.0d0*dl*(-4.0d0*m*dm) + l*(-4.0d0*(dm**2 + m*d2m)) )
         endselect

      !---------------
      case (4) ! l=1: d_xy p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         case (1) ! pi
            Cmnj = 2.0d0*(d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         endselect

      !---------------
      case (5) ! l=2: d_xy d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 6.0d0 * ( (dl**2 + l*d2l)*m**2 + 4.0d0*l*dl*m*dm + l**2*(dm**2 + m*d2m) )
         case (1) ! pi
            !alpha = l**2 + m**2
            !d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
            Cmnj = d2_alpha - 8.0d0*((dl**2 + l*d2l)*m*dm + l*dl*(dm**2 + m*d2m))
         case (2) ! delta
            Cmnj = 2.0d0*( (dn**2 + n*d2n) + (dl**2 + l*d2l)*m**2 + 4.0d0*l*dl*m*dm + l**2*(dm**2 + m*d2m) )
         endselect

      !---------------
      case (6) ! l=2: d_xy d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( d2l*n*m**2 + 2.0d0*dl*dn*m**2 + 4.0d0*dl*n*m*dm + &
                             l*d2n*m**2 + 4.0d0*l*dn*m*dm + &
                             2.0d0*(l*n*dm**2 + l*n*m*d2m) )
         case (1) ! pi
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)*(1.0d0-4.0d0*m**2) + 16.0d0*(dl*n + l*dn)*(-m*dm) + l*n*(-8.0d0*(dm**2 + m*d2m))
         case (2) ! delta
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)*(m**2 - 1.0d0) + (dl*n + l*dn)*(4.0d0*m*dm) + 2.0d0*l*n*(dm**2 + m*d2m)
         endselect

      !---------------
      case (7) ! l=2: d_xy d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*l**2 + (dm*n + m*dn)*2.0d0*l*dl + &
                              2.0d0*(dm*n*l*dl + m*dn*l*dl + m*n*dl*dl + m*n*l*d2l) )
         case (1) ! pi
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)*(1.0d0-4.0d0*l**2) + (dm*n + m*dn)*(-16.0d0*l*dl) + m*n*(-8.0d0*(dl**2 + l*d2l))
         case (2) ! delta
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)*(l**2 - 1.0d0) + (dm*n + m*dn)*(4.0d0*l*dl) + 2.0d0*m*n*(dl**2 + l*d2l)
         endselect

      !---------------
      case (8) ! l=3: d_xy d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ((d2l*m + 2.0d0*dl*dm + l*d2m)*beta + 2.0d0*(dl*m + l*dm)*d_beta + l*m*d_beta + l*m*d2_beta)
         case (1) ! pi
            Cmnj = -2.0d0 * ((d2l*m + 2.0d0*dl*dm + l*d2m)*beta + 2.0d0*(dl*m + l*dm)*d_beta + l*m*d_beta + l*m*d2_beta)
         case (2) ! delta
            Cmnj = 0.5d0 * ((d2l*m + 2.0d0*dl*dm + l*d2m)*beta + 2.0d0*(dl*m + l*dm)*d_beta + l*m*d_beta + l*m*d2_beta)
         endselect

      !---------------
      case (9) ! l=4: d_xy d_3z^2-r^2
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(l*d2l + dl**2 + m*d2m + dm**2)
            Cmnj = m_sqrt3 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*(n**2 - 0.5d0*alpha) + &
                               (dl*m + l*dm)*(2.0d0*n*dn - 0.5d0*d_alpha) + l*m*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha))
         case (1) ! pi
            Cmnj = -2.0d0*m_sqrt3 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*n**2 + (dl*m + l*dm)*2.0d0*n*dn + 2.0d0*l*m*(dn**2 + n*d2n))
         case (2) ! delta
            Cmnj = m_sqrt3_half * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*(1.0d0 + n**2) + (dl*m + l*dm)*(4.0d0*n*dn) + l*m*(dn**2 + n*d2n) )
         endselect
      endselect ! select case (l2)


   !===============
   case (6) ! l=2: d_yz

      select case (l2)
      case (1) ! l=0: d_yz s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (d2m*n + 2.0d0*dm*dn + m*d2n)
         endselect

      !---------------
      case (2) ! l=1: d_yz p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * ( d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         case (1) ! pi
            Cmnj = 2.0d0*( d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         endselect

      !---------------
      case (3) ! l=1: d_yz p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * ( 2.0d0*((dm**2 + m*d2m)*n + m*dm*dn) + 2.0d0*m*dm*dn + m**2*d2n )
         case (1) ! pi
            Cmnj = -( d2n*(1.0d0 - 2.0d0*m**2) + dn*(-4.0d0*m*dm) + n*(-4.0d0*(dm**2 + m*d2m)) )
         endselect

      !---------------
      case (4) ! l=1: d_yz p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*((dn**2 + n*d2n)*m + n*dn*dm) + 2.0d0*n*dn*dm + n**2*d2m)
         case (1) ! pi
            Cmnj = -( d2m*(1.0d0 - 2.0d0*n**2) + 2.0d0*dm*(-4.0d0*n*dn) + m*(-4.0d0*(dn**2 + n*d2n)) )
         endselect

      !---------------
      case (5) ! l=2: d_yz d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*m**2 + 2.0d0*(dl*n + l*dn)*m*dm + 2.0d0*l*n*(dm**2 + m*d2m) )
         case (1) ! pi
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)*(1.0d0-4.0d0*m**2) + (dl*n + l*dn)*(-16.0d0*m*dm) + l*n*(-8.0d0*(dm**2 + m*d2m))
         case (2) ! delta
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)*(m**2 - 1.0d0) + 2.0d0*(dl*n + l*dn)*(2.0d0*m*dm) + l*n*(2.0d0*(dm**2 + m*d2m))
         endselect

      !---------------
      case (6) ! l=2: d_yz d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 6.0d0 * (dm**2 + m*d2m * dn**2 + n*d2n)
         case (1) ! pi
            !alpha = n**2 + m**2
            !d_alpha = 2.0d0*(n*dn + m*dm)
            d2_alpha = 2.0d0*(dn**2 + n*d2n + dm**2 + m*d2m)
            Cmnj = d2_alpha - 8.0d0*((dn**2 + n*d2n)*m**2 + 4.0d0*n*dn*m*dm + n**2*(dm**2 + m*d2m))
         case (2) ! delta
            Cmnj = 2.0d0*((dl**2+l*d2l) + (dn**2+n*d2n)*m**2 + 2.0d0*n*dn*m*dm + n**2*(dm**2+m*d2m))
         endselect

      !---------------
      case (7) ! l=2: d_yz d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (d2m*l + 2.0d0*dm*dl + m*d2l)*n**2 + 4.0d0*(dm*l + m*dl)*n*dn + 2.0d0*m*l*(dn**2 + n*d2n) )
         case (1) ! pi
            Cmnj = (d2m*l + 2.0d0*dm*dl + m*d2l)*(1.0d0-4.0d0*n**2) + (dm*l + m*dl)*(-16.0d0*n*dn) + m*l*(-8.0d0*(dn**2 + n*d2n))
         case (2) ! delta
            Cmnj = (d2m*l + 2.0d0*dm*dl + m*d2l)*(n**2 - 1.0d0) + 4.0d0*(dm*l + m*dl)*(2.0d0*n*dn) + m*l*(2.0d0*(dn**2 + n*d2n))
         endselect

      !---------------
      case (8) ! l=3: d_yz d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(l*d2l + dl**2 - m*d2m - dm**2)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*beta + 2.0d0*(dm*n + m*dn)*d_beta + m*n*d2_beta )
         case (1) ! pi
            Cmnj = -( (d2m*n + 2.0d0*dm*dn + m*d2n)*(1.0d0 + 2.0d0*beta) + (dm*n + m*dn)*2.0d0*d_beta + m*n*2.0d0*d2_beta )
         case (2) ! delta
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)*(1.0d0 + 0.5d0*beta) + (dm*n + m*dn)*d_beta + m*n*0.5d0*d2_beta
         endselect

      !---------------
      case (9) ! l=4: d_yz d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*(n**2 - 0.5d0*alpha) + &
                               2.0d0*(dm*n + m*dn)*(2.0d0*n*dn - 0.5d0*d_alpha) + m*n*(2.0d0*(dn**2+n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*(alpha - n**2) + 2.0d0*(dm*n + m*dn)*(d_alpha - 2.0d0*n*dn) + &
                                m*n*(d2_alpha - 2.0d0*(dn**2 + n*d2n)) )
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ((d2m*n + 2.0d0*dm*dn + m*d2n)*alpha + 2.0d0*(dm*n + m*dn)*d_alpha + m*n*d2_alpha)
         endselect
      endselect ! select case (l2)


   !===============
   case (7) ! l=2: d_xz

      select case (l2)
      case (1) ! l=0: d_xz s
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * (d2l*n + 2.0d0*dl*dn + l*d2n)
         endselect

      !---------------
      case (2) ! l=1: d_xz p_x
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*((dl**2 + l*d2l)*n + l*dl*dn) + 2.0d0*l*dl*dn + l**2*d2n)
         case (1) ! pi
            Cmnj = -( d2n*(1.0d0 - 2.0d0*l**2) + dn*(-8.0d0*l*dl) + n*(-4.0d0*(dl**2 + l*d2l)) )
         endselect

      !---------------
      case (3) ! l=1: d_xz p_y
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         case (1) ! pi
            Cmnj = 2.0d0*(d2l*m*n + l*d2m*n + l*m*d2n + 2.0d0*(dl*dm*n + l*dm*dn + dl*m*dn))
         endselect

      !---------------
      case (4) ! l=1: d_xz p_z
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3 * (2.0d0*((dn**2 + n*d2n)*l + n*dn*dl) + 2.0d0*n*dn*dl + n**2*d2l)
         case (1) ! pi
            Cmnj = -( d2l*(1.0d0 - 2.0d0*n**2) + dl*(-8.0d0*n*dn) + l*(-4.0d0*(dn**2 + n*d2n)) )
         endselect

      !---------------
      case (5) ! l=2: d_xz d_xy
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*l**2 + 4.0d0*(dm*n + m*dn)*l*dl + 2.0d0*m*n*(dl**2 + l*d2l) )
         case (1) ! pi
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)*(1.0d0-4.0d0*l**2) + (dm*n + m*dn)*(-16.0d0*l*dl) + m*n*(-8.0d0*(dl**2 + l*d2l))
         case (2) ! delta
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)*(l**2 - 1.0d0) + 4.0d0*(dm*n + m*dn)*l*dl + 2.0d0*m*n*(dl**2 + l*d2l)
         endselect

      !---------------
      case (6) ! l=2: d_xz d_yz
         select case (j)
         case (0) ! sigma
            Cmnj = 3.0d0 * ( (d2m*l + 2.0d0*dm*dl + m*d2l)*n**2 + 4.0d0*(dm*l + m*dl)*n*dn + 2.0d0*m*l*(dn**2 + n*d2n) )
         case (1) ! pi
            Cmnj = (d2m*l + 2.0d0*dm*dl + m*d2l)*(1.0d0-4.0d0*n**2) + (dm*l + m*dl)*(-16.0d0*n*dn) + m*l*(-8.0d0*(dn**2 + n*d2n))
         case (2) ! delta
            Cmnj = (d2m*l + 2.0d0*dm*dl + m*d2l)*(n**2 - 1.0d0) + 4.0d0*(dm*l + m*dl)*n*dn + 2.0d0*m*l*(dn**2 + n*d2n)
         endselect

      !---------------
      case (7) ! l=2: d_xz d_xz
         select case (j)
         case (0) ! sigma
            Cmnj = 6.0d0*( (dl**2 + l*d2l)*n**2 + 4.0d0*l*dl*n*dn + l**2*(dn**2 + n*d2n) )
         case (1) ! pi
            !alpha = n**2 + l**2
            !d_alpha = 2.0d0*(n*dn + l*dl)
            d2_alpha = 2.0d0*(dn**2 + n*d2n + dl**2 + l*d2l)
            Cmnj = d2_alpha - 8.0d0*((dn**2 + n*d2n)*l**2 + 4.0d0*n*dn*l*dl + n**2*(dl**2 + l*d2l))
         case (2) ! delta
            Cmnj = 2.0d0*( (dm**2 + m*d2m) + (dn**2 + n*d2n)*l**2 + 4.0d0*n*dn*l*dl + n**2*(dl**2 + l*d2l) )
         endselect

      !---------------
      case (8) ! l=3: d_xz d_x^2-y^2
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*beta + 2.0d0*(dl*n + l*dn)*d_beta + l*n*d2_beta )
         case (1) ! pi
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)*(1.0d0 - 2.0d0*beta) + (dl*n + l*dn)*(-4.0d0*d_beta) + l*n*(-2.0d0*d2_beta)
         case (2) ! delta
            Cmnj = -( (d2l*n + 2.0d0*dl*dn + l*d2n)*(1.0d0 - 0.5d0*beta) + (dl*n + l*dn)*(-d_beta) + l*n*(-0.5d0*d2_beta) )
         endselect

      !---------------
      case (9) ! l=4: d_xz d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*(n**2 - 0.5d0*alpha) + 2.0d0*(dl*n + l*dn)*(2.0d0*n*dn - 0.5d0*d_alpha) + &
                               l*n*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*(alpha - n**2) + 2.0d0*(dl*n + l*dn)*(d_alpha - 2.0d0*n*dn) + &
                                l*n*(d2_alpha - 2.0d0*(dn**2 + n*d2n)) )
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*alpha + 2.0d0*(dl*n + l*dn)*d_alpha + l*n*d2_alpha )
         endselect
      endselect ! select case (l2)

   !===============
   case (8) ! l=3: d_x^2-y^2

     select case (l2)
      case (1) ! l=0: d_x^2-y^2 s
         select case (j)
         case (0) ! sigma
            !beta = d1**2 - d2**2
            !d_beta = 2.0d0*(l*dl - m*dm)
            d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
            Cmnj = m_sqrt3_half * d2_beta
         endselect

      !---------------
      case (2) ! l=1: d_x^2-y^2 p_x
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * ( d2l*beta + 2.0d0*dl*d_beta + l*d2_beta )
         case (1) ! pi
            Cmnj = -( d2l*(1.0d0 - beta) + 2.0d0*dl*(-d_beta) + l*(-d2_beta) )
         endselect

      !---------------
      case (3) ! l=1: d_x^2-y^2 p_y
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * ( d2m*beta + 2.0d0*dm*d_beta + m*d2_beta )
         case (1) ! pi
            Cmnj = d2m*(1.0d0 + beta) + 2.0d0*dm*d_beta + m*d2_beta
         endselect

      !---------------
      case (4) ! l=1: d_x^2-y^2 p_z
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = -m_sqrt3_half * (d2n*beta + 2.0d0*dn*d_beta + n*d2_beta)
         case (1) ! pi
            Cmnj = d2n*beta + 2.0d0*dn*d_beta + n*d2_beta
         endselect

      !---------------
      case (5) ! l=2: d_x^2-y^2 d_xy
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*beta + 2.0d0*(dl*m + l*dm)*d_beta + l*m*d2_beta )
         case (1) ! pi
            Cmnj = -2.0d0 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*beta + 2.0d0*(dl*m + l*dm)*d_beta + l*m*d2_beta )
         case (2) ! delta
            Cmnj = 0.5d0 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*beta + 2.0d0*(dl*m + l*dm)*d_beta + l*m*d2_beta )
         endselect

      !---------------
      case (6) ! l=2: d_x^2-y^2 d_yz
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*beta + (dm*n + m*dn)*d_beta + m*n*d2_beta )
         case (1) ! pi
            Cmnj = -( (d2m*n + 2.0d0*dm*dn + m*d2n)*(1.0d0 + 2.0d0*beta) + (dm*n + m*dn)*4.0d0*d_beta + m*n*2.0d0*d2_beta )
         case (2) ! delta
            Cmnj = (d2m*n + 2.0d0*dm*dn + m*d2n)*(1.0d0 + 0.5d0*beta) + (dm*n + m*dn)*d_beta + m*n*0.5d0*d2_beta
         endselect

      !---------------
      case (7) ! l=2: d_x^2-y^2 d_xz
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*beta + 2.0d0*(dl*n + l*dn)*d_beta + l*n*d2_beta )
         case (1) ! pi
            Cmnj = (d2l*n + 2.0d0*dl*dn + l*d2n)*(1.0d0 - 2.0d0*beta) + (dl*n + l*dn)*(-4.0d0*d_beta) + l*n*(-2.0d0*d2_beta)
         case (2) ! delta
            Cmnj = -(d2l*n + 2.0d0*dl*dn + l*d2n)*(1.0d0 - 0.5d0*beta) - (dl*n + l*dn)*d_beta + l*n*0.5d0*d2_beta
         endselect

      !---------------
      case (8) ! l=3: d_x^2-y^2 d_x^2-y^2
         !d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 1.5d0*(d_beta**2 + beta*d2_beta)
         case (1) ! pi
            Cmnj = d2_alpha - 2.0d0*(d_beta**2 + beta*d2_beta)
         case (2) ! delta
            Cmnj = 2.0d0*(dn**2 + n*d2n) + 0.5d0*(d_beta**2 + beta*d2_beta)
         endselect

      !---------------
      case (9) ! l=4: d_x^2-y^2 d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * ( d2_beta*(n**2 - 0.5d0*alpha) + 2.0d0*d_beta*(2.0d0*n*dn - 0.5d0*d_alpha) + &
                                    beta*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = -m_sqrt3*( 2.0d0*((dn**2 + n*d2n)*beta + n*dn*d_beta) + 2.0d0*n*dn*d_beta + n**2*d2_beta )
         case (2) ! delta
            Cmnj = 0.25d0*m_sqrt3*( 2.0d0*((dn**2 + n*d2n)*beta + n*dn*d_beta) + 2.0d0*n*dn*d_beta + (1.0d0 + n**2)*d2_beta )
         endselect
      endselect ! select case (l2)
   !===============
   case (9) ! l=4: d_3z^2-r^2

      select case (l2)
      case (1) ! l=0: d_3z^2-r^2 s
         select case (j)
         case (0) ! sigma
            !alpha = d1**2 + d2**2
            !d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
            Cmnj = 2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha
         endselect

      !---------------
      case (2) ! l=1: d_3z^2-r^2 p_x
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
            Cmnj = -( d2l*(n**2 - 0.5d0*alpha) + 2.0d0*dl*(2.0d0*n*dn - 0.5d0*d_alpha) + l*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3*( d2l*n**2 + 4.0d0*dl*n*dn + 2.0d0*l*(dn**2 + n*d2n) )
         endselect

      !---------------
      case (3) ! l=1: d_3z^2-r^2 p_y
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            d_alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
            Cmnj = -( d2m*(n**2 - 0.5d0*alpha) + 2.0d0*dm*(2.0d0*n*dn - 0.5d0*d_alpha) + m*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3*( d2m*n**2 + 4.0d0*dm*n*dn + 2.0d0*m*(dn**2 + n*d2n) )
         endselect

      !---------------
      case (4) ! l=1: d_3z^2-r^2 p_z
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = d2n*(n**2 - 0.5d0*alpha) + 2.0d0*dn*(2.0d0*n*dn - 0.5d0*d_alpha) + n*(2.0d0*(dn**2+n*d2n) - 0.5d0*d2_alpha)
         case (1) ! pi
            Cmnj = m_sqrt3*( d2n*alpha + 2.0d0*dn*d_alpha + n*d2_alpha )
         endselect

      !---------------
      case (5) ! l=2: d_3z^2-r^2 d_xy
         select case (j)
         case (0) ! sigma
            alpha = l**2 + m**2
            alpha = 2.0d0*(l*dl + m*dm)
            d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
            Cmnj = m_sqrt3 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*(n**2 - 0.5d0*alpha) + 2.0d0*(dl*m + l*dm)*(2.0d0*n*dn - 0.5d0*d_alpha) + &
                                l*m*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = -2.0d0*m_sqrt3 * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*n**2 + 2.0d0*(dl*m + l*dm)*n*dn + 2.0d0*l*m*(dn**2 + n*d2n) )
         case (2) ! delta
            Cmnj = m_sqrt3_half * ( (d2l*m + 2.0d0*dl*dm + l*d2m)*(1.0d0 + n**2) + 2.0d0*(dl*m + l*dm)*(2.0d0*n*dn) + &
                                     l*m*(2.0d0*(dn**2 + n*d2n)) )
         endselect

      !---------------
      case (6) ! l=2: d_3z^2-r^2 d_yz
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*(n**2 - 0.5d0*alpha) + 2.0d0*(dm*n + m*dn)*(2.0d0*n*dn - 0.5d0*d_alpha) + &
                                m*n*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (d2m*n + 2.0d0*dm*dn + m*d2n)*(alpha - n**2) + 2.0d0*(dm*n + m*dn)*(d_alpha - 2.0d0*n*dn) + &
                                m*n*(d2_alpha - 2.0d0*(dn**2 + n*d2n)) )
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ((d2m*n + 2.0d0*dm*dn + m*d2n)*alpha + 2.0d0*(dm*n + m*dn)*d_alpha + m*n*d2_alpha)
         endselect

      !---------------
      case (7) ! l=2: d_3z^2-r^2 d_xz
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*(n**2 - 0.5d0*alpha) + &
                               2.0d0*(dl*n + l*dn)*(2.0d0*n*dn - 0.5d0*d_alpha) + l*n*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = m_sqrt3 * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*(alpha - n**2) + 2.0d0*(dl*n + l*dn)*(d_alpha - 2.0d0*n*dn) + &
                                l*n*(d2_alpha - 2.0d0*(dn**2 + n*d2n)) )
         case (2) ! delta
            Cmnj = -m_sqrt3_half * ( (d2l*n + 2.0d0*dl*dn + l*d2n)*alpha + 2.0d0*(dl*n + l*dn)*d_alpha + l*n*d2_alpha )
         endselect

      !---------------
      case (8) ! l=3: d_3z^2-r^2 d_x^2-y^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         beta = l**2 - m**2
         d_beta = 2.0d0*(l*dl - m*dm)
         d2_beta = 2.0d0*(dl**2 + l*d2l - dm**2 - m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = m_sqrt3_half * ( d2_beta*(n**2 - 0.5d0*alpha) + 2.0d0*d_beta*(2.0d0*n*dn - 0.5d0*d_alpha) + &
                                    beta*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = -m_sqrt3*( 2.0d0*((dn**2 + n*d2n)*beta + n*dn*d_beta) + 2.0d0*n*dn*d_beta + n**2*d2_beta )
         case (2) ! delta
            Cmnj = 0.25d0*m_sqrt3*( 2.0d0*((dn**2 + n*d2n)*beta + n*dn*d_beta)+ 2.0d0*n*dn*d_beta + (1.0d0 + n**2)*d2_beta )
         endselect

      !---------------
      case (9) ! l=4: d_3z^2-r^2 d_3z^2-r^2
         alpha = l**2 + m**2
         d_alpha = 2.0d0*(l*dl + m*dm)
         d2_alpha = 2.0d0*(dl**2 + l*d2l + dm**2 + m*d2m)
         select case (j)
         case (0) ! sigma
            Cmnj = 2.0d0*( (2.0d0*n*dn - 0.5d0*d_alpha)**2 + (n**2 - 0.5d0*alpha)*(2.0d0*(dn**2 + n*d2n) - 0.5d0*d2_alpha) )
         case (1) ! pi
            Cmnj = 3.0d0*( 2.0d0*((dm**2 + m*d2m)*alpha + m*dm*d_alpha) + 2.0d0*m*dm*d_alpha + m**2*d2_alpha )
         case (2) ! delta
            Cmnj = 1.5d0*(d_alpha**2 + alpha*d2_alpha)
         endselect
      endselect ! select case (l2)
   endselect ! select case (l1)

end function d2_KS_Cmnj_orbital





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
pure subroutine KS_sp3d5_hetero_TEST(Vsr12, Vsr21, l, m, n, ts)
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
   ts(1,2) = t_s_px(l, Vsr12(2))	! s px
   ts(1,3) = t_s_px(m, Vsr12(2))	! s py
   ts(1,4) = t_s_px(n, Vsr12(2))	! s pz
   ts(1,5) = t_s_dab(l, m, Vsr12(3))	! s dxy
   ts(1,6) = t_s_dab(l, n, Vsr12(3))	! s dxz
   ts(1,7) = t_s_dab(m, n, Vsr12(3))	! s dyz
   ts(1,8) = t_s_dx2_y2(l, m, Vsr12(3))	! s dx2-y2
   ts(1,9) = t_s_dz2_r2(l, m, n, Vsr12(3))	! s d3z2-r2


   ts(2,1) = -t_s_px(l, Vsr21(2))		! px s
   ts(2,2) = t_pa_pa(l, Vsr12(4), Vsr12(5))		! px px
   ts(2,3) = t_pa_pb(l, m, Vsr12(4), Vsr12(5))	! px py
   ts(2,4) = t_pa_pb(l, n, Vsr12(4), Vsr12(5))	! px pz
   ts(2,5) = t_px_dxy(l, m, Vsr12(6), Vsr12(7))	! px dxy
   ts(2,6) = t_px_dxy(l, n, Vsr12(6), Vsr12(7))	! px dxz
   ts(2,7) = t_px_dyz(l, m, n, Vsr12(6), Vsr12(7))	! px dyz
   ts(2,8) = t_px_dx2_y2(l, m, Vsr12(6), Vsr12(7))	! px dx2-y2
   ts(2,9) = t_pa_d3z2_r2(l, l, m, n, Vsr12(6), Vsr12(7))	! px d3z2-r2


   ts(3,1) = -t_s_px(m, Vsr21(2))	! py s
   ts(3,2) = t_pa_pb(l, m, Vsr21(4), Vsr21(5))	! py px
   ts(3,3) = t_pa_pa(m, Vsr12(4), Vsr12(5))		! py py
   ts(3,4) = t_pa_pb(m, n, Vsr12(4), Vsr12(5))	! py pz
   ts(3,5) = t_px_dxy(m, l, Vsr12(6), Vsr12(7))	! py dxy
   ts(3,6) = t_px_dyz(m, l, n, Vsr12(6), Vsr12(7))	! py dxz
   ts(3,7) = t_px_dxy(m, n,  Vsr12(6), Vsr12(7))	! py dyz
   ts(3,8) = t_py_dx2_y2(l, m, Vsr12(6), Vsr12(7))	! py dx2-y2
   ts(3,9) = t_pa_d3z2_r2(m, l, m, n, Vsr12(6), Vsr12(7))	! py d3z2-r2


   ts(4,1) = -t_s_px(n, Vsr21(2))	! pz s
   ts(4,2) = t_pa_pb(l, n, Vsr21(4), Vsr21(5))	! pz px
   ts(4,3) = t_pa_pb(m, n, Vsr21(4), Vsr21(5))	! pz py
   ts(4,4) = t_pa_pa(n, Vsr12(4), Vsr12(5))		! pz pz
   ts(4,5) = t_px_dyz(n, l, m, Vsr12(6), Vsr12(7))	! pz dxy
   ts(4,6) = t_px_dxy(n, l, Vsr12(6), Vsr12(7))	! pz dxz
   ts(4,7) = t_px_dxy(n, m, Vsr12(6), Vsr12(7))	! pz dyz
   ts(4,8) = t_pz_dx2_y2(l, m, n, Vsr12(6), Vsr12(7))	! pz dx2-y2
   ts(4,9) = t_pz_d3z2_r2(l, m, n, Vsr12(6), Vsr12(7))	! pz d3z2-r2


   ts(5,1) = t_s_dab(l, m, Vsr21(3))	! dxy s
   ts(5,2) = -t_px_dxy(l, m, Vsr21(6), Vsr21(7))	! dxy px
   ts(5,3) = -t_px_dxy(m, l, Vsr21(6), Vsr21(7))	! dxy py
   ts(5,4) = -t_px_dyz(n, l, m, Vsr21(6), Vsr21(7))	! dxy pz
   ts(5,5) = t_dab_dab(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dxy
   ts(5,6) = t_dab_dbg(m, l, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dxz
   ts(5,7) = t_dab_dbg(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dyz
   ts(5,8) = t_dxy_dx2_y2(l, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxy dx2-y2
   ts(5,9) = t_dxy_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dxy d3z2-r2


   ts(6,1) = t_s_dab(l, n, Vsr21(3))	! dxz s
   ts(6,2) = -t_px_dxy(l, n, Vsr21(6), Vsr21(7))	! dxz px
   ts(6,3) = -t_px_dyz(m, l, n, Vsr21(6), Vsr21(7))	! dxz py
   ts(6,4) = -t_px_dxy(n, l, Vsr21(6), Vsr21(7))	! dxz pz
   ts(6,5) = t_dab_dbg(m, l, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dxz dxy
   ts(6,6) = t_dab_dab(l, n, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dxz
   ts(6,7) = t_dab_dbg(l, n, m, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dyz
   ts(6,8) = t_dxz_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dxz dx2-y2
   ts(6,9) = t_dxz_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dxz d3z2-r2


   ts(7,1) = t_s_dab(m, n, Vsr21(3))	! dyz s
   ts(7,2) = -t_px_dyz(l, m, n, Vsr21(6), Vsr21(7))	! dyz px
   ts(7,3) = -t_px_dxy(m, n,  Vsr21(6), Vsr21(7))	! dyz py
   ts(7,4) = -t_px_dxy(n, m, Vsr21(6), Vsr21(7))	! dyz pz
   ts(7,5) = t_dab_dbg(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dyz dxy
   ts(7,6) = t_dab_dbg(l, n, m, Vsr21(8), Vsr21(9), Vsr21(10))	! dyz dxz
   ts(7,7) = t_dab_dab(m, n, l, Vsr12(8), Vsr12(9), Vsr12(10))	! dyz dyz
   ts(7,8) = t_dyz_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dyz dx2-y2
   ts(7,9) = t_dyz_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dyz d3z2-r2


   ts(8,1) = t_s_dx2_y2(l, m, Vsr21(3))	! dx2-y2 s
   ts(8,2) = -t_px_dx2_y2(l, m, Vsr21(6), Vsr21(7))	! dx2-y2 px
   ts(8,3) = -t_py_dx2_y2(l, m, Vsr21(6), Vsr21(7))	! dx2-y2 py
   ts(8,4) = -t_pz_dx2_y2(l, m, n, Vsr21(6), Vsr21(7))	! dx2-y2 pz
   ts(8,5) = t_dxy_dx2_y2(l, m, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dxy
   ts(8,6) = t_dxz_dx2_y2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dxz
   ts(8,7) = t_dyz_dx2_y2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10))	! dx2-y2 dyz
   ts(8,8) = t_dx2_y2_dx2_y2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10))	! dx2-y2 dx2-y2
   ts(8,9) = t_dx2_y2_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! dx2-y2 d3z2-r2


   ts(9,1) = t_s_dz2_r2(l, m, n, Vsr21(3))	! d3z2-r2 s
   ts(9,2) = -t_pa_d3z2_r2(l, l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 px
   ts(9,3) = -t_pa_d3z2_r2(m, l, m, n, Vsr21(6), Vsr21(7))	! d3z2-r2 py
   ts(9,4) = -t_pz_d3z2_r2(l, m, n, Vsr21(6), Vsr21(7))	    ! d3z2-r2 pz
   ts(9,5) = t_dxy_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dxy
   ts(9,6) = t_dxz_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dxz
   ts(9,7) = t_dyz_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dyz
   ts(9,8) = t_dx2_y2_d3z2_r2(l, m, n, Vsr21(8), Vsr21(9), Vsr21(10)) 	! d3z2-r2 dx2-y2
   ts(9,9) = t_d3z2_r2_d3z2_r2(l, m, n, Vsr12(8), Vsr12(9), Vsr12(10)) 	! d3z2-r2 d3z2-r2
end subroutine KS_sp3d5_hetero_TEST


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
   
   ts(2,1) = -t_s_px(l, Vsr21(2))   ! (px s) = -(s px)
   ts(2,2) = t_pa_pa(l, Vsr12(3), Vsr12(4))
   ts(2,3) = t_pa_pb(l, m, Vsr12(3), Vsr12(4))
   ts(2,4) = t_pa_pb(l, n, Vsr12(3), Vsr12(4))
   
   ts(3,1) = -t_s_px(m, Vsr21(2))   ! (py s) = -(s py)
   ts(3,2) = t_pa_pb(l, m, Vsr21(3), Vsr21(4))
   ts(3,3) = t_pa_pa(m, Vsr12(3), Vsr12(4))
   ts(3,4) = t_pa_pb(m, n, Vsr12(3), Vsr12(4))
   
   ts(4,1) = -t_s_px(n, Vsr21(2))   ! (pz s) = -(s pz)
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
