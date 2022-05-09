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
! This module contains subroutines to deal with van der Waals forces within TB

MODULE Exponential_wall
use Universal_constants
use Objects
use Variables
use Little_subroutines
use Atomic_tools
use Coulomb, only : f_cut_L_C, d_f_cut_L_C

implicit none

 contains


subroutine get_Exp_wall_s(TB_Expwall, Scell, NSC, numpar, a)   ! Exponential wall energy
   type(TB_Exp_wall_simple), dimension(:,:), intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a	! [eV] exponential "wall" energy
   !=====================================================
   real(8), pointer :: a_r
   real(8) :: sum_a, Expwall_pot
   INTEGER(4) i, atom_2
   integer, pointer :: j, KOA1, KOA2, m
   
   sum_a = 0.0d0
   !$omp PARALLEL private(i, m, KOA1, atom_2, j, KOA2, a_r, Expwall_pot)
   !$omp do reduction( + : sum_a)
   do i = 1, Scell(NSC)%Na ! all atoms
      m => Scell(NSC)%Near_neighbor_size(i)
      KOA1 => Scell(NSC)%MDatoms(i)%KOA
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            KOA2 => Scell(NSC)%MDatoms(j)%KOA
            a_r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)	! at this distance, R
            Expwall_pot = Exp_wall_pot(TB_Expwall(KOA1,KOA2), a_r)    ! function below
            sum_a = sum_a + Expwall_pot
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp end parallel
   a = sum_a*0.5d0	! [eV], factor to compensate for double-counting
   nullify(a_r, j, m, KOA1, KOA2)
end subroutine get_Exp_wall_s


function Exp_wall_pot(TB_Expwall, a_r) result(Pot)
   real(8) Pot	! [eV] Exponential wall potential
   type(TB_Exp_wall_simple), intent(in), target :: TB_Expwall	! Exponential wall parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !====================================
   real(8) :: f_cut_large
   real(8), pointer :: C	! [eV]
   real(8), pointer :: r0	! [A]
   real(8), pointer :: d0	! [A]
   real(8), pointer :: dd	! [A]
   C => TB_Expwall%C 	! [eV]
   r0 => TB_Expwall%r0	! [A]
   d0 => TB_Expwall%d0	! [A]
   dd => TB_Expwall%dd	! [A]
   if (a_r > d0+dd*10.0d0) then ! anything beyond cut-offs is zero:
      Pot = 0.0d0
   else ! at shorter distances we use proper potential:
      f_cut_large = f_cut_L_C(a_r, d0, dd)	! module "Coulomb"
      Pot = Wall_potential(C, r0, a_r) * f_cut_large	! ! [eV]
   endif
   nullify(C, r0, d0, dd)
end function Exp_wall_pot


pure function Wall_potential(C, r0, r) ! Exponential wall potential
   real(8) :: Wall_potential
   real(8), intent(in) :: C, r0, r
   real(8) arg
   arg = r - r0
   if (abs(arg) <= log(TINY(r))) then	! atoms are too close
      Wall_potential = sign(1.0d25,C)	! "infinity"
   else	! proper potential:
      Wall_potential = C*exp(1.0d0/arg) ! [eV]
   endif
end function Wall_potential


pure function d_Wall_potential(C, r0, r)  ! derivative of exponential wall potential
   real(8) :: d_Wall_potential
   real(8), intent(in) :: C, r0, r
   real(8) :: rr0, arg
   arg = r - r0
   if (abs(arg) <= log(TINY(r))) then	! atoms are too close
      d_Wall_potential = 0.0d0	! no contribution to force at "infinite" potential
   else	! proper potential:
      rr0 = 1.0d0/arg
      d_Wall_potential = -C*exp(rr0)*rr0*rr0
   endif
end function d_Wall_potential




function d_Exp_wall_pot(TB_Expwall, a_r) result(dPot)
   real(8) :: dPot		! derivative of the exponential wall potential
   type(TB_Exp_wall_simple), intent(in), target :: TB_Expwall	! Exponential wall parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !====================================
   real(8) :: E_C, f_cut_large, d_f_large, d_E_C
   real(8), pointer :: C	! [eV]
   real(8), pointer :: r0	! [A]
   real(8), pointer :: d0	! [A]
   real(8), pointer :: dd	! [A]
   C => TB_Expwall%C 	! [eV]
   r0 => TB_Expwall%r0	! [A]
   d0 => TB_Expwall%d0	! [A]
   dd => TB_Expwall%dd	! [A]

   if (a_r > d0+dd*10.0d0) then ! anything beyond cut-offs is zero:
      dPot = 0.0d0
   else ! at shorter distances we use proper potential:
      f_cut_large = f_cut_L_C(a_r, d0, dd)	! module "Coulomb"
      E_C = Wall_potential(C, r0, a_r)	! function above
      d_E_C = d_Wall_potential(C, r0, a_r)		! derivative of  exponential wall part of the potential
      d_f_large = d_f_cut_L_C(a_r, d0, dd)	! module "Coulomb"
      dPot = d_E_C*f_cut_large + E_C*d_f_large
   endif
   nullify(C, r0, d0, dd)
end function d_Exp_wall_pot


! Derivatives of the exponential "wall" energy by s:
subroutine d_Exp_wall_pot_s(Scell, NSC, TB_Expwall, numpar) 
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Exp_wall_simple), dimension(:,:),intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !---------------------------------------
   real(8), dimension(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, ian, dik, djk, n, atom_2
   real(8), dimension(:,:), allocatable :: Erx_s
   integer, pointer :: KOA1, KOA2, m, j1
   real(8), pointer ::  x, y, z
   n = Scell(NSC)%Na ! number of atoms
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

   !$omp PARALLEL private(ian, i1, dik, dpsi, m, KOA1, atom_2, j1, djk, KOA2, x,y,z, x1, b, a_r,ddlta,b_delta)
   !$omp DO
   do ian = 1, n	! Forces for all atoms
      do i1 = 1, n	! contribution from all atoms
         if (ian == i1) then	! Kroniker delta
            dik = 1
         else
            dik = 0
         endif
         dpsi = 0.0d0
         m => Scell(NSC)%Near_neighbor_size(i1)
         KOA1 => Scell(NSC)%MDatoms(i1)%KOA
         do atom_2 = 1,m		! do only for atoms close to that one  
            j1 => Scell(NSC)%Near_neighbor_list(i1, atom_2)	! this is the list of such close atoms
            if (j1 > 0) then
               if (ian == j1) then	! Kroniker delta
                  djk = 1
               else
                  djk = 0
               endif
               cos_if:if ((dik-djk) /= 0) then ! without it, it gives ERROR
                  KOA2 => Scell(NSC)%MDatoms(j1)%KOA
                  x => Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X, Y, Z
                  y => Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z => Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z
                  
                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
                  x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                  x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
                  
                  a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R
                  b = d_Exp_wall_pot(TB_Expwall(KOA1,KOA2), a_r)	! function above

                  ddlta = dble(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(:) = dpsi(:) + b_delta*x1(:)
               endif cos_if
            endif ! j1 > 0
         enddo ! j1
         
         Erx_s(:,ian) = Erx_s(:,ian) + dpsi(:) ! potential part in X-coordinate
      enddo ! i1
      ! Add exponential wall force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)*0.5d0	! factor 0.5 to compensate for double-counting
   enddo ! ian
   !$omp end do
   !$omp end parallel

   deallocate(Erx_s)
   nullify(j1, m, KOA1, KOA2, x, y, z)
END subroutine d_Exp_wall_pot_s



! Derivatives of the exponential wall energy by h:
subroutine d_Exp_wall_Pressure_s(Scell, NSC, TB_Expwall, numpar)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Exp_wall_simple), dimension(:,:),intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !===============================================
   REAL(8), DIMENSION(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i, k, l, n, atom_2
   integer, pointer :: KOA1, KOA2, m, j
   real(8) r, rcur(3), scur(3), PForce(3,3)
   real(8) df_psy, psi, dpsy

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = Scell(NSC)%Na ! number of atoms

      PForce = 0.0d0 ! to start with
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         m => Scell(NSC)%Near_neighbor_size(i)
         KOA1 => Scell(NSC)%MDatoms(i)%KOA
         do atom_2 = 1,m		! do only for atoms close to that one  
            j => Scell(NSC)%Near_neighbor_list(i, atom_2)	! this is the list of such close atoms
            if (j > 0) then
               KOA2 => Scell(NSC)%MDatoms(j)%KOA
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
               dpsy = d_Exp_wall_pot(TB_Expwall(KOA1,KOA2), r)	! function above

               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! (j > 0)
         enddo ! atom_2

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)*0.5d0	! factor 0.5 to compensate for double-counting
            enddo ! l
         enddo ! k
      enddo ! i
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add "exponential wall" part to existing TB part
   endif
   nullify(KOA1, KOA2, m, j)
end subroutine d_Exp_wall_Pressure_s


END MODULE Exponential_wall
