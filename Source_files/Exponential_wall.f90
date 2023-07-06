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
! This module contains subroutines to deal with short-range classical potential ("exponential wall")

MODULE Exponential_wall
use Universal_constants
use Objects
use Coulomb, only : f_cut_L_C, d_f_cut_L_C
use ZBL_potential, only : ZBL_pot, d_ZBL_pot

implicit none
PRIVATE

public :: get_Exp_wall_s, d_Exp_wall_Pressure_s, d_Exp_wall_pot_s
public :: get_short_range_rep_s, d_Short_range_pot_s, d_Short_range_Pressure_s

 contains

! New, general repulsive potential:
subroutine get_short_range_rep_s(TB_Expwall, Scell, NSC, matter, numpar, a)   ! Short-range potential energy
   type(TB_Short_Rep), dimension(:,:), intent(in) :: TB_Expwall    ! General short-range repulsive parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in), target :: matter   ! all material parameters
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   real(8), intent(out) :: a  ! [eV] short-range energy
   !------------------------
   real(8) :: sum_a, Short_range_pot
   integer(4) :: i, atom_2
   real(8), pointer :: a_r, Z1, Z2
   integer, pointer :: j, KOA1, KOA2, m

   sum_a = 0.0d0
   !$omp PARALLEL private(i, m, KOA1, Z1, atom_2, j, KOA2, Z2, a_r, Short_range_pot)
   !$omp do reduction( + : sum_a)
   do i = 1, Scell(NSC)%Na ! all atoms
      m => Scell(NSC)%Near_neighbor_size(i)
      KOA1 => Scell(NSC)%MDatoms(i)%KOA
      Z1 => matter%Atoms(KOA1)%Z
      do atom_2 = 1,m ! do only for atoms close to that one
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            KOA2 => Scell(NSC)%MDatoms(j)%KOA
            Z2 => matter%Atoms(KOA2)%Z
            a_r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)	! at this distance, R
            Short_range_pot = Shortrange_pot(TB_Expwall(KOA1,KOA2), a_r, Z1, Z2)    ! function below
            sum_a = sum_a + Short_range_pot
!             print*, 'get_short_range_rep_s', sum_a, Short_range_pot, a_r, KOA1, KOA2
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp end parallel
   a = sum_a*0.5d0   ! [eV], factor to compensate for double-counting
   nullify(a_r, j, m, KOA1, KOA2, Z1, Z2)
end subroutine get_short_range_rep_s


function Shortrange_pot(TB_Expwall, a_r, Z1, Z2) result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   type(TB_Short_Rep), intent(in), target :: TB_Expwall  ! Exponential wall parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   real(8), intent(in) :: Z1, Z2 ! atomic numbers of elements 1 and 2 (for ZBL)
   !------------------------
   real(8) :: f_cut_large, f_pow, f_exp, f_invexp, f_ZBL

   Pot = 0.0d0 ! to start with
   if (a_r < TB_Expwall%f_cut%d0 + TB_Expwall%f_cut%dd*10.0d0) then ! only at close range

      ! Cut-off function:
      f_cut_large = f_cut_L_C(a_r, TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd)   ! module "Coulomb"

      ! Contribution of the exponential function:
      f_exp = exp_function(a_r, TB_Expwall%f_exp%use_it, TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a)  ! below

      ! Contribution of the inverse exponential function:
      f_invexp = inv_exp_function(a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0)  ! below

      ! Contribution of the power functions:
      f_pow = power_function(a_r, TB_Expwall%f_pow)   ! below

      ! Contribution of the ZBL potential:
      if (TB_Expwall%f_ZBL%use_it) then
         f_ZBL = ZBL_pot(Z1, Z2, a_r)   ! module "ZBL_potential"
      else
         f_ZBL = 0.0d0
      endif

      !print*, 'Shortrange_pot-1:', a_r, f_invexp, f_exp, f_pow, f_ZBL, f_cut_large
      !print*, 'Shortrange_pot-2:', a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0

      ! Combine all:
      Pot = f_invexp + f_exp + f_pow + f_ZBL
      ! Augment with the cut-off function:
      Pot = Pot * f_cut_large  ! [eV]
   endif
end function Shortrange_pot




function d_Short_range_pot(TB_Expwall, a_r, Z1, Z2) result(dPot)
   real(8) :: dPot   ! derivative of the short-range potential
   type(TB_Short_Rep), intent(in), target :: TB_Expwall  ! short-range parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   real(8), intent(in) :: Z1, Z2 ! atomic numbers of elements 1 and 2 (for ZBL)
   !------------------------
   real(8) :: f_cut_large, d_f_large, f_exp, d_f_exp, f_invexp, d_f_invexp, f_pow, d_f_pow, Pot, d_Pot, f_ZBL, d_f_ZBL

   dPot = 0.0d0   ! to start with
   if (a_r < TB_Expwall%f_cut%d0 + TB_Expwall%f_cut%dd*10.0d0) then ! only at close range
      Pot = 0.0d0 ! to start with
      d_Pot = 0.0d0  ! to start with

      !------------------
      ! Cut-off function:
      f_cut_large = f_cut_L_C(a_r, TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd)   ! module "Coulomb"
      ! And its derivative:
      d_f_large = d_f_cut_L_C(a_r, TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd)   ! module "Coulomb"

      !------------------
      ! Contribution of the exponential function:
      f_exp = exp_function(a_r, TB_Expwall%f_exp%use_it, TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a)  ! below
      ! And its derivative:
      d_f_exp = d_exp_function(a_r, TB_Expwall%f_exp%use_it, TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a)  ! below

      !------------------
      ! Contribution of the inverse exponential function:
      f_invexp = inv_exp_function(a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0)  ! below
      ! And its derivative:
      d_f_invexp = d_inv_exp_function(a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0)  ! below

      !------------------
      ! Contribution of the power functions:
      f_pow = power_function(a_r, TB_Expwall%f_pow)   ! below
      ! And its derivative:
      d_f_pow = d_power_function(a_r, TB_Expwall%f_pow)   ! below

      !------------------
      ! Contribution of the ZBL potential:
      if (TB_Expwall%f_ZBL%use_it) then
         f_ZBL = ZBL_pot(Z1, Z2, a_r)   ! module "ZBL_potential"
         d_f_ZBL = d_ZBL_pot(Z1, Z2, a_r)   ! module "ZBL_potential"
      else
         f_ZBL = 0.0d0
         d_f_ZBL = 0.0d0
      endif

      !------------------
      ! Combine all:
      Pot = f_invexp + f_exp + f_pow + f_ZBL
      d_Pot = d_f_invexp + d_f_exp + d_f_pow + d_f_ZBL

      !------------------
      ! Augment the potential with the cut-off function (and its derivative):
      dPot = Pot*d_f_large + d_Pot*f_cut_large
   endif
end function d_Short_range_pot


! Derivatives of the Short-range energy by s:
subroutine d_Short_range_pot_s(Scell, NSC, matter, TB_Expwall, numpar)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in), target :: matter   ! all material parameters
   type(TB_Short_Rep), dimension(:,:),intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !------------------------
   real(8), dimension(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, ian, dik, djk, n, atom_2
   real(8), dimension(:,:), allocatable :: Erx_s
   real(8), pointer ::  x, y, z, Z1, Z2
   integer, pointer :: KOA1, KOA2, m, j1

   n = Scell(NSC)%Na ! number of atoms
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

   !$omp PARALLEL private(ian, i1, dik, dpsi, m, KOA1, Z1, atom_2, j1, djk, KOA2, Z2, x,y,z, x1, b, a_r,ddlta,b_delta)
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
         Z1 => matter%Atoms(KOA1)%Z ! atomic number
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
                  Z2 => matter%Atoms(KOA2)%Z ! atomic number
                  x => Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X, Y, Z
                  y => Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z => Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z

                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
                  x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                  x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

                  a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)   ! at this distance, R
                  b = d_Short_range_pot(TB_Expwall(KOA1,KOA2), a_r, Z1, Z2)  ! function above
                  !print*, KOA1, KOA2, TB_Expwall(KOA1,KOA2)%Param, TB_Expwall(KOA1,KOA2)%f_ZBL%use_it, TB_Expwall(KOA1,KOA2)%f_pow%use_it

                  ddlta = dble(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(:) = dpsi(:) + b_delta*x1(:)
               endif cos_if
            endif ! j1 > 0
         enddo ! j1

         Erx_s(:,ian) = Erx_s(:,ian) + dpsi(:) ! potential part in X-coordinate
      enddo ! i1
      ! Add extra short-range force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)*0.5d0
      ! factor 0.5 to compensate for double-counting
   enddo ! ian
   !$omp end do
   !$omp end parallel

   !pause 'd_Short_range_pot_s'

   deallocate(Erx_s)
   nullify(j1, m, KOA1, KOA2, x, y, z, Z1, Z2)
END subroutine d_Short_range_pot_s



! Derivatives of the short-range energy by h:
subroutine d_Short_range_Pressure_s(Scell, NSC, TB_Expwall, matter, numpar)
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Short_Rep), dimension(:,:),intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Solid), intent(in), target :: matter   ! all material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !------------------------
   real(8), dimension(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i, k, l, n, atom_2
   integer, pointer :: KOA1, KOA2, m, j
   real(8), pointer :: Z1, Z2
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
         Z1 => matter%Atoms(KOA1)%Z ! atomic number
         do atom_2 = 1,m		! do only for atoms close to that one
            j => Scell(NSC)%Near_neighbor_list(i, atom_2)	! this is the list of such close atoms
            if (j > 0) then
               KOA2 => Scell(NSC)%MDatoms(j)%KOA
               Z2 => matter%Atoms(KOA2)%Z ! atomic number
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
               dpsy = d_Short_range_pot(TB_Expwall(KOA1,KOA2), r, Z1, Z2)	! function above

               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! (j > 0)
         enddo ! atom_2

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)*0.5d0   ! factor 0.5 to compensate for double-counting
            enddo ! l
         enddo ! k
      enddo ! i
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add extra short-range part to existing TB part
   endif
   nullify(KOA1, KOA2, m, j, Z1, Z2)
end subroutine d_Short_range_Pressure_s





pure function power_function(r, f_pow) result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   real(8), intent(in) :: r   ! [A] distance
   type(Rep_pow), dimension(:), allocatable, intent(in) :: f_pow
   integer :: i, Nsiz
   Pot = 0.0d0
   if (allocated(f_pow)) then ! and only then
      Nsiz = size(f_pow)
      do i = 1, Nsiz
         Pot = Pot + f_pow(i)%Phi * (r/f_pow(i)%r0)**f_pow(i)%m
      enddo
   endif
end function power_function

pure function d_power_function(r, f_pow) result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   real(8), intent(in) :: r   ! [A] distance
   type(Rep_pow), dimension(:), allocatable, intent(in) :: f_pow
   integer :: i, Nsiz
   Pot = 0.0d0
   if (allocated(f_pow)) then ! and only then
      Nsiz = size(f_pow)
      do i = 1, Nsiz
         Pot = Pot + f_pow(i)%Phi * (r/f_pow(i)%r0)**f_pow(i)%m * f_pow(i)%m / r
      enddo
   endif
end function d_power_function


pure function inv_exp_function(r, use_it, C, r0)  result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   logical, intent(in) :: use_it ! flag to use this function or not
   real(8), intent(in) :: r   ! [A] distance
   real(8), intent(in) :: C    ! [eV] energy of the "wall"
   real(8), intent(in) :: r0     ! [A] "wall" position
   real(8) :: eps
   Pot = 0.0d0
   if (use_it) then  ! and only then
      Pot = Wall_potential(C, r0, r)  ! below
   endif
end function inv_exp_function

pure function d_inv_exp_function(r, use_it, C, r0)  result(Pot)
   real(8) Pot ! [eV] Derivative of the repulsive potential
   logical, intent(in) :: use_it ! flag to use this function or not
   real(8), intent(in) :: r   ! [A] distance
   real(8), intent(in) :: C    ! [eV] energy of the "wall"
   real(8), intent(in) :: r0     ! [A] "wall" position
   real(8) :: eps
   Pot = 0.0d0
   if (use_it) then  ! and only then
      Pot = d_Wall_potential(C, r0, r)  ! below
   endif
end function d_inv_exp_function


pure function exp_function(r, use_it, Phi, r0, a)  result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   logical, intent(in) :: use_it ! flag to use this function or not
   real(8), intent(in) :: r   ! [A] distance
   real(8), intent(in) :: Phi    ! [eV] energy of the "wall"
   real(8), intent(in) :: r0     ! [A] "wall" position
   real(8), intent(in) :: a  ! [A] width
   real(8) :: eps
   Pot = 0.0d0
   if (use_it) then  ! and only then
      eps = 1.0d-5
      if (abs(a) > eps) then
         Pot = Phi * exp(-(r - r0)/a)
      endif
   endif
end function exp_function

pure function d_exp_function(r, use_it, Phi, r0, a)  result(Pot)
   real(8) Pot ! [eV] Derivative of the repulsive potential
   logical, intent(in) :: use_it ! flag to use this function or not
   real(8), intent(in) :: r   ! [A] distance
   real(8), intent(in) :: Phi    ! [eV] energy of the "wall"
   real(8), intent(in) :: r0     ! [A] "wall" position
   real(8), intent(in) :: a  ! [A] width
   real(8) :: eps
   Pot = 0.0d0
   if (use_it) then  ! and only then
      eps = 1.0d-5
      if (abs(a) > eps) then
         Pot = -Phi * exp(-(r - r0)/a) / a
      endif
   endif
end function d_exp_function


!------------------------------------------------------
! Obsolete simple-wall (inverse exponential) potential:
subroutine get_Exp_wall_s(TB_Expwall, Scell, NSC, numpar, a)   ! Exponential wall energy
   type(TB_Exp_wall_simple), dimension(:,:), intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a	! [eV] exponential "wall" energy
   !----------------------
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
   !-------------------
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
   !---------------------
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
