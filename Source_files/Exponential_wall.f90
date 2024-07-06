! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT-3
! available at: https://zenodo.org/doi/10.5281/zenodo.8392569
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
! This module contains subroutines to deal with short-range classical potential ("exponential wall")

MODULE Exponential_wall
use Universal_constants
use Objects
use Coulomb, only : f_cut_L_C, d_f_cut_L_C, d2_f_cut_L_C, ddija_dria
use ZBL_potential, only : ZBL_pot, d_ZBL_pot, d2_ZBL_pot
use Little_subroutines, only : Find_monotonous_LE, linear_interpolation
use Algebra_tools, only : cubic_function, d_cubic_function, d2_cubic_function

#ifdef MPI_USED
use MPI_subroutines, only : do_MPI_Allreduce
#endif


implicit none
PRIVATE

public :: get_Exp_wall_s, d_Exp_wall_Pressure_s, d_Exp_wall_pot_s, d_Exponential_wall_forces
public :: get_short_range_rep_s, d_Short_range_pot_s, d_Short_range_Pressure_s

 contains

! New, general repulsive potential:
subroutine get_short_range_rep_s(TB_Expwall, Scell, NSC, matter, numpar, a)   ! Short-range potential energy
   type(TB_Short_Rep), dimension(:,:), intent(in) :: TB_Expwall    ! General short-range repulsive parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in), target :: matter   ! all material parameters
   type(Numerics_param), intent(inout) :: numpar   ! all numerical parameters
   real(8), intent(out) :: a  ! [eV] short-range energy
   !------------------------
   real(8) :: sum_a, Short_range_pot, E_pot
   integer(4) :: i, atom_2
   real(8), pointer :: a_r, Z1, Z2
   integer, pointer :: j, KOA1, KOA2, m
   real(8), dimension(Scell(1)%Na) :: E_pot_array
   integer :: Nstart, Nend, N_incr
   character(100) :: error_part

   sum_a = 0.0d0


#ifdef MPI_USED   ! only does anything if the code is compiled with MPI

   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank
   Nend = Scell(NSC)%Na
   E_pot_array = 0.0d0  ! initialize
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, Scell(NSC)%Na ! all atoms
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
            sum_a = sum_a + Short_range_pot  ! total
            E_pot_array(i) = E_pot_array(i) + Short_range_pot  ! for one atom
!             print*, 'get_short_range_rep_s', sum_a, Short_range_pot, a_r, KOA1, KOA2
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_short_range_rep_s'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{sum_a}', sum_a) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{E_pot_array}', E_pot_array) ! module "MPI_subroutines"
   ! And save for each atom:
   Scell(NSC)%MDAtoms(:)%Epot = Scell(NSC)%MDAtoms(:)%Epot + E_pot_array(:)*0.5d0 ! to exclude double-counting

#else ! use OpenMP instead
   !$omp PARALLEL private(i, m, KOA1, Z1, atom_2, j, KOA2, Z2, a_r, Short_range_pot, E_pot)
   !$omp do reduction( + : sum_a)
   do i = 1, Scell(NSC)%Na ! all atoms
      E_pot = 0.0d0 ! to start
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
            sum_a = sum_a + Short_range_pot  ! total
            E_pot = E_pot + Short_range_pot  ! for one atom
!             print*, 'get_short_range_rep_s', sum_a, Short_range_pot, a_r, KOA1, KOA2
         endif ! (j .GT. 0)
      enddo ! j
      Scell(NSC)%MDAtoms(i)%Epot = Scell(NSC)%MDAtoms(i)%Epot + E_pot*0.5d0 ! to exclude double-counting
   enddo ! i
   !$omp end do
   !$omp end parallel
#endif
   a = sum_a*0.5d0   ! [eV], factor to compensate for double-counting
   nullify(a_r, j, m, KOA1, KOA2, Z1, Z2)
end subroutine get_short_range_rep_s


function Shortrange_pot(TB_Expwall, a_r, Z1, Z2) result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   type(TB_Short_Rep), intent(in), target :: TB_Expwall  ! Exponential wall parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   real(8), intent(in) :: Z1, Z2 ! atomic numbers of elements 1 and 2 (for ZBL)
   !------------------------
   real(8) :: f_cut_large, f_pow, f_exp, f_invexp, f_ZBL, f_tab

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

      ! Contribution of tabulated potential:
      f_tab = tabulated_potential(a_r, TB_Expwall%f_tab) ! below


      !print*, 'Shortrange_pot-1:', a_r, f_invexp, f_exp, f_pow, f_ZBL, f_cut_large
      !print*, 'Shortrange_pot-2:', a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0

      ! Combine all:
      Pot = f_invexp + f_exp + f_pow + f_ZBL + f_tab
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
   real(8) :: f_tab, d_f_tab

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
      ! Contribution of tabulated potential:
      f_tab = tabulated_potential(a_r, TB_Expwall%f_tab) ! below
      d_f_tab = d_tabulated_potential(a_r, TB_Expwall%f_tab) ! below


      !------------------
      ! Combine all:
      Pot   = f_invexp + f_exp + f_pow + f_ZBL + f_tab
      d_Pot = d_f_invexp + d_f_exp + d_f_pow + d_f_ZBL + d_f_tab

      !------------------
      ! Augment the potential with the cut-off function (and its derivative):
      dPot = Pot*d_f_large + d_Pot*f_cut_large
   endif
end function d_Short_range_pot



function d2_Short_range_pot(TB_Expwall, a_r, Z1, Z2) result(dPot)
   real(8) :: dPot   ! second derivative of the short-range potential
   type(TB_Short_Rep), intent(in), target :: TB_Expwall  ! short-range parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   real(8), intent(in) :: Z1, Z2 ! atomic numbers of elements 1 and 2 (for ZBL)
   !------------------------
   real(8) :: f_cut_large, d_f_large, f_exp, d_f_exp, f_invexp, d_f_invexp, f_pow, d_f_pow, Pot, d_Pot, f_ZBL, d_f_ZBL
   real(8) :: f_tab, d_f_tab, d2_f_large, d2_Pot, d2_f_exp, d2_f_invexp, d2_f_pow, d2_f_ZBL, d2_f_tab

   dPot = 0.0d0   ! to start with
   if (a_r < TB_Expwall%f_cut%d0 + TB_Expwall%f_cut%dd*10.0d0) then ! only at close range
      Pot = 0.0d0 ! to start with
      d_Pot = 0.0d0  ! to start with
      d2_Pot = 0.0d0  ! to start with

      !------------------
      ! Cut-off function:
      f_cut_large = f_cut_L_C(a_r, TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd)   ! module "Coulomb"
      ! its derivative:
      d_f_large = d_f_cut_L_C(a_r, TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd)   ! module "Coulomb"
      ! and second derivative:
      d2_f_large = d2_f_cut_L_C(a_r, TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd)   ! module "Coulomb"

      !------------------
      ! Contribution of the exponential function:
      f_exp = exp_function(a_r, TB_Expwall%f_exp%use_it, TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a)  ! below
      ! its derivative:
      d_f_exp = d_exp_function(a_r, TB_Expwall%f_exp%use_it, TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a)  ! below
      ! and second derivative:
      d2_f_exp = d2_exp_function(a_r, TB_Expwall%f_exp%use_it, TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a)  ! below

      !------------------
      ! Contribution of the inverse exponential function:
      f_invexp = inv_exp_function(a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0)  ! below
      ! its derivative:
      d_f_invexp = d_inv_exp_function(a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0)  ! below
      ! and second derivative:
      d2_f_invexp = d2_inv_exp_function(a_r, TB_Expwall%f_inv_exp%use_it, TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0)  ! below

      !------------------
      ! Contribution of the power functions:
      f_pow = power_function(a_r, TB_Expwall%f_pow)   ! below
      ! its derivative:
      d_f_pow = d_power_function(a_r, TB_Expwall%f_pow)   ! below
      ! and second derivative:
      d2_f_pow = d2_power_function(a_r, TB_Expwall%f_pow)   ! below

      !------------------
      ! Contribution of the ZBL potential:
      if (TB_Expwall%f_ZBL%use_it) then
         f_ZBL = ZBL_pot(Z1, Z2, a_r)   ! module "ZBL_potential"
         ! its derivative:
         d_f_ZBL = d_ZBL_pot(Z1, Z2, a_r)   ! module "ZBL_potential"
         ! and second derivative:
         d2_f_ZBL = d2_ZBL_pot(Z1, Z2, a_r)   ! module "ZBL_potential"
      else
         f_ZBL = 0.0d0
         d_f_ZBL = 0.0d0
         d2_f_ZBL = 0.0d0
      endif

      !------------------
      ! Contribution of tabulated potential:
      f_tab = tabulated_potential(a_r, TB_Expwall%f_tab) ! below
      ! its derivative:
      d_f_tab = d_tabulated_potential(a_r, TB_Expwall%f_tab) ! below
      ! and second derivative:
      d2_f_tab = d2_tabulated_potential(a_r, TB_Expwall%f_tab) ! below

      !------------------
      ! Combine all:
      Pot   = f_invexp + f_exp + f_pow + f_ZBL + f_tab
      d_Pot = d_f_invexp + d_f_exp + d_f_pow + d_f_ZBL + d_f_tab
      d2_Pot = d2_f_invexp + d2_f_exp + d2_f_pow + d2_f_ZBL + d2_f_tab

      !------------------
      ! Construct the second derivative:
      dPot = Pot*d2_f_large + 2.0d0*d_Pot*d_f_large + d2_Pot*f_cut_large
   endif
end function d2_Short_range_pot




subroutine d_Exponential_wall_forces(Scell, NSC, matter, numpar, F_wall, dF_wall) ! get derivatives of the exponential wall forces
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of super-cell
   type(solid), intent(in), target :: matter   ! materil parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), dimension(:,:), allocatable, intent(out) :: F_wall, dF_wall ! force and its derivative
   !----------------------------------
   integer :: i1, atom_2
   real(8) :: drdrx, drdry, drdrz, d2rdr2x, d2rdr2y, d2rdr2z, F_r, dF_r, F(3), dF(3), a_r
   integer, pointer :: m, KOA1, KOA2, j1
   real(8), pointer :: x, y, z, Z1, Z2


   ! Make sure the forces are allocated:
   if (.not.allocated(F_wall)) allocate(F_wall(3,Scell(NSC)%Na))
   if (.not.allocated(dF_wall)) allocate(dF_wall(3,Scell(NSC)%Na))
   F_wall(:,:) = 0.0d0	! just to start with, forces
   dF_wall(:,:) = 0.0d0	! just to start with, derivatives of forces


   ! Check if there is any vdW forces in this parameterization:
   if (.not.allocated(Scell(NSC)%TB_Expwall)) return ! nothing to do


   !$omp PARALLEL private(i1, m, KOA1, Z1, atom_2, j1, KOA2, Z2, x,y,z, a_r, drdrx, drdry, drdrz, d2rdr2x, d2rdr2y, d2rdr2z, F_r, dF_r, F, dF)
   !$omp DO
   do i1 = 1, Scell(NSC)%Na	! contribution from all atoms
      m => Scell(NSC)%Near_neighbor_size(i1)
      KOA1 => Scell(NSC)%MDatoms(i1)%KOA
      Z1 => matter%Atoms(KOA1)%Z ! atomic number
      F = 0.0d0   ! to restart
      dF = 0.0d0  ! to restart
      do atom_2 = 1,m		! do only for atoms close to that one
         j1 => Scell(NSC)%Near_neighbor_list(i1, atom_2)	! this is the list of such close atoms
         if (j1 > 0) then
            KOA2 => Scell(NSC)%MDatoms(j1)%KOA
            Z2 => matter%Atoms(KOA2)%Z ! atomic number

            x => Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X, Y, Z
            y => Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
            z => Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z
            a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R

            ! Derivatives d r_{i,j} / d r_{i,alpha}:
            drdrx = x/a_r
            drdry = y/a_r
            drdrz = z/a_r

            ! Second derivatives d2 r_{ij} / d r2_{i,alpha}:
            d2rdr2x = ddija_dria(x, a_r)  ! module "Coulomb"
            d2rdr2y = ddija_dria(y, a_r)  ! module "Coulomb"
            d2rdr2z = ddija_dria(z, a_r)  ! module "Coulomb"

            call get_exponential_wall_F_dF(Scell(NSC)%TB_Expwall(KOA1, KOA2), a_r, Z1, Z2, F_r, dF_r)   ! below

            ! Construct the force and derivative:
            F(1) = F(1) + F_r*drdrx
            F(2) = F(2) + F_r*drdry
            F(3) = F(3) + F_r*drdrz
            dF(1) = dF(1) + dF_r*drdrx + F_r*d2rdr2x
            dF(2) = dF(2) + dF_r*drdry + F_r*d2rdr2y
            dF(3) = dF(3) + dF_r*drdrz + F_r*d2rdr2z

            !if (abs(F_r) > 1.0d-6) write(*,'(a,i0,i0,f,f,f,f,f)') 'a', i1, j1, F_r, drdrx, dF_r, d2rdr2x, F(1)
         endif ! j1 > 0
      enddo ! j1
      ! And save for each atom:
      F_wall(:,i1) = F_wall(:,i1) + F
      dF_wall(:,i1) = dF_wall(:,i1) + dF
   enddo ! i1
   !$omp end do
   !$omp end parallel

   nullify(m, KOA1, KOA2, j1, x, y, z, Z1, Z2)
end subroutine d_Exponential_wall_forces


subroutine get_exponential_wall_F_dF(TB_Expwall, a_r, Z1, Z2, F_r, dF_r) ! wrapper around select_type; ASSOCIATE does not work inside OMP region
   class(TB_Exp_wall), intent(in) :: TB_Expwall
   real(8), intent(in) :: a_r, Z1, Z2
   real(8), intent(out) :: F_r, dF_r   ! force and its derivative
   !-----------------
   F_r = 0.0d0 ! to start with
   dF_r = 0.0d0 ! to start with

   select type (TB_Expwall)
   type is (TB_Exp_wall_simple)
      F_r = d_Exp_wall_pot(TB_Expwall, a_r) ! below
      dF_r = d2_Exp_wall_pot(TB_Expwall, a_r) ! below
      !if (abs(F_r) > 1.0d-6) print*, 'd_Exp_wall_pot', F_r, dF_r

   type is (TB_Short_Rep)
      F_r = d_Short_range_pot(TB_Expwall, a_r, Z1, Z2)  ! below
      dF_r = d2_Short_range_pot(TB_Expwall, a_r, Z1, Z2)  ! below
      !if (abs(F_r) > 1.0d-6) print*, 'd_Short_range_pot', F_r, dF_r

   end select
end subroutine get_exponential_wall_F_dF



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



!----------------------------------
! Functions:

function tabulated_potential(r, f_tab) result(Pot)
   real(8) Pot ! [eV] Repulsive potential
   real(8), intent(in) :: r   ! [A] distance
   type(Rep_tab), intent(in) :: f_tab
   !---------------------------
   integer :: i_closest, Nsiz, i
   real(8) :: E, x

   if (f_tab%use_it) then  ! and only then

      Nsiz = size(f_tab%R)
      if (r > f_tab%R(Nsiz)) then ! nullify potential beyond the grid:
         E = 0.0d0
      else  ! at short distances
         ! Find the value of radius in the array:
         call Find_monotonous_LE(f_tab%R, r, i_closest) ! module "Little_subroutines"

         if (f_tab%use_spline) then  ! use cubic spline instead of the finite difference

            ! For spline, only the interval N-1 is defined, so use it for too large distances:
            if (i_closest == Nsiz) i_closest = Nsiz - 1

            E = cubic_function(r-f_tab%R(i_closest), f_tab%a(i_closest), f_tab%b(i_closest), &
                               f_tab%c(i_closest), f_tab%d(i_closest))  ! module "Algebra_tools"

         else ! use finite diference
            i_closest = i_closest + 1  ! set upper point for interpolation between i-1 and i
            ! Interpolate between the grid points:
            call linear_interpolation(f_tab%R, f_tab%E, r, E, i_closest)   ! module "Little_subroutines"

         endif ! (f_tab%use_spline)
      endif ! (r > f_tab%R(Nsiz))
   else  ! no potential to use
      E = 0
   endif ! (f_tab%use_it)

   Pot = E ! [eV] potential
end function tabulated_potential


function d_tabulated_potential(r, f_tab) result(Pot)
   real(8) Pot ! [eV/A] Derivative of the repulsive potential
   real(8), intent(in) :: r   ! [A] distance
   type(Rep_tab), intent(in) :: f_tab
   !---------------------------
   integer :: i_closest, Nsiz
   real(8) :: E

   if (f_tab%use_it) then  ! and only then
      Nsiz = size(f_tab%R)

      ! Find the value of radius in the array:
      if (r > f_tab%R(Nsiz)) then ! nullify potential beyond the grid:
         E = 0.0d0
      else  ! (r > f_tab%R(Nsiz))
         call Find_monotonous_LE(f_tab%R, r, i_closest) ! module "Little_subroutines"
         if (i_closest == 1) then
            E = (f_tab%E(i_closest+1) - f_tab%E(i_closest)) / (f_tab%R(i_closest+1) - f_tab%R(i_closest))
         else
            if (f_tab%use_spline) then  ! use cubic spline instead of the finite difference

               ! For spline, only the interval N-1 is defined, so use it for too large distances:
               if (i_closest == Nsiz) i_closest = Nsiz - 1

               E = d_cubic_function(r-f_tab%R(i_closest), f_tab%b(i_closest), &
                               f_tab%c(i_closest), f_tab%d(i_closest))  ! module "Algebra_tools"

            else ! use finite diference
               i_closest = i_closest + 1  ! set upper point for interpolation between i-1 and i
               ! Numerical derivative of lineaar function:
               E = (f_tab%E(i_closest) - f_tab%E(i_closest-1)) / (f_tab%R(i_closest) - f_tab%R(i_closest-1))

            endif ! (f_tab%use_spline)
         endif ! (i_closest == 1)
      endif ! (f_tab%use_spline)
   else  ! no potenital => no force
      E = 0.0d0
   endif ! (f_tab%use_it)

   Pot = E ! [eV/A] derivative of the potential
end function d_tabulated_potential


function d2_tabulated_potential(r, f_tab) result(Pot)
   real(8) Pot ! [eV/A^2] second derivative of the repulsive potential
   real(8), intent(in) :: r   ! [A] distance
   type(Rep_tab), intent(in) :: f_tab
   !---------------------------
   integer :: i_closest, Nsiz
   real(8) :: E

   if (f_tab%use_it) then  ! and only then
      Nsiz = size(f_tab%R)

      ! Find the value of radius in the array:
      if (r > f_tab%R(Nsiz)) then ! nullify potential beyond the grid:
         E = 0.0d0
      else  ! (r > f_tab%R(Nsiz))
         call Find_monotonous_LE(f_tab%R, r, i_closest) ! module "Little_subroutines"
         if (i_closest == 1) then
            E = (f_tab%E(i_closest+1) - f_tab%E(i_closest)) / (f_tab%R(i_closest+1) - f_tab%R(i_closest))
         else
            if (i_closest == Nsiz) i_closest = Nsiz - 1
            if (f_tab%use_spline) then  ! use cubic spline instead of the finite difference
               ! For spline, only the interval N-1 is defined, so use it for too large distances:
               E = d2_cubic_function(r-f_tab%R(i_closest), f_tab%c(i_closest), f_tab%d(i_closest))  ! module "Algebra_tools"

            else ! use finite diference
               i_closest = i_closest + 1  ! set upper point for interpolation between i-1 and i
               ! Numerical second derivative:
               E = (f_tab%E(i_closest) - 2.0d0*f_tab%E(i_closest-1)+f_tab%E(i_closest-2)) / &
                   (f_tab%R(i_closest) - f_tab%R(i_closest-1))*(f_tab%R(i_closest-1) - f_tab%R(i_closest-2))

            endif ! (f_tab%use_spline)
         endif ! (i_closest == 1)
      endif ! (f_tab%use_spline)
   else  ! no potenital => no force
      E = 0.0d0
   endif ! (f_tab%use_it)

   Pot = E ! [eV/A^2] second derivative of the potential
end function d2_tabulated_potential


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
   real(8) Pot ! [eV/A] Repulsive potential
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

pure function d2_power_function(r, f_pow) result(Pot)
   real(8) Pot ! [eV/A^2] Repulsive potential
   real(8), intent(in) :: r   ! [A] distance
   type(Rep_pow), dimension(:), allocatable, intent(in) :: f_pow
   integer :: i, Nsiz
   Pot = 0.0d0
   if (allocated(f_pow)) then ! and only then
      Nsiz = size(f_pow)
      do i = 1, Nsiz
         Pot = Pot + f_pow(i)%Phi * (r/f_pow(i)%r0)**f_pow(i)%m * f_pow(i)%m*(f_pow(i)%m-1.0d0) / r**2
      enddo
   endif
end function d2_power_function



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
   real(8) Pot ! [eV/A] Derivative of the repulsive potential
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

pure function d2_inv_exp_function(r, use_it, C, r0)  result(Pot)
   real(8) Pot ! [eV/A^2] Derivative of the repulsive potential
   logical, intent(in) :: use_it ! flag to use this function or not
   real(8), intent(in) :: r   ! [A] distance
   real(8), intent(in) :: C    ! [eV] energy of the "wall"
   real(8), intent(in) :: r0     ! [A] "wall" position
   real(8) :: eps
   Pot = 0.0d0
   if (use_it) then  ! and only then
      Pot = d2_Wall_potential(C, r0, r)  ! below
   endif
end function d2_inv_exp_function



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
   real(8) Pot ! [eV/A] Derivative of the repulsive potential
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


pure function d2_exp_function(r, use_it, Phi, r0, a)  result(Pot)
   real(8) Pot ! [eV/A^2] Second derivative of the repulsive potential
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
         Pot = Phi * exp(-(r - r0)/a) / a**2
      endif
   endif
end function d2_exp_function


!------------------------------------------------------
! Obsolete simple-wall (inverse exponential) potential:
subroutine get_Exp_wall_s(TB_Expwall, Scell, NSC, numpar, a)   ! Exponential wall energy
   type(TB_Exp_wall_simple), dimension(:,:), intent(in) :: TB_Expwall	! Exponential wall parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a	! [eV] exponential "wall" energy
   !----------------------
   real(8), pointer :: a_r
   real(8) :: sum_a, Expwall_pot, E_pot
   real(8), dimension(Scell(1)%Na) :: E_pot_array
   INTEGER(4) i, atom_2
   integer, pointer :: j, KOA1, KOA2, m
   integer :: Nstart, Nend, N_incr
   character(100) :: error_part
   
   sum_a = 0.0d0


#ifdef MPI_USED   ! only does anything if the code is compiled with MPI

   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank
   Nend = Scell(NSC)%Na
   E_pot_array = 0.0d0  ! initializing

   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, Scell(NSC)%Na ! all atoms
      m => Scell(NSC)%Near_neighbor_size(i)
      KOA1 => Scell(NSC)%MDatoms(i)%KOA
      do atom_2 = 1,m ! do only for atoms close to that one
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            KOA2 => Scell(NSC)%MDatoms(j)%KOA
            a_r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)	! at this distance, R
            Expwall_pot = Exp_wall_pot(TB_Expwall(KOA1,KOA2), a_r)    ! function below
            sum_a = sum_a + Expwall_pot   ! total
            E_pot_array(i) = E_pot_array(i) + Expwall_pot   ! for one atom
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i

   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_Exp_wall_s'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{sum_a}', sum_a) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'{E_pot_array}', E_pot_array) ! module "MPI_subroutines"

   ! And save for each atom:
   Scell(NSC)%MDAtoms(:)%Epot = Scell(NSC)%MDAtoms(:)%Epot + E_pot_array(:) * 0.5d0   ! exclude double-counting

#else ! use OpenMP instead
   !$omp PARALLEL private(i, m, KOA1, atom_2, j, KOA2, a_r, Expwall_pot, E_pot)
   !$omp do reduction( + : sum_a)
   do i = 1, Scell(NSC)%Na ! all atoms
      E_pot = 0.0d0
      m => Scell(NSC)%Near_neighbor_size(i)
      KOA1 => Scell(NSC)%MDatoms(i)%KOA
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            KOA2 => Scell(NSC)%MDatoms(j)%KOA
            a_r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4)	! at this distance, R
            Expwall_pot = Exp_wall_pot(TB_Expwall(KOA1,KOA2), a_r)    ! function below
            sum_a = sum_a + Expwall_pot   ! total
            E_pot = E_pot + Expwall_pot   ! for one atom
         endif ! (j .GT. 0)
      enddo ! j
      ! And save for each atom:
      Scell(NSC)%MDAtoms(i)%Epot = Scell(NSC)%MDAtoms(i)%Epot + E_pot * 0.5d0   ! exclude double-counting
   enddo ! i
   !$omp end do
   !$omp end parallel
#endif

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


pure function d2_Wall_potential(C, r0, r) result(dF) ! second derivative of exponential wall potential
   real(8) :: dF
   real(8), intent(in) :: C, r0, r
   real(8) :: rr0, arg
   arg = r - r0
   if (abs(arg) <= log(TINY(r))) then	! atoms are too close
      dF = 0.0d0	! no contribution to force at "infinite" potential
   else	! proper potential:
      rr0 = 1.0d0/arg
      dF = C*exp(rr0)*rr0**3 * (rr0 + 2.0d0)
   endif
end function d2_Wall_potential




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
      E_C   = Wall_potential(C, r0, a_r)	! function above
      d_E_C = d_Wall_potential(C, r0, a_r)		! derivative of  exponential wall part of the potential

      f_cut_large = f_cut_L_C(a_r, d0, dd)	! module "Coulomb"
      d_f_large   = d_f_cut_L_C(a_r, d0, dd)	! module "Coulomb"

      dPot = d_E_C*f_cut_large + E_C*d_f_large
   endif
   nullify(C, r0, d0, dd)
end function d_Exp_wall_pot


function d2_Exp_wall_pot(TB_Expwall, a_r) result(dPot)
   real(8) :: dPot		! derivative of the exponential wall potential
   type(TB_Exp_wall_simple), intent(in), target :: TB_Expwall	! Exponential wall parameters
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !---------------------
   real(8) :: E_C, f_cut_large, d_f_large, d2_f_large, d_E_C, d2_E_C
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
      E_C    = Wall_potential(C, r0, a_r)       ! above
      d_E_C  = d_Wall_potential(C, r0, a_r)   ! derivative of exponential wall part of the potential
      d2_E_C = d2_Wall_potential(C, r0, a_r) ! second derivative of exponential wall part of the potential

      f_cut_large = f_cut_L_C(a_r, d0, dd)   ! module "Coulomb"
      d_f_large   = d_f_cut_L_C(a_r, d0, dd)   ! module "Coulomb"
      d2_f_large  = d2_f_cut_L_C(a_r, d0, dd) ! module "Coulomb"

      dPot = d2_E_C*f_cut_large + 2.0d0*d_E_C*d_f_large + E_C*d2_f_large
   endif
   nullify(C, r0, d0, dd)
end function d2_Exp_wall_pot


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
