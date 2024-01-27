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
! References used in the module:
! [1]  Jepps et al., PRE 62, 4757 (2000)
! This module includes subroutines to calculate thermodynamic properties from MD simulation:
MODULE Atomic_thermodynamics
use Universal_constants
use Objects

implicit none
PRIVATE

public :: partial_temperatures, get_temperature_from_equipartition, temperature_from_moments, get_temperature_from_distribution, &
get_temperature_from_entropy, atomic_entropy, set_Gibbs_x_powerDOS, set_Maxwell_distribution, temperature_from_moments_pot, &
Maxwell_entropy

real(8), parameter :: m_two_third = 2.0d0 / 3.0d0


 contains


subroutine partial_temperatures(Scell, matter, numpar)
   type(Super_cell), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in), target :: matter	! materil parameters
   type(Numerics_param), intent(in) :: numpar   ! numerical parameters
   !--------------------
   integer :: Nat, i
   real(8) :: prefac, Ekin(3)
   real(8), pointer :: Mass

   Nat = size(Scell%MDAtoms)  ! number of atoms

   ! Kinetic temperatures:
   Ekin(:) = 0.0d0   ! to start with
   prefac = 1d10/g_e ! to get [eV]
   do i = 1, Nat
      Mass => matter%Atoms(Scell%MDatoms(i)%KOA)%Ma
      Ekin(:) = Ekin(:) + 0.5d0*Mass*Scell%MDatoms(i)%V(:)*Scell%MDatoms(i)%V(:)
   enddo
   Ekin(:) = Ekin(:) * prefac / dble(Nat)    ! [eV]
   Scell%Ta_r_var(1:3) = 2.0d0 * Ekin(1:3) * g_kb     ! [K]

   ! Configurational temperatures:
   prefac = (Scell%V * 1e-30) / dble(Nat) / g_e * g_kb  ! to get temperature in [K]
   Scell%Ta_r_var(4) = -Scell%Pot_Stress(1,1) * prefac   ! X
   Scell%Ta_r_var(5) = -Scell%Pot_Stress(2,2) * prefac   ! Y
   Scell%Ta_r_var(6) = -Scell%Pot_Stress(3,3) * prefac   ! Z

   Scell%Ta_r_var(:) = abs(Scell%Ta_r_var(:))   ! ensure it is non-negative

   nullify(Mass)
end subroutine partial_temperatures


function get_temperature_from_equipartition(Scell, matter, numpar, non_periodic) result(Ta) ! works for non-periodic boundaries
   real(8) :: Ta  ! [K] configurational temperature
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in), target :: matter	! materil parameters
   type(Numerics_param), intent(in) :: numpar   ! numerical parameters
   logical, intent(in), optional :: non_periodic   ! if we want to use nonperiodic expression
   !------------------------------
   integer :: Nat, i
   real(8) :: F(3), acc(3), r(3), Pot, Pot_r(3), Pot_tot
   real(8), pointer :: Mass
   logical :: do_nonper

   if (present(non_periodic)) then
      do_nonper = non_periodic
   else
      do_nonper = .false.  ! by default, use periodic definition
   endif

   Nat = size(Scell%MDAtoms)  ! number of atoms

   if ( .not.do_nonper ) then ! periodic boundaries are used
     ! Get it from the pressure, calculated for the periodic boundaries:
     Ta = -Scell%Pot_Pressure * (Scell%V * 1e-30) / dble(Nat) / g_e   ! [eV]

   else ! for nonperiodic systems (it is more straightforward):
      Pot_tot = 0.0d0   ! to start with
      do i = 1, Nat  ! for all atoms
         ! Convert acceleration into SI units:
         acc(:) = Scell%MDAtoms(i)%accel(:) * 1.0d20 ! [A/fs^2] -> [m/s^2]
         Mass => matter%Atoms(Scell%MDatoms(i)%KOA)%Ma ! atomic mass [kg]

         ! Get the force:
         F(:) = Mass * acc(:) ! [N]
         ! Get the coordinate relative to the center of the supercell:
         r(:) = position_relative_to_center(Scell, i) ! below
         r(:) = r(:) * 1.0d-10   ! [A] -> [m]

         ! Construct the potential energy contribution:
         Pot = SUM(F(:) * r(:)) / g_e ! [eV]

         ! Total potential contribution to get the temperature
         Pot_tot = Pot_tot + Pot
      enddo
      ! Configurational temperature from the equipartition theorem as potential energy per atom per degree of freedom:
      Ta = -Pot_tot / (3.0d0 * dble(Nat))   ! [eV]
   endif

   ! Convert [eV] -> {K}:
   Ta = abs(Ta) * g_kb  ! ensure it is non-negative, even though pressure may be
end function get_temperature_from_equipartition


function position_relative_to_center(Scell, i_at) result(Rrc)
   real(8), dimension(3) :: Rrc  ! [A] position relative to the center of the supercell
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: i_at   ! atom index
   !-----------------
   real(8) :: Sj(3), a_r

   ! relative coords of the center of the supercell:
   Sj(:) = 0.5d0
   ! Shortest distance to the center in a cell with periodic boundaries:
   call shortest_distance_to_point(Scell, i_at, Sj, a_r, x1=Rrc(1), y1=Rrc(2), z1=Rrc(3))   ! below
end function position_relative_to_center


subroutine temperature_from_moments(Scell, Ta, E0)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   real(8), intent(out) :: Ta, E0   ! [K] and [eV] temperature and shift
   !---------------
   real(8) :: E1, E2, one_Nat, T_test
   integer :: Nat

   ! Number of atoms:
   Nat = size(Scell%MDAtoms(:))
   one_Nat = 1.0d0 / dble(Nat)

   ! First moment of the distribution:
   E1 = SUM( Scell%MDAtoms(:)%Ekin ) * one_Nat

   ! Second moment of the distribution:
   E2 = SUM( Scell%MDAtoms(:)%Ekin**2 ) * one_Nat

   ! Define temperature assuming Maxwell distribution:
   Ta = sqrt( m_two_third * (E2 - E1**2) )   ! [eV]

   ! Define the shift of the generalized maxwell distribution:
   E0 = E1 - 1.5d0*Ta   ! [eV]

   ! Convert [eV] -> [K]:
   Ta = Ta * g_kb ! [K]

   ! For maxwellian distribution, the result is iudentical to:
   ! T_test = get_T_from_fluctuation(E1, E2) ! below
   ! print*, 'temperature_from_moments', Ta, T_test* g_kb
end subroutine temperature_from_moments



subroutine temperature_from_n_moment(Scell, n, Ta, from_abs)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: n      ! moment of distribution to be used
   real(8), intent(out) :: Ta    ! [K] temperature from n-th moment of distribution
   logical, optional, intent(in) :: from_abs ! get temperature from the n-th moment instead of n-th order fluctuation
   !---------------
   real(8) :: E1, En, one_Nat, den
   integer :: Nat
   logical :: use_fluc

   if (present(from_abs)) then
      use_fluc = .not.from_abs
   else  ! by default, use fluctuation instead of absolute value
      use_fluc = .true.
   endif

   ! Number of atoms:
   Nat = size(Scell%MDAtoms(:))
   one_Nat = 1.0d0 / dble(Nat)

   ! First moment of the distribution:
   E1 = SUM( Scell%MDAtoms(:)%Ekin ) * one_Nat

   ! n-th moment of the distribution:
   En = SUM( Scell%MDAtoms(:)%Ekin**n ) * one_Nat

   ! Define temperature from n-th moment:
   if (use_fluc) then   ! n-th-order fluctuation:
      den = 2.0d0/g_sqrt_Pi * ( GAMMA(n+3.0d0/2.0d0) - (3.0d0/2.0d0)**n )
      Ta = ( (En - E1**n) / den ) ** (1.0d0/dble(n))  ! [eV]
   else ! n-th moment:
      den = 2.0d0/g_sqrt_Pi * ( GAMMA(n+3.0d0/2.0d0))
      Ta = (En / den) ** (1.0d0/dble(n))  ! [eV]
   endif

   ! Convert [eV] -> [K]:
   Ta = Ta * g_kb ! [K]
end subroutine temperature_from_n_moment



subroutine temperature_from_moments_pot(Scell, Ta, E0)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   real(8), intent(out) :: Ta, E0   ! [K] and [eV] potential temperature and shift
   !---------------
   real(8), dimension(size(Scell%MDAtoms)) :: E_pot
   real(8) :: E1, E2, one_Nat
   integer :: Nat

   ! Number of atoms:
   Nat = size(Scell%MDAtoms(:))
   one_Nat = 1.0d0 / dble(Nat)

   ! Set the potential energy for each atom:
   E_pot(:) = Scell%MDAtoms(:)%Epot !*0.5d0

   ! First moment of the distribution:
   E1 = SUM( E_pot(:) ) * one_Nat

   ! Second moment of the distribution:
   E2 = SUM( E_pot(:)**2 ) * one_Nat

   ! Define temperature assuming Maxwell distribution:
   Ta = sqrt( m_two_third * (E2 - E1**2) )   ! [eV]

   ! Define the shift of the generalized maxwell distribution:
   E0 = E1 - 1.5d0*Ta   ! [eV]

   ! Convert [eV] -> [K]:
   Ta = Ta * g_kb ! [K]
end subroutine temperature_from_moments_pot



pure function Maxwell_entropy(Ta) result(Sa) ! for equilibrium Maxwell distribution
   real(8) Sa
   real(8), intent(in) :: Ta  ! [eV]
   !------------------
   if (Ta > 0.0d0) then ! possible to get entropy
      !Sa = g_kb_EV * (log(g_sqrt_Pi * Ta) + (g_Eulers_gamma + 1.0d0)*0.5d0)   ! [eV/K]
      Sa = g_kb_EV * ( 2.5d0 - log(2.0d0/(g_sqrt_Pi * Ta**(1.5d0))) )   ! [eV/K]
   else  ! undefined
      Sa = 0.0d0
   endif
end function Maxwell_entropy


pure function get_temperature_from_entropy(Sa) result(Ta)
   real(8) Ta  ! [K]
   real(8), intent(in) :: Sa
   !--------------------
   !Ta = 1.0d0/g_sqrt_Pi * exp(Sa / g_kb_EV - 0.5d0*(g_Eulers_gamma + 1.0d0))  ! [eV]
   Ta = (4.0d0/g_Pi)**(1.0d0/3.0d0) * exp( Sa / g_kb_EV * 2.0d0/3.0d0 - 5.0d0/3.0d0 )  ! [eV]
   Ta = Ta * g_kb ! [K]
end function get_temperature_from_entropy



pure function get_temperature_from_distribution(E_grid, fa) result(Ta)
   real(8) Ta  ! [K]
   real(8), dimension(:), intent(in) :: E_grid, fa
   !--------------------
   real(8), dimension(size(fa)) :: dE
   real(8) :: Ekin
   integer :: i, j, Nsiz

   Nsiz = size(E_grid)
   Ekin = 0.0d0   ! to start with
   do j = 1, Nsiz-1
      if (j == 1) then
         dE(j) = E_grid(j+1) - E_grid(j)
      else
         dE(j) = E_grid(j) - E_grid(j-1)
      endif
      Ekin = Ekin + (E_grid(j+1) + E_grid(j))*0.5d0 * fa(j) * dE(j)
   enddo

   Ta = 2.0d0/3.0d0 * Ekin * g_kb ! [K]
end function get_temperature_from_distribution




subroutine atomic_entropy(E_grid, fa, Sa, i_start, i_end)
   real(8), dimension(:), intent(in) :: E_grid, fa ! atomic distribution function
   real(8), intent(out) :: Sa ! atomic entropy
   integer, intent(in), optional :: i_start, i_end  ! starting and ending levels to include
   ! Boltzmann H-function: Se = -kB * int [ ( f * (ln(f/sqrt(E)) -1) ) ]
   !----------------------------
   real(8), dimension(size(fa)) :: f_lnf, dE
   real(8) :: eps, E0
   integer :: i, Nsiz, i_low, i_high, j
   !============================
   eps = 1.0d-12  ! precision
   Nsiz = size(fa)

   if (present(i_start)) then
      i_low = i_start
   else  ! default, start from 1
      i_low = 1
   endif

   if (present(i_end)) then
      i_high = i_end
   else  ! default, end at the end
      i_high = Nsiz
   endif

   ! To start with:
   Sa = 0.0d0
   f_lnf = 0.0d0

   ! Set integration step as an array:
   do j = 1, Nsiz
      if (j == 1) then
         dE(j) = E_grid(j+1) - E_grid(j)
      else
         dE(j) = E_grid(j) - E_grid(j-1)
      endif
   enddo ! j

   ! Shift, if any:
   E0 = min(E_grid(1),0.0d0)

   ! Entropy (via Boltzmann H-theorem definition):
   !where (fa(i_low:i_high) > eps) f_lnf(i_low:i_high) = fa(i_low:i_high)*log(fa(i_low:i_high))
   where ( (fa(i_low:i_high) > eps) .and. (E_grid(i_low:i_high) > E0) ) f_lnf(i_low:i_high) = &
            fa(i_low:i_high)*(log(fa(i_low:i_high)/sqrt(E_grid(i_low:i_high)-E0)) - 1.0d0)
   Sa = SUM(f_lnf(i_low:i_high) * dE(i_low:i_high))

   ! Make proper units:
   Sa = -g_kb_EV*Sa  ! [eV/K]
end subroutine atomic_entropy


subroutine set_Maxwell_distribution(numpar, Scell, NSC)
   type(Numerics_param), intent(in) :: numpar   ! numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   !----------------------------------
   integer :: j, Nsiz
   real(8) :: arg, Tfact

   Nsiz = size(Scell(NSC)%Ea_grid) ! size of the energy grid
   if (.not.allocated(Scell(NSC)%fa_eq)) allocate(Scell(NSC)%fa_eq(Nsiz), source = 0.0d0)

   if (Scell(NSC)%TaeV > 0.0d0) then
      Tfact = (1.0d0 / Scell(NSC)%TaeV)**1.5d0
      do j = 1, Nsiz
         arg = Scell(NSC)%Ea_grid(j) / Scell(NSC)%TaeV
         Scell(NSC)%fa_eq(j) = 2.0d0 * sqrt(Scell(NSC)%Ea_grid(j) / g_Pi) * Tfact * exp(-arg)
      enddo
   else  ! zero-temperature distribution
      Scell(NSC)%fa_eq(:) = 0.0d0
      Scell(NSC)%fa_eq(1) = 1.0d0
   endif

   ! For printout:
   if (numpar%save_fa) then
      Nsiz = size(Scell(NSC)%Ea_grid_out) ! size of the energy grid
      if (.not.allocated(Scell(NSC)%fa_eq_out)) allocate(Scell(NSC)%fa_eq_out(Nsiz), source = 0.0d0)
      if (Scell(NSC)%TaeV > 0.0d0) then
         Tfact = (1.0d0 / Scell(NSC)%TaeV)**1.5d0
         do j = 1, Nsiz
            arg = Scell(NSC)%Ea_grid_out(j) / Scell(NSC)%TaeV
            Scell(NSC)%fa_eq_out(j) = 2.0d0 * sqrt(Scell(NSC)%Ea_grid_out(j) / g_Pi) * Tfact * exp(-arg)
         enddo
      else  ! zero-temperature distribution
         Scell(NSC)%fa_eq_out(:) = 0.0d0
         Scell(NSC)%fa_eq_out(1) = 1.0d0
      endif
   endif
end subroutine set_Maxwell_distribution



subroutine set_Gibbs_x_powerDOS(numpar, Scell, NSC) ! below
   type(Numerics_param), intent(in) :: numpar   ! numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   !----------------------------------
   integer :: j, Nsiz, Nat, Niter
   real(8) :: A, U0, T, b  ! Generealized distribution parameters
   real(8) :: arg, Tfact, E_shift, Ta, dE, U1, U2, U3, T_cur, eps, b_cur, U0_cur, T0
   real(8) :: alpha_T, alpha_U0, alpha_b

   ! Find the parameters iteratively:
   Nat = size(Scell(NSC)%MDAtoms)   ! number of atoms
   ! Define mean energy, square, and power-3:
   U1 = SUM(Scell(NSC)%MDAtoms(:)%Epot) / dble(Nat)
   U2 = SUM(Scell(NSC)%MDAtoms(:)%Epot**2) / dble(Nat)
   U3 = SUM(Scell(NSC)%MDAtoms(:)%Epot**3) / dble(Nat)
   !print*, 'Us:', U1, U2, U3

   ! Define parametres of the (Gibbs * DOS) distribution:
   U0 = minval(Scell(NSC)%MDAtoms(:)%Epot)   ! to start with
   ! Assume fixed DOS of harmonic oscillator just for comparisons and tests:
   b = numpar%power_b   ! Fixed, or starting, b for distribution: n(U)=A*exp(-(U-U0)/T)*U^b
   call find_Gibbs_x_powerDOS_parameters(U1, U2, U3, A, T, U0, b, 0)   ! below

   ! Save potential temperature:
   Scell(NSC)%Ta_var(5) = T * g_kb  ! potential temperature [K]

   Nsiz = size(Scell(NSC)%Ea_grid) ! size of the energy grid
   if (.not.allocated(Scell(NSC)%fa_eq_pot)) allocate(Scell(NSC)%fa_eq_pot(Nsiz), source = 0.0d0)
   if (T > 0.0d0) then
      do j = 1, Nsiz ! for all grid points:
         Scell(NSC)%fa_eq_pot(j) = Gibbs_x_powerDOS(Scell(NSC)%Ea_pot_grid(j), T, U0, A, b)  ! below
      enddo
   else
      Scell(NSC)%fa_eq_pot(:) = 0.0d0
      Scell(NSC)%fa_eq_pot(1) = 1.0d0
   endif

   ! For printout:
   if (numpar%save_fa) then
      Nsiz = size(Scell(NSC)%Ea_pot_grid_out) ! size of the energy grid
      if (.not.allocated(Scell(NSC)%fa_eq_pot_out)) allocate(Scell(NSC)%fa_eq_pot_out(Nsiz), source = 0.0d0)
      if (T > 0.0d0) then
         do j = 1, Nsiz
            Scell(NSC)%fa_eq_pot_out(j) = Gibbs_x_powerDOS(Scell(NSC)%Ea_pot_grid_out(j), T, U0, A, b)  ! below
         enddo
      else
         Scell(NSC)%fa_eq_pot_out(:) = 0.0d0
         Scell(NSC)%fa_eq_pot_out(1) = 1.0d0
      endif
   endif
end subroutine set_Gibbs_x_powerDOS



subroutine find_Gibbs_x_powerDOS_parameters(U1, U2, U3, A, T, U0, b, inx)
   ! Initial parameters MUST be defined:
   ! b (if inx = 0)
   ! b and U0 (if inx /= 0)
   real(8), intent(in) :: U1, U2, U3   ! first 3 moments of the distribution
   real(8), intent(inout) :: A, T, U0, b ! parameters of the distribution: A*exp(-(U-U0)/T)*U^b
   integer, intent(in) :: inx ! index, which version of this subroutine to use
   !-----------------
   real(8) :: T0, eps, alpha_b, alpha_U0, alpha_T, U0_cur, T_cur, b_cur
   integer :: Niter


   select case (inx)
   case (0) ! assume b is given, don't vary it:
      ! Define temperature from fluctuations:
      T = get_T_from_fluctuation(U1, U2, b) ! below
      ! Define shift parameter from the mean and temperature:
      U0 = get_U0_for_Gibbs_x_powerDOS(T, U1, b)  ! below

   case default   ! use b as variable:
      ! 1) Estimate the starting T and b:
      T = get_T_from_fluctuation(U1, U2, b) ! below

      ! Find corresponding b:
      b = find_b_for_Gibbs_x_powerDOS(U1, U2, U3, T, U0)   ! below
      T0 = 1.1d10    ! to start with

      ! Redefine T and U0:
      Niter = 0   ! to start with
      ! Iterate until converges:
      eps = 1.0d-5   ! precision
      ! convergence factors:
      alpha_T = 0.15d0
      alpha_U0 = 0.25d0
      alpha_b = 0.05d0
      do while ( abs(T - T0) > eps*alpha_T )
         T0 = T ! save for the next iteration
         Niter = Niter + 1 ! number of iterations

         T_cur = get_T_from_fluctuation(U1, U2, b) ! below
         T = T*(1.0d0 - alpha_T) + T_cur*alpha_T

         ! Find corresponding shift:
         U0_cur = get_U0_for_Gibbs_x_powerDOS(T, U1, b)  ! below
         U0 = U0*(1.0d0 - alpha_U0) + U0_cur*alpha_U0

         ! Find corresponding b:
         b_cur = find_b_for_Gibbs_x_powerDOS(U1, U2, U3, T, U0)   ! below
         b = b*(1.0d0 - alpha_b) + b_cur*alpha_b

         if (Niter > 10000) then
            print*, 'set_Gibbs_x_powerDOS did not converge:', Niter
            exit ! did not converge
         endif
         !write(*,'(a,i0,f,f,f,f)') 'TUb', Niter, T*g_kb, U0, b
      enddo
   end select

   ! Define A:
   A = define_A_for_GIbbs_x_powerDOS(T, b)   ! below

!    print*, 'find_Gibbs_x_powerDOS_parameters:', A, T*g_kb, U0, b
end subroutine find_Gibbs_x_powerDOS_parameters



function Gibbs_x_powerDOS(U, T, U0, A, b) result(Distrib)
   real(8) Distrib
   real(8), intent(in) :: U, T, U0, A, b
   !---------------------
   real(8) :: dU, eps
   dU = (U-U0)
   eps = 1.0d-12

   ! Gibbs distribution * power-function DOS:
   ! Note: reduces to Maxwell for b=1/2 and U0=0
   if (dU < 0.0d0) then ! no distribution at negative energies
      Distrib = 0.0d0
   else  ! there is some distribution
      Distrib = A * exp(-dU/T) * (dU**b)
   endif ! (dU < 0.0d0)
end function Gibbs_x_powerDOS


function define_A_for_GIbbs_x_powerDOS(T, b) result(A)
   real(8) A
   real(8), intent(in) :: T, b
   real(8) Gb
   Gb = GAMMA(b+1.0d0) ! intrinsic
   if ((T > 1.0d-12) .and. (abs(Gb) > 0.0d0)) then
      A = 1.0d0 / (T**(b+1.0d0) * Gb)
   else
      A = 1.0d20
   endif
end function define_A_for_GIbbs_x_powerDOS


function find_b_for_Gibbs_x_powerDOS(U1, U2, U3, T, U0) result(b)
   real(8) b   ! Gibbs x powerDOS parameter: power of the energy in DOS
   !type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   real(8), intent(in) :: U1, U2, U3   ! 1st, 2d, 3d moments of energy distribution
   real(8), intent(in) :: T, U0  ! temperature [eV], and energy shift [eV]
   !------------------
   real(8) :: LHS, b_min, b_max, eps, Fb
   integer :: i, Nat, Niter, Niter_max

   ! Define left-hand-side of the equation defining b:
   !LHS = -1.0d0/T**3 * ( 3.0d0*U0*(U2 - U1**2) - (U3 - U1**3) )
   if (T > 0.0d0) then
      LHS = 1.0d0/T**3 * (U3 - 3.0d0*U2*U0 + 3.0d0*U1*U0**2 - U0**3)
   else  ! undefined
      LHS = 1.0d10
   endif

   !print*, 'LHS:', LHS, T, U3, U2, U1, U0

   ! Find b by bisection:
   eps = 1.0d-5    ! precision

   b_min = 0.0d0   ! to start with
   b_max = 35.0d0  ! to start with

   b = (b_max + b_min)/2.0d0
   Niter = 0   ! to start with
   Niter_max = 1000
   Find_b:do while ( abs(b_max-b_min) > eps )
      Niter = Niter + 1 ! next iteration
      if (Niter > Niter_max) then
         print*, 'find_b_for_Gibbs_x_powerDOS: too many iterations:', Niter
         exit Find_b
      endif
      Fb = get_RHS_for_b(b)   ! below
      if (Fb > LHS) then
         b_max = b
      else
         b_min = b
      endif
      b = (b_max + b_min)/2.0d0
      !print*, Niter, b, Fb, LHS
   enddo Find_b

   !pause 'find_b_for_Gibbs_x_powerDOS'
end function find_b_for_Gibbs_x_powerDOS


pure function get_RHS_for_b(b) result(RHS)
   real(8) RHS
   real(8), intent(in) :: b
   !-------------
   real(8) :: Gb1, Gb2, Gb4
   Gb1 = GAMMA(b+1.0d0) ! intrinsic
   !Gb2 = GAMMA(b+2.0d0) ! intrinsic
   Gb4 = GAMMA(b+4.0d0) ! intrinsic
   !RHS = Gb4/Gb1 - (Gb2/Gb1)**3
   RHS = Gb4/Gb1
end function get_RHS_for_b



function get_T_from_fluctuation(E1, E2, b) result(T)
   real(8) T   ! [eV] temperature
   real(8), intent(in) :: E1, E2 ! [eV] Mean square, and mean, energy of atoms
   real(8), intent(in), optional :: b  ! power of DOS
   !-----------------
   integer :: Nat
   real(8) :: power_b, Fact, arg
   real(8) :: Gb1, Gb2, Gb3

   if (present(b)) then ! given power function for DOS: ~E^b
      power_b = b
   else  ! default: free-particle DOS: ~sqrt(E)
      power_b = 0.5d0
   endif

   ! Define the factor associated with DOS:
   Gb1 = GAMMA(power_b+1.0d0) ! intrinsic
   Gb2 = GAMMA(power_b+2.0d0) ! intrinsic
   Gb3 = GAMMA(power_b+3.0d0) ! intrinsic
   arg = Gb3/Gb1 - (Gb2/Gb1)**2
   Fact = 1.0d0 / sqrt( abs(arg) )

   ! Define temeprature:
   T = sqrt(E2 - E1**2) * Fact
end function get_T_from_fluctuation


pure function get_U0_for_Gibbs_x_powerDOS(T, U1, b) result(U0)
   real(8) U0  ! [eV]
   real(8), intent(in) :: T, U1, b  ! temperature [eV], mean energy [eV], power of DOS
   !--------------
   real(8) :: Gb1, Gb2
   Gb1 = GAMMA(b+1.0d0) ! intrinsic
   Gb2 = GAMMA(b+2.0d0) ! intrinsic

   U0 = U1 - T*Gb2/Gb1
end function get_U0_for_Gibbs_x_powerDOS




subroutine set_Maxwell_distribution_pot(numpar, Scell, NSC) ! Obsolete
   type(Numerics_param), intent(in) :: numpar   ! numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   !----------------------------------
   integer :: j, Nsiz
   real(8) :: arg, Tfact, E_shift, Ta, dE

   E_shift = Scell(NSC)%Pot_distr_E_shift ! shift of the distribution
   Ta = Scell(NSC)%Ta_var(5) / g_kb  ! potential temperature

   Nsiz = size(Scell(NSC)%Ea_grid) ! size of the energy grid
   if (.not.allocated(Scell(NSC)%fa_eq_pot)) allocate(Scell(NSC)%fa_eq_pot(Nsiz), source = 0.0d0)
   if (Ta > 0.0d0) then
      Tfact = (1.0d0 / Ta)**1.5d0
      do j = 1, Nsiz
         dE = Scell(NSC)%Ea_grid(j) - E_shift
         if (dE > 0.0d0) then
            arg = dE / Ta
            Scell(NSC)%fa_eq_pot(j) = 2.0d0 * sqrt( dE / g_Pi) * Tfact * exp(-arg)
         else
            Scell(NSC)%fa_eq_pot(j) = 0.0d0
         endif
         !arg = ( Scell(NSC)%Ea_grid(j) ) / Ta
         !Scell(NSC)%fa_eq_pot(j) = 2.0d0 * sqrt( (Scell(NSC)%Ea_grid(j) ) / g_Pi) * Tfact * exp(-arg)
      enddo
   else
      Scell(NSC)%fa_eq_pot(:) = 0.0d0
      Scell(NSC)%fa_eq_pot(1) = 1.0d0
   endif

   ! For printout:
   if (numpar%save_fa) then
      Nsiz = size(Scell(NSC)%Ea_pot_grid_out) ! size of the energy grid
      if (.not.allocated(Scell(NSC)%fa_eq_pot_out)) allocate(Scell(NSC)%fa_eq_pot_out(Nsiz), source = 0.0d0)
      if (Ta > 0.0d0) then
         do j = 1, Nsiz
            dE = Scell(NSC)%Ea_pot_grid_out(j) - E_shift
            if (dE > 0.0d0) then
               arg = dE / Ta
               Scell(NSC)%fa_eq_pot_out(j) = 2.0d0 * sqrt( dE / g_Pi) * Tfact * exp(-arg)
            else
               Scell(NSC)%fa_eq_pot_out(j) = 0.0d0
            endif
            !arg = ( Scell(NSC)%Ea_grid_out(j) ) / Ta
            !Scell(NSC)%fa_eq_pot_out(j) = 2.0d0 * sqrt( (Scell(NSC)%Ea_grid_out(j) ) / g_Pi) * Tfact * exp(-arg)
         enddo
      else
         Scell(NSC)%fa_eq_pot_out(:) = 0.0d0
         Scell(NSC)%fa_eq_pot_out(1) = 1.0d0
      endif
   endif
end subroutine set_Maxwell_distribution_pot




subroutine shortest_distance_to_point(Scell, i1, Sj, a_r, x1, y1, z1, sx1, sy1, sz1, cell_x, cell_y, cell_z)
   type(Super_cell), intent(in), target :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: i1 ! atom index
   real(8), dimension(3), intent(in) :: Sj   ! relative cooerds of the point, distance to which we seek
   real(8), intent(out) ::  a_r	! [A] shortest distance between the two atoms within supercell with periodic boundaries
   real(8), intent(out), optional :: x1, y1, z1		! [A] projections of the shortest distance
   real(8), intent(out), optional :: sx1, sy1, sz1 	! relative projections of the shortest distance
   integer, intent(out), optional :: cell_x, cell_y, cell_z ! cell numbers
   real(8) x, y, z, zb(3), r, x0, y0, z0, r1
   integer i, j, k, ik
   type(Atom), dimension(:), pointer :: atoms	! array of atoms in the supercell

   atoms => Scell%MDAtoms
   x = 0.0d0
   y = 0.0d0
   z = 0.0d0

   ! For the case of periodic boundaries:
   do ik = 1,3
      x = x + (atoms(i1)%S(ik) - Sj(ik))*Scell%supce(ik,1)
      y = y + (atoms(i1)%S(ik) - Sj(ik))*Scell%supce(ik,2)
      z = z + (atoms(i1)%S(ik) - Sj(ik))*Scell%supce(ik,3)
   enddo ! ik
   a_r = DSQRT(x*x + y*y + z*z)
   if (present(x1)) x1 = x
   if (present(y1)) y1 = y
   if (present(z1)) z1 = z
   if (present(sx1)) sx1 = atoms(i1)%S(1) - Sj(1)
   if (present(sy1)) sy1 = atoms(i1)%S(2) - Sj(2)
   if (present(sz1)) sz1 = atoms(i1)%S(3) - Sj(3)
   if (present(cell_x)) cell_x = 0
   if (present(cell_y)) cell_y = 0
   if (present(cell_z)) cell_z = 0

   do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
      ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell:
      ! periodic boundary conditions
      zb(1) = dble(i)
      do j =-1,1
         zb(2) = dble(j)
         do k = -1,1
            zb(3) = dble(k)
            x0 = 0.0d0
            y0 = 0.0d0
            z0 = 0.0d0
            do ik = 1,3
               x0 = x0 + (atoms(i1)%S(ik) - Sj(ik) + zb(ik))*Scell%supce(ik,1)
               y0 = y0 + (atoms(i1)%S(ik) - Sj(ik) + zb(ik))*Scell%supce(ik,2)
               z0 = z0 + (atoms(i1)%S(ik) - Sj(ik) + zb(ik))*Scell%supce(ik,3)
            enddo ! ik
            r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
            if (r1 <= a_r) then
               x = x0
               y = y0
               z = z0
               a_r = r1
               if (present(x1)) x1 = x
               if (present(y1)) y1 = y
               if (present(z1)) z1 = z
               if (present(sx1)) sx1 = atoms(i1)%S(1) - Sj(1) + zb(1)
               if (present(sy1)) sy1 = atoms(i1)%S(2) - Sj(2) + zb(2)
               if (present(sz1)) sz1 = atoms(i1)%S(3) - Sj(3) + zb(3)
               if (present(cell_x)) cell_x = i
               if (present(cell_y)) cell_y = j
               if (present(cell_z)) cell_z = k
            endif
         enddo ! k
      enddo ! j
   enddo ! i
end subroutine shortest_distance_to_point


END MODULE Atomic_thermodynamics
