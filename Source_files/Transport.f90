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
! This module contains transport equations for electron particle and energy,
! and for heat transport in the atomic system out of the simulation box,
! using Berendsen thermostat:
! https://en.wikipedia.org/wiki/Berendsen_thermostat

MODULE Transport
use Universal_constants
use Electron_tools
use Atomic_tools
use TB

implicit none
 
 contains


subroutine Electron_transport(trans_mod, time, Scell, numpar, matter, dt, tau, Err)
   integer, intent(in) :: trans_mod ! model for the transport to use
   real(8), intent(in) :: time ! current timestep [fs]
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter ! material parameters
   real(8), intent(in) :: dt, tau ! [fs]: time-step, characteristic time of relaxation
   type(Error_handling), intent(inout) :: Err ! error save
   !-------------------------------
   integer NSC
   select case (trans_mod) ! in case we later add different models for the electron transport
   case (0) ! excluded transport
      print*, 'No electron transport out of the simulation box is included'
   case default ! simple rate equation (Berendsen thermostat)
      do NSC = 1, size(Scell) ! for all supercells
         !print*, 'Electron_transport', matter%T_bath_e, dt, tau
         call rate_equation(Scell(NSC)%TeeV, matter%T_bath_e, dt, tau, 2) ! below
         select case (numpar%el_ion_scheme)
         case (4)    ! relaxation-time approximation
            call rate_equation_for_fe(Scell(NSC), matter%T_bath_e, dt, tau)   ! below
         case default
            call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution, module "Electron_tools"
         endselect
         call get_new_energies(Scell, matter, numpar, time, Err) ! module "TB"
      enddo
   endselect
end subroutine Electron_transport



subroutine rate_equation_for_fe(Scell, T_bath_e, dt, tau)   ! below
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   real(8), intent(in) :: T_bath_e  ! [eV] electronic temperature of the bath
   real(8), intent(in) :: dt, tau   ! MD timestep [fs], and characteristic cooling time [fs]
   !-----------------------------
   real(8), dimension(size(Scell%Ei)) :: fe_eq  ! distribution function of electronic bath
   real(8) :: Ntot, mu, exp_dttau
   integer :: i_fe, i

   ! 1) Find the chemical potential of the bath (corresponding to the constant number of electrons):
   !Ntot = dble(Scell%Ne)   ! number of electrons to be conserved
   Ntot = Scell%Ne_low  ! transient number of low-energy electrons (excluding high-energy ones!)
   call find_mu_from_N_T(Scell%Ei, Ntot, mu, T_bath_e) ! module "Electron_tools"

   ! 2) Set the disitrbution function of the electornic bath:
   call set_Fermi(Scell%Ei, T_bath_e, mu, fe_eq)   ! module "Electron_tools"
   !print*, 'rate_equation_for_fe 1:', Ntot, SUM(fe_eq(:))

   ! 3) Do the rate equation for the electron distribution:
   i_fe = size(Scell%fe)   ! number of grid points in distribution function
   exp_dttau = dexp(-dt / tau)   ! finite-time relaxation
   do i = 1, i_fe ! for all grid points (MO energy levels)
      Scell%fe(i) = fe_eq(i) + (Scell%fe(i) - fe_eq(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
   enddo

   ! 4) Update the equivalent distribution:
   ! Update the electronic band energy:
   call set_total_el_energy(Scell%Ei, Scell%fe, Scell%nrg%El_low) ! module "Electron_tools"
   ! Get the equivalent (kinetic) temperature and chemical potential:
   call Electron_Fixed_Etot(Scell%Ei, Scell%Ne_low, Scell%nrg%El_low, Scell%mu, Scell%TeeV, .true.) ! module "Electron_tools"
   Scell%Te = Scell%TeeV*g_kb ! save also in [K]
   ! Construct Fermi function with the given transient parameters:
   if (.not.allocated(Scell%fe_eq)) allocate(Scell%fe_eq(i_fe))
   call set_Fermi(Scell%Ei, Scell%TeeV, Scell%mu, Scell%fe_eq)   ! below

   !print*, 'rate_equation_for_fe 2:', Ntot, Scell%Ne_low, SUM(fe_eq(:)), SUM(Scell%fe(:))
end subroutine rate_equation_for_fe



subroutine rate_equation(X, X0, dt, tau, ind) ! solving dX/dt = -(X - X0)/tau
   real(8), intent(inout) :: X ! the value that changes according to the rate equation
   real(8), intent(in) :: X0   ! the equilibrium value that X approaches
   real(8), intent(in) :: dt   ! the time-step
   real(8), intent(in) :: tau  ! the characteristic time of the rate equation
   integer, intent(in) :: ind  ! which scheme to use
   real(8) dttau
   dttau = dt/tau
   select case (ind)
   case (1) ! implicit
      X = (X + dttau*X0)/(1.0d0 + dttau)
   case (2) ! analytical
      X = X0 + (X - X0)*dexp(-dttau)
   case default ! explicit
      X = X - dttau*(X - X0)
   end select
end subroutine rate_equation



subroutine Atomic_heat_transport(trans_mod, Scell, matter, dt, tau)
   integer, intent(in) :: trans_mod ! model for the transport to use
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter ! material parameters
   real(8), intent(in) :: dt, tau ! [fs]: time-step, characteristic time of relaxation
   !-------------------------------
   real(8) Ekin, dE
   integer NSC
   select case (trans_mod) ! in case we later add different models for the electron transport
   case (0) ! exclude heat transport within atomic subsystem
      print*, 'No heat transport in atomic system out of simulation box is included'
   case default ! simple rate equation (Berendsen thermostat)
      do NSC = 1, size(Scell) ! for all supercells
         call rate_equation(Scell(NSC)%TaeV, matter%T_bath, dt, tau, 2) ! above
         Scell(NSC)%Ta = Scell(NSC)%TaeV*g_kb	! [K]
         call get_energy_from_temperature(dble(Scell(NSC)%Na), Scell(NSC)%TaeV, Ekin) ! get new total kinetic energy, module "Atomic_tools"
         dE = Ekin - SUM(Scell(NSC)%MDatoms(:)%Ekin) ! change of the total kinetic energy due to change of temperature [eV]
         ! change velocities according to energy gain/loss:
         call Rescale_atomic_velocities(dE, matter, Scell, NSC, Scell(NSC)%nrg) ! module "Atomic tools"
         call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg) ! update energy value, module "Atomic_tools"
         call save_last_timestep(Scell) ! Update the last time-step data accordingly, module "Atomic_tools"
      enddo
   endselect
end subroutine Atomic_heat_transport



subroutine Change_affected_layer(trans_mod, dd, Scell, dt, tau) ! change affected layer for optical consts calculation
   real(8), intent(inout) :: dd ! Scell(NSC)%eps%dd : affected layer thickness [A]
   integer, intent(in) :: trans_mod ! model for the transport to use
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   real(8), intent(in) :: dt, tau ! [fs]: time-step, characteristic time of relaxation
   real(8) Ekin, dE
   integer NSC
   select case (trans_mod) ! in case we later add different models for the electron transport
   case (0) ! exclude heat transport within atomic subsystem
      print*, 'No change of the affected layer is included'
   case default ! simple rate equation
      do NSC = 1, size(Scell) ! for all supercells
         ! Change the layer thickness to the given FIXED value (used for testing!!!)
         call rate_equation(dd, 0.0d0, dt, tau, 2) ! above
      enddo
   endselect
end subroutine Change_affected_layer


END MODULE Transport
