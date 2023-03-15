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
MODULE Electron_tools
use Universal_constants
use Objects
!use Variables
use Algebra_tools
use Atomic_tools
use Little_subroutines

implicit none
 contains


subroutine update_cross_section(Scell, matter)
   type(Super_cell), intent(in) :: Scell ! supercell with all the atoms as one object
   type(solid), intent(inout) :: matter	! materil parameters
   integer :: Nshl, i, N_Te
   real(8) :: dT, T_left
   ! Get the mean free paths vs Te:
   Nshl = size(matter%Atoms(1)%Ip)
   select case (matter%Atoms(1)%TOCS(Nshl)) ! Valence band and CDF only
   case (1) ! CDF
      if (Scell%Te > 100.0d0) then  ! recalculate:
         ! Temperature (grid defined in subroutine get_MFPs as Te_temp = dble((i-1)*1000)):
         dT = 1000.0d0 ! [K] grid step
         T_left = FLOOR(Scell%Te/1000)*1000
         N_Te = CEILING(Scell%Te/1000)
         if (N_Te > size(matter%Atoms(1)%El_MFP_vs_T)) N_Te = size(matter%Atoms(1)%El_MFP_vs_T)   ! maximal energy set
         ! Interpolate valence band MFP for the given temperature:
         if (N_Te == 1) then
            matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) + &
          (matter%Atoms(1)%El_MFP_vs_T(N_Te+1)%L(:) - matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:))/dT * (Scell%Te-T_left)
         else
            matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) + &
          (matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) - matter%Atoms(1)%El_MFP_vs_T(N_Te-1)%L(:))/dT * (Scell%Te-T_left)
         endif
      else ! no need to recalculate, the tempereature is too small:
         N_Te = 1
         matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:)
      endif
   endselect

!    print*, 'Te=', Scell%Te, matter%Atoms(1)%El_MFP(Nshl)%L(1), matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(1)
end subroutine update_cross_section


subroutine find_band_gap(wr, Scell, matter, numpar)
   REAL(8), DIMENSION(:), INTENT(in) ::  wr	! [eV] energy levels
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(solid), intent(inout) :: matter	! materil parameters
   type(Numerics_param), intent(inout) :: numpar ! numerical parameters, including MC energy cut-off
   integer i, sumNe, siz, j

   siz = size(wr)
   
   if (Scell%N_Egap <= 0) then
      sumNe = 0
      !i = 1
      do i = 1, siz
      !do while (sumNe <= Scell%Ne)
         sumNe = sumNe + 2
         if (sumNe >= Scell%Ne) exit ! HOMO level
      enddo
      Scell%N_Egap = i ! save, this is the energy level where bandgap starts
   else
      i = Scell%N_Egap
   endif

   Scell%E_gap = ABS(wr(i+1) - wr(i))	! bandgap [eV]
   Scell%E_bottom = wr(i+1)	! bottom of the CB [eV]
   ! If top energy level is excluded (in parameterization, or just too high due to convergence issues)
   j = siz  ! start from the last one
   Scell%E_top = wr(j)		! [eV] current top of the (meaningful) conduction band
   do while (wr(j)>=79.0d0)
      j = j - 1
      Scell%E_top = wr(j)		! [eV] current top of the (meaningful) conduction band
   enddo
   Scell%E_VB_bottom = wr(1)	! [eV] current bottom of the valence band
   Scell%E_VB_top = wr(i)		! [eV] current top of the valence band
   
   ! Set MC high-energy electron cut-off energy equal to the uppermost level of CB:
   !if (numpar%E_cut_dynamic) numpar%E_cut = wr(siz) - Scell%E_bottom ! [eV]
   if (numpar%E_cut_dynamic) numpar%E_cut = Scell%E_top - Scell%E_bottom ! [eV]
   ! For noneuqilibrium distributions (BO or relaxation time), threshold cannot be higher than
   ! the topmost level of CB, otherwise, there is no way to place an incomming electron:
   select case (numpar%el_ion_scheme)
   case (3:4)
      numpar%E_cut = min(numpar%E_cut, Scell%E_top-Scell%E_bottom) ! [eV]
   endselect

   ! In case we have a very wide gap material, cut-off cannot be smaller than the gap:
   if (numpar%E_cut < Scell%E_gap) numpar%E_cut = Scell%E_gap

   select case (matter%Atoms(1)%TOCS(size(matter%Atoms(1)%TOCS))) ! which inelastic cross section to use (BEB vs CDF):
      case (1) ! CDF cross section
         matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)) = Scell%E_gap ! [eV] ionization potential of the valence band
      case default  ! BEB:
         ! Renormalization is optional:
         !if (numpar%E_cut < matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)) ) numpar%E_cut = matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip))
   end select
   
!    print*, 'numpar%E_cut =', numpar%E_cut , Scell%E_gap , matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip))
!    print*, 'Ne:', Scell%Ne_low, Scell%Ne, Scell%Na
!    print*, 'E :', Scell%E_gap, Scell%E_VB_top, Scell%E_bottom
end subroutine find_band_gap


subroutine get_number_of_CB_electrons(Scell, NSC)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   Scell(NSC)%Ne_CB = SUM(Scell(NSC)%fe(:), MASK = (Scell(NSC)%Ei >= Scell(NSC)%E_bottom))
   Scell(NSC)%Ne_CB = Scell(NSC)%Ne_CB + Scell(NSC)%Ne_high ! total number of conduction band electrons
end subroutine get_number_of_CB_electrons



subroutine update_fe(Scell, matter, numpar, t, Err, do_E_tot)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(solid), intent(inout) :: matter ! materil parameters
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   real(8), intent(in) :: t ! current timestep [fs]
   type(Error_handling), intent(inout) :: Err ! error save
   logical, intent(in), optional :: do_E_tot  ! total energy is given or temperature?
   !=========================================
   real(8) :: E_tot, mu_cur, Te_cur
   integer :: NSC, i_fe
   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
      do NSC = 1, size(Scell)
         ! Which scheme to use:
         ! 0=decoupled electrons; 1=enforced energy conservation; 2=T=const; 3=BO
         select case (numpar%el_ion_scheme)
         case (1) ! Enforced energy conservation (Etot = Ee + Eat = const):
            if (t .GT. numpar%t_Te_Ee) then ! Total energy is fixed:
               !call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! (SLOW) below
               call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! (FAST) below
            else ! electron temperature is fixed:
               call Electron_Fixed_Te(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! below
            endif
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]
            call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution
            Scell(NSC)%fe_eq = Scell(NSC)%fe ! instanteneous thermalization means both functions are the same

         case (2) ! Fixed temperature (Te=const):
!             if (numpar%scc) then ! SCC, so the total energy is defined by the part H_0 without charge energy:
!                call Electron_Fixed_Te(Scell(NSC)%Ei_scc_part, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! below
!             else
               call Electron_Fixed_Te(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! below
!             endif
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]
            call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution
            Scell(NSC)%fe_eq = Scell(NSC)%fe ! instanteneous thermalization means both functions are the same

         case (3) ! Born-Oppenheimer:
            ! Do nothing with fe!
            ! Only get the kinetic temperature of electrons (out-of-equilibrium):
            call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, &
                                          Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! below (FAST)
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]

            ! Construct Fermi function with the given transient parameters (equivalent Te and mu):
            i_fe = size(Scell(NSC)%fe)   ! number of grid points in distribution function
            if (.not.allocated(Scell(NSC)%fe_eq)) allocate(Scell(NSC)%fe_eq(i_fe))
            call set_Fermi(Scell(NSC)%Ei, Scell(NSC)%TeeV, Scell(NSC)%mu, Scell(NSC)%fe_eq)   ! below

         case (4) ! Relaxation-time approximation [ df/dt=(f-f0)/tau ]:
            ! Relaxing electrons via rate equation with given characteristic time:
            !call Do_relaxation_time(Scell(NSC), numpar)  ! below
            ! We only update it once per simulation step, not every time this subroutine called!

         case (50) ! Boltzmann electron-electron collision integral (NOT READY, DO NOT USE!):
            if (t > -8.5d0) then ! testing, unfnished
               call test_evolution_of_fe(Scell(NSC)%Ei, Scell(NSC)%fe, t) ! see below
            endif

            ! Only get the kinetic temperature of electrons (out-of-equilibrium):
            call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, &
                                          Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! below (FAST)
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]

            ! Construct Fermi function with the given transient parameters (equivalent Te and mu):
            i_fe = size(Scell(NSC)%fe)   ! number of grid points in distribution function
            if (.not.allocated(Scell(NSC)%fe_eq)) allocate(Scell(NSC)%fe_eq(i_fe))
            call set_Fermi(Scell(NSC)%Ei, Scell(NSC)%TeeV, Scell(NSC)%mu, Scell(NSC)%fe_eq)   ! below


         case default ! Decoupled electrons and ions (Ee = const; instant thermalization of electrons):
            !call set_total_el_energy(Scell(NSC)%Ei, Scell(NSC)%fe, Scell(NSC)%nrg%E_tot) ! get the total electron energy
!             if (numpar%scc) then ! SCC, so the total energy is defined by the part H_0 without charge energy:
!                ! get the total electron energy:
!                call set_total_el_energy(Scell(NSC)%Ei_scc_part, Scell(NSC)%fe, Scell(NSC)%nrg%El_low)
!                ! thermalize the distribution function over the SCC-band levels:
!                call Electron_Fixed_Etot(Scell(NSC)%Ei_scc_part, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, &
!                                           Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! below (FAST)
!             else  ! non-SCC, reause the precalculated energy:

               if (.not.present(do_E_tot)) then ! if we do not have total energy given, get it from distribution:
                  call set_total_el_energy(Scell(NSC)%Ei, Scell(NSC)%fe, Scell(NSC)%nrg%El_low) ! get the total electron energy
               endif
               call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, &
                                          Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! below (FAST)
!             endif
!

            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]
            call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution
            if (.not.allocated(Scell(NSC)%fe_eq)) allocate(Scell(NSC)%fe_eq(i_fe))
            Scell(NSC)%fe_eq = Scell(NSC)%fe ! instanteneous thermalization means both functions are the same

         end select

         ! Update the number of CB electrons:
         call get_number_of_CB_electrons (Scell, NSC)
      enddo
   endif DO_TB
end subroutine update_fe


!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! Nonequilibrium electron kinetics

subroutine Electron_thermalization(Scell, numpar, skip_thermalization)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   logical, intent(in), optional :: skip_thermalization
   !----------------------
   select case (numpar%el_ion_scheme)
   case (4)    ! relaxation-time approximation
      if (present(skip_thermalization)) then
         call Do_relaxation_time(Scell(1), numpar, skip_thermalization)  ! below
      else
         call Do_relaxation_time(Scell(1), numpar)  ! below
      endif
   endselect
end subroutine Electron_thermalization



! Relaxation time approximation:
subroutine Do_relaxation_time(Scell, numpar, skip_thermalization)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   logical, intent(in), optional :: skip_thermalization
   !----------------------
   real(8) :: exp_dttau, extra_dt, extra_tau
   integer :: i_fe, i, i_cycle, N_cycle
   logical :: skip_step, extra_cycle

   if (present(skip_thermalization)) then
      skip_step = skip_thermalization
   else
      skip_step = .false.
   endif

   ! Get the equivalent (kinetic) temperature and chemical potential:
   call Electron_Fixed_Etot(Scell%Ei, Scell%Ne_low, Scell%nrg%El_low, Scell%mu, Scell%TeeV, .true.) ! below (FAST)
   Scell%Te = Scell%TeeV*g_kb ! save also in [K]

   ! Construct Fermi function with the given transient parameters:
   i_fe = size(Scell%fe)   ! number of grid points in distribution function
   if (.not.allocated(Scell%fe_eq)) allocate(Scell%fe_eq(i_fe))
   call set_Fermi(Scell%Ei, Scell%TeeV, Scell%mu, Scell%fe_eq)   ! below

   ! Solve rate equation:
   if (.not.skip_step) then ! do the thermalization step:
      if (numpar%tau_fe < numpar%dt/30.0d0) then ! it's basically instantaneous
         exp_dttau = 0.0d0
      else  ! finite time relaxation
         exp_dttau = dexp(-numpar%dt / numpar%tau_fe)
      endif
      do i = 1, i_fe ! for all grid points (MO energy levels)
         Scell%fe(i) = Scell%fe_eq(i) + (Scell%fe(i) - Scell%fe_eq(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
      enddo

      !--------------------------
      ! Extra check for smoothening unphysical artefacts that may be present after MC:
      extra_cycle = .false.   ! by default, assume no artifact
      do i = 1, i_fe ! for all grid points (MO energy levels)
         if ((Scell%fe(i) > 2.0d0) .or. (Scell%fe(i) < 0.0d0)) then
            extra_cycle = .true.   ! some artefacts present, do extra thermalization to get rid of them
            print*, 'Extra thermalization step needed:', i, Scell%fe(i)
            exit
         endif
      enddo
      if (extra_cycle) then   ! do extra thermalization
         extra_dt = numpar%dt*0.1d0 ! use this small step to minimize the effect of extra smoothing
         extra_tau = min(numpar%tau_fe, 10.0d0) ! artificial thermalization steps
         N_cycle = 10000 ! limit for the cycles
         i_cycle = 0 ! to start
         do while (extra_cycle)
            i_cycle = i_cycle + 1
            extra_cycle = .false.   ! assume the problem is solved
            if (numpar%tau_fe < extra_dt/30.0d0) then ! it's basically instantaneous
               exp_dttau = 0.0d0
            else  ! finite time relaxation
               exp_dttau = dexp(-extra_dt / extra_tau)
            endif
            do i = 1, i_fe ! for all grid points (MO energy levels)
               Scell%fe(i) = Scell%fe_eq(i) + (Scell%fe(i) - Scell%fe_eq(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
               if ((Scell%fe(i) > 2.0d0) .or. (Scell%fe(i) < 0.0d0)) then
                  !print*, 'Step:', i, Scell%fe(i)
                  if (i_cycle < N_cycle) extra_cycle = .true.   ! artefacts still present, do another cycle of extra thermalization
               endif
            enddo ! i = 1, i_fe
         enddo ! while (extra_cycle)
      endif ! (extra_cycle)
   endif ! (.not.skip_step)
end subroutine Do_relaxation_time



! Electron-electron collision integral (UNFINISHED DUE TO PROBLEMS WITH ENERGY CONSERVATION):
subroutine Boltzmann_e_e_IN(Ev, fe, M_ee, dt) ! calculates change of distribution function via Boltzmann collision integral
! See examples of Boltzmann equation in energy space e.g. in 
! [B. Rethfeld, A. Kaiser, M. Vicanek and G. Simon, Phys.Rev.B 65, 214303, (2002)]
! (although we use here only electron-electron integral, and with very different matrix element)
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   real(8), dimension(:), intent(inout) :: fe  ! electron distribution function (ocupation numbers of the energy levels Ei)
   real(8), dimension(:,:), intent(in) :: M_ee ! electron-electron scattering matrix elements
   real(8), intent(in) :: dt ! [fs] time-step
   !-----------------------------------------
   real(8) :: Ei, Ef, Ei2, Ef2 ! energies of initial and final states for #1 and #2 electron
   real(8), dimension(size(fe)) :: fe_temp ! temporary distribution to work with
   real(8), dimension(size(fe)) :: two_minus_fe ! temporary distribution to work with
   real(8), dimension(size(fe),size(fe)) :: fe_2_fe! array of multiplied fe(i)*(2-fe(j))

   real(8) :: fe_f ! average final state
   real(8) :: E_tot, E_tot2 ! total energy used for testing [eV]
   real(8) :: prefactor, N_tot, N_tot2
   integer :: i, i2, j2, N, i_cur, i_f1, i_f2
   logical :: within_range ! check if final energy level is within possible range

   N = size(fe) ! total number of energy levels available

   prefactor =  2.0d0*g_Pi/g_h*g_e*1.0d-15*dt ! prefactor of 2Pi/h_bar*dt -> in eV and fs
   
   fe_temp = fe ! just to start
   two_minus_fe = 2.0d0 - fe ! number of free places in the final states

   ! Test the energy conservation:
   call set_total_el_energy(Ev,fe_temp,E_tot)
!    print*, 'FIRST :', E_tot/64.0d0, SUM(fe_temp)
   
   ! Construct an array that will be reused multiple times: fe(i)*(2-fe(j)):
   call Two_Vect_Matr(fe,two_minus_fe,fe_2_fe) ! module "Algebra_tools"
   
!$omp PARALLEL private(i)
!$omp do schedule(dynamic) reduction( + : fe_temp)
   do i = 1, N ! all energy levels on which we looking for changes (left-hand side of the equation)
      call Boltzmann_solution(i, Ev, M_ee, fe_temp, fe, two_minus_fe, fe_2_fe, prefactor, -1) ! see below
   enddo
!$omp end do
!$omp end parallel

   ! Test the energy conservation:
   call set_total_el_energy(Ev,fe_temp,E_tot2)
!    print*, 'SECOND:', E_tot2/64.0d0, SUM(fe_temp), ABS(E_tot2 - E_tot)/ABS(E_tot)*100.0d0
   N_tot = SUM(fe)
   N_tot2 = SUM(fe_temp)
   print*, 'Boltzmann_e_e N1:', N_tot, N_tot2, ABS(N_tot2-N_tot)/N_tot*100.0d0
   print*, 'Boltzmann_e_e E1:', E_tot, E_tot2, ABS(E_tot2 - E_tot)/ABS(E_tot)*100.0d0
 
   
   fe = fe_temp ! output: updated distribution function
   
!     PAUSE "Boltzmann_e_e PUASE"
!     ! Test:
!     do i = 1,N
!        print*, 'fe=', i, fe(i)
!     enddo
   
   
end subroutine Boltzmann_e_e_IN


subroutine Boltzmann_solution(i, Ev, M_ee, fe_temp, fe, two_minus_fe, fe_2_fe, prefactor, scheme)
   integer, intent(in) :: i ! current energy level
   real(8), dimension(:), intent(in), target :: Ev     ! [eV] electron energy levels
   real(8), dimension(:), intent(inout) :: fe_temp  ! electron distribution function on the curent time-step
   real(8), dimension(:), intent(in) :: fe  ! electron distribution function on the last time-step
   real(8), dimension(:), intent(in) :: two_minus_fe ! (2-fe)
   real(8), dimension(:,:), intent(in) :: M_ee    ! electron-electron scattering matrix elements
   real(8), dimension(:,:), intent(in) :: fe_2_fe ! array of multiplied fe(i)*(2-fe(j))
   real(8), intent(in) :: prefactor ! multiplier in the Boltzmann collision integral
   integer, intent(in) :: scheme ! which integration scheme in the Boltzmann equation to use: 0=explicit, 1=implicit
   !--------------------------
   integer :: i2, j2, N, i_f1, i_f2
   real(8), pointer :: Ei, Ei2, Ef2
   real(8) :: Ef, fe_f, F1, F2, F_temp, a, F2F1
   logical :: within_range ! check if final energy level is within possible range
   
   N = size(fe)
   Ei => Ev(i) ! [eV]
   F1 = 0.0d0
   F2 = 0.0d0
   
   E_FROM:do i2 = 1, N ! all levels from where electron can scatter off
      Ei2 => Ev(i2) ! [eV]
         
      E_TO:do j2 = 1, N ! all levels where electron can come to
         if (j2 /= i2) then ! transition must be between 2 different levels
            Ef2 => Ev(j2) ! [eV]
            Ef = Ei + (Ei2 - Ef2) ! [eV]
            within_range = .true. ! model it by default, check later
            ! Exclude transitions outside of energy levels given within TB:
            if (Ef < Ev(1) .or. Ef > Ev(N)) then 
               within_range = .false.
            endif

            if (within_range) then ! if transition is possible
               ! Find the value of distribution in-between the energy levels:
               call mean_distribution(Ev, fe, Ef, fe_f=fe_f, i_f1=i_f1, method=0) ! see below

               select case (scheme) ! scheme of integration
               case (0) ! explicit
!                 F1 = F1 + M_ee(i,i_f1)*( fe_f*two_minus_fe(i)*fe(j2)*two_minus_fe(i2) - fe(i)*(2.0d0 - fe_f)*fe(i2)*two_minus_fe(j2) )
                  F1 = F1 + M_ee(i,i_f1)*( fe_f*two_minus_fe(i)*fe_2_fe(j2,i2) - fe(i)*(2.0d0 - fe_f)*fe_2_fe(i2,j2) )
               case default ! implicit: much more stable!
!                   F_temp = M_ee(i,i_f1)*( fe_f*fe(j2)*two_minus_fe(i2) )
                  F_temp = M_ee(i,i_f1)*fe_f*fe_2_fe(j2,i2)
                  F1 = F1 + F_temp
!                   F2 = F2 + F_temp + M_ee(i,i_f1)*( fe(i2)*two_minus_fe(j2)*(2.0d0 - fe_f) )
                  F2 = F2 + F_temp + M_ee(i,i_f1)*( fe_2_fe(i2,j2)*(2.0d0 - fe_f) )
               end select
            endif ! within_range
         endif ! (j2 /= i2)
      enddo E_TO
   enddo E_FROM
   
   select case (scheme) ! scheme of integration
   case (0) ! explicit
      fe_temp(i) = fe_temp(i) + prefactor*F1
   case (1) ! implicit: much more stable (unconditionally stable?)
      fe_temp(i) = (fe_temp(i) + 2.0d0*prefactor*F1)/(1.0d0 + prefactor*F2)
   case default ! exact solution of dx/dt = -F2*x + F1
      F2F1 = 2.0d0*F1/F2
      fe_temp(i) = (fe_temp(i) - F2F1)*exp(-F2*prefactor) + F2F1
   end select
   
   nullify(Ei, Ei2, Ef2) ! free pointers
end subroutine Boltzmann_solution


subroutine mean_distribution(Ev, fe, Ef, fe_f, i_f1, method)
   real(8), dimension(:), intent(in) :: Ev ! [eV] electron energy levels
   real(8), dimension(:), intent(in) :: fe ! electron distribution function
   real(8), intent(in) :: Ef    ! given energy (in between energy levels)
   real(8), intent(out) :: fe_f  ! distribution function at the given point between the energy levels
   integer, intent(out), optional :: i_f1  ! number of energy level the closest one
   integer, intent(in), optional :: method ! what kind of average do we use
   real(8) Fm, sigma, mu, Gaus, three_sigma
   integer i_f, i_cur, N, i_start, i_end
   
   N = size(eV) ! number of energy levels
   sigma = 0.1d0 ! [eV] width of each energy level assigned
   three_sigma = 3.0d0*sigma
   
   call Find_in_array_monoton(Ev, Ef, i_f) ! module "Little_subroutines"
   if (present(i_f1)) i_f1 = i_f ! for the output

   select case (method)
   case (1) ! try Fermi interpolation :: POOR CONSERVATION PROPERTIES
      call Fermi_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
   case (2) ! try "geometric" average :: DOESN'T WORK
      Fm = fe(i_f-1)*(Ef - Ev(i_f-1)) + fe(i_f)*(Ev(i_f) - Ef)
      fe_f = fe(i_f-1)*fe(i_f)*(Ev(i_f) - Ev(i_f-1))/Fm
   case (3) ! do gaussian width of energy levels :: DOESN'T WORK
      fe_f = 0.0d0
      call Find_in_array_monoton(Ev, Ef-three_sigma, i_start) ! module "Little_subroutines"
      call Find_in_array_monoton(Ev, Ef+three_sigma, i_end) ! module "Little_subroutines"
      do i_cur = i_start, i_end ! within three sigma range
         call Gaussian(Ev(i_cur), sigma, Ef, Gaus)
         fe_f = fe_f + Gaus !*(2.0d0*sqrt(2.0d0*g_Pi)*sigma) ! we use this kind of nuormalization to make 2 electrons exctly at a level
         !print*, 'Gaus', i_cur, fe_f, Gaus, (2.0d0*sqrt(2.0d0*g_Pi)*sigma)
      enddo
   case (4) ! try linear interpolation + gaussian width :: DOESN'T WORK
      if (i_f > 1) then
         call linear_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
         if ((Ef - Ev(i_f-1)) > (Ev(i_f) - Ef)) then
            call Gaussian(Ev(i_f), sigma, Ef, Gaus, 2.0d0)
         else 
            call Gaussian(Ev(i_f-1), sigma, Ef, Gaus, 2.0d0)
         endif
         fe_f = fe_f*Gaus
      else
         fe_f = fe(1)
      endif
   case (5) ! linear interpolation with exclusion of too-far-lying levels
      if (i_f > 1) then
         if ((Ef - Ev(i_f-1)) < three_sigma .or. (Ev(i_f) - Ef) < three_sigma) then
            call linear_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
         endif
      else
         fe_f = fe(1)
      endif
   case default ! linear interpolation
      if (i_f > 1) then
         call linear_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
      else
         fe_f = fe(1)
      endif
   endselect

end subroutine mean_distribution


subroutine find_average_level(fe, Ev, i_f1, E_mean, fe_mean)
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   real(8), dimension(:), intent(in) :: fe  ! electron distribution function (ocupation numbers of the energy levels Ei)
   integer, intent(in) :: i_f1 ! number of the level
   real(8), intent(in) :: E_mean ! [eV] exact energy we want to find corresponding occupation to
   real(8), intent(out) :: fe_mean ! average occupation number
   integer :: i_plus
   if (i_f1 == size(Ev)) then ! borderline:
      fe_mean = fe(i_f1)
   else ! within the borders:
      i_plus = i_f1+1
      fe_mean = (fe(i_f1)*(Ev(i_plus) - E_mean) + fe(i_plus)*(E_mean - Ev(i_f1))) / (Ev(i_plus) - Ev(i_f1))
   endif
end subroutine find_average_level


subroutine test_evolution_of_fe(Ev, fe, t)
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   real(8), dimension(:), intent(inout) :: fe  ! electron distribution function (ocupation numbers of the energy levels Ei)
   real(8), intent(in) :: t ! [fs] current time-step
   !------------------------------
   real(8), dimension(size(fe)) :: fe_test  ! electron distribution function (ocupation numbers of the energy levels Ei)
   real(8), dimension(size(fe),size(fe)) :: M_ee ! electron-electron scattering matrix elements
   real(8) :: dt ! [fs] time-step
   real(8) :: Se ! electronic entropy
   integer i

   if (t < -8.4) call test_change_of_fe(fe, fe_test) ! just to test!!

   M_ee = 0.1d0  ! just to test!!
   dt = 0.001d0  ! [fs]

   do i = 1, 100
!       call Boltzmann_e_e(Ev, fe_test, M_ee, dt) ! calculates change of distribution function via Boltzmann collision integral
      call Boltzmann_e_e_IN(Ev, fe_test, M_ee, dt) ! calculates change of distribution function via Boltzmann collision integral
      print*, 'Boltzmann_e_e is done', i
      call electronic_entropy(fe_test, Se) ! subroutine above
      print*, 'Entropy:', i, Se
   enddo

end subroutine test_evolution_of_fe


subroutine test_change_of_fe(fe, fe_test)
   real(8), dimension(:), intent(in) :: fe
   real(8), dimension(:), intent(out) :: fe_test
   real(8) :: RN
   integer :: i, N
   N = size(fe)
   do i = 1, N
      call random_number(RN)
      fe_test(i) = fe(i) + 0.001d0*(RN-0.5d0)*fe(i)
   enddo

   where (fe_test(:) > 2.0d0) ! cannot be
      fe_test(:) = 2.0d0
   end where
   where (fe_test(:) < 0.0d0) ! cannot be
      fe_test(:) = 0.0d0
   end where
end subroutine test_change_of_fe


subroutine share_energy(Ev, i, j, Eprime, a, b) 
! Shares energy between levels i and j so that total energy Eprime is conserved
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   integer, intent(in) :: i, j ! levels to which we give parts of the energy Eprime
   real(8), intent(in) :: Eprime ! [eV] energy that is somewhere in between levels Ev(i) and Ev(j)
   real(8), intent(out) :: a, b  ! coefficients, how much energy must be in the level i and j to conserve total energy Eprime
   real(8) :: Ediff
   Ediff = (Ev(i)-Ev(j))
   a = (Eprime - Ev(j))/Ediff
   b = (Ev(i) - Eprime)/Ediff
!    print*, 'a,b', a, b, Eprime, Ev(i), Ev(j)
!    pause 
end subroutine share_energy



!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! General subroutines that are often used:

subroutine electronic_entropy(fe, Se, norm_fe)
   real(8), dimension(:), intent(in) :: fe ! electron distribution function
   real(8), intent(out) :: Se ! self-explanatory
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   ! Se = -kB * int [ DOS*( f * ln(f) + (1-f) * ln(1-f) ) ]
   ! E.G. [https://doi.org/10.1103/PhysRevB.50.14686]
   !----------------------------
   real(8), dimension(size(fe)) :: f_lnf
   real(8) :: f_norm, eps
   integer :: i

   eps = 1.0d-12  ! precision

   if (present(norm_fe)) then   ! user provided
      f_norm = norm_fe
   else ! by default, not spin resolved
      f_norm = 2.0d0
   endif
   Se = 0.0d0
   f_lnf = 0.0d0
   ! First term of the total entropy:
   where (fe(:) > eps) f_lnf(:) = fe(:)*log(fe(:)/f_norm) ! our f is normalized to f_norm, which means it includes DOS in it, so divide by f_norm where needed
!    do i = 1, size(fe)
!       if (fe(i) > 0.0d0) then
!          if (fe(i)/f_norm < 1.0d-10) print*, i, fe(i), f_norm
!          f_lnf(i) = fe(i)*log(fe(i)/f_norm)
!          if (fe(i)/f_norm < 1.0d-10) print*, i, 'done'
!       endif
!    enddo
   Se = SUM(f_lnf(:))
   f_lnf = 0.0d0
   ! Second term of the total entropy:
   where (fe(:) < f_norm-eps) f_lnf(:) = (f_norm - fe(:))*log((f_norm - fe(:))/f_norm)
   Se = Se + SUM(f_lnf(:))
   !Se = -g_kb*Se
   Se = -g_kb_EV*Se  ! [eV/K]
end subroutine  electronic_entropy



subroutine get_DOS_sort(Ei, DOS, smearing, partial_DOS, masks_DOS, Hij, CHij)
   real(8), dimension(:), intent(in) :: Ei	! [eV] energy levels
   real(8), dimension(:,:), intent(inout) :: DOS	! [eV] grid; [a.u.] DOS
   real(8), intent(in) :: smearing	! [eV] smearing used for DOS calculations
   real(8), dimension(:,:,:), intent(inout), optional :: partial_DOS    ! partial DOS made of each orbital type, if required to be constructed
   logical, dimension(:,:,:), intent(in), optional :: masks_DOS   ! partial DOS made of each orbital type, if required to be constructed
   real(8), dimension(:,:), intent(in), optional :: Hij      ! real eigenvectors
   complex, dimension(:,:), intent(in), optional :: CHij ! complex eigenvectors
   !-----------------------------------------
   integer :: i, Nsiz, Ngridsiz, j_center, j, N_at, N_types, i_at, i_types
   real(8) :: Gaus, epsylon, sigma, temp
   real(8), dimension(:), allocatable :: DOS_sum
   real(8), dimension(:,:,:), allocatable :: partial_DOS_sum
   logical :: do_partial
   
!    print*, 'get_DOS_sort test 0'
   
   epsylon = 1.0d-12	! precision
   sigma = smearing	! [eV] gaussian smearing
   Ngridsiz = size(DOS,2)	! number of grid points
   Nsiz = size(Ei)	! number of energy levels
   DOS(2,:) = 0.0d0	! to start from
   allocate(DOS_sum(Ngridsiz))
   DOS_sum = 0.0d0	! to start from
   do_partial = (present(partial_DOS) .and. present(masks_DOS) .and. (present(Hij) .or. present(CHij)))
   if (do_partial) then
      N_at = size(partial_DOS,1)
      N_types = size(partial_DOS,2)
      allocate(partial_DOS_sum(N_at, N_types, Ngridsiz))
      partial_DOS_sum = 0.0d0
   endif
   
!     print*, 'get_DOS_sort test 1', Ngridsiz
   
   !$omp PARALLEL private(i, j_center, j, Gaus, i_at, i_types, temp)
   !$omp do schedule(dynamic) reduction( + : DOS_sum, partial_DOS_sum)
   do i = 1, Ngridsiz	! for all grid points
      ! Do the summation in two parts:
!        print*, 'get_DOS_sort test 1.5'
      
      call Find_in_monotonous_1D_array(Ei, DOS(1,i), j_center)	! module "Little_subroutines"
      
!        print*, 'get_DOS_sort test 2'
      
      ! 1) Contribution from the levels above the chosen point:
      if (j_center <= Nsiz) then
         EL:do j = j_center, Nsiz	! for all energy levels above, up to the last one
            call Gaussian(Ei(j), sigma, DOS(1,i), Gaus)	! module "Little_subroutines"
            if (Gaus < epsylon) exit EL	! no need to continue, the contribution from higher levels is negligible
            DOS_sum(i) = DOS_sum(i) + Gaus
            if (do_partial) then
               if (present(Hij)) then
                  temp = SUM( Hij(:,j) * Hij(:,j) )
                  if (abs(temp) > 1.0d-12) then
                     do i_at = 1, N_at
                        do i_types = 1, N_types
                           !partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(Hij(:,j)*Hij(:,j)/temp, MASK = masks_DOS(i_at, i_types, :))
!                             print*, 'get_DOS_sort test 3a'
                           partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(Hij(:,j)*Hij(:,j), MASK = masks_DOS(i_at, i_types, :))/temp
!                             print*, 'get_DOS_sort test 4a'
                        enddo
                     enddo
                  endif
               elseif (present(CHij)) then
                  !temp = SUM( dconjg(CHij(:,j)) * CHij(:,j) )
                  temp = SUM( conjg(CHij(:,j)) * CHij(:,j) )
                  if (abs(temp) > 1.0d-12) then
                     do i_at = 1, N_at
                        do i_types = 1, N_types
                           !partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(dconjg(CHij(:,j))*CHij(:,j)/temp, MASK = masks_DOS(i_at, i_types, :))
!                            partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(dconjg(CHij(:,j))*CHij(:,j), MASK = masks_DOS(i_at, i_types, :))/temp
!                             print*, 'get_DOS_sort test 3b'
                           partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(conjg(CHij(:,j))*CHij(:,j), MASK = masks_DOS(i_at, i_types, :))/temp
!                             print*, 'get_DOS_sort test 4b'
                        enddo
                     enddo
                  endif
               endif
            endif
         enddo EL
      endif ! (j_center <= Nsiz)
      ! 2) Contribution from the levels below the given point:
      if (j_center > 1) then
         EL2:do j = (j_center-1), 1, -1	! for all energy levels below, down to the first one
            call Gaussian(Ei(j), sigma, DOS(1,i), Gaus)	! module "Little_subroutines"
            if (Gaus < epsylon) exit EL2	! no need to continue, the contribution from higher levels is negligible
            DOS_sum(i) = DOS_sum(i) + Gaus
            if (do_partial) then
               if (present(Hij)) then
                  temp = SUM( Hij(:,j) * Hij(:,j) )
                  if (abs(temp) > 1.0d-12) then
                     do i_at = 1, N_at
                        do i_types = 1, N_types
                           !partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(Hij(:,j)*Hij(:,j)/temp, MASK = masks_DOS(i_at, i_types, j))
!                            print*, 'get_DOS_sort test 3c'
                           partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(Hij(:,j)*Hij(:,j), MASK = masks_DOS(i_at, i_types, j))/temp
!                            print*, 'get_DOS_sort test 4c'
                        enddo
                     enddo
                  endif
               elseif (present(CHij)) then
!                   temp = SUM( dconjg(CHij(:,j)) * CHij(:,j) )
                  temp = SUM( conjg(CHij(:,j)) * CHij(:,j) )
                  if (abs(temp) > 1.0d-12) then
                     do i_at = 1, N_at
                        do i_types = 1, N_types
                           !partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(dconjg(CHij(:,j))*CHij(:,j)/temp, MASK = masks_DOS(i_at, i_types, :))
!                            partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(dconjg(CHij(:,j))*CHij(:,j), MASK = masks_DOS(i_at, i_types, :))/temp
!                            print*, 'get_DOS_sort test 3d'
                           partial_DOS_sum(i_at, i_types, i) = partial_DOS_sum(i_at, i_types, i) + Gaus*SUM(conjg(CHij(:,j))*CHij(:,j), MASK = masks_DOS(i_at, i_types, :))/temp
!                            print*, 'get_DOS_sort test 4d'
                        enddo
                     enddo
                  endif
               endif
            endif
         enddo EL2
      endif ! (j_center > 1)
   enddo !  i = 1, Ngridsiz
   !$omp end do
   !$omp end parallel
   
!     print*, 'get_DOS_sort test 5'
   
   DOS(2,:) = DOS_sum(:)
   if (do_partial) partial_DOS(:,:,:) = partial_DOS_sum(:,:,:)
   
   deallocate(DOS_sum)
   if (allocated(partial_DOS_sum)) deallocate(partial_DOS_sum)
!     print*, 'get_DOS_sort test 6'
end subroutine get_DOS_sort



subroutine get_total_el_energy(Scell, matter, numpar, t, Err) ! get total electron energy, from module "Electron_tools"
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), intent(in) :: t ! current timestep [fs]
   type(Error_handling), intent(inout) :: Err ! error save
   !=========================================
   real(8) nat
   integer NSC
   do NSC = 1, size(Scell)
      nat = dble(Scell(NSC)%Na) ! number of atoms
      select case (numpar%el_ion_scheme)
      case (1) ! Enforced energy conservation:
         Scell(NSC)%nrg%E_tot = Scell(NSC)%nrg%E_glob - Scell(NSC)%nrg%At_kin*nat - Scell(NSC)%nrg%E_rep - Scell(NSC)%nrg%E_supce*nat	! total electron energy
         Scell(NSC)%nrg%El_low = Scell(NSC)%nrg%E_tot - Scell(NSC)%nrg%El_high*nat ! energy of low-energy electrons [eV]
         ! Correspondingly update distribution function:
         call update_fe(Scell, matter, numpar, t, Err) ! module "Electron_tools"
         !call get_low_e_energy(Scell, matter) ! get the total electron energy
         ! Get the energy associated with the electrons populating band structure:
         if (numpar%scc) then ! energy for SCC is defined by for H_0:
            !call get_low_e_energy(Scell, matter, numpar) ! below
            call get_low_e_energy(Scell, matter) ! below
         else  ! no SCC, the usual expression then
            call get_low_e_energy(Scell, matter) ! below
         endif
         Scell(NSC)%nrg%Total = Scell(NSC)%nrg%At_pot + Scell(NSC)%nrg%At_kin + Scell(NSC)%nrg%E_vdW + &
               Scell(NSC)%nrg%E_coul + Scell(NSC)%nrg%E_expwall + Scell(NSC)%nrg%E_coul_scc/nat ! [eV/atom] initial total energy
      case default   ! Separate electronic and atomic energies:

         call update_fe(Scell, matter, numpar, t, Err) ! module "Electron_tools"

         !call get_low_e_energy(Scell, matter) ! get the total electron energy
         ! Get the energy associated with the electrons populating band structure:
         if (numpar%scc) then ! energy for SCC is defined by for H_0:
            !call get_low_e_energy(Scell, matter, numpar) ! below
            call get_low_e_energy(Scell, matter) ! below
         else  ! no SCC, the usual expression then
            call get_low_e_energy(Scell, matter) ! below
         endif

         Scell(NSC)%nrg%E_tot = Scell(NSC)%nrg%El_low + Scell(NSC)%nrg%El_high*nat ! energy of all electrons (low + high energies)
         Scell(NSC)%nrg%Total = Scell(NSC)%nrg%At_pot + Scell(NSC)%nrg%At_kin + Scell(NSC)%nrg%E_vdW + &
               Scell(NSC)%nrg%E_coul + Scell(NSC)%nrg%E_expwall + Scell(NSC)%nrg%E_coul_scc/nat ! [eV/atom] initial total energy
      end select
   enddo
end subroutine get_total_el_energy



subroutine get_glob_energy(Scell, matter)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter ! material parameters
   integer NSC
   do NSC = 1, size(Scell)
      call get_new_global_energy(Scell(NSC), Scell(NSC)%nrg)
   enddo
end subroutine get_glob_energy

subroutine get_new_global_energy(Scell, nrg)
   type(Super_cell), intent(in) :: Scell  ! supercell with all the atoms as one object
   type(Energies), intent(inout) :: nrg	! energies in the material
   real(8) :: nat
   nat = dble(Scell%Na)
   nrg%Total = nrg%At_pot + nrg%At_kin + nrg%E_vdW + nrg%E_coul + nrg%E_expwall + nrg%E_coul_scc/nat ! [eV/atom] initial total energy
   nrg%E_glob = (nrg%Total + nrg%E_supce)*nat ! [eV] total energy in the super-cell, save it
end subroutine get_new_global_energy



subroutine get_low_e_energy(Scell, matter, numpar)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter ! material parameters
   type(Numerics_param), intent(in), optional :: numpar  ! numerical parameters
   integer NSC
   do NSC = 1, size(Scell)
      if (present(numpar)) then  ! there may be SCC calculations involved
         if (numpar%scc) then ! SCC, so the total energy is defined by the part H_0 without charge energy:
            call set_total_el_energy(Scell(NSC)%Ei_scc_part, Scell(NSC)%fe, Scell(NSC)%nrg%El_low) ! below
         else ! non-SCC:
            call set_total_el_energy(Scell(NSC)%Ei, Scell(NSC)%fe, Scell(NSC)%nrg%El_low) ! below
         endif
      else
         call set_total_el_energy(Scell(NSC)%Ei, Scell(NSC)%fe, Scell(NSC)%nrg%El_low) ! below
      endif
   enddo
end subroutine get_low_e_energy



subroutine set_total_el_energy(Ei,fe,E_tot)
   real(8), dimension(:), intent(in) :: Ei	! energy levels, eigenvalues of the hamiltonian matrix [eV]
   real(8), dimension(:), intent(in) :: fe	! electron distribution
   real(8), intent(out) :: E_tot
   integer n
!    real(8) sumNe, sumEe
   n = size(Ei)
   E_tot = 0.0d0
   !E_tot = SUM(Ei(:)*fe(:))	! total electron energy, or do it with mkl subroutine (faster):
   call dgemm ('T','N', 1, 1, n, 1.0d0, Ei, n, fe, n, 0.0d0, E_tot, 1) ! mkl

!    sumNe = 0.0d0
!    sumEe = 0.0d0
!    todo:do n = 1, size(Ei)
!       sumNe = sumNe + fe(n)
!       sumEe = sumEe + Ei(n)*fe(n)
!       if (fe(n) <= 1.0d-14) then
!          print*, n, fe(n-1), fe(n), Ei(n), sumNe, (sumEe)/64.0d0, E_tot/64.0d0
!          exit todo
!       endif
!    enddo todo
end subroutine set_total_el_energy



!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
! Instant equilibration of electron distribution towards Fermi-Dirac distribution:

subroutine set_initial_fe(Scell, matter, Err, norm_fe)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter ! material parameters
   type(Error_handling), intent(inout) :: Err	! error save
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   !=========================================
   integer :: NSC, i
   real(8) :: sumNe
   real(8) :: f_norm
   character(200) :: Error_descript
   DO_TB:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
      if (present(norm_fe)) then    ! user provided
         f_norm = norm_fe
      else  ! by default, spin degenerate
         f_norm = 2.0d0
      endif
      Error_descript = ''
      do NSC = 1, size(Scell)
         ! Allocate electron distribution function (Fermi-function):
         if (.not. allocated(Scell(NSC)%fe)) allocate(Scell(NSC)%fe(size(Scell(NSC)%Ei)))
         ! Define the distribution function:
         TE:if (Scell(NSC)%Te <= 1.0d-12) then ! zero temperature
            sumNe = 0.0d0
            Scell(NSC)%mu = -1d15   ! just to start with
            do i = 1,size(Scell(NSC)%fe) ! all energy levels
               sumNe = sumNe + f_norm
               if (sumNe <= Scell(NSC)%Ne) then
                  Scell(NSC)%fe(i) = f_norm
               else
                  if (Scell(NSC)%mu <= -1d14) Scell(NSC)%mu = Scell(NSC)%Ei(i-1) ! chem.potential [eV]
                  Scell(NSC)%fe(i) = 0.0d0
               endif
!                print*, i, Scell(NSC)%Ei(i), Scell(NSC)%fe(i)
            enddo
!             pause 'set_initial_fe'
         else TE ! finite temperature
            ! Find chem.potential mu for given temperature Te:
            call Electron_fixed_Te(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! module "Electron_tools"
            ! Knowing the chem.potential, set the Fermi-distribution:
            if (present(norm_fe)) then
               call set_Fermi(Scell(NSC)%Ei, Scell(NSC)%TeeV, Scell(NSC)%mu, Scell(NSC)%fe, Error_descript, norm_fe)	! module "Electron_tools"
            else
               call set_Fermi(Scell(NSC)%Ei, Scell(NSC)%TeeV, Scell(NSC)%mu, Scell(NSC)%fe, Error_descript)	! module "Electron_tools"
            endif
            if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
               call Save_error_details(Err, 7, Error_descript)
               print*, trim(adjustl(Error_descript))
            endif
            ! Also set temperature in [K]:
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! [eV] -> [K]
         endif TE
      enddo
   endif DO_TB
end subroutine set_initial_fe


pure function Fermi_distribution(E, mu, T, norm_fe) result(fe)
   real(8) :: fe    ! fermi-function
   real(8), intent(in) :: E, mu, T
   real(8), intent(in), optional :: norm_fe ! with spin degeneracy or without
   real(8) :: f_norm
   
   if (present(norm_fe)) then    ! user provided
      f_norm = norm_fe
   else  ! by default, spin degenerate
      f_norm = 2.0d0
   endif
   
   fe = f_norm/(1.0d0 + exp((E - mu)/T))
end function Fermi_distribution 



subroutine set_Fermi(Ei,Te,mu,fe, Error_desript, norm_fe)
   real(8), dimension(:), intent(in) :: Ei	! energy levels, eigenvalues of the hamiltonian matrix
   real(8), intent(in) :: Te, mu	! temperature [eV], and chemical potential [eV]
   real(8), dimension(:), intent(out) :: fe
   character(*), intent(out), optional :: Error_desript ! description of error if occured
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer i, n
   real(8) :: f_norm
   
   if (present(norm_fe)) then    ! user provided
      f_norm = norm_fe
   else  ! by default, spin degenerate
      f_norm = 2.0d0
   endif
      
   if (Te > 0.0d0) then
      where ((Ei(:) - mu)/Te >= log(HUGE(mu))) ! exp(x) -> infinity
         fe(:) = 0.0d0
      else where
         fe(:) = f_norm/(1.0d0 + exp((Ei(:) - mu)/Te)) ! fermi-function
      end where
   elseif (Te == 0.0d0) then
      where (Ei(:) <= mu)
         fe(:) = f_norm
      else where
         fe(:) = 0.0d0
      end where
   else  !Te < 0
      fe(:) = 0.0d0
      if (present(Error_desript)) then
         Error_desript = 'Electron temperature is negative in set_Fermi-subroutine!'
      else
         !print*, 'Electron temperature is negative in set_Fermi-subroutine!'
      endif
   endif
end subroutine set_Fermi



subroutine set_Erf_distribution(Ei,Te,mu,fe, Error_desript, norm_fe)
   real(8), dimension(:), intent(in) :: Ei	! energy levels, eigenvalues of the hamiltonian matrix
   real(8), intent(in) :: Te, mu	! temperature [eV], and chemical potential [eV]
   real(8), dimension(:), intent(out) :: fe
   character(*), intent(out), optional :: Error_desript ! description of error if occured
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer i, n
   real(8) :: f_norm
   
   if (present(norm_fe)) then    ! user provided
      f_norm = norm_fe
   else  ! by default, spin degenerate
      f_norm = 2.0d0
   endif
   
   print*, 'set_Erf_distribution'
   
   
   if (Te > 0.0d0) then
      do i = 1,size(Ei)
         if ((Ei(i) - mu)/Te >= log(HUGE(mu))) then ! exp(x) -> infinity
            fe(i) = 0.0d0
         else
            fe(i) = (1.0d0-erf((Ei(i)-mu)/Te))   ! erf-function
!             print*, i, Ei(i), fe(i)
         endif
      enddo
   elseif (Te == 0.0d0) then
      where (Ei(:) <= mu)
         fe(:) = f_norm
      else where
         fe(:) = 0.0d0
      end where
   else  !Te < 0
      fe(:) = 0.0d0
      if (present(Error_desript)) then
         Error_desript = 'Electron temperature is negative in set_Erf_distribution subroutine!'
      else
         !print*, 'Electron temperature is negative in set_Erf_distribution subroutine!'
      endif
   endif
end subroutine set_Erf_distribution



subroutine get_electronic_heat_capacity(Scell, NSC, Ce, do_kappa, DOS_weights, Ce_partial, norm_fe)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), intent(out) :: Ce ! current electron heat capacity [J/(m^3 K)]
   logical, intent(in) :: do_kappa  ! if kappa calculations are requested
   real(8), dimension(:,:,:), intent(in), optional :: DOS_weights ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:), intent(out), allocatable, optional :: Ce_partial ! band-resolved Ce [J/(m^3 K)]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8) :: Ntot	! number of electrons
   real(8) :: nat   ! number of atoms
   real(8) :: Te	! current electron temperature [eV]
   real(8) :: mu	! current electron chemical potential [eV]
   real(8) :: dmu, dTe, mu0	! electron differential chemical potential [eV], temperature [eV]
   real(8) :: Dens	 ! atomic density
   real(8) :: coef   ! conversion coefficients with units
   real(8) :: C1, C2
   logical :: do_partial
   real(8), dimension(size(Scell(NSC)%Ei)) :: Ce_i

    if (present(DOS_weights) .and. present(Ce_partial)) then ! partial contributions required:
       do_partial = .true.
       if (.not.allocated(Ce_partial)) allocate(Ce_partial(size(Scell(NSC)%G_ei_partial,1)))
    else
       do_partial = .false.
    endif

   dTe = 10.0d0/g_kb	! [eV] -> [K]
   Ntot = dble(Scell(NSC)%Ne)
   nat = dble(Scell(NSC)%Na) ! number of atoms

   ! 1) Low-energy electrons populating TB-band structure:
   Te = Scell(NSC)%TeeV
   mu = Scell(NSC)%mu
   if (present(norm_fe)) then
      call Electron_Fixed_Te(Scell(NSC)%Ei, Ntot, mu0, Te+dTe, norm_fe) ! in case if the electron temperature is given
   else
      call Electron_Fixed_Te(Scell(NSC)%Ei, Ntot, mu0, Te+dTe) ! in case if the electron temperature is given
   endif
   dmu = (mu0 - mu)/dTe
   if (present(norm_fe)) then ! normalization of fe provided:
      if (do_partial) then
         call Get_Ce(Ce_i, Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, DOS_weights, Ce_partial, norm_fe)
      else  ! no partial contributions required:
         call Get_Ce(Ce_i, Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, norm_fe=norm_fe)
      endif
   else  ! default normalization of fe:
      if (do_partial) then
         call Get_Ce(Ce_i, Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, DOS_weights, Ce_partial)
      else  ! no partial contributions required:
         call Get_Ce(Ce_i, Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce)
      endif
   endif

   !Dens = Scell(NSC)%Ne_low/(Scell(NSC)%V)*1d24 ! [1/cm^3]
   !Dens = Dens*1d6   ! -> [1/m^3]
   !Ce = Ce/Scell(NSC)%Ne_low * Dens*g_e/g_kb	! [J/(m^3 K)]

   coef = 1.0d30*g_e/g_kb  ! [eV/A^3] -> [J/m^3/K]
   Dens = 1.0d0/(Scell(NSC)%V) ! [1/A^3]
   Ce = Ce * Dens*coef  ! [J/(m^3 K)]
   if (do_partial) then
      Ce_partial = Ce_partial * Dens*coef  ! [J/(m^3 K)]
   endif

   ! Save energy-level-resolved heat capacities, if needed:
   if (do_kappa) then
      Scell(NSC)%Ce_i = Ce_i/g_kb   ! [eV/K]
   endif

   ! 2) High-energy electrons from MC:
   C2 = Scell(NSC)%Ne_high * Dens*coef
   Ce = Ce + C2

   if (do_partial) then ! For simplicity, distribute the high-energy electrons equally:
      Ce_partial = Ce_partial + C2/dble(size(Ce_partial))
   endif

   if (isnan(Ce) .or. abs(Ce) >= 1d30) Ce = 0.0d0 ! if undefined or infinite

end subroutine get_electronic_heat_capacity


subroutine Get_Ce(Ce_i, Ei, Te, mu, dmu, C, DOS_weights, Ce_partial, norm_fe)
    real(8), dimension(:), intent(out) :: Ce_i
    real(8), dimension(:), intent(in) :: Ei
    real(8), intent(in) :: Te, mu, dmu
    real(8), intent(out) :: C ! heat capacity
    real(8), dimension(:,:,:), intent(in), optional :: DOS_weights ! weigths of the particular type of orbital on each energy level
    real(8), dimension(:), intent(out), optional :: Ce_partial ! band-resolved Ce [J/(m^3 K)]
    real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
    !------------------------
    real(8) :: dfdT, E, C_temp
    integer :: i, N, N_at, N_types, i_at, i_types, i_G1
    logical :: do_partial

    if (present(DOS_weights) .and. present(Ce_partial)) then ! partial contributions required:
       do_partial = .true.
       Ce_partial = 0.0d0 ! to start with
    else
       do_partial = .false.
    endif

    if (do_partial) then
      N_at = size(DOS_weights,1)    ! number of kinds of atoms
      N_types = size(DOS_weights,2) ! number of atomic shells (basis set size)
    endif
    N = size(Ei)
    C = 0.0d0
    do i = 1, N   ! all energy levels
        E = Ei(i)
        if (present(norm_fe)) then
           dfdT = Diff_Fermi_Te(Te, mu, dmu, E, norm_fe)
        else
           dfdT = Diff_Fermi_Te(Te, mu, dmu, E)
        endif
        C_temp = dfdT*(E-mu)  ! correct definition from Cv = T*dS/dT; S=entropy
        C = C + C_temp     ! total
        Ce_i(i) = C_temp   ! save for output

        if (do_partial) then ! partial contributions required:
           do i_at = 1, N_at ! all elements
              do i_types = 1, N_types  ! all shells of each element
                 i_G1 = (i_at-1) * N_types + i_types
                 Ce_partial(i_G1) = Ce_partial(i_G1) + C_temp * DOS_weights(i_at, i_types, i)
              enddo ! i_types
           enddo ! i_at
        endif ! do_partial
    enddo
end subroutine


pure function Diff_Fermi_Te(Te, mu, dmu, E, norm_fe)
   real(8), intent(in) :: Te, mu, dmu, E   ! [eV], temperature, chem.potential, and energy
   real(8) :: Diff_Fermi_Te   ! Derivative of the Fermi-function by temperature Te
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8) F, buf
   real(8) :: f_norm
   
   if (present(norm_fe)) then    ! user provided
      f_norm = norm_fe
   else  ! by default, spin degenerate
      f_norm = 2.0d0
   endif
   
   if (((E - mu)/Te) >= log(HUGE(mu))) then ! dealing with the problem of large and small numbers
      Diff_Fermi_Te = 0.0d0
   else
      buf = dexp((E - mu)/Te)
      if ( buf .gt. 1.0d30) then ! dealing with the problem of large and small numbers
         F = f_norm/(buf)
         Diff_Fermi_Te = F*(E - mu + Te*dmu)/(Te*Te)
      else	! in case everything is ok
         F = 1.0d0/(1.0d0 + buf)
         Diff_Fermi_Te = f_norm*buf*F*F*(E - mu + Te*dmu)/(Te*Te)
      endif
   endif
end function Diff_Fermi_Te



pure function Diff_Fermi_E(Te, mu, E, norm_fe) result(dfdE)
   real(8), intent(in) :: Te, mu, E   ! [eV], temperature, chem.potential, and energy
   real(8) :: dfdE   ! Derivative of the Fermi-function by energy
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8) F, buf
   real(8) :: f_norm

   if (present(norm_fe)) then    ! user provided
      f_norm = norm_fe
   else  ! by default, spin degenerate
      f_norm = 2.0d0
   endif

   if (((E - mu)/Te) >= log(HUGE(mu))) then ! dealing with the problem of large and small numbers
      dfdE = 0.0d0
   else
      buf = dexp((E - mu)/Te)
      if ( buf > 1.0d30) then ! dealing with the problem of large and small numbers
         dfdE = -f_norm/(buf*Te)
      else  ! in case everything is ok
         F = 1.0d0/(1.0d0 + buf)
         dfdE = -f_norm*buf*F*F/Te
      endif
   endif
end function Diff_Fermi_E



subroutine Electron_Fixed_Te(wrD, Netot, mu, Te, norm_fe) ! in case if the electron temperature is given
! but the total energy can change, then we have to find the chem.potential that conserves 
! the given total number of electons Netot
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/atom to normalize the distribution function
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(out) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8) a, b, Ncur
   integer i

   a = wrD(1) - 10.0d0
   b = wrD(size(wrD)) + 5.0d0
   Ncur = 0.0d0

   do while (ABS(Ncur-Netot) .GT. Netot*1d-12) !
      mu = (a+b)/2.0d0
      Ncur = get_N_tot(wrD, mu, Te) ! function from below
      if (Ncur .GT. Netot) then
         b = mu
      else
         a = mu
      endif
      if (ABS(a-b) .LT. 1d-12) exit ! it's too close anyway...
   enddo ! while
end subroutine Electron_Fixed_Te



function get_N_partial(fe, i_start, i_end) result(Ne)
   real(8) Ne  ! number of electrons in the given inteval according to distribution function
   real(8), dimension(:), intent(in) :: fe   ! given distribution function
   integer, intent(in) :: i_start, i_end  ! starting and ending levels to include
   !----------------------
   Ne = SUM(fe(i_start:i_end))   ! total number of electrons in the givel interval
end function get_N_partial

function get_E_partial(wr, fe, i_start, i_end) result(Ee)
   real(8) Ee  ! energy of electrons in the given inteval according to distribution function
   real(8), dimension(:), intent(in) :: wr, fe   ! given energy levels and distribution function
   integer, intent(in) :: i_start, i_end  ! starting and ending levels to include
   !----------------------
   Ee = SUM(wr(i_start:i_end) * fe(i_start:i_end))   ! total energy of electrons in the givel interval
end function get_E_partial


function get_N_partial_Fermi(wr, mu, Te, i_start, i_end) result(Ne)
   real(8) Ne  ! number of electrons in the given inteval according to Fermi distribution
   real(8), dimension(:), intent(in) :: wr   ! given energy levels
   real(8), intent(in) :: Te ! electron temperature [eV]
   real(8), intent(in) ::  mu ! chem.potential [eV]
   integer, intent(in) :: i_start, i_end  ! starting and ending levels to include
   !----------------------
   if (Te > 1.0d-6) then
      Ne = 2.0d0 * SUM(1.0d0/(1.0d0 + dexp((wr(:) - mu)/Te))) ! Fermi-function
   else
      Ne = 2.0d0 * dble(COUNT(wr(i_start:i_end) <= mu)) ! Fermi-function at T=0
   endif
end function get_N_partial_Fermi

function get_E_partial_Fermi(wr, mu, Te, i_start, i_end) result(Ee)
   real(8) Ee  ! energy of electrons in the given inteval according to Fermi distribution
   real(8), dimension(:), intent(in) :: wr   ! given energy levels
   real(8), intent(in) :: Te ! electron temperature [eV]
   real(8), intent(in) ::  mu ! chem.potential [eV]
   integer, intent(in) :: i_start, i_end  ! starting and ending levels to include
   !----------------------
   if (Te > 1.0d-6) then
      Ee = 2.0d0 * SUM(wr(i_start:i_end)/(1.0d0 + dexp((wr(i_start:i_end) - mu)/Te))) ! Fermi-function
   else
      Ee = 2.0d0 * SUM(wr(i_start:i_end)) ! Fermi-function at T=0
   endif
end function get_E_partial_Fermi




function get_N_tot(wrD, mu, Te, norm_fe)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(in) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8) :: get_N_tot ! number of electrons for given band structure and Fermi-function
   integer i
   real(8) :: f_norm
   
   if (Te .GT. 1.0d-12) then
      !get_N_tot = SUM(2.0d0/(1.0d0 + dexp((wrD(:) - mu)/Te))) ! fermi-function
      !get_N_tot = SUM(2.0d0/(1.0d0 + dexp((wrD(:) - mu)/Te)), mask = (wrD(:) - mu)/Te .LT. log(HUGE(Te))) ! fermi-function
      get_N_tot = 0.0d0
      DO_SUM:do i = 1, size(wrD)
         if ((wrD(i) - mu)/Te .LT. log(HUGE(mu))) then ! exp(x) -> infinity
            !get_N_tot = get_N_tot + (2.0d0/(1.0d0 + dexp((wrD(i) - mu)/Te))) ! fermi-function
            if (present(norm_fe)) then    ! user provided
               get_N_tot = get_N_tot + Fermi_distribution (wrD(i), mu, Te, norm_fe)  ! function above
            else    ! default
               get_N_tot = get_N_tot + Fermi_distribution (wrD(i), mu, Te)  ! function above
            endif
         endif
      enddo DO_SUM
   else
      if (present(norm_fe)) then    ! user provided
         f_norm = norm_fe
      else  ! by default, spin degenerate
         f_norm = 2.0d0
      endif
      get_N_tot = f_norm*COUNT(wrD(:) .LT. mu)
   endif
end function get_N_tot

function get_E_tot(wrD, mu, Te, norm_fe)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(in) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8) :: get_E_tot ! number of electrons for given band structure and Fermi-function
   integer i
   real(8) :: f_norm
   if (Te .GT. 1.0d-12) then
      !get_E_tot = 2.0d0*SUM(wrD(:)/(1.0d0 + exp((wrD(:) - mu)/Te)))
      get_E_tot = 0.0d0
      DO_SUM:do i = 1, size(wrD)
         if ((wrD(i) - mu)/Te .LT. log(HUGE(mu))) then ! exp(x) -> infinity
            !get_E_tot = get_E_tot + 2.0d0*wrD(i)/(1.0d0 + dexp((wrD(i) - mu)/Te)) ! fermi-function
            if (present(norm_fe)) then    ! user provided
               get_E_tot = get_E_tot + wrD(i)*Fermi_distribution (wrD(i), mu, Te, norm_fe)  ! function above
            else    ! default
               get_E_tot = get_E_tot + wrD(i)*Fermi_distribution (wrD(i), mu, Te)  ! function above
            endif
         endif
      enddo DO_SUM
   else
      if (present(norm_fe)) then    ! user provided
         f_norm = norm_fe
      else  ! by default, spin degenerate
         f_norm = 2.0d0
      endif
!       get_E_tot = 2.0d0*SUM(wrD(:), mask = wrD(:) .LT. mu)
      get_E_tot = f_norm*SUM(wrD(:), mask = wrD(:) .LT. mu)
   endif
end function get_E_tot



subroutine find_mu_from_N_T(wrD, Netot, mu, Te) ! in case if the electron temperature is given
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/atom to normalize the distribution function
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(out) ::  mu ! chem.potential to be found [eV]
   call Electron_Fixed_Te(wrD, Netot, mu, Te)
end subroutine find_mu_from_N_T

subroutine find_T_from_N_mu(wrD, Netot, mu, Te) ! in case if the electron temperature is given
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/atom to normalize the distribution function
   REAL(8), INTENT(out) :: Te ! electron temperature [eV]
   REAL(8), INTENT(in) ::  mu ! chem.potential to be found [eV]
   real(8) a, b, Ncur
   integer i
   a = 0.0d0
   b = 1.0d5
   Ncur = 0.0d0
   do while (ABS(Ncur-Netot)/Netot .GT. 1d-12) !
      Te = (a+b)/2.0d0
      Ncur = get_N_tot(wrD, mu, Te)
      if (Ncur .GT. Netot) then
         b = Te
      else
         a = Te
      endif
      if (ABS(a-b) .LT. 1d-12) exit ! it's too close anyway...
   enddo ! while
end subroutine find_T_from_N_mu



subroutine Electron_Fixed_Etot(wrD, Netot, Eetot, mu, Te, subnum)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(in) :: Eetot  ! energy of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]
   logical, intent(in), optional :: subnum ! which subroutine to use?
   if (present(subnum)) then ! do the new one
!       call Electron_Fixed_Etot_1(wrD, Netot, Eetot, mu, Te) ! default choice
!       call Electron_Fixed_Etot_2(wrD, Netot, Eetot, mu, Te) ! DO NOT USE, DOES NOT CONVERGE
      call Electron_Fixed_Etot_3(wrD, Netot, Eetot, mu, Te)   ! faster version
   else ! do the old one
      call Electron_Fixed_Etot_1(wrD, Netot, Eetot, mu, Te)
   endif
end subroutine Electron_Fixed_Etot



subroutine Electron_Fixed_Etot_2(wrD, Netot, Eetot, mu, Te) ! DOES NOT CONVERGE
   ! Newton method of iterative solution of system of equations:
   ! https://www.cmu.edu/math/undergrad/suami/pdfs/2014_newton_method.pdf
   !--------------------------------------------------------------------
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(in) :: Eetot  ! energy of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]
   !--------------------------------------------------------------------
   real(8) :: muN, TeN, dN_dmu, dN_dTe, dE_dmu, dE_dTe, mu_cur, Te_cur
   real(8) :: f1, f2, fermi, Efermi, exp_func(size(wrD)), exp_fermi(size(wrD))
   integer :: i
   logical :: presision_reached
   
   ! To start iterations from:
   muN = mu 
   TeN = Te
   presision_reached = .false.
   
   ! Do the iterations:
   i = 0
!    print*, 'i', i, mu_cur, muN, Te_cur, TeN
   do while ( .not.presision_reached )
      i = i + 1 ! count iterations
      
      ! get integrals of Fermi function:
      if (TeN > 1.0d-10) then
         where ((wrD(:) - muN) < TeN*log(HUGE(muN)))
            exp_func(:) = exp((wrD(:) - muN)/TeN)
         elsewhere
            exp_func(:) = 0.0d0
         endwhere
      else
         exp_func(:) = 0.0d0
      endif
      exp_fermi(:) = 2.0d0 / (1.0d0 + exp_func(:))
      ! Get Newton functions:
      f1 = SUM(exp_fermi(:)) - Netot 
      f2 = SUM(exp_fermi(:) * wrD(:)) - Eetot
      ! get the values of dN/dmu, dN/dTe, dE/dmu, dE/dTe for given mu and Te:
      call get_dN_dE_dmu(wrD, muN, TeN, dN_dmu, dE_dmu) ! below
      call get_dN_dE_dTe(wrD, muN, TeN, dN_dTe, dE_dTe) ! below
      ! Iterative solution:
      mu_cur = muN + (dN_dmu * f1 + dN_dTe * f2)
      Te_cur = TeN + (dE_dmu * f1 + dE_dTe * f2)
!       print*, i, dN_dmu, dE_dmu
!       print*, 'i', dN_dTe, dE_dTe
!       print*, 'f', f1, f2
!       print*, 'n', Netot, Eetot
!       print*, 'l', muN, TeN
!       print*, 'c', mu_cur, Te_cur
      
      ! Estimate if presision is reached:
      if ( (abs(mu_cur - muN) < 1.0-2 * abs(mu_cur)) .and. (abs(Te_cur - TeN) < 1.0-3 * abs(Te_cur)) ) then
         presision_reached = .true.
      endif
      ! Update old data:
      muN = mu_cur
      TeN = Te_cur
      if (i >= 100) exit ! enough iterations, it either converged enough or diverged
   enddo
   ! Save output:
   Te = TeN
   mu = muN
end subroutine Electron_Fixed_Etot_2


pure subroutine get_dN_dE_dmu(wrD, mu, Te, dN_dmu, dE_dmu)
   real(8), intent(out) :: dN_dmu, dE_dmu
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]
   real(8), dimension(size(wrD)) :: exp_Fermi, exp_func
   if (Te .GT. 1.0d-12) then    ! nonzero temperature
      where ((wrD(:) - mu) < Te*log(HUGE(mu)))
         exp_Fermi(:) = exp((wrD(:) - mu)/Te)
      elsewhere
         exp_Fermi(:) = 0.0d0
      endwhere
      where (exp_Fermi(:) < sqrt(HUGE(mu)))
         exp_func(:) = exp_Fermi(:) / (1.0d0 + exp_Fermi(:))**2
      elsewhere
         exp_func(:) = 0.0d0
      endwhere
      dN_dmu = 2.0d0 * SUM( exp_func(:) ) / Te  ! 2 from spin
      dE_dmu = 2.0d0 * SUM( exp_func(:) * wrD(:) ) / Te
   else     ! zero temperature
      dN_dmu = 0.0d0    ! analytical limit dN/dmu at Te->0
      dE_dmu = 0.0d0
   endif   
end subroutine get_dN_dE_dmu



pure subroutine get_dN_dE_dTe(wrD, mu, Te, dN_dTe, dE_dTe)
   real(8), intent(out) :: dN_dTe, dE_dTe
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]
   real(8), dimension(size(wrD)) :: exp_Fermi, exp_func
   real(8) :: Te2
   if (Te .GT. 1.0d-10) then    ! nonzero temperature
      where ((wrD(:) - mu) < Te*log(HUGE(mu)))
         exp_Fermi(:) = exp((wrD(:) - mu)/Te)
      elsewhere
         exp_Fermi(:) = 0.0d0
      endwhere
      where (exp_Fermi(:) < sqrt(HUGE(mu)))
         exp_func(:) = exp_Fermi(:) * (wrD(:) - mu) / (1.0d0 + exp_Fermi(:))**2
      elsewhere
         exp_func(:) = 0.0d0
      endwhere
      Te2 = 1.0d0 / Te**2
      dN_dTe = 2.0d0 * SUM( exp_func(:) ) * Te2 ! 2 from spin
      dE_dTe = 2.0d0 * SUM( exp_func(:) * wrD(:) ) * Te2
   else     ! zero temperature
      dN_dTe = 0.0d0    ! analytical limit dN/dTe at Te->0
      dE_dTe = 0.0d0
   endif
end subroutine get_dN_dE_dTe




subroutine Electron_Fixed_Etot_1(wrD, Netot, Eetot, mu, Te)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(in) :: Eetot  ! energy of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]

!    real(8), dimension(size(wrD)) :: fe	! transient Fermi-distribution
   real(8) a, b, Ncur, Ecur, Elast, muN, muE, Telast, Telast2, popul, che1, che2
   real(8) curMu, minmunmue, dT
   integer i, j, cou, countr
   !parameter (kb = 11604.0d0)   ! Boltzmann constant    [K/eV]

   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
   ! precalculations before finding mu and Te
   dT = 20.0d0	! temperature step [K]
   Te = dT/g_kb ! [eV] electron temperature to be calculated
   call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te nad Ne
   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

   !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
   ! finding mu and Te
   countr = 0
437 cou = 0 
   Telast = Te
   Telast2 = Telast + dT/g_kb
   Te = (Telast+Telast2)/2.0d0
   Ecur = -1d12
   che1 = (Ecur-Eetot)  ! checker for finding the crossing point
   do while ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-9)) .GT. 1d-9) .AND. (ABS((Te - Telast)/MAX(Te,Telast,1d-12)) .GT. 1d-12)) ! main while
!       write(*,'(a,f,f,f,f)') 'T:', Ecur, Eetot, Te, muN
      Te = (Telast+Telast2)/2.0d0 ! [eV]
      if (ABS((Te - Telast)/MAX(Te,Telast,1d-12)) .LT. 1d-12) exit

      ! Finding chem.potential for the neq temperature:
      call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te nad Ne

      ! Checking, whether this chem.potential and temperature are also solutions for energy: 
      Ecur = 0.0d0
      Ecur = get_E_tot(wrD, muN, Te) ! function from above

      ! if they are, we can stop the calculations:
      if (ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-9)) .LT. 1d-9) exit

      ! Otherwise, we have to continue for a new temperature:
      che2 = (Ecur-Eetot)  ! checker for finding the crossing point

      if (((che1*che2) .GE. 0.0d0) .AND. (cou .EQ. 0)) then ! they are of the same sign => no crossing point in between
         Telast = Telast2
         Telast2 = Telast2 + dT/g_kb ! here step must be reduced for fast-oscillating functions muN(Te) and/or muE(Te)
      else ! if there is a crossing point in between, they change the sign:
         if (cou .EQ. 0) then 
            Telast = Telast - dT/g_kb
            cou = 1 ! counter, to never repeat this action more than one time
         else ! lets find then the crossing point precisely with besiction method:
            if (che2 .GE. 0.0d0) then
               Telast2 = Te
            else
               Telast = Te
            endif
         endif
      endif
      Te = (Telast+Telast2)/2.0d0 ! [eV] save for output
      mu = muN ! [eV] save for output
      che1 = che2
   enddo ! main while

   if ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-5)) .GT. 1d-5) .AND. (countr .LT. 100))then
      countr = countr + 1
      Te = Te + 2.0d0*dT/g_kb
      goto 437 ! just try it over again if it didn't work well the first time...
   endif
!    if (countr .GE. 100) pause 'PAUSE HERE'
   !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

end subroutine Electron_Fixed_Etot_1




subroutine Electron_Fixed_Etot_3(wrD, Netot, Eetot, mu, Te)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(in) :: Eetot  ! energy of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]

!    real(8), dimension(size(wrD)) :: fe	! transient Fermi-distribution
   real(8) a, b, Ncur, Ecur, Elast, muN, muE, Telast, Telast2, popul, che1, che2
   real(8) curMu, minmunmue, dT
   integer i, j, cou, countr
   !parameter (kb = 11604.0d0)   ! Boltzmann constant    [K/eV]

   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
   ! precalculations before finding mu and Te
   dT = max(Te/10.0d0, 20.0d0)	! temperature step [K]
   Te = max(Te/2.0d0, dT/g_kb) ! [eV] electron temperature to be calculated
   call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te nad Ne
   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

   !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
   ! finding mu and Te
   countr = 0
438 cou = 0 
   Telast = Te
   Telast2 = Telast + dT/g_kb
   Te = (Telast+Telast2)/2.0d0
   Ecur = -1d12
   che1 = (Ecur-Eetot)  ! checker for finding the crossing point
   do while ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-9)) .GT. 1d-9) .AND. (ABS((Te - Telast)/MAX(Te,Telast,1d-12)) .GT. 1d-12)) ! main while
!       write(*,'(a,f,f,f,f)') 'T:', Ecur, Eetot, Te, muN
      Te = (Telast+Telast2)/2.0d0 ! [eV]
      if (ABS((Te - Telast)/MAX(Te,Telast,1d-12)) .LT. 1d-12) exit

      ! Finding chem.potential for the neq temperature:
      call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te nad Ne

      ! Checking, whether this chem.potential and temperature are also solutions for energy: 
      Ecur = 0.0d0
      Ecur = get_E_tot(wrD, muN, Te) ! function from above

      ! if they are, we can stop the calculations:
      if (ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-9)) .LT. 1d-9) exit

      ! Otherwise, we have to continue for a new temperature:
      che2 = (Ecur-Eetot)  ! checker for finding the crossing point

      if (((che1*che2) .GE. 0.0d0) .AND. (cou .EQ. 0)) then ! they are of the same sign => no crossing point in between
         Telast = Telast2
         Telast2 = Telast2 + dT/g_kb ! here step must be reduced for fast-oscillating functions muN(Te) and/or muE(Te)
      else ! if there is a crossing point in between, they change the sign:
         if (cou .EQ. 0) then 
            Telast = Telast - dT/g_kb
            cou = 1 ! counter, to never repeat this action more than one time
         else ! lets find then the crossing point precisely with besiction method:
            if (che2 .GE. 0.0d0) then
               Telast2 = Te
            else
               Telast = Te
            endif
         endif
      endif
      Te = (Telast+Telast2)/2.0d0 ! [eV] save for output
      mu = muN ! [eV] save for output
      che1 = che2
   enddo ! main while

   if ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-5)) .GT. 1d-5) .AND. (countr .LT. 10))then
      countr = countr + 1
      Te = Te + 2.0d0*dT/g_kb
      goto 438 ! just try it over again if it didn't work well the first time...
   endif
end subroutine Electron_Fixed_Etot_3




END MODULE Electron_tools
