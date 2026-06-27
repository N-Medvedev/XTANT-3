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
MODULE Electron_tools
use Universal_constants
use Objects
use Algebra_tools, only : Two_Vect_Matr
use Little_subroutines, only : Find_in_array_monoton, Fermi_interpolation, linear_interpolation, print_progress
use MC_cross_sections, only : TotIMFP, Mean_free_path
use Electron_electron_scattering, only: get_Boltzmann_alpha_beta, Boltzmann_solution, test_change_of_fe, &
                                        Electron_electron_scattering_Kij


#ifdef MPI_USED
   use MPI_subroutines, only : MPI_barrier_wrapper, do_MPI_Allreduce
#endif




implicit none
PRIVATE

public :: get_low_e_energy, find_band_gap, Diff_Fermi_E, get_number_of_CB_electrons, set_Fermi
public :: set_Erf_distribution, update_fe, Electron_thermalization, get_glob_energy, update_cross_section
public :: Do_relaxation_time, set_initial_fe, find_mu_from_N_T, set_total_el_energy, Electron_Fixed_Etot
public :: get_new_global_energy, get_electronic_heat_capacity, get_total_el_energy, electronic_entropy
public :: get_low_energy_distribution, set_high_DOS, get_Ce_and_mu, Diff_Fermi_Te, get_orbital_resolved_data
public :: patch_distribution, identify_fragment_for_orbital, get_fragments_data_for_electrons

 contains



subroutine get_low_energy_distribution(Scell, numpar) ! On the grid, if requested
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar   ! numerical parameters, including lists of earest neighbors
   !--------------------
   integer :: Nsiz, Nei, i, j, Nlev

   if (numpar%save_fe_grid) then  ! only user requested
      ! Count number of steps, over which the distribution is averaged:
      numpar%fe_aver_num = numpar%fe_aver_num + 1

      ! Recalculate the low-energy part of the distribution:
      if (numpar%save_fe_grid) then  ! only user requested
         Nsiz = size(Scell%E_fe_grid)
         Nei = size(Scell%Ei)
         ! Fill in the distribution function:
         i = 1 ! to start with
         SRTL:do j = 1, Nsiz ! all energy intervals
            Nlev = 0 ! count levels within this grid interval
            do while (Scell%Ei(i) < Scell%E_fe_grid(j)) ! all energy levels below the given one in this step
               Scell%fe_on_grid(j) = Scell%fe_on_grid(j) + Scell%fe(i)
               !print*, j, i, Scell%E_fe_grid(j), Scell%fe_on_grid(j), Scell%fe(i)
               Nlev = Nlev + 1
               Scell%fe_norm_on_grid(j) = Scell%fe_on_grid(j)/Nlev
               i = i + 1
               if (i > Nei) exit SRTL
            enddo
         enddo SRTL
      endif
   endif
end subroutine get_low_energy_distribution


subroutine set_high_DOS(Scell, numpar)
   type(Super_cell), intent(in) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar ! numerical parameters, including MC energy cut-off
   !-----------------------
   integer :: Nsiz, i
   real(8) :: eps, k, coef, prefac

   eps = 1.0d-12
   Nsiz = size(numpar%high_DOS)

   if (numpar%high_DOS(Nsiz) <= eps) then ! DOS is not set, so set it:
      k = 2.0d0*g_me/g_h**2
      coef = g_e*sqrt(g_e)/(1.0d30*2.0d0*g_Pi**2)
      prefac = coef * k**(1.5d0) * Scell%V
      do i = 1, Nsiz
         numpar%high_DOS(i) = prefac * sqrt( abs(Scell%E_fe_grid(i)) )  ! Free-electron DOS [1/(A^3 * eV)]
         !print*, i, Scell%E_fe_grid(i), numpar%high_DOS(i)
      enddo
   endif
end subroutine set_high_DOS



subroutine update_cross_section(Scell, matter, Te_in)
   type(Super_cell), intent(in) :: Scell  ! supercell with all the atoms as one object
   type(solid), intent(inout) :: matter   ! materil parameters
   real(8), intent(in), optional :: Te_in ! user-defined temperature [K]
   !------------------------
   integer :: Nshl, i, N_Te, j, N_Tmax
   real(8) :: dT, T_left, Te, T_right
   ! Get the mean free paths vs Te:
   Nshl = size(matter%Atoms(1)%Ip)
   select case (matter%Atoms(1)%TOCS(Nshl)) ! Valence band and CDF only
   case (1) ! CDF

      if (present(Te_in)) then   ! provided
         Te = Te_in
      else  ! default
         Te = Scell%Te
      endif

      if (Te > 100.0d0) then  ! recalculate the VB/CB scattering:
         ! Temperature (grid defined in subroutine get_MFPs as Te_temp = dble((i-1)*1000)):
         dT = 1000.0d0 ! [K] grid step
         T_left = floor((Te+1.0d-8)/1000.0d0)*1000.0d0
         N_Te = ceiling((Te+1.0d-8)/1000.0d0)
         T_right = dble(N_Te)*1000.0d0
         N_Tmax = size(matter%Atoms(1)%El_MFP_vs_T)
         ! Interpolate valence band MFP for the given temperature:
         if (N_Te == 1) then ! extrapolate to lower Te
            matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) - &
               (matter%Atoms(1)%El_MFP_vs_T(N_Te+1)%L(:) - matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:))/dT * (T_right-Te) ! [1/A]
         elseif (N_Te > N_Tmax) then  ! etrapolate to higher Te
            matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Tmax)%L(:) + &
               (matter%Atoms(1)%El_MFP_vs_T(N_Tmax)%L(:) - matter%Atoms(1)%El_MFP_vs_T(N_Tmax-1)%L(:))/dT * (Te-T_left) ! [1/A]
         else ! within the array of data
            if (Te > T_left + 0.5d0*dT) then ! from above
               matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) - &
                  (matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) - matter%Atoms(1)%El_MFP_vs_T(N_Te-1)%L(:))/dT * (T_right-Te) ! [1/A]
            else ! from below
               matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te-1)%L(:) + &
                  (matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) - matter%Atoms(1)%El_MFP_vs_T(N_Te-1)%L(:))/dT * (Te-T_left) ! [1/A]
            endif ! (Te > T_left + 0.5d0*dT)

!             print*, 'N_Te', Te, N_Te, T_left, matter%Atoms(1)%El_MFP(Nshl)%L(10), matter%Atoms(1)%El_MFP_vs_T(N_Te-1)%L(10), matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(10)
!             pause

         endif ! (N_Te == 1)

      else ! no need to recalculate, the tempereature is too small:
         N_Te = 1
         matter%Atoms(1)%El_MFP(Nshl)%L(:) = matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) ! [1/A]
      endif ! (Te > 100.0d0)


      !-------------------------------------------
      ! Update the total inelastic mean free path:
      matter%El_MFP_tot%L(:) = 0.0d0   ! to start with
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         do j = 1, Nshl ! for all shells of this atom
            matter%El_MFP_tot%L(:) = matter%El_MFP_tot%L(:) + matter%Atoms(i)%El_MFP(j)%L(:) ! [1/A] inverse MFP
         enddo
      enddo
   endselect

!    print*, 'Te=', Scell%Te, matter%Atoms(1)%El_MFP(Nshl)%L(1), matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(1)
end subroutine update_cross_section


subroutine find_band_gap(wr, Scell, matter, numpar)
   REAL(8), DIMENSION(:), INTENT(in) ::  wr	! [eV] energy levels
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(solid), intent(inout) :: matter	! materil parameters
   type(Numerics_param), intent(inout) :: numpar ! numerical parameters, including MC energy cut-off
   !-----------
   real(8) :: L, Ele, Egap_old
   integer :: i, sumNe, siz, j, Nshl, N_grid, k

   siz = size(wr)
   Nshl = size(matter%Atoms(1)%Ip)  ! last shell, corresponding to bandgap
   
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

   Egap_old = Scell%E_gap  ! store old gap to use below for MC cross section

   if (i < size(wr)) then
      Scell%E_gap = ABS(wr(i+1) - wr(i))	! bandgap [eV]
      Scell%E_bottom = wr(i+1)	! bottom of the CB [eV]
   else ! undefined, so set the closest guess
      Scell%E_gap = 0.0d0
      Scell%E_bottom = wr(i)
   endif
   ! If top energy level is excluded (in parameterization, or just too high due to convergence issues)
   j = siz  ! start from the last one
   Scell%E_top = wr(j)		! [eV] current top of the (meaningful) conduction band
   do while (wr(j)>=79.0d0)
      j = j - 1
      Scell%E_top = wr(j)		! [eV] current top of the (meaningful) conduction band
   enddo
   Scell%E_VB_bottom = wr(1)	! [eV] current bottom of the valence band
   Scell%E_VB_top = wr(i)		! [eV] current top of the valence band
   
   !--------------------------
   ! Set MC high-energy electron cut-off energy equal to the uppermost level of CB:
   if (numpar%E_cut_dynamic) numpar%E_cut = Scell%E_top - Scell%E_bottom ! [eV]

   ! Ionization potential of the valence band:
   matter%Atoms(1)%Ip(Nshl) = Scell%E_gap ! [eV]

   ! Update VB impact-ionization scattering CS:
   call update_ionization_CS(matter, Nshl, Scell%E_gap - Egap_old, numpar)  ! below

   ! For noneuqilibrium distributions (BO or relaxation time), threshold cannot be higher than
   ! the topmost level of CB, otherwise, there is no way to place an incomming electron:
   select case (numpar%el_ion_scheme)
   case (3:5)
      if ( (numpar%verbose) .and. (numpar%E_cut > (Scell%E_top-Scell%E_bottom) ) ) then
         if (numpar%MPI_param%process_rank == 0) then ! only master process does it
            print*, 'E_cut > E_CB, which is impossible in nonequilibrium simulation,', &
            ' resetting it to E_cut=', Scell%E_top-Scell%E_bottom
         endif
      endif
      numpar%E_cut = min(numpar%E_cut, Scell%E_top-Scell%E_bottom) ! [eV]

      ! Also ionization potential and the cross-section may need to be adjusted:
      if (matter%Atoms(1)%Ip(Nshl) > (Scell%E_top-Scell%E_bottom) ) then
         if (numpar%verbose) then
            if (numpar%MPI_param%process_rank == 0) then ! only master process does it
               print*, 'Ionization potential is smaller than CB width,', ' it will be reset to Ip=', Scell%E_top-Scell%E_bottom
               print*, 'Note that it will affect the cross-section'
            endif
         endif
         Egap_old = matter%Atoms(1)%Ip(Nshl)
         matter%Atoms(1)%Ip(Nshl) = Scell%E_top-Scell%E_bottom ! change it to the width of CB
         N_grid = size( matter%Atoms(1)%El_MFP(Nshl)%E )
         do k = 1, N_grid  ! recalculate the cross-section of scattering on CB (UNFINISHED, not combined with changes due to temperature!)
            ! Total recalculation (very time-consuming!):
            !Ele = matter%Atoms(1)%El_MFP(Nshl)%E(k) ! [eV] energy
            !call TotIMFP(Ele, matter, Scell%TeeV, 1, Nshl, L)  ! module "MC_cross_sections"
            !matter%Atoms(1)%El_MFP(Nshl)%L(k) = L ! [A] MFP

            ! Alternatively, use a simple shift of the energy scale:
            !call update_ionization_CS(matter, Nshl, (matter%Atoms(1)%Ip(Nshl) - Egap_old) )  ! below

            if (numpar%verbose) call print_progress('Progress:', k, N_grid)    ! module "Little_subroutines"
         enddo
      endif
   endselect

   ! In case we have a very wide gap material, cut-off cannot be smaller than the gap:
   if (numpar%E_cut < Scell%E_gap) then
      if (numpar%MPI_param%process_rank == 0) then ! only master process does it
         if (numpar%verbose) print*, 'Potential problem: E_cut < E_gap, resetting it'
      endif
      numpar%E_cut = Scell%E_gap
   endif

   select case (matter%Atoms(1)%TOCS(size(matter%Atoms(1)%TOCS))) ! which inelastic cross section to use (BEB vs CDF):
      case (1) ! CDF cross section
         matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)) = Scell%E_gap ! [eV] ionization potential of the valence band
      case default  ! BEB:
         ! Renormalization is optional:
         if (numpar%E_cut < matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)) ) then
            if (numpar%MPI_param%process_rank == 0) then ! only master process does it
               if (numpar%verbose) print*, 'Electron cut-off energy cannot be smaller than the highest ionization otential,', &
               'resetting it to E_cut=', matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip))
            endif
            numpar%E_cut = matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip))
         endif
   end select
   
!    print*, 'numpar%E_cut =', numpar%E_cut , Scell%E_gap , matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip))
!    print*, 'Ne:', Scell%Ne_low, Scell%Ne, Scell%Na
!    print*, 'E :', Scell%E_gap, Scell%E_VB_top, Scell%E_bottom
end subroutine find_band_gap



subroutine update_ionization_CS(matter, Nshl, dE, numpar)
   type(solid), intent(inout) :: matter ! materil parameters
   integer, intent(in) :: Nshl   ! VB index
   real(8), intent(in) :: dE  ! energy shift
   type(Numerics_param), intent(inout) :: numpar ! numerical parameters, including MC energy cut-off
   !----------------------------------------
   real(8) :: IMFP
   integer :: i, j, iE, N_sh

   matter%Atoms(1)%El_MFP(Nshl)%E(:) = matter%Atoms(1)%El_MFP(Nshl)%E(:) + dE

   ! And total cross-section:
   do i = 1, size(matter%Atoms) ! for all atoms
      N_sh = size(matter%Atoms(i)%Ip)
      do j = 1, N_sh ! for all shells of this atom
         if ((i == 1) .and. (j==N_sh)) then ! VB
            ! skip VB, to be added separately
         else ! core shells
            matter%El_MFP_tot%L(:) = matter%El_MFP_tot%L(:) + matter%Atoms(i)%El_MFP(j)%L(:) ! [1/A] inverse MFP
         endif
      enddo
   enddo
   ! Now, add the shifted VB:
   do iE = 1, size(matter%Atoms(1)%El_MFP(Nshl)%E)
      call Mean_free_path(matter%El_MFP_tot%E(iE), matter%Atoms(1)%El_MFP(Nshl), IMFP, inversed=.true.) ! [1/A], module "MC_cross_section"
      matter%El_MFP_tot%L(:) = matter%El_MFP_tot%L(:) + IMFP ! [1/A] inverse MFP
   enddo
end subroutine update_ionization_CS


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
         case (1) ! Enforced energy conservation (Etot = Ee + Eat = const; OUTDATED, NOT NEEDED):
            if (t .GT. numpar%t_Te_Ee) then ! Total energy is fixed:
               !call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! (SLOW) below
               call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! (FAST) below
            else ! electron temperature is fixed:
               call Electron_Fixed_Te(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! below
            endif
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]
            call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution
            Scell(NSC)%fe_eq = Scell(NSC)%fe ! instanteneous thermalization means both functions are the same

         case (2) ! Fixed temperature (Te=const; EXTERNALLY ENFORCED):
!             if (numpar%scc) then ! SCC, so the total energy is defined by the part H_0 without charge energy:
!                call Electron_Fixed_Te(Scell(NSC)%Ei_scc_part, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! below
!             else
               call Electron_Fixed_Te(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%mu, Scell(NSC)%TeeV) ! below
!             endif
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]
            call set_initial_fe(Scell, matter, Err) ! recalculate new electron distribution
            Scell(NSC)%fe_eq = Scell(NSC)%fe ! instanteneous thermalization means both functions are the same

         case (3) ! Born-Oppenheimer (constant populations):
            ! Do nothing with fe!
            ! Only get the kinetic temperature of electrons (out-of-equilibrium):
            call Electron_Fixed_Etot(Scell(NSC)%Ei, Scell(NSC)%Ne_low, Scell(NSC)%nrg%El_low, &
                                          Scell(NSC)%mu, Scell(NSC)%TeeV, .true.) ! below (FAST)
            Scell(NSC)%Te = Scell(NSC)%TeeV*g_kb ! save also in [K]

            ! Construct Fermi function with the given transient parameters (equivalent Te and mu):
            i_fe = size(Scell(NSC)%fe)   ! number of grid points in distribution function
            if (.not.allocated(Scell(NSC)%fe_eq)) allocate(Scell(NSC)%fe_eq(i_fe))
            call set_Fermi(Scell(NSC)%Ei, Scell(NSC)%TeeV, Scell(NSC)%mu, Scell(NSC)%fe_eq)   ! below

         case (4:5) ! Relaxation-time approximation; electron-electron collision integral
            ! Relaxing electrons via rate equation with given characteristic time:
            !call Do_relaxation_time(Scell(NSC), numpar)  ! below
            ! We only update it once per simulation step, not every time this subroutine called!

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
   case (5)    ! electron-electron collision integral
      if (present(skip_thermalization)) then
         call Do_e_e_collision(Scell(1), numpar, skip_thermalization)  ! below
      else
         call Do_e_e_collision(Scell(1), numpar)  ! below
      endif

      ! For output, get the partial equivalent Fermi distributions of CB and VB:
      call construct_CB_and_VB_equivalent_Fermi(Scell(1), numpar)    ! below
   endselect
end subroutine Electron_thermalization


! Electron-electron collision integral:
subroutine Do_e_e_collision(Scell, numpar, skip_thermalization)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   logical, intent(in), optional :: skip_thermalization
   !----------------------
   logical :: skip_step
   integer :: i_fe
   !real(8), dimension(size(Scell%fe),size(Scell%fe)) :: M_ee ! electron-electron scattering matrix elements
   real(8), dimension(:,:), allocatable :: M_ee ! electron-electron scattering matrix elements

   if (present(skip_thermalization)) then
      skip_step = skip_thermalization
   else
      skip_step = .false.
   endif

   allocate(M_ee(size(Scell%fe),size(Scell%fe)))

   ! Get the equivalent (kinetic) temperature and chemical potential:
   call Electron_Fixed_Etot(Scell%Ei, Scell%Ne_low, Scell%nrg%El_low, Scell%mu, Scell%TeeV, .true.) ! below (FAST)
   Scell%Te = Scell%TeeV*g_kb ! save also in [K]

   ! Construct Fermi function with the given transient parameters:
   i_fe = size(Scell%fe)   ! number of grid points in distribution function
   if (.not.allocated(Scell%fe_eq)) allocate(Scell%fe_eq(i_fe), source = 0.0d0)
   call set_Fermi(Scell%Ei, Scell%TeeV, Scell%mu, Scell%fe_eq)   ! below

   M_ee = 0.1d0  ! just to test!!

   ! Solve Boltzmann collision integral:
   if (.not.skip_step) then ! do the e-e collisions step:
      ! If we don't use substeps:s
      !call Boltzmann_e_e_IN(Scell, numpar, Scell%Ei, Scell%fe, numpar%dt) ! below

      ! If we want to use substeps:
      call Boltzmann_e_e_IN(Scell, numpar, Scell%Ei, Scell%fe, numpar%dt, Npoints=int(numpar%tau_fe)) ! below
   endif ! (.not.skip_step)


   !--------------------------
   ! Extra check for smoothening unphysical artefacts that may be present after MC:
   call smoothening_step(Scell, numpar, eps_precision=1.0d-3) ! below

   ! clean up:
   deallocate(M_ee)
end subroutine Do_e_e_collision



! Electron-electron collision integral (UNFINISHED DUE TO PROBLEMS WITH ENERGY AND PARTICLE NUMBER CONSERVATION):
subroutine Boltzmann_e_e_IN(Scell, numpar, Ev, fe, dt, Npoints) ! calculates change of distribution function via Boltzmann collision integral
! See examples of Boltzmann equation in energy space e.g. in
! [B. Rethfeld, A. Kaiser, M. Vicanek and G. Simon, Phys.Rev.B 65, 214303, (2002)]
! (although we use here only electron-electron integral, and with very different matrix element)
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   real(8), dimension(:), intent(inout) :: fe  ! electron distribution function (ocupation numbers of the energy levels Ei)
   real(8), intent(in) :: dt ! [fs] time-step
   integer, intent(in), optional :: Npoints     ! number of additional time-points (substeps)
   !-----------------------------------------
   real(8) :: Ei, Ef, Ei2, Ef2 ! energies of initial and final states for #1 and #2 electron
   !real(8), dimension(size(fe)) :: fe_temp, f_cur, f_cur0 ! temporary distribution to work with
   !real(8), dimension(size(fe),size(fe)) :: alpha_ij, beta_ij     ! coefficients in Boltzmann collision integral
   !real(8), dimension(size(fe),size(fe)) :: alpha_ij0, beta_ij0     ! coefficients in Boltzmann collision integral
   !real(8), dimension(size(fe),size(fe)) :: M_ee ! electron-electron scattering matrix elements
   real(8), dimension(:), allocatable :: fe_temp, f_cur, f_cur0 ! temporary distribution to work with
   real(8), dimension(:,:), allocatable :: alpha_ij, beta_ij     ! coefficients in Boltzmann collision integral
   real(8), dimension(:,:), allocatable :: alpha_ij0, beta_ij0     ! coefficients in Boltzmann collision integral
   real(8), dimension(:,:), allocatable :: M_ee ! electron-electron scattering matrix elements
   real(8) :: dt_small, split_coef
   real(8) :: fe_f ! average final state
   real(8) :: E_tot, E_tot2 ! total energy used for testing [eV]
   real(8) :: N_tot, N_tot2
   integer :: i, j, i2, j2, N, i_cur, i_f1, i_f2, N_steps, i_step, N_steps_default, i_PC, N_PC
   logical :: within_range ! check if final energy level is within possible range

   N = size(fe)   ! total number of energy levels available

   allocate(fe_temp(N))
   allocate(f_cur(N))
   allocate(f_cur0(N))
   allocate(alpha_ij(N,N))
   allocate(beta_ij(N,N))
   allocate(alpha_ij0(N,N))
   allocate(beta_ij0(N,N))

   fe_temp = fe   ! just to start
   f_cur(:) = fe_temp(:)   ! to start with
   f_cur0(:) = fe_temp(:)   ! to start with
   alpha_ij = 0.0d0     ! to start with
   beta_ij = 0.0d0      ! to start with
   alpha_ij0 = 0.0d0    ! to start with
   beta_ij0 = 0.0d0     ! to start with


   ! Test the energy conservation:
   call set_total_el_energy(Ev,fe_temp,E_tot)


   ! Precalculate the overlap-matrix coefficients:
   call get_electron_electron_overlap(Scell, M_ee)    ! below


#ifdef MPI_USED   ! use the MPI version
   ! NOT DONE YET

#else ! use OpenMP instead
   ! Do with a number of smaller time-steps:
   N_steps_default = 1
   if (present(Npoints)) then
      if (Npoints > 0) then
         N_steps = Npoints
      else
         N_steps = N_steps_default   ! number of time-substeps
      endif
   else
      N_steps = N_steps_default   ! number of time-substeps
   endif

   ! Define the number of Predictor-Correction iterations:
   N_PC = 1      ! max (N=1 means no predictor-correction used)
   split_coef = 0.3d0   ! splitting coefficient in PC

   dt_small = dt / dble(N_steps)    ! [fs] size of the substep
   do i_step = 1, N_steps    ! time substeps
      do i_PC  = 1, N_PC     ! Predictor-corrector

      ! Get the upper triangle of coefficients alpha and beta:
!$omp PARALLEL shared(alpha_ij, beta_ij, f_cur, fe_temp, alpha_ij0, beta_ij0)
!$omp do schedule(dynamic) private(i)
         do i = 1, N
            ! Without PC:
            !call get_Boltzmann_alpha_beta(i, Ev, M_ee, f_cur, dt_small, alpha_ij, beta_ij)  ! module "Electron_electron_scattering"
            ! With PC:
            call get_Boltzmann_alpha_beta(Scell, i, Ev, M_ee, f_cur, dt_small, alpha_ij0, beta_ij0)  ! module "Electron_electron_scattering"
         enddo ! i
!$omp end do
!$omp barrier

      ! Get the lower triangle of coefficients alpha and beta:
!$omp do schedule(dynamic) private(i,j)
         do i = 2, N
            do j = 1, i-1 ! all levels from where second electron can scatter off
               ! using symmetry coef(i,j) = coef(j,i)
               ! Without PC:
               !alpha_ij(i,j) = alpha_ij(j,i)
               !beta_ij(i,j) = beta_ij(j,i)
               ! With PC:
               alpha_ij0(i,j) = alpha_ij0(j,i)
               beta_ij0(i,j) = beta_ij0(j,i)
            enddo ! j = 1, i-1
         enddo ! i
!$omp end do
!$omp barrier


         ! Predictor corrector algorithm:
         if (i_PC == 1) then
            alpha_ij = alpha_ij0
            beta_ij = beta_ij0
         else
            ! Upper triangle:
!$omp do schedule(dynamic) private(i,j)
            do i = 1, N
               do j = i, N ! all levels from where second electron can scatter off
                  alpha_ij(i,j) = split_coef * alpha_ij0(i,j) + (1.0d0-split_coef)*alpha_ij(i,j)
                  beta_ij(i,j)  = split_coef * beta_ij0(i,j)  + (1.0d0-split_coef)*beta_ij(i,j)

                  !if (j==i) print*, i,j, alpha_ij(i,j), beta_ij(i,j)

               enddo ! j = 1, i-1
            enddo ! i
!$omp end do
!$omp barrier
            ! Lower triangle:
!$omp do schedule(dynamic) private(i,j)
            do i = 2, N
               do j = 1, i-1 ! all levels from where second electron can scatter off
                  alpha_ij(i,j) = alpha_ij(j,i)
                  beta_ij(i,j) = beta_ij(j,i)
                  ! Test (wrong)
                  !alpha_ij(i,j) = beta_ij(j,i)
                  !beta_ij(i,j) = alpha_ij(j,i)
               enddo ! j
            enddo ! i
!$omp end do
!$omp barrier
         endif

!$omp do schedule(dynamic) private(i)
         do i = 1, N ! all energy levels on which we looking for changes (left-hand side of the equation)
            ! Without substeps:
            !call Boltzmann_solution(i, fe_temp, fe, dt, 0) ! module "Electron_electron_scattering"
            ! With substeps:
            call Boltzmann_solution(i, f_cur, fe_temp, dt_small, alpha_ij, beta_ij, 0) ! module "Electron_electron_scattering"
         enddo ! i
!$omp end do
!$omp end parallel

         ! PC update of the distribution:
         f_cur(:) = split_coef*f_cur(:) + (1.0d0-split_coef)*f_cur0(:)

         ! Explicit solution may make f<0 or f>2, correct such cases (shouldn't happen, but just in case):
         call patch_distribution(f_cur, Ev, Scell, numpar)    ! below

         ! Correction to ensure long-term particle and energy conservation:
         call correction_distribution(fe, Ev, f_cur)   ! below

         ! Save for the next step
         f_cur0(:) = f_cur(:)

         ! Tests:
         !N_tot = SUM(fe)
         !N_tot2 = SUM(f_cur)
         !write(*,'(a, i0, f, f, f, a)') 'Boltzmann_e_e N1:', i_PC, N_tot, N_tot2, ABS(N_tot2 - N_tot)/N_tot*100.0d0, '%'
         ! Test the energy conservation:
         !call set_total_el_energy(Ev, f_cur, E_tot2)      ! below
         !write(*,'(a, i0, f, f, f, a, f)') 'Boltzmann_e_e E1:', i_PC, E_tot, E_tot2, ABS(E_tot2 - E_tot)/ABS(E_tot)*100.0d0, '%', &
         !       (E_tot2 - E_tot)/(N_tot2 - N_tot)

      enddo ! i_PC

      ! Time-propagation:
      fe_temp(:) = f_cur(:)   ! update for the next substep

   enddo ! i_step
#endif

   ! Test the particle and energy conservation:
   !call set_total_el_energy(Ev, fe_temp, E_tot2)      ! below
   !N_tot = SUM(fe)
   !N_tot2 = SUM(fe_temp)
   !print*, 'Boltzmann_e_e N1:', N_tot, N_tot2, ABS(N_tot2-N_tot)/N_tot*100.0d0
   !print*, 'Boltzmann_e_e E1:', E_tot, E_tot2, ABS(E_tot2 - E_tot)/ABS(E_tot)*100.0d0

   fe = fe_temp ! output: updated distribution function

   ! clean up:
   deallocate(fe_temp, f_cur, f_cur0, alpha_ij, beta_ij, alpha_ij0, beta_ij0)
end subroutine Boltzmann_e_e_IN


subroutine get_electron_electron_overlap(Scell, M_ee)    ! below
   type(Super_cell), intent(in) :: Scell ! supercell with all the atoms as one object
   real(8), dimension(:,:), intent(inout) :: M_ee ! matrix element for electron-ion coupling
   !=====================================

   ! Calculate nonadiabatic-coupling matrix element:
   ASSOCIATE (ARRAY => Scell%TB_Hamil(:,:)) ! this is the sintax we have to use to check the class of defined types
      ! Different expressions for orthogonal and non-orthogonal bases:
      select type(ARRAY)
      type is (TB_H_Pettifor) ! TB parametrization according to Pettifor: orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee)      ! module "Electron_electron_scattering"

      type is (TB_H_Molteni)  ! TB parametrization accroding to Molteni: orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee)      ! module "Electron_electron_scattering"

      type is (TB_H_Fu)  ! TB parametrization accroding to Fu: orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee)      ! module "Electron_electron_scattering"

      type is (TB_H_NRL)  ! TB parametrization accroding to NRL method: non-orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee, Scell%Sij)      ! module "Electron_electron_scattering"

      type is (TB_H_DFTB)  ! TB parametrization accroding to DFTB: non-orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee, Scell%Sij)      ! module "Electron_electron_scattering"

      type is (TB_H_3TB)  ! TB parametrization accroding to 3TB: non-orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee, Scell%Sij)      ! module "Electron_electron_scattering"

      type is (TB_H_xTB)  ! TB parametrization accroding to xTB: non-orthogonal
         call Electron_electron_scattering_Kij(Scell, Scell%Ha, M_ee, Scell%Sij)      ! module "Electron_electron_scattering"

      end select
   END ASSOCIATE
end subroutine get_electron_electron_overlap




! Two-levels corection to the electron-distribution function to ensure particle and energy conservation:
subroutine correction_distribution(fe, Ei, f_cur)
   real(8), dimension(:),intent(in) :: fe ! original distribution (before e-e scattering)
   real(8), dimension(:),intent(in) :: Ei ! energy levels
   real(8), dimension(:),intent(inout) :: f_cur ! distribution to be corrected (after e-e scattering)
   !-------------------------------------
   integer :: i, j, Nsiz, i1, i2, Npart, coun, i_low, i_high, cou1, cou2
   real(8) :: E_in, E_fin, N_in, N_fin, dE, dN, RN, N_cur, N_rand, df1, df2, temp, df1_min, df2_min
   !logical, dimension(size(fe)) :: fe_mask
   !real(8), dimension(size(fe)) :: f_temp
   logical, dimension(:), allocatable :: fe_mask
   real(8), dimension(:), allocatable :: f_temp
   logical :: found_yet

   Nsiz = size(fe)

   allocate(fe_mask(Nsiz))
   allocate(f_temp(Nsiz))

   fe_mask = .false.    ! to start with
   f_temp = f_cur       ! to start with

   ! Define initial number of particles and energy:
   N_in = sum(fe) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, fe, E_in)      ! below

   ! Define final number of particles and energy:
   N_fin = sum(f_cur) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, f_cur, E_fin)      ! below

   ! Number of particles and Energy difference:
   dN = N_in - N_fin
   dE = E_in - E_fin

   temp = dN/N_in
   if (abs(temp) < 1.0d-12) return ! no correction meeded

   ! Go through all the pairs of levels, finding the minimal deviation:
   df1_min = 1.0d10     ! to start with
   df2_min = 1.0d10     ! to start with
   found_yet = .false.  ! to start with
   do i = 1, Nsiz !
      do j = i, Nsiz
         if (j == i) cycle     ! skip the same level

         ! check all the levels, until we find the one that works:
         df1 = (dN*Ei(j) - dE) / (Ei(j) - Ei(i))
         df2 = (dN*Ei(i) - dE) / (Ei(i) - Ei(j))

         ! Correct the distribution:
         f_temp(i) = f_cur(i) + df1
         f_temp(j) = f_cur(j) + df2

         if ( (f_temp(i) >= 0.0d0) .and. (f_temp(i) <= 2.0d0) .and. &
              (f_temp(j) >= 0.0d0) .and. (f_temp(j) <= 2.0d0) ) then ! it is appropriate
            ! Save these energy levels, if needed:
            if ( (abs(df1) < abs(df1_min)) .and. (abs(df2) < abs(df2_min)) ) then
               found_yet = .true.      ! found at least something...
               i_low = i   ! save for later
               i_high = j  ! save for later
               df1_min = df1
               df2_min = df2
            endif
         endif
      enddo ! j
   enddo ! j

   ! Update distribution:
   if (found_yet) then
      f_cur(i_low) = f_cur(i_low) + df1_min
      f_cur(i_high) = f_cur(i_high) + df2_min
   else
      !print*, 'Found nothing...'
   endif

   ! Test:
   ! Define final number of particles and energy:
   N_fin = sum(f_cur) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, f_cur, E_fin)      ! below
   ! Test the energy conservation:
   dE = -(E_fin - E_in)
   dN = -(N_fin - N_in)
   !write(*,'(a, f, f, f, a)') 'Boltzmann_e_e N2:', N_in, N_fin, ABS(dN)/N_in*100.0d0, '%'
   !rite(*,'(a, f, f, f, a)') 'Boltzmann_e_e E2:', E_in, E_fin, ABS(dE)/ABS(E_in)*100.0d0, '%'

   ! clean up:
   deallocate(fe_mask, f_temp)
end subroutine correction_distribution


! Two-levels corection to the electron-distribution function to ensure particle and energy conservation:
subroutine correction_distribution_first(fe, Ei, f_cur)
   real(8), dimension(:),intent(in) :: fe ! original distribution (before e-e scattering)
   real(8), dimension(:),intent(in) :: Ei ! energy levels
   real(8), dimension(:),intent(inout) :: f_cur ! distribution to be corrected (after e-e scattering)
   !-------------------------------------
   integer :: i, j, Nsiz, i1, i2, Npart, coun, i_low, i_high, cou1, cou2
   real(8) :: E_in, E_fin, N_in, N_fin, dE, dN, RN, N_cur, N_rand, df1, df2, temp
   !logical, dimension(size(fe)) :: fe_mask
   !real(8), dimension(size(fe)) :: f_temp
   logical, dimension(:), allocatable :: fe_mask
   real(8), dimension(:), allocatable :: f_temp
   logical :: found_yet

   Nsiz = size(fe)

   allocate(fe_mask(Nsiz))
   allocate(f_temp(Nsiz))

   fe_mask = .false.    ! to start with
   f_temp = f_cur       ! to start with

   ! Define initial number of particles and energy:
   N_in = sum(fe) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, fe, E_in)      ! below

   ! Define final number of particles and energy:
   N_fin = sum(f_cur) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, f_cur, E_fin)      ! below

   ! Number of particles and Energy difference:
   dN = N_in - N_fin
   dE = E_in - E_fin

   temp = dN/N_in
   if (abs(temp) < 1.0d-12) return ! no correction meeded

   ! Chose levels we can use for correction:
   do i = 1, Nsiz
      if ( (f_cur(i) > 0.1d0) .and. (f_cur(i) < 1.9d0) ) then ! let's use those levels
         fe_mask(i) =.true.
      endif
   enddo
   Npart = COUNT(fe_mask)

   !print*, 'n', dN, dE, Npart

   if (Npart < 2) return ! no correction possible

   ! Select random levels from there:
   ! Level 1:
   call random_number(RN)
   N_rand = RN * dble(Npart)
   ! Find the randomly selected level:
   cou1 = 0 ! to start with
   do i = 1, Nsiz
      if (fe_mask(i)) then ! only for the allowed levels:
         cou1 = cou1 + 1
         i1 = i
         if (cou1 >= int(N_rand)) exit ! found the level
      endif
   enddo
   i_low = i1 ! save for later

   ! Level 2:
   do i2 = 1, Nsiz
      if (i2 == i1) cycle     ! skip the same level
      i_high = i2 ! save for later
      ! check all the levels, until we find the one that works:
      df1 = (dN*Ei(i2) - dE) / (Ei(i2) - Ei(i1))
      df2 = (dN*Ei(i1) - dE) / (Ei(i1) - Ei(i2))

      ! Correct the distribution:
      f_temp(i1) = f_cur(i1) + df1
      f_temp(i2) = f_cur(i2) + df2

      if ( (f_temp(i1) >= 0.0d0) .and. (f_temp(i1) <= 2.0d0) .and. &
           (f_temp(i2) >= 0.0d0) .and. (f_temp(i2) <= 2.0d0) ) then ! it is appropriate
         exit     ! done, no need to continue
      endif
   enddo
   ! Update distribution:
   f_cur(i_low) = f_cur(i_low) + df1
   f_cur(i_high) = f_cur(i_high) + df2

   ! Test:
   ! Define final number of particles and energy:
   N_fin = sum(f_cur) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, f_cur, E_fin)      ! below
   ! Test the energy conservation:
   dE = -(E_fin - E_in)
   dN = -(N_fin - N_in)
   !write(*,'(a, f, f, f, a)') 'Boltzmann_e_e N2:', N_in, N_fin, ABS(dN)/N_in*100.0d0, '%'
   !rite(*,'(a, f, f, f, a)') 'Boltzmann_e_e E2:', E_in, E_fin, ABS(dE)/ABS(E_in)*100.0d0, '%'


   ! clean up:
   deallocate(fe_mask, f_temp)
end subroutine correction_distribution_first


! Linear function corection to the electron-distribution function to ensure particle and energy conservation:
! ERROR: THIS VERSION OF CORRECTION PRODUCES f>2 or f<0
subroutine lilnear_correction_distribution(fe, Ei, f_cur)
   real(8), dimension(:),intent(in) :: fe ! original distribution (before e-e scattering)
   real(8), dimension(:),intent(in) :: Ei ! energy levels
   real(8), dimension(:),intent(inout) :: f_cur ! distribution to be corrected (after e-e scattering)
   !-------------------------------------
   integer :: i, Nsiz, Npart
   real(8) :: E_in, E_fin, N_in, N_fin, dE, dN, E2, Epart, temp
   real(8) :: Alpha, Beta     ! linear coefficients for distribution rescaling
   !logical, dimension(size(fe)) :: fe_mask
   !real(8), dimension(size(fe)) :: df ! distribution to be corrected (after e-e scattering)
   logical, dimension(:), allocatable :: fe_mask
   real(8), dimension(:), allocatable :: df ! distribution to be corrected (after e-e scattering)

   Nsiz = size(fe)

   allocate(fe_mask(Nsiz))
   allocate(df(Nsiz))

   fe_mask = .false.    ! to start with

   ! Define initial number of particles and energy:
   N_in = sum(fe) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, fe, E_in)      ! below

   ! Define final number of particles and energy:
   N_fin = sum(f_cur) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, f_cur, E_fin)      ! below

   ! Number of particles and Energy difference:
   dN = N_in - N_fin
   dE = E_in - E_fin

   ! Tests:
   !write(*,'(a, f, f, f, a)') 'Boltzmann_e_e N1:', N_in, N_fin, ABS(dN)/N_in*100.0d0, '%'
   ! Test the energy conservation:
   !write(*,'(a, f, f, f, a)') 'Boltzmann_e_e E1:', E_in, E_fin, ABS(dE)/ABS(E_in)*100.0d0, '%'


   ! Make a correction with a linear function, redistributing the excessive (or missing) particles and energy:
   ! If there is excessive number of particles, remove them from where there are plenty:
   ! 1) Define the levels to be used to correct:

   ! This version of the mask produces ERRORS: f>2 or f<0
   if (dN < 0.0d0) then
      do i = 1, Nsiz
         if (fe(i) <= 1.0d0) then ! half-empty
            fe_mask(i) = .true.     ! mark pessimistic levels
         endif
      enddo
   else ! if there are missing particles, add them where there are few:
      do i = 1, Nsiz
         if (fe(i) >= 1.0d0) then ! half-full
            fe_mask(i) = .true.     ! mark optimistic levels
         endif
      enddo
   endif ! (dN > 0.0d0)

   ! 2) Calculate properties to be used in correction:
   Npart = COUNT(fe_mask)
   if (Npart < 1) return ! no correction possible

   Epart = sum(Ei(:), MASK = fe_mask)
   E2 = sum(Ei(:)**2, MASK = fe_mask)
   temp = E2 - Epart**2/dble(Npart)
   if (abs(temp) < 1.0d-12) return ! no correction possible

   Alpha = ( dE - dN/dble(Npart)*Epart ) / temp
   Beta = (dN - Epart*Alpha)/dble(Npart)

   ! 3) Define the correction to the distribution:
   df = 0.0d0     ! to start with
   do i = 1, Nsiz
      if (fe_mask(i)) then ! only for those levels:
         df(i) = Alpha * Ei(i) + Beta
      endif
      print*, i, df(i), f_cur(i), f_cur(i)+df(i)
   enddo

   ! 4) Correct the distribution:
   f_cur(:) = f_cur(:) + df(:)


   ! Test:
   ! Define final number of particles and energy:
   N_fin = sum(f_cur) ! total number of particles
   ! and total energy:
   call set_total_el_energy(Ei, f_cur, E_fin)      ! below
   ! Test the energy conservation:
   dE = -(E_fin - E_in)
   dN = -(N_fin - N_in)
   !write(*,'(a, f, f, f, a)') 'Boltzmann_e_e N2:', N_in, N_fin, ABS(dN)/N_in*100.0d0, '%'
   !write(*,'(a, f, f, f, a)') 'Boltzmann_e_e E2:', E_in, E_fin, ABS(dE)/ABS(E_in)*100.0d0, '%'

   !pause 'correction_distribution'

   ! clean up:
   deallocate(fe_mask, df)
end subroutine lilnear_correction_distribution



! Correct distribution for f>2 or f<0:
subroutine patch_distribution(fe, Ei, Scell, numpar)
   real(8), dimension(:),intent(in) :: Ei ! energy levels
   real(8), dimension(:),intent(inout) :: fe ! distribution
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !-----------------------------
   integer :: i, N_siz, i_below, i_above, counter, j
   real(8) :: eps, df, dE, Ee
   logical :: trouble_present, trouble_in_this_point

   eps = 1.0d-10     ! how close energy levels are allowed to be
   N_siz = size(fe)  ! number of energy levels
   trouble_present = .true.   ! just to start
   counter = 0 ! to start with

   TP:do while (trouble_present) ! do until trouble is solved
      trouble_present = .false.   ! assume no problem remains
      counter = counter + 1   ! count iterations
      ! Now, check if there is a problem:
      do i = 1, N_siz ! check that there is no problem in distribution function change
         trouble_in_this_point = .false.   ! assume no problem at this point "i"
         if (fe(i) > 2.0d0+eps) then
            trouble_present = .true.   ! there is an unphysical value, correct it and check again
            trouble_in_this_point = .true.   ! there is an unphysical value, correct it and check again
            df = fe(i) - 2.0d0   ! excessive part to be removed
            !fe(i) = 2.0d0        ! distribution adjusted to accceptable
         elseif (fe(i) > 2.0d0) then   ! it's within [2; 2+eps]
            fe(i) = 2.0d0        ! distribution adjusted to accceptable
         elseif (fe(i) < 0.0d0-eps) then
            trouble_present = .true.   ! there is an unphysical value, correct it and check again
            trouble_in_this_point = .true.   ! there is an unphysical value, correct it and check again
            df = fe(i)      ! missing part to be added
            !fe(i) = 0.0d0   ! distribution adjusted to accceptable
         elseif (fe(i) < 0.0d0) then  ! it's within [0-eps;0]
            fe(i) = 0.0d0        ! distribution adjusted to accceptable
         endif

         if (trouble_in_this_point) then ! try to solve it:
            !print*, 'Error in patch_distribution:', counter, 'i=', i, fe(i)+df, Ei(i)
            !if (i>1)print*, fe(i-1), Ei(i-1)
            !if (i<size(fe)) print*, fe(i+1), Ei(i+1)

            if (i >= N_siz .or. (i <= 1)) print*, 'Potential problem in patch_distribution #1:', i, fe(i), Ei(i), N_siz
            if (counter > 2000) then
               print*, 'Too many iterations patch_distribution:', counter, i, fe(i), Ei(i), N_siz
               ! Just try to cut out the problematic points:
               do j = 1, N_siz
                  if (fe(j) > 2.0d0) then
                     fe(j) = 2.0d0
                  elseif (fe(j) < 0.0d0) then
                     fe(j) = 0.0d0
                  endif
               enddo
               ! And we are done here:
               exit TP
            endif

            Ee = Ei(i)  ! this energy level contains problematic distribution point

            ! Find the levels to redistribute electrons to:
            call choose_level(Ei, fe, df, i, i_below, i_above)   ! below

            if ( (i_above > N_siz .or. (i_above < 1)) .or. (i_below > N_siz .or. (i_below < 1)) ) then
               !if (i_above > N_siz .or. (i_above < 1)) print*, 'Potential problem in patch_distribution #2a:', i, i_above, fe(i), Ei(i)
               !if (i_below > N_siz .or. (i_below < 1)) print*, 'Potential problem in patch_distribution #2b:', i, i_below, fe(i), Ei(i)
               trouble_present = .true.
               exit TP
            else
               ! Fractions of electron distributed between two levels,
               ! ensuring conservation of particles and energy:
               dE = Ei(i_above)-Ei(i_below)   ! energy levels difference
               fe(i_below) = fe(i_below) + (Ei(i_above) - Ee)/dE * df
               fe(i_above) = fe(i_above) + (Ee - Ei(i_below))/dE * df
               fe(i) = fe(i) - df   ! distribution changed
            endif

            !print*, 'i=', i, i_above, i_below, df, fe(i_above), fe(i_above) - (Ee - Ei(i_below))/dE * df , fe(i_below), fe(i_below) - (Ei(i_above) - Ee)/dE * df

         endif
      enddo ! i = 1, N_siz
   enddo TP ! trouble_present

   ! Check if there was a situation that electrons could not be redistributed:
   if (trouble_present) then  ! Do thermalization instead, adjust all levels:
      call Do_relaxation_time(Scell, numpar, skip_partial=.true.) ! module "Electron_tools"
      trouble_present = .false.
      ! And the final check:
      do i = 1, N_siz ! check that there is no problem in distribution function change
         if (fe(i) > 2.0d0+eps) then   ! it's outside of [2; 2+eps]
            print*, 'Problem in patch_distribution #3a:', i, fe(i)
            fe(i) = 2.0d0
            trouble_present = .true.   ! there still is a problem
         elseif (fe(i) > 2.0d0) then   ! it's within [2; 2+eps]
            fe(i) = 2.0d0        ! distribution adjusted to accceptable
         elseif (fe(i) < -eps) then  ! it's outside of [0-eps;0]
            print*, 'Problem in patch_distribution #3b:', i, fe(i)
            fe(i) = 0.0d0        ! distribution adjusted to accceptable
            trouble_present = .true.   ! there still is a problem
         elseif (fe(i) < 0.0d0) then  ! it's within [0-eps;0]
            fe(i) = 0.0d0        ! distribution adjusted to accceptable
         endif
      enddo

      if (trouble_present) then  ! Oh well, I am out of ideas...
         print*, 'Problem persists, we may have lost some particles and energy...'
      else  ! Yey, it worked!
         print*, 'Problem avoided, nothing to worry about.'
      endif
   endif
end subroutine patch_distribution



subroutine choose_level(Ei, fe, df, i, i_below, i_above)
   real(8), dimension(:),intent(in) :: Ei, fe ! distribution
   real(8), intent(in) :: df  ! change of fe
   integer, intent(in) :: i   ! the level given
   integer, intent(out) :: i_below, i_above  ! two levels chosen
   !--------------------
   integer :: j, k, N_siz
   real(8) :: dE, eps, fe1, fe2
   logical :: found_it

   eps = 1.0d-10
   N_siz = size(fe)
   found_it = .false.   ! to start with
   ! Check below the level i:
   FV:do j = i-1, 1, -1
      do k = j+1, N_siz
         if (k /= i) then
            dE = Ei(k) - Ei(j)   ! energy levels difference
            if (dE > eps) then
               fe1 = fe(j) + (Ei(k) - Ei(i))/dE * df
               fe2 = fe(k) + (Ei(i) - Ei(j))/dE * df
               ! Check if we found an acceptable pair of levels to redistribute electrons into:
               if ( ((fe1 >= 0.0d0) .and. (fe1 <= 2.0d0)) .and. ((fe2 >= 0.0d0) .and. (fe2 <= 2.0d0)) ) then
                  found_it = .true.
                  exit FV
               endif
            endif
         endif ! k/=i
      enddo ! k
   enddo FV ! j
   ! Check above the level i:
   if (.not. found_it) then
      FV2:do j = i+1, N_siz-1
         do k = j+1, N_siz
            if (k /= i) then
               dE = Ei(k) - Ei(j)   ! energy levels difference
               if (dE > eps) then
                  fe1 = fe(j) + (Ei(k) - Ei(i))/dE * df
                  fe2 = fe(k) + (Ei(i) - Ei(j))/dE * df
                  ! Check if we found an acceptable pair of levels to redistribute electrons into:
                  if ( ((fe1 >= 0.0d0) .and. (fe1 <= 2.0d0)) .and. ((fe2 >= 0.0d0) .and. (fe2 <= 2.0d0)) ) then
                     found_it = .true.
                     exit FV2
                  endif
               endif
            endif ! k/=i
         enddo ! k
      enddo FV2 ! j
   endif
   ! Get the output values:
   i_below = j
   i_above = k

   ! Check if there is a situation when it's impossible to find the levels:
!   if (.not. found_it) then
!      print*, 'Problem in choose_level:', i, i_below, i_above
!      print*, 'fe=', fe(i), df, 'Ee=', Ei(i)
!   endif

end subroutine choose_level



subroutine choose_level_OLD(Ei, fe, df, i, i_below, i_above)
   real(8), dimension(:),intent(in) :: Ei, fe ! distribution
   real(8), intent(in) :: df  ! change of fe
   integer, intent(in) :: i   ! the level given
   integer, intent(out) :: i_below, i_above  ! two levels chosen
   !--------------------
   integer :: j, k, N_siz
   real(8) :: dE, eps, fe1, fe2
   logical :: found_it

   eps = 1.0d-10
   N_siz = size(fe)
   found_it = .false.   ! to start with

   FV:do j = 1, N_siz-1
      do k = j, N_siz
         dE = Ei(k) - Ei(j)   ! energy levels difference
         if (dE > eps) then
            fe1 = fe(j) + (Ei(k) - Ei(i))/dE * df
            fe2 = fe(k) + (Ei(i) - Ei(j))/dE * df
            ! Check if we found an acceptable pair of levels to redistribute electrons into:
            if ( ((fe1 >= 0.0d0) .and. (fe1 <= 2.0d0)) .and. ((fe2 >= 0.0d0) .and. (fe2 <= 2.0d0)) ) then
               found_it = .true.
               exit FV
            endif
         endif
      enddo ! k
   enddo FV ! j
   ! Get the output values:
   i_below = j
   i_above = k
end subroutine choose_level_OLD



subroutine construct_CB_and_VB_equivalent_Fermi(Scell, numpar)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !----------------------------------
   real(8) :: Te_VB, Te_CB
   integer :: i_fe

   if (.not. numpar%do_partial_thermal) return  ! no need to do anything

   i_fe = size(Scell%fe)   ! number of grid points in distribution function

   ! Get the Fermi function for VB and CB separately:
   ! VB:
   ! Get the number of particles and energy in the band:
   Scell%Ne_low_VB = get_N_partial(Scell%fe, 1, Scell%N_Egap) ! below
   Scell%El_low_VB = get_E_partial(Scell%Ei, Scell%fe, 1, Scell%N_Egap) ! below

   ! Get the equivalent temperature and chem.potential of the band:
   call Electron_Fixed_Etot_partial(Scell%Ei, Scell%Ne_low_VB, Scell%El_low_VB, Scell%mu_VB, Te_VB, i_end_in=Scell%N_Egap) ! below
   Scell%Te_VB = Te_VB*g_kb ! save in [K]

   ! Construct equivalent Fermi distribution:
   call set_Fermi(Scell%Ei, Te_VB, Scell%mu_VB, Scell%fe_eq_VB, i_end=Scell%N_Egap)  ! below

   ! CB:
   ! Get the number of particles and energy in the band:
   Scell%Ne_low_CB = get_N_partial(Scell%fe, Scell%N_Egap+1, i_fe) ! below
   Scell%El_low_CB = get_E_partial(Scell%Ei, Scell%fe,  Scell%N_Egap+1, i_fe) ! below

   ! Get the equivalent temperature and chem.potential of the band:
   call Electron_Fixed_Etot_partial(Scell%Ei, Scell%Ne_low_CB, Scell%El_low_CB, Scell%mu_CB, Te_CB, i_start_in=Scell%N_Egap+1) ! below
   Scell%Te_CB = Te_CB*g_kb ! save in [K]

   ! Construct equivalent Fermi distribution:
   call set_Fermi(Scell%Ei, Te_CB, Scell%mu_CB, Scell%fe_eq_CB, i_start=Scell%N_Egap+1)  ! below
end subroutine construct_CB_and_VB_equivalent_Fermi




! Relaxation time approximation:
subroutine Do_relaxation_time(Scell, numpar, skip_thermalization, skip_partial, do_fast)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   logical, intent(in), optional :: skip_thermalization, skip_partial, do_fast
   !----------------------
   real(8) :: exp_dttau, extra_dt, extra_tau, eps, Te_CB, Te_VB, Ne, Ne_eq
   integer :: i_fe, i, i_cycle, N_cycle
   logical :: skip_step, extra_cycle, skip_part, do_fast_search

   if (present(skip_thermalization)) then
      skip_step = skip_thermalization
   else
      skip_step = .false.
   endif

   if (present(skip_partial)) then
      skip_part = skip_partial
   else
      skip_part = .false.
   endif

   if (present(do_fast)) then
      do_fast_search = do_fast
   else
      do_fast_search = .true.
   endif

   ! Get the equivalent (kinetic) temperature and chemical potential:
   call Electron_Fixed_Etot(Scell%Ei, Scell%Ne_low, Scell%nrg%El_low, Scell%mu, Scell%TeeV, do_fast_search) ! below (FAST)
   Scell%Te = Scell%TeeV*g_kb ! save also in [K]

   ! Construct Fermi function with the given transient parameters:
   i_fe = size(Scell%fe)   ! number of grid points in distribution function
   if (.not.allocated(Scell%fe_eq)) allocate(Scell%fe_eq(i_fe), source = 0.0d0)
   call set_Fermi(Scell%Ei, Scell%TeeV, Scell%mu, Scell%fe_eq)   ! below

   ! Solve rate equation:
   if (.not.skip_step) then ! do the thermalization step:
      !--------------------------
      ! I. Check partial-band thermalization (separate for VB and CB)
      BNDS:if (numpar%do_partial_thermal .and. (.not.skip_part)) then ! Partial thermalization is on:
         ! Get the Fermi function for VB and CB separately:
         ! VB:
         ! Get the number of particles and energy in the band:
         Scell%Ne_low_VB = get_N_partial(Scell%fe, 1, Scell%N_Egap) ! below
         Scell%El_low_VB = get_E_partial(Scell%Ei, Scell%fe, 1, Scell%N_Egap) ! below

         ! Get the equivalent temperature and chem.potential of the band:
         call Electron_Fixed_Etot_partial(Scell%Ei, Scell%Ne_low_VB, Scell%El_low_VB, Scell%mu_VB, &
                                          Te_VB, Te_start=Scell%TeeV, mu_start=Scell%mu, i_end_in=Scell%N_Egap) ! below
         Scell%Te_VB = Te_VB*g_kb ! save in [K]

         ! Construct equivalent Fermi distribution:
         call set_Fermi(Scell%Ei, Te_VB, Scell%mu_VB, Scell%fe_eq_VB, i_end=Scell%N_Egap)  ! below

         ! Make the thermalization for this band:
         if (numpar%tau_fe_VB < numpar%dt/30.0d0) then ! it's basically instantaneous
            exp_dttau = 0.0d0
         else  ! finite time relaxation
            exp_dttau = dexp(-numpar%dt / numpar%tau_fe_VB)
         endif
         do i = 1, Scell%N_Egap  ! for VB grid points (MO energy levels)
            Scell%fe(i) = Scell%fe_eq_VB(i) + (Scell%fe(i) - Scell%fe_eq_VB(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
         enddo

         ! CB:
         ! Get the number of particles and energy in the band:
         Scell%Ne_low_CB = get_N_partial(Scell%fe, Scell%N_Egap+1, i_fe) ! below
         Scell%El_low_CB = get_E_partial(Scell%Ei, Scell%fe,  Scell%N_Egap+1, i_fe) ! below

         ! Get the equivalent temperature and chem.potential of the band:
         call Electron_Fixed_Etot_partial(Scell%Ei, Scell%Ne_low_CB, Scell%El_low_CB, Scell%mu_CB, &
                                          Te_CB, Te_start=Scell%TeeV, mu_start=Scell%mu, i_start_in=Scell%N_Egap+1) ! below
         Scell%Te_CB = Te_CB*g_kb ! save in [K]

         ! Construct equivalent Fermi distribution:
         call set_Fermi(Scell%Ei, Te_CB, Scell%mu_CB, Scell%fe_eq_CB, i_start=Scell%N_Egap+1)  ! below

         ! Make the thermalization for this band:
         if (numpar%tau_fe_CB < numpar%dt/30.0d0) then ! it's basically instantaneous
            exp_dttau = 0.0d0
         else  ! finite time relaxation
            exp_dttau = dexp(-numpar%dt / numpar%tau_fe_CB)
         endif
         do i = Scell%N_Egap+1, i_fe ! for CB grid points (MO energy levels)
            Scell%fe(i) = Scell%fe_eq_CB(i) + (Scell%fe(i) - Scell%fe_eq_CB(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
         enddo

!          print*, 'Ne=', Scell%Ne_low_VB, Scell%Ne_low_CB, Scell%Ne_low_VB+Scell%Ne_low_CB, Scell%Ne_low
!          print*, 'Ne2=', get_N_partial(Scell%fe, 1, Scell%N_Egap), get_N_partial(Scell%fe_eq_VB, 1, Scell%N_Egap)
!          print*, 'Ne3=', get_N_partial(Scell%fe, Scell%N_Egap+1, i_fe), get_N_partial(Scell%fe_eq_CB, Scell%N_Egap+1, i_fe), get_N_partial(Scell%fe, 1, i_fe)
!          print*, 'Ee=', Scell%El_low_VB, Scell%El_low_CB, Scell%El_low_VB+Scell%El_low_CB, Scell%nrg%El_low
!          print*, 'Ee2=', get_E_partial(Scell%Ei, Scell%fe, 1, Scell%N_Egap), get_E_partial(Scell%Ei, Scell%fe,  Scell%N_Egap+1, i_fe), get_E_partial(Scell%Ei, Scell%fe, 1, i_fe)
!          print*, 'Te=', Te_VB, Te_CB, Scell%mu_VB, Scell%mu_CB
      endif BNDS

      !--------------------------
      ! II. Do complete thermalization
      if (numpar%tau_fe < numpar%dt/30.0d0) then ! it's basically instantaneous
         exp_dttau = 0.0d0
      else  ! finite time relaxation
         exp_dttau = dexp(-numpar%dt / numpar%tau_fe)
      endif
      do i = 1, i_fe ! for all grid points (MO energy levels)
         Scell%fe(i) = Scell%fe_eq(i) + (Scell%fe(i) - Scell%fe_eq(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
         !print*, i, Scell%fe(i), Scell%fe_eq(i), exp_dttau
      enddo

      !--------------------------
      ! III. Extra check for smoothening unphysical artefacts that may be present after MC:
      call smoothening_step(Scell, numpar) ! below

   endif ! (.not.skip_step)
end subroutine Do_relaxation_time


subroutine smoothening_step(Scell, numpar, eps_precision)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   real(8), intent(in), optional :: eps_precision     ! given precision
   !-------------------------------
   real(8) :: eps, Ne, Ne_eq, extra_dt, extra_tau, exp_dttau
   integer :: i, i_fe, i_cycle, N_cycle
   logical :: extra_cycle

   if (present(eps_precision)) then
      eps = eps_precision
   else ! use default value:
      eps = 1.0d-6
   endif

   i_fe = size(Scell%fe)   ! number of grid points in distribution function
   extra_cycle = .false.   ! by default, assume no artifact
   do i = 1, i_fe ! for all grid points (MO energy levels)
      if ((Scell%fe(i) > 2.0d0) .or. (Scell%fe(i) < 0.0d0)) then
         extra_cycle = .true.   ! some artefacts present, do extra thermalization to get rid of them
         print*, 'Extra thermalization step needed (1):', i, Scell%fe(i), SUM(Scell%fe)
         exit
      endif
   enddo

   ! Or, if the number of electrons does not coincide with the initial one:
   Ne = get_N_partial(Scell%fe, 1, i_fe)  ! current number of electrons
   Ne_eq = get_N_partial(Scell%fe_eq, 1, i_fe)  ! initial number of electrons
   if (ABS(Ne - Ne_eq) > eps*Ne_eq) then
      extra_cycle = .true. ! shouldn't be too different
      print*, 'Extra thermalization step needed (2):', Ne, Ne_eq
   endif

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

         !print*, 'Ne_extra=', i_cycle, SUM(Scell%fe), SUM(Scell%fe_eq), Scell%Ne_low, extra_cycle

         ! Also a consistency check for the number of electrons:
         Ne = get_N_partial(Scell%fe, 1, i_fe)
         Ne_eq = get_N_partial(Scell%fe_eq, 1, i_fe)
         if (ABS(Ne - Ne_eq) > eps*Ne_eq) extra_cycle = .true.
         if (i_cycle > N_cycle) exit
      enddo ! while (extra_cycle)
   endif ! (extra_cycle)
end subroutine smoothening_step



subroutine identify_fragment_for_orbital(Scell)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   !--------------
   integer :: Nat, Nlev, Norb, i, i_max
   real(8) :: Max_D_i

   Nat = size(Scell%MDAtoms)     ! number of atoms
   Nlev = size(Scell%Ei)         ! number of energy levels
   Norb = Nlev / Nat             ! number of orbitals per atom

   if (.not.allocated(Scell%orb_in_atom)) then
      allocate(Scell%orb_in_atom(Nlev), source = 0)
   endif

   ! For all energy levels, identify which fragment they belong by the maximal contribution of the atomic orbitals into this level:
   do i = 1, Nlev
      ! Which atom contributes to this orbital the most:
      !Max_D_i = maxval(Scell%Dmatrix(i,:))      ! atomic orbital with maximal contribution to this level
      !i_max = (transfer( maxloc(Scell%Dmatrix(i,:)), i_max) - 1)/ Norb + 1 ! find which atom this orbital belongs to
      i_max = (transfer( maxloc(Scell%Dmatrix(:,i)), i_max) - 1)/ Norb + 1 ! find which atom this orbital belongs to

      ! Which atom it belongs to:
      Scell%orb_in_atom(i) = i_max

      ! reminder:
      !Scell%fragments%indices(i_max)     ! this identifies the fragment
   enddo
end subroutine identify_fragment_for_orbital



subroutine get_fragments_data_for_electrons(Scell, NSC, numpar, matter)
   type(Super_cell), dimension(:), intent(inout) :: Scell   ! suoer-cell with all the atoms inside
   integer, intent(in) :: NSC                   ! number of the super-cell
   type(Numerics_param), intent(inout) :: numpar   ! numerical parameters
   type(solid), intent(inout) :: matter         ! material parameters
   !-----------------------------------
   real(8) :: Te_frag, mu_frag, Ne, Ee
   integer :: i, Nsiz, j, Nlev, i_at, i_orb, N_at, N_orb, N_orb_tot
   logical :: needs_allocation
   logical, dimension(:), allocatable :: orbital_fragments
   !logical, dimension(size(Scell(NSC)%MDAtoms)) :: fragment_mask
   logical, dimension(:), allocatable :: fragment_mask

   !real(8), dimension(:), allocatable :: fe_eq_frag

   ! Number of different fragments:
   Nsiz = maxval(Scell(NSC)%fragments%indices)  ! what's the maximal fragment index among all the atoms
   N_at = size(Scell(NSC)%MDAtoms)       ! total number of atoms
   N_orb_tot = size(Scell(NSC)%Ha,1)     ! total number of orbitals
   N_orb = N_orb_tot/N_at           ! orbitals per atom

   allocate(fragment_mask(N_at))

   if (.not.allocated(Scell(NSC)%orb_in_atom)) then ! it's a first call
      call identify_fragment_for_orbital(Scell(NSC))  ! above
   endif

   !-----------------------------------
   ! Redefine the arrays of electronic parameters for each fragment:
   needs_allocation = .false. ! to start with
   ! 1) Array of number of electrons:
   if (.not.allocated(Scell(NSC)%fragments%N_e)) then
      needs_allocation = .true.
   else ! may need reallocation, if the number of fragments changed:
      if (size(Scell(NSC)%fragments%N_e) /= Nsiz) then
         deallocate(Scell(NSC)%fragments%N_e)  ! wrong size, reallocate
         needs_allocation = .true.
      endif
   endif
   if (needs_allocation) then ! make sure it is allocated with the correct size
      allocate(Scell(NSC)%fragments%N_e(Nsiz))
   endif
   ! 2) Array of electron energy:
   if (.not.allocated(Scell(NSC)%fragments%E_e)) then
      needs_allocation = .true.
   else ! may need reallocation, if the number of fragments changed:
      if (size(Scell(NSC)%fragments%E_e) /= Nsiz) then
         deallocate(Scell(NSC)%fragments%E_e)  ! wrong size, reallocate
         needs_allocation = .true.
      endif
   endif
   if (needs_allocation) then ! make sure it is allocated with the correct size
      allocate(Scell(NSC)%fragments%E_e(Nsiz))
   endif
   ! 3) Array of electronic chemical potential:
   if (.not.allocated(Scell(NSC)%fragments%mu)) then
      needs_allocation = .true.
   else ! may need reallocation, if the number of fragments changed:
      if (size(Scell(NSC)%fragments%mu) /= Nsiz) then
         deallocate(Scell(NSC)%fragments%mu)  ! wrong size, reallocate
         needs_allocation = .true.
      endif
   endif
   if (needs_allocation) then ! make sure it is allocated with the correct size
      allocate(Scell(NSC)%fragments%mu(Nsiz))
   endif
   ! 4) Array of electronic temperatures:
   if (.not.allocated(Scell(NSC)%fragments%T_e)) then
      needs_allocation = .true.
   else ! may need reallocation, if the number of fragments changed:
      if (size(Scell(NSC)%fragments%T_e) /= Nsiz) then
         deallocate(Scell(NSC)%fragments%T_e)  ! wrong size, reallocate
         needs_allocation = .true.
      endif
   endif
   if (needs_allocation) then ! make sure it is allocated with the correct size
      allocate(Scell(NSC)%fragments%T_e(Nsiz))
   endif

   !-----------------------------------
   ! Identify the electronic parameters for each fragment,
   ! use Mulliken analysis to determinre electronic populations on particular atoms:

   Scell(NSC)%fragments%N_e(:) = 0.0d0 ! restart counting
   Scell(NSC)%fragments%E_e(:) = 0.0d0 ! restart counting

   do i = 1, Nsiz
      !Scell(NSC)%fragments%N_e(i) = SUM( Scell(NSC)%fe(:), MASK = orbital_fragments(:) ) ! TOO CRUDE, DESON'T WORK WELL

      ! Update mask for atoms - which ones belong to this fragment "i":
      fragment_mask(:) = (Scell(NSC)%fragments%indices(:) == i)
      Ne = 0.0d0  ! to start with
      Ee = 0.0d0  ! to start with

      !$omp PARALLEL private(i_at, i_orb, j)
      !$omp do reduction(+ : Ne, Ee)
      do i_at = 1, N_at ! all atoms
         do i_orb = 1, N_orb  ! all orbitals of each atom
            j = (i_at-1)*N_orb + i_orb ! current orbital among all

            if (fragment_mask(i_at)) then ! this atom belongs to this fragment
               ! 1) Number of electrons in the fragment "i":
               !Ne = Ne + SUM(Scell(NSC)%fe(:) * Scell(NSC)%Dmatrix(j,:))     ! incorrect
               Ne = Ne + SUM(Scell(NSC)%fe(:) * Scell(NSC)%Dmatrix(:,j))      ! correct, tested

               ! 2) Electrons energy in the fragment "i"
               !Ee = Ee + SUM(Scell(NSC)%fe(:) * Scell(NSC)%Ei(:) * Scell(NSC)%Dmatrix(j,:))    ! incorrect
               Ee = Ee + SUM(Scell(NSC)%fe(:) * Scell(NSC)%Ei(:) * Scell(NSC)%Dmatrix(:,j))     ! correct, tested
            endif
         enddo   ! i_orb
      enddo ! i_at
      !$omp end do
      !$omp end parallel

      Scell(NSC)%fragments%N_e(i) = Ne
      ! Normalizatin per atom in the fragment:
      Scell(NSC)%fragments%E_e(i) = Ee / dble(Scell(NSC)%fragments%N_at(i))  ! [eV/atom]

      ! 3) Electronic chemical potential and temperature in the fragment "i"
      ! Get the equivalent temperature and chem.potential of electrons in this fragment:
      !call Electron_Fixed_Etot_partial(Scell(NSC)%Ei, Scell(NSC)%fragments%N_e(i), Scell(NSC)%fragments%E_e(i), mu_frag, &
      !      Te_frag, Te_start=Scell(NSC)%TeeV, mu_start=Scell(NSC)%mu, orbital_fragments=orbital_fragments) ! below
      !Scell(NSC)%fragments%T_e(i) = Te_frag*g_kb      ! save in [K]
      !Scell(NSC)%fragments%mu(i) = mu_frag            ! [eV]

      !print*, i, Scell(NSC)%fragments%N_at(i), Scell(NSC)%fragments%N_e(i), Scell(NSC)%fragments%E_e(i), SUM(matter%Atoms(Scell(NSC)%MDAtoms(:)%KOA )%NVB - Scell(NSC)%MDAtoms(:)%q, MASK=fragment_mask(:))
   enddo ! i


   ! Consistency check:
   if ( abs(sum(Scell(NSC)%fragments%N_e(:)) - Scell(NSC)%Ne) > 1.0d-4*Scell(NSC)%Ne ) then
      print*, 'PROBLEM in get_fragments_data_for_electrons:'
      print*, 'electrons in fragments do not add up to total Ne:', sum(Scell(NSC)%fragments%N_e(:)), Scell(NSC)%Ne
   endif

   if (allocated(orbital_fragments)) deallocate(orbital_fragments)
   if (allocated(fragment_mask)) deallocate(fragment_mask)
end subroutine get_fragments_data_for_electrons


!    allocate(fe_eq_frag(size(Scell(NSC)%Ei)), source = 0.0d0)
!       ! Construct equivalent Fermi distribution:
!       call set_Fermi(Scell%Ei, Te_frag, mu_frag, fe_eq_frag, orbital_fragments=orbital_fragments)  ! below
!
!       ! Make the thermalization for this band:
!       if (numpar%tau_fe_VB < numpar%dt/30.0d0) then ! it's basically instantaneous
!          exp_dttau = 0.0d0
!       else  ! finite time relaxation
!          exp_dttau = dexp(-numpar%dt / numpar%tau_fe_VB)
!       endif
!       do i = 1, Scell%N_Egap  ! for VB grid points (MO energy levels)
!          Scell%fe(i) = Scell%fe_eq_VB(i) + (Scell%fe(i) - Scell%fe_eq_VB(i))*exp_dttau   ! exact solution of df/dt=-(f-f0)/tau
!       enddo


!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! General subroutines that are often used:

subroutine electronic_entropy(fe, Se, norm_fe, i_start, i_end)
   real(8), dimension(:), intent(in) :: fe ! electron distribution function
   real(8), intent(out) :: Se ! self-explanatory
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer, intent(in), optional :: i_start, i_end  ! starting and ending levels to include
   ! Se = -kB * int [ DOS*( f * ln(f) + (1-f) * ln(1-f) ) ]
   ! E.G. [https://doi.org/10.1103/PhysRevB.50.14686]
   !----------------------------
   !real(8), dimension(size(fe)) :: f_lnf
   real(8), dimension(:), allocatable :: f_lnf
   real(8) :: f_norm, eps
   integer :: i, Nsiz, i_low, i_high
   !============================
   eps = 1.0d-12  ! precision
   Nsiz = size(fe)

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

   if (present(norm_fe)) then   ! user provided
      f_norm = norm_fe
   else ! by default, not spin resolved
      f_norm = 2.0d0
   endif

   ! To start with:
   allocate(f_lnf(Nsiz))
   Se = 0.0d0
   f_lnf = 0.0d0

   ! First term of the total entropy:
   ! (our f is normalized to f_norm, which means it includes DOS in it, so divide by f_norm where needed)
   where (fe(i_low:i_high) > eps) f_lnf(i_low:i_high) = fe(i_low:i_high)*log(fe(i_low:i_high)/f_norm)
   ! First term:
   Se = SUM(f_lnf(i_low:i_high))
   f_lnf = 0.0d0  ! fill the empty points
   ! Second term of the total entropy:
   where (fe(i_low:i_high) < f_norm-eps) f_lnf(i_low:i_high) = (f_norm - fe(i_low:i_high))*log((f_norm - fe(i_low:i_high))/f_norm)
   ! Assemble the terms:
   Se = Se + SUM(f_lnf(i_low:i_high))

   ! Make proper units:
   Se = -g_kb_EV*Se  ! [eV/K]

   ! clean up:
   deallocate(f_lnf)
end subroutine  electronic_entropy



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



subroutine set_Fermi(Ei,Te,mu,fe, Error_desript, norm_fe, i_start, i_end)
   real(8), dimension(:), intent(in) :: Ei	! energy levels, eigenvalues of the hamiltonian matrix
   real(8), intent(in) :: Te, mu	! temperature [eV], and chemical potential [eV]
   real(8), dimension(:), intent(out) :: fe
   character(*), intent(out), optional :: Error_desript ! description of error if occured
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer, intent(in), optional :: i_start, i_end
   integer :: Nsiz, i, n, i_low, i_high
   real(8) :: f_norm, arg
   

   if (present(i_start)) then
      i_low = i_start
   else  ! default, start from 1
      i_low = 1
   endif
   Nsiz = size(Ei)
   if (present(i_end)) then
      i_high = i_end
   else  ! default, end at the end
      i_high = Nsiz
   endif

   if (present(norm_fe)) then    ! user provided
      f_norm = norm_fe
   else  ! by default, spin degenerate
      f_norm = 2.0d0
   endif


   fe(:) = 0.0d0  ! to start with
   if (Te > 0.0d0) then
      !where ((Ei(i_low:i_high) - mu)/Te >= log(HUGE(mu))) ! exp(x) -> infinity
      !   fe(:) = 0.0d0
      !else where
      !   fe(:) = f_norm/(1.0d0 + exp((Ei(i_low:i_high) - mu)/Te)) ! fermi-function
      !end where
      do i = i_low, i_high
         arg = (Ei(i) - mu)/Te
         if (arg <= log(huge(0.0d0))) then
            fe(i) = f_norm/(1.0d0 + exp(arg))
         endif
      enddo
   elseif (Te == 0.0d0) then
      where (Ei(i_low:i_high) <= mu)
         fe(:) = f_norm
      else where
         fe(:) = 0.0d0
      end where
   else  !Te < 0
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




subroutine get_orbital_resolved_data(Scell, matter, DOS_weights, numpar)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(solid), intent(in) :: matter   ! materil parameters
   real(8), dimension(:,:,:), intent(in) :: DOS_weights ! weigths of the particular type of orbital on each energy level
   type(Numerics_param), intent(in) :: numpar  ! numerical parameters
   !-----------------
   integer :: Nat, Nsh, i, j, Nsiz
   real(8) :: N_MDat

   ! Number of kinds of atoms:
   Nat = size(matter%Atoms)
   ! Numer of atoms in the supercell:
   N_MDat = dble(Scell%Na)

   !N_at = size(DOS_weights,1)    ! number of kinds of atoms
   ! consistency test:
   if (Nat /= size(DOS_weights,1)) then
      print*, 'Problem in get_orbital_resolved_data, inconsistent sizes:', Nat, size(DOS_weights,1)
      print*, 'Cannot get orbital-resolved electron data.'
      return
   endif
   ! Make sure output data are allocated:
   if (.not.allocated(Scell%Orb_data)) allocate(Scell%Orb_data(Nat))

   ! For all kinds of atoms
   do i = 1, Nat
      Nsh = size(DOS_weights,2) ! number of atomic shells (basis set size)
      ! Make sure output data are allocated:
      if (.not.allocated(Scell%Orb_data(i)%Ne)) allocate(Scell%Orb_data(i)%Ne(Nsh), source = 0.0d0)
      if (.not.allocated(Scell%Orb_data(i)%Ee)) allocate(Scell%Orb_data(i)%Ee(Nsh), source = 0.0d0)
      if (numpar%save_fe_orb) then
         Nsiz = size(Scell%fe)
         if (.not.allocated(Scell%Orb_data(i)%fe)) allocate(Scell%Orb_data(i)%fe(Nsh, Nsiz), source = 0.0d0)
      endif
      ! For all orbitals:
      do j = 1, Nsh
         ! Number of electrons in each orbital:
         Scell%Orb_data(i)%Ne(j) = SUM( Scell%fe(:) * DOS_weights(i, j, :) )
         ! Energy of electrons in each orbital:
         Scell%Orb_data(i)%Ee(j) = SUM( Scell%Ei(:) * Scell%fe(:) * DOS_weights(i, j, :) )
         ! Normalize the data per atom:
         Scell%Orb_data(i)%Ne(j) = Scell%Orb_data(i)%Ne(j) / N_MDat   ! [electron/atom]
         Scell%Orb_data(i)%Ee(j) = Scell%Orb_data(i)%Ee(j) / N_MDat   ! [eV/atom]
         ! Orbital-resolved distribution function, if requested:
         if (numpar%save_fe_orb) then
            Scell%Orb_data(i)%fe(j,:) = Scell%fe(:) * DOS_weights(i, j, :)
         endif
      enddo ! j
   enddo ! i

end subroutine get_orbital_resolved_data




subroutine get_electronic_heat_capacity(Scell, NSC, Ce, do_kappa, DOS_weights, Ce_partial, norm_fe, mu_in, Te_in)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), intent(out) :: Ce ! current electron heat capacity [J/(m^3 K)]
   logical, intent(in) :: do_kappa  ! if kappa calculations are requested
   real(8), dimension(:,:,:), intent(in), optional :: DOS_weights ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:), intent(out), allocatable, optional :: Ce_partial ! band-resolved Ce [J/(m^3 K)]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   real(8), intent(in), optional :: mu_in, Te_in   ! [eV], [K], provided chemical potential and temperature
   !-----------------
   real(8) :: Ntot	! number of electrons
   real(8) :: nat   ! number of atoms
   real(8) :: Te	! current electron temperature [eV]
   real(8) :: mu	! current electron chemical potential [eV]
   real(8) :: dmu, dTe, mu0	! electron differential chemical potential [eV], temperature [eV]
   real(8) :: Dens	 ! atomic density
   real(8) :: coef   ! conversion coefficients with units
   real(8) :: C1, C2
   logical :: do_partial
   !real(8), dimension(size(Scell(NSC)%Ei)) :: Ce_i
   real(8), dimension(:), allocatable :: Ce_i

   if (present(DOS_weights) .and. present(Ce_partial)) then ! partial contributions required:
      do_partial = .true.
      if (.not.allocated(Ce_partial)) allocate(Ce_partial(size(Scell(NSC)%G_ei_partial,1)))
   else
      do_partial = .false.
   endif

   allocate(Ce_i(size(Scell(NSC)%Ei)))

   dTe = 10.0d0/g_kb    ! [eV] -> [K]
   Ntot = dble(Scell(NSC)%Ne)
   nat = dble(Scell(NSC)%Na) ! number of atoms

   ! 1) Low-energy electrons populating TB-band structure:
   if (present(Te_in)) then   ! use the provided value
      Te = Te_in/g_kb    ! [eV] -> [K]
   else ! use the default one from the supercell
      Te = Scell(NSC)%TeeV
   endif
   if (present(mu_in)) then   ! use the provided value
      mu = mu_in
   else ! use the default one from the supercell
      mu = Scell(NSC)%mu
   endif

   if (present(norm_fe)) then
      call Electron_Fixed_Te(Scell(NSC)%Ei, Ntot, mu0, Te+dTe, norm_fe) ! in case if the electron temperature is given
   else
      call Electron_Fixed_Te(Scell(NSC)%Ei, Ntot, mu0, Te+dTe) ! in case if the electron temperature is given
   endif
   dmu = (mu0 - mu)/dTe
   if (present(norm_fe)) then ! normalization of fe provided:
      if (do_partial) then
         call Get_Ce(Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, Ce_i, DOS_weights, Ce_partial, norm_fe)   ! below
      else  ! no partial contributions required:
         call Get_Ce(Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, Ce_i, norm_fe=norm_fe)   ! below
      endif
   else  ! default normalization of fe:
      if (do_partial) then
         call Get_Ce(Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, Ce_i, DOS_weights, Ce_partial)   ! below
      else  ! no partial contributions required:
         call Get_Ce(Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce, Ce_i)   ! below
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

   ! clean up:
   deallocate(Ce_i)
end subroutine get_electronic_heat_capacity




subroutine get_Ce_and_mu(Scell, NSC, Te_in, Ei, Ce, mu, mu_on_gamma)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), intent(in) :: Te_in   ! [K] electronic temperature
   real(8), dimension(:), intent(in) :: Ei   ! [eV] electronic energy levels
   real(8), intent(out) :: Ce ! electron heat capacity [J/(m^3 K)]
   real(8), intent(out) :: mu ! electron chemical potential [eV]
   logical, intent(in), optional :: mu_on_gamma ! calculate chem.pot. in gamma point (or on given k-point)
   !-----------------
   real(8) :: Ntot   ! number of electrons
   real(8) :: nat    ! number of atoms
   real(8) :: Te     ! current electron temperature [eV]
   real(8) :: dmu, dTe, mu0   ! electron differential chemical potential [eV], temperature [eV]
   real(8) :: Dens   ! atomic density
   real(8) :: coef   ! conversion coefficients with units
   logical :: mu_gamma

   if (present(mu_on_gamma)) then
      mu_gamma = mu_on_gamma
   else  ! by default, get mu in gamma point only
      mu_gamma = .true.
   endif

   Te = Te_in/g_kb      ! [eV] <- [K]
   ! Step in temperature for (d mu/ d Te):
   dTe = 10.0d0/g_kb    ! [eV] <- [K]
   Ntot = dble(Scell(NSC)%Ne) ! number of electrons
   nat = dble(Scell(NSC)%Na)  ! number of atoms

   ! Numerical derivative of chem.pot.:
   if (mu_gamma) then ! on gamma-point:
      ! Get mu:
      call Electron_Fixed_Te(Scell(NSC)%Ei, Ntot, mu, Te) ! below
      ! Get mu0:
      call Electron_Fixed_Te(Scell(NSC)%Ei, Ntot, mu0, Te+dTe) ! below
      dmu = (mu0 - mu)/dTe    ! (d mu/ d Te)
!       ! Get Ce in arb.units:
!       call Get_Ce(Scell(NSC)%Ei, Te+dTe/2.0d0, mu, dmu, Ce) ! below
   else ! on given k-point
      ! Get mu:
      call Electron_Fixed_Te(Ei, Ntot, mu, Te) ! below
      ! Get mu0:
      call Electron_Fixed_Te(Ei, Ntot, mu0, Te+dTe) ! below
      dmu = (mu0 - mu)/dTe    ! (d mu/ d Te)
   endif

   ! Get Ce in arb.units:
   call Get_Ce(Ei, Te+dTe/2.0d0, mu, dmu, Ce) ! below

   coef = 1.0d30*g_e/g_kb  ! [eV/A^3] -> [J/m^3/K]
   Dens = 1.0d0/(Scell(NSC)%V) ! [1/A^3]
   Ce = Ce * Dens*coef  ! [J/(m^3 K)]

   ! If undefined or infinite:
   if (isnan(Ce) .or. abs(Ce) >= 1d30) Ce = 0.0d0
end subroutine get_Ce_and_mu



subroutine Get_Ce(Ei, Te, mu, dmu, C, Ce_i, DOS_weights, Ce_partial, norm_fe)
    real(8), dimension(:), intent(in) :: Ei
    real(8), intent(in) :: Te, mu, dmu
    real(8), intent(out) :: C ! heat capacity [arb.units]
    real(8), dimension(:), intent(out), optional :: Ce_i
    real(8), dimension(:,:,:), intent(in), optional :: DOS_weights ! weigths of the particular type of orbital on each energy level
    real(8), dimension(:), intent(out), optional :: Ce_partial ! band-resolved Ce [arb.units]
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
        if (present(Ce_i)) Ce_i(i) = C_temp   ! save for output

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
   real(8) :: dfdE   ! Derivative of the Fermi-function by energy [1/eV]
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



subroutine Electron_Fixed_Etot_partial(Ei, Netot, Eetot, mu, Te, mu_start, Te_start, i_start_in, i_end_in, orbital_fragments)
   real(8), dimension(:), intent(in) ::  Ei  ! eigenvalues of TB-Hamiltonian for electrons
   real(8), intent(in) :: Netot  ! number of electrons/supercell to normalize the distribution function
   real(8), intent(in) :: Eetot  ! energy of electrons/supercell to normalize the distribution function
   real(8), intent(inout) :: Te ! electron temperature [eV]
   real(8), intent(inout) :: mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: mu_start, Te_start    ! initial guess for the electron chem.potential and temperature [eV]
   integer, intent(in), optional :: i_start_in, i_end_in
   logical, dimension(:), intent(in), optional :: orbital_fragments
   !-------------------------
   real(8) :: mu_0, mu_1, Te_0, Te_1, mix_fact, eps, mu_diff, Te_diff, eps_Te, eps_mu
   integer :: Nsiz, i_start, i_end, coun
   logical :: cycle_continue, do_steps

   if (present(i_start_in)) then
      i_start = i_start_in
   else  ! default, start from 1
      i_start = 1
   endif
   Nsiz = size(Ei)
   if (present(i_end_in)) then
      i_end = i_end_in
   else  ! default, end at the end
      i_end = Nsiz
   endif

   mix_fact = 0.35d0          ! empirically chosen mixing parameter for self-consistent cycle
   eps = 1.0d-7               ! relative precision for chem.potential
   eps_Te = 10.0d0/g_kb       ! [eV] 10 K, absolite allowed precision for the temperatures
   eps_mu = 1.0d-3            ! [eV] absolute allowed precision for chem.potential

   ! Starting cycle:
   if (present(Te_start)) then
      Te_0 = Te_start   ! [eV] electron temperature to be calculated
   else ! use default starting value
      Te_0 = 1.0d0/g_kb   ! [eV] electron temperature to be calculated
   endif
   if (present(mu_start)) then
      mu_0 = mu_start   ! [eV] electron temperature to be calculated
   else ! find it from the given electron temperature
      call Electron_Fixed_Te(Ei, Netot, mu_0, Te_0, i_start=i_start, i_end=i_end) ! find partial mu for given Te and Ne
   endif

   ! Get Te for this mu:
   call Electron_Fixed_mu(Ei, Eetot, mu_0, Te_1, i_start=i_start, i_end=i_end) ! find partial Te for given mu and Ee

   coun = 0 ! to start with
   do_steps =.false. ! to start with

   !mu_1 = mu_0 ! to start with
   ! Find chem.potential for this temperature:
   call Electron_Fixed_Te(Ei, Netot, mu_1, Te_1, i_start=i_start, i_end=i_end) ! find partial mu for given Te and Ne

   !print*, 'SC1', coun, Te_0*g_kb, Te_1*g_kb, mu_0, mu_1, cycle_continue, Eetot

   !Te_1 = Te_0 * (1.0d0-mix_fact) + mix_fact * Te_1   ! mixing for the next step
   Te_1 = min( Te_0*(1.0d0-mix_fact) + mix_fact*Te_1 , 2.0d0*Te_0 )  ! mixing for the next step
   mu_diff = abs(mu_0 - mu_1)
   Te_diff = abs(Te_0 - Te_1)
   cycle_continue = ( (Te_diff > eps*(max(Te_0,Te_1))) .or. (mu_diff > eps*(max(mu_0,mu_1))) )

   do while ( cycle_continue )
      coun = coun + 1
      mu_0 = mu_1
      Te_0 = Te_1

      ! Get mu:
      call Electron_Fixed_Te(Ei, Netot, mu_1, Te_1, i_start=i_start, i_end=i_end) ! find partial mu for given Te and Ne
      mu_1 = mu_0 * (1.0d0-mix_fact) + mix_fact * mu_1   ! mixing for the next step
      ! Limit the chem.potential step by 50% of the current value:
      !mu_1 = min( mu_0 * (1.0d0-mix_fact) + mix_fact * mu_1, 2.0d0*mu_0 )   ! mixing for the next step

      ! Get Te:
      call Electron_Fixed_mu(Ei, Eetot, mu_1, Te_1, i_start=i_start, i_end=i_end) ! find partial Te for given mu and Ee
      !Te_1 = Te_0 * (1.0d0-mix_fact) + mix_fact * Te_1   ! mixing for the next step
      ! Limit the temperature step by 50% of the current temeprature:
      Te_1 = min( Te_0*(1.0d0-mix_fact) + mix_fact*Te_1 , 2.0d0*Te_0 )  ! mixing for the next step

      ! Check if mixing should be reduced:
      if ( (abs(Te_0 - Te_1) > max(Te_diff, Te_0*0.1d0) ) .or. &
           (abs(mu_0 - mu_1) > max(mu_diff, mu_0*0.1d0) ) ) then      ! difference not smaller than 10% of the value itself
         if (mix_fact > 0.1d0) mix_fact = mix_fact*0.9d0
         !mu_1 = (mu_0+mu_1)*0.5d0
         !Te_1 = (Te_0+Te_1)*0.5d0
      endif

      ! Check if we need to continue the cycle:
      mu_diff = abs(mu_0 - mu_1)
      Te_diff = abs(Te_0 - Te_1)
      ! Relative precision:
      cycle_continue = ( (Te_diff > eps*(max(Te_0,Te_1))) .or. (mu_diff > eps*(max(mu_0,mu_1))) )

      ! Absolute precision:
      cycle_continue = ( (Te_diff > eps_Te) .or. (mu_diff > eps_mu) )

      !print*, 'SC2', coun, Te_0*g_kb, Te_1*g_kb, mu_0, mu_1, cycle_continue, mix_fact

      ! It's enough of iterations:
      if (coun > 100) then
         !if ( Te_diff > 10.0d0 * eps*(max(Te_0,Te_1)) ) then ! printout only if the difference is noticeable
         !   write(*,'(a,i0)') 'Too many iterations in Electron_Fixed_Etot_partial: ', coun
         !   write(*,'(f,f,f,f)') Te_0, Te_1, mu_0, mu_1
         !endif
         ! If still didn't find the temperature and chem.potential, do it differently:
         do_steps = .true.
         exit
      endif
      ! And make sure there is no runaway:
      if (abs(mu_1) > 1.0d10) then
         write(*,'(a)') 'Runaway condition reached in Electron_Fixed_Etot_partial: '
         write(*,'(f,f,f,f)') Te_0, Te_1, mu_0, mu_1
         exit
      endif
   enddo

   ! Output:
   Te = (Te_0 + Te_1)*0.5d0
   mu = (mu_0 + mu_1)*0.5d0

   !print*, 'one', mu, Te

   ! Check if another method of finding Te and mu is required:
   if (do_steps) then
      call Electron_Fixed_Etot_partial_1(Ei, Netot, Eetot, mu, Te, i_start_in=i_start, i_end_in=i_end) ! below
   endif

   !print*, 'two', mu, Te
end subroutine Electron_Fixed_Etot_partial




subroutine Electron_Fixed_Etot_partial_1(wrD, Netot, Eetot, mu, Te, i_start_in, i_end_in)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(in) :: Eetot  ! energy of electrons/supercell to normalize the distribution function
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]
   integer, intent(in), optional :: i_start_in, i_end_in
   real(8) :: a, b, Ncur, Ecur, Elast, muN, muE, Telast, Telast2, popul, che1, che2
   real(8) :: curMu, minmunmue, dT
   integer :: Nsiz, i, j, cou, countr, i_start, i_end

   if (present(i_start_in)) then
      i_start = i_start_in
   else  ! default, start from 1
      i_start = 1
   endif
   Nsiz = size(wrD)
   if (present(i_end_in)) then
      i_end = i_end_in
   else  ! default, end at the end
      i_end = Nsiz
   endif

   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
   ! precalculations before finding mu and Te
   dT = 50.0d0    ! temperature step [K]
   Te = dT/g_kb   ! [eV] electron temperature to be calculated
   call Electron_Fixed_Te(wrD, Netot, muN, Te, i_start=i_start, i_end=i_end) ! find partial mu for given Te and Ne
   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

   !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
   ! finding mu and Te
   countr = 0
4370 cou = 0
   Telast = Te
   Telast2 = Telast + dT/g_kb
   Te = (Telast+Telast2)/2.0d0
   Ecur = -1d12
   che1 = (Ecur-Eetot)  ! checker for finding the crossing point
   do while ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-9)) .GT. 1d-9) .AND. (ABS((Te - Telast)/MAX(Te,Telast,1d-12)) .GT. 1d-12)) ! main while
      !write(*,'(a,es,es,f,f)') 'T:', Ecur, Eetot, Te, muN
      Te = (Telast+Telast2)/2.0d0 ! [eV]
      if (ABS((Te - Telast)/MAX(Te,Telast,1d-12)) .LT. 1d-12) exit

      ! Finding chem.potential for the neq temperature:
      call Electron_Fixed_Te(wrD, Netot, muN, Te, i_start=i_start, i_end=i_end) ! find mu for given Te and Ne

      ! Checking, whether this chem.potential and temperature are also solutions for energy:
      Ecur = get_E_tot(wrD, muN, Te, i_start=i_start, i_end=i_end) ! function from above

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

   if ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-9)) .GT. 1d-9) .AND. (countr .LT. 100))then
      countr = countr + 1
      Te = Te + 2.0d0*dT/g_kb
      goto 4370 ! just try it over again if it didn't work well the first time...
   endif
end subroutine Electron_Fixed_Etot_partial_1




subroutine Electron_Fixed_mu(wrD, Eetot, mu, Te, norm_fe, i_start, i_end) ! For given chemical potential and total energy, find Te
   real(8), dimension(:), intent(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Eetot  ! number of electrons/atom to normalize the distribution function
   REAL(8), INTENT(out) :: Te ! electron temperature [eV]
   REAL(8), INTENT(in) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer, intent(in), optional :: i_start, i_end
   !-----------------------
   real(8) :: a, b, Ecur, Tcur, Elast, Emin, Emax, E_sign, dT
   integer :: Nsiz, i, i_low, i_high, N_steps

   if (present(i_start)) then
      i_low = i_start
   else  ! default, start from 1
      i_low = 1
   endif
   Nsiz = size(wrD)
   if (present(i_end)) then
      i_high = i_end
   else  ! default, end at the end
      i_high = Nsiz
   endif

   a = 0.1d0/g_kb ! [eV]
   b = 20.0d0 ! [eV}
   Ecur = 1.0d10   ! to start with
   Te = 0.0d0  ! to start with
   ! Energy on the borders of the search interval:
   Emin = get_E_tot(wrD, mu, a, i_start=i_low, i_end=i_high) ! function from below
   Emax = get_E_tot(wrD, mu, b, i_start=i_low, i_end=i_high) ! function from below


   ! Step one: find the interval where the function crosses the Eetot:
   N_steps = 1000  ! not more than this number of steps for interval search
   dT = (b - a)/dble(N_steps)   ! step to use for this search
   do i = 1, N_steps
      Te = a + dble(i) * dT
      Emax = get_E_tot(wrD, mu, Te, i_start=i_low, i_end=i_high) ! function from below
      ! check if both border values are on the same side of the value we search,
      ! both greater or both smaller. If it is the case, keep searching.
      ! If on is smaller, the other one is bigger, than this is the interval
      ! where the function crosses the value we are looking for.
      E_sign = (Emin-Eetot) * (Emax-Eetot)

      !print*, 's --->', i, Emin, Emax, Eetot, E_sign, mu, Te

      if (E_sign > 0.0d0) then
         ! both on the same size, keep searching
      else
         a = Te - dT
         b = Te   ! found the upper bound
         exit
      endif
   enddo
   ! use found boundary of the interval:
   Emin = get_E_tot(wrD, mu, a, i_start=i_low, i_end=i_high) ! function from below
   Elast = Emin   ! value from the last step

   ! Second step: bisection method to find the exact value from the defined interval:
   do while (ABS(Ecur-Eetot) .GT. Eetot*1d-12)
      Te = (a+b)/2.0d0
      Ecur = get_E_tot(wrD, mu, Te, i_start=i_low, i_end=i_high) ! function from below

      !print*, 'p --->', Ecur, Eetot, a, b

      if (Ecur > Emin) then ! increasing function:
         if (Ecur .GT. Eetot) then
            b = Te
         else
            a = Te
         endif
      else ! decareasing function
         if (Ecur .GT. Eetot) then
            a = Te
         else
            b = Te
         endif
      endif

      !if (ABS(a-b) .LT. 1d-12) exit ! it's too close anyway...
      if (ABS(a-b) .LT. 1d-10) exit ! it's too close anyway...
   enddo ! while

   !pause 'Electron_Fixed_mu'
end subroutine Electron_Fixed_mu





subroutine Electron_Fixed_Te(wrD, Netot, mu, Te, norm_fe, i_start, i_end) ! in case if the electron temperature is given
! but the total energy can change, then we have to find the chem.potential that conserves 
! the given total number of electons Netot
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Netot  ! number of electrons/atom to normalize the distribution function
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(out) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer, intent(in), optional :: i_start, i_end
   real(8) :: a, b, Ncur
   integer :: Nsiz, i, i_low, i_high

   Nsiz = size(wrD)

   if (present(i_start)) then
      i_low = i_start
      if (i_low > Nsiz) i_low = Nsiz
   else  ! default, start from 1
      i_low = 1
   endif

   if (present(i_end)) then
      i_high = i_end
      if (i_high > Nsiz) i_high = Nsiz
   else  ! default, end at the end
      i_high = Nsiz
   endif

   !a = wrD(i_low) - 100.0d0
   !b = wrD(i_high) + 50.0d0
   a = wrD(i_low) - 50.0d0
   b = wrD(i_high) + 50.0d0
   Ncur = 0.0d0

   do while (ABS(Ncur-Netot) .GT. Netot*1d-12) !
      mu = (a+b)/2.0d0
      !Ncur = get_N_tot(wrD, mu, Te) ! function from below
      Ncur = get_N_tot(wrD, mu, Te, i_start=i_low, i_end=i_high) ! function from below
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
      Ne = 2.0d0 * SUM(1.0d0/(1.0d0 + dexp((wr(i_start:i_end) - mu)/Te))) ! Fermi-function
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
      Ee = 2.0d0 * SUM(wr(i_start:i_end), mask = wr(:) <= mu) ! Fermi-function at T=0
   endif
end function get_E_partial_Fermi




function get_N_tot(wrD, mu, Te, norm_fe, i_start, i_end)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(in) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer, intent(in), optional :: i_start, i_end  ! starting and ending levels to include
   real(8) :: get_N_tot ! number of electrons for given band structure and Fermi-function
   !----------------
   real(8) :: f_norm
   integer :: i, Nsiz, i_low, i_high

   if (present(i_start)) then
      i_low = i_start
   else  ! default, start from 1
      i_low = 1
   endif
   Nsiz = size(wrD)
   if (present(i_end)) then
      i_high = i_end
   else  ! default, end at the end
      i_high = Nsiz
   endif

   
   if (Te .GT. 1.0d-12) then
      get_N_tot = 0.0d0
      DO_SUM:do i = i_low, i_high
         if ((wrD(i) - mu)/Te .LT. log(HUGE(mu))) then ! exp(x) -> infinity
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
      get_N_tot = f_norm*COUNT(wrD(i_low:i_high) <= mu)
   endif
end function get_N_tot


function get_E_tot(wrD, mu, Te, norm_fe, i_start, i_end)
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(in) :: Te ! electron temperature [eV]
   REAL(8), INTENT(in) ::  mu ! chem.potential to be found [eV]
   real(8), intent(in), optional :: norm_fe ! normalization of distribution: spin resolved or not
   integer, intent(in), optional :: i_start, i_end  ! starting and ending levels to include
   real(8) :: get_E_tot ! number of electrons for given band structure and Fermi-function
   !------------------------
   integer i, Nsiz, i_low, i_high
   real(8) :: f_norm

   if (present(i_start)) then
      i_low = i_start
   else  ! default, start from 1
      i_low = 1
   endif
   Nsiz = size(wrD)
   if (present(i_end)) then
      i_high = i_end
   else  ! default, end at the end
      i_high = Nsiz
   endif


   if (Te .GT. 1.0d-12) then
      !get_E_tot = 2.0d0*SUM(wrD(:)/(1.0d0 + exp((wrD(:) - mu)/Te)))
      Nsiz = size(wrD)
      get_E_tot = 0.0d0
      DO_SUM:do i = i_low, i_high
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
      get_E_tot = f_norm*SUM(wrD(i_low:i_high), mask = wrD(:) <= mu)
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
   !---------------------
   !real(8), dimension(size(wrD)) :: exp_Fermi, exp_func
   real(8), dimension(:), allocatable :: exp_Fermi, exp_func

   allocate(exp_Fermi(size(wrD)))
   allocate(exp_func(size(wrD)))

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

   ! clean up:
   deallocate(exp_Fermi, exp_func)
end subroutine get_dN_dE_dmu



pure subroutine get_dN_dE_dTe(wrD, mu, Te, dN_dTe, dE_dTe)
   real(8), intent(out) :: dN_dTe, dE_dTe
   REAL(8), DIMENSION(:), INTENT(in) ::  wrD  ! eigenvalues of TB-Hamiltonian for electrons
   REAL(8), INTENT(inout) :: Te ! electron temperature [eV]
   REAL(8), INTENT(inout) :: mu ! chem.potential to be found [eV]
   !---------------------
   !real(8), dimension(size(wrD)) :: exp_Fermi, exp_func
   real(8), dimension(:), allocatable :: exp_Fermi, exp_func
   real(8) :: Te2

   allocate(exp_Fermi(size(wrD)))
   allocate(exp_func(size(wrD)))

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

   ! clean up:
   deallocate(exp_Fermi, exp_func)
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
      call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te and Ne

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

   if ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-5)) .GT. 1d-5) .AND. (countr .LT. 100)) then
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
   dT = max(Te/10.0d0, 10.0d0)	! temperature step [K]
   Te = max(Te/2.0d0, dT/g_kb) ! [eV] electron temperature to be calculated
   !call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te nad Ne
   !PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

   !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
   ! finding mu and Te
   countr = 0
438 cou = 0 
   Telast = Te
   Telast2 = Telast + dT/g_kb
   Te = (Telast+Telast2)/2.0d0
   call Electron_Fixed_Te(wrD, Netot, muN, Te) ! find mu for given Te nad Ne
   mu = muN ! [eV] save for output
   !print*, 'muN=', muN, Netot

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
      !print*, 'muN=', muN, Telast, Telast2, Ecur, Eetot

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

   if ((ABS((Ecur - Eetot)/MAX(ABS(Ecur),ABS(Eetot),1d-5)) .GT. 1d-5) .AND. (countr .LT. 10)) then
      countr = countr + 1
      Te = Te + 2.0d0*dT/g_kb
      goto 438 ! just try it over again if it didn't work well the first time...
   endif
end subroutine Electron_Fixed_Etot_3




END MODULE Electron_tools
