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
! This module contains subroutines to deal with cross-sections:

MODULE Monte_carlo
use Objects
use Universal_constants
use Little_subroutines, only : sample_gaussian, Find_in_array_monoton
use MC_cross_sections, only : which_shell, which_atom, NRG_transfer_elastic_atomic, Electron_energy_transfer_inelastic, &
                        Mean_free_path
use Electron_tools, only : update_cross_section, Do_relaxation_time, set_high_DOS

implicit none
PRIVATE

public :: MC_Propagate


 contains

subroutine MC_Propagate(MC, numpar, matter, Scell, laser, tim, Err) ! The entire MC procidure is here:
   type(MC_data), dimension(:), allocatable, intent(inout) :: MC    ! all MC arrays for photons, electrons and holes
   type(Numerics_param), intent(inout) :: numpar     ! all numerical parameters
   type(Solid), intent(inout) :: matter           ! all material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Pulse), dimension(:), intent(in) :: laser ! laser pulse parameters
   real(8), intent(in) :: tim ! [fs] current time step
   type(Error_handling), intent(inout) :: Err     ! errors to save
   !==============================================================
   real(8) Nphot, NMC_real, Eetot_stat, noeVB_stat, Eetot_cur, noeVB_cur, Ee_HE, Ne_high, Ne_emit, Ne_holes, N_ph, Eh
   real(8) :: E_atoms_heating, E_atoms_cur, min_df
   integer NSC, stat, i, j, N, Nph
   real(8), dimension(:,:), allocatable :: MChole
   real(8), dimension(size(Scell(1)%fe)) :: d_fe, d_fe_cur
   
   NSC = 1 ! for the moment, we have only one supercell

   if (numpar%NMC > 0) then ! if there are MC iterations at all
      call total_photons(laser, numpar, tim, Nphot) ! estimate the total number of absorbed photons
      ! Do MC run only if we have particles to trace (photons, high-energy-electrons, or core-holes):
      IF_MC:if ( (Nphot .GT. 1d-10) .OR. ( (Scell(NSC)%Ne_high-Scell(NSC)%Ne_emit) .GT. 1d-10) .OR. (Scell(NSC)%Nh .GT. 1d-10) ) then

         ! Update inelastic scattering cross section depending on Te:
         call update_cross_section(Scell(NSC), matter)  ! module "Electron_tools"

         NMC_real = dble(numpar%NMC)
         min_df = 1.0d0 / NMC_real  ! minimal allowed change of the distribution function
         Eetot_stat = 0.0d0   ! Total energy of low-energy electrons
         noeVB_stat = 0.0d0   ! Total number of low-energy electrons
         d_fe(:) = 0.0d0   ! the same as the distribution function itself
         if (.not.allocated(MChole)) then
            N = maxval(matter%Atoms(:)%sh)
            allocate(MChole(size(matter%Atoms),N)) ! different kinds of atoms
         endif
         MChole = 0.0d0
         Ee_HE = 0.0d0
         Ne_high = 0.0d0
         Ne_emit = 0.0d0
         Ne_holes = 0.0d0
         N_ph = 0.0d0
         Eh = 0.0d0
         E_atoms_heating = 0.0d0
         ! The iteration in MC are largely independent, so they can be parallelized with openmp:
!$omp parallel &
!$omp private (stat, Eetot_cur, noeVB_cur, Nph, E_atoms_cur, d_fe_cur, i)
!$omp do schedule(dynamic) reduction( + : Eetot_stat, noeVB_stat, Ee_HE, Ne_high, Ne_emit, &
!$omp                                     Eh, Ne_holes, MChole, N_ph, E_atoms_heating, d_fe)
         DO_STAT:do stat = 1,numpar%NMC ! Statistics in MC, iterate the same thing and average
            ! Set initial data:
            Eetot_cur = Scell(NSC)%nrg%El_low	! [eV] starting total energy of low-energy electrons
            noeVB_cur = Scell(NSC)%Ne_low		! number of low-energy electrons
            d_fe_cur = 0.0d0  ! change of the distribution in each iteration

            ! Perform the MC run:
            call MC_run(tim, MC(stat), Scell(NSC), laser, matter, numpar, Eetot_cur, noeVB_cur, Nph, E_atoms_cur, d_fe_cur, min_df)
            ! Add up the data from this run to total data:
            Eetot_stat = Eetot_stat + Eetot_cur  ! [eV] VB-electrons
            noeVB_stat = noeVB_stat + noeVB_cur  ! number VB-electrons
            d_fe = d_fe + d_fe_cur  ! change in distribution of low-energy electrons
            E_atoms_heating = E_atoms_heating + E_atoms_cur ! [eV] heating of atoms by electron elastic scattering
            call sort_out_electrons(MC(stat), MC(stat)%electrons, Ne_high, Ne_emit, Ee_HE) ! below
            call sort_out_holes(MC(stat), MC(stat)%holes, MChole, Ne_holes, Eh)   ! below
            N_ph = N_ph + Nph  ! number absorbed photons

         enddo DO_STAT
!$omp end do
!$omp end parallel
         ! Average all the outcome over the statistics:
         Scell(NSC)%nrg%E_glob = Scell(NSC)%nrg%E_glob + (Eetot_stat/NMC_real - Scell(NSC)%nrg%El_low)
         Scell(NSC)%nrg%El_low = Eetot_stat/NMC_real
         Scell(NSC)%Ne_low = noeVB_stat/NMC_real
         Scell(NSC)%nrg%El_high = Ee_HE/NMC_real/size(Scell(NSC)%MDatoms)
         Scell(NSC)%Ne_high = Ne_high/NMC_real
         Scell(NSC)%Ne_emit = Ne_emit/NMC_real
         Scell(NSC)%Nh = Ne_holes/NMC_real
         Scell(NSC)%nrg%Eh_tot = Eh/NMC_real/size(Scell(NSC)%MDatoms)
         Scell(NSC)%nrg%E_high_heating = E_atoms_heating/NMC_real
         do i = 1, size(matter%Atoms) ! for each kind of atoms:
            do j = 1, matter%Atoms(i)%sh ! each shell:
               Scell(NSC)%MChole(i)%Noh(j) = MChole(i,j)/NMC_real ! holes
            enddo
         enddo
         Scell(NSC)%Nph = N_ph/NMC_real
         ! Update the distribution function:
         Scell(NSC)%fe(:) = Scell(NSC)%fe(:) + d_fe(:)/NMC_real

         ! And the distribution on the grid:
         call get_high_energy_distribution(Scell(NSC), MC, numpar)  ! below

         ! Consistency checks:
         ! Make sure there are no unphysical values (fe<0 or fe>2):
         call patch_distribution(Scell(NSC)%fe, Scell(NSC)%Ei, Scell(NSC), numpar) ! below
         ! Also check that the total number of particles is conserved:
         if ( abs(Scell(NSC)%Ne_low - SUM(Scell(NSC)%fe(:))) > 1.0d-5*Scell(NSC)%Ne_low ) then
            print*, 'Error noticed in MC_E1:', Scell(NSC)%Ne_low, SUM(Scell(NSC)%fe(:))
            print*, 'Avoiding it by updating number of low-energy electron (may lose some!)'
            Scell(NSC)%Ne_low = SUM(Scell(NSC)%fe(:))
         endif
         ! And the total energy is conserved:
         if ( abs(Scell(NSC)%nrg%El_low - SUM(Scell(NSC)%fe(:) * Scell(NSC)%Ei(:))) > 1.0d-6*abs(Scell(NSC)%nrg%El_low) ) then
            print*, 'Error in MC_E2:', Scell(NSC)%nrg%El_low, SUM(Scell(NSC)%fe(:) * Scell(NSC)%Ei(:))
         endif

      else IF_MC ! if there is no more particles in MC, don't do MC at all:
         Scell(NSC)%Nph = 0.0d0 ! no photons here
         Scell(NSC)%Ne_high = Scell(NSC)%Ne_emit ! no high-energy electrons
         Scell(NSC)%Nh = 0.0d0 ! no holes
         do i = 1, size(matter%Atoms) ! for each kind of atoms:
            do j = 1, matter%Atoms(i)%sh ! and each shell:
               Scell(NSC)%MChole(i)%Noh(j) = 0.0d0
            enddo
         enddo
         Scell(NSC)%nrg%E_high_heating = 0.0d0 ! no energy transfer to atoms
         ! No changes in the electron distribution function
      endif IF_MC
      Scell(NSC)%Q = Scell(NSC)%Ne_emit/Scell(NSC)%Na ! mean unballanced charge
   else ! No MC iteratins, no changes in the system:
      Scell(NSC)%Q = 0.0d0  ! mean unballanced charge
      Scell(NSC)%nrg%E_high_heating = 0.0d0 ! no energy transfer to atoms
      ! No changes in the electron distribution function
   endif

!    print*, 'Charge:', Scell(NSC)%Q
end subroutine MC_Propagate



subroutine get_high_energy_distribution(Scell, MC, numpar)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(MC_data), dimension(:), intent(in) :: MC   ! all MC arrays for photons, electrons and holes
   type(Numerics_param), intent(inout) :: numpar   ! numerical parameters, including lists of earest neighbors
   !--------------------
   integer :: stat, i_el, j
   real(8) :: NMC

   if (numpar%save_fe_grid) then  ! only user requested
      ! Check if high-energy DOS is set, and set it if not:
      call set_high_DOS(Scell, numpar)  ! module "Electron_tools"

      ! Number of MC iterations to normalize to:
      NMC = 1.0d0/dble(numpar%NMC)

      ! Sort electrons from each iteration onto the given grid for distribution:
      do stat = 1,numpar%NMC ! Statistics in MC, iterate the same thing and average
         do i_el = 1, MC(stat)%noe  ! all active electrons
            ! Find the level, closest to where electron is incomming into:
            call Find_in_array_monoton(Scell%E_fe_grid, MC(stat)%electrons(i_el)%E, j) ! module "Little_subroutine"
            j = j - 1   ! one level below
            ! Add electron to this energy grid point:
            Scell%fe_high_on_grid(j) = Scell%fe_high_on_grid(j) + NMC ! density of electrons per iteration (without DOS)
            Scell%fe_norm_high_on_grid(j) = Scell%fe_norm_high_on_grid(j) + NMC / numpar%high_DOS(j) ! electron per iteration per DOS
         enddo ! i_el
      enddo ! stat

   endif
end subroutine get_high_energy_distribution




subroutine MC_run(tim, MC, Scell, laser, matter, numpar, Eetot_cur, noeVB_cur, Nph, E_atoms_cur, d_fe, min_df)
   real(8), intent(in) :: tim	! [fs] current time
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   type(solid), intent(in) :: matter	! materil parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), intent(inout) :: Eetot_cur, noeVB_cur	! [eV] CB electrons energy; and number
   integer, intent(out) :: Nph	! number of absorbed photons
   real(8), intent(out) :: E_atoms_cur ! [eV] heating of atoms by electrons via elastic scattering
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), intent(in) :: min_df ! minimal allowed change in the distribution function
   !========================================================
   real(8) :: RN
   integer :: i

   Nph = 0 ! just to start
   ! How many photons are absorbed in this run over time-step dt:
   call absorbed_photons_dt(laser, numpar, tim, Nph)
   E_atoms_cur = 0.0d0 ! to start with no energy transfer
   
   ! Yet again, skip those iterations where there is no particles:
   !if ( (Nph .GT. 0) .OR. (MC%noe .GT. 0) .OR. (MC%noh_tot .GT. 0) ) then
   if ( (Nph .GT. 0) .OR. ((MC%noe-MC%noe_emit) .GT. 0) .OR. (MC%noh_tot .GT. 0) ) then
      allocate(MC%photons(Nph))
      ! If we have more then one pulse, select from which we have this photon:
      call choose_photon_energy(numpar, laser, MC, tim) ! see below
      ! MC modeling of photoabsorbtion:
      call MC_for_photon(tim, MC, Nph, numpar, matter, laser, Scell, Eetot_cur, noeVB_cur, d_fe, min_df) ! below

      ! MC modeling of electron propagation:
      call MC_for_electron(tim, MC, matter, numpar, Scell, Eetot_cur, noeVB_cur, E_atoms_cur, d_fe, min_df) ! below
      
      ! MC modeling of hole decay:
      call MC_for_hole(tim, MC, matter, numpar, Scell, Eetot_cur, noeVB_cur, d_fe, min_df) ! below
      
      ! At the end, clean up:
      deallocate(MC%photons)	! free it for the next time-step
   endif
!    print*, 'MC_run finished'
end subroutine MC_run



!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! For electrons:
subroutine MC_for_electron(tim, MC, matter, numpar, Scell, Eetot_cur, noeVB_cur, E_atoms_cur, d_fe, min_df)
   real(8), intent(in) :: tim	! [fs] current time
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(solid), intent(in) :: matter	! materil parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   real(8), intent(inout) :: Eetot_cur, noeVB_cur	! [eV] CB electrons energy; and number
   real(8), intent(out) :: E_atoms_cur ! [eV] heating of atoms by electrons via elastic scattering
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), intent(in) :: min_df
   !=============================================
   real(8) hw	! transferred energy [eV]
   real(8) kind_of_coll ! elasctic vs inelastic
   real(8) RN, Ekin, IMFP, EMFP
   integer j, shl, KOA
   j = 1
   do while ((MC%noe - j) .GE. 0)  ! all electrons
      do while (MC%electrons(j)%ti .LT. tim) ! until time of this electron becomes larger than the current timestep
         Ekin = MC%electrons(j)%E - Scell%E_bottom ! [eV] kinetic energy of an electron is counted from the bottom of CB
         call random_number(RN)
         ! Get MFPs for inelastic and elastic scattering to compare:
         call Mean_free_path(Ekin, matter%El_MFP_tot, IMFP, inversed=.true.) ! inelastic MFP [1/A], module "MC_cross_section"
         call Mean_free_path(Ekin, matter%El_EMFP_tot, EMFP, inversed=.true.) ! elastic MFP [1/A], module "MC_cross_section"
         
         elast_vs_inelast: if (RN <= IMFP/(IMFP + EMFP)) then ! it is inelastic scattering
!              print*, 'Inelastic scattering event:', RN, j, Ekin, IMFP/(IMFP + EMFP), IMFP, EMFP
            
            ! Which shell can be impact-ionized:
            !call which_shell(MC%electrons(j)%E, matter%Atoms, matter, 1, KOA, SHL) ! module 'MC_cross_sections'
            call which_shell(Ekin, matter%Atoms, matter, 1, KOA, SHL) ! module 'MC_cross_sections'

!             if (MC%electrons(j)%E .LT. matter%Atoms(KOA)%Ip(SHL)) then
            if (Ekin .LT. matter%Atoms(KOA)%Ip(SHL)) then
               print*, 'Attention! In subroutine MC_for_electron:'
               print*, 'Electron energy is lower than the ionization potential!'
               write(*,'(f25.16, f25.16, i2, i2)') MC%electrons(j)%E, matter%Atoms(KOA)%Ip(SHL), KOA, SHL
!                write(*,'(f25.16,$)') matter%Atoms(KOA)%Ip(SHL)
!                write(*,'(a)') ''
!                pause 'MC_for_electron'
            endif

            ! how much energy it loses:
            !call Electron_energy_transfer_inelastic(matter, MC%electrons(j)%E, KOA, SHL, matter%Atoms(KOA)%El_MFP(shl), hw) ! module "Cross_sections"
            call Electron_energy_transfer_inelastic(matter, Scell%TeeV, Ekin, KOA, SHL, matter%Atoms(KOA)%El_MFP(shl), hw) ! module "Cross_sections"

            ! new electron (and may be hole) is created:
            call New_born_electron_n_hole(MC, KOA, shl, numpar, Scell, matter, hw, MC%electrons(j)%ti, &
                                          Eetot_cur, noeVB_cur, d_fe, min_df)
         else  elast_vs_inelast ! it is elastic scattering
            call which_atom(Ekin, matter%Atoms, EMFP, KOA) ! get which kind of atoms we scatter on, module "MC_cross_sections"
            call NRG_transfer_elastic_atomic(matter%Atoms(KOA)%Ma, matter%Atoms(KOA)%Z, Ekin, hw) ! module "MC_cross_sections"
            E_atoms_cur = E_atoms_cur + hw ! [eV] energy transferred to atoms by high-energy electrons' elastic scattering
         endif elast_vs_inelast

         ! new parameters for the incident electrons:
         call update_electron_data(MC, Scell, matter, numpar, j, hw, Eetot_cur, noeVB_cur, d_fe, min_df) ! below
      enddo ! time
      j = j + 1 ! next electron
   enddo ! number 
end subroutine MC_for_electron



subroutine update_electron_data(MC, Scell, matter, numpar, j, hw, Eetot_cur, noeVB_cur, d_fe, min_df)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(solid), intent(in) :: matter	! materil parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
!    type(MFP), intent(in) :: mfps	! total electron mean free paths
   integer, intent(in) :: j	! number of electron in the aray
   real(8), intent(in) :: hw	! energy transferred to the electron
   real(8), intent(inout) :: Eetot_cur, noeVB_cur ! [eV] CB electrons energy; and number
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), intent(in) :: min_df ! minimal allowed change in the distribution function (from MC)
   !-----------------------------------------------
   logical :: emitted
   emitted = .false.
   MC%electrons(j)%colls = MC%electrons(j)%colls + 1 ! new collision occured for this electron
   MC%electrons(j)%E = MC%electrons(j)%E - hw	! it lost this much energy
   
   ! Check if the electron is emitted out of the material:
   if (numpar%E_work <= 0.0d0) then ! emission after a number of collisions:
      ! exclude electrons that made more collisions than allowed
      ! also account for the accumulating positive charge inside of the material
      if ((MC%electrons(j)%colls > ABS(numpar%E_work)) .and. (MC%electrons(j)%E > numpar%E_Coulomb)) then
         call emit_electron(MC, j) ! see below
         emitted = .true.
      endif
   else ! emission in case if energy is above 'work-function':
      ! also account for the accumulating positive charge inside of the material
      if ((MC%electrons(j)%E > numpar%E_work) .and. (MC%electrons(j)%E > numpar%E_Coulomb)) then 
         ! exclude electrons with the energy above the work functinon:
         call emit_electron(MC, j) ! see below
         emitted = .true.
      endif
   endif

   if (.not.emitted) then
      ! Find where this electron goes: high-energy part or low-energy part:
      if (MC%electrons(j)%E .LT. (Scell%E_bottom + numpar%E_cut)) then ! Electron joins VB:
         ! Corresponding change in the distribution function:
         call d_distribution_between_levels(d_fe, Scell%Ei, MC%electrons(j)%E, Scell%E_bottom, Scell, min_df)  ! below
         Eetot_cur = Eetot_cur + MC%electrons(j)%E
         noeVB_cur = noeVB_cur + 1
         call electron_disappears(MC, j) ! one electron less
      else ! normal electron propagation:
         call get_electron_MFP(MC, matter%El_MFP_tot, matter%El_EMFP_tot, j, MC%electrons(j)%E-Scell%E_bottom) ! next scattering event
      endif
   endif
end subroutine update_electron_data



subroutine emit_electron(MC, j)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   integer, intent(in) :: j	! number of electron in the aray
   MC%electrons(j)%ti = 1d33 ! [fs] it will never return back and make the next scattering event
   MC%noe_emit = MC%noe_emit + 1 ! one more emitted electron
!    print*, 'An electron has been emitted!',  MC%electrons(j)%ti, MC%noe_emit
end subroutine emit_electron


subroutine electron_disappears(MC, j)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   integer, intent(in) :: j	! number of electron in the aray
   integer k
   do k = j, MC%noe-1 ! one electron from Hight Energy states disappeared 
      MC%electrons(k)%E = MC%electrons(k+1)%E	! [eV]
      MC%electrons(k)%ti = MC%electrons(k+1)%ti	! [fs]
   enddo ! k
   MC%electrons(MC%noe)%E = 0.0d0	! [eV]
   MC%electrons(MC%noe)%ti = 1d20	! [fs]
   MC%noe = MC%noe - 1 ! electron has fallen into the sea of CB
end subroutine electron_disappears


subroutine New_born_electron_n_hole(MC, KOA, SHL, numpar, Scell, matter, hw, t_cur, Eetot_cur, noeVB_cur, d_fe, min_df)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   integer, INTENT(in) ::  KOA, SHL   ! number of atom and its shell that is being ionized
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(solid), intent(in) :: matter	! materil parameters
   real(8), intent(in) :: hw	! energy transferred to the electron
   real(8), intent(in) :: t_cur ! [fs] time of electron creation
   real(8), intent(inout) :: Eetot_cur, noeVB_cur	! [eV] CB electrons energy; and number
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), intent(in) :: min_df ! minimal allowed change in the distribution function
   !--------------------------------------
   REAL(8) Ee ! electron energy after ionization [eV]
   real(8) RN, norm_sum, t_dec
   integer i, IONIZ ! which level is ionized

   ! Find from which shell an electron is excited:
   if ((KOA .EQ. 1) .and. (SHL .EQ. matter%Atoms(KOA)%sh)) then ! VB:
      !call sample_VB_level(Scell%Ne_low, Scell%fe, i, wr=Scell%Ei, Ee=hw, min_df=min_df) ! NEW
      ! Include the influence of dynamically changing distribution:
      call sample_VB_level(Scell%Ne_low, (Scell%fe+d_fe*min_df), i, wr=Scell%Ei, Ee=hw, min_df=min_df) ! NEW
      if (i == 0) i = Scell%N_Egap ! very top of VB
      Ee = hw + Scell%Ei(i) ! [eV] electron energy
      IONIZ = i ! from this level
      noeVB_cur = noeVB_cur - 1 ! one electron has left VB going up
      Eetot_cur = Eetot_cur - Scell%Ei(i) ! and brought energy with it
      d_fe(i) = d_fe(i) - 1.0d0  ! change in the electron distribution
      !if ((Scell%fe(i)+d_fe(i)*min_df) < 0.0d0) then
      if ((Scell%fe(i)+d_fe(i)*min_df) < 0.0d0) then
         print*, 'Error New_born_electron_n_hole:', i, Scell%fe(i), d_fe(i)*min_df, Ee, Scell%E_gap, numpar%E_cut
      endif
!       if (i > 129) then ! testing
!          print*, 'Potential New_born_electron_n_hole:', i, Scell%fe(i), d_fe(i)*min_df
!       endif
   else ! deep shell:
      Ee = hw - matter%Atoms(KOA)%Ip(SHL) ! [eV] electron energy
      IONIZ = 0	! deep shell, which one is determined by "SHL"
      ! a new hole is created:
      MC%noh_tot = MC%noh_tot + 1    ! one more hole is created
      MC%holes(MC%noh_tot)%KOA = KOA ! this atom is excited
      MC%holes(MC%noh_tot)%Sh = shl  ! this shell is excited
      MC%holes(MC%noh_tot)%E = matter%Atoms(KOA)%Ip(shl)	  ! [eV] energy of this hole
      call Hole_decay_time_sampled(matter%Atoms, KOA, shl, t_dec) ! below
      MC%holes(MC%noh_tot)%ti = t_cur + t_dec ! [fs] decay time
      if (Ee < 0.0d0) then
         print*, 'Subroutine New_born_electron_n_hole', hw
         print*, 'Produced negative electron energy:', Ee, IONIZ
         print*, 'deep shell:', matter%Atoms(KOA)%Ip(SHL), KOA, SHL
         pause
      endif
   endif

   ! Find where this electron goes: high-energy part or low-energy part
   if (Ee > (Scell%E_bottom + numpar%E_cut)) then ! high-energy, MC part
      MC%noe =  MC%noe + 1 	! new electron was born!
      if (MC%noe > size(MC%electrons)) then
         call extend_MC_array(MC%electrons)  ! below
      endif
      MC%electrons(MC%noe)%E = Ee ! [eV] with this energy
      MC%electrons(MC%noe)%ti = t_cur ! [eV] with this energy
      ! Get the new electron's next-scattering time:
      call get_electron_MFP(MC, matter%El_MFP_tot, matter%El_EMFP_tot, MC%noe, MC%electrons(MC%noe)%E-Scell%E_bottom) ! below
   else	! Electron joins VB:
      Eetot_cur = Eetot_cur + Ee
      noeVB_cur = noeVB_cur + 1
!       print*, 'New_born_electron_n_hole:', Ee, (Scell%E_bottom + numpar%E_cut)
      call d_distribution_between_levels(d_fe, Scell%Ei, Ee, Scell%E_bottom, Scell, min_df)  ! below
   endif
end subroutine New_born_electron_n_hole


subroutine d_distribution_between_levels(d_fe, Ei, Ee, E_bottom, Scell, min_df, delta_f)
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), dimension(:), intent(in) :: Ei ! energy levels [eV]
   real(8), intent(in) :: Ee  ! [eV] incoming electron's energy
   real(8), intent(in) :: E_bottom  ! [eV] bottom of conduction band, where Ee is counted from
   type(Super_cell), intent(in) :: Scell ! supercell with all the atoms as one object
   real(8), intent(in) :: min_df ! minimal allowed change in the distribution function
   real(8), intent(in), optional :: delta_f  ! given change of the distribution (number of particle)
   !--------------
   integer :: j
   real(8) :: eps, E_abs, dE, df

   if (present(delta_f)) then ! given change in the number of electrons
      df = delta_f
   else  ! assume one electron
      df = 1.0d0
   endif

   eps = 1.0d-8  ! precision

   ! Find the level, closest to where electron is incomming into:
   call Find_in_array_monoton(Ei, Ee, j) ! module "Little_subroutine"
   j = j - 1   ! one level below

   if (j >= size(Ei) .or. (j <= 1)) print*, 'd_distribution_between_levels trouble', j, Ee, Ei(size(Ei))

   dE = Ei(j+1)-Ei(j)   ! energy levels difference

   ! change in the electron distribution:
   if (dE < eps) then   ! degenerate levels
      d_fe(j) = d_fe(j) + df  !1.0d0
   else  ! different levels
      ! Fractions of electron distributed between two closest levels,
      ! ensuring conservation of particles and energy:
      d_fe(j)   = d_fe(j) + (Ei(j+1) - Ee)/dE * df
      d_fe(j+1) = d_fe(j+1) + (Ee - Ei(j))/dE * df
   endif

   if ( ((Scell%fe(j)+d_fe(j)*min_df) > 2.0d0) .or. (Scell%fe(j+1)+d_fe(j+1)*min_df) > 2.0d0 ) then
      print*, 'Possible trouble d_distribution_between_levels:', j, Scell%fe(j), d_fe(j)*min_df, Scell%Ei(j), Scell%fe(j+1), d_fe(j+1)*min_df, Scell%Ei(j+1)
   endif

end subroutine d_distribution_between_levels



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
         if (fe(i) > 2.0d0+eps) then   ! it's within [2; 2+eps]
            print*, 'Problem in patch_distribution #3a:', i, fe(i)
            fe(i) = 2.0d0
            trouble_present = .true.   ! there still is a problem
         elseif (fe(i) > 2.0d0) then   ! it's within [2; 2+eps]
            fe(i) = 2.0d0        ! distribution adjusted to accceptable
         elseif (fe(i) < -eps) then  ! it's within [0-eps;0]
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



subroutine get_electron_MFP(MC, Imfps, Emfps, j, Ekin) ! calculate total electron mean free path and scattering time
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(MFP), intent(in) :: Imfps	! total inelastic mean free path
   type(MFP), intent(in) :: Emfps	! total elastic mean free path
   integer, INTENT(in) :: j   ! number of electron
   real(8), INTENT(in) :: Ekin ! [eV] electron kinetic energy
   !-----------------
   real(8) RN, MFPath, IMFPath, EMFPath
   call random_number(RN)
   ! find mean free paths from the array given [A]:
   !call Mean_free_path(MC%electrons(j)%E, mfps, MFPath) !module "MC_cross_section"
   call Mean_free_path(Ekin, Imfps, IMFPath, inversed=.true.) ! inelastic MFP [1/A], module "MC_cross_section"
   call Mean_free_path(Ekin, Emfps, EMFPath, inversed=.true.) ! elastic MFP [1/A], module "MC_cross_section"
   MFPath = 1.0d0/(IMFPath+EMFPath) ! total mean free path [A]
   
   !MC%electrons(j)%ti = MC%electrons(j)%ti - log(RN)*(MFPath*1d-10)/sqrt(2.0d0*MC%electrons(j)%E*g_e/g_me)*1d15 ! [fs]
   MC%electrons(j)%ti = MC%electrons(j)%ti - log(RN)*(MFPath*1d-10)/sqrt(2.0d0*Ekin*g_e/g_me)*1d15 ! [fs]
!    print*, 'Ti:', MC%electrons(j)%E, MC%electrons(j)%ti, MFPath, (MFPath*1d-10)/sqrt(2.0d0*MC%electrons(j)%E*g_e/g_me)*1d15
end subroutine get_electron_MFP



subroutine sample_VB_level(Ne_low, fe, i, wr, Ee, min_df)
   real(8), intent(in) :: Ne_low ! number of low-energy electrons
   real(8), dimension(:), INTENT(in) :: fe  ! electron distribution function [eV, number]
   integer, intent(out) :: i  ! sampled VB level, from where electron is emitted
   REAL(8), DIMENSION(:), INTENT(in), optional ::  wr	! [eV] energy levels
   real(8), intent(in), optional :: Ee ! energy transfer [eV]
   real(8), intent(in), optional :: min_df   ! minimal allowed change in the distribution function
   !===================================
   real(8) :: RN, norm_sum, cur_sum, sampled_sum, min_df_used, eps
   real(8), dimension(size(fe)) :: E_weight
   real(8), dimension(size(fe)) :: fe_final ! construct an array of final states populations
   integer j, N_levels, j_fin
   logical :: close_levels

   eps = 1.0d0 ! [eV] acceptence interval, an electron can come into in between the levels

   ! Set the minimal allowed change in the distribution function:
   if (present(min_df)) then
      min_df_used = min_df ! only the change defined by MC iterations that does not make unphysical values in fe
   else
      min_df_used = 0.0d0  ! assume any change is allowed
   endif

   call random_number(RN) ! sample random number
   norm_sum = 0.0d0
   PART_VB:if ((present(Ee)) .and. (present(wr))) then
      N_levels = size(wr) ! total number of electron energy levels
      ! Construct probabilities of scattering events for all the levels:
      do j = 1, N_levels
         ! find what is the final state number:
         if (wr(j)+Ee > wr(N_levels)) then ! transferred energy brings electron out of low-energy domain:
            fe_final(j) = 0.0d0 ! these levels are totally free in our model
         else ! final state is within low-energy domain, find it:
            ! SIDENOTE: it is the closest level existing, but it is not exactly equal to (wr(j)+Ee)!
            call Find_in_array_monoton(wr, wr(j)+Ee, j_fin) ! module "Little_subroutine"
            !print*, 'Test:', wr(j)+Ee, wr(j_fin-1), wr(j_fin)
            j_fin = j_fin - 1 ! one level below
            fe_final(j) = max(fe(j_fin),fe(j_fin+1)) ! that's the transient population on the final level
            ! Check that levels are not >too far apart:
            if ((wr(j_fin+1)-wr(j_fin)) > eps) then
               fe_final(j) = 2.0d0  ! exclude levels that are too far apart
            endif
         endif
         
!          ! Weight is calculated according to the probability of ionization of a level,
!          ! assumming 1/E dependence of the ionization cross section:
!          if (abs(wr(j)) <= 1d-2) then
!             E_weight(j) = 1d2
!          else
!             E_weight(j) = abs(1.0d0/wr(j))
!          endif
         E_weight(j) = 1.0d0 ! Exclude weights
         
         ! that is the relative probability of the scattering event, accounting for 
         ! number of electrons in the inital level an electron can scatter off, and
         ! the number of free places in the final state an electron can come to:
         if ( (fe(j) > min_df_used) .and. (2.0d0-fe_final(j) > min_df_used) ) then ! only if allowed
            norm_sum = norm_sum + fe(j)*(2.0d0 - fe_final(j))*E_weight(j)
         endif
         !if (Ee < 6.5) print*, 'sample_VB_level #0:', Ee, j, fe(j), fe_final(j), wr(j), wr(j+1), norm_sum
         !write(*,'(i4,f,f,f,f,f)') j, Ee, wr(j), fe(j), fe_final(j), norm_sum
      enddo

      if (norm_sum <= 0.0d0) then ! no ionization of VB is possible for some reason...
         i = 0
         !print*, 'Error sample_VB_level #1:', i, norm_sum, Ee
         !print*, 'sample_VB_level', Ee, j_fin, fe(j_fin), fe(j_fin+1), wr(j_fin), wr(j_fin+1)
      else ! something is possible:
         cur_sum = 0.0d0
         i = 0 ! start from bottom of VB
         sampled_sum = RN*norm_sum ! this is the sampled value of the probability we need to reach
         do while ((cur_sum - sampled_sum) < 0.0d0) ! find which scattering event it is
            i = i+1 ! next level
            if (i > N_levels) exit ! the final state is reached
            !if ( (fe(i) > min_df_used) .and. (2.0d0-fe_final(i) > min_df_used) .and. ((wr(i+1)-wr(i)) < eps)) then ! only if allowed
            if ( (fe(i) > min_df_used) .and. (2.0d0-fe_final(i) > min_df_used) ) then ! only if allowed
               cur_sum = cur_sum + fe(i)*(2.0d0 - fe_final(i))*E_weight(i) ! go through the probabilities until reach the given value
            endif
!             if (fe(i) < min_df_used) then  ! testing
!                write(*,'(a,i4,f,f,f,f,f)') 'sample_VB_level', i, wr(i), fe(i), fe_final(i), sampled_sum, cur_sum, RN
!             endif
         enddo ! while
         !if (i == 0) print*, 'Error sample_VB_level #2:', i, norm_sum, cur_sum

      endif
      
   else PART_VB ! Full VB is possible to ionize:
      i = 0
      sampled_sum = RN*Ne_low ! this is the sampled value of the probability we need to reach
      do while ((norm_sum - sampled_sum) <= 0.0d0) ! find which electron it is among all of them
         i = i+1
         norm_sum = norm_sum + fe(i)
      enddo ! while
      if (i == 0) print*, 'Error sample_VB_level #3:', i, norm_sum, sampled_sum
   endif PART_VB
end subroutine sample_VB_level


subroutine sample_VB_level_OLD(Ne_low, fe, i, wr, Ee)
   real(8), intent(in) :: Ne_low ! number of low-energy electrons
   real(8), dimension(:), INTENT(in) :: fe  ! electron distribution function [eV, number]
   integer, intent(out) :: i
   REAL(8), DIMENSION(:), INTENT(in), optional ::  wr	! [eV] energy levels
   real(8), intent(in), optional :: Ee ! energy transfer [eV]
   !===================================
   real(8) RN, norm_sum, cur_sum
   integer j
   call random_number(RN)
   norm_sum = 0.0d0
   PART_VB:if ((present(Ee)) .and. (present(wr))) then
      j = size(wr) ! start from top
      do while ((wr(j) + Ee) >= 0.0d0) ! find how many levels are possible to ionize with this energy
         norm_sum = norm_sum + fe(j) ! that's how many electrons are on this shell, add it
         j = j - 1 ! go deeper in the VB
         if (j == 0) exit
      enddo
      if (norm_sum <= 0.0d0) then ! no ionization of VB is possible!
         i = 0
      else ! something is possible:
         cur_sum = 0.0d0
         i = size(wr) ! start from top
         do while ((cur_sum - RN*norm_sum) <= 0.0d0) ! find which electron it is
            i = i-1
            cur_sum = cur_sum + fe(i)
         enddo ! while
      endif
   else PART_VB ! Full VB is possible to ionize:
      i = 0
      do while ((norm_sum - RN*Ne_low) <= 0.0d0) ! find which electron it is among all of them
         i = i+1
         norm_sum = norm_sum + fe(i)
      enddo ! while
   endif PART_VB
end subroutine sample_VB_level_OLD

!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE




!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
! MC for Holes:
!subroutine hole_decay(stat, MC, matter, numpar, mfps, tim, fe, Ei, Eetot_cur, noeVB_cur)
subroutine MC_for_hole(tim, MC, matter, numpar, Scell, Eetot_cur, noeVB_cur, d_fe, min_df)
   real(8), intent(in) :: tim	! current timestep
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(solid), intent(in) :: matter	! materil parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   real(8), intent(inout) :: Eetot_cur, noeVB_cur	! [eV] CB electrons energy; and number
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), intent(in) :: min_df
   !=============================================
   real(8) t_cur, hw, cur
   integer i, j, KOA, shl, Nsh1, Nsh2, N_val1, KOA1, KOA2
   j = 1
   do while ((MC%noh_tot - j) .GE. 0)  ! all holes
      KOA = MC%holes(j)%KOA     ! this atom has a hole
      shl = MC%holes(j)%Sh	! this hole is in this shell
      ! Trace until time of this hole becomes larger than the current timestep (a few decays may be possible):
      do while (MC%holes(j)%ti .LT. tim)
         ! Find to which shell the first hole pops after Auger:
         call which_shell_Auger(matter, MC, Scell, KOA, shl, KOA1, Nsh1, N_val1, hw, min_df)   ! below
         t_cur = MC%holes(j)%ti	! [fs] time of hole deacy, save for future
         ! First (old) hole:
         call update_this_hole(KOA1, Nsh1, N_val1, Scell%Ei, matter, MC, j, Eetot_cur, noeVB_cur, d_fe)   ! below
         ! Second hole -- for simplicity, assume it occurs via virtual photon:
         call which_shell(hw, matter%Atoms, matter, 0, KOA2, Nsh2) ! module 'MC_cross_sections'
         ! Second (new) electron and hole:
         call New_born_electron_n_hole(MC, KOA2, Nsh2, numpar, Scell, matter, hw, t_cur, Eetot_cur, noeVB_cur, d_fe, min_df)  ! below
         KOA = MC%holes(j)%KOA     ! this atom has a hole here now
         shl = MC%holes(j)%Sh      ! this hole is in this shell now
      enddo	! time
      j = j + 1 ! next hole
   enddo	! all holes
end subroutine MC_for_hole



subroutine update_this_hole(KOA, Nsh, N_val, Ei, matter, MC, j, Eetot_cur, noeVB_cur, d_fe)
   integer, intent(in) :: KOA, Nsh, N_val ! shells that participate in Auger-decay (final states), states in VB if it's VB
   real(8), dimension(:), intent(in) :: Ei	! [eV] energy level, eigenvalues of H_TB
   type(solid), intent(in) :: matter	! materil parameters
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   integer, intent(in) :: j	! number of hole
   real(8), intent(inout) :: Eetot_cur, noeVB_cur	! [eV] CB electrons energy; and number
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   !-----------------------
   real(8) t_dec

   if ((KOA == 1) .and. (Nsh == matter%Atoms(1)%sh)) then ! VB:
      Eetot_cur = Eetot_cur - Ei(N_val)   ! the hole is in VB now, at this level
      noeVB_cur = noeVB_cur - 1           ! meaning, an electron from here disappeared
      call hole_disappears(MC, j) ! VB holes are considered in other domain, not MC
      d_fe(N_val) = d_fe(N_val) - 1.0d0   ! corresponding change in the electron distribution
   else  ! core shell
      MC%holes(j)%E = matter%Atoms(KOA)%Ip(Nsh)	! the hole is in this shell now
      call Hole_decay_time_sampled(matter%Atoms, KOA, Nsh, t_dec) ! below
      MC%holes(j)%ti = MC%holes(j)%ti + t_dec ! [fs] next decay time
      MC%holes(j)%KOA = KOA	! it moved to this atom
      MC%holes(j)%Sh = Nsh	! it moved to this shell
   endif
end subroutine update_this_hole


subroutine hole_disappears(MC, j)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   integer, intent(in) :: j	! this hole disappears
   integer k
   do k = j, MC%noh_tot-1 ! one hole disappeared
      MC%holes(k)%E = MC%holes(k+1)%E	! [eV]
      MC%holes(k)%ti = MC%holes(k+1)%ti	! [fs]
      MC%holes(k)%Sh = MC%holes(k+1)%Sh ! this shell
      MC%holes(k)%KOA = MC%holes(k+1)%KOA ! this atom
   enddo ! k
   MC%holes(MC%noh_tot)%E = 0.0d0	! [eV]
   MC%holes(MC%noh_tot)%ti = 1d20	! [fs]
   MC%holes(MC%noh_tot)%Sh = 0
   MC%holes(MC%noh_tot)%KOA = 0
   MC%noh_tot = MC%noh_tot - 1	! holes disappears
end subroutine hole_disappears


subroutine which_shell_Auger(matter, MC, Scell, KOA, shl, KOA1, Nsh1, N_val1, hw, min_df)
   type(solid), intent(in) :: matter	! materil parameters
   type(MC_data), intent(in) :: MC	! all MC arrays for photons, electrons and holes
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: KOA, shl	! shell atom and which decays
   integer, intent(out) :: KOA1, Nsh1 ! first atom and shell participating in Auger (final state of the hole)
   integer, intent(out) :: N_val1 ! states in VB if it's VB
   real(8), intent(out) :: hw	! energy given to the ionized electron in this Auger [eV]
   real(8), intent(in) :: min_df
   !----------------
   real(8) RN, Ecur
   integer i, Ne
   N_val1 = 0
   ! For simplicity, assume it's always the next one of the same atom (Nsh1 = shl + 1):
   if (((KOA == 1) .and. (shl == matter%Atoms(1)%sh-1)) .or. (shl == matter%Atoms(KOA)%sh)) then ! Next is VB:
      KOA1 = 1
      Nsh1 = matter%Atoms(1)%sh
      call sample_VB_level(Scell%Ne_low, Scell%fe, N_val1, wr=Scell%Ei, Ee=100.0d0, min_df=min_df)  ! above
      if (N_val1 == 0) N_val1 = Scell%N_Egap ! very top of VB
      hw = matter%Atoms(KOA)%Ip(shl) + Scell%Ei(N_val1)
   else ! it can be atomic shell, not only VB:
      KOA1 = KOA     ! assume the same atom for simplicity
      Nsh1 = shl + 1 ! try next shell
      N_val1 = 0      
      hw = matter%Atoms(KOA)%Ip(shl) - matter%Atoms(KOA1)%Ip(Nsh1)
      SHL_CHECK:do while (matter%Atoms(KOA)%Ip(shl) - matter%Atoms(KOA1)%Ip(Nsh1) < Scell%E_gap) ! if shells are too close, no Auger is possible
         Nsh1 = Nsh1 + 1 ! try next shell
         if ((KOA1 == 1) .and. (Nsh1 >= size(matter%Atoms(KOA1)%Ip))) then ! it's VB:
            Nsh1 = matter%Atoms(1)%sh
            call sample_VB_level(Scell%Ne_low, Scell%fe, N_val1, wr=Scell%Ei, Ee=100.0d0, min_df=min_df)  ! above
            if (N_val1 == 0) N_val1 = Scell%N_Egap ! very top of VB
            hw = matter%Atoms(KOA)%Ip(shl) + Scell%Ei(N_val1)
            exit SHL_CHECK
         elseif (Nsh1 > size(matter%Atoms(KOA1)%Ip)) then ! other atoms don't have VB in this description, so make it:
            KOA1 = 1
            Nsh1 = matter%Atoms(1)%sh
            call sample_VB_level(Scell%Ne_low, Scell%fe, N_val1, wr=Scell%Ei, Ee=100.0d0, min_df=min_df)  ! above
            if (N_val1 == 0) N_val1 = Scell%N_Egap ! very top of VB
            hw = matter%Atoms(KOA)%Ip(shl) + Scell%Ei(N_val1)
            exit SHL_CHECK
         else ! it is still within the array:
            hw = matter%Atoms(KOA)%Ip(shl) - matter%Atoms(KOA1)%Ip(Nsh1)
         endif
      enddo SHL_CHECK
   endif
end subroutine which_shell_Auger


subroutine Hole_decay_time_sampled(Atoms, KOA, shl, t_dec)
   type(At_data), dimension(:), intent(in) :: Atoms	! all kinds of atoms of the compound
   integer, INTENT(in) :: KOA, shl	! atoms and shell
   REAL(8), INTENT(out) :: t_dec	! hole decay
   real(8) RN
   call random_number(RN)
   t_dec = -log(RN)*Atoms(KOA)%Auger(shl)  ! [fs]
end subroutine Hole_decay_time_sampled

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH



!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! MC for photons:
! Photon absorption and corresponding creation of electron and hole:
subroutine MC_for_photon(tim, MC, Nph, numpar, matter, laser, Scell, Eetot_cur, noeVB_cur, d_fe, min_df)
   real(8), intent(in) :: tim	! [fs] current timestep
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   integer, intent(in) :: Nph	! number of absorbed photons
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(solid), intent(in) :: matter	! materil parameters
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   real(8), intent(inout) :: Eetot_cur, noeVB_cur	! [eV] CB electrons energy; and number
   real(8), dimension(:), intent(inout) :: d_fe ! change in the electron distribution function
   real(8), intent(in) :: min_df ! minimal allowed change in the distribution function
   !=============================================
   real(8) RN, t_cur
   integer KOA, SHL, j, PN
   ! If some photons were absorbed at this iteration at this timestep
   if (Nph .GT. 0) then
      do j = 1, Nph ! all photons absorbed during this timestep
         call random_number(RN)
         MC%photons(j)%ti = tim + RN*numpar%dt ! [fs] time of photoabsorbtion within this timestep
         t_cur = MC%photons(j)%ti	! [fs] time of photoabsorbtion
         !! find number of pulse to which this photon belongs:
         !call Find_in_1D_array(laser(:)%hw, MC%photons(j)%E, PN) ! module 'Little_subroutines'
         !call Photon_which_shell(PN, matter, shl) ! finds shell from which electron is ionized
         ! Find by which shell this photon is absorbed:
         call which_shell(MC%photons(j)%E, matter%Atoms, matter, 0, KOA, SHL) ! module 'MC_cross_sections'
         ! Photon is absorbed, electron-hole pair is created:
         call New_born_electron_n_hole(MC, KOA, SHL, numpar, Scell, matter, MC%photons(j)%E, t_cur, &
                                       Eetot_cur, noeVB_cur, d_fe, min_df) ! above
       enddo ! j
   endif
end subroutine MC_for_photon
!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP


!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
subroutine sort_out_electrons(MC, electrons, Ne_high, Ne_emit, Ee_HE)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(Electron), dimension(:), intent(in) :: electrons ! all electrons in MC
   real(8), intent(inout) :: Ee_HE, Ne_high, Ne_emit ! number and energy of high-energy electrons, and number of emitted electrons
   integer i
!    do i = 1, MC%noe
!       Ee_HE = Ee_HE + electrons(i)%E ! [eV] high-energy electrons
!    enddo
   Ee_HE = Ee_HE + SUM(electrons(1:MC%noe)%E) ! [eV] high-energy electrons
   Ne_high = Ne_high + MC%noe ! number high-energy electrons
   Ne_emit = Ne_emit + MC%noe_emit ! number emitted electrons
end subroutine sort_out_electrons


subroutine sort_out_holes(MC, holes, MChole, Ne_holes, Eh)
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   type(Hole), dimension(:), intent(in), target :: holes ! all holes in MC
   real(8), dimension(:,:), intent(inout) :: MChole ! array of holes in each shell of eahc atom
   real(8), intent(inout) :: Ne_holes, Eh ! number and energy of holes
   integer i
   integer, pointer :: KOA, SHL ! kind of atom; number of shell
   ! Total energy of all holes:
   !Eh = Eh + SUM(holes(:)%E, MASK=(holes(:)%E .GT. 1.0d-10))	! [eV] holes
   ! Total number of holes:
   !Ne_holes = Ne_holes + COUNT(holes(:)%E, MASK=(holes(:)%E .GT. 1.0d-10)) ! number core holes
   ! Holes in different shells:
   do i = 1, MC%noh_tot
      KOA => holes(i)%KOA
      SHL => holes(i)%Sh
      MChole(KOA,SHL) = MChole(KOA,SHL) + 1 ! holes in this shell
   enddo
   ! Sum them up from all iterations:
   Eh = Eh + SUM(holes(1:MC%noh_tot)%E)	! [eV] holes
   Ne_holes = Ne_holes + MC%noh_tot ! number of holes
   nullify(KOA, SHL)
end subroutine sort_out_holes


subroutine extend_MC_array(array1)
   type(Electron), dimension(:), allocatable, intent(inout) :: array1
   integer N
   integer, dimension(:), allocatable :: arrayE, arrayti, arrayt
   N = size(array1)
   allocate(arrayE(N))
   allocate(arrayti(N))
   allocate(arrayt(N))
   arrayE = array1%E
   arrayti = array1%ti
   arrayt = array1%t
   deallocate(array1)
   allocate(array1(2*N))
   array1(1:N)%E = arrayE(1:N)
   array1(1:N)%ti = arrayti(1:N)
   array1(1:N)%t = arrayt(1:N)
   deallocate(arrayE)
   deallocate(arrayti)
   deallocate(arrayt)
end subroutine extend_MC_array


!LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
! Laser pulse and photons:
subroutine choose_photon_energy(numpar, laser, MC, tim) ! see below
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   type(MC_data), intent(inout) :: MC	! all MC arrays for photons, electrons and holes
   real(8), intent(in) :: tim  ! current time-step [fs]
   real(8), dimension(:), allocatable :: Nphot_pulses
   real(8) RN, Tot, sum_phot
   integer coun, i, N
   N = size(MC%photons) ! total number of photons (in all pulses)
   call photons_per_pulse(laser, numpar, tim, Nphot_pulses)
   sum_phot = SUM(Nphot_pulses(:))      ! total number of photons from all pulses
   if (sum_phot <= 0) then
      Nphot_pulses = 0
   else
      Nphot_pulses = Nphot_pulses/sum_phot ! normalize fractions of different pulses to 1
   endif
   sum_phot = Nphot_pulses(1) ! now, to summation
   do i = 1, N
      call random_number(RN)
      coun = 1 ! start counting pulses
      do while (sum_phot .LT. RN)
         coun = coun + 1 ! next pulse
         sum_phot = sum_phot + Nphot_pulses(coun) ! total number of photons
      enddo

      ! Distribution of the photon energy:
      if (laser(coun)%FWHM_hw > 0.0d0) then  ! if it is a distribution
         MC%photons(i)%E = sample_gaussian(laser(coun)%hw, laser(coun)%FWHM_hw, .true.) ! module "Little_subroutine"
      else  ! single photon energy
         MC%photons(i)%E = laser(coun)%hw ! [eV] photon energy in this pulse #coun
      endif
!       print*, MC%photons(i)%E
   enddo
end subroutine choose_photon_energy


subroutine total_photons(laser, numpar, tim, Nphot)
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar     ! numerical parameters, including lists of earest neighbors
   real(8), intent(in) :: tim      ! [fs] current time
   real(8), intent(inout) :: Nphot ! number of absorbed photons
   real(8), dimension(:), allocatable :: Nphot_pulses
   integer i, N
   N = size(laser) ! how many pulses
   call photons_per_pulse(laser, numpar, tim, Nphot_pulses)
   Nphot = SUM(Nphot_pulses) ! total number in all pulses
end subroutine total_photons


subroutine photons_per_pulse(laser, numpar, tim, Nphot)
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar     ! numerical parameters, including lists of earest neighbors
   real(8), intent(in) :: tim ! [fs] current time
   real(8), dimension(:), allocatable, intent(inout) :: Nphot ! number of absorbed photons
   integer i, N
   N = size(laser) ! how many pulses
   if (.not.allocated(Nphot)) allocate(Nphot(N)) ! absorbed photons for each pulse
   Nphot = 0.0d0
   do i = 1, N ! number of photons in each pulse at this timestep
      SELECT CASE (laser(i)%KOP)   ! Laser pulse shape:
         CASE (0)  ! rectangular pulse
            Nphot(i) = laser(i)%Nph*Rect_pulse(laser(i)%t0, laser(i)%t, tim)*numpar%dt
         CASE (2)     ! SASE pulse
            Nphot(i) = laser(i)%Nph*SASE_pulse(laser(i)%t0, laser(i)%t, tim)*numpar%dt
         CASE default ! gaussian pulse
            Nphot(i) = laser(i)%Nph*Gaus_pulse(laser(i)%t0, laser(i)%t, tim)*numpar%dt
      END SELECT
   enddo
end subroutine photons_per_pulse


subroutine absorbed_photons_dt(laser, numpar, tim, Nph, Nphot)
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar     ! numerical parameters, including lists of earest neighbors
   real(8), intent(in) :: tim	! [fs] current time
   integer, intent(inout) :: Nph	! number of absorbed photons
   real(8), intent(inout), optional :: Nphot
   real(8) Nphot0
   real(8) RN
   integer Nph_cur
   if (.not.present(Nphot)) then
      call total_photons(laser, numpar, tim, Nphot0)
   else
      Nphot0 = Nphot
   endif
   call random_number(RN)
   if ((FLOOR(Nphot0)+RN) .GT. Nphot0) then
      Nph_cur = FLOOR(Nphot0)   ! realized 'Number of photons'
   else
      Nph_cur = CEILING(Nphot0) ! realized 'Number of photons'
   endif
   Nph = Nph + Nph_cur

!    if (Nph > 0) write(*,'(a, i2, i2, e, e)') 'A:', Nph, Nph_cur, Nphot0, RN
end subroutine absorbed_photons_dt


function Rect_pulse(mu, sigma, x) ! number of photons at timestep x, according to flat-top pulse
   real(8) mu, sigma, x, Rect_pulse
   real(8) Gaus_sigma
   Gaus_sigma = sigma !*2.35482d0
   if ((x .GT. (mu - Gaus_sigma/2.0d0)) .AND. (x .LE. (mu + Gaus_sigma/2.0d0)) ) then
      Rect_pulse = 1.0d0/Gaus_sigma
   else
      Rect_pulse = 0.0d0
   endif
end function Rect_pulse

function Gaus_pulse(mu, sigma, x) ! number of photons at the time x according to Gaussian shape
   real(8) mu, sigma, x
   real(8) Gaus_pulse
   !if ((x .GT. (mu - 2.5d0*sigma)) .AND. (x .LE. (mu + 2.5d0*sigma)) ) then ! cutting > 2.5 sigma loses ~1.2% of energy; adjust as needed
   if ((x .GT. (mu - 3.0d0*sigma)) .AND. (x .LE. (mu + 3.0d0*sigma)) ) then
      Gaus_pulse = 1.0d0/(sqrt(2.0d0*g_Pi)*sigma)*dexp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))
   else
      Gaus_pulse = 0.0d0
   endif
end function Gaus_pulse

function SASE_pulse(mu, sigma, x) ! number of photons at the time x according to schematic SASE shape
   real(8) mu, sigma, x
   real(8) SASE_pulse
   integer i
   real(8) RN(4), Co(4), y, f, SASE
   if ((x .GT. (mu - sigma)) .AND. (x .LE. (mu + sigma)) ) then
      RN(1) = 0.4563449303d0
      RN(2) = 0.1271999433d0
      RN(3) = 1.0d0 - RN(1)
      RN(4) = 1.0d0 - RN(2)
      Co(1) = 2.0d0
      Co(2) = 4.0d0
      Co(3) = 5.0d0
      Co(4) = 9.0d0
      SASE = 0.0d0
      do i = 1,4
         y = Co(i)*g_Pi*(x+sigma-mu)/sigma
         f = sin(y)
         SASE = SASE + f*f/sigma*RN(i)
      enddo
      SASE_pulse = SASE/2.0d0 ! normalization to 1
   else
      SASE_pulse = 0.0d0
   endif
end function SASE_pulse
!LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


END MODULE Monte_carlo
