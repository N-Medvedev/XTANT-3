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
module Plots_python

#ifndef __GFORTRAN__
USE IFLPORT, only : system, chdir   ! library, allowing to operate with directories in intel fortran
#endif

use Universal_constants
use Objects
use Dealing_with_files, only : Count_lines_in_file
use Gnuplotting, only : write_gnuplot_script_header_new, write_gnuplot_script_ending_new
use Little_subroutines, only : set_starting_time, order_of_time, name_of_orbitals, number_of_types_of_orbitals
use Dealing_with_EADL, only : define_PQN


implicit none 
PRIVATE


character(50), parameter :: m_Python_plot_all = 'OUTPUT_Python_PLOT_ALL.py'


public :: collect_python_plots, execute_all_pyplots, create_python_plot_scripts, &
            Plot_electron_MFP_python, Plot_photon_MFP_python, Plot_laser_spectrum_python


 contains




subroutine create_python_plot_scripts(Scell,matter,numpar,laser, file_path, file_temperatures, file_pressure, file_energies, &
file_supercell, file_electron_properties, file_heat_capacity, file_heat_capacity_dyn, &
file_numbers, file_orb, file_deep_holes, file_optics, file_Ei, file_PCF, file_NN, file_element_NN, file_electron_entropy, file_Te, file_mu, &
file_atomic_entropy, file_atomic_temperatures, file_atomic_temperatures_part, file_sect_displ, &
file_diffraction_peaks, file_diffraction_peaks_part, file_diffraction_powder, file_diffraction_peaks_DW, file_Debye_temperature, &
file_testmode, file_coupling, file_high_e, file_fragments)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   character(*), intent(in) :: file_path
   character(*), intent(in) :: file_temperatures	! time [fs], Te [K], Ta [K]
   character(*), intent(in) :: file_pressure	! pressure and stress tensore
   character(*), intent(in) :: file_energies	! energies [eV]
   character(*), intent(in) :: file_supercell	! supercell vectors
   character(*), intent(in) :: file_electron_properties	! electron properties
   character(*), intent(in) :: file_heat_capacity  ! electronic heat capacity
   character(*), intent(in) :: file_heat_capacity_dyn  ! electronic heat conductivity dynamical
   character(*), intent(in) :: file_numbers	! total numbers of electrons and holes
   character(*), intent(in) :: file_orb   ! orbital-resolved electron parameters
   character(*), intent(in) :: file_deep_holes	! deep-shell holes
   character(*), intent(in) :: file_optics		! optical coefficients
   character(*), intent(in) :: file_Ei		! energy levels
   character(*), intent(in) :: file_PCF		! pair correlation function
   character(*), intent(in) :: file_NN      ! nearest neighbors
   character(*), dimension(:), allocatable, intent(in) :: file_element_NN      ! element-specific nearest neighbors
   character(*), intent(in) :: file_electron_entropy  ! electron entropy
   character(*), intent(in) :: file_Te ! electron temperatures
   character(*), intent(in) :: file_mu ! electron chem.potentials
   character(*), intent(in) :: file_atomic_entropy ! atomic entropy
   character(*), intent(in) :: file_atomic_temperatures ! atomic temperatures (various definitions)
   character(*), intent(in) :: file_atomic_temperatures_part  ! partial atomic temperatures (X, Y, Z)
   character(*), dimension(:), allocatable, intent(in) :: file_sect_displ
   character(*), intent(in) :: file_diffraction_peaks, file_diffraction_powder, file_diffraction_peaks_DW  ! diffraction peaks
   character(*), dimension(:), allocatable, intent(in) :: file_diffraction_peaks_part
   character(*), intent(in) :: file_Debye_temperature ! Debye temperatures
   character(*), intent(in) :: file_testmode    ! testmode data
   character(*), intent(in) :: file_coupling    ! electron-ion coupling, including partial
   character(*), intent(in) :: file_high_e      ! high-energy electrons by type
   character(*), intent(in) :: file_fragments       ! fragments of the target
   !----------------

   character(300) :: File_name, File_name2
   real(8) :: t0, t_last, x_tics, E_temp
   integer FN, i, j, Nshl, counter, iret, Nsiz
   character(300) :: chtemp, command
   character(11) :: chtemp11, sh_cmd, call_slash
   character(8) :: temp, time_order

   ! Starting time, to give enough time for system to thermalize before the pulse:
   call set_starting_time(laser, t0, numpar%t_start) ! module "Little_subroutines"
   ! Finishing time:
   t_last = numpar%t_total ! total duration of simulation [fs]



   ! Energy levels (molecular orbitals):
   if (numpar%save_Ei) then
      call Python_plot_energy_levels(numpar, t0, t_last, Scell, file_Ei, 'OUTPUT_energy_levels.py')  ! below
   endif


   ! Energies:
   call Python_plot_energies(numpar, file_energies, t0, t_last, 'OUTPUT_energies.py') ! below


   ! Temepratures:
   call Python_plot_temperatures(numpar, matter, file_temperatures, t0, t_last, 'OUTPUT_temperatures.py') ! below


   ! Mean displacement:
   if (abs(numpar%MSD_power) > 1.0e-6) then ! only plot it if it's not zero
      call Python_plot_MSD(matter, numpar, file_temperatures, t0, t_last, &
            'OUTPUT_mean_displacement.py', &
            numpar%MSD_power) ! below
   endif


   ! Atomic masks for sectional displacements:
   if (allocated(Scell(1)%Displ)) then
      Nsiz = size(Scell(1)%Displ)   ! how many masks
      do j = 1, Nsiz    ! for all masks
         call Python_plot_displacements(matter, numpar, file_sect_displ(j), t0, t_last, &
            'OUTPUT_mean_displacements_'//trim(adjustl(Scell(1)%Displ(j)%mask_name))//'.py', &
            Scell(1)%Displ(j)%MSD_power, trim(adjustl(Scell(1)%Displ(j)%mask_name)) ) ! below

         ! Partial by elements, if there is more than one:
         if (matter%N_KAO > 1) then
            call Python_plot_displacements_partial(matter, numpar, file_sect_displ(j), t0, t_last, 'OUTPUT_mean_displacement_'// &
               trim(adjustl(Scell(1)%Displ(j)%mask_name))//'_partial.py', &
               Scell(1)%Displ(j)%MSD_power, trim(adjustl(Scell(1)%Displ(j)%mask_name)) ) ! below
         endif
      enddo ! j
   endif


   ! Diffraction:
   if (numpar%save_diff_peaks) then
      ! Diffraction peaks for selected Miller indices:
      call Python_plot_diffraction_peaks(Scell(1), numpar, file_diffraction_peaks, t0, t_last, &
                                    'OUTPUT_diffraction_peaks.py', &
                                    'OUTPUT_diffraction_peaks', 'Diffraction peaks', &
                                    .false.) ! below

      ! For element-specific diffraction data:
      if (size(matter%Atoms) > 1) then
          do j = 1, size(matter%Atoms)    ! for all elements
             call Python_plot_diffraction_peaks(Scell(1), numpar, file_diffraction_peaks_part(j), t0, t_last, &
                                    'OUTPUT_diffraction_peaks_'//trim(adjustl(matter%Atoms(j)%Name))//'.py', &
                                    'OUTPUT_diffraction_peaks_'//trim(adjustl(matter%Atoms(j)%Name)), &
                                    'Diffraction peaks in '//trim(adjustl(matter%Atoms(j)%Name)), &
                                    .false.) ! below
          enddo
      endif

      ! Powder diffraction:
      call Python_plot_powder_diffraction(Scell(1), matter, numpar, file_diffraction_powder, 'OUTPUT_diffraction_powder.py', &
                                          trim(adjustl(numpar%vid_extention)))      ! below

      ! Check if Debye-Waller analysis is required:
      if ( abs(numpar%DW_theta) > 1.0d-6 ) then
         call Python_plot_diffraction_peaks(Scell(1), numpar, file_diffraction_peaks_DW, t0, t_last, &
                                    'OUTPUT_diffraction_peaks_DW.py', &
                                    'OUTPUT_diffraction_peaks_DW', 'Debye Waller peaks',  &
                                    .true.) ! below

         call Python_plot_Debye_temperatures(Scell(1), numpar, file_Debye_temperature, t0, t_last, &
                                    'OUTPUT_Debye_temperatures.py') ! below
      endif
   endif


   ! Pressure:
   call Python_plot_pressure(numpar, file_pressure, t0, t_last, 'OUTPUT_pressure.py') ! below


   ! Stress tensor:
   call Python_plot_stress(numpar, file_pressure, t0, t_last, 'OUTPUT_stress_tensor.py') ! below


   ! Numbers of particles:
   call Python_plot_numbers(numpar, file_numbers, t0, t_last, 'OUTPUT_electrons_holes.py') ! below

   ! High-energy electrons_by_type:
   call Python_plot_high_energy_el(numpar, file_high_e, t0, t_last, 'OUTPUT_electrons_high_energy.py') ! below


   ! Orbital-resolved electron parameters:
   call Python_plot_orbital_resolved(Scell(1), matter, numpar, file_orb, t0, t_last, 'OUTPUT_orbital_resolved_Ne.py') ! below


   ! Numbers of CB electrons:
   call Python_plot_CB_electrons(numpar, file_numbers, t0, t_last, 'OUTPUT_CB_electrons.py')     ! below


   ! Numbers of deep-shell holes:
   call Python_plot_holes(matter, numpar, file_deep_holes, t0, t_last, 'OUTPUT_deep_shell_holes.py') ! below


   ! Band gap:
   call Python_plot_Egap(numpar, file_electron_properties, t0, t_last, 'OUTPUT_Egap.py') ! below


   ! Chemical potential and Ne:
   call Python_plot_mu(numpar, file_electron_properties, t0, t_last, 'OUTPUT_mu_and_Ne.py') ! below


   ! Boundaries of the bands:
   call Python_plot_Ebands(numpar, file_electron_properties, t0, t_last, 'OUTPUT_bands.py') ! below


   ! Electron heat capacity:
   call Python_plot_capacity(numpar, file_electron_properties, t0, t_last, 'OUTPUT_electron_Ce.py') ! below


   ! Electron heat conductivity:
   if (numpar%do_kappa) then
      call Python_plot_heat_conductivity(numpar, file_heat_capacity, numpar%kappa_Te_min, numpar%kappa_Te_max, &
            'OUTPUT_electron_heat_conductivity.py') ! below
   endif
   if (numpar%do_kappa_dyn) then
      call Python_plot_heat_conductivity_dyn(numpar, file_heat_capacity_dyn, numpar%kappa_Te_min, numpar%kappa_Te_max, &
            'OUTPUT_electron_heat_conductivity_dyn.py') ! below
   endif


   ! Electron entropy:
   call Python_plot_entropy(numpar, file_electron_entropy, t0, t_last, 'OUTPUT_electron_entropy.py') ! below


   ! Electron temperatures and chemical potential (for band-resolved calculations):
   if (numpar%do_partial_thermal) then
      call Python_plot_el_temperatures(numpar, file_Te, t0, t_last, 'OUTPUT_electron_temperatures.py') ! below
      call Python_plot_chempots(numpar, file_mu, t0, t_last, 'OUTPUT_electron_chempotentials.py') ! below
   endif


   ! Atomic temperatures (various definitions):
   if (numpar%print_Ta) then
      ! Atomic entropy:
      call Python_plot_entropy_atomic(numpar, file_atomic_entropy, t0, t_last, 'OUTPUT_atomic_entropy.py') ! below

      call Python_plot_at_temperatures(numpar, file_atomic_temperatures, t0, t_last, 'OUTPUT_atomic_temperatures.py') ! below

      call Python_plot_at_temperatures_part(numpar, file_atomic_temperatures_part, t0, t_last, 'OUTPUT_atomic_temperatures_partial.py') ! below
   endif


   ! Data for fragments (if required):
   if (numpar%print_fragments) then
      call Python_plot_fragments_data(Scell(1), numpar, file_fragments, t0, t_last, &
            'OUTPUT_fragments_Nat.py', 'OUTPUT_fragments_Ta.py', 'OUTPUT_fragments_Ne.py') ! below
   endif


   ! Electron-ion coupling parameter:
   call Python_plot_coupling(numpar, file_electron_properties, t0, t_last, 'OUTPUT_coupling.py') ! below

   ! Partial electron-ion coupling:
   call Python_plot_partial_coupling(Scell(1), matter, numpar, file_coupling, t0, t_last, &
            'OUTPUT_coupling_by_element.py', 'OUTPUT_coupling_by_orbital.py') ! below


   ! Volume:
   call Python_plot_volume(numpar, file_supercell, t0, t_last, 'OUTPUT_volume.py') ! below


   ! Mulliken charges:
   if (numpar%Mulliken_model >= 1) then
      call Python_plot_Mulliken_charges(matter, numpar, file_electron_properties, t0, t_last, 'OUTPUT_Mulliken_charges.py') ! below
   endif


   ! Nearest neighbors:
   if (numpar%save_NN) then
      call Python_plot_nearest_neighbors(numpar, file_NN, t0, t_last, 'OUTPUT_nearest_neighbors.py') ! below
   endif


   ! Element-specific nearest neighbors:
   if (allocated(numpar%NN_radii)) then
      do i = 1, size(numpar%NN_radii) ! for all requested elements
         call Python_plot_nearest_neighbors_elements(matter, numpar, file_element_NN(i), trim(adjustl(numpar%NN_radii(i)%Name)), &
              t0, t_last, &
              'OUTPUT_nearest_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'.py') ! below
      enddo ! i
   endif


   ! Pair correlation function:
   if (numpar%save_PCF) then
      ! Pair correlation function can only be plotted as animated plot:
      call Python_plot_pair_correlation(Scell(1), matter, numpar, 'OUTPUT_pair_correlation_function.dat', &
                                          'OUTPUT_pair_correlation.py', trim(adjustl(numpar%vid_extention)))   ! below
   endif


   ! Distribution function of electrons:
   if (numpar%save_fe) then
      ! Distribution function can only be plotted as animated plot:
      call Python_plot_distribution(numpar, 'OUTPUT_electron_distribution.dat', &
                                    'OUTPUT_electron_distribution.py', trim(adjustl(numpar%vid_extention)))   ! below
   endif


   ! Orbital-resoloved distribution function of electrons:
   if (numpar%save_fe_orb) then
      ! Distribution function can only be plotted as animated plot:
      call Python_plot_orb_distribution(Scell(1), matter, numpar, 'OUTPUT_electron_distribution.dat', &
                                    'OUTPUT_orbital_resolved_fe.py', trim(adjustl(numpar%vid_extention)))   ! below
   endif


   ! Distribution function of all electrons on the grid:
   if (numpar%save_fe_grid) then
      ! Distribution function can only be plotted as animated plot:
      call Python_plot_distribution_on_grid(Scell(1), numpar, 'OUTPUT_electron_distribution_on_grid.dat', &
                                    'OUTPUT_electron_distribution_on_grid.py', trim(adjustl(numpar%vid_extention)))   ! below
   endif



   ! Distribution function of atoms:
   if (numpar%save_fa) then
      ! 1) Distribution of kinetic energies:
      call Python_plot_atomic_distribution(Scell(1), numpar, 'OUTPUT_atomic_distribution.dat', 'OUTPUT_atomic_distribution_kin', &
           'Distribution of kinetic energies', 'OUTPUT_atoms_distribution.py', trim(adjustl(numpar%vid_extention)))   ! below

      ! 2) Distribution of potential energies:
      call Python_plot_atomic_distribution(Scell(1), numpar, 'OUTPUT_atomic_distribution_pot.dat', 'OUTPUT_atomic_distribution_pot', &
           'Distribution of potential energies', 'OUTPUT_atomic_distribution_pot.py', trim(adjustl(numpar%vid_extention)))   ! below

      ! 3) Distribution of total energies:
      call Python_plot_atomic_distribution(Scell(1), numpar, 'OUTPUT_atomic_distribution_tot.dat', 'OUTPUT_atomic_distribution_tot', &
       'Distribution of total energies', 'OUTPUT_atomic_distribution_tot.py', trim(adjustl(numpar%vid_extention)), skip_eq=.true.)   ! below
   endif



   ! DOS of electrons:
   if (numpar%save_DOS) then  ! Material DOS
      select case (numpar%DOS_splitting)
      case (1) ! with partial DOS
         ! DOS can only be plotted as animated plot:
         call Python_plot_DOS(Scell(1), matter, numpar, 'OUTPUT_DOS.dat', 'OUTPUT_DOS.py', trim(adjustl(numpar%vid_extention)))   ! below
      endselect
   endif



   ! Optical coefficients
   if (numpar%do_drude) then
      call Python_plot_optical_coefficients(numpar, file_optics, t0, t_last, 'OUTPUT_optical_coefficients.py') ! below
      ! also n and k:
      call Python_plot_n_and_k(numpar, file_optics, t0, t_last, 'OUTPUT_optical_n_and_k.py') ! below
   endif


   !cccccccccccccccccccccccccccccccccccccccccccccc
   ! Create also convolved plots:
   CONV:if (Scell(1)%eps%tau > 0.0d0) then ! convolved files too:
      ! Energies:
      call Python_plot_energies(numpar, file_energies, t0, t_last, 'OUTPUT_energies_CONVOLVED.py', convolved=.true.) ! below


      ! Temepratures:
      call Python_plot_temperatures(numpar, matter, file_temperatures, t0, t_last, 'OUTPUT_temepratures_CONVOLVED.py', convolved=.true.) ! below


      ! Mean displacement:
      if (abs(numpar%MSD_power) > 1.0e-6) then ! only plot it if it's not zero
            call Python_plot_MSD(matter, numpar, file_temperatures, t0, t_last, &
                  'OUTPUT_mean_displacement_CONVOLVED.py', &
                  numpar%MSD_power, convolved=.true.) ! below
      endif


      ! Atomic masks for sectional displacements:
      if (allocated(Scell(1)%Displ)) then
            Nsiz = size(Scell(1)%Displ)   ! how many masks
            do j = 1, Nsiz    ! for all masks
            call Python_plot_displacements(matter, numpar, file_sect_displ(j), t0, t_last, &
                  'OUTPUT_mean_displacements_'//trim(adjustl(Scell(1)%Displ(j)%mask_name))//'CONVOLVED.py', &
                  Scell(1)%Displ(j)%MSD_power, trim(adjustl(Scell(1)%Displ(j)%mask_name)) , convolved=.true.) ! below

            ! Partial by elements, if there is more than one:
            if (matter%N_KAO > 1) then
                  call Python_plot_displacements_partial(matter, numpar, file_sect_displ(j), t0, t_last, 'OUTPUT_mean_displacement_'// &
                  trim(adjustl(Scell(1)%Displ(j)%mask_name))//'_partial_CONVOLVED.py', &
                  Scell(1)%Displ(j)%MSD_power, trim(adjustl(Scell(1)%Displ(j)%mask_name)) , convolved=.true.) ! below
            endif
            enddo ! j
      endif


      ! Diffraction:
      if (numpar%save_diff_peaks) then
            ! Diffraction peaks for selected Miller indices:
            call Python_plot_diffraction_peaks(Scell(1), numpar, file_diffraction_peaks, t0, t_last, &
                                          'OUTPUT_diffraction_peaks_CONVOLVED.py', &
                                          'OUTPUT_diffraction_peaks', 'Diffraction peaks', &
                                          .false., convolved=.true.) ! below

            ! For element-specific diffraction data:
            if (size(matter%Atoms) > 1) then
               do j = 1, size(matter%Atoms)    ! for all elements
                  call Python_plot_diffraction_peaks(Scell(1), numpar, file_diffraction_peaks_part(j), t0, t_last, &
                                          'OUTPUT_diffraction_peaks_'//trim(adjustl(matter%Atoms(j)%Name))//'_CONVOLVED.py', &
                                          'OUTPUT_diffraction_peaks_'//trim(adjustl(matter%Atoms(j)%Name)), &
                                          'Diffraction peaks in '//trim(adjustl(matter%Atoms(j)%Name)), &
                                          .false., convolved=.true.) ! below
               enddo
            endif

            ! Check if Debye-Waller analysis is required:
            if ( abs(numpar%DW_theta) > 1.0d-6 ) then
               call Python_plot_diffraction_peaks(Scell(1), numpar, file_diffraction_peaks_DW, t0, t_last, &
                                          'OUTPUT_diffraction_peaks_DW_CONVOLVED.py', &
                                          'OUTPUT_diffraction_peaks_DW', 'Debye Waller peaks',  &
                                          .true., convolved=.true.) ! below

               call Python_plot_Debye_temperatures(Scell(1), numpar, file_Debye_temperature, t0, t_last, &
                                          'OUTPUT_Debye_temperatures_CONVOLVED.py', convolved=.true.) ! below
            endif
      endif


      ! Pressure:
      call Python_plot_pressure(numpar, file_pressure, t0, t_last, 'OUTPUT_pressure_CONVOLVED.py', convolved=.true.) ! below


      ! Stress tensor:
      call Python_plot_stress(numpar, file_pressure, t0, t_last, 'OUTPUT_stress_tensor_CONVOLVED.py', convolved=.true.) ! below


      ! Numbers of particles:
      call Python_plot_numbers(numpar, file_numbers, t0, t_last, 'OUTPUT_electrons_holes_CONVOLVED.py', convolved=.true.) ! below


      ! High-energy electrons_by_type:
      call Python_plot_high_energy_el(numpar, file_high_e, t0, t_last, 'OUTPUT_electrons_high_energy_CONVOLVED.py', convolved=.true.) ! below


      ! Orbital-resolved electron parameters:
      call Python_plot_orbital_resolved(Scell(1), matter, numpar, file_orb, t0, t_last, 'OUTPUT_orbital_resolved_Ne_CONVOLVED.py', convolved=.true.) ! below


      ! Numbers of CB electrons:
      call Python_plot_CB_electrons(numpar, file_numbers, t0, t_last, 'OUTPUT_CB_electrons_CONVOLVED.py', convolved=.true.)     ! below


      ! Numbers of deep-shell holes:
      call Python_plot_holes(matter, numpar, file_deep_holes, t0, t_last, 'OUTPUT_deep_shell_holes_CONVOLVED.py', convolved=.true.) ! below


      ! Band gap:
      call Python_plot_Egap(numpar, file_electron_properties, t0, t_last, 'OUTPUT_Egap_CONVOLVED.py', convolved=.true.) ! below


      ! Chemical potential and Ne:
      call Python_plot_mu(numpar, file_electron_properties, t0, t_last, 'OUTPUT_mu_and_Ne_CONVOLVED.py', convolved=.true.) ! below


      ! Boundaries of the bands:
      call Python_plot_Ebands(numpar, file_electron_properties, t0, t_last, 'OUTPUT_bands_CONVOLVED.py', convolved=.true.) ! below


      ! Electron heat capacity:
      call Python_plot_capacity(numpar, file_electron_properties, t0, t_last, 'OUTPUT_electron_Ce_CONVOLVED.py', convolved=.true.) ! below


      ! Electron entropy:
      call Python_plot_entropy(numpar, file_electron_entropy, t0, t_last, 'OUTPUT_electron_entropy.py', convolved=.true.) ! below


      ! Electron temperatures and chemical potential (for band-resolved calculations):
      if (numpar%do_partial_thermal) then
            call Python_plot_el_temperatures(numpar, file_Te, t0, t_last, 'OUTPUT_electron_temperatures_CONVOLVED.py', convolved=.true.) ! below
            call Python_plot_chempots(numpar, file_mu, t0, t_last, 'OUTPUT_electron_chempotentials_CONVOLVED.py', convolved=.true.) ! below
      endif


      ! Atomic temperatures (various definitions):
      if (numpar%print_Ta) then
            ! Atomic entropy:
            call Python_plot_entropy_atomic(numpar, file_atomic_entropy, t0, t_last, 'OUTPUT_atomic_entropy_CONVOLVED.py', convolved=.true.) ! below

            call Python_plot_at_temperatures(numpar, file_atomic_temperatures, t0, t_last, 'OUTPUT_atomic_temperatures_CONVOLVED.py', convolved=.true.) ! below

            call Python_plot_at_temperatures_part(numpar, file_atomic_temperatures_part, t0, t_last, 'OUTPUT_atomic_temperatures_partial_CONVOLVED.py', convolved=.true.) ! below
      endif


      ! Electron-ion coupling parameter:
      call Python_plot_coupling(numpar, file_electron_properties, t0, t_last, 'OUTPUT_coupling_CONVOLVED.py', convolved=.true.) ! below

      ! Partial electron-ion coupling:
      call Python_plot_partial_coupling(Scell(1), matter, numpar, file_coupling, t0, t_last, 'OUTPUT_coupling_by_element.py', &
            'OUTPUT_coupling_by_orbital.py', convolved=.true.) ! below


      ! Volume:
      call Python_plot_volume(numpar, file_supercell, t0, t_last, 'OUTPUT_volume_CONVOLVED.py', convolved=.true.) ! below


      ! Mulliken charges:
      if (numpar%Mulliken_model >= 1) then
            call Python_plot_Mulliken_charges(matter, numpar, file_electron_properties, t0, t_last, 'OUTPUT_Mulliken_charges_CONVOLVED.py', convolved=.true.) ! below
      endif


      ! Nearest neighbors:
      if (numpar%save_NN) then
            call Python_plot_nearest_neighbors(numpar, file_NN, t0, t_last, 'OUTPUT_nearest_neighbors_CONVOLVED.py', convolved=.true.) ! below
      endif


      ! Element-specific nearest neighbors:
      if (allocated(numpar%NN_radii)) then
            do i = 1, size(numpar%NN_radii) ! for all requested elements
            call Python_plot_nearest_neighbors_elements(matter, numpar, file_element_NN(i), trim(adjustl(numpar%NN_radii(i)%Name)), &
                  t0, t_last, &
                  'OUTPUT_nearest_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'_CONVOLVED.py', convolved=.true.) ! below
            enddo ! i
      endif


      ! Optical coefficients
      if (numpar%do_drude) then
            call Python_plot_optical_coefficients(numpar, file_optics, t0, t_last, 'OUTPUT_optical_coefficients_CONVOLVED.py', convolved=.true.) ! below
            ! also n and k:
            call Python_plot_n_and_k(numpar, file_optics, t0, t_last, 'OUTPUT_optical_n_and_k_CONVOLVED.py', convolved=.true.) ! below
      endif

   endif CONV
end subroutine create_python_plot_scripts





subroutine Python_plot_n_and_k(numpar, file_optics, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_optics, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 2

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 4
   col_labels(1) = '"n"'
   call select_linestyle(1,linestyle(1))
   col_nums(2) = 5
   col_labels(2) = '"k"'
   call select_linestyle(2,linestyle(2))


   Plot_name = 'OUTPUT_optical_n_and_k'
   Data_file_name = file_optics
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_optics(1:len(trim(adjustl(file_optics)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Optical n and k', 'Complex refractive index', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_n_and_k



subroutine Python_plot_optical_coefficients(numpar, file_optics, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_optics, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 3

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 1
   col_labels(1) = '"Reflectivity"'
   call select_linestyle(1,linestyle(1))
   col_nums(2) = 2
   col_labels(2) = '"Transmission"'
   call select_linestyle(2,linestyle(2))
   col_nums(3) = 3
   col_labels(3) = '"Absorption"'
   call select_linestyle(3,linestyle(3))


   Plot_name = 'OUTPUT_optical_coefficients'
   Data_file_name = file_optics
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_optics(1:len(trim(adjustl(file_optics)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Optical coefficients', 'Optical coefficients', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, y_max = 1.0d0)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_optical_coefficients





subroutine Python_plot_DOS(Scell, matter, numpar, file_DOS, script_name, video_format)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_DOS, script_name ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   !----------------
   integer :: FN, i, i_start, ind, j, k, Nsiz, i_at, N_types, NKOA, N_at, norb, i_col, i_cur, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle, symbols
   character(300) :: File_name
   character(30) :: chtemp1
   real(8) :: x_start, x_end

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Plot limits:
   x_start = -25.0d0
   x_end = 25.0d0

   ! Find number of orbitals per atom:
   NKOA = matter%N_KAO      ! number of kinds of atoms
   N_at = size(Scell%MDatoms) ! number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   norb = Nsiz/N_at ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

   ! allocate the arrays:
   Nsiz = NKOA * N_types + 1     ! assuming basis set is the same for all elements
   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 1
   col_labels(1) = 'Total'
   call select_linestyle(1, linestyle(1)) ! below

   i_col = 1      ! column number after which orbital-resolved data start
   i_cur = 1      ! to start with
   do i_at = 1, NKOA
      do i_types = 1, N_types
         i_cur = i_cur + 1    ! count columns
         ! Set the arrays:
         col_nums(i_cur) = (i_col - 1) + i_cur    ! index for python
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
         col_labels(i_cur) = trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
         call select_linestyle(i_cur, linestyle(i_cur)) ! below
      enddo ! i_types
   enddo ! i_at

   call Create_Python_animation(FN, file_DOS, col_nums, col_labels, &
      'Energy (eV)', 'DOS (electrons/eV)', 'Electronic density of states', &
      "best", 'OUTPUT_DOS', trim(adjustl(video_format)), &
      x_min=x_start, x_max=x_end, y_min=0.0d0, &
      t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total, &
      l_style=linestyle, &
      first_line_in_front =.false.)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_DOS




subroutine Python_plot_atomic_distribution(Scell, numpar, file_distribution, plot_name, plot_title, script_name, video_format, skip_eq)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_distribution, plot_name, script_name, plot_title ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   logical, intent(in), optional :: skip_eq     ! skip plotting equilibrium (equivalent) distribution
   !----------------
   integer :: FN, i, i_start, ind, j, k, Nsiz, i_at, N_types, NKOA, N_at, norb, i_col, i_cur, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle, symbols
   character(300) :: File_name
   character(30) :: chtemp1
   real(8) :: x_start, x_end
   logical :: skip_eq_plot


   skip_eq_plot = .false. ! default
   if (present(skip_eq)) then
      if (skip_eq) skip_eq_plot = skip_eq
   endif

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! allocate the arrays:
   if (skip_eq_plot) then
      Nsiz = 1   ! one less
   else
      Nsiz = 2 ! default
   endif

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))
   allocate(symbols(Nsiz))

   if (.not.skip_eq_plot) then
      col_nums(1) = 2
      col_labels(1) = 'Equilibrium'
      call select_linestyle(1, linestyle(1)) ! below
      symbols(1) = '"None"'
      col_nums(2) = 1
      col_labels(2) = 'Nonequilibrium'
      linestyle(2) = '"None"'
      call select_symbols(1, symbols(2))     ! below
   else
      col_nums(1) = 1
      col_labels(1) = 'Nonequilibrium'
      linestyle(1) = '"None"'
      call select_symbols(1, symbols(1))     ! below
   endif

   call Create_Python_animation(FN, file_distribution, col_nums, col_labels, &
      'Energy (eV)', 'Atomic distribution (a.u.)', trim(adjustl(plot_title)), &
      "best", trim(adjustl(plot_name)), trim(adjustl(video_format)), &
      t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total, &
      l_style=linestyle, symbols=symbols, &
      first_line_in_front =.false.)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle, symbols)
end subroutine Python_plot_atomic_distribution





subroutine Python_plot_distribution_on_grid(Scell, numpar, file_distribution, script_name, video_format)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_distribution, script_name ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   !----------------
   integer :: FN, i, i_start, ind, j, k, Nsiz, i_at, N_types, NKOA, N_at, norb, i_col, i_cur, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle, symbols
   character(300) :: File_name
   character(30) :: chtemp1
   real(8) :: x_start, x_end

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Plot limits:
   x_start = -25.0d0
   x_end = Scell%E_fe_grid(size(Scell%E_fe_grid))

   ! allocate the arrays:
   Nsiz = 4
   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))
   allocate(symbols(Nsiz))

   col_nums(1) =1
   col_labels(1) = 'Total'
   linestyle(1) = '"None"'
   call select_symbols(1, symbols(1))     ! below

   col_nums(2) =3
   col_labels(2) = 'Photoelectrons'
   linestyle(2) = '"None"'
   call select_symbols(2, symbols(2))     ! below

   col_nums(3) =4
   col_labels(3) = 'Impact electrons'
   linestyle(3) = '"None"'
   call select_symbols(3, symbols(3))     ! below

   col_nums(4) =5
   col_labels(4) = 'Auger electrons'
   linestyle(4) = '"None"'
   call select_symbols(4, symbols(4))     ! below

   call Create_Python_animation(FN, file_distribution, col_nums, col_labels, &
      'Energy (eV)', 'Electron distribution (1/(V*E))', 'Electron spectrum', &
      "best", 'OUTPUT_electron_distribution_on_grid', trim(adjustl(video_format)), &
      x_min=x_start, x_max=x_end, y_min=1.0d-6, set_y_log=.true., &
      t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total, &
      l_style=linestyle, symbols=symbols, &
      first_line_in_front =.false.)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle, symbols)
end subroutine Python_plot_distribution_on_grid





subroutine Python_plot_orb_distribution(Scell, matter, numpar, file_distribution, script_name, video_format)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_distribution, script_name ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   !----------------
   integer :: FN, i, i_start, ind, j, k, Nsiz, i_at, N_types, NKOA, N_at, norb, i_col, i_cur, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle, symbols
   character(300) :: File_name
   character(30) :: chtemp1
   real(8) :: x_start, x_end

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Plot limits:
   x_start = -25.0d0
   x_end = 25.0d0

   ! Find number of orbitals per atom:
   NKOA = matter%N_KAO      ! number of kinds of atoms
   N_at = size(Scell%MDatoms) ! number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   norb = Nsiz/N_at ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

   ! allocate the arrays:
   Nsiz = NKOA * N_types      ! assuming basis set is the same for all elements
   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))
   allocate(symbols(Nsiz))

   if (numpar%do_partial_thermal) then ! includes band-resolved equivalent distributions
      i_col = 5  ! column number after which orbital-resolved data start
   else
      select case (numpar%el_ion_scheme)
      case (3:5)
         i_col = 3  ! column number after which orbital-resolved data start
      case default
         i_col = 2  ! column number after which orbital-resolved data start
      endselect
   endif

   i_cur = 0      ! to start with
   do i_at = 1, NKOA
      do i_types = 1, N_types
         i_cur = i_cur + 1    ! count columns
         ! Set the arrays:
         col_nums(i_cur) = (i_col - 1) + i_cur    ! index for python
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
         col_labels(i_cur) = trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
         linestyle(i_cur) = '"None"'
         call select_symbols(i_cur, symbols(i_cur))     ! below
      enddo ! i_types
   enddo ! i_at

   call Create_Python_animation(FN, file_distribution, col_nums, col_labels, &
      'Energy (eV)', 'Electron distribution (1/eV)', 'Orbital-resolved electron distribution', &
      "best", 'OUTPUT_orbital_resolved_fe', trim(adjustl(video_format)), &
      x_min=x_start, x_max=x_end, y_min=0.0d0, &
      t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total, &
      l_style=linestyle, symbols=symbols, &
      first_line_in_front =.false.)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle, symbols)
end subroutine Python_plot_orb_distribution



subroutine Python_plot_distribution(numpar, file_distribution, script_name, video_format)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_distribution, script_name ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   !----------------
   integer :: FN, i, i_start, ind, j, k, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle, symbols
   character(300) :: File_name
   real(8) :: x_start, x_end

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Plot limits:
   x_start = -25.0d0
   x_end = 25.0d0

   if (numpar%do_partial_thermal) then ! with band-resolved equivalent distributions
      Nsiz = 4
   else ! only full distribution
      Nsiz = 2
   endif

   ! allocate the arrays:
   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))
   allocate(symbols(Nsiz))

   ! Set the arrays:
   col_nums(1) = 2
   col_labels(1) = 'Equivalent Fermi'
   linestyle(1) = '"-"'
   symbols(1) = '"None"'

   if (numpar%do_partial_thermal) then ! with band-resolved equivalent distributions
      col_nums(2) = 3
      col_labels(2) = 'f$_{eq}$ in valence band'
      call select_linestyle(2, linestyle(2))      ! below
      symbols(2) = '"None"'

      col_nums(3) = 4
      col_labels(3) = 'f$_{eq}$ in conduction band'
      call select_linestyle(3, linestyle(3))      ! below
      symbols(3) = '"None"'

      ind = 4     ! for the last line
   else
      ind = 2     ! for the last line
   endif

   col_nums(ind) = 1
   col_labels(ind) = 'Nonequilibrium'
   linestyle(ind) = '"None"'
   symbols(ind) = '"o"'


   call Create_Python_animation(FN, file_distribution, col_nums, col_labels, &
      'Energy (eV)', 'Electron distribution (1/eV)', 'Electron distribution', &
      "best", 'OUTPUT_electron_distribution', trim(adjustl(video_format)), &
      x_min=x_start, x_max=x_end, y_min=0.0d0, &
      t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total, &
      l_style=linestyle, symbols=symbols, &
      first_line_in_front =.false.)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle, symbols)
end subroutine Python_plot_distribution




subroutine Python_plot_pair_correlation(Scell, matter, numpar, file_pair_correlation, script_name, video_format)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_pair_correlation, script_name ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   !----------------
   integer :: FN, i, i_start, ind, j, k, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name
   real(8) :: x_start, x_end

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Plot limits:
   x_start = 1.0d0
   x_end = 1+ceiling( 0.5d0* abs (min (Scell%supce(1,1), Scell%supce(2,2), Scell%supce(3,3) ) ) ) ! half of the supercell

   ! Get the number of columns to print:
   ind = 1     ! top start counting columns
   if (matter%N_KAO > 1) then ! no need to do for 1 element
      do j = 1, matter%N_KAO     ! for all elements
         do k = j, matter%N_KAO  ! for all different pairs
            ind = ind + 1
         enddo
      enddo
   endif
   Nsiz = ind
   ! allocate the arrays:
   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))
   ! Set the arrays:
   col_nums(1) = 1
   col_labels(1) = 'Total'
   linestyle(1) = '"-"'
   ind = 1  ! restart
   if (matter%N_KAO > 1) then ! no need to do for 1 element
      do j = 1, matter%N_KAO     ! for all elements
         do k = j, matter%N_KAO  ! for all different pairs
            ind = ind + 1
            col_nums(ind) = ind
            col_labels(ind) = trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))
            call select_linestyle(ind, linestyle(ind))      ! below
         enddo
      enddo
   endif


   call Create_Python_animation(FN, file_pair_correlation, col_nums, col_labels, &
      'Radius (A)', 'Pair correlation function (a.u.)', 'Pair correlation function', &
      "best", 'OUTPUT_pair_correlation', trim(adjustl(video_format)), &
      x_min=x_start, x_max=x_end, y_min=0.0d0, &
      t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total, l_style=linestyle, &
      first_line_in_front =.false.)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_pair_correlation





subroutine Python_plot_nearest_neighbors_elements(matter, numpar, file_element_NN, element_name, t0, t_last, script_name, convolved)
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: element_name, file_element_NN, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz, i_cur, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = matter%N_KAO+1

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 1
   col_labels(1) = '"Total"'
   linestyle(1) = '"-"'

   do i = 1, matter%N_KAO
      col_nums(i+1) = 1+i
      col_labels(i+1) = '"'//trim(adjustl(matter%Atoms(i)%Name))//'"'
      call select_linestyle(i+1, linestyle(i+1)) ! below
   enddo


   Plot_name = 'OUTPUT_nearest_neighbors_'//trim(adjustl(element_name))
   Data_file_name = file_element_NN
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_element_NN(1:len(trim(adjustl(file_element_NN)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Nearest neighbors fraction', 'Nearest neighbors of '//trim(adjustl(element_name)), &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_nearest_neighbors_elements





subroutine Python_plot_nearest_neighbors(numpar, file_NN, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_NN, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz, i_cur, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 7

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 2
   col_labels(1) = '"Single atom"'
   linestyle(1) = '"-"'

   col_nums(2) = 3
   col_labels(2) = '"One neighbor"'
   linestyle(2) = '"--"'

   col_nums(3) = 4
   col_labels(3) = '"Two neighbors"'
   linestyle(3) = '"-."'

   col_nums(4) = 5
   col_labels(4) = '"Three neighbors"'
   linestyle(4) = '(0,(5,2,1,2,1,2))'

   col_nums(5) = 6
   col_labels(5) = '"Four neighbors"'
   linestyle(5) = '(0,(10,3,4,3,1,2))'

   col_nums(6) = 7
   col_labels(6) = '"Five neighbors"'
   linestyle(6) = '":"'

   col_nums(7) = 8
   col_labels(7) = '"Six neighbors"'
   linestyle(7) = '(0,(10,3,1,2,4,3,1,2))'


   Plot_name = 'OUTPUT_nearest_neighbors'
   Data_file_name = file_NN
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_NN(1:len(trim(adjustl(file_NN)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Nearest neighbors fraction', 'Number of nearest neighbors', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_nearest_neighbors



subroutine Python_plot_Mulliken_charges(matter, numpar, file_electron_properties, t0, t_last, script_name, convolved)
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_properties, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz, i_cur, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = size(matter%Atoms)

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   i_start = 9
   do i = 1, matter%N_KAO   ! intermediate elements
      i_cur = i_start + i
      col_nums(i) = i_cur
      write(col_labels(i),'(a)') '"'//matter%Atoms(i)%Name//'"'
      !set different style for different elements:
      call select_linestyle(i, linestyle(i))      ! below
   enddo



   Plot_name = 'OUTPUT_Mulliken_charges'
   Data_file_name = file_electron_properties
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat'
      endif
   endif


   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Mulliken charge (e)', 'Mulliken charge', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_Mulliken_charges





subroutine Python_plot_volume(numpar, file_supercell, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_supercell, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 1

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))

   col_nums(1) = 1
   col_labels(1) = '"Volume"'


   Plot_name = 'OUTPUT_volume'
   Data_file_name = file_supercell
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_supercell(1:len(trim(adjustl(file_supercell)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Volume (A$^3$)', 'Volume', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_volume



subroutine Python_plot_partial_coupling(Scell, matter, numpar, file_coupling, t0, t_last, script_name_element, script_name_shell, convolved)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter      ! material properties
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_coupling, script_name_element, script_name_shell ! file with energy levels, scripts
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz, Nat, j, col1, col2, N_types, norb, i_start, i_orb, j_orb, N_at, N_cols
   character(50) :: x_min, x_max, temp_col1, temp_col2, linestyle, orb_name1, orb_name2
   character(300) :: File_name, Plot_name, Data_file_name, temp_txt

   ! Only do it if the data-file exists:
   select case (numpar%DOS_splitting)  ! orbital-resolved data
   case (1) ! with partial DOS
      Nat = size(matter%Atoms)      ! number of different kinds of atoms

      ! 1) Create element-resolved script, only if there are more than one element:
      if (Nat > 1) then
         ! Py script file:
         File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name_element))
         open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

         Plot_name = 'OUTPUT_coupling_parameter_by_element'
         Data_file_name = file_coupling
         if (present(convolved)) then
            if (convolved) then
               Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
               Data_file_name  = trim(adjustl(file_coupling(1:len(trim(adjustl(file_coupling)))-4)))//'_CONVOLVED.dat'
            endif
         endif

         ! Create the python script:
         write(FN, '(a)')  'import pandas as pd'
         write(FN, '(a)')  'import numpy as np'
         write(FN, '(a)')  'import matplotlib.pyplot as plt'
         write(FN, '(a)')  'import matplotlib.cm as cm'
         write(FN, '(a)')  '# Set the format:'
         write(FN, '(a)')  'plt.xlabel("Time (fs)", fontsize=14)'
         write(FN, '(a)')  'plt.ylabel("Coupling parameter (W/(m$^3$ K))", fontsize=14)'
         write(FN, '(a)')  'plt.xticks(fontsize=12)'
         write(FN, '(a)')  'plt.yticks(fontsize=12)'
         write(FN, '(a)')  'plt.grid(False)'
         write(FN, '(a)')  '# Set the axes:'
         write(FN, '(a)')  '# Read the output file:'
         write(FN, '(a)')  'df = pd.read_csv(r"'// trim(adjustl(Data_file_name))// &
                        '", sep=r'//"'\s+',"//'header=None, comment="#", skipinitialspace=True)'
         write(FN, '(a)')  '# Prepare the plot:'
         write(FN, '(a)')  '# -----------------------------'
         write(FN, '(a)')  '# Time axis'
         write(FN, '(a)')  '# -----------------------------'
         write(FN, '(a)')  't = df.iloc[:, 0]'

         ! To add up off-diagonal contributions:
         write(FN, '(a)')  'pairs = ['
         ! First ,add the total coupling:
         write(FN,'(a)') '(1,1, "Total", "-"),'
         ! Now, all element-resolved ones:
         do i = 1, Nat ! Atoms resolved:
            do j = i, Nat
               ! Get the column numberes in the file corresponding to these two elements:
               col1 = 2+(i-1)*Nat+(j-1)
               col2 = 2+(j-1)*Nat+(i-1)
               write(temp_col1,'(i0)') col1
               write(temp_col2,'(i0)') col2
               ! Also add the label and the linestyle:
               call select_linestyle(col1, linestyle) ! below

               if ((i /= Nat) .or. (j/=Nat)) then ! add coma
                  write(FN,'(a)') '('//trim(adjustl(temp_col1))//','//trim(adjustl(temp_col2))//', "'// &
                     trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//'", '// &
                     trim(adjustl(linestyle))//'),'
               else ! no coma
                  write(FN,'(a)') '('//trim(adjustl(temp_col1))// ','//trim(adjustl(temp_col2))//', "'// &
                     trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//'", '// &
                     trim(adjustl(linestyle))//')'
               endif
            enddo
         enddo
         write(FN, '(a)')  ']'

         write(FN, '(a)')  '# Create N distinct colors from a colormap'
         write(FN, '(a)')  'colors = plt.cm.tab20(np.linspace(0, 1, 20))'
         write(FN, '(a)')  'for idx, (c1, c2, label, ls) in enumerate(pairs):'
         write(FN, '(a)')  '     # Skip reversed duplicates (keep only c1 >= c2)'
         write(FN, '(a)')  '     if c1 > c2:'
         write(FN, '(a)')  '           continue'
         write(FN, '(a)')  '     if c1 == c2:'
         write(FN, '(a)')  '           # Plot column i alone'
         write(FN, '(a)')  '           summed = df.iloc[:, c1]'
         write(FN, '(a)')  '     else:'
         write(FN, '(a)')  '           # Plot sum of columns i + j'
         write(FN, '(a)')  '           summed = df.iloc[:, c1] + df.iloc[:, c2]'
         write(FN, '(a)')  '     plt.plot(t, summed, color=colors[idx], linestyle=ls, label=label)'

         write(x_min, '(f12.3)') t0
         write(x_max, '(f12.3)') t_last
         ! Make sure the maximum value is adjusted but not larger than the given one:
         write(FN,'(a)') 'xmax = df.iloc[:, 0].max()'
         x_max = 'min(xmax, '//trim(adjustl(x_max))//')'


         write(FN, '(a)')  'plt.xlim('// trim(adjustl(x_min)) //','// trim(adjustl(x_max)) //')'
         write(FN, '(a)')  'plt.ylim(None,None)'
         write(FN, '(a)')  'plt.title("Electron-ion coupling by element")'
         N_cols = 1+Nat*(Nat+1)/2     ! estimate the number of columns in the legend
         if (N_cols < 8) then ! normal:
            write(temp_txt,'(a)') 'fontsize=12'
         elseif (N_cols < 17) then ! smaller
            write(temp_txt,'(a)') 'fontsize=10'
         elseif (N_cols < 21) then
            write(temp_txt,'(a)') 'fontsize=9'
         elseif (N_cols < 33) then
            write(temp_txt,'(a)') 'fontsize=10, ncol=2'
         elseif (N_cols < 41) then
            write(temp_txt,'(a)') 'fontsize=9, ncol=2'
         else
            write(temp_txt,'(a)') 'fontsize=9, ncol=3'
         endif
         write(FN, '(a)')  'plt.legend(loc="best", '//trim(adjustl(temp_txt))//')'
         write(FN, '(a)')  '# Save the plot in this format:'
         write(FN, '(a)')  'plt.savefig("'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'", dpi=300, bbox_inches="tight")'
         write(FN, '(a)')  'plt.close()'

         close(FN)
      endif


      ! 2) Create orbital-resolved script, only if there are more than one orbital:
      ! Find number of different orbital types:
      N_at = size(Scell%MDatoms) ! total number of atoms
      Nsiz = size(Scell%Ha,1) ! total number of orbitals
      norb = Nsiz/N_at ! orbitals per atom
      ! Find number of different orbital types:
      N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

      if (N_types > 1) then
         ! Py script file:
         File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name_shell))
         open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

         Plot_name = 'OUTPUT_coupling_parameter_by_orbital'
         Data_file_name = file_coupling
         if (present(convolved)) then
            if (convolved) then
               Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
               Data_file_name  = trim(adjustl(file_coupling(1:len(trim(adjustl(file_coupling)))-4)))//'_CONVOLVED.dat'
            endif
         endif

         ! Create the python script:
         write(FN, '(a)')  'import pandas as pd'
         write(FN, '(a)')  'import numpy as np'
         write(FN, '(a)')  'import matplotlib.pyplot as plt'
         write(FN, '(a)')  'import matplotlib.cm as cm'
         write(FN, '(a)')  '# Set the format:'
         write(FN, '(a)')  'plt.xlabel("Time (fs)", fontsize=14)'
         write(FN, '(a)')  'plt.ylabel("Coupling parameter (W/(m$^3$ K))", fontsize=14)'
         write(FN, '(a)')  'plt.xticks(fontsize=12)'
         write(FN, '(a)')  'plt.yticks(fontsize=12)'
         write(FN, '(a)')  'plt.grid(False)'
         write(FN, '(a)')  '# Set the axes:'
         write(FN, '(a)')  '# Read the output file:'
         write(FN, '(a)')  'df = pd.read_csv(r"'// trim(adjustl(Data_file_name))// &
                        '", sep=r'//"'\s+',"//'header=None, comment="#", skipinitialspace=True)'
         write(FN, '(a)')  '# Prepare the plot:'
         write(FN, '(a)')  '# -----------------------------'
         write(FN, '(a)')  '# Time axis'
         write(FN, '(a)')  '# -----------------------------'
         write(FN, '(a)')  't = df.iloc[:, 0]'

         ! To add up off-diagonal contributions:
         write(FN, '(a)')  'pairs = ['
         ! First ,add the total coupling:
         write(FN,'(a)') '(1,1, "Total", "-"),'
         i_start = Nat**2     ! to skip the lines with elements
         ! Now, all element-resolved ones:
         do i = 1, Nat        ! for all atoms i
            do i_orb = 1, N_types      ! and all orbital types of atom i
               orb_name1 = name_of_orbitals(norb, i_orb) ! module "Little_subroutines"

               do j = i, Nat     ! for all atoms j
                  do j_orb = 1, N_types   ! and all orbital types of atom j
                     orb_name2 = name_of_orbitals(norb, j_orb) ! module "Little_subroutines"

                     ! Get the column numberes in the file corresponding to these two elements:
                     col1 = (2+i_start) + (i-1)*Nat*N_types**2 + (i_orb-1)*Nat*N_types + (j-1)*N_types + (j_orb-1)
                     col2 = (2+i_start) + (j-1)*Nat*N_types**2 + (j_orb-1)*Nat*N_types + (i-1)*N_types + (i_orb-1)
                     write(temp_col1,'(i0)') col1
                     write(temp_col2,'(i0)') col2
                     ! Also add the label and the linestyle:
                     call select_linestyle((col1-i_start), linestyle) ! below

                     if ((i /= Nat) .or. (j/=Nat) .or. (i_orb/=N_types) .or. (j_orb/=N_types)) then ! add coma
                        write(FN,'(a)') '('//trim(adjustl(temp_col1))//','//trim(adjustl(temp_col2))//', "'// &
                           trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                           trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//'", '// &
                           trim(adjustl(linestyle))//'),'
                     else ! no coma
                        write(FN,'(a)') '('//trim(adjustl(temp_col1))// ','//trim(adjustl(temp_col2))//', "'// &
                           trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                           trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//'", '// &
                           trim(adjustl(linestyle))//')'
                     endif
                  enddo
               enddo
            enddo
         enddo
         write(FN, '(a)')  ']'

         write(FN, '(a)')  '# Create N distinct colors from a colormap'
         write(FN, '(a)')  'colors = plt.cm.tab20(np.linspace(0, 1, 20))'

         write(FN, '(a)')  'for idx, (c1, c2, label, ls) in enumerate(pairs):'
         write(FN, '(a)')  '     # Skip reversed duplicates (keep only c1 >= c2)'
         write(FN, '(a)')  '     if c1 > c2:'
         write(FN, '(a)')  '           continue'
         write(FN, '(a)')  '     if c1 == c2:'
         write(FN, '(a)')  '           # Plot column i alone'
         write(FN, '(a)')  '           summed = df.iloc[:, c1]'
         write(FN, '(a)')  '     else:'
         write(FN, '(a)')  '           # Plot sum of columns i + j'
         write(FN, '(a)')  '           summed = df.iloc[:, c1] + df.iloc[:, c2]'
         write(FN, '(a)')  '     plt.plot(t, summed, color=colors[idx % len(colors)], linestyle=ls, label=label)'
         write(x_min, '(f12.3)') t0
         write(x_max, '(f12.3)') t_last
         ! Make sure the maximum value is adjusted but not larger than the given one:
         write(FN,'(a)') 'xmax = df.iloc[:, 0].max()'
         x_max = 'min(xmax, '//trim(adjustl(x_max))//')'


         write(FN, '(a)')  'plt.xlim('// trim(adjustl(x_min)) //','// trim(adjustl(x_max)) //')'
         write(FN, '(a)')  'plt.ylim(None,None)'
         write(FN, '(a)')  'plt.title("Electron-ion coupling by orbital")'

         ! Make appropriate font size:
         N_cols = 1+Nat*N_types*(Nat+1)*(N_types+1)/4     ! estimate the number of columns in the legend
         if (N_cols < 8) then ! normal:
            write(temp_txt,'(a)') 'fontsize=12'
         elseif (N_cols < 17) then ! smaller
            write(temp_txt,'(a)') 'fontsize=10'
         elseif (N_cols < 21) then
            write(temp_txt,'(a)') 'fontsize=9'
         elseif (N_cols < 33) then
            write(temp_txt,'(a)') 'fontsize=10, ncol=2'
         elseif (N_cols < 41) then
            write(temp_txt,'(a)') 'fontsize=9, ncol=2'
         else
            write(temp_txt,'(a)') 'fontsize=9, ncol=3'
         endif
         write(FN, '(a)')  'plt.legend(loc="best", '//trim(adjustl(temp_txt))//')'
         write(FN, '(a)')  '# Save the plot in this format:'
         write(FN, '(a)')  'plt.savefig("'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'", dpi=300, bbox_inches="tight")'
         write(FN, '(a)')  'plt.close()'

         close(FN)
      endif

   end select
end subroutine Python_plot_partial_coupling



subroutine Python_plot_fragments_data(Scell, numpar, file_fragments, t0, t_last, script_name_Nat, script_name_Ta, script_name_Ne) ! below
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_fragments, script_name_Nat, script_name_Ta ! file with atomic numbers and temperatures
   character(*), intent(in) :: script_name_Ne   ! file with electronic numbers
   !logical, intent(in), optional :: convolved   ! is it a convolved copy of files (not implemented)
   !----------------
   character(300) :: File_name, Plot_name, Data_file_name, temp_txt
   integer :: FN, N_cols
   logical :: hide_legend

   !----------------
   ! 1) Number of atoms in each fragment:
   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name_Nat))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Plot_name = 'OUTPUT_fragments_Nat'
   Data_file_name = trim(adjustl(file_fragments))


   write(FN, '(a)') 'import glob'
   write(FN, '(a)') 'import re'
   write(FN, '(a)') 'import matplotlib.pyplot as plt'

   write(FN, '(a)') 'def numeric_key(filename):'
   write(FN, '(a)') '    m = re.search(r"'// trim(adjustl(Data_file_name))//'(\d+)\.dat", filename)'
   write(FN, '(a)') '    return int(m.group(1)) if m else 0'
   write(FN, '(a)') 'files = sorted(glob.glob("'//trim(adjustl(Data_file_name))//'*.dat"), key=numeric_key)'

   write(FN, '(a)') 'for i, fname in enumerate(files, start=1):'
   write(FN, '(a)') '    times = []'
   write(FN, '(a)') '    natoms = []'

   write(FN, '(a)') '    with open(fname, "r") as f:'
   write(FN, '(a)') '        for line in f:'
   write(FN, '(a)') '            if line.lstrip().startswith("#"):'
   write(FN, '(a)') '                continue'
   write(FN, '(a)') '            if not line.strip():'
   write(FN, '(a)') '                continue'

   write(FN, '(a)') '            cols = line.split()'

   write(FN, '(a)') '            # col1 = time, col2 = natoms'
   write(FN, '(a)') '            t = float(cols[0])'
   write(FN, '(a)') '            n = float(cols[1])'

   write(FN, '(a)') '            times.append(t)'
   write(FN, '(a)') '            natoms.append(n)'

   write(FN, '(a)') '    plt.plot(times, natoms, label=f"Fragment #{i}")'

   write(FN, '(a)') 'plt.xlabel("Time (fs)", fontsize=14)'
   write(FN, '(a)') 'plt.ylabel("Number of atoms", fontsize=14)'

   write(FN, '(a)') 'plt.xlim(None,None)'
   write(FN, '(a)') 'plt.ylim(0,None)'

   write(FN, '(a)') 'plt.title("Fragment evolution", fontsize=14)'

   write(FN, '(a)') 'plt.xticks(fontsize=12)'
   write(FN, '(a)') 'plt.yticks(fontsize=12)'

   ! Make appropriate font size:
   hide_legend = .false.      ! to start with
   N_cols = Scell%fragments%N_frag_max
   if (N_cols < 8) then ! normal:
      write(temp_txt,'(a)') 'fontsize=12'
   elseif (N_cols < 17) then ! smaller
      write(temp_txt,'(a)') 'fontsize=10'
   elseif (N_cols < 21) then
      write(temp_txt,'(a)') 'fontsize=9'
   elseif (N_cols < 33) then
      write(temp_txt,'(a)') 'fontsize=10, ncol=2'
   elseif (N_cols < 41) then
      write(temp_txt,'(a)') 'fontsize=9, ncol=2'
   elseif (N_cols < 65) then
      write(temp_txt,'(a)') 'fontsize=9, ncol=3'
   else ! too many curves, no need for a legend at all
      hide_legend = .true.
   endif

   if (.not.hide_legend) then
      write(FN, '(a)') 'plt.legend(loc="best", '//trim(adjustl(temp_txt))//')'
   else ! hide the legend
      write(FN, '(a)') 'plt.legend().set_visible(False)'
   endif
   write(FN, '(a)') 'plt.tight_layout()'

   write(FN, '(a)') 'plt.savefig("'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'", dpi=300, bbox_inches="tight")'
   close(FN)

   !----------------
   ! 2) Temperatures of atoms in each fragment:
   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name_Ta))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Plot_name = 'OUTPUT_fragments_Ta'
   Data_file_name = trim(adjustl(file_fragments))

   write(FN, '(a)') 'import glob'
   write(FN, '(a)') 'import re'
   write(FN, '(a)') 'import matplotlib.pyplot as plt'

   write(FN, '(a)') 'def numeric_key(filename):'
   write(FN, '(a)') '    m = re.search(r"'// trim(adjustl(Data_file_name))//'(\d+)\.dat", filename)'
   write(FN, '(a)') '    return int(m.group(1)) if m else 0'
   write(FN, '(a)') 'files = sorted(glob.glob("'//trim(adjustl(Data_file_name))//'*.dat"), key=numeric_key)'

   write(FN, '(a)') '# Define a set of linestyles to cycle through'
   write(FN, '(a)') 'linestyles = ['
   write(FN, '(a)') '    "-",'
   write(FN, '(a)') '    "--"'
   write(FN, '(a)') ']'

   write(FN, '(a)') "colors = plt.rcParams['axes.prop_cycle'].by_key()['color']"

   write(FN, '(a)') 'for i, fname in enumerate(files, start=1):'
   write(FN, '(a)') '    times = []'
   write(FN, '(a)') '    Tkin = []'
   write(FN, '(a)') '    Tfluc = []'

   write(FN, '(a)') '    with open(fname, "r") as f:'
   write(FN, '(a)') '        for line in f:'
   write(FN, '(a)') '            if line.lstrip().startswith("#"):'
   write(FN, '(a)') '                continue'
   write(FN, '(a)') '            if not line.strip():'
   write(FN, '(a)') '                continue'

   write(FN, '(a)') '            cols = line.split()'

   write(FN, '(a)') '            t = float(cols[0])'
   write(FN, '(a)') '            n = float(cols[2])'
   write(FN, '(a)') '            m = float(cols[3])'

   write(FN, '(a)') '            times.append(t)'
   write(FN, '(a)') '            Tkin.append(n)'
   write(FN, '(a)') '            Tfluc.append(m)'


   write(FN, '(a)') '    # Pick one color per fragment, cycle if needed'
   write(FN, '(a)') '    color = colors[(i-1) % len(colors)]'

   write(FN, '(a)') '    # Assign two different linestyles per fragment'
   write(FN, '(a)') '    ls1 = linestyles[(2*(i-1)) % len(linestyles)]'
   write(FN, '(a)') '    ls2 = linestyles[(2*(i-1)+1) % len(linestyles)]'

   write(FN, '(a)') '    plt.plot(times, Tkin,  linestyle=ls1, color=color, label=f"Fragment #{i} T$_{{kin}}$")'
   write(FN, '(a)') '    plt.plot(times, Tfluc, linestyle=ls2, color=color, label=f"Fragment #{i} T$_{{fluc}}$")'

   write(FN, '(a)') 'plt.xlabel("Time (fs)", fontsize=14)'
   write(FN, '(a)') 'plt.ylabel("Temperature (K)", fontsize=14)'

   write(FN, '(a)') 'plt.xlim(None,None)'
   write(FN, '(a)') 'plt.ylim(0,None)'

   write(FN, '(a)') 'plt.title("Fragment temperatures", fontsize=14)'

   write(FN, '(a)') 'plt.xticks(fontsize=12)'
   write(FN, '(a)') 'plt.yticks(fontsize=12)'

   ! Make appropriate font size:
   hide_legend = .false.      ! to start with
   N_cols = Scell%fragments%N_frag_max*2
   if (N_cols < 8) then ! normal:
      write(temp_txt,'(a)') 'fontsize=12'
   elseif (N_cols < 17) then ! smaller
      write(temp_txt,'(a)') 'fontsize=10'
   elseif (N_cols < 21) then
      write(temp_txt,'(a)') 'fontsize=9'
   elseif (N_cols < 33) then
      write(temp_txt,'(a)') 'fontsize=10, ncol=2'
   elseif (N_cols < 41) then
      write(temp_txt,'(a)') 'fontsize=9, ncol=2'
   elseif (N_cols < 65) then
      write(temp_txt,'(a)') 'fontsize=9, ncol=3'
   else ! too many curves, no need for a legend at all
      hide_legend = .true.
   endif

   if (.not.hide_legend) then
      write(FN, '(a)') 'plt.legend(loc="best", '//trim(adjustl(temp_txt))//')'
   else ! hide the legend
      write(FN, '(a)') 'plt.legend().set_visible(False)'
   endif
   write(FN, '(a)') 'plt.tight_layout()'

   write(FN, '(a)') 'plt.savefig("'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'", dpi=300, bbox_inches="tight")'

   close(FN)



   !----------------
   ! 3) Number of electrons in each fragment:
   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name_Ne))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Plot_name = 'OUTPUT_fragments_Ne'
   Data_file_name = trim(adjustl(file_fragments))


   write(FN, '(a)') 'import glob'
   write(FN, '(a)') 'import re'
   write(FN, '(a)') 'import matplotlib.pyplot as plt'

   write(FN, '(a)') 'def numeric_key(filename):'
   write(FN, '(a)') '    m = re.search(r"'// trim(adjustl(Data_file_name))//'(\d+)\.dat", filename)'
   write(FN, '(a)') '    return int(m.group(1)) if m else 0'
   write(FN, '(a)') 'files = sorted(glob.glob("'//trim(adjustl(Data_file_name))//'*.dat"), key=numeric_key)'

   write(FN, '(a)') 'for i, fname in enumerate(files, start=1):'
   write(FN, '(a)') '    times = []'
   write(FN, '(a)') '    natoms = []'

   write(FN, '(a)') '    with open(fname, "r") as f:'
   write(FN, '(a)') '        for line in f:'
   write(FN, '(a)') '            if line.lstrip().startswith("#"):'
   write(FN, '(a)') '                continue'
   write(FN, '(a)') '            if not line.strip():'
   write(FN, '(a)') '                continue'

   write(FN, '(a)') '            cols = line.split()'

   write(FN, '(a)') '            # col1 = time, col2 = natoms'
   write(FN, '(a)') '            t = float(cols[0])'
   write(FN, '(a)') '            n = float(cols[4])'

   write(FN, '(a)') '            times.append(t)'
   write(FN, '(a)') '            natoms.append(n)'

   write(FN, '(a)') '    plt.plot(times, natoms, label=f"Fragment #{i}")'

   write(FN, '(a)') 'plt.xlabel("Time (fs)", fontsize=14)'
   write(FN, '(a)') 'plt.ylabel("Number of electrons", fontsize=14)'

   write(FN, '(a)') 'plt.xlim(None,None)'
   write(FN, '(a)') 'plt.ylim(0,None)'

   write(FN, '(a)') 'plt.title("Fragment electrons", fontsize=14)'

   write(FN, '(a)') 'plt.xticks(fontsize=12)'
   write(FN, '(a)') 'plt.yticks(fontsize=12)'

   ! Make appropriate font size:
   hide_legend = .false.      ! to start with
   N_cols = Scell%fragments%N_frag_max
   if (N_cols < 8) then ! normal:
      write(temp_txt,'(a)') 'fontsize=12'
   elseif (N_cols < 17) then ! smaller
      write(temp_txt,'(a)') 'fontsize=10'
   elseif (N_cols < 21) then
      write(temp_txt,'(a)') 'fontsize=9'
   elseif (N_cols < 33) then
      write(temp_txt,'(a)') 'fontsize=10, ncol=2'
   elseif (N_cols < 41) then
      write(temp_txt,'(a)') 'fontsize=9, ncol=2'
   elseif (N_cols < 65) then
      write(temp_txt,'(a)') 'fontsize=9, ncol=3'
   else ! too many curves, no need for a legend at all
      hide_legend = .true.
   endif

   if (.not.hide_legend) then
      write(FN, '(a)') 'plt.legend(loc="best", '//trim(adjustl(temp_txt))//')'
   else ! hide the legend
      write(FN, '(a)') 'plt.legend().set_visible(False)'
   endif
   write(FN, '(a)') 'plt.tight_layout()'

   write(FN, '(a)') 'plt.savefig("'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'", dpi=300, bbox_inches="tight")'
   close(FN)


end subroutine Python_plot_fragments_data



subroutine Python_plot_coupling(numpar, file_electron_properties, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_properties, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 1

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))


   col_nums(1) = 5
   col_labels(1) = '"Coupling"'


   Plot_name = 'OUTPUT_coupling_parameter'
   Data_file_name = file_electron_properties
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat'
      endif
   endif


   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Coupling parameter (W/(m$^3$ K))', 'Electron-ion coupling', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_coupling






subroutine Python_plot_at_temperatures_part(numpar, file_atomic_temperatures_part, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_atomic_temperatures_part, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 6

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))


   col_nums(1) = 1
   col_labels(1) = '"Kinetic: X"'
   linestyle(1) = '"--"'
   col_nums(2) = 2
   col_labels(2) = '"Kinetic: Y"'
   linestyle(2) = '(0,(5,2,1,2,1,2))'
   col_nums(3) = 3
   col_labels(3) = '"Kinetic: Z"'
   linestyle(3) = '"-"'

   col_nums(4) = 4
   col_labels(4) = '"Virial: X"'
   linestyle(4) = '"--"'
   col_nums(5) = 5
   col_labels(5) = '"Virial: Y"'
   linestyle(5) = '(0,(5,2,1,2,1,2))'
   col_nums(6) = 6
   col_labels(6) = '"Virial: Z"'
   linestyle(6) = '"-"'

   Plot_name = 'OUTPUT_atomic_temperatures_partial'
   Data_file_name = file_atomic_temperatures_part
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_atomic_temperatures_part(1:len(trim(adjustl(file_atomic_temperatures_part)))-4)))//'_CONVOLVED.dat'
      endif
   endif


   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Atomic temperature (K)', 'Partial atomic temperatures', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_at_temperatures_part




subroutine Python_plot_at_temperatures(numpar, file_atomic_temperatures, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_atomic_temperatures, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 4

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))


   col_nums(1) = 10
   col_labels(1) = '"Configurational"'
   linestyle(1) = '"--"'
   col_nums(2) = 11
   col_labels(2) = '"Hyperconfig"'
   linestyle(2) = '(0,(5,2,1,2,1,2))'
   col_nums(3) = 1
   col_labels(3) = '"Kinetic"'
   linestyle(3) = '"-"'
   col_nums(4) = 4
   col_labels(4) = '"Fluctuational"'
   linestyle(4) = '"-."'

   Plot_name = 'OUTPUT_atomic_temperatures'
   Data_file_name = file_atomic_temperatures
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_atomic_temperatures(1:len(trim(adjustl(file_atomic_temperatures)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Atomic temperature (K)', 'Atomic temperatures', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_at_temperatures





subroutine Python_plot_entropy_atomic(numpar, file_atomic_entropy, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_atomic_entropy, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 5

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))


   col_nums(1) = 2
   col_labels(1) = '"Equilibrium"'
   linestyle(1) = '"--"'
   col_nums(2) = 3
   col_labels(2) = '"Equilibrium (num)"'
   linestyle(2) = '":"'
   col_nums(3) = 4
   col_labels(3) = '"Configurational"'
   linestyle(3) = '"-."'
   col_nums(4) = 5
   col_labels(4) = '"Total"'
   linestyle(4) = '"-"'
   col_nums(5) = 1
   col_labels(5) = '"Nonequilibrium"'
   linestyle(5) = '(0,(5,2,1,2,1,2))'


   Plot_name = 'OUTPUT_atomic_entropy'
   Data_file_name = file_atomic_entropy
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_atomic_entropy(1:len(trim(adjustl(file_atomic_entropy)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Atomic entropy (eV/K)', 'Atomic entropy', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_entropy_atomic




subroutine Python_plot_chempots(numpar, file_mu, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_mu, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 3

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 1
   col_labels(1) = '"Total"'
   linestyle(1) = '"-"'
   col_nums(2) = 2
   col_labels(2) = '"Valence"'
   linestyle(2) = '"--"'
   col_nums(3) = 3
   col_labels(3) = '"Conduction"'
   linestyle(3) = '"-."'


   Plot_name = 'OUTPUT_electron_chempotentials'
   Data_file_name = file_mu
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_mu(1:len(trim(adjustl(file_mu)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Chemical potentials (eV)', 'Electron chemical potentials', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_chempots



subroutine Python_plot_el_temperatures(numpar, file_Te, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_Te, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   Nsiz = 3

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))

   col_nums(1) = 1
   col_labels(1) = '"Total (kinetic)"'
   linestyle(1) = '"-"'
   col_nums(2) = 2
   col_labels(2) = '"Valence"'
   linestyle(2) = '"--"'
   col_nums(3) = 3
   col_labels(3) = '"Conduction"'
   linestyle(3) = '"-."'

   Plot_name = 'OUTPUT_electron_temperatures'
   Data_file_name = file_Te
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_Te(1:len(trim(adjustl(file_Te)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Electron tempereature (K)', 'Electron tempereatures', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_el_temperatures




subroutine Python_plot_entropy(numpar, file_electron_entropy, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_entropy, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, Nsiz
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   if (numpar%do_partial_thermal) then ! partial entropies too
      Nsiz = 6
   else
      Nsiz = 2
   endif

   allocate(col_nums(Nsiz), source = 0)
   allocate(col_labels(Nsiz))
   allocate(linestyle(Nsiz))


   col_nums(1) = 2
   col_labels(1) = '"Equilibrium"'
   linestyle(1) = '"--"'
   col_nums(2) = 1
   col_labels(2) = '"Nonequilibrium"'
   linestyle(2) = '"-"'
   if (numpar%do_partial_thermal) then ! add partial entropies
      col_nums(3) = 4
      col_labels(3) = '"Equilibrium VB"'
      linestyle(3) = '"-."'
      col_nums(4) = 3
      col_labels(4) = '"Nonequilibrium VB"'
      linestyle(4) = '"-."'
      col_nums(5) = 6
      col_labels(5) = '"Equilibrium CB"'
      linestyle(5) = '":"'
      col_nums(6) = 5
      col_labels(6) = '"Nonequilibrium CB"'
      linestyle(6) = '":"'
   endif


   Plot_name = 'OUTPUT_electron_entropy'
   Data_file_name = file_electron_entropy
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_entropy(1:len(trim(adjustl(file_electron_entropy)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Electron entropy (eV/K)', 'Electron entropy', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_entropy



subroutine Python_plot_heat_conductivity_dyn(numpar, file_heat_conductivity_dyn, t0, t_last, script_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_heat_conductivity_dyn, script_name ! file with energy levels, script
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(3), source = 0)
   allocate(col_labels(3))
   allocate(linestyle(3))

   col_nums(1) = 1
   col_labels(1) = 'r"$\kappa$"'
   linestyle(1) = '"-"'
   col_nums(2) = 2
   col_labels(2) = 'r"$\kappa_{e-ph}$"'
   linestyle(2) = '"--"'
   col_nums(3) = 3
   col_labels(3) = 'r"$\kappa_{e-e}$"'
   linestyle(3) = '"-."'

   call Create_python_plot(FN, file_heat_conductivity_dyn, col_nums, col_labels, &
      'Time (fs)', 'Heat conductivity (W/(m K))', 'Electron heat conductivity', &
      "best", 'OUTPUT_electron_heat_conductivity_dyn', trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_heat_conductivity_dyn



subroutine Python_plot_heat_conductivity(numpar, file_heat_conductivity, t0, t_last, script_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_heat_conductivity, script_name ! file with energy levels, script
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(3), source = 0)
   allocate(col_labels(3))
   allocate(linestyle(3))

   col_nums(1) = 1
   col_labels(1) = 'r"$\kappa$"'
   linestyle(1) = '"-"'
   col_nums(2) = 2
   col_labels(2) = 'r"$\kappa_{e-ph}$"'
   linestyle(2) = '"--"'
   col_nums(3) = 3
   col_labels(3) = 'r"$\kappa_{e-e}$"'
   linestyle(3) = '"-."'

   call Create_python_plot(FN, file_heat_conductivity, col_nums, col_labels, &
      'Time (fs)', 'Heat conductivity (W/(m K))', 'Electron heat conductivity', &
      "best", 'OUTPUT_electron_heat_conductivity', trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_heat_conductivity




subroutine Python_plot_capacity(numpar, file_electron_properties, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_properties, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(1), source = 0)
   allocate(col_labels(1))

   col_nums(1) = 4
   col_labels(1) = '"C$_e$"'


   Plot_name = 'OUTPUT_electron_Ce'
   Data_file_name = file_electron_properties
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Heat capacity (J/(m$^3$ K))', 'Electron heat capacity', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_capacity



subroutine Python_plot_Ebands(numpar, file_electron_properties, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_properties, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name


   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(5), source = 0)
   allocate(col_labels(5))
   allocate(linestyle(5))

   col_nums(1) = 9
   col_labels(1) = '"Top of CB"'
   col_nums(2) = 8
   col_labels(2) = '"Bottom of CB"'
   col_nums(3) = 7
   col_labels(3) = '"Top of VB"'
   col_nums(4) = 6
   col_labels(4) = '"Bottom of VB"'
   col_nums(5) = 2
   col_labels(5) = '"Chemical potential"'

   linestyle = '"-"'          ! all solid
   linestyle(5) = '"--"'      ! chem.pot dashed

   Plot_name = 'OUTPUT_bands'
   Data_file_name = file_electron_properties
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Energy (eV)', 'Band borders', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_Ebands




subroutine Python_plot_mu(numpar, file_electron_properties, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_properties, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !-------------
   character(300) :: File_name, Plot_name, Data_file_name
   character(30) :: x_min_txt, x_max_txt
   integer :: FN
   !-------------

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   write(x_min_txt,'(f24.8)') t0
   write(x_max_txt,'(f24.8)') t_last


   Plot_name = 'OUTPUT_mu_and_Ne'
   Data_file_name = file_electron_properties
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   ! Two Y-axes require special treatement here:
   write(FN, '(a)')  'import numpy as np'
   write(FN, '(a)')  'import pandas as pd'
   write(FN, '(a)')  'import matplotlib.pyplot as plt'

   write(FN, '(a)')  '# --- Load data ---'
   write(FN, '(a)')  'df = pd.read_csv('
   write(FN, '(a)')  '"'//trim(adjustl(Data_file_name))//'",'
   write(FN, '(a)')  'sep=r"\s+",'
   write(FN, '(a)')  'header=None,'
   write(FN, '(a)')  'comment="#"'
   write(FN, '(a)')  ')'

   write(FN, '(a)') 'time = df.iloc[:, 0]'
   write(FN, '(a)') 'Ne   = df.iloc[:, 1]   # electron density (y2 axis)'
   write(FN, '(a)') 'mu   = df.iloc[:, 2]   # chemical potential (left axis)'

   ! Make sure the maximum X value is adjusted but not larger than the given one:
   write(FN,'(a)') 'xmax = df.iloc[:, 0].max()'
   x_max_txt = 'min(xmax, '//trim(adjustl(x_max_txt))//')'

   write(FN, '(a)') '# --- Create figure ---'
   write(FN, '(a)') 'fig, ax1 = plt.subplots(figsize=(8, 6), dpi=150)'
   write(FN, '(a)') 'ax1.set_xlim('//trim(adjustl(x_min_txt))//','//trim(adjustl(x_max_txt))//')'

   write(FN, '(a)') 'LW = 2.0'

   write(FN, '(a)') '# --- Left axis: chemical potential ---'
   write(FN, '(a)') 'line1, = ax1.plot('
   write(FN, '(a)') '   time, mu,'
   write(FN, '(a)') '   lw=LW,'
   write(FN, '(a)') '   color="black",'
   write(FN, '(a)') '   label="Chemical potential"'
   write(FN, '(a)') ')'

   write(FN, '(a)') 'ax1.set_xlabel("Time (fs)", fontsize=16)'
   write(FN, '(a)') 'ax1.set_ylabel("Energy (eV)", fontsize=16)'
   write(FN, '(a)') 'ax1.tick_params(axis="both", labelsize=12)'

   write(FN, '(a)') '# --- Right axis: electron density ---'
   write(FN, '(a)') 'ax2 = ax1.twinx()'

   write(FN, '(a)') '# --- Compute padded limits for right axis ---'
   write(FN, '(a)') 'y2_min = Ne.min()'
   write(FN, '(a)') 'y2_max = Ne.max()'
   write(FN, '(a)') 'padding = max(0.1 * (y2_max - y2_min),0.1)'
   write(FN, '(a)') 'ax2.set_ylim(y2_min - padding, y2_max + padding)'

   write(FN, '(a)') 'line2, = ax2.plot('
   write(FN, '(a)') '   time, Ne,'
   write(FN, '(a)') '   lw=LW,'
   write(FN, '(a)') '   linestyle="--",'
   write(FN, '(a)') '   color="tab:blue",'
   write(FN, '(a)') '   label="Electron density"'
   write(FN, '(a)') ')'

   write(FN, '(a)') 'ax2.set_ylabel("Electron density (1/atom)", fontsize=16, color="tab:blue")'
   write(FN, '(a)') 'ax2.tick_params(axis="y", labelcolor="tab:blue", labelsize=12)'

   write(FN, '(a)') '# y2 ticks every 0.1 (as in gnuplot)'
   write(FN, '(a)') 'ax2.set_yticks(np.arange(ax2.get_ylim()[0], ax2.get_ylim()[1], 0.1))'

   write(FN, '(a)') '# --- Legend (combined) ---'
   write(FN, '(a)') 'lines = [line1, line2]'
   write(FN, '(a)') 'labels = [l.get_label() for l in lines]'
   write(FN, '(a)') 'ax1.legend(lines, labels, loc="best", fontsize=14)'

   write(FN, '(a)') 'plt.title("Chemical potential and electron density", fontsize=14)'
   write(FN, '(a)') 'plt.tight_layout()'

   write(FN, '(a)') 'plt.savefig("'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'", dpi=300, bbox_inches="tight")'
   write(FN, '(a)') 'plt.close()'

   close(FN)
end subroutine Python_plot_mu




subroutine Python_plot_Egap(numpar, file_electron_properties, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_electron_properties, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name
   character(20) :: chtemp1

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(1), source = 0)
   allocate(col_labels(1))

   ! All shells resolved:
   col_nums(1) = 3
   col_labels(1) = '"Band gap"'

   Plot_name = 'OUTPUT_Egap'
   Data_file_name = file_electron_properties
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Band gap (eV)', 'Band gap', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_Egap




subroutine Python_plot_holes(matter, numpar, file_deep_holes, t0, t_last, script_name, convolved)
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_deep_holes, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, j, Na, N_sh_max, Nshl, N_sh_tot, i_cur
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name
   character(20) :: chtemp1, chtemp11

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find how many shells we have for plotting:
   Na = size(matter%Atoms)
   i_cur = 0 ! to start with
   do i = 1, Na !size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl ! for all shells of this atom
         i_cur = i_cur + 1    ! next shell
      enddo
   enddo
   N_sh_tot = i_cur - 1 ! estimate the total amount of core shells (exclude valence band)

   allocate(col_nums(N_sh_tot), source = 0)
   allocate(col_labels(N_sh_tot))
   allocate(linestyle(N_sh_tot))

   ! All shells:
   i_cur = 1      ! to start with
   ATOMS:do i = 1, Na ! size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)     ! number of shells in this element
      SHELLS:do j = 1, Nshl ! for all shells of this atom
         if ((i .NE. 1) .or. (j .NE. Nshl)) then ! atomic shell:
            call define_PQN(matter%Atoms(i)%Shl_dsgnr(j), chtemp11) ! module "Dealing_with_EADL"
            write(col_labels(i_cur),'(a)') '"'//trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(chtemp11))//'"'

            col_nums(i_cur) = i_cur    ! column number

            !set different style for different elements:
            call select_linestyle(i, linestyle(i_cur))      ! below

            i_cur = i_cur + 1    ! next shell
         else ! Valence band
            ! Skip the valence band
         endif
      enddo SHELLS
   enddo ATOMS


   Plot_name = 'OUTPUT_deep_shell_holes'
   Data_file_name = file_deep_holes
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_deep_holes(1:len(trim(adjustl(file_deep_holes)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Number of holes (total)', 'Core holes', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_holes



subroutine Python_plot_CB_electrons(numpar, file_numbers, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_numbers, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(2), source = 0)
   allocate(col_labels(2))
   allocate(linestyle(2))

   ! All shells resolved:
   col_nums(1) = 2
   col_labels(1) = '"Electrons"'
   linestyle(1) = '"-"'
   col_nums(2) = 6
   col_labels(2) = '"Photons"'
   linestyle(2) = '"--"'


   Plot_name = 'OUTPUT_CB_electrons'
   Data_file_name = file_numbers
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Electrons (1/atom)', 'Conduction-band electrons', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_CB_electrons




subroutine Python_plot_orbital_resolved(Scell, matter, numpar, file_numbers, t0, t_last, script_name, convolved)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_numbers, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_col, NKOA, N_at, Nsiz, norb, N_types, N_col, i_at, i_types
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name
   character(20) :: chtemp1

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find number of orbitals per atom:
   NKOA = matter%N_KAO      ! number of kinds of atoms
   N_at = size(Scell%MDatoms) ! number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   norb = Nsiz/N_at ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

   N_col = NKOA * N_types

   allocate(col_nums(N_col+1), source = 0)
   allocate(col_labels(N_col+1))
   allocate(linestyle(N_col+1))

   ! All shells resolved:
   col_nums(1) = 1
   col_labels(1) = '"Total"'
   linestyle(1) = '"-"'
   i_col = 1   ! to start with
   do i_at = 1, NKOA
      do i_types = 1, N_types
         ! Get name of the orbital:
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
         i_col = i_col + 1 ! column number

         col_nums(i_col) = i_col
         col_labels(i_col) = '"'//trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))//'"'

         !set different style for different elements:
         select case (i_at)   ! currently, supports 7 different elements (may be added if needed)
         case (1)    ! dashed
            linestyle(i_col) = '"--"'
         case (2)    ! dash-dot
            linestyle(i_col) = '"-."'
         case (3)    ! dash-dot-dot
            linestyle(i_col) = '(0,(5,2,1,2,1,2))'
         case (4)    ! long dash-short dash-dot
            linestyle(i_col) = '(0,(10,3,4,3,1,2))'
         case (5)    ! dot
            linestyle(i_col) = '":"'
         case default ! Long dash - dot - short dash - dot
            linestyle(i_col) = '(0,(10,3,1,2,4,3,1,2))'
         end select
      enddo   ! i_types
   enddo ! i_at


   Plot_name = 'OUTPUT_orbital_resolved_Ne'
   Data_file_name = file_numbers
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Electrons (1/atom)', 'Orbital-resolved electrons', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_orbital_resolved




subroutine Python_plot_numbers(numpar, file_numbers, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_numbers, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(4), source = 0)
   allocate(col_labels(4))
   allocate(linestyle(4))

   col_nums(1) = 3
   col_labels(1) = '"High-energy electrons"'
   linestyle(1) = '"-"'
   col_nums(2) = 4
   col_labels(2) = '"All deep-shell holes"'
   linestyle(2) = '"-."'
   col_nums(3) = 6
   col_labels(3) = '"Photons"'
   linestyle(3) = '"--"'
   col_nums(4) = 5
   col_labels(4) = '"Conservation error"'
   linestyle(4) = '":"'


   Plot_name = 'OUTPUT_electrons_and_holes'
   Data_file_name = file_numbers
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Particles (1/atom)', 'Number of particles', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_numbers




subroutine Python_plot_high_energy_el(numpar, file_high_e, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_high_e, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(5), source = 0)
   allocate(col_labels(5))
   allocate(linestyle(5))

   ! Set columns description:
   do i = 1, 5
      col_nums(i) = i
      call select_linestyle(i, linestyle(i))    ! below
   enddo
   col_labels(1) = '"Total"'
   col_labels(2) = '"Photo"'
   col_labels(3) = '"Auger"'
   col_labels(4) = '"Impact-ionized"'
   col_labels(5) = '"Error"'


   Plot_name = 'OUTPUT_electrons_high_energy'
   Data_file_name = file_high_e
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_high_e(1:len(trim(adjustl(file_high_e)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Electrons (1/atom)', 'High-energy electrons', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_high_energy_el



subroutine Python_plot_stress(numpar, file_pressure, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_pressure, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(9), source = 0)
   allocate(col_labels(9))
   allocate(linestyle(9))

   do i = 1, 9
      col_nums(i) = 1+i
   enddo
   col_labels(1) = '"P$_{xx}$"'
   col_labels(2) = '"P$_{xy}$"'
   col_labels(3) = '"P$_{xz}$"'
   col_labels(4) = '"P$_{yx}$"'
   col_labels(5) = '"P$_{yy}$"'
   col_labels(6) = '"P$_{yz}$"'
   col_labels(7) = '"P$_{zx}$"'
   col_labels(8) = '"P$_{zy}$"'
   col_labels(9) = '"P$_{zz}$"'
   ! Dashed lines for off diagonal:
   linestyle = '"--"'
   ! Solid lines for diagonal:
   linestyle(1) = '"-"'
   linestyle(5) = '"-"'
   linestyle(9) = '"-"'


   Plot_name = 'OUTPUT_pressure_tensor'
   Data_file_name = file_pressure
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Pressure tensor (GPa)', 'Pressure tensor', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_stress




subroutine Python_plot_pressure(numpar, file_pressure, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_pressure, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(1), source = 0)
   allocate(col_labels(1))

   col_nums(1) = 1
   col_labels(1) = '"Pressure"'


   Plot_name = 'OUTPUT_pressure'
   Data_file_name = file_pressure
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Pressure (GPa)', 'Pressure', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_pressure





subroutine Python_plot_powder_diffraction(Scell, matter, numpar, file_powder_diffraction, script_name, video_format)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_powder_diffraction, script_name ! file with energy levels, script
   character(*), intent(in) :: video_format     ! which video format to print it out in
   !----------------
   integer :: FN, i, i_start, ind, j, k
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")


   ! Get the number of columns to print:
   ind = 1     ! top start counting columns

   if (matter%N_KAO > 1) then ! no need to do for 1 element
      do j = 1, matter%N_KAO     ! for all elements
         do k = j, matter%N_KAO  ! for all different pairs
            ind = ind + 1
         enddo
      enddo
   endif
   ! allocate the arrays:
   allocate(col_nums(ind), source = 0)
   allocate(col_labels(ind))
   !allocate(linestyle(size(Scell%diff_peaks%I_diff_peak)))
   ! Set the arrays:
   col_nums(1) = 1
   col_labels(1) = 'Total'
   ind = 1  ! restart
   if (matter%N_KAO > 1) then ! no need to do for 1 element
      do j = 1, matter%N_KAO     ! for all elements
         do k = j, matter%N_KAO  ! for all different pairs
            ind = ind + 1
            col_nums(ind) = ind
            col_labels(ind) = trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))
         enddo
      enddo
   endif

   call Create_Python_animation(FN, file_powder_diffraction, col_nums, col_labels, &
      '2theta (deg)', 'Peak Intensity (a.u.)', 'Powder diffraction', &
      "best", 'OUTPUT_diffraction_powder', trim(adjustl(video_format)), &
      x_min=10.0d0, x_max=180.0d0, t_start=numpar%t_start, dt=numpar%dt_save, t_end=numpar%t_total)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_powder_diffraction



subroutine Python_plot_Debye_temperatures(Scell, numpar, file_Debye_temperature, t0, t_last, script_name, convolved)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_Debye_temperature, script_name ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   allocate(col_nums(size(Scell%diff_peaks%I_diff_peak)), source = 0)
   allocate(col_labels(size(Scell%diff_peaks%I_diff_peak)))
   !allocate(linestyle(size(Scell%diff_peaks%I_diff_peak)))

   do i = 1, size(Scell%diff_peaks%I_diff_peak)
      col_nums(i) = i
      col_labels(i) = make_diff_peak_name(Scell, i) ! below
      col_labels(i) = '"'//trim(adjustl(col_labels(i)))//'"'
   enddo

   Plot_name = 'OUTPUT_Debye_temperatures'
   Data_file_name = file_Debye_temperature
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_Debye_temperature(1:len(trim(adjustl(file_Debye_temperature)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Debye temperature (K)', 'Debye temperatures', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_Debye_temperatures



subroutine Python_plot_diffraction_peaks(Scell, numpar, file_diffraction, t0, t_last, script_name, plot_name, plot_label, DW, convolved)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_diffraction, script_name, plot_name, plot_label ! file with energy levels, script, plot
   logical, intent(in) :: DW  ! if it's DW, change the axis label
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, y_axis_label, Plot_name_used, Data_file_name
   character(10) :: units, temp

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")


   if (DW) then
      y_axis_label = 'DW peak intensity (arb. units)'
   else
      y_axis_label = 'Peak intensity (arb. units)'
   endif

   allocate(col_nums(size(Scell%diff_peaks%I_diff_peak)), source = 0)
   allocate(col_labels(size(Scell%diff_peaks%I_diff_peak)))
   !allocate(linestyle(size(Scell%diff_peaks%I_diff_peak)))

   do i = 1, size(Scell%diff_peaks%I_diff_peak)
      col_nums(i) = i
      col_labels(i) = make_diff_peak_name(Scell, i) ! below
      col_labels(i) = '"'//trim(adjustl(col_labels(i)))//'"'
   enddo


   Plot_name_used = trim(adjustl(plot_name))
   Data_file_name = file_diffraction
   if (present(convolved)) then
      if (convolved) then
         Plot_name_used = trim(adjustl(Plot_name_used))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_diffraction(1:len(trim(adjustl(file_diffraction)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', trim(adjustl(y_axis_label)), trim(adjustl(plot_label)), &
      "best", trim(adjustl(Plot_name_used)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0)     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Python_plot_diffraction_peaks

function make_diff_peak_name(Scell, i) result(peak_name)
   character(30) :: peak_name
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: i   ! number of the peak
   !---------------
   character(5) :: text1, text2, text3
   write(text1, '(i0)') Scell%diff_peaks%ijk_diff_peak(1,i)
   write(text2, '(i0)') Scell%diff_peaks%ijk_diff_peak(2,i)
   write(text3, '(i0)') Scell%diff_peaks%ijk_diff_peak(3,i)
   peak_name = '('//trim(adjustl(text1))//trim(adjustl(text2))//trim(adjustl(text3))//')'
end function make_diff_peak_name



subroutine Python_plot_displacements_partial(matter, numpar, file_MSD, t0, t_last, script_name, MSD_power, mask_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_MSD, script_name  ! file with energy levels, script
   real(8), intent(in) :: MSD_power ! power of MSD
   character(*), intent(in) :: mask_name
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name
   character(10) :: units, temp

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")


   if (abs(MSD_power) > 1.0) then
      write(temp, '(i0)') int(MSD_power)
      write(units, '(a)') '(A$^'//trim(adjustl(temp))//'$)'
   else
      units = '(A)'
   endif

   allocate(col_nums(4*matter%N_KAO), source = 0)
   allocate(col_labels(4*matter%N_KAO))
   allocate(linestyle(4*matter%N_KAO))

   i_start = 5
   do i = 1, matter%N_KAO
      col_nums(1 + (i-1)*4 ) = i_start + (i-1)*4
      col_labels(1 + (i-1)*4) = '"'//trim(adjustl(matter%Atoms(i)%Name))//'"'
      col_nums(2+ (i-1)*4) = i_start + (i-1)*4 + 1
      col_labels(2+ (i-1)*4) = '"'//trim(adjustl(matter%Atoms(i)%Name))//':X"'
      col_nums(3+ (i-1)*4) = i_start + (i-1)*4 + 2
      col_labels(3+ (i-1)*4) = '"'//trim(adjustl(matter%Atoms(i)%Name))//':Y"'
      col_nums(4+ (i-1)*4) = i_start + (i-1)*4 + 3
      col_labels(4+ (i-1)*4) = '"'//trim(adjustl(matter%Atoms(i)%Name))//':Z"'
      linestyle(1+ (i-1)*4) = '"-"'
      linestyle(2+ (i-1)*4:4+ (i-1)*4) = '"--"'
   enddo


   Plot_name = 'OUTPUT_mean_displacement_'//trim(adjustl(mask_name))//'_partial'
   Data_file_name = file_MSD
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_MSD(1:len(trim(adjustl(file_MSD)))-4)))//'_CONVOLVED.dat'
      endif
   endif


   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Mean displacement '//trim(adjustl(units)), 'Mean displacement', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_displacements_partial



subroutine Python_plot_displacements(matter, numpar, file_MSD, t0, t_last, script_name, MSD_power, mask_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_MSD, script_name  ! file with energy levels, script
   real(8), intent(in) :: MSD_power ! power of MSD
   character(*), intent(in) :: mask_name
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name
   character(10) :: units, temp

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")


   if (abs(MSD_power) > 1) then
      write(temp, '(i0)') int(MSD_power)
      write(units, '(a)') '(A$^'//trim(adjustl(temp))//'$)'
   else
      units = '(A)'
   endif

   allocate(col_nums(4), source = 0)
   allocate(col_labels(4))
   allocate(linestyle(4))
   col_nums(1) = 1
   col_labels(1) = '"Average"'
   col_nums(2) = 2
   col_labels(2) = '"X"'
   col_nums(3) = 3
   col_labels(3) = '"Y"'
   col_nums(4) = 4
   col_labels(4) = '"Z"'
   linestyle(1) = '"-"'
   linestyle(2:4) = '"--"'


   Plot_name = 'OUTPUT_mean_displacement_'//trim(adjustl(mask_name))
   Data_file_name = file_MSD
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_MSD(1:len(trim(adjustl(file_MSD)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Mean displacement '//trim(adjustl(units)), 'Mean displacement', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_displacements



subroutine Python_plot_MSD(matter, numpar, file_MSD, t0, t_last, script_name, MSD_power, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_MSD, script_name  ! file with energy levels, script
   integer, intent(in) :: MSD_power ! power of MSD
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i, i_start
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name
   character(10) :: units, temp

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")


   if (abs(MSD_power) > 1) then
      write(temp, '(i)') MSD_power
      write(units, '(a)') '(A$^'//trim(adjustl(temp))//'$)'
   else
      units = '(A)'
   endif


   if (matter%N_KAO == 1) then      ! single spieces
      ! Prepare column indices and names:
      allocate(col_nums(1), source = 0)
      allocate(col_labels(1))
      allocate(linestyle(1))
      col_nums(1) = 4
      col_labels(1) = '"Displacement"'
      linestyle(1) = '"-"'
   else ! more than one element:
      allocate(col_nums(matter%N_KAO+1), source = 0)
      allocate(col_labels(matter%N_KAO+1))
      allocate(linestyle(matter%N_KAO+1))
      i_start = 3 + matter%N_KAO
      col_nums(1) = i_start
      col_labels(1) = '"Average"'
      linestyle(1) = '"-"'
      ! All elements:
      do i = 1, matter%N_KAO
         col_nums(1+i) = i_start+i
         col_labels(1+i) = '"'//trim(adjustl(matter%Atoms(i)%Name))//' atoms"'
         linestyle(1+i) = '"--"'
      enddo
   endif


   Plot_name = 'OUTPUT_mean_displacement'
   Data_file_name = file_MSD
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_MSD(1:len(trim(adjustl(file_MSD)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Mean displacement '//trim(adjustl(units)), 'Mean displacement', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_MSD




subroutine Python_plot_temperatures(numpar, matter, file_temperatures, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_temperatures, script_name  ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN, i
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")


   if (matter%N_KAO == 1) then      ! single spieces
      ! Prepare column indices and names:
      allocate(col_nums(2), source = 0)
      allocate(col_labels(2))
      allocate(linestyle(2))
      col_nums(1) = 1
      col_labels(1) = '"Electrons"'
      linestyle(1) = '"--"'
      col_nums(2) = 2
      col_labels(2) = '"Atoms"'
      linestyle(2) = '"-"'
   else ! more than one element:
      allocate(col_nums(matter%N_KAO+2), source = 0)
      allocate(col_labels(matter%N_KAO+2))
      allocate(linestyle(matter%N_KAO+2))
      col_nums(1) = 1
      col_labels(1) = '"Electrons"'
      linestyle(1) = '"--"'
      col_nums(2) = 2
      col_labels(2) = '"Average atoms"'
      linestyle(2) = '"-"'
      ! All elements:
      do i = 1, matter%N_KAO
         col_nums(i+2) = 2+i
         col_labels(i+2) = '"'//trim(adjustl(matter%Atoms(i)%Name))//' atoms"'
         linestyle(i+2) = '"-."'
      enddo
   endif


   Plot_name = 'OUTPUT_temepratures'
   Data_file_name = file_temperatures
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_temperatures(1:len(trim(adjustl(file_temperatures)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Temperature (K)', 'Temperatures', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, y_min=0.0d0, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_temperatures




subroutine Python_plot_energies(numpar, file_energies, t0, t_last, script_name, convolved)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_Energies, script_name  ! file with energy levels, script
   logical, intent(in), optional :: convolved   ! is it a convolved copy of files
   !----------------
   integer :: FN
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels, linestyle
   character(300) :: File_name, Plot_name, Data_file_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Prepare column indices and names:
   allocate(col_nums(4), source = 0)
   allocate(col_labels(4))
   allocate(linestyle(4))

   col_nums(1) = 7
   col_labels(1) = '"Total energy"'
   linestyle(1) = '"-"'
   col_nums(2) = 6
   col_labels(2) = '"Atoms and electrons"'
   linestyle(2) = '"-."'
   col_nums(3) = 5
   col_labels(3) = '"Atomic energy"'
   linestyle(3) = '":"'
   col_nums(4) = 3
   col_labels(4) = '"Potential energy"'
   linestyle(4) = '"--"'


   Plot_name = 'OUTPUT_energies'
   Data_file_name = file_energies
   if (present(convolved)) then
      if (convolved) then
         Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
         Data_file_name  = trim(adjustl(file_energies(1:len(trim(adjustl(file_energies)))-4)))//'_CONVOLVED.dat'
      endif
   endif

   call Create_python_plot(FN, trim(adjustl(Data_file_name)), col_nums, col_labels, &
      'Time (fs)', 'Energy (eV/atom)', 'Energies', &
      "best", trim(adjustl(Plot_name)), trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last, l_style=linestyle)     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Python_plot_energies




subroutine Python_plot_energy_levels(numpar, t0, t_last, Scell, file_Ei, script_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   character(*), intent(in) :: file_Ei, script_name  ! file with energy levels, script
   !----------------
   integer i, M, NSC, FN
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels
   character(300) :: File_name

   do NSC = 1, size(Scell)    ! for all supercells (which is one)

      ! Py script file:
      File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

      ! Prepare column indices:
      M = size(Scell(NSC)%Ei) ! number of levels
      ! Prepare column indices:
      allocate(col_nums(M), source = 0)
      do i = 1, M
         col_nums(i) = i
      enddo ! i

      ! Use them to plot the data:
      call Create_python_plot(FN, file_Ei, col_nums, col_labels, &
            'Time (fs)', 'Energy levels (eV)', 'Electron energy levels', &
            "best", 'OUTPUT_energy_levels', trim(adjustl(numpar%fig_extention)), &
            x_min=t0, x_max=t_last, y_max=25.0d0 )     ! below

      close(FN)
      deallocate(col_nums)
   enddo    ! NSC
end subroutine Python_plot_energy_levels





subroutine Plot_electron_MFP_python(numpar, matter, file_electron_IMFP, file_electron_EMFP)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   character(*), intent(in) :: file_electron_IMFP, file_electron_EMFP     ! file with data
   !-------------
   character(300) :: py_electron_MFP
   character(11) :: call_slash
   character(8) :: col
   integer :: N_grid, FN, count_col, i, Nshl, j_start, j, N_tot_col
   real(8) :: t0, t_last, x_tics, x_min, x_max, y_min, y_max
   integer, dimension(:), allocatable :: col_nums, col_nums2
   character(30), dimension(:), allocatable :: col_labels, col_labels2
   !-------------

   N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
   x_min = 1.0d0  ! [eV] starting energy point
   x_max = matter%Atoms(1)%Ph_MFP(1)%E(N_grid) ! ending energy point
   y_min = 1.0d0
   y_max = 1.0d5
   x_tics = 10.0d0

   ! Electron MFPs plotting:
   py_electron_MFP = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//'OUTPUT_MFP_electron.py'
   open(NEWUNIT=FN, FILE = trim(adjustl(py_electron_MFP)), action="write", status="replace")

   ! Prepare column indices and names:
   call count_shells(matter, N_tot_col)     ! below
   allocate(col_nums(N_tot_col+1), source = 0)
   allocate(col_labels(N_tot_col+1))
   ! Fill the column numbers and labels:
   count_col = 0  ! to start with
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl    ! for all shells of this atom
         if ((i == 1) .and. (j == Nshl)) then ! valence
            ! skip the column, it is not here but at the end
         else
            count_col = count_col + 1  ! number of columns
            col_nums(count_col) = count_col
            ! Names of the columns (skipping the "Time"):
            ! core shell
            col_labels(count_col) = '"'//trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))//'"'
         endif
      enddo
   enddo
   ! Add the valence band MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_labels(count_col) = '"Valence"'
   ! And the total MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_labels(count_col) = '"Total inelastic"'

   ! Add the elastic part:
   allocate(col_nums2(1), source = 0)
   allocate(col_labels2(1))
   col_nums2(1) = 1+size(matter%Atoms)
   col_labels2(1) = '"Elastic"'

   ! Use them to plot the data:
   call Create_python_plot(FN, file_electron_IMFP, col_nums, col_labels, &
      'Electron energy (eV)', 'Mean free path (A)', 'Electron mean free paths', &
      "best", 'OUTPUT_MFP_electron', trim(adjustl(numpar%fig_extention)), &
      x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, &
      set_x_log=.true., set_y_log=.true., &
      Data_file2=file_electron_EMFP, col_nums2=col_nums2, col_labels2=col_labels2, &
      colors_inverted=.true. )     ! below

   close(FN)
   deallocate(col_nums, col_labels, col_nums2, col_labels2)
end subroutine Plot_electron_MFP_python




subroutine Plot_photon_MFP_python(numpar, matter, file_photon_MFP)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   character(*), intent(in) :: file_photon_MFP     ! file with data
   !-------------
   character(300) :: py_photon_MFP
   character(11) :: call_slash
   character(8) :: col
   integer :: N_grid, FN, count_col, i, Nshl, j_start, j, N_tot_col
   real(8) :: t0, t_last, x_tics, x_min, x_max, y_min, y_max
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_labels
   !-------------

   N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
   x_min = 1.0d0  ! [eV] starting energy point
   x_max = matter%Atoms(1)%Ph_MFP(1)%E(N_grid) ! ending energy point
   y_min = 10.0d0
   y_max = 1.0d7
   x_tics = 10.0d0

   ! Photon MFPs plotting:
   py_photon_MFP = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//'OUTPUT_MFP_photon.py'
   open(NEWUNIT=FN, FILE = trim(adjustl(py_photon_MFP)), action="write", status="replace")

   ! Prepare column indices and names:
   call count_shells(matter, N_tot_col)     ! below
   allocate(col_nums(N_tot_col+1), source = 0)
   allocate(col_labels(N_tot_col+1))
   ! Fill the column numbers and labels:
   count_col = 0  ! to start with
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl    ! for all shells of this atom
         if ((i == 1) .and. (j == Nshl)) then ! valence
            ! skip the column, it is not here but at the end
         else
            count_col = count_col + 1  ! number of columns
            col_nums(count_col) = count_col
            ! Names of the columns (skipping the "Time"):
            ! core shell
            col_labels(count_col) = '"'//trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))//'"'
         endif
      enddo
   enddo
   ! Add the valence band MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_labels(count_col) = '"Valence"'
   ! And the total MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_labels(count_col) = '"Total"'

   ! Use them to plot the data:
   call Create_python_plot(FN, file_photon_MFP, col_nums, col_labels, &
      'Photon energy (eV)', 'Attenuation length (A)', 'Photon attenuation length', &
      "best", 'OUTPUT_MFP_photon', trim(adjustl(numpar%fig_extention)), &
      x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, &
      set_x_log=.true., set_y_log=.true., colors_inverted=.true. )     ! below

   close(FN)
   deallocate(col_nums, col_labels)
end subroutine Plot_photon_MFP_python



subroutine Plot_laser_spectrum_python(numpar, laser, file_spectrum, i_pulse)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   character(*), intent(in) :: file_spectrum     ! file with data
   integer, intent(in) :: i_pulse   ! index of the pulse
   !-------------
   character(300) :: py_photon_spectrum
   character(11) :: call_slash
   character(8) :: col
   integer :: N_grid, FN, count_col, i, Nshl, j_start, j, N_tot_col
   real(8) :: t0, t_last, x_tics, x_min, x_max, y_min, y_max
   integer, dimension(:), allocatable :: col_nums, col_nums2
   character(30), dimension(:), allocatable :: col_labels, linestyle
   !-------------

   x_min = laser(i_pulse)%Spectrum(1,1)   ! starting point
   N_grid = size(laser(i_pulse)%Spectrum,2)
   x_max = laser(i_pulse)%Spectrum(1,N_grid)    ! ending point
   y_min = 0.0d0

   ! Photon MFPs plotting:
   py_photon_spectrum = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//'OUTPUT_photon_spectrum.py'
   open(NEWUNIT=FN, FILE = trim(adjustl(py_photon_spectrum)), action="write", status="replace")

   ! Prepare column indices and names:
   allocate(col_nums(3), source = 0)
   allocate(col_labels(3))
   allocate(linestyle(3))
   ! Fill the column numbers and labels:
   col_nums(1) = 1
   col_nums(2) = 2
   col_nums(3) = 3
   ! Names of the columns (skipping the "Energy"):
   col_labels(1) = '"Incoming"'
   col_labels(2) = '"Absorbed"'
   col_labels(3) = '"MC Sampled"'
   ! Style of curves:
   linestyle(1)='"--"'
   linestyle(2)='"-"'
   linestyle(3)='"-."'

   ! Use them to plot the data:
   call Create_python_plot(FN, file_spectrum, col_nums, col_labels, &
      'Photon energy (eV)', 'Photon spectrum (arb. units)', 'Photon spectrum', &
      "best", 'OUTPUT_photon_spectrum', trim(adjustl(numpar%fig_extention)), &
      x_min=x_min, x_max=x_max, y_min=y_min, &
      l_style=linestyle, &
      colors_inverted=.false. )     ! below

   close(FN)
   deallocate(col_nums, col_labels, linestyle)
end subroutine Plot_laser_spectrum_python







!===================================================================
! General routines to use python for plotting:

subroutine Create_Python_animation(FN, Data_file, col_nums, col_labels, &
      x_axis_label, y_axis_label, title, legend_position, out_name, out_format, &
      x_min, x_max, y_min, y_max, t_start, dt, t_end, &
      x_tics, y_tics, &
      set_x_log, set_y_log, &
      l_style, symbols, &
      Data_file2, col_nums2, col_labels2, &
      colors_inverted, first_line_in_front &
      )
   integer, intent(in) :: FN  ! file number (must be opened)
   character(*), intent(in) :: Data_file  ! data file to plot data from
   integer, dimension(:), intent(in) :: col_nums      ! array of columns to plot
   character(*), dimension(:), allocatable, intent(in) :: col_labels    ! array of column labels
   character(*), intent(in) :: x_axis_label, y_axis_label, title
   character(*), intent(in) :: out_name, out_format   ! plot file name and extension
   character(*), intent(in) :: legend_position  ! e.g. 'upper left', 'right bottom' etc.
   real(8), intent(in), optional :: x_min, x_max, y_min, y_max, dt, t_start, t_end
   real(8), intent(in), optional :: x_tics, y_tics
   logical, intent(in), optional :: set_x_log, set_y_log
   character(*), dimension(:), allocatable, intent(in), optional :: l_style    ! array of line-styles
   character(*), dimension(:), allocatable, intent(in), optional :: symbols     ! array of symbols
   character(*), intent(in), optional :: Data_file2   ! data file #2 to plot data from
   integer, dimension(:), intent(in), optional :: col_nums2      ! array of columns to plot #2
   character(*), dimension(:), allocatable, intent(in), optional :: col_labels2    ! array of column labels #2
   logical, intent(in), optional :: colors_inverted, first_line_in_front
   !------------
   integer :: i
   character(10) :: i_ch, i_ch1
   character(40) :: x_min_txt, x_max_txt, y_min_txt, y_max_txt, dt_txt, t_start_txt, t_max_txt, col_txt
   character(10000) :: linestyle_txt, symbols_txt
   logical :: first_in_front


   if (present(first_line_in_front)) then
      first_in_front = first_line_in_front
   else ! default
      first_in_front = .true.
   endif

   ! Set axes:
   if (present(x_min)) then
      write(x_min_txt,'(f24.8)') x_min
   else
      write(x_min_txt,'(a)') 'x.min()'
   endif
   if (present(x_max)) then
      write(x_max_txt,'(f24.8)') x_max
   else
      write(x_max_txt,'(a)') 'x.max()'
   endif
   if (present(y_min)) then
      write(y_min_txt,'(f24.8)') y_min
   else
      write(y_min_txt,'(a)') 'None'
   endif
   if (present(y_max)) then
      write(y_max_txt,'(f24.8)') y_max
   else
      write(y_max_txt,'(a)') 'None'
   endif
   ! Time parameters:
   if (present(dt)) then
      write(dt_txt,'(f24.8)') dt
   else
      write(dt_txt,'(a)') '1.0'
   endif
   if (present(dt)) then
      write(t_start_txt,'(f24.8)') t_start
   else
      write(t_start_txt,'(a)') '0.0'
   endif
   write(t_max_txt,'(f24.1)') max(t_start, abs(t_end))


   ! Line styles, if required:
   if (present(l_style)) then
      linestyle_txt = '['
      do i = 1, size(col_nums)
         ! Line styles:
         if (i > 1) then   ! add come in between
            linestyle_txt = trim(adjustl(linestyle_txt))//','
         endif
         linestyle_txt = trim(adjustl(linestyle_txt))//' '//trim(adjustl(l_style(i)))
      enddo
      linestyle_txt = trim(adjustl(linestyle_txt))//']'
   else ! use default - solid lines:
      linestyle_txt = '['
      do i = 1, size(col_nums)
         ! Line styles:
         if (i > 1) then   ! add come in between
            linestyle_txt = trim(adjustl(linestyle_txt))//','
         endif
         linestyle_txt = trim(adjustl(linestyle_txt))//' "-"'
      enddo
      linestyle_txt = trim(adjustl(linestyle_txt))//']'
   endif

   ! Symbols styles, if required:
   if (present(symbols)) then
      symbols_txt = '['
      do i = 1, size(col_nums)
         ! Line styles:
         if (i > 1) then   ! add come in between
            symbols_txt = trim(adjustl(symbols_txt))//','
         endif
         symbols_txt = trim(adjustl(symbols_txt))//' '//trim(adjustl(symbols(i)))
      enddo
      symbols_txt = trim(adjustl(symbols_txt))//']'
   else ! use default - no symbol, pure lines:
      symbols_txt = '['
      do i = 1, size(col_nums)
         ! Line styles:
         if (i > 1) then   ! add come in between
            symbols_txt = trim(adjustl(symbols_txt))//','
         endif
         symbols_txt = trim(adjustl(symbols_txt))//' "None"'
      enddo
      symbols_txt = trim(adjustl(symbols_txt))//']'
   endif


   !-----------------------------
   write(FN,'(a)') 'import numpy as np'
   write(FN,'(a)') 'import matplotlib.pyplot as plt'

   select case (trim(adjustl(out_format)))
   case ('GIF', 'Gif', 'gif')
      write(FN,'(a)') 'from PIL import Image'
   case default ! assume avi:
      write(FN,'(a)') 'from matplotlib.animation import FFMpegWriter'
   endselect

   !-----------------------------
   write(FN,'(a)') '# 1. Read the multi-block file'
   write(FN,'(a)') 'blocks = []'
   write(FN,'(a)') 'current = []'

   write(FN,'(a)') 'with open("' // trim(adjustl(Data_file)) // '") as f:'
   write(FN,'(a)') '    for line in f:'
   write(FN,'(a)') '          line = line.strip()'
   write(FN,'(a)') '          # Skip comment lines'
   write(FN,'(a)') '          if line.startswith("#"):'
   write(FN,'(a)') '                continue'
   write(FN,'(a)') '          # Block separator'
   write(FN,'(a)') '          if line == "":'
   write(FN,'(a)') '                if current:'
   write(FN,'(a)') '                      blocks.append(np.array(current, float))'
   write(FN,'(a)') '                      current = []'
   write(FN,'(a)') '          else:'
   write(FN,'(a)') '                current.append(line.split())'
   write(FN,'(a)') '# Add last block if not empty'
   write(FN,'(a)') 'if current:'
   write(FN,'(a)') '    blocks.append(np.array(current, float))'

   !-----------------------------
   write(FN,'(a)') '# 2. Prepare animation frames'

   write(FN,'(a)') '# Find the optimal Y-axis:'
   write(FN,'(a)') 'first_block = blocks[0]'
   write(FN,'(a)') 'x = first_block[:, 0]'
   ! Set axes:
   if (present(x_min) .and. present(x_max)) then
      write(FN,'(a)') 'mask = (x >= '//trim(adjustl(x_min_txt))//') & (x <= '//trim(adjustl(x_max_txt))//')'
      write(FN,'(a)') '# Compute y-max only from filtered region'
      write(FN,'(a)') 'ymax = first_block[mask, 1:].max()'
      write(FN,'(a)') 'ymin = first_block[mask, 1:].min()'
   else
      write(FN,'(a)') 'ymax = first_block[0, 1:].max()'
      write(FN,'(a)') 'ymin = first_block[0, 1:].min()'
   endif
   write(FN,'(a)') 'padding = 0.1 * (ymax - ymin)'
   if (present(y_min)) then
      write(FN,'(a)') 'y_lower = '//trim(adjustl(y_min_txt))
   else ! use padding:
      write(FN,'(a)') 'y_lower = min(ymin - padding, 0)'
   endif
   if (present(y_max)) then
      write(FN,'(a)') 'y_upper = '//trim(adjustl(y_max_txt))
   else ! use padding
      write(FN,'(a)') 'y_upper = ymax + padding'
   endif

   !-----------------------------
   write(FN,'(a)') '# 3. Prepare figure ONCE (no jitter)'

   write(FN,'(a)') 'fig, ax = plt.subplots(figsize=(8, 6), dpi=150)'
   write(FN,'(a)') '# Create placeholder lines ONCE'
   write(FN,'(a)') 'linestyles = '//trim(adjustl(linestyle_txt))
   write(FN,'(a)') 'markers = '//trim(adjustl(symbols_txt))

   do i = 1, size(col_nums)   ! for all lines to plot
      write(i_ch, '(i0)') i
      write(i_ch1, '(i0)') i-1
      write(col_txt, '(i0)') col_nums(i)
      if (i == 1) then ! first line
         if (first_in_front) then ! enforce the first line to be in front of others:
            write(FN,'(a)') 'line'//trim(adjustl(i_ch))// &
               ', = ax.plot([], [], lw=2, color="black", linestyle=linestyles[0], marker=markers[0], label="'// &
               trim(adjustl(t_max_txt))//' fs '//trim(adjustl(col_labels(i)))//'", zorder=10)'
         else ! first line last:
            write(FN,'(a)') 'line'//trim(adjustl(i_ch))// &
               ', = ax.plot([], [], lw=2, color="black", linestyle=linestyles[0], marker=markers[0], label="'// &
               trim(adjustl(t_max_txt))//' fs '//trim(adjustl(col_labels(i)))//'")'
         endif
      else ! the rest
         write(FN,'(a)') 'line'//trim(adjustl(i_ch))// &
         ', = ax.plot([], [], lw=1.5, linestyle=linestyles['//trim(adjustl(i_ch1))// &
         '], marker=markers['//trim(adjustl(i_ch1))//'],  label="    '// &
         trim(adjustl(col_labels(i)))//'")'
      endif
   enddo

   write(FN,'(a)') '# Create legend once to prevent jitter:'
   write(FN,'(a)') 'legend = ax.legend('
   write(FN,'(a)') '    loc="upper right",'
   write(FN,'(a)') '    bbox_to_anchor=(1.0, 1.0),'
   write(FN,'(a)') '    fontsize=14,'
   write(FN,'(a)') '    frameon=True,'
   write(FN,'(a)') "    prop={'family': 'monospace', 'size': 14}"
   write(FN,'(a)') ')'


   write(FN,'(a)') '# Fix axes once:'
   write(FN,'(a)') 'ax.set_xlim('//trim(adjustl(x_min_txt))//', '//trim(adjustl(x_max_txt))//')'

   write(FN,'(a)') 'ax.set_ylim(y_lower, y_upper)'

   if (present(set_y_log)) then
      if (set_y_log) write(FN,'(a)') 'ax.set_yscale("log")'
   endif

   write(FN,'(a)') 'ax.set_xlabel("'//trim(adjustl(x_axis_label))//'", fontsize=14)'
   write(FN,'(a)') 'ax.set_ylabel("'//trim(adjustl(y_axis_label))//'", fontsize=14)'
   ! Set tics:
   write(FN,'(a)') 'plt.xticks(fontsize=12)'
   write(FN,'(a)') 'plt.yticks(fontsize=12)'
   write(FN,'(a)') 'ax.set_title("'//trim(adjustl(title))//'", fontsize=14, pad=12)'
   write(FN,'(a)') 'fig.tight_layout()'


   !-----------------------------
   ! Choose which format of animation to use:
   select case (trim(adjustl(out_format)))      ! animated gif:
   case ('gif')
      write(FN,'(a)') '# 4. Animation loop (update only data + legend text):'

      write(FN,'(a)') 'frames = []'
      write(FN,'(a)') 'for i, block in enumerate(blocks):'
      write(FN,'(a)') '    x = block[:, 0]'

      ! For all columns that we want to plot:
      do i = 1, size(col_nums)
         write(i_ch, '(i0)') i
         write(col_txt, '(i0)') col_nums(i)
         write(FN,'(a)') '    y'//trim(adjustl(i_ch))//' = block[:, '// trim(adjustl(col_txt)) //']'
      enddo

      write(FN,'(a)') '   # Update line data:'
      do i = 1, size(col_nums)
          write(i_ch, '(i0)') i
          write(FN,'(a)') '    line'//trim(adjustl(i_ch))//'.set_data(x, y'//trim(adjustl(i_ch))//')'
      enddo

      write(FN,'(a)') '    # Update legend text (first entry only):'
      write(FN,'(a)') '    time_fs = i*'//trim(adjustl(dt_txt))//' + ('//trim(adjustl(t_start_txt))//')'
      write(FN,'(a)') '    legend.get_texts()[0].set_text(f"{time_fs:6.1f} fs '//trim(adjustl(col_labels(1)))//'")'

      write(FN,'(a)') '    # Convert figure to image'
      write(FN,'(a)') '    fig.canvas.draw()'
      ! New Matplotlib 3.8+ method:
      write(FN,'(a)') '    buf = fig.canvas.buffer_rgba()'
      write(FN,'(a)') '    w, h = fig.canvas.get_width_height()'
      write(FN,'(a)') '    frame = Image.frombuffer("RGBA", (w, h), buf, "raw", "RGBA", 0, 1).copy()'
      write(FN,'(a)') '    frames.append(frame)'

      write(FN,'(a)') '# Save as animated GIF'
      write(FN,'(a)') 'frames[0].save('
      write(FN,'(a)') '    "'//trim(adjustl(out_name))//'.gif",'
      write(FN,'(a)') '    save_all=True,'
      write(FN,'(a)') '    append_images=frames[1:],'
      write(FN,'(a)') '    duration=100,'   ! 100 ms = gnuplot delay 10
      write(FN,'(a)') '    loop=0'
      write(FN,'(a)') ')'
      return      ! this part is done, for gif there is nothing else to do
   case ('mp4')
      write(FN,'(a)') 'writer = FFMpegWriter('
      write(FN,'(a)') '       fps=10,'
      !write(FN,'(a)') '       codec="libx264",'      ! this codec is often missing on HPC
      write(FN,'(a)') '       codec="mpeg4",'
      write(FN,'(a)') '       bitrate=5000,'
      write(FN,'(a)') '       extra_args=["-pix_fmt", "yuv420p"]'
      write(FN,'(a)') ' )'
      !write(FN,'(a)') 'with writer.saving(fig, "'//trim(adjustl(out_name))//'.mp4", dpi=150):'
   case ('webm')
      write(FN,'(a)') 'writer = FFMpegWriter('
      write(FN,'(a)') 'fps=10,'
      write(FN,'(a)') 'codec="libvpx-vp9",'
      write(FN,'(a)') 'extra_args=["-crf", "30", "-b:v", "0"]'
      write(FN,'(a)') ')'
   case default ! all other format of video: avi, mkv, mov
      write(FN,'(a)') '# Save as avi'
      write(FN,'(a)') 'writer = FFMpegWriter(fps=10, bitrate=1800)'
      !write(FN,'(a)') 'with writer.saving(fig, "'//trim(adjustl(out_name))//'.avi", dpi=150):'
   end select

   write(FN,'(a)') 'with writer.saving(fig, "'//trim(adjustl(out_name))//'.'//trim(adjustl(out_format))//'", dpi=150):'

   write(FN,'(a)') '    for i, block in enumerate(blocks):'
   write(FN,'(a)') '        x = block[:, 0]'
   do i = 1, size(col_nums)      ! for all curves
      write(i_ch, '(i0)') i
      write(col_txt, '(i0)') col_nums(i)
      write(FN,'(a)') '        y'//trim(adjustl(i_ch))//' = block[:, '// trim(adjustl(col_txt)) //']'
   enddo

   write(FN,'(a)') '        # Update line data:'
   do i = 1, size(col_nums)
      write(i_ch, '(i0)') i
      write(FN,'(a)') '        line'//trim(adjustl(i_ch))//'.set_data(x, y'//trim(adjustl(i_ch))//')'
   enddo

   write(FN,'(a)') '        # Update legend text (first entry only):'
   write(FN,'(a)') '        time_fs = i*'//trim(adjustl(dt_txt))//' + ('//trim(adjustl(t_start_txt))//')'
   write(FN,'(a)') '        legend.get_texts()[0].set_text(f"{time_fs:6.1f} fs '//trim(adjustl(col_labels(1)))//'")'

   write(FN,'(a)') '        # Draw and save frame'
   write(FN,'(a)') '        fig.canvas.draw()'
   write(FN,'(a)') '        writer.grab_frame()'
   write(FN,'(a)') '    # Force final frame flush'
   write(FN,'(a)') '    fig.canvas.draw()'
   write(FN,'(a)') '    writer.grab_frame()'
end subroutine Create_Python_animation



subroutine Create_python_plot(FN, Data_file, col_nums, col_labels, &
      x_axis_label, y_axis_label, title, legend_position, out_name, out_format, &
      x_min, x_max, y_min, y_max, &
      x_tics, y_tics, &
      set_x_log, set_y_log, &
      l_style, &
      Data_file2, col_nums2, col_labels2, &
      colors_inverted &
      )
   integer, intent(in) :: FN  ! file number (must be opened)
   character(*), intent(in) :: Data_file  ! data file to plot data from
   integer, dimension(:), intent(in) :: col_nums      ! array of columns to plot
   character(*), dimension(:), allocatable, intent(in) :: col_labels    ! array of column labels
   character(*), intent(in) :: x_axis_label, y_axis_label, title
   character(*), intent(in) :: out_name, out_format   ! plot file name and extension
   character(*), intent(in) :: legend_position  ! e.g. 'upper left', 'right bottom' etc.
   real(8), intent(in), optional :: x_min, x_max, y_min, y_max
   real(8), intent(in), optional :: x_tics, y_tics
   logical, intent(in), optional :: set_x_log, set_y_log
   character(*), dimension(:), allocatable, intent(in), optional :: l_style    ! array of line-styles
   character(*), intent(in), optional :: Data_file2   ! data file #2 to plot data from
   integer, dimension(:), intent(in), optional :: col_nums2      ! array of columns to plot #2
   character(*), dimension(:), allocatable, intent(in), optional :: col_labels2    ! array of column labels #2
   logical, intent(in), optional :: colors_inverted
   !------------
   character(100) :: x_min_txt, x_max_txt, y_min_txt, y_max_txt, x_tics_text, y_tics_text, temp_txt
   character(25) :: ch_var
   character(3) :: font_size
   logical :: tics_present, color_inverse, ls_present
   integer :: i, N_cols
   !------------

   if (present(x_tics) .or. present(x_tics)) then
      tics_present = .true.
   else
      tics_present = .false.
   endif

   if (present(colors_inverted)) then
      color_inverse = colors_inverted
   else ! use default
      color_inverse = .false. ! by default, don't invert
   endif

   ! Total number of curves to be plotted:
   N_cols = size(col_nums)
   if (present(col_nums2)) N_cols = N_cols + size(col_nums2)   ! include the second one, if present


   write(FN,'(a)') 'import pandas as pd'
   write(FN,'(a)') 'import numpy as np'
   write(FN,'(a)') 'import matplotlib.pyplot as plt'
   write(FN,'(a)') 'import matplotlib.cm as cm'
   if (tics_present) then
      write(FN,'(a)') 'from matplotlib.ticker import MultipleLocator'
   endif

   write(FN,'(a)') '# Set the format:'
   write(FN,'(a)') 'plt.xlabel("'// trim(adjustl(x_axis_label)) //'", fontsize=14)'
   write(FN,'(a)') 'plt.ylabel("'// trim(adjustl(y_axis_label)) //'", fontsize=14)'

   ! Set tics:
   write(FN,'(a)') 'plt.xticks(fontsize=12)'
   write(FN,'(a)') 'plt.yticks(fontsize=12)'
   if (present(x_tics)) then
      write(x_tics_text,'(f24.8)') x_tics
      write(FN,'(a)') 'plt.gca().xaxis.set_major_locator(MultipleLocator('//trim(adjustl(x_tics_text))//'))'
   endif
   if (present(y_tics)) then
      write(y_tics_text,'(f24.8)') y_tics
      write(FN,'(a)') 'plt.gca().yaxis.set_major_locator(MultipleLocator('//trim(adjustl(y_tics_text))//'))'
   endif

   write(FN,'(a)') 'plt.grid(False)'
   write(FN,'(a)') '# Set the axes:'


   ! Process the data file:
   write(FN,'(a)') '# Read the output file:'
   write(FN,'(a)') 'df = pd.read_csv(r"'//trim(adjustl(Data_file))//'", '//"sep=r'\s+', "//'header=None, comment="#", skipinitialspace=True)'

   ! Log-scale, if required:
   if (present(set_x_log)) then
      if (set_x_log) write(FN,'(a)') 'plt.xscale("log")'
   endif
   if (present(set_y_log)) then
      if (set_y_log) write(FN,'(a)')'plt.yscale("log")'
   endif

   ! Prepare the plot:
   write(FN,'(a)') '# Prepare the plot:'
   ! Set a list of which columns to plot:
   write(FN,'(a)', advance='no') 'columns_to_plot = ['
   do i = 1, size(col_nums)
      ! Column numbers:
      if (i > 1) then   ! add come in between
         write(FN,'(a)', advance='no') ', '
      endif
      write(temp_txt,'(i)') col_nums(i)
      write(FN,'(a)', advance='no') trim(adjustl(temp_txt))
   enddo
   write(FN,'(a)') ']'

   ! Column labels, if required:
   if (allocated(col_labels)) then ! the legend is required
      write(FN,'(a)', advance='no') 'labels = ['
      do i = 1, size(col_labels)
         ! Column titles:
         if (i > 1) then   ! add come in between
            write(FN,'(a)', advance='no') ', '
         endif
         write(FN,'(a)', advance='no') trim(adjustl(col_labels(i)))
      enddo
      write(FN,'(a)') ']'
   endif

   ! Colormap:
   write(FN,'(a)') '# Create N distinct colors from a colormap'
   write(FN,'(a)') 'N = len(columns_to_plot)'

   if (N_cols <=10) then ! use standard pallette
      if (color_inverse) then
         write(FN,'(a)') 'colors = plt.cm.tab10(np.linspace(1, 0, 10))'
      else ! default color range
         write(FN,'(a)') 'colors = plt.cm.tab10(np.linspace(0, 1, 10))'
      endif
   elseif (N_cols <= 20) then ! use extended one
      if (color_inverse) then
         write(FN,'(a)') 'colors = plt.cm.tab20(np.linspace(1, 0, 20))'
      else ! default color range
         write(FN,'(a)') 'colors = plt.cm.tab20(np.linspace(0, 1, 20))'
      endif
   else     ! use combined one
      if (color_inverse) then
         write(FN,'(a)') 'colors = np.vstack([plt.cm.tab20(np.linspace(1, 0, 20)), plt.cm.tab20b(np.linspace(1, 0, 20))])'
      else ! default color range
         write(FN,'(a)') 'colors = np.vstack([plt.cm.tab20(np.linspace(0, 1, 20)), plt.cm.tab20b(np.linspace(0, 1, 20))])'
      endif
   endif

   ! Line-styles:
   write(FN,'(a)') '# Set line styles:'
   write(FN,'(a)', advance='no') 'linestyles = '
      ! Line styles, if required:
   if (present(l_style)) then
      write(FN,'(a)', advance='no') '['
      do i = 1, size(col_nums)
         ! Line styles:
         if (i > 1) then   ! add come in between
            write(FN,'(a)', advance='no') ', '
         endif
         write(FN,'(a)', advance='no') trim(adjustl(l_style(i)))
      enddo
      write(FN,'(a)') ']'
   else ! use default - solid lines:
      write(FN,'(a)', advance='no') '['
      do i = 1, size(col_nums)
         ! Line styles:
         if (i > 1) then   ! add come in between
            write(FN,'(a)', advance='no') ', '
         endif
         write(FN,'(a)', advance='no') '"-"'
      enddo
      write(FN,'(a)') ']'
   endif

   ! Make the plot:
   write(FN,'(a)') 'for i, col in enumerate(columns_to_plot):'          ! cycle for all curves
   write(FN,'(a)') '    color = colors[i % len(colors)]'                ! repeat colors
   write(FN,'(a)') '    ls    = linestyles[i % len(linestyles)]'        ! repeat line styles
   if (allocated(col_labels)) then
      write(FN,'(a)') '    label = labels[i % len(labels)]'                ! labels
   endif
   write(FN,'(a)') '    plt.plot(df.iloc[:, 0], df.iloc[:, col],'       ! plot columns
   write(FN,'(a)') '    color=color,'                                   ! set color
   if (allocated(col_labels)) then                                      ! the legend is required
      write(FN,'(a)') '    label=label,'
   endif
   write(FN,'(a)') '    linestyle=ls)'                                  ! line style


   !-----------------------
   ! If we want to add data from another file on the same plot:
   if (present(Data_file2) .and. present(col_nums2) .and. present(col_labels2)) then !# Add a curve from the second file
      write(FN,'(a)') '# Add a curve from the second file'
      write(FN,'(a)') 'df2 = pd.read_csv(r"'//trim(adjustl(Data_file2))//'", '//"sep=r'\s+', "// &
                        'header=None, comment="#", skipinitialspace=True)'

      ! Set a list of which columns to plot:
      write(FN,'(a)', advance='no') 'columns_to_plot2 = ['
      do i = 1, size(col_nums2)
         ! Column numbers:
         if (i > 1) then   ! add come in between
            write(FN,'(a)', advance='no') ', '
         endif
         write(temp_txt,'(i)') col_nums2(i)
         write(FN,'(a)', advance='no') trim(adjustl(temp_txt))
      enddo
      write(FN,'(a)')  ']'

      ! Column labels, if required:
      if (allocated(col_labels2)) then ! the legend is required:
         write(FN,'(a)', advance='no') 'labels2 = ['
         do i = 1, size(col_labels2)
            ! Column numbers:
            if (i > 1) then   ! add come in between
               write(FN,'(a)', advance='no') ', '
            endif
            write(FN,'(a)', advance='no') trim(adjustl(col_labels2(i)))
         enddo
         write(FN,'(a)') ']'
      endif

      write(FN,'(a)') '# Create second part of the plot:'
      if (allocated(col_labels2)) then ! the legend is required:
         write(FN,'(a)') 'for col, label in zip(columns_to_plot2, labels2):'
         write(FN,'(a)') '    plt.plot(df2.iloc[:, 0], df2.iloc[:, col], label=label, linestyle="--")'
      else ! No legend:
         write(FN,'(a)') 'for col in columns_to_plot2:'
         write(FN,'(a)') '    plt.plot(df2.iloc[:, 0], df2.iloc[:, col])'
      endif
   endif
   !-----------------------


   ! Now add the details to the plot:
   ! Set axes:
   if (present(x_min)) then
      write(x_min_txt,'(f24.8)') x_min
   else
      write(x_min_txt,'(a)') 'None'
   endif
   if (present(x_max)) then
      write(x_max_txt,'(f24.8)') x_max
      ! Make sure the maximum value is adjusted but not larger than the given one:
      write(FN,'(a)') 'xmax = df.iloc[:, 0].max()'
      x_max_txt = 'min(xmax, '//trim(adjustl(x_max_txt))//')'
   else
      write(x_max_txt,'(a)') 'None'
   endif
   if (present(y_min)) then
      write(y_min_txt,'(f24.8)') y_min
   else
      write(y_min_txt,'(a)') 'None'
   endif
   if (present(y_max)) then
      write(y_max_txt,'(f24.8)') y_max
   else
      write(y_max_txt,'(a)') 'None'
   endif
   write(FN,'(a)') 'plt.xlim('//trim(adjustl(x_min_txt))//','//trim(adjustl(x_max_txt))//')'
   write(FN,'(a)') 'plt.ylim('//trim(adjustl(y_min_txt))//','//trim(adjustl(y_max_txt))//')'

   write(FN,'(a)') 'plt.title("'//trim(adjustl(title))//'")'

   ! Set the legend:
   if (allocated(col_labels)) then ! the legend is required
      ! Make appropriate font size:
      if (N_cols < 8) then ! normal:
         write(temp_txt,'(a)') 'fontsize=12'
      elseif (N_cols < 17) then ! smaller
         write(temp_txt,'(a)') 'fontsize=10'
      elseif (N_cols < 21) then
         write(temp_txt,'(a)') 'fontsize=9'
      elseif (N_cols < 33) then
         write(temp_txt,'(a)') 'fontsize=10, ncol=2'
      elseif (N_cols < 41) then
         write(temp_txt,'(a)') 'fontsize=9, ncol=2'
      else
         write(temp_txt,'(a)') 'fontsize=9, ncol=3'
      endif
      write(FN,'(a)') 'plt.legend(loc="'//trim(adjustl(legend_position))//'", '//trim(adjustl(temp_txt))//')'
   endif


   write(FN,'(a)') '# Save the plot in this format:'
   write(FN,'(a)') 'plt.savefig("'//trim(adjustl(out_name))//'.'//trim(adjustl(out_format))//'", dpi=300, bbox_inches="tight")'
end subroutine Create_python_plot




subroutine define_python_call(path_sep, py_call)
   character(*), intent(in) :: path_sep
   character(*), intent(out) :: py_call

   if (path_sep .EQ. '\') then	! if it is Windows
      py_call = 'python '
   else ! It is linux
      py_call = 'python3 '
   endif
end subroutine define_python_call


subroutine count_shells(matter, Tot_shl)
   type(Solid), intent(in) :: matter ! parameters of the material
   integer, intent(out) :: Tot_shl
   !----------------------
   integer :: count_col, i, j, Nshl

   count_col = 0 ! to start with
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl    ! for all shells of this atom
         count_col = count_col + 1  ! number of columns
      enddo
   enddo
   Tot_shl = count_col
end subroutine count_shells



subroutine select_linestyle(i, line_style)
   integer, intent(in) :: i   ! index
   character(*), intent(out) :: line_style      ! linestyle for python plotting
   !------------------------
   integer :: i_allowed

   ! Currently, there are 7 stiles, so cycle them:
   i_allowed = MOD(i, 7)

   select case (i_allowed)   ! currently, supports 7 different elements (may be added if needed)
   case (1)    ! solid
      line_style = '"-"'
   case (2)    ! dashed
      line_style = '"--"'
   case (3)    ! dash-dot
      line_style = '"-."'
   case (4)    ! dash-dot-dot
      line_style = '(0,(5,2,1,2,1,2))'
   case (5)    ! long dash-short dash-dot
      line_style = '(0,(10,3,4,3,1,2))'
   case (6)    ! dot
      line_style = '":"'
   case default ! Long dash - dot - short dash - dot
      line_style = '(0,(10,3,1,2,4,3,1,2))'
   end select
end subroutine select_linestyle




subroutine select_symbols(i, symbol)
   integer, intent(in) :: i   ! index
   character(*), intent(out) :: symbol      ! symbol for python plotting
   !------------------------
   integer :: i_allowed
   ! Currently, there are 15 stiles, so cycle them:
   i_allowed = MOD(i, 15)

   select case (i_allowed)
   case (1)
      symbol = '"o"' ! circle
   case (2)
      symbol = '"s"'	! square
   case (3)
      symbol = '"D"'	! diamond
   case (4)
      symbol = '"p"'	! pentagon
   case (5)
      symbol = '"h"'	! hexagon 1
   case (6)
      symbol = '"H"'	! hexagon 2
   case (7)
      symbol = '"X"'	! x-filled (thick)
   case (8)
      symbol = '"8"'	! octagon
   case (9)
      symbol = '"^"'	! triangle up
   case (10)
      symbol = '"v"'	! triangle down
   case (11)
      symbol = '">"'	! triangle right
   case (12)
      symbol = '"<"'	! triangle left
   case (13)
      symbol = '"P"'	! plus‑filled (thick)
   case (14)
      symbol = '"d"'	! thin diamond
   case default
      symbol = '"*"'	! star
   end select
end subroutine select_symbols



!===================================================
! Collecting all python plotting all the scripts:

subroutine execute_all_pyplots(numpar, file_path)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_path
   !----------------
   character(300) :: command
   integer :: iret

   !call chdir(trim(adjustl(file_path)))
   command = trim(adjustl(file_path))
   iret = chdir(command)

   if (trim(adjustl(numpar%path_sep)) == '\') then	! if it is Windows
      command = 'python '//trim(adjustl(m_Python_plot_all))
      iret = system(command)
   else ! linux:
      command = 'python3 '//trim(adjustl(m_Python_plot_all))
      iret = system(command)
   endif

   iret = chdir("../")
end subroutine execute_all_pyplots


subroutine collect_python_plots(path_sep, out_path, skip_execution)
   character(*), intent(in) :: path_sep, out_path    ! folder with the cmd-files
   logical, intent(in), optional :: skip_execution  ! if you don't want to execute all gnuplots
   !------------------------
   character(300) :: File_name, command, Py_plot_all_files, File_name_withquotes
   integer :: FN, N_f, i, n_slash
   integer :: open_status, iret, idir, leng
   character(300), dimension(:), allocatable :: All_files
   character(300) :: output_path
   character(5) ::  call_slash, sh_cmd
   logical :: skip_exec

   if (present(skip_execution)) then   ! if requested, you may skip execution of gnuplot scripts
      skip_exec = skip_execution
   else  ! be default, execute gnuplot scripts
      skip_exec = .false.
   endif

   ! In which directory to collect all gnu scripts:
   output_path = out_path

   ! Create a temporary file:
   Py_plot_all_files = trim(adjustl(m_Python_plot_all))

   ! Include the path to the directory:
   File_name_withquotes = '"'//trim(adjustl(output_path))//'"'//trim(adjustl(path_sep))//trim(adjustl(Py_plot_all_files))
   File_name = trim(adjustl(output_path))//trim(adjustl(path_sep))//trim(adjustl(Py_plot_all_files))

   ! Save the names of all gnuplot scripts into this file:
   if (trim(adjustl(path_sep)) == '\') then  ! if it is Windows
      command = 'dir /o-d "'//trim(adjustl(output_path))//'"\*.py /b >'//trim(adjustl(File_name_withquotes))
   else ! linux:
      command = "ls -t "//trim(adjustl(output_path))//" | grep '.py' >"//trim(adjustl(File_name))
   endif

   iret = system(trim(adjustl(command)))   ! execute the command to save file names in the temp file

   ! Open the files with gnuplot script names:
   open(NEWUNIT=FN, file=trim(adjustl(File_name)), iostat=open_status)
   if ( open_status /= 0 ) then
      print *, 'Could not open ',trim(adjustl(File_name)),' for py-plotting. Unit = ', FN
   endif

   ! Find out how many there are:
   call Count_lines_in_file(FN, N_f) ! module "Dealing_with_files"

   ! Allocate array with them:
   allocate(All_files(N_f)) ! array with all relevant file names
   !print*, 'N_f:', N_f

   ! Read file names:
   do i = 1,N_f
      read(FN,*) All_files(i)
      !print*, i, trim(adjustl(All_files(i)))
   enddo

   ! Rewind file to overwrite including the calls:
   rewind(FN)

   ! Make the master py-script:
   write(FN,'(a)') 'import subprocess'
   write(FN,'(a)') 'import sys'

   ! Collect all the py-names into this file:
   do i = 1,N_f
      if (trim(adjustl(All_files(i))) /= trim(adjustl(Py_plot_all_files))) then ! exclude the file itself
         write(FN,'(a)') 'subprocess.run([sys.executable, "'//trim(adjustl(All_files(i)))//'"])'
      endif
   enddo
   close (FN)

   !--------------
   if (skip_exec) return   ! If execution of py-plot scripts is not requested, we are done;
   ! otherwise, execute all the gnuplot scripts, if requested:
   idir = chdir(trim(adjustl(output_path))) ! go into the directory with output files

   if (trim(adjustl(path_sep)) == '\') then	! if it is Windows
      iret = system('python '//trim(adjustl(Py_plot_all_files)))
   else ! linux:
      iret = system('python3 '//trim(adjustl(Py_plot_all_files)))
   endif

   ! Count how many times the system has to go out of the directory to get back into the original directory:
   ! Defined by the number of slashes in the path given:
   n_slash = count( (/ (trim(adjustl(output_path(i:i))), i=1,len_trim(output_path)) /) == trim(adjustl(path_sep)) )
   do i = 1, n_slash+1  ! go up in the directory tree as many times as needed
      idir = chdir("../")    ! exit the directory with output files
   enddo
end subroutine collect_python_plots


end module Plots_python
