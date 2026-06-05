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
module Plots_gnuplot

#ifndef __GFORTRAN__
USE IFLPORT, only : system, chdir   ! library, allowing to operate with directories in intel fortran
#endif

use Universal_constants
use Objects
use Dealing_with_files, only : Count_lines_in_file
use Gnuplotting, only : write_gnuplot_script_header_new, write_gnuplot_script_ending_new, select_linestyle_gnu
use Little_subroutines, only : set_starting_time, order_of_time, name_of_orbitals, number_of_types_of_orbitals
use Dealing_with_EADL, only : define_PQN


implicit none 
PRIVATE


character(50), parameter :: m_Gnuplot_all = 'OUTPUT_Gnuplot_all'


public :: m_Gnuplot_all, Plot_electron_MFP_gunplot, Plot_photon_MFP_gunplot, create_gnuplot_scripts, &
            execute_all_gnuplots, Plot_laser_spectrum_gnuplot


 contains




subroutine create_gnuplot_scripts(Scell,matter,numpar,laser, file_path, file_temperatures, file_pressure, file_energies, &
file_supercell, file_electron_properties, file_heat_capacity, file_heat_capacity_dyn, &
file_numbers, file_orb, file_deep_holes, file_optics, file_Ei, file_PCF, file_NN, file_element_NN, file_electron_entropy, file_Te, file_mu, &
file_atomic_entropy, file_atomic_temperatures, file_atomic_temperatures_part, file_sect_displ, &
file_diffraction_peaks, file_diffraction_peaks_part, file_diffraction_powder, file_diffraction_peaks_DW, file_Debye_temperature, &
file_testmode, file_coupling, file_high_e)
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

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      call_slash = 'call '
      sh_cmd = '.cmd'
   else ! It is linux
      call_slash = './'
      sh_cmd = '.sh'
   endif


   ! Energy levels:
   if (numpar%save_Ei) then
      ! Find order of the number, and set number of tics as tenth of it:
      call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

      File_name  = trim(adjustl(file_path))//'OUTPUT_energy_levels_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 0.2d0, x_tics, 'Energy levels', 'Time (fs)', 'Energy levels (eV)', 'OUTPUT_energy_levels.'//trim(adjustl(numpar%fig_extention)), numpar%path_sep, setkey=4)

      call write_energy_levels_gnuplot(numpar, FN, t0, Scell, 'OUTPUT_energy_levels.dat')
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)
   endif


   ! Energies:
   File_name  = trim(adjustl(file_path))//'OUTPUT_energies_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_energies(numpar, File_name, file_energies, t0, t_last, 'OUTPUT_energies.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Temepratures:
   File_name  = trim(adjustl(file_path))//'OUTPUT_temperatures_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_temperatures(numpar, matter, File_name, file_temperatures, t0, t_last, 'OUTPUT_temepratures.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Mean square displacement:
   File_name  = trim(adjustl(file_path))//'OUTPUT_mean_displacement_Gnuplot'//trim(adjustl(sh_cmd))
   if (abs(numpar%MSD_power) > 1.0e-6) then ! only plot it if it's not zero
      call gnu_MSD(matter, numpar, File_name, file_temperatures, t0, t_last, &
            'OUTPUT_mean_displacement.'//trim(adjustl(numpar%fig_extention)), &
            numpar%MSD_power) ! below
   endif

   ! Atomic masks for sectional displacements:
   if (allocated(Scell(1)%Displ)) then
      Nsiz = size(Scell(1)%Displ)   ! how many masks
      do j = 1, Nsiz    ! for all masks
         File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(j)%mask_name))// &
                  '_Gnuplot'//trim(adjustl(sh_cmd))
         call gnu_displacements(numpar, File_name, file_sect_displ(j), t0, t_last, 'OUTPUT_mean_displacement_'// &
               trim(adjustl(Scell(1)%Displ(j)%mask_name))//'.'//trim(adjustl(numpar%fig_extention)), &
               Scell(1)%Displ(j)%MSD_power) ! below
         ! Partial by elements, if there is more than one:
         File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(j)%mask_name))// &
                  '_partial_Gnuplot'//trim(adjustl(sh_cmd))
         if (matter%N_KAO > 1) then
            call gnu_displacements_partial(File_name, file_sect_displ(j), t0, t_last, 'OUTPUT_mean_displacement_'// &
               trim(adjustl(Scell(1)%Displ(j)%mask_name))//'_partial.'//trim(adjustl(numpar%fig_extention)), &
               Scell(1)%Displ(j)%MSD_power, matter, numpar) ! below
         endif
      enddo ! j
   endif


   ! Diffraction:
   if (numpar%save_diff_peaks) then
      ! Diffraction peaks:
      File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_diffraction_peaks(Scell(1), numpar, File_name, file_diffraction_peaks, t0, t_last, &
                                    'OUTPUT_diffraction_peaks.'//trim(adjustl(numpar%fig_extention)), .false.) ! below

      ! For element-specific diffraction data:
      if (size(matter%Atoms) > 1) then
          do j = 1, size(matter%Atoms)    ! for all elements
             File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_'//trim(adjustl(matter%Atoms(j)%Name))// &
                          '_Gnuplot'//trim(adjustl(sh_cmd))
             call gnu_diffraction_peaks(Scell(1), numpar, File_name, file_diffraction_peaks_part(j), t0, t_last, &
                                    'OUTPUT_diffraction_peaks_'//trim(adjustl(matter%Atoms(j)%Name))// &
                                    '.'//trim(adjustl(numpar%fig_extention)), .false.) ! below
          enddo
      endif


      ! Powder diffraction:
      File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_powder_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 10.0d0, 'Powder', '2theta (deg)', 'Peak Intensity (a.u.)', 'OUTPUT_diffraction_powder.gif', numpar%path_sep, setkey=0)

      call write_diffraction_powder_gnuplot(FN, Scell, matter, numpar, trim(adjustl(file_diffraction_powder)))      ! below

      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)

      ! Check if Debye-Waller analysis is required:
      if ( abs(numpar%DW_theta) > 1.0d-6 ) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_DW_Gnuplot'//trim(adjustl(sh_cmd))
         call gnu_diffraction_peaks(Scell(1), numpar, File_name, file_diffraction_peaks_DW, t0, t_last, &
                                    'OUTPUT_diffraction_peaks_DW.'//trim(adjustl(numpar%fig_extention)), .true.) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_Debye_temperatures_Gnuplot'//trim(adjustl(sh_cmd))
         call gnu_Debye_temperatures(Scell(1), numpar, File_name, file_Debye_temperature, t0, t_last, &
                                    'OUTPUT_Debye_temperatures.'//trim(adjustl(numpar%fig_extention))) ! below
      endif
   endif

   ! Pressure:
   File_name  = trim(adjustl(file_path))//'OUTPUT_pressure_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_pressure(numpar, File_name, file_pressure, t0, t_last, 'OUTPUT_pressure.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Stress tensor:
   File_name  = trim(adjustl(file_path))//'OUTPUT_stress_tensor_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_stress(numpar, File_name, file_pressure, t0, t_last, 'OUTPUT_pressure_tensor.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Numbers of particles:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_and_holes_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_numbers(numpar, File_name, file_numbers, t0, t_last, 'OUTPUT_electrons_holes.'//trim(adjustl(numpar%fig_extention))) ! below

   ! High-energy electrons:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_high_energy'//trim(adjustl(sh_cmd))
   call gnu_high_energy_el(numpar, File_name, file_high_e, t0, t_last, 'OUTPUT_electrons_high_energy.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Orbital-resolved electron parameters:
   File_name  = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_orbital_resolved(Scell(1), matter, numpar, File_name, file_orb, t0, t_last, 'OUTPUT_orbital_resolved_Ne.'// &
               trim(adjustl(numpar%fig_extention))) ! below

   ! Numbers of CB electrons:
   File_name  = trim(adjustl(file_path))//'OUTPUT_CB_electrons_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_CB_electrons(numpar, File_name, file_numbers, t0, t_last, 'OUTPUT_CB_electrons.'//trim(adjustl(numpar%fig_extention)))

   ! Numbers of deep-shell holes:
   File_name  = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_holes(numpar, File_name, file_deep_holes, t0, t_last, matter, 'OUTPUT_deep_shell_holes.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Band gap:
   File_name  = trim(adjustl(file_path))//'OUTPUT_Egap_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_Egap(numpar, File_name, file_electron_properties, t0, t_last, 'OUTPUT_Egap.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Chemical potential and Ne:
   File_name  = trim(adjustl(file_path))//'OUTPUT_mu_and_Ne_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_mu(numpar, File_name, file_electron_properties, t0, t_last, 'OUTPUT_mu_and_Ne.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Boundaries of the bands:
   File_name  = trim(adjustl(file_path))//'OUTPUT_bands_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_Ebands(numpar, File_name, file_electron_properties, t0, t_last, 'OUTPUT_bands.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron heat capacity:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electron_Ce'//trim(adjustl(sh_cmd))
   call gnu_capacity(numpar, File_name, file_electron_properties, t0, t_last, 'OUTPUT_electron_Ce.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron heat conductivity:
   if (numpar%do_kappa) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity'//trim(adjustl(sh_cmd))
      call gnu_heat_conductivity(numpar, File_name, file_heat_capacity, &
            numpar%kappa_Te_min, numpar%kappa_Te_max, &
            'OUTPUT_electron_heat_conductivity.'//trim(adjustl(numpar%fig_extention))) ! below
   endif
   if (numpar%do_kappa_dyn) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_dyn'//trim(adjustl(sh_cmd))
      call gnu_heat_conductivity_dyn(numpar, File_name, file_heat_capacity_dyn, &
            numpar%kappa_Te_min, numpar%kappa_Te_max, &
            'OUTPUT_electron_heat_conductivity_dyn.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Electron entropy:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electron_entropy'//trim(adjustl(sh_cmd))
   call gnu_entropy(numpar, File_name, file_electron_entropy, t0, t_last, 'OUTPUT_electron_entropy.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron temperatures and chemical potential (for band-resolved calculations):
   if (numpar%do_partial_thermal) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_temperatures'//trim(adjustl(sh_cmd))
      call gnu_el_temperatures(numpar, File_name, file_Te, t0, t_last, &
               'OUTPUT_electron_temperatures.'//trim(adjustl(numpar%fig_extention))) ! below

      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_chempotentials'//trim(adjustl(sh_cmd))
      call gnu_chempots(numpar, File_name, file_mu, t0, t_last, &
               'OUTPUT_electron_chempotentials.'//trim(adjustl(numpar%fig_extention))) ! below
   endif


   ! Atomic temperatures (various definitions):
   if (numpar%print_Ta) then
      ! Atomic entropy:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_entropy'//trim(adjustl(sh_cmd))
      call gnu_entropy_atomic(numpar, File_name, file_atomic_entropy, t0, t_last, &
            'OUTPUT_atomic_entropy.'//trim(adjustl(numpar%fig_extention))) ! below

      File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures'//trim(adjustl(sh_cmd))
      call gnu_at_temperatures(numpar, File_name, file_atomic_temperatures, t0, t_last, &
               'OUTPUT_atomic_temperatures.'//trim(adjustl(numpar%fig_extention))) ! below

      File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_partial'//trim(adjustl(sh_cmd))
      call gnu_at_temperatures_part(numpar, File_name, file_atomic_temperatures_part, t0, t_last, &
               'OUTPUT_atomic_temperatures_partial.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Testmode additional data:
   if (numpar%save_testmode) then   ! if we want to gnuplot testmode data...
      !File_name  = trim(adjustl(file_path))//'OUTPUT_center_of_mass'//trim(adjustl(sh_cmd))
      !call gnu_center_of_mass(File_name, file_testmode, t0, t_last, 'OUTPUT_center_of_mass.'//trim(adjustl(numpar%fig_extention))) ! below
   endif


   ! Electron-ion coupling parameter:
   File_name  = trim(adjustl(file_path))//'OUTPUT_coupling_parameter'//trim(adjustl(sh_cmd))
   call gnu_coupling(numpar, File_name, file_electron_properties, t0, t_last, 'OUTPUT_coupling.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Partial electron-ion coupling:
   call gnu_partial_coupling(Scell(1), matter, numpar, trim(adjustl(file_coupling)), t0, t_last, &
            'OUTPUT_coupling_by_element', 'OUTPUT_coupling_by_orbital')      ! below

   ! Volume:
   File_name  = trim(adjustl(file_path))//'OUTPUT_volume_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_volume(numpar, File_name, file_supercell, t0, t_last, 'OUTPUT_volume.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Mulliken charges:
   if (numpar%Mulliken_model >= 1) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_Mulliken_charges_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_Mulliken_charges(matter, numpar, File_name, file_electron_properties, t0, t_last, &
      'OUTPUT_Mulliken_charges.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Nearest neighbors:
   if (numpar%save_NN) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_neighbors_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_nearest_neighbors(numpar, File_name, file_NN, t0, t_last, 'OUTPUT_nearest_neighbors.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Element-specific nearest neighbors:
   if (allocated(numpar%NN_radii)) then
      do i = 1, size(numpar%NN_radii) ! for all requested elements
         File_name  = trim(adjustl(file_path))//'OUTPUT_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'_Gnuplot'//trim(adjustl(sh_cmd))
         call gnu_nearest_neighbors_elements(numpar, File_name, file_element_NN(i), trim(adjustl(numpar%NN_radii(i)%Name)), &
              matter, t0, t_last, &
              'OUTPUT_nearest_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'.'//trim(adjustl(numpar%fig_extention))) ! below
      enddo ! i
   endif


   ! Pair correlation function:
   if (numpar%save_PCF) then
      ! Pair correlation function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_pair_correlation_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 1.0d0, 'PCF', 'Radius (A)', 'Pair correlation function (a.u.)', 'OUTPUT_pair_correlation.gif', numpar%path_sep, setkey=0)
      call write_pair_correlation_gnuplot(FN, Scell, numpar, matter, 'OUTPUT_pair_correlation_function.dat')   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)
   endif


   ! Distribution function of electrons:
   if (numpar%save_fe) then
      !! Find order of the number, and set number of tics as tenth of it:
      !call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

      ! Distribution function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_distribution_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 5.0d0, 'Distribution', 'Energy (eV)', 'Electron distribution (a.u.)', 'OUTPUT_electron_distribution.gif', numpar%path_sep, setkey=0)
      call write_distribution_gnuplot(FN, Scell, numpar, 'OUTPUT_electron_distribution.dat')   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)
   endif

   ! Orbital-resoloved distribution function of electrons:
   if (numpar%save_fe_orb) then
      ! Distribution function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_fe_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 5.0d0, 'Distribution', 'Energy (eV)', 'Electron distribution (a.u.)', 'OUTPUT_orbital_resolved_fe.gif', numpar%path_sep, setkey=0)
      call write_orb_distribution_gnuplot(FN, Scell, numpar, matter, 'OUTPUT_electron_distribution.dat')   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)
   endif

   ! Distribution function of all electrons on the grid:
   if (numpar%save_fe_grid) then
      ! Find order of the max energy grid, and set number of tics as tenth of it:
      call order_of_time((Scell(1)%E_fe_grid(size(Scell(1)%E_fe_grid)) - 30.0), time_order, temp, x_tics)	! module "Little_subroutines"

      ! Distribution function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_distribution_on_grid_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, x_tics, 'Distribution', 'Energy (eV)', 'Electron density (1/(V*E))', 'OUTPUT_electron_distribution_on_grid.gif', numpar%path_sep, setkey=0)
      call write_distribution_on_grid_gnuplot(FN, Scell, numpar, 'OUTPUT_electron_distribution_on_grid.dat')   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)
   endif


   ! Distribution function of atoms:
   if (numpar%save_fa) then
      call order_of_time((Scell(1)%Ea_grid_out(size(Scell(1)%Ea_grid_out))), time_order, temp, x_tics)	! module "Little_subroutines"

      ! Distribution function can only be plotted as animated gif:

      ! 1) Distribution of kinetic energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atoms_distribution_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, x_tics, 'Distribution', 'Energy (eV)', 'Atomic distribution (a.u.)', 'OUTPUT_atomic_distribution.gif', numpar%path_sep, setkey=0)
      call write_atomic_distribution_gnuplot(FN, Scell, numpar, 'OUTPUT_atomic_distribution.dat')   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)

      ! 2) Distribution of potential energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atoms_distribution_pot_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, x_tics, 'Distribution', 'Energy (eV)', 'Atomic distribution (a.u.)', 'OUTPUT_atomic_distribution_pot.gif', numpar%path_sep, setkey=0)
      call write_atomic_distribution_gnuplot(FN, Scell, numpar, &
            'OUTPUT_atomic_distribution_pot.dat', its_pot=.true. )   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)

      ! 3) Distribution of total energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atoms_distribution_tot_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, x_tics, 'Distribution', 'Energy (eV)', 'Atomic distribution (a.u.)', 'OUTPUT_atomic_distribution_tot.gif', numpar%path_sep, setkey=0)
      call write_atomic_distribution_gnuplot(FN, Scell, numpar, &
            'OUTPUT_atomic_distribution_tot.dat', its_pot=.true., no_maxwell=.true.)   ! below
      call write_gnuplot_script_ending(numpar, FN, File_name, 1)
      close(FN)
   endif


   ! DOS of electrons:
   if (numpar%save_DOS) then  ! Material DOS
      select case (numpar%DOS_splitting)
      case (1) ! with partial DOS
         ! DOS can only be plotted as animated gif:
         File_name  = trim(adjustl(file_path))//'OUTPUT_DOS_Gnuplot'//trim(adjustl(sh_cmd))
         open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
         call write_gnuplot_script_header_new(FN, 6, 1.0d0, 5.0d0, 'DOS', 'Energy (eV)', 'DOS (electrons/eV)', 'OUTPUT_DOS.gif', numpar%path_sep, setkey=0)
         call write_DOS_gnuplot(FN, Scell, numpar, matter, 'OUTPUT_DOS.dat')   ! below
         call write_gnuplot_script_ending(numpar, FN, File_name, 1)
         close(FN)
      endselect
   endif

   ! Optical coefficients
   if (numpar%do_drude) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_optical_coefficients'//trim(adjustl(sh_cmd))
      call gnu_optical_coefficients(numpar, File_name, file_optics, t0, t_last, &
            'OUTPUT_optical_coefficients.'//trim(adjustl(numpar%fig_extention))) ! below
      ! also n and k:
      File_name  = trim(adjustl(file_path))//'OUTPUT_optical_n_and_k'//trim(adjustl(sh_cmd))
      call gnu_n_and_k(numpar, File_name, file_optics, t0, t_last, 'OUTPUT_optical_n_and_k.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   !ccccccccccccccccccccccccccccc
   ! Create also convolved plots:
   CONV:if (Scell(1)%eps%tau > 0.0d0) then ! convolved files too:
      ! Energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_energies_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_energies(numpar, File_name, trim(adjustl(file_energies(1:len(trim(adjustl(file_energies)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_energies_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Temepratures:
      File_name  = trim(adjustl(file_path))//'OUTPUT_temperatures_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_temperatures(numpar, matter, File_name, trim(adjustl(file_temperatures(1:len(trim(adjustl(file_temperatures)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_temepratures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Mean displacement:
      File_name  = trim(adjustl(file_path))//'OUTPUT_mean_displacement_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      if (abs(numpar%MSD_power) > 1.0e-6) then ! only plot it if it's not zero
         call gnu_MSD(matter, numpar, File_name, &
            trim(adjustl(file_temperatures(1:len(trim(adjustl(file_temperatures)))-4)))//'_CONVOLVED.dat', t0, t_last, &
            'OUTPUT_mean_displacement_CONVOLVED.'//trim(adjustl(numpar%fig_extention)), numpar%MSD_power) ! below
      endif

      ! Atomic masks for sectional displacements:
      if (allocated(Scell(1)%Displ)) then
         Nsiz = size(Scell(1)%Displ)   ! how many masks
         do i = 1, Nsiz    ! for all masks
            File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(i)%mask_name))// &
                  '_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
            call gnu_displacements(numpar, File_name, trim(adjustl(file_sect_displ(i)(1:len(trim(adjustl(file_sect_displ(i))))-4) )), &
                  t0, t_last, 'OUTPUT_mean_displacement_'// &
                  trim(adjustl(Scell(1)%Displ(i)%mask_name))//'_CONVOLVED.' &
                  //trim(adjustl(numpar%fig_extention)), Scell(1)%Displ(i)%MSD_power) ! below
            ! Partial by elements, if there is more than one:
            File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(i)%mask_name))// &
                  '_partial_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
            if (matter%N_KAO > 1) then
               call gnu_displacements_partial(File_name, trim(adjustl(file_sect_displ(i)(1:len(trim(adjustl(file_sect_displ(i))))-4) )), &
                  t0, t_last, 'OUTPUT_mean_displacement_'// &
                  trim(adjustl(Scell(1)%Displ(i)%mask_name))//'_partial_CONVOLVED.' &
                  //trim(adjustl(numpar%fig_extention)), Scell(1)%Displ(i)%MSD_power, matter, numpar) ! below
            endif
         enddo ! i
      endif

      ! Diffraction:
      if (numpar%save_diff_peaks) then
         ! Diffraction peaks:
         File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_diffraction_peaks(Scell(1), numpar, File_name, &
            trim(adjustl(file_diffraction_peaks(1:len(trim(adjustl(file_diffraction_peaks)))-4)))//'_CONVOLVED.dat' , &
            t0, t_last, 'OUTPUT_diffraction_peaks_CONVOLVED.'//trim(adjustl(numpar%fig_extention)), .false.) ! below

         ! Check if Debye-Waller analysis is required:
         if ( abs(numpar%DW_theta) > 1.0d-6 ) then
            File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_DW_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
            call gnu_diffraction_peaks(Scell(1), numpar, File_name, &
               trim(adjustl(file_diffraction_peaks_DW(1:len(trim(adjustl(file_diffraction_peaks_DW)))-4)))//'_CONVOLVED.dat' , &
               t0, t_last, 'OUTPUT_diffraction_peaks_DW_CONVOLVED.'//trim(adjustl(numpar%fig_extention)), .true.) ! below

            File_name  = trim(adjustl(file_path))//'OUTPUT_Debye_temperatures_Gnuplot_COLVOLVED'//trim(adjustl(sh_cmd))
            call gnu_Debye_temperatures(Scell(1), numpar, File_name, &
               trim(adjustl(file_Debye_temperature(1:len(trim(adjustl(file_Debye_temperature)))-4)))//'_CONVOLVED.dat' , &
               t0, t_last, 'OUTPUT_Debye_temperatures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
         endif
      endif



      ! Pressure:
      File_name  = trim(adjustl(file_path))//'OUTPUT_pressure_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_pressure(numpar, File_name, trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_pressure_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Stress tensor:
      File_name  = trim(adjustl(file_path))//'OUTPUT_stress_tensor_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_stress(numpar, File_name, trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat', &
           t0, t_last, 'OUTPUT_pressure_tensor_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Numbers of particles:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_and_holes_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_numbers(numpar, File_name, trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electrons_holes_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! High-energy electrons:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_high_energy_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_high_energy_el(numpar, File_name, trim(adjustl(file_high_e(1:len(trim(adjustl(file_high_e)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_electrons_high_energy_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Orbital-resolved electron parameters:
      File_name  = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_orbital_resolved(Scell(1), matter, numpar, File_name, &
               trim(adjustl(file_orb(1:len(trim(adjustl(file_orb)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_orbital_resolved_Ne_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Numbers of CB electrons:
      File_name  = trim(adjustl(file_path))//'OUTPUT_CB_electrons_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_CB_electrons(numpar, File_name, trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat', &
           t0, t_last, 'OUTPUT_CB_electrons_CONVOLVED.'//trim(adjustl(numpar%fig_extention)))

      ! Numbers of deep-shell holes:
      File_name  = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_holes(numpar, File_name, trim(adjustl(file_deep_holes(1:len(trim(adjustl(file_deep_holes)))-4)))//'_CONVOLVED.dat', &
           t0, t_last, matter, 'OUTPUT_deep_shell_holes_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Band gap:
      File_name  = trim(adjustl(file_path))//'OUTPUT_Egap_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_Egap(numpar, File_name, &
           trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', &
           t0, t_last, 'OUTPUT_Egap_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Chemical potential and Ne:
      File_name  = trim(adjustl(file_path))//'OUTPUT_mu_and_Ne_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_mu(numpar, File_name, &
            trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_mu_and_Ne_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Boundaries of the bands:
      File_name  = trim(adjustl(file_path))//'OUTPUT_bands_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_Ebands(numpar, File_name, &
            trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_bands_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron heat capacity:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_Ce_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_capacity(numpar, File_name, &
            trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_electron_Ce_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron heat conductivity:
      if (numpar%do_kappa) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_heat_conductivity(numpar, File_name, &
            trim(adjustl(file_heat_capacity(1:len(trim(adjustl(file_heat_capacity)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_electron_heat_conductivity.'//trim(adjustl(numpar%fig_extention))) ! below
      endif
      if (numpar%do_kappa_dyn) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_dyn_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_heat_conductivity_dyn(numpar, File_name, &
            trim(adjustl(file_heat_capacity_dyn(1:len(trim(adjustl(file_heat_capacity_dyn)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_electron_heat_conductivity_dyn.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Electron entropy:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_entropy(numpar, File_name, &
            trim(adjustl(file_electron_entropy(1:len(trim(adjustl(file_electron_entropy)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_electron_entropy_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron temperatures and chemical potential (for band-resolved calculations):
      if (numpar%do_partial_thermal) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_temperatures_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_el_temperatures(numpar, File_name, trim(adjustl(file_Te(1:len(trim(adjustl(file_Te)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_electron_temperatures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_chempotentials_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_chempots(numpar, File_name, trim(adjustl(file_mu(1:len(trim(adjustl(file_mu)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_electron_chempotentials_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Atomic temperatures (various definitions):
      if (numpar%print_Ta) then
          ! Atomic entropy:
         File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_entropy_atomic(numpar, File_name, &
            trim(adjustl(file_atomic_entropy(1:len(trim(adjustl(file_atomic_entropy)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_atomic_entropy_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_at_temperatures(numpar, File_name, &
            trim(adjustl(file_atomic_temperatures(1:len(trim(adjustl(file_atomic_temperatures)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_atomic_temperatures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_partial_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_at_temperatures_part(numpar, File_name, &
            trim(adjustl(file_atomic_temperatures_part(1:len(trim(adjustl(file_atomic_temperatures_part)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_atomic_temperatures_partial.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Electron-ion coupling parameter:
      File_name  = trim(adjustl(file_path))//'OUTPUT_coupling_parameter_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_coupling(numpar, File_name, &
            trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_coupling_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Partial coupling:
      call gnu_partial_coupling(Scell(1), matter, numpar, trim(adjustl(file_coupling)), t0, t_last, &
            'OUTPUT_coupling_by_element', 'OUTPUT_coupling_by_orbital', convolved=.true.)      ! below

      ! Volume:
      File_name  = trim(adjustl(file_path))//'OUTPUT_volume_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_volume(numpar, File_name, trim(adjustl(file_supercell(1:len(trim(adjustl(file_supercell)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_volume_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Mulliken charges:
      if (numpar%Mulliken_model >= 1) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_Mulliken_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_Mulliken_charges(matter, numpar, File_name, &
         trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', &
         t0, t_last, 'OUTPUT_Mulliken_charges_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Nearest neighbors:
      if (numpar%save_NN) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_neighbors_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_nearest_neighbors(numpar, File_name, trim(adjustl(file_NN(1:len(trim(adjustl(file_NN)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_nearest_neighbors_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      if (numpar%do_drude) then
         ! optical coefficients:
         File_name  = trim(adjustl(file_path))//'OUTPUT_optical_coefficients_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_optical_coefficients(numpar, File_name, &
            trim(adjustl(file_optics(1:len(trim(adjustl(file_optics)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_optical_coefficients_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
         ! also n and k:
         File_name  = trim(adjustl(file_path))//'OUTPUT_optical_n_and_k_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_n_and_k(numpar, File_name, trim(adjustl(file_optics(1:len(trim(adjustl(file_optics)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_optical_n_and_k_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif
   endif CONV
end subroutine create_gnuplot_scripts



subroutine call_vs_slash(path_sep, text, FN) ! how a file is executed: "call " vs "./" (windows vs linux)
   character(*), intent(in) :: path_sep ! 0=windows, 1=linux
   character(*), intent(in) :: text
   integer, intent(in) :: FN ! file into which we write the text
   if (path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,a,a)') 'call ', trim(adjustl(text)), '.cmd'
   else ! It is linux
      write(FN, '(a,a,a)') './', trim(adjustl(text)), '.sh'
   endif
end subroutine call_vs_slash





subroutine gnu_energies(numpar, File_name, file_energies, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_energies ! input file
   real(8), intent(in) :: t0, t_last	! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   real(8) :: x_tics
   character(8) :: temp, time_order
   integer :: FN
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics,  'Energies', 'Time (fs)', 'Energy (eV/atom)',  trim(adjustl(eps_name)), numpar%path_sep, 1)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_energies)), ' "u 1:4 w l lw LW title "Potential energy" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_energies)), ' "u 1:6 w l lw LW title "Atomic energy" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_energies)), ' "u 1:7 w l lw LW title "Energy of atoms and electrons" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_energies)), ' "u 1:8 w l lw LW title "Total energy" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_energies)), '\"u 1:4 w l lw \"$LW\" title \"Potential energy\" ,\'
      write(FN, '(a,a,a)') '\"', trim(adjustl(file_energies)), '\"u 1:6 w l lw \"$LW\" title \"Atomic energy\" ,\'
      write(FN, '(a,a,a)') '\"', trim(adjustl(file_energies)), '\"u 1:7 w l lw \"$LW\" title \"Energy of atoms and electrons\" ,\'
      write(FN, '(a,a,a)') '\"', trim(adjustl(file_energies)), '\"u 1:8 w l lw \"$LW\" title \"Total energy\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_energies


subroutine gnu_temperatures(numpar, matter, File_name, file_temperatures, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter ! parameters of the material
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_temperatures ! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN, i
   real(8) :: x_tics
   character(8) :: temp, time_order
   character(25) :: chtemp

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   !call write_gnuplot_script_header(FN, 1, 3, 'Temperatures','Time (fs)', 'Temperatures (K)', trim(adjustl(file_path))//'OUTPUT_temepratures.'//trim(adjustl(numpar%fig_extention)))
   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Temperatures','Time (fs)', 'Temperatures (K)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Temperatures', 'Time (fs)', 'Temperature (K)', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (matter%N_KAO == 1) then
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_temperatures)), ' "u 1:2 w l lw LW title "Electrons" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(file_temperatures)), ' "u 1:3 w l lw LW title "Atoms" '
!          write(FN, '(a,a,a)') ' "', trim(adjustl(file_temperatures)), ' "u 1:3 w l lw LW title "Atoms (Tkin)" ,\'
!          write(FN, '(a,a,a)') ' "', trim(adjustl(file_temperatures)), ' "u 1:4 w l lt rgb "#0000FF" lw LW title "Atoms (Tconfig)" '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_temperatures)), '\"u 1:2 w l lw \"$LW\" title \"Electrons\" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(file_temperatures)), '\"u 1:3 w l lw \"$LW\" title \"Atoms\" '
!          write(FN, '(a,a,a)') '\"', trim(adjustl(file_temperatures)), '\"u 1:3 w l lw \"$LW\" title \"Atoms (Tkin)\" ,\'
!          write(FN, '(a,a,a)') '\"', trim(adjustl(file_temperatures)), '\"u 1:4 w l lt rgb \"#0000FF\" lw  \"$LW\" title \"Atoms (Tconfig)\" '
      endif
   else ! more than one element:
       if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_temperatures)), ' "u 1:2 w l lw LW title "Electrons" ,\'
         write(FN, '(a,a,a)') ' "', trim(adjustl(file_temperatures)), ' "u 1:3 w l lw LW title "Atoms average" ,\'
         do i = 1, matter%N_KAO - 1
            write(chtemp,'(a)') matter%Atoms(i)%Name//' atoms'
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', 3+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         enddo
         write(chtemp,'(a)') matter%Atoms(matter%N_KAO)%Name//' atoms'
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', 3+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_temperatures)), '\"u 1:2 w l lw \"$LW\" title \"Electrons\" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(file_temperatures)), '\"u 1:3 w l lw \"$LW\" title \"Atoms average\" ,\'
         do i = 1, matter%N_KAO - 1
            write(chtemp,'(a)') matter%Atoms(i)%Name//' atoms'
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', 3+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) ,'\" ,\'
         enddo
         write(chtemp,'(a)') matter%Atoms(matter%N_KAO)%Name//' atoms'
         write(FN, '(a,i3,a,a)') '\"\" u 1:', 3+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//'\"'
      endif
   endif

   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_temperatures


subroutine gnu_MSD(matter, numpar, File_name, file_MSD, t0, t_last, eps_name, MSD_power)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_MSD	! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer, intent(in) :: MSD_power ! power of MSD
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp, MSD_text

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   if (MSD_power > 1) then
      ! Power of MSD:
      write(MSD_text,'(i2)') MSD_power

      !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Mean square displacement','Time (fs)', 'Mean square displacement (A^2)', trim(adjustl(eps_name)))
      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement','Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), numpar%path_sep, 1)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement','Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), numpar%path_sep, 1)	! module "Gnuplotting"
   endif

   if (matter%N_KAO == 1) then
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:5 w l lw LW title "Displacement" '
!          write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:5 w l lw LW title "Displacement" '	! if included Tconf
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:5 w l lw \"$LW\" title \"Displacement\" '
      endif
   else ! more than one element:
      i_start = 4 + matter%N_KAO
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:', i_start ,' w l lw LW title "Average" ,\'
         do i = 1, matter%N_KAO - 1
            write(chtemp,'(a)') matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         enddo
         write(chtemp,'(a)') matter%Atoms(matter%N_KAO)%Name
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " '
      else
         write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:', i_start, ' w l lw \"$LW\" title \"Average\" ,\'
         do i = 1, matter%N_KAO - 1
            write(chtemp,'(a)') matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) ,'\" ,\'
         enddo
         write(chtemp,'(a)') matter%Atoms(matter%N_KAO)%Name
         write(FN, '(a,i3,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//'\"'
      endif
   endif ! (matter%N_KAO == 1)


   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_MSD



subroutine gnu_displacements(numpar, File_name, file_MSD, t0, t_last, eps_name, MSD_power)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_MSD	! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   real(8), intent(in) :: MSD_power ! power of MSD
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp, MSD_text

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   if (MSD_power > 1) then
      ! Power of MSD:
      write(MSD_text,'(i2)') int(MSD_power)

      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), numpar%path_sep, 2)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), numpar%path_sep, 2)	! module "Gnuplotting"
   endif

   i_start = 2
   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:', i_start ,' w l lw LW title "Average" ,\'
      write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+1 ,' w l lw LW title " ', 'X'  ,' " ,\'
      write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+2 ,' w l lw LW title " ', 'Y'  ,' " ,\'
      write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+3 ,' w l lw LW title " ', 'Z'  ,' " '
   else
      write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:', i_start, &
                                          ' w l lw \"$LW\" title \"Average\" ,\'
      write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+1, ' w l lw \"$LW\" title \" ', 'X' ,'\" ,\'
      write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+2, ' w l lw \"$LW\" title \" ', 'Y' ,'\" ,\'
      write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+3, ' w l lw \"$LW\" title \" ', 'Z' ,'\"'
   endif

   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_displacements



subroutine gnu_displacements_partial(File_name, file_MSD, t0, t_last, eps_name, MSD_power, matter, numpar)
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_MSD	! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   real(8), intent(in) :: MSD_power ! power of MSD
   type(Solid), intent(in) :: matter     ! material parameters
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   !------------------------
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp, MSD_text

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   if (MSD_power > 1) then
      ! Power of MSD:
      write(MSD_text,'(i2)') int(MSD_power)

      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), numpar%path_sep, 2)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), numpar%path_sep, 2)	! module "Gnuplotting"
   endif

   i_start = 6
   do i = 1, matter%N_KAO
      chtemp = trim(adjustl(matter%Atoms(i)%Name))
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         if (i == 1) then
            write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:', i_start , &
                     ' w l lw LW title "'//trim(adjustl(chtemp))//'" ,\'
         else
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start + (i-1)*4,' w l lw LW title "'//trim(adjustl(chtemp))//'" ,\'
         endif
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+1+(i-1)*4 ,' w l lw LW title " ', trim(adjustl(chtemp))//':X'  ,' " ,\'
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+2+(i-1)*4,' w l lw LW title " ', trim(adjustl(chtemp))//':Y'  ,' " ,\'
         if (i /= matter%N_KAO) then
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+3+(i-1)*4,' w l lw LW title " ', trim(adjustl(chtemp))//':Z'  ,' " ,\'
         else
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+3+(i-1)*4,' w l lw LW title " ', trim(adjustl(chtemp))//':Z'  ,' " '
         endif

      else
         if (i == 1) then
            write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:', i_start, &
                                          ' w l lw \"$LW\" title \"'//trim(adjustl(chtemp))//'\" ,\'
         else
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+(i-1)*4, ' w l lw \"$LW\" title \"'//trim(adjustl(chtemp))//'\" ,\'
         endif
         write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+1+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':X' ,'\" ,\'
         write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+2+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':Y' ,'\" ,\'

         if (i /= matter%N_KAO) then
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+3+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':Z' ,'\" ,\'
         else
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+3+(i-1)*4, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//':Z' ,'\"'
         endif
      endif
   enddo
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_displacements_partial


subroutine gnu_Debye_temperatures(Scell, numpar, File_name, file_Debye_temperature, t0, t_last, fig_name)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_Debye_temperature ! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: fig_name ! name of the figure
   !------------------------
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp
   character(20) :: peak_name

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Debye temperature', &
            'Time (fs)', 'Debye temperature (K)', trim(adjustl(fig_name)), numpar%path_sep, 0)      ! module "Gnuplotting"

   write(FN, '(a)') 'set encoding utf8'   ! unicode to make teta


   if (size(Scell%diff_peaks%I_diff_peak) == 1) then  ! only one peak to plot
      peak_name = make_diff_peak_name(Scell, 1) ! below
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Debye_temperature)), ' "u 1:2 w l lw LW title "\U+03B8_D'//trim(adjustl(peak_name))//'" '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Debye_temperature)), '\"u 1:2 w l lw \"$LW\" title \"\U+03B8_D'//trim(adjustl(peak_name))//'\" '
      endif
   else ! more than one peak:

      peak_name = make_diff_peak_name(Scell, 1) ! below
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         ! First peak:
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Debye_temperature)), ' "u 1:2 w l lw LW title "\U+03B8_D'//trim(adjustl(peak_name))//'" ,\'
         ! Next peaks:
         do i = 2, size(Scell%diff_peaks%I_diff_peak)-1
            peak_name = make_diff_peak_name(Scell, i) ! below
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', 1+i ,' w l lw LW title "\U+03B8_D', trim(adjustl(peak_name))  ,'" ,\'
         enddo
         ! Last peak:
         peak_name = make_diff_peak_name(Scell, i) ! below
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', 1+i ,' w l lw LW title "\U+03B8_D', trim(adjustl(peak_name))  ,'" '
      else  ! Linux:
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Debye_temperature)), &
               '\"u 1:2 w l lw \"$LW\" title \"\U+03B8_D'//trim(adjustl(peak_name))//'\" ,\'
         do i = 2, size(Scell%diff_peaks%I_diff_peak)-1
            peak_name = make_diff_peak_name(Scell, i) ! below
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', 1+i ,' w l lw \"$LW\" title \"\U+03B8_D', trim(adjustl(peak_name)), '\" ,\'
         enddo
         ! Last peak:
         peak_name = make_diff_peak_name(Scell, i) ! below
         write(FN, '(a,i3,a,a,a)') '\"\" u 1:', 1+i ,' w l lw \"$LW\" title \"\U+03B8_D', trim(adjustl(peak_name)), '\" '
      endif

   endif

   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_Debye_temperatures



subroutine gnu_diffraction_peaks(Scell, numpar, File_name, file_diffraction_peaks, t0, t_last, fig_name, DW)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_diffraction_peaks ! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: fig_name ! name of the figure
   logical, intent(in) :: DW  ! if it is DW, change the axis label
   !------------------------
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp
   character(20) :: peak_name
   character(50) :: y_axis_label

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   if (DW) then
      y_axis_label = 'DW peak intensity (arb. units)'
   else
      y_axis_label = 'Peak intensity (arb. units)'
   endif

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Diffraction peak', &
            'Time (fs)', trim(adjustl(y_axis_label)), trim(adjustl(fig_name)), numpar%path_sep, 0)      ! module "Gnuplotting"
            !'Time (fs)', 'Peak intensity (arb. units)', trim(adjustl(fig_name)), numpar%path_sep, 0)      ! module "Gnuplotting"


   if (size(Scell%diff_peaks%I_diff_peak) == 1) then  ! only one peak to plot
      peak_name = make_diff_peak_name(Scell, 1) ! below
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_diffraction_peaks)), ' "u 1:2 w l lw LW title "'//trim(adjustl(peak_name))//'" '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_diffraction_peaks)), '\"u 1:2 w l lw \"$LW\" title \"'//trim(adjustl(peak_name))//'\" '
      endif
   else ! more than one peak:

      peak_name = make_diff_peak_name(Scell, 1) ! below
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         ! First peak:
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_diffraction_peaks)), ' "u 1:2 w l lw LW title "'//trim(adjustl(peak_name))//'" ,\'
         ! Next peaks:
         do i = 2, size(Scell%diff_peaks%I_diff_peak)-1
            peak_name = make_diff_peak_name(Scell, i) ! below
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', 1+i ,' w l lw LW title "', trim(adjustl(peak_name))  ,'" ,\'
         enddo
         ! Last peak:
         peak_name = make_diff_peak_name(Scell, i) ! below
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', 1+i ,' w l lw LW title "', trim(adjustl(peak_name))  ,'" '
      else  ! Linux:
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_diffraction_peaks)), '\"u 1:2 w l lw \"$LW\" title \"'//trim(adjustl(peak_name))//'\" ,\'
         do i = 2, size(Scell%diff_peaks%I_diff_peak)-1
            peak_name = make_diff_peak_name(Scell, i) ! below
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', 1+i ,' w l lw \"$LW\" title \"', trim(adjustl(peak_name)), '\" ,\'
         enddo
         ! Last peak:
         peak_name = make_diff_peak_name(Scell, i) ! below
         write(FN, '(a,i3,a,a,a)') '\"\" u 1:', 1+i ,' w l lw \"$LW\" title \"', trim(adjustl(peak_name)), '\" '
      endif

   endif

   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_diffraction_peaks


function make_diff_peak_name(Scell, i) result(peak_name)
   character(20) :: peak_name
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: i   ! number of the peak
   !---------------
   character(5) :: text1, text2, text3
   write(text1, '(i0)') Scell%diff_peaks%ijk_diff_peak(1,i)
   write(text2, '(i0)') Scell%diff_peaks%ijk_diff_peak(2,i)
   write(text3, '(i0)') Scell%diff_peaks%ijk_diff_peak(3,i)
   peak_name = '('//trim(adjustl(text1))//trim(adjustl(text2))//trim(adjustl(text3))//')'
end function make_diff_peak_name



subroutine gnu_Mulliken_charges(matter, numpar, File_name, file_electron_properties, t0, t_last, eps_name)
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_properties	! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Mulliken charge','Time (fs)', 'Mulliken charge', trim(adjustl(eps_name)), numpar%path_sep, 1)	! module "Gnuplotting"

   if (matter%N_KAO == 1) then
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:11 w l lw LW title "Charge" '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:11 w l lw \"$LW\" title \"Charge\" '
      endif
   else ! more than one element:
      i_start = 10
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(chtemp,'(a)') matter%Atoms(1)%Name ! first element
         write(FN, '(a,es25.16,a,a,a,i3,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), '" u 1:', i_start+1 ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         do i = 2, matter%N_KAO - 1   ! intermediate elements
            write(chtemp,'(a)') matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         enddo
         write(chtemp,'(a)') matter%Atoms(matter%N_KAO)%Name    ! last element
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " '
      else
         write(chtemp,'(a)') matter%Atoms(1)%Name ! first element
         write(FN, '(a,es25.16,a,a,a,i3,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:', i_start+1, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) ,'\" ,\'
         do i = 2, matter%N_KAO - 1   ! intermediate elements
            write(chtemp,'(a)') matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) , '\" ,\'
         enddo
         write(chtemp,'(a)') matter%Atoms(matter%N_KAO)%Name  ! last element
         write(FN, '(a,i3,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//'\"'
      endif
   endif ! (matter%N_KAO == 1)

   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_Mulliken_charges


subroutine gnu_nearest_neighbors(numpar, File_name, file_NN, t0, t_last, fig_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name    ! file to create
   character(*), intent(in) :: file_NN      ! data file
   real(8), intent(in) :: t0, t_last    ! time instance [fs]
   character(*), intent(in) :: fig_name ! name of the figure
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Nearest neighbors','Time (fs)', 'Nearest neighbors fraction', trim(adjustl(fig_name)), numpar%path_sep, 1)	! module "Gnuplotting"

   if (numpar%path_sep == '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_NN)), ' "u 1:3 w l lw LW title "Single atom"  ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:4 w l lw LW title "One neighbor" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:5 w l lw LW title "Two neighbors" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:6 w l lw LW title "Three neighbors" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:7 w l lw LW title "Four neighbors" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:8 w l lw LW title "Five neighbors" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:9 w l lw LW title "Six neighbors" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_NN)), '\"u 1:3 w l lw \"$LW\" title \"Single atom\"  ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:4 w l lw \"$LW\" title \"One neighbor\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:5 w l lw \"$LW\" title \"Two neighbors\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:6 w l lw \"$LW\" title \"Three neighbors\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:7 w l lw \"$LW\" title \"Four neighbors\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:8 w l lw \"$LW\" title \"Five neighbors\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:9 w l lw \"$LW\" title \"Six neighbors\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_nearest_neighbors




subroutine gnu_nearest_neighbors_elements(numpar, File_name, file_NN, Name, matter, t0, t_last, fig_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name    ! file to create
   character(*), intent(in) :: file_NN      ! data file
   character(*), intent(in) :: Name       ! element name
   type(solid), intent(in) :: matter      ! material parameters
   real(8), intent(in) :: t0, t_last      ! time instance [fs]
   character(*), intent(in) :: fig_name   ! name of the figure
   !------------------------
   integer :: FN, i, i_start
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Nearest neighbors', 'Time (fs)', &
      'Nearest neighbors of '//trim(adjustl(Name)), &
      trim(adjustl(fig_name)), numpar%path_sep, 1)	! module "Gnuplotting"

   if (numpar%path_sep == '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_NN)), ' " u 1:2 w l lw LW title "Total"  ,\'

      do i = 1, matter%N_KAO-1
         write(temp, '(i0)') 2+i    ! column with data
         write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' " u 1:'//trim(adjustl(temp))// &
                              ' w l lw LW title "'//trim(adjustl(matter%Atoms(i)%Name))//'" ,\'
      enddo
      i = matter%N_KAO
      write(temp, '(i0)') matter%N_KAO+2    ! last column with data
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_NN)), ' "u 1:'//trim(adjustl(temp))// &
                              ' w l lw LW title "'//trim(adjustl(matter%Atoms(i)%Name))//'"'
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_NN)), '\" u 1:2 w l lw \"$LW\" title \"Total\"  ,\'
      do i = 1, matter%N_KAO-1
         write(temp, '(i0)') 2+i    ! column with data
         write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\" u 1:'//trim(adjustl(temp))// &
                              ' w l lw \"$LW\" title \"'//trim(adjustl(matter%Atoms(i)%Name))//'\" ,\'
      enddo
      i = matter%N_KAO
      write(temp, '(i0)') matter%N_KAO+2    ! last column with data
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_NN)), '\"u 1:'//trim(adjustl(temp))// &
                              ' w l lw \"$LW\" title \"'//trim(adjustl(matter%Atoms(i)%Name))//'\"'
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_nearest_neighbors_elements


subroutine gnu_pressure(numpar, File_name, file_pressure, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_pressure ! input file
   real(8), intent(in) :: t0, t_last	! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Pressure','Time (fs)', 'Pressure (GPa)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Pressure','Time (fs)', 'Pressure (GPa)', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_pressure)), ' "u 1:2 w l lw LW title "Pressure" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_pressure)), '\"u 1:2 w l lw \"$LW\" title \"Pressure\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_pressure


subroutine gnu_stress(numpar, File_name, file_pressure, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_pressure ! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Stress tensor','Time (fs)', 'Stress tensor (GPa)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Pressure tensor','Time (fs)', 'Pressure tensor (GPa)', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_pressure)), ' "u 1:3 w l lw LW title "Pressure (x,x)"  ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:4 w l lw LW title "Pressure (x,y)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:5 w l lw LW title "Pressure (x,z)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:6 w l lw LW title "Pressure (y,x)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:7 w l lw LW title "Pressure (y,y)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:8 w l lw LW title "Pressure (y,z)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:9 w l lw LW title "Pressure (z,x)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:10 w l lw LW title "Pressure (z,y)" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_pressure)), ' "u 1:11 w l lt rgb "#545454" lw LW title "Pressure (z,z)" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_pressure)), '\"u 1:3 w l lw \"$LW\" title \"Pressure (x,x)\"  ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:4 w l lw \"$LW\" title \"Pressure (x,y)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:5 w l lw \"$LW\" title \"Pressure (x,z)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:6 w l lw \"$LW\" title \"Pressure (y,x)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:7 w l lw \"$LW\" title \"Pressure (y,y)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:8 w l lw \"$LW\" title \"Pressure (y,z)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:9 w l lw \"$LW\" title \"Pressure (z,x)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:10 w l lw \"$LW\" title \"Pressure (z,y)\" ,\'
      write(FN, '(a,a,a)') ' \"', trim(adjustl(file_pressure)), '\"u 1:11 w l lt rgb \"#545454\" lw \"$LW\" title \"Pressure (z,z)\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_stress



subroutine gnu_numbers(numpar, File_name, file_numbers, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_numbers ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Numbers', 'Time (fs)', &
                        'Particles (1/atom)', trim(adjustl(eps_name)), numpar%path_sep, 0) ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_numbers)), ' "u 1:4 w l lw LW title "High-energy electrons" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_numbers)), ' "u 1:5 w l lw LW title "All deep-shell holes"  ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_numbers)), ' "u 1:7 w l lw LW title "Photons" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_numbers)), ' "u 1:6 w l lw LW title "Error in particle conservation" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_numbers)), '\"u 1:4 w l lw \"$LW\" title \"High-energy electrons\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_numbers)), '\"u 1:5 w l lw \"$LW\" title \"All deep-shell holes\"  ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_numbers)), '\"u 1:7 w l lw \"$LW\" title \"Photons\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_numbers)), '\"u 1:6 w l lw \"$LW\" title \"Error in particle conservation\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_numbers



subroutine gnu_high_energy_el(numpar, File_name, file_high_e, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_high_e ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Numbers', 'Time (fs)', &
                        'Electrons (1/atom)', trim(adjustl(eps_name)), numpar%path_sep, 0) ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_high_e)), ' "u 1:2 w l lw LW title "Total" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_high_e)), ' "u 1:3 w l lw LW title "Photo"  ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_high_e)), ' "u 1:4 w l lw LW title "Auger" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_high_e)), ' "u 1:5 w l lw LW title "Impact-ionized" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_high_e)), ' "u 1:6 w l lw LW title "Error" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_high_e)), '\"u 1:2 w l lw \"$LW\" title \"Total\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_high_e)), '\"u 1:3 w l lw \"$LW\" title \"Photo\"  ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_high_e)), '\"u 1:4 w l lw \"$LW\" title \"Auger\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_high_e)), '\"u 1:5 w l lw \"$LW\" title \"Impact-ionization\" '
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_high_e)), '\"u 1:6 w l lw \"$LW\" title \"Error\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_high_energy_el


! Orbital-resolved electron parameters:
subroutine gnu_orbital_resolved(Scell, matter, numpar, File_name, file_orb, t0, t_last, eps_name)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(in) :: numpar   ! numerical parameters, including lists of earest neighbors
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_orb ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   !------------------------
   integer :: FN, i_at, i_types, N_at, N_types, norb, i_col, Nsiz, NKOA
   real(8) :: x_tics
   character(8) :: temp, time_order, chtemp1, chtemp, ch_col

   ! Find number of orbitals per atom:
   NKOA = matter%N_KAO      ! number of kinds of atoms
   N_at = size(Scell%MDatoms) ! number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   norb = Nsiz/N_at ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"


   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Numbers', 'Time (fs)', &
                        'Electrons (1/atom)', trim(adjustl(eps_name)), numpar%path_sep, 1) ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows

      ! All shells resolved:
      i_col = 2   ! to start with
      do i_at = 1, NKOA
         do i_types = 1, N_types
            chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
            write(chtemp,'(a)') trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
            i_col = i_col + 1 ! column number
            write(ch_col, '(i4)') i_col

            if ( (i_at == 1) .and. (i_types == 1) ) then ! first line
               write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_orb)), '" u 1:'//trim(adjustl(ch_col))// &
                        ' w l lw LW title "'//trim(adjustl(chtemp))//'" ,\'
            else ! regular orrbital
               write(FN, '(a,a,a)') ' "', trim(adjustl(file_orb)), '" u 1:'//trim(adjustl(ch_col))//' w l lw LW title "'// &
                        trim(adjustl(chtemp))//'" ,\'
            endif
         enddo   ! i_types
      enddo ! i_at
      ! Last: total Ne:
      write(FN, '(a)') ' "'//trim(adjustl(file_orb))//'" u 1:2 w l lw LW title "Total" '

   else
      ! All shells resolved:
      i_col = 2   ! to start with
      do i_at = 1, NKOA
         do i_types = 1, N_types
            chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
            write(chtemp,'(a)') trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
            i_col = i_col + 1 ! column number
            write(ch_col, '(i4)') i_col

            if ( (i_at == 1) .and. (i_types == 1) ) then ! first line
               write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_orb)), '\" u 1:'//trim(adjustl(ch_col))// &
                        ' w l lw \"$LW\" title \"'//trim(adjustl(chtemp))//'\" ,\'
            else ! regular orrbital
               write(FN, '(a,a,a)') '\"', trim(adjustl(file_orb)), '\" u 1:'//trim(adjustl(ch_col))//' w l lw \"$LW\" title \"'// &
                        trim(adjustl(chtemp))//'\" ,\'
            endif
         enddo   ! i_types
      enddo ! i_at
      ! Last: total Ne:
      write(FN, '(a)') ' \"'//trim(adjustl(file_orb))//'\" u 1:2 w l lw \"$LW\" title \"Total\" '
   endif

   call write_gnuplot_script_ending(numpar, FN, trim(adjustl(File_name)), 1) ! below
   close(FN)
end subroutine gnu_orbital_resolved



subroutine gnu_CB_electrons(numpar, File_name, file_numbers, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_numbers ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'CB_electrons','Time (fs)', 'Electrons (per atom)', trim(adjustl(eps_name)), numpar%path_sep, 1)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_numbers)), ' "u 1:($3/4*100) w l lw LW title "CB electrons" ,\'
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_numbers)), ' "u 1:($3) w l lw LW title "CB electrons" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_numbers)), ' "u 1:($7) w l lw LW title "Photons"'
   else
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_numbers)), '\"u 1:(\$3/4*100) w l lw \"$LW\" title \"CB electrons\" ,\'
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_numbers)), '\"u 1:(\$3) w l lw \"$LW\" title \"CB electrons\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_numbers)), '\"u 1:(\$7) w l lw \"$LW\" title \"Photons\"'
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_CB_electrons


subroutine gnu_holes(numpar, File_name, file_deep_holes, t0, t_last, matter, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_deep_holes ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   type(Solid), intent(in) :: matter
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN, counter, Nshl, i, j, Na, N_sh_max, font_siz, N_sh_tot
   character(100) :: chtemp, ch_temp
   character(11) :: chtemp11
   real(8) :: x_tics
   character(8) :: temp, time_order
   logical :: first_line, last_atom_no_shells

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

    ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   ! Find how many shells we have for plotting:
   Na = size(matter%Atoms)
   N_sh_max = size(matter%Atoms(1)%Ip)
   do i = 1, Na !size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      if (N_sh_max < Nshl) N_sh_max = Nshl ! to find the maximal value
   enddo
   N_sh_tot = Na*N_sh_max  ! estimate the total amount of shells
   if (N_sh_tot > 70) then ! make tiny font
      font_siz = 8
   elseif (N_sh_tot > 45) then ! make very small font
      font_siz = 10
   elseif (N_sh_tot > 26) then ! make small font
      font_siz = 12
   else  ! standard font
      font_siz = 14
   endif

   ! prepare gnuplot header:
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Holes','Time (fs)', 'Particles (total)', &
         trim(adjustl(eps_name)), numpar%path_sep, setkey=0, fontsize=font_siz)  ! module "Gnuplotting"

   counter = 0 ! to start with
   first_line = .true.  ! to start from the first line
   last_atom_no_shells = .false.   ! to start with
   W_vs_L:if (numpar%path_sep .EQ. '\') then	! if it is Windows
     ATOMS0:do i = 1, Na ! size(matter%Atoms) ! for all atoms
         write(ch_temp,'(a,i0)') "dashtype ", i ! to set line type
         Nshl = size(matter%Atoms(i)%Ip)
         SHELLS0:do j = 1, Nshl ! for all shells of this atom
            if ((i .NE. 1) .or. (j .NE. Nshl)) then ! atomic shell:
               counter = counter + 1
               call define_PQN(matter%Atoms(i)%Shl_dsgnr(j), chtemp11) ! module "Dealing_with_EADL"
               write(chtemp,'(a,a,a)') trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(chtemp11))
               select case(Na)
               case (1)
                  !if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] "' , trim(adjustl(file_deep_holes)), &
                        '" u 1:', 1+j ,' w l lw LW '//trim(adjustl(ch_temp))//' title "', trim(adjustl(chtemp))  ,'"'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") ' "', trim(adjustl(file_deep_holes)), '" u 1:', 1+j ,' w l lw LW '// &
                        trim(adjustl(ch_temp))//' title "', trim(adjustl(chtemp))  ,'"'
                  endif
                  if ((i .NE. Na) .OR. (j .LT. Nshl-1)) then
                     write(FN, '(a)') ',\'
                  else
                     write(FN, '(a)') ''
                  endif
               case default
                  !if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] "' , trim(adjustl(file_deep_holes)), &
                        ' "u 1:', 1+counter,' w l lw LW '//trim(adjustl(ch_temp))//' title " ', trim(adjustl(chtemp))  ,' "'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") ' "', trim(adjustl(file_deep_holes)), '" u 1:', 1+counter, &
                        ' w l lw LW '//trim(adjustl(ch_temp))//' title "', trim(adjustl(chtemp))  ,'"'
                  endif

                  ! Check if the next element is last, and if it has any shells:
                  if (i+1 .EQ. Na) then
                     ! check if it has any shells:
                     if (matter%Atoms(Na)%Z < 3) then ! element H and He have no core shells, exclude from plotting
                        last_atom_no_shells = .true.
                     else
                        last_atom_no_shells = .false.
                     endif
                  else
                     last_atom_no_shells = .false.
                  endif

                  if ( ((i .NE. Na) .OR. (j .LT. Nshl)) .and. .not.last_atom_no_shells ) then
                     write(FN, '(a)') ',\'
                  else
                     write(FN, '(a)') ''
                  endif
               end select
            else ! VB:
               !write(chtemp,'(a)') 'Valence Band'
               !skip it
            endif
         enddo SHELLS0
      enddo ATOMS0
   else W_vs_L ! It is linux
      ATOMS:do i = 1, Na   ! size(matter%Atoms) ! for all atoms
         write(ch_temp,'(a,i0)') "dashtype ", i ! to set line type
         Nshl = size(matter%Atoms(i)%Ip)
         SHELLS:do j = 1, Nshl ! for all shells of this atom
            if ((i .NE. 1) .or. (j .NE. Nshl)) then ! atomic shell:
               counter = counter + 1
               call define_PQN(matter%Atoms(i)%Shl_dsgnr(j), chtemp11) ! module "Dealing_with_EADL"
               write(chtemp,'(a,a,a)') trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(chtemp11))
               select case(Na)
               case (1)
!                   if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] \"' , trim(adjustl(file_deep_holes)), &
                        '\"u 1:', 1+j ,' w l lw \"$LW\" '//trim(adjustl(ch_temp))//' title \" ', trim(adjustl(chtemp))  ,'\"'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") '\"', trim(adjustl(file_deep_holes)), '\"u 1:', 1+j , &
                        ' w l lw \"$LW\" '//trim(adjustl(ch_temp))//' title \" ', trim(adjustl(chtemp))  ,'\"'
                  endif
                  if ((i .NE. Na) .OR. (j .LT. Nshl-1)) then
                     write(FN, '(a)') ',\'
                  else
                     write(FN, '(a)') ''
                  endif
               case default
!                   if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] \"' , trim(adjustl(file_deep_holes)), &
                        '\"u 1:', 1+counter,' w l lw \"$LW\" '//trim(adjustl(ch_temp))//' title \" ', trim(adjustl(chtemp))  ,'\"'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") '\"', trim(adjustl(file_deep_holes)), '\"u 1:', 1+counter, &
                        ' w l lw \"$LW\" '//trim(adjustl(ch_temp))//' title \" ', trim(adjustl(chtemp))  ,'\"'
                  endif

                  ! Check if the next element is last, and if it has any shells:
                  if (i+1 .EQ. Na) then
                     ! check if it has any shells:
                     if (matter%Atoms(Na)%Z < 3) then ! element H and He have no core shells, exclude from plotting
                        last_atom_no_shells = .true.
                     else
                        last_atom_no_shells = .false.
                     endif
                  else
                     last_atom_no_shells = .false.
                  endif

                  !if ((i .NE. Na) .OR. (j .LT. Nshl)) then
                  if ( ((i .NE. Na) .OR. (j .LT. Nshl)) .and. .not.last_atom_no_shells ) then
                     write(FN, '(a)') ',\'
                  else
                     write(FN, '(a)') ''
                  endif
               end select
            else ! VB:
               !write(chtemp,'(a)') 'Valence Band'
               !skip it
            endif
         enddo SHELLS
      enddo ATOMS
   endif W_vs_L
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_holes


subroutine gnu_Egap(numpar, File_name, file_electron_properties, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_properties ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Egap','Time (fs)', 'Bandgap (eV)', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), '"u 1:4 w l lw LW title "Band gap" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:4 w l lw \"$LW\" title \"Band gap \" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_Egap



subroutine gnu_mu(numpar, File_name, file_electron_properties, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_properties ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'mu and Ne', &
            'Time (fs)', 'Energy (eV)', trim(adjustl(eps_name)), numpar%path_sep, 0) ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      ! Make the right axis with difference scale and ticks for Ne:
      write(FN, '(a)') 'set y2tics 0.1 nomirror tc lt 2'
      write(FN, '(a)') 'set y2label "Electron density (1/atom)" font "arial,18" '

      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), &
                  ' "u 1:3 w l lw LW title "Chemical potential" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), &
                  ' "u 1:2 w l lw LW axes x1y2 title "Electron density" '
   else ! It is linux
      ! Make the right axis with difference scale and ticks for Ne:
      write(FN, '(a)') 'set y2tics 0.1 nomirror tc lt 2'
      write(FN, '(a)') 'set y2label \"Electron density (1/atom)\" font \"arial,18\" '

      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), &
                  '\"u 1:3 w l lw \"$LW\" title \"Chemical potential\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), &
                  '\"u 1:2 w l lw \"$LW\" axes x1y2 title \"Electron density\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_mu


subroutine gnu_Ebands(numpar, File_name, file_electron_properties, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_properties ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'bands boundaries','Time (fs)', 'Energy (eV)', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:10 w l lw LW title "Top of CB",\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), ' "u 1:9 w l lw LW title "Bottom of CB",\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), ' "u 1:8 w l lw LW title "Top of VB",\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), ' "u 1:7 w l lw LW title "Bottom of VB",\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), ' "u 1:3 w l lt rgb "#000000" lw LW title "Chemical potential" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:10 w l lw \"$LW\" title \"Top of CB \" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:9 w l lw \"$LW\" title \"Bottom of CB\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:8 w l lw \"$LW\" title \"Top of VB\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:7 w l lw \"$LW\" title \"Bottom of VB\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:3 w l lt rgb \"#000000\"  lw \"$LW\" title \"Chemical potential\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_Ebands



subroutine gnu_capacity(numpar, File_name, file_electron_properties, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_properties ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Electron Ce','Time (fs)', 'Heat capacity (J/(m^3 K))', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron Ce','Time (fs)', 'Heat capacity (J/(m^3 K))', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:5 w l lw LW title "Electron heat capacity"  '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:5 w l lw \"$LW\" title \"Electron heat capacity\"  '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_capacity


subroutine gnu_heat_conductivity_dyn(numpar, File_name, file_heat_capacity, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_heat_capacity ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron Kappa','Time (fs)', 'Heat conductivity (W/(m K))', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_heat_capacity)), ' "u 1:2 w l lw LW title "kappa total",\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_heat_capacity)), '" u 1:3 w l lw LW title "kappa e-ph" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_heat_capacity)), '" u 1:4 w l lw LW title "kappa e-e" '

   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_heat_capacity)), '\"u 1:2 w l lw \"$LW\" title \"kappa total\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_heat_capacity)), '\" u 1:3 w l lw \"$LW\" title \"kappa e-ph\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_heat_capacity)), '\" u 1:4 w l lw \"$LW\" title \"kappa e-e\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_heat_conductivity_dyn


subroutine gnu_heat_conductivity(numpar, File_name, file_heat_capacity, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_heat_capacity ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron K','Electron temperature (K)', &
            'Heat conductivity (W/(m K))', trim(adjustl(eps_name)), numpar%path_sep, 0)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_heat_capacity)), &
               ' "u 1:2 w l lw LW title "Electron heat conductivity"  '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_heat_capacity)), &
               '\"u 1:2 w l lw \"$LW\" title \"Electron heat conductivity\"  '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_heat_conductivity


subroutine gnu_entropy(numpar, File_name, file_electron_entropy, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_entropy ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron entropy','Time (fs)', 'Electron entropy (eV/K)', trim(adjustl(eps_name)), numpar%path_sep, 0)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_entropy)), '" u 1:3 w l lw LW title "Equilibrium" ,\'

      if (numpar%do_partial_thermal) then ! partial entropies too
         write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:2 w l lw LW title "Nonequilibrium" ,\'
         write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:5 w l lw LW title "Equilibrium VB" ,\'
         write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:4 w l lw LW title "Nonequilibrium VB" ,\'
         write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:7 w l lw LW title "Equilibrium CB" ,\'
         write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:6 w l lw LW title "Nonequilibrium CB" '
      else  ! only total, no partial
         write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:2 w l lw LW title "Nonequilibrium" '
      endif
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_entropy)), '\" u 1:3 w l lw \"$LW\" title \"Equilibrium\" ,\'

      if (numpar%do_partial_thermal) then ! partial entropies too
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:2 w l lw \"$LW\" title \"Nonequilibrium\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:5 w l lw \"$LW\" title \"Equilibrium VB\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:4 w l lw \"$LW\" title \"Nonequilibrium VB\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:7 w l lw \"$LW\" title \"Equilibrium CB\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:6 w l lw \"$LW\" title \"Nonequilibrium CB\" '
      else  ! only total, no partial
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:2 w l lw \"$LW\" title \"Nonequilibrium\" '
      endif
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_entropy



subroutine gnu_entropy_atomic(numpar, File_name, file_atomic_entropy, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_atomic_entropy ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Atomic entropy', 'Time (fs)', &
         'Atomic entropy (eV/K)', trim(adjustl(eps_name)), numpar%path_sep, 1)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_atomic_entropy)), '" u 1:3 w l lw LW title "Equilibrium" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_atomic_entropy)), '" u 1:4 w l lw LW dashtype 2 title "Equilibrium (num)" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_atomic_entropy)), '" u 1:5 w l lw LW dashtype 4 title "Configurational" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_atomic_entropy)), '" u 1:6 w l lw LW dashtype 5 title "Total" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_atomic_entropy)), '" u 1:2 w l lw LW title "Nonequilibrium" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_atomic_entropy)), '\" u 1:3 w l lw \"$LW\" title \"Equilibrium\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_atomic_entropy)), '\" u 1:4 w l lw \"$LW\" dashtype 2 title \"Equilibrium (num)\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_atomic_entropy)), '\" u 1:5 w l lw \"$LW\" dashtype 4 title \"Configurational\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_atomic_entropy)), '\" u 1:6 w l lw \"$LW\" dashtype 5 title \"Total\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_atomic_entropy)), '\" u 1:2 w l lw \"$LW\" title \"Nonequilibrium\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_entropy_atomic


subroutine gnu_el_temperatures(numpar, File_name, file_Te, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_Te ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron tempereature', 'Time (fs)', 'Electron temperature (K)', trim(adjustl(eps_name)), numpar%path_sep, 0)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Te)), '" u 1:2 w l lw LW title "Total (kinetic)" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Te)), '" u 1:3 w l lw LW title "Valence" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Te)), '" u 1:4 w l lw LW title "Conduction" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Te)), '\" u 1:2 w l lw \"$LW\" title \"Total (kinetic)\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Te)), '\" u 1:3 w l lw \"$LW\" title \"Valence\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Te)), '\" u 1:4 w l lw \"$LW\" title \"Conduction\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_el_temperatures



subroutine gnu_at_temperatures(numpar, File_name, file_Ta, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_Ta ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Atomic tempereature', &
         'Time (fs)', 'Atomic temperature (K)', trim(adjustl(eps_name)), numpar%path_sep, 0)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      ! Not needed old bits:
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Ta)), '" u 1:4 w l lw 1 dashtype 2 title "Distributional" ,\'
      !write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:3 w l lw 1 dashtype 4 title "Entropic" ,\'
      !write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:6 w l lw 1.5 dashtype 5 title "Potential" ,\'
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Ta)), '" u 1:8 w l lw 2 dashtype 2 title "Sine^2" ,\'
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "', trim(adjustl(file_Ta)), '" u 1:11 w l lt rgb "blue" lw LW title "Configurational" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:12 w l lt rgb "green" lw 2 dashtype "__" title "Hyperconfig" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:2 w l lt rgb "black" lw LW title "Kinetic" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:5 w l lt rgb "red" dashtype "_." lw 2 title "Fluctuational" '
   else ! It is linux
      ! Not needed old bits:
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Ta)), '\" u 1:4 w l lw 1 dashtype 2 title \"Distributional\" ,\'
      !write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:3 w l lw 1 dashtype 4 title \"Entropic\" ,\'
      !write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:6 w l lw 1.5 dashtype 5 title \"Potential\" ,\'
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Ta)), '\" u 1:8 w l lw 2 dashtype 2 title \"Sine^2\" ,\'
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Ta)), '\" u 1:11 w l lt rgb \"blue\" lw \"$LW\" title \"Configurational\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:12 w l lt rgb \"green\" lw 2 dashtype \"__\" title \"Hyperconfig\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:2 w l lt rgb \"black\" lw \"$LW\" title \"Kinetic\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:5 w l lt rgb \"red\" dashtype \"_.\" lw 2 title \"Fluctuational\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_at_temperatures



subroutine gnu_at_temperatures_part(numpar, File_name, file_Ta, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_Ta ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Atomic tempereature', &
         'Time (fs)', 'Atomic temperature (K)', trim(adjustl(eps_name)), numpar%path_sep, 0)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Ta)), '" u 1:2 w l lw LW dashtype "__" title "Kinetic: X" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:3 w l lw LW dashtype "_." title "Kinetic: Y" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:4 w l lw LW title "Kinetic: Z" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:5 w l lw LW dashtype "__" title "Virial: X" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:6 w l lt rgb "green" lw LW dashtype "_." title "Virial: Y" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Ta)), '" u 1:7 w l lt rgb "black" lw LW title "Virial: Z" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Ta)), '\" u 1:2 w l lw \"$LW\" dashtype \"__\" title \"Kinetic: X\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:3 w l lw \"$LW\" dashtype \"_.\" title \"Kinetic: Y\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:4 w l lw \"$LW\" title \"Kinetic: Z\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:5 w l lw \"$LW\" dashtype \"__\" title \"Virial: X\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:6 w l lt rgb \"green\" lw \"$LW\" dashtype \"_.\" title \"Virial: Y\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Ta)), '\" u 1:7 w l lt rgb \"black\" lw \"$LW\" title \"Virial: Z\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_at_temperatures_part



subroutine gnu_chempots(numpar, File_name, file_mu, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_mu ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Chemical potential', 'Time (fs)', 'Electron chemical potential (eV)', trim(adjustl(eps_name)), numpar%path_sep, 0)   ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_mu)), '" u 1:2 w l lw LW title "Total" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_mu)), '" u 1:3 w l lw LW title "Valence" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_mu)), '" u 1:4 w l lw LW title "Conduction" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_mu)), '\" u 1:2 w l lw \"$LW\" title \"Total\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_mu)), '\" u 1:3 w l lw \"$LW\" title \"Valence\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_mu)), '\" u 1:4 w l lw \"$LW\" title \"Conduction\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_chempots



subroutine gnu_coupling(numpar, File_name, file_electron_properties, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_electron_properties ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Coupling parameter','Time (fs)', &
      'Coupling parameter (W/(m^3 K))', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), &
         ' "u 1:6 w l lw LW title "Electron-ion coupling" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), &
         '\"u 1:6 w l lw \"$LW\" title \"Electron-ion coupling\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_coupling



subroutine gnu_partial_coupling(Scell, matter, numpar, file_coupling, t0, t_last, script_name_element, script_name_shell, convolved)
   type(Super_cell), intent(in) :: Scell
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(in) :: numpar
   real(8), intent(in) :: t0, t_last
   character(*), intent(in) :: file_coupling, script_name_element, script_name_shell
   logical, intent(in), optional :: convolved

   integer :: FN, i, j, Nat, col1, col2, N_types, norb, i_start, i_orb, j_orb, N_at, N_cols
   character(300) :: File_name, Plot_name, Data_file_name, temp_txt
   character(30) :: x_min, x_max, temp_col1, temp_col2, linestyle, orb_name1, orb_name2
   character(10) :: sh_cmd, call_slash

   ! -----------------------------------------------------------
   ! OS detection: Windows vs Linux
   ! -----------------------------------------------------------
   if (numpar%path_sep .EQ. '\') then
      sh_cmd = '.cmd'
      call_slash = 'call '
   else
      sh_cmd = '.sh'
      call_slash = './'
   endif

   select case (numpar%DOS_splitting)

   case (1)
      Nat = size(matter%Atoms)

      ! ============================================================
      ! 1) ELEMENT-RESOLVED PLOT
      ! ============================================================
      if (Nat > 1) then

         File_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))// &
                     trim(adjustl(script_name_element))//trim(adjustl(sh_cmd))
         open(NEWUNIT=FN, FILE=trim(adjustl(File_name)), action="write", status="replace")

         Plot_name = 'OUTPUT_coupling_by_element'
         Data_file_name = file_coupling

         if (present(convolved)) then
            if (convolved) then
               Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
               Data_file_name = trim(adjustl(file_coupling(1:len(trim(adjustl(file_coupling)))-4)))//'_CONVOLVED.dat'
            endif
         endif

         write(x_min,'(f12.3)') t0
         write(x_max,'(f12.3)') t_last

         ! ============================================================
         ! LINUX VERSION (.sh)
         ! ============================================================
         if (sh_cmd .EQ. '.sh') then
            write(FN,'(a)') '#!/bin/bash'
            write(FN,'(a)') 'gnuplot << EOF'
            write(FN,'(a)') 'set terminal pngcairo enhanced font "Arial,14"'
            write(FN,'(a)') 'set output "'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'"'
            write(FN,'(a)') 'set xlabel "Time (fs)" font "arial,18"'
            write(FN,'(a)') 'set ylabel "Coupling parameter (W/(m^3 K))" font "arial,18"'
            write(FN,'(a)') 'set xrange ['//trim(adjustl(x_min))//':'//trim(adjustl(x_max))//']'
            write(FN,'(a)') 'set key top left'
            write(FN,'(a)') 'unset grid'
            write(FN,'(a)') 'plot \'

         ! ============================================================
         ! WINDOWS VERSION (.cmd)
         ! ============================================================
         else
            write(FN,'(a)') '@echo off & call gnuplot.exe -e "echo='//"'#'"//';set macros" "%~f0" & goto :eof'
            !write(FN,'(a)') 'gnuplot -persist -e "'
            write(FN,'(a)') 'set terminal pngcairo enhanced font ''Arial,14'';'
            write(FN,'(a)') 'set output '''//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//''';'
            write(FN,'(a)') 'set xlabel ''Time (fs)'' font "arial,18";'
            write(FN,'(a)') 'set ylabel ''Coupling parameter (W/(m^3 K))'' font "arial,18";'
            write(FN,'(a)') 'set xrange ['//trim(adjustl(x_min))//':'//trim(adjustl(x_max))//'];'
            write(FN,'(a)') 'set key top left;'
            write(FN,'(a)') 'unset grid;'
            write(FN,'(a)') 'plot \'
         endif

         ! ============================================================
         ! TOTAL COUPLING
         ! ============================================================
         if (sh_cmd .EQ. '.sh') then
            write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:2 with lines lw 3 lt 1 title "Total", \'
         else
            write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:2 with lines lw 3 lt 1 title ''Total'', \'
         endif

         ! ============================================================
         ! ELEMENT-RESOLVED PAIRS
         ! ============================================================
         do i = 1, Nat
            do j = i, Nat

               col1 = 3 + (i-1)*Nat + (j-1)
               col2 = 3 + (j-1)*Nat + (i-1)

               write(temp_col1,'(i0)') col1
               write(temp_col2,'(i0)') col2

               call select_linestyle_gnu(col1, linestyle)   ! module "Gnuplotting"

               if (col1 == col2) then
                  if (sh_cmd .EQ. '.sh') then
                     if ((i == Nat) .and. (j==Nat)) then
                        write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:'//trim(adjustl(temp_col1))// &
                        ' with lines lw 3 '//trim(adjustl(linestyle))//' title "'// &
                        trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//'"'
                     else
                        write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:'//trim(adjustl(temp_col1))// &
                        ' with lines lw 3 '//trim(adjustl(linestyle))//' title "'// &
                        trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//'", \'
                     endif
                  else
                     write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:'//trim(adjustl(temp_col1))// &
                        ' with lines lw 3 '//trim(adjustl(linestyle))//' title '''// &
                        trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//''', \'
                  endif
               elseif (col1 < col2) then ! only print the sum once
                  if (sh_cmd .EQ. '.sh') then
                     write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:(column('//trim(adjustl(temp_col1))// &
                        ')+column('//trim(adjustl(temp_col2))//')) with lines lw 3 '// &
                        trim(adjustl(linestyle))//' title "'// &
                        trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//'", \'
                  else
                     write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:(column('//trim(adjustl(temp_col1))// &
                        ')+column('//trim(adjustl(temp_col2))//')) with lines lw 3 '// &
                        trim(adjustl(linestyle))//' title '''// &
                        trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(j)%Name))//''', \'
                  endif
               endif

            enddo
         enddo

         if (sh_cmd .EQ. '.sh') then
            write(FN,'(a)') 'EOF'
         else
            !write(FN,'(a)') '"'
         endif

         call write_gnuplot_script_ending(numpar, FN, File_name, 1)
         close(FN)
      endif

      ! ============================================================
      ! 2) ORBITAL-RESOLVED PLOT
      ! ============================================================
      N_at = size(Scell%MDatoms)
      norb = size(Scell%Ha,1) / N_at
      N_types = number_of_types_of_orbitals(norb)

      if (N_types > 1) then

         File_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))// &
                     trim(adjustl(script_name_shell))//trim(adjustl(sh_cmd))
         open(NEWUNIT=FN, FILE=trim(adjustl(File_name)), action="write", status="replace")

         Plot_name = 'OUTPUT_coupling_by_orbital'
         Data_file_name = file_coupling

         if (present(convolved)) then
            if (convolved) then
               Plot_name = trim(adjustl(Plot_name))//'_CONVOLVED'
               Data_file_name = trim(adjustl(file_coupling(1:len(trim(adjustl(file_coupling)))-4)))//'_CONVOLVED.dat'
            endif
         endif

         write(x_min,'(f12.3)') t0
         write(x_max,'(f12.3)') t_last

         ! Linux or Windows header
         if (sh_cmd .EQ. '.sh') then
            write(FN,'(a)') '#!/bin/bash'
            write(FN,'(a)') 'gnuplot << EOF'
            write(FN,'(a)') 'set terminal pngcairo enhanced font "Arial,14"'
            write(FN,'(a)') 'set output "'//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//'"'
            write(FN,'(a)') 'set xlabel "Time (fs)" font "arial,18"'
            write(FN,'(a)') 'set ylabel "Coupling parameter (W/(m^3 K))" font "arial,18"'
            write(FN,'(a)') 'set xrange ['//trim(adjustl(x_min))//':'//trim(adjustl(x_max))//']'
            write(FN,'(a)') 'set key top left'
            write(FN,'(a)') 'unset grid'
            write(FN,'(a)') 'plot \'
         else
            write(FN,'(a)') '@echo off & call gnuplot.exe -e "echo='//"'#'"//';set macros" "%~f0" & goto :eof'
            !write(FN,'(a)') 'gnuplot -persist -e "'
            write(FN,'(a)') 'set terminal pngcairo enhanced font ''Arial,14'';'
            write(FN,'(a)') 'set output '''//trim(adjustl(Plot_name))//'.'//trim(adjustl(numpar%fig_extention))//''';'
            write(FN,'(a)') 'set xlabel ''Time (fs)'' font "arial,18";'
            write(FN,'(a)') 'set ylabel ''Coupling parameter (W/(m^3 K))'' font "arial,18";'
            write(FN,'(a)') 'set xrange ['//trim(adjustl(x_min))//':'//trim(adjustl(x_max))//'];'
            write(FN,'(a)') 'set key top left;'
            write(FN,'(a)') 'unset grid;'
            write(FN,'(a)') 'plot \'
         endif

         ! Total
         if (sh_cmd .EQ. '.sh') then
            write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:2 with lines lw 3 lt 1 title "Total", \'
         else
            write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:2 with lines lw 3 lt 1 title ''Total'', \'
         endif

         i_start = Nat**2

         do i = 1, Nat
            do i_orb = 1, N_types
               orb_name1 = name_of_orbitals(norb, i_orb)

               do j = i, Nat
                  do j_orb = 1, N_types
                     orb_name2 = name_of_orbitals(norb, j_orb)

                     col1 = (3+i_start) + (i-1)*Nat*N_types**2 + (i_orb-1)*Nat*N_types + (j-1)*N_types + (j_orb-1)
                     col2 = (3+i_start) + (j-1)*Nat*N_types**2 + (j_orb-1)*Nat*N_types + (i-1)*N_types + (i_orb-1)

                     write(temp_col1,'(i0)') col1
                     write(temp_col2,'(i0)') col2

                     call select_linestyle_gnu(col1-i_start, linestyle) ! module "Gnuplotting"

                     if (col1 == col2) then
                        if (sh_cmd .EQ. '.sh') then
                           if ((i == Nat) .and. (j==Nat)) then
                              write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:'//trim(adjustl(temp_col1))// &
                              ' with lines lw 3 '//trim(adjustl(linestyle))//' title "'// &
                              trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                              trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//'"'
                           else
                              write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:'//trim(adjustl(temp_col1))// &
                              ' with lines lw 3 '//trim(adjustl(linestyle))//' title "'// &
                              trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                              trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//'", \'
                           endif
                        else
                           write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:'//trim(adjustl(temp_col1))// &
                              ' with lines lw 3 '//trim(adjustl(linestyle))//' title '''// &
                              trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                              trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//''', \'
                        endif
                     elseif (col1 < col2) then ! only print the sum once
                        if (sh_cmd .EQ. '.sh') then
                           write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:(column('//trim(adjustl(temp_col1))// &
                              ')+column('//trim(adjustl(temp_col2))//')) with lines lw 3 '// &
                              trim(adjustl(linestyle))//' title "'// &
                              trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                              trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//'", \'
                        else
                           write(FN,'(a)') '"'//trim(adjustl(Data_file_name))//'" using 1:(column('//trim(adjustl(temp_col1))// &
                              ')+column('//trim(adjustl(temp_col2))//')) with lines lw 3 '// &
                              trim(adjustl(linestyle))//' title '''// &
                              trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(orb_name1))//'-'// &
                              trim(adjustl(matter%Atoms(j)%Name))//' '//trim(adjustl(orb_name2))//''', \'
                        endif
                     endif

                  enddo
               enddo
            enddo
         enddo

         if (sh_cmd .EQ. '.sh') then
            write(FN,'(a)') 'EOF'
         else
            !write(FN,'(a)') '"'
         endif

         call write_gnuplot_script_ending(numpar, FN, File_name, 1)
         close(FN)
      endif

   end select
end subroutine gnu_partial_coupling





subroutine gnu_volume(numpar, File_name, file_supercell, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_supercell ! input file
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Volume','Time (fs)', 'Volume (A^3)', &
      trim(adjustl(eps_name)), numpar%path_sep, 1)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_supercell)), ' "u 1:2 w l lw LW title "Supercell volume" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"', trim(adjustl(file_supercell)), '\"u 1:2 w l lw \"$LW\" title \"Supercell volume\" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_volume


subroutine gnu_optical_coefficients(numpar, File_name, file_optics, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_optics ! optical coefficients
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Optical coefficients','Time (fs)', &
      'Optical coefficients', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_optics)), ' "u 1:2 w l lw LW title "Reflectivity" ,\'
      write(FN, '(a,a,a,a,a)') ' "', trim(adjustl(file_optics)), ' "u 1:3 w l lw LW title "Transmission" ,\'
      write(FN, '(a,a,a,a,a)') ' "', trim(adjustl(file_optics)), ' "u 1:4 w l lw LW title "Absorption" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_optics)), '\"u 1:2 w l lw \"$LW\" title \"Reflectivity\" ,\'
      write(FN, '(a,a,a,a,a)') '\"', trim(adjustl(file_optics)), '\"u 1:3 w l lw \"$LW\" title \"Transmission\" ,\'
      write(FN, '(a,a,a,a,a)') '\"', trim(adjustl(file_optics)), '\"u 1:4 w l lw \"$LW\" title \"Absorption \" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_optical_coefficients


subroutine gnu_n_and_k(numpar, File_name, file_optics, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_optics ! optical coefficients
   real(8), intent(in) :: t0, t_last	 ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   integer :: FN
   real(8) :: x_tics
   character(8) :: temp, time_order

   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics, 'Optical n and k','Time (fs)', &
      'Optical parameters', trim(adjustl(eps_name)), numpar%path_sep, 0)	! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_optics)), ' "u 1:5 w l lw LW title "n" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_optics)), ' "u 1:6 w l lw LW title "k" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_optics)), '\"u 1:5 w l lw \"$LW\" title \"n\" ,\'
      write(FN, '(a,a,a)') '\"', trim(adjustl(file_optics)), '\"u 1:6 w l lw \"$LW\" title \"k \" '
   endif
   call write_gnuplot_script_ending(numpar, FN, File_name, 1)
   close(FN)
end subroutine gnu_n_and_k





subroutine write_gnuplot_script_ending(numpar, FN, File_name, ind)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   integer, intent(in) :: FN, ind
   character(*), intent(in) :: File_name
   character(300) :: command
   integer :: iret

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      select case (ind)
      case (:1)
         write(FN, '(a)') 'reset'
         write(FN, '(a)') '" | gnuplot '
         !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
         command = 'chmod +x '//trim(adjustl(File_name))
         iret = system(command)

         !call system(trim(adjustl(File_name))) ! execute the prepared script
      case (2:)
         write(FN, '(a)') 'reset'
         write(FN, '(a)') '" | gnuplot '
         !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
         command = 'chmod +x '//trim(adjustl(File_name))
         iret = system(command)
         !call system(trim(adjustl(File_name))) ! execute the prepared script
      endselect
   endif
end subroutine write_gnuplot_script_ending


subroutine execute_all_gnuplots(numpar, file_path)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   character(*), intent(in) :: file_path
   character(300) :: command
   integer :: iret

   !call chdir(trim(adjustl(file_path)))
   command = trim(adjustl(file_path))
   iret = chdir(command)

   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      !call system("OUTPUT_Gnuplot_all.cmd")
      !command = "OUTPUT_Gnuplot_all.cmd"
      command = trim(adjustl(m_Gnuplot_all))//".cmd"
      iret = system(command)
   else ! linux:
      !call system("./OUTPUT_Gnuplot_all.sh")
      !command = "./OUTPUT_Gnuplot_all.sh"
      command = "./"//trim(adjustl(m_Gnuplot_all))//".sh"
      iret = system(command)
   endif
   !call chdir("../")
   iret = chdir("../")
end subroutine execute_all_gnuplots


subroutine execute_gnuplot(File_name)
   character(*), intent(in) :: File_name
   character(100) :: command
   integer :: iret
   !call system(trim(adjustl(File_name))) ! execute the prepared script
   command = "./OUTPUT_Gnuplot_all.sh"
   iret = system(command)
end subroutine execute_gnuplot


subroutine write_energy_levels_gnuplot(numpar, FN, t0, Scell, file_Ei)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   integer, intent(in) :: FN            ! file to write into
   real(8), intent(in) :: t0  ! starting time
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   character(*), intent(in) :: file_Ei  ! file with energy levels
   integer i, M, NSC
   character(30) :: ch_temp, ch_temp2

   write(ch_temp2,'(f)') t0

   do NSC = 1, size(Scell)
      M = size(Scell(NSC)%Ei)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') 25.0d0      ! Scell(NSC)%E_top

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,a,a,i5,a)') 'p ['//trim(adjustl(ch_temp2))//':][:'//trim(adjustl(ch_temp))// &
                                      '] "', trim(adjustl(file_Ei)), '"u 1:', 2, ' pt 7 ps 0.2 ,\'
         do i = 3, M
            write(FN, '(a,a,a,i5,a)') ' "', trim(adjustl(file_Ei)), '"u 1:', i, ' pt 7 ps 0.2 ,\'
         enddo
         write(FN, '(a,a,a,i5,a)') ' "', trim(adjustl(file_Ei)), '"u 1:', M+1, ' pt 7 ps 0.2'
      else
         write(FN, '(a,a,a,i5,a)') 'p ['//trim(adjustl(ch_temp2))//':][:'//trim(adjustl(ch_temp))// &
                                      '] \"', trim(adjustl(file_Ei)), &
                                       '\"u 1:', 2, ' pt 7 ps \"$LW\" ,\'
         do i = 3, M
            write(FN, '(a,a,a,i5,a)') '\"', trim(adjustl(file_Ei)), '\"u 1:', i, ' pt 7 ps \"$LW\" ,\'
         enddo
         write(FN, '(a,a,a,i5,a)') '\"', trim(adjustl(file_Ei)), '\"u 1:', M+1, ' pt 7 ps \"$LW\"'
      endif
   enddo
end subroutine write_energy_levels_gnuplot




subroutine write_diffraction_powder_gnuplot(FN, Scell, matter, numpar, file_powder, min_x, max_x)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter            ! material parameters
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_powder  ! file with powder diffraction spectrum
   real(8), intent(in), optional :: min_x, max_x      ! start and end of x-grid
   !-----------------------
   real(8) :: x_start, x_end
   integer :: i, M, NSC, ind, j, k
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4, ch_i
   logical :: do_fe_eq

   if (present(min_x)) then
      x_start = min_x
   else ! default
      x_start = 10.0d0
   endif

   if (present(max_x)) then
      x_end = max_x
   else ! default
      x_end = 180.0d0
   endif

   do NSC = 1, size(Scell)
      write(ch_temp,'(f)') x_end
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f)') numpar%dt_save
      ! grid start:
      write(ch_temp4,'(f)') x_start

      if (numpar%path_sep .EQ. '\') then      ! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_powder))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         ! Add partial ones:
         if (matter%N_KAO > 1) then ! only makes sence if there is more then one element in the compound:
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][:] "'//trim(adjustl(file_powder))// &
                  '" index (i-1) u 1:2 w l lw 2 lt rgb "black" title sprintf("%i fs total",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) ,\'

            ind = 2     ! top start counting columns
            do j = 1, matter%N_KAO     ! for all elements
               do k = j, matter%N_KAO  ! for all different pairs
                  ind = ind + 1
                  write(ch_i,'(i0)') ind
                  if ( (j==matter%N_KAO) .and. (k==matter%N_KAO) ) then ! last line
                     write(FN, '(a)') ' "'//trim(adjustl(file_powder))// &
                        '" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title "'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'" '
                  else
                     write(FN, '(a)') ' "'//trim(adjustl(file_powder))// &
                        '" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title "'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'" ,\'
                  endif
               enddo ! k
            enddo ! j
         else ! there is only one column
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][:] "'//trim(adjustl(file_powder))// &
                  '" index (i-1) u 1:2 w l lw 2 lt rgb "black" title sprintf("%i fs",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif

      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_powder))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         if (matter%N_KAO > 1) then ! if there are more then 1 column
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][:] \"'//trim(adjustl(file_powder))// &
                  '\" index (i-1) u 1:2 w l lw 2 lt rgb \"black\" title sprintf(\"%i fs total\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) ,\'
            ind = 2     ! top start counting columns
            do j = 1, matter%N_KAO     ! for all elements
               do k = j, matter%N_KAO  ! for all different pairs
                  ind = ind + 1
                  write(ch_i,'(i0)') ind
                  if ( (j==matter%N_KAO) .and. (k==matter%N_KAO) ) then ! last line
                     write(FN, '(a)') ' \"'//trim(adjustl(file_powder))// &
                        '\" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title \"'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'\" '
                  else
                     write(FN, '(a)') ' \"'//trim(adjustl(file_powder))// &
                        '\" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title \"'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'\" ,\'
                  endif
               enddo ! k
            enddo ! j
         else ! there is only one column
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][:] \"'//trim(adjustl(file_powder))// &
                  '\" index (i-1) u 1:2 w l lw 2 lt rgb \"black\" title sprintf(\"%i fs\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_diffraction_powder_gnuplot



subroutine write_atomic_distribution_gnuplot(FN, Scell, numpar, file_fe, its_pot, no_maxwell)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_fe  ! file with atomic distribution function
   logical, optional :: its_pot, no_maxwell
   !-----------------------
   integer :: i, M, NSC
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4, ch_titel
   logical :: do_fe_eq, no_maxw, poten

   if (present(no_maxwell)) then
      no_maxw = no_maxwell
   else  ! assume there is maxwellian distribution
      no_maxw = .false.
   endif


   if (present(its_pot)) then
      poten = its_pot
   else
      poten = .false.
   endif


   do NSC = 1, size(Scell)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)')  Scell(NSC)%Ea_pot_grid_out(size(Scell(NSC)%Ea_pot_grid_out))
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f)') numpar%dt_save

      select case (numpar%el_ion_scheme)
         case (3:5)
            do_fe_eq = .true.
         case default
            do_fe_eq = .false.
      endselect
      ! minimal energy grid:
      write(ch_temp4,'(f)') 0.0d0
      if (poten) then
         write(ch_temp4,'(f)') Scell(NSC)%Ea_pot_grid_out(1)
         ch_titel = 'Equivalent Gibbs'
      else
         ch_titel = 'Equivalent Maxwell'
      endif

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         !write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:5] "'//trim(adjustl(file_fe))// &

         if (no_maxw) then ! just distribution, without equivalent Maxwell:
            write(FN, '(a)') 'p [:][0:5] "'//trim(adjustl(file_fe))// &
                  '" index (i) u 1:2 pt 7 ps 1 title sprintf("%i fs",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         else ! with equivalent Maxwell
            write(FN, '(a)') 'p [:][0:5] "'//trim(adjustl(file_fe))// &
                  '" index (i) u 1:3 w l lw 2 lt rgb "grey" title "'// trim(adjustl(ch_titel)) //'" ,\'

            write(FN, '(a)') ' "'//trim(adjustl(file_fe))// &
                  '" index (i) u 1:2 pt 7 ps 1 title sprintf("%i fs",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif ! no_maxw
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         !write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:5] \"'//trim(adjustl(file_fe))// &
         if (no_maxw) then ! just distribution, without equivalent Maxwell:
            write(FN, '(a)') 'p [:][0:5] \"'//trim(adjustl(file_fe))// &
            '\" index (i) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         else
            write(FN, '(a)') 'p [:][0:5] \"'//trim(adjustl(file_fe))// &
                  '\" index (i) u 1:3 w l lw 2 lt rgb \"grey\" title \"'// trim(adjustl(ch_titel)) //'\" ,\'

            write(FN, '(a)') ' \"'//trim(adjustl(file_fe))// &
                  '\" index (i) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif ! no_maxw
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_atomic_distribution_gnuplot



subroutine write_pair_correlation_gnuplot(FN, Scell, numpar, matter, file_PCF, min_x, max_x)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   type(solid), intent(in) :: matter	! materil parameters
   character(*), intent(in) :: file_PCF  ! file with atomic pair correlation function
   real(8), intent(in), optional :: min_x, max_x      ! start and end of x-grid
   !-----------------------------------------
   real(8) :: x_start, x_end
   integer :: i, M, NSC, j, k, ind
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4, ch_i
   logical :: do_fe_eq

   if (present(min_x)) then
      x_start = min_x
   else ! default
      x_start = 1.0d0
   endif

   if (present(max_x)) then
      x_end = max_x
   else ! default
      x_end = 1+ceiling( 0.5d0* abs(min (Scell(1)%supce(1,1), Scell(1)%supce(2,2), Scell(1)%supce(3,3) ) ) ) ! half of the supercell
   endif

   do NSC = 1, size(Scell)
      ! Choose the maximal distance [A]:
      write(ch_temp,'(f)') x_end
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f)') numpar%dt_save

      ! minimal energy grid:
      write(ch_temp4,'(f)') x_start

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_PCF))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         ! Add partial ones:
         if (matter%N_KAO > 1) then ! only makes sence if there is more then one element in the compound:
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:10] "'//trim(adjustl(file_PCF))// &
                  '" index (i-1) u 1:2 w l lw 2 lt rgb "black" title sprintf("%i fs total",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) ,\'

            ind = 2     ! top start counting columns
            do j = 1, matter%N_KAO     ! for all elements
               do k = j, matter%N_KAO  ! for all different pairs
                  ind = ind + 1
                  write(ch_i,'(i0)') ind
                  if ( (j==matter%N_KAO) .and. (k==matter%N_KAO) ) then ! last line
                     write(FN, '(a)') ' "'//trim(adjustl(file_PCF))// &
                        '" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title "'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'" '
                  else
                     write(FN, '(a)') ' "'//trim(adjustl(file_PCF))// &
                        '" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title "'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'" ,\'
                  endif
               enddo ! k
            enddo ! j
         else ! there is only one column
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:10] "'//trim(adjustl(file_PCF))// &
                  '" index (i-1) u 1:2 w l lw 2 lt rgb "black" title sprintf("%i fs total",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '

         endif
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_PCF))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         if (matter%N_KAO > 1) then ! if there are more then 1 column
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:10] \"'//trim(adjustl(file_PCF))// &
                  '\" index (i-1) u 1:2 w l lw 2 lt rgb \"black\" title sprintf(\"%i fs total\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) ,\'
            ind = 2     ! top start counting columns
            do j = 1, matter%N_KAO     ! for all elements
               do k = j, matter%N_KAO  ! for all different pairs
                  ind = ind + 1
                  write(ch_i,'(i0)') ind
                  if ( (j==matter%N_KAO) .and. (k==matter%N_KAO) ) then ! last line
                     write(FN, '(a)') ' \"'//trim(adjustl(file_PCF))// &
                        '\" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title \"'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'\" '
                  else
                     write(FN, '(a)') ' \"'//trim(adjustl(file_PCF))// &
                        '\" index (i-1) u 1:'//trim(adjustl(ch_i))//' w l lw 2 title \"'// &
                        trim(adjustl(matter%Atoms(j)%Name))//'-'//trim(adjustl(matter%Atoms(k)%Name))//'\" ,\'
                  endif
               enddo ! k
            enddo ! j
         else ! there is only one column
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:10] \"'//trim(adjustl(file_PCF))// &
                  '\" index (i-1) u 1:2 w l lw 2 lt rgb \"black\" title sprintf(\"%i fs total\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif
      endif
      write(FN, '(a)') '}'
   enddo ! NSC
end subroutine write_pair_correlation_gnuplot


subroutine write_distribution_gnuplot(FN, Scell, numpar, file_fe, min_x, max_x)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_fe  ! file with electronic distribution function
   real(8), intent(in), optional :: min_x, max_x      ! start and end of x-grid
   !-----------------------
   real(8) :: x_start, x_end
   integer :: i, M, NSC
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4
   logical :: do_fe_eq

   if (present(min_x)) then
      x_start = min_x
   else ! default
      x_start = -25.0d0
   endif

   if (present(max_x)) then
      x_end = max_x
   else ! default
      x_end = 25.0d0
   endif


   do NSC = 1, size(Scell)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') x_end      ! Scell(NSC)%E_top
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f)') numpar%dt_save

      select case (numpar%el_ion_scheme)
         case (3:4)
            do_fe_eq = .true.
         case default
            do_fe_eq = .false.
      endselect
      ! minimal energy grid:
      write(ch_temp4,'(f)') x_start  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         !if (do_fe_eq) then  ! plot also equivalent Fermi distribution

            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:2] "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:3 w l lw 2 lt rgb "grey" title "Equivalent Fermi" ,\'

            if (numpar%do_partial_thermal) then ! add band-resolved equivalent distributions
               write(FN, '(a)') ' "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:4 w l lw 2 lt rgb "blue" title "Equivalent Fermi in VB" ,\'
               write(FN, '(a)') ' "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:5 w l lw 2 lt rgb "red" title "Equivalent Fermi in CB" ,\'
            endif

            write(FN, '(a)') ' "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         !else
         !   write(FN, '(a)') 'p [:'//trim(adjustl(ch_temp))//'][0:2] "'//trim(adjustl(file_fe))// &
         !         '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i-1+'// &
         !         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
         !endif
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         !if (do_fe_eq) then ! plot also equivalent Fermi distribution
            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:2] \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:3 w l lw 2 lt rgb \"grey\" title \"Equivalent Fermi\" ,\'

            if (numpar%do_partial_thermal) then ! add band-resolved equivalent distributions
               write(FN, '(a)') ' \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:4 w l lw 2 lt rgb \"blue\" title \"Equivalent Fermi in VB\" ,\'
               write(FN, '(a)') ' \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:5 w l lw 2 lt rgb \"red\" title \"Equivalent Fermi in CB\" ,\'
            endif

            write(FN, '(a)') ' \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",((i-1)*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         !else
         !   write(FN, '(a)') 'p [:'//trim(adjustl(ch_temp))//'][0:2] \"'//trim(adjustl(file_fe))// &
         !         '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i-1+'// &
         !         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
         !endif
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_distribution_gnuplot


subroutine write_orb_distribution_gnuplot(FN, Scell, numpar, matter, file_fe)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   type(Solid), intent(in) :: matter	! Material parameters
   character(*), intent(in) :: file_fe  ! file with electronic distribution function
   !-----------------------
   integer :: N_at, N_types, Nsiz, i_at, i_types, nat, norb
   integer :: i, M, NSC, col, i_col, NKOA
   character(2) :: chtemp1
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4, ch_name, ch_col
   logical :: do_fe_eq

   do NSC = 1, size(Scell)

      NKOA = matter%N_KAO    ! number of kinds of atoms
      ! Find number of orbitals per atom:
      nat = size(Scell(NSC)%MDatoms) ! number of atoms
      Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
      norb =  Nsiz/nat ! orbitals per atom
      ! Find number of different orbital types:
      N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"


      ! minimal energy grid:
      write(ch_temp4,'(f)') -20.0d0  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') 25.0d0      ! Scell(NSC)%E_top
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f10.4)') numpar%dt_save

      if (numpar%do_partial_thermal) then ! includes band-resolved equivalent distributions
         col = 5  ! column number after which orbital-resolved data start
      else
         select case (numpar%el_ion_scheme)
            case (3:5)
               col = 3  ! column number after which orbital-resolved data start
            case default
               col = 2  ! column number after which orbital-resolved data start
         endselect
      endif

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         i_col = col ! to start with
         do i_at = 1, NKOA
            do i_types = 1, N_types
               chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
               write(ch_name,'(a)') trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
               i_col = i_col + 1 ! column number
               write(ch_col, '(i4)') i_col

               if ( (i_at == 1) .and. (i_types == 1) ) then ! first column
                  write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:2] "'//trim(adjustl(file_fe))// &
                      '" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' pt 7 ps 1 title sprintf("%i fs ' // trim(adjustl(ch_name)) // &
                      '",((i-1)*' // trim(adjustl(ch_temp3)) //'+'// trim(adjustl(ch_temp2))// ')) ,\'
               elseif ( (i_at == NKOA) .and. (i_types == N_types) ) then ! last column
                  write(FN, '(a)') '"'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' pt 7 ps 1 title "'// trim(adjustl(ch_name))//'"'
               else  ! normal column
                  write(FN, '(a)') '"'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' pt 7 ps 1 title "'// trim(adjustl(ch_name))//'" ,\'
               endif
            enddo ! i_types
         enddo ! i_at

      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         i_col = col ! to start with
         do i_at = 1, NKOA
            do i_types = 1, N_types
               chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
               write(ch_name,'(a)') trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
               i_col = i_col + 1 ! column number
               write(ch_col, '(i4)') i_col

               if ( (i_at == 1) .and. (i_types == 1) ) then ! first column
                  write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:2] \"'//trim(adjustl(file_fe))// &
                      '\" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' pt 7 ps 1 title sprintf(\"\%i fs ' // trim(adjustl(ch_name)) // &
                      '\",((i-1)*'// trim(adjustl(ch_temp3)) //'+'// trim(adjustl(ch_temp2))// ')) ,\'
               elseif ( (i_at == NKOA) .and. (i_types == N_types) ) then ! last column
                  write(FN, '(a)') '\"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' pt 7 ps 1 title \"'// trim(adjustl(ch_name))//'\"'
               else  ! normal column
                  write(FN, '(a)') '\"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' pt 7 ps 1 title \"'// trim(adjustl(ch_name))//'\" ,\'
               endif
            enddo ! i_types
         enddo ! i_at

      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_orb_distribution_gnuplot


subroutine write_distribution_on_grid_gnuplot(FN, Scell, numpar, file_fe)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_fe  ! file with electronic distribution function
   !-----------------------
   real(8) :: E_max
   integer :: i, M, NSC
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4

   do NSC = 1, size(Scell)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      E_max = Scell(1)%E_fe_grid(size(Scell(1)%E_fe_grid))
      write(ch_temp,'(f16.2)') E_max
      write(ch_temp2,'(f16.2)') numpar%t_start
      if (numpar%t_start > 0.0d0) then ! positive, add plus
         ch_temp2 = '+'//trim(adjustl(ch_temp2))
      endif
      write(ch_temp3,'(f13.6)') numpar%dt_save

      ! minimal energy grid:
      write(ch_temp4,'(f)') -25.0d0  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'set logscale y'
         write(FN, '(a)') 'set format y "10^{%L}"'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:200] "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs total",(i-1)'// &
                  '*' // trim(adjustl(ch_temp3)) // trim(adjustl(ch_temp2))// ') ,\'

         write(FN, '(a)') '"" index (i-1) u 1:4 pt 6 ps 1 title "Photoelectrons" ,\'
         write(FN, '(a)') '"" index (i-1) u 1:5 pt 5 ps 1 title "Impact electrons" ,\'
         write(FN, '(a)') '"" index (i-1) u 1:6 pt 4 ps 1 title "Auger electrons"'
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'set logscale y'
         write(FN, '(a)') 'set format y \"10^{%L}\"'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:200] \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i-1)'// &
                  '*' // trim(adjustl(ch_temp3)) // trim(adjustl(ch_temp2))// ') ,\'
         write(FN, '(a)') '\"\" index (i-1) u 1:4 pt 6 ps 1 title \"Photoelectrons\" ,\'
         write(FN, '(a)') '\"\" index (i-1) u 1:5 pt 5 ps 1 title \"Impact electrons\" ,\'
         write(FN, '(a)') '\"\" index (i-1) u 1:6 pt 4 ps 1 title \"Auger electrons\"'
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_distribution_on_grid_gnuplot




subroutine write_DOS_gnuplot(FN, Scell, numpar, matter, file_fe)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   type(Solid), intent(in) :: matter	! Material parameters
   character(*), intent(in) :: file_fe  ! file with electronic distribution function
   !-----------------------
   integer :: N_at, N_types, Nsiz, i_at, i_types, nat, norb
   integer :: i, M, NSC, col, i_col, NKOA
   character(2) :: chtemp1
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4, ch_name, ch_col
   logical :: do_fe_eq

   do NSC = 1, size(Scell)

      NKOA = matter%N_KAO    ! number of kinds of atoms
      ! Find number of orbitals per atom:
      nat = size(Scell(NSC)%MDatoms) ! number of atoms
      Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
      norb =  Nsiz/nat ! orbitals per atom
      ! Find number of different orbital types:
      N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"


      ! minimal energy grid:
      write(ch_temp4,'(f)') -20.0d0  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') 25.0d0      ! Scell(NSC)%E_top
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f10.4)') numpar%dt_save

      col = 2  ! column number after which orbital-resolved data start

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:] "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:2 w l lw 2 lt rgb "black" title sprintf("%i fs Total",(i-1)'// &
                  '*' // trim(adjustl(ch_temp3)) // '+' // trim(adjustl(ch_temp2))// ') ,\'

         i_col = col ! to start with
         do i_at = 1, NKOA
            do i_types = 1, N_types
               chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
               write(ch_name,'(a)') trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
               i_col = i_col + 1 ! column number
               write(ch_col, '(i4)') i_col

               if ( (i_at == NKOA) .and. (i_types == N_types) ) then ! last column
                  write(FN, '(a)') '"'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' w l lw 2 title "'// trim(adjustl(ch_name))//'"'
               else  ! normal column
                  write(FN, '(a)') '"'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' w l lw 2 title "'// trim(adjustl(ch_name))//'" ,\'
               endif
            enddo ! i_types
         enddo ! i_at

      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:] \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:2 w l lw 2 lt rgb \"black\" title sprintf(\"%i fs Total\",(i-1)'// &
                  '*' // trim(adjustl(ch_temp3)) // '+' // trim(adjustl(ch_temp2))// ') ,\'

         i_col = col ! to start with
         do i_at = 1, NKOA
            do i_types = 1, N_types
               chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
               write(ch_name,'(a)') trim(adjustl(matter%Atoms(i_at)%Name))//' '//trim(adjustl(chtemp1))
               i_col = i_col + 1 ! column number
               write(ch_col, '(i4)') i_col

               if ( (i_at == NKOA) .and. (i_types == N_types) ) then ! last column
                  write(FN, '(a)') '\"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' w l lw 2 title \"'// trim(adjustl(ch_name))//'\"'
               else  ! normal column
                  write(FN, '(a)') '\"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:' // trim(adjustl(ch_col)) // ' w l lw 2 title \"'// trim(adjustl(ch_name))//'\" ,\'
               endif
            enddo ! i_types
         enddo ! i_at

      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_DOS_gnuplot




subroutine Plot_electron_MFP_gunplot(numpar, matter, file_electron_IMFP, file_electron_EMFP)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   character(*), intent(in) :: file_electron_IMFP, file_electron_EMFP     ! file with data
   !-------------
   character(300) :: gnu_electron_MFP
   character(11) :: sh_cmd, call_slash
   character(8) :: col
   integer :: N_grid, FN, count_col, i, Nshl, j_start, j
   real(8) :: t0, t_last, x_tics
   !-------------

   ! Define if sh or cmd script required:
   call define_sh_vs_cmd(numpar%path_sep, call_slash, sh_cmd)      ! below


   N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
   t0 = 1.0d0  ! [eV] starting energy point
   t_last = matter%Atoms(1)%Ph_MFP(1)%E(N_grid) ! ending energy point

   ! 1) Electron MFPs gnuplotting:
   gnu_electron_MFP = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//'OUTPUT_MFP_electron'//trim(adjustl(sh_cmd))
   open(NEWUNIT=FN, FILE = trim(adjustl(gnu_electron_MFP)), action="write", status="replace")

   x_tics = 10.0d0
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics,  'MFPs', &
      'Electron energy (eV)', 'Mean free path (A)', 'OUTPUT_MFP_electron.'//trim(adjustl(numpar%fig_extention)), &
      numpar%path_sep, setkey=1, set_x_log=.true., set_y_log=.true.)  ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then  ! if it is Windows
      count_col = 2  ! to start with
      write(FN, '(a,es15.6,a,es15.6,a,a,a)') 'p [', t0, ':', t_last, '][1:1e5] "' , trim(adjustl(file_electron_IMFP)), &
         ' "u 1:2 w l lw LW title "'//trim(adjustl(matter%Atoms(1)%Name))//' '//trim(adjustl(matter%Atoms(1)%Shell_name(1))) &
         //'" ,\'
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         if (i == 1) then
            j_start = 2
         else
            j_start = 1
         endif

         do j = j_start, Nshl    ! for all shells of this atom
            if ((i == 1) .and. (j == Nshl)) then ! valence
               ! skip here, plot later
            else ! core shell
               count_col = count_col + 1  ! number of columns
               write(col, '(i4)') count_col
               write(FN, '(a)') '"'//trim(adjustl(file_electron_IMFP))// &
                  '" u 1:'//trim(adjustl(col))//' w l lw LW title "' &
                  //trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j))) &
                  //'" ,\'
            endif
         enddo
      enddo
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Valence
      write(FN, '(a)') '"'//trim(adjustl(file_electron_IMFP))// &
                  '" u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Total
      write(FN, '(a)') '"'//trim(adjustl(file_electron_IMFP))// &
                     '" u 1:'//trim(adjustl(col))//' w l lw LW title "Total inelastic" ,\'
      ! Elastic MFP:
      write(col, '(i4)') 1+size(matter%Atoms)+1
      write(FN, '(a)') '"'//trim(adjustl(file_electron_EMFP))// &
                     '" u 1:'//trim(adjustl(col))//' w l lw LW title "Elastic" '
   else  ! Linux
      count_col = 2  ! to start with
      write(FN, '(a,es15.6,a,es15.6,a,a,a)') 'p [', t0, ':', t_last, '][1:1e5] \"' , trim(adjustl(file_electron_IMFP)), &
         '\" u 1:2 w l lw \"$LW\" title \"'//trim(adjustl(matter%Atoms(1)%Name))//' '//trim(adjustl(matter%Atoms(1)%Shell_name(1))) &
         //'\" ,\'
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         if (i == 1) then
            j_start = 2
         else
            j_start = 1
         endif

         do j = j_start, Nshl    ! for all shells of this atom
            if ((i == 1) .and. (j == Nshl)) then ! valence
               ! skip here, plot later
            else ! core shell
               count_col = count_col + 1  ! number of columns
               write(col, '(i4)') count_col
               write(FN, '(a)') '\"'//trim(adjustl(file_electron_IMFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' &
                  //trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j))) &
                  //'\" ,\'
            endif
         enddo
      enddo
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Valence
      write(FN, '(a)') '\"'//trim(adjustl(file_electron_IMFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Total
      write(FN, '(a)') '\"'//trim(adjustl(file_electron_IMFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total inelastic\" ,\'
      ! Elastic MFP:
      write(col, '(i4)') 1+size(matter%Atoms)+1
      write(FN, '(a)') '\"'//trim(adjustl(file_electron_EMFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Elastic\" '
   endif
   !call write_gnuplot_script_ending(numpar, FN, gnu_electron_MFP, 1)   ! below
   call write_gnuplot_script_ending_new(FN, gnu_electron_MFP, numpar%path_sep) ! module "Gnuplotting"
   close(FN)
end subroutine Plot_electron_MFP_gunplot



subroutine Plot_photon_MFP_gunplot(numpar, matter, file_photon_MFP)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   character(*), intent(in) :: file_photon_MFP     ! file with data
   !-------------
   character(300) :: gnu_photon_MFP
   character(11) :: sh_cmd, call_slash
   character(8) :: col
   integer :: N_grid, FN, count_col, i, Nshl, j_start, j
   real(8) :: t0, t_last, x_tics
   !-------------

   ! Define if sh or cmd script required:
   call define_sh_vs_cmd(numpar%path_sep, call_slash, sh_cmd)      ! below

   N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
   t0 = 1.0d0  ! [eV] starting energy point
   t_last = matter%Atoms(1)%Ph_MFP(1)%E(N_grid) ! ending energy point

   gnu_photon_MFP = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//'OUTPUT_MFP_photon'//trim(adjustl(sh_cmd))
   open(NEWUNIT=FN, FILE = trim(adjustl(gnu_photon_MFP)), action="write", status="replace")

   x_tics = 10.0d0
   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics,  'MFPs', &
      'Photon energy (eV)', 'Attenuation length (A)', 'OUTPUT_MFP_photon.'//trim(adjustl(numpar%fig_extention)), &
      numpar%path_sep, setkey=1, set_x_log=.true., set_y_log=.true.)  ! module "Gnuplotting"

   if (numpar%path_sep .EQ. '\') then  ! if it is Windows
      count_col = 2  ! to start with
      write(FN, '(a,es15.6,a,es15.6,a,a,a)') 'p [',t0, ':', t_last, '][10:1e7] "' , trim(adjustl(file_photon_MFP)), &
         '" u 1:2 w l lw LW title "'//trim(adjustl(matter%Atoms(1)%Name))//' '//trim(adjustl(matter%Atoms(1)%Shell_name(1))) &
         //'" ,\'
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         if (i == 1) then
            j_start = 2
         else
            j_start = 1
         endif

         do j = j_start, Nshl    ! for all shells of this atom
            if ((i == 1) .and. (j == Nshl)) then ! valence
               ! skip here, plot later
            else ! core shell
               count_col = count_col + 1  ! number of columns
               write(col, '(i4)') count_col
               write(FN, '(a)') '"'//trim(adjustl(file_photon_MFP))// &
                  '" u 1:'//trim(adjustl(col))//' w l lw LW title "' &
                  //trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j))) &
                  //'" ,\'
            endif
         enddo
      enddo
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Valence
      write(FN, '(a)') '"'//trim(adjustl(file_photon_MFP))// &
                  '" u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Total
      write(FN, '(a)') '"'//trim(adjustl(file_photon_MFP))// &
                     '" u 1:'//trim(adjustl(col))//' w l lw LW title "Total" '
   else  ! Linux
      count_col = 2  ! to start with
      write(FN, '(a,es15.6,a,es15.6,a,a,a)') 'p [', t0, ':', t_last, '][10:1e7] \"' , trim(adjustl(file_photon_MFP)), &
         '\" u 1:2 w l lw \"$LW\" title \"'//trim(adjustl(matter%Atoms(1)%Name))//' '//trim(adjustl(matter%Atoms(1)%Shell_name(1))) &
         //'\" ,\'
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         if (i == 1) then
            j_start = 2
         else
            j_start = 1
         endif

         do j = j_start, Nshl    ! for all shells of this atom
            if ((i == 1) .and. (j == Nshl)) then ! valence
               ! skip here, plot later
            else ! core shell
               count_col = count_col + 1  ! number of columns
               write(col, '(i4)') count_col
               write(FN, '(a)') '\"'//trim(adjustl(file_photon_MFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' &
                  //trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j))) &
                  //'\" ,\'
            endif
         enddo
      enddo
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Valence
      write(FN, '(a)') '\"'//trim(adjustl(file_photon_MFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'
      count_col = count_col + 1  ! number of columns
      write(col, '(i4)') count_col  ! Total
      write(FN, '(a)') '\"'//trim(adjustl(file_photon_MFP))// &
                  '\" u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total\"'

   endif
   !call write_gnuplot_script_ending(numpar, FN, gnu_photon_MFP, 1)  ! below
   call write_gnuplot_script_ending_new(FN, gnu_photon_MFP, numpar%path_sep) ! module "Gnuplotting"
   close(FN)
end subroutine Plot_photon_MFP_gunplot



subroutine Plot_laser_spectrum_gnuplot(numpar, laser, file_spectrum, i_pulse)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   character(*), intent(in) :: file_spectrum     ! file with data
   integer, intent(in) :: i_pulse   ! index of the pulse
   !-------------
   character(300) :: gnu_photon_spectrum
   character(11) :: sh_cmd, call_slash
   character(8) :: col
   character(18) :: temp, ch_temp
   integer :: N_grid, FN, count_col, i, Nshl, j_start, j, N_pulse
   real(8) :: t0, t_last, x_tics
   !-------------

   ! Define if sh or cmd script required:
   call define_sh_vs_cmd(numpar%path_sep, call_slash, sh_cmd)      ! below

   N_pulse = size(laser)      ! number if laser pulses
   if (N_pulse > 1) then
      write(temp,'(i0)') i_pulse ! number of the pulse
      temp = '_pulse_'//trim(adjustl(temp))
   else
      temp = ''   ! empty
   endif

   ! Grid size:
   N_grid = size(laser(i_pulse)%Spectrum,2)
   t0 = laser(i_pulse)%Spectrum(1,1)
   t_last = laser(i_pulse)%Spectrum(1,N_grid)


   ! Photon spectrum gnuplotting:
   gnu_photon_spectrum = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))// &
            'OUTPUT_photon_spectrum'//trim(adjustl(temp))//trim(adjustl(sh_cmd))
   open(NEWUNIT=FN, FILE = trim(adjustl(gnu_photon_spectrum)), action="write", status="replace")

   ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), ch_temp, x_tics=x_tics)	! module "Little_subroutines"

   call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics,  'Spectrum', &
      'Photon energy (eV)', 'Photon spectrum (arb. units)', 'OUTPUT_photon_spectrum.'//trim(adjustl(numpar%fig_extention)), &
      numpar%path_sep, setkey=0, set_x_log=.false., set_y_log=.false.)  ! module "Gnuplotting"


   if (numpar%path_sep .EQ. '\') then  ! if it is Windows
      write(FN, '(a,es15.6,a,es15.6,a,a,a)') 'p [', t0, ':', t_last, '][0:] "' , trim(adjustl(file_spectrum)), &
            ' "u 1:2 w l lw LW dashtype "_" title "Incoming" ,\'

      write(FN, '(a)') '"'//trim(adjustl(file_spectrum))// &
                  '" u 1:3 w l lw LW title "Absorbed" ,\'

      write(FN, '(a)') '"'//trim(adjustl(file_spectrum))// &
                  '" u 1:4 w l lt rgb "#000000" lw LW dashtype "_." title "MC Sampled" '
   else  ! Linux
      write(FN, '(a,es15.6,a,es15.6,a,a,a)') 'p [', t0, ':', t_last, '][0:] \"' , trim(adjustl(file_spectrum)), &
            '\" u 1:2 w l lw \"$LW\" dashtype \"_\" title \"Incomming\" ,\'

      write(FN, '(a)') '\"'//trim(adjustl(file_spectrum))// &
                  '\" u 1:3 w l lw \"$LW\" title \"Absorbed\" ,\'

      write(FN, '(a)') '\"'//trim(adjustl(file_spectrum))// &
                  '\" u 1:4 w l lt rgb \"#000000\" lw \"$LW\" dashtype \"_.\" title \"MC Sampled\" '
   endif

   call write_gnuplot_script_ending_new(FN, gnu_photon_spectrum, numpar%path_sep) ! module "Gnuplotting"
   close(FN)
end subroutine Plot_laser_spectrum_gnuplot


subroutine define_sh_vs_cmd(path_sep, call_slash, sh_cmd)
   character(*), intent(in) :: path_sep
   character(*), intent(out) :: call_slash, sh_cmd

   if (path_sep .EQ. '\') then	! if it is Windows
      call_slash = 'call '
      sh_cmd = '.cmd'
   else ! It is linux
      call_slash = './'
      sh_cmd = '.sh'
   endif
end subroutine define_sh_vs_cmd


end module Plots_gnuplot
