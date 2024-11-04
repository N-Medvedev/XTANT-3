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
! This file contains all use-statements for the main file "XTANT_MAIN_FILE.f90"

use Universal_constants
use Objects
use Variables
use Algebra_tools, only : DET_3X3
use Little_subroutines, only : print_time, parse_yes_no, set_starting_time, print_time_step, parse_time
use Dealing_with_files, only : close_file
use Read_input_data, only : get_add_data, Read_Input_Files, m_starline, m_INPUT_MATERIAL, m_NUMERICAL_PARAMETERS, &
        m_INPUT_directory, m_INPUT_ALL, check_all_warnings
use Dealing_with_output_files, only : reset_dt, Print_title, prepare_output_files, write_output_files, convolve_output, &
        communicate, close_save_files, close_output_files, save_duration, execute_all_gnuplots, write_energies, &
        XTANT_label, m_Error_log_file, printout_CDF_file, print_a_comforting_message, printout_MFP_file
use Gnuplotting, only : collect_gnuplots
use Initial_configuration, only : set_initial_configuration
use Atomic_thermodynamics, only : get_atomic_distribution, update_Ta_config_running_average
use Atomic_tools, only : get_mean_square_displacement, save_last_timestep, Cooling_atoms, &
        Coordinates_rel_to_abs, velocities_abs_to_rel, shortest_distance, &
        make_time_step_atoms, make_time_step_supercell, numerical_acceleration, Get_testmode_add_data, &
        get_diffraction_peaks
use Electron_tools, only : Electron_thermalization, get_glob_energy, get_low_energy_distribution, update_fe
use Transport, only : Electron_transport, Atomic_heat_transport, Change_affected_layer
use TB, only : get_new_energies, get_DOS, Get_pressure, get_electronic_thermal_parameters, &
        vdW_interplane, Electron_ion_coupling, update_nrg_after_change,  get_Hamilonian_and_E, MD_step, get_Mullikens_all, &
        Get_configurational_temperature
use Optical_parameters, only : get_optical_parameters
use MC_cross_sections, only : get_mfps, get_photon_attenuation
use Monte_carlo, only : MC_Propagate, process_laser_parameters
use ZBL_potential, only : get_total_ZBL
use TB_complex, only : use_complex_Hamiltonian

use MPI_subroutines, only : initialize_MPI, Initialize_ScaLAPACK, initialize_random_seed, MPI_share_add_data, MPI_error_wrapper, MPI_barrier_wrapper, &
        MPI_share_Read_Input_Files, MPI_share_initial_configuration, MPI_share_electron_MFPs, MPI_share_photon_attenuation

#ifdef MPI_USED
use mpi
#endif
