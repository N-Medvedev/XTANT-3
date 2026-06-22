!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! XTANT-3: X-ray-induced Thermal And Nonthermal Transitions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! 000000000000000000000000000000000000000000000000000000000000
! The hybrid code is written by
!
! Dr. Nikita Medvedev
!
! The model is described in: 
! https://arxiv.org/abs/1805.07524
!
! Should you have any questions, contact the author: nikita.medvedev@fzu.cz
! Or by private email: n.a.medvedev@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONVENTIONS OF PROGRAMMING:
! 1) All global variables start with "g_", e.g. g_numpar, and all defined in the module "Variables"
! 2) All modular variable names are defined starting as "m_", e.g. "m_number"
! 3) All local variables used within subrounies should NOT start with "g_" or "m_"
! 4) Add a comment after each subroutine and function specifying in which module it can be found
! 5) Leave comments describing what EACH LINE of the code is doing
! 6) Each end(smth) statement should be commented to which block it belongs, e.g.: if (i<k) then ... endif ! (i<k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM XTANT
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Initiate modules with all the 'use' statements collected in a separate file:
#include "Use_statements.f90"  ! include part of the code from an external file

implicit none


!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! MPI initialization (only if MPI is present, checked automatically via preprocessing):
call initialize_MPI(g_numpar%MPI_param, g_Err)   ! module "MPI_subroutines"
! Initialize ScaLAPACK:
call Initialize_ScaLAPACK(g_numpar%MPI_param, g_Err)   ! module "MPI_subroutines"
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Print XTANT label on the screen

if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
#ifdef _OPENMP
   call XTANT_label(6, 1)   ! module "Dealing_with_output_files"
#else
#ifdef MPI_USED
   call XTANT_label(6, 10)   ! module "Dealing_with_output_files"
#else ! no parallelzarion of any kind
   call XTANT_label(6, 4)   ! module "Dealing_with_output_files"
#endif
#endif
endif

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Set some starting default values:
g_numpar%which_input = 0 ! starting with the default input files
g_numpar%allow_rotate = .false. ! do not allow rotation of the target, remove angular momentum from initial conditions
1984 call Save_error_details(g_Err, 0, '', empty=.true.) ! module "Objects"
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   open(UNIT = g_Err%File_Num, FILE = trim(adjustl(m_Error_log_file)))
endif

! Check if the user needs any additional info (by setting the flags):
call get_add_data(g_numpar%path_sep, change_size=g_numpar%change_size, contin=g_Err%Stopsignal, &
      allow_rotate=g_numpar%allow_rotate, verbose=g_numpar%verbose, nonverbose=g_numpar%nonverbose) ! module "Read_input_data"
!--------------------------------------------------------------
! Master thread shares read info with all the other MPI-processes:
call MPI_share_add_data(g_numpar, g_Err)  ! module "MPI_subroutines"
!--------------------------------------------------------------
if (g_Err%Err) goto 2016     ! if something went wrong, cannot proceed
if (g_Err%Stopsignal) goto 2016     ! if the USER does not want to run the calculations, stop
! Otherwise, run the calculations:

! Make sure each MPI process is using a different random seed:
call initialize_random_seed(g_numpar%MPI_param)   ! module "MPI_subroutines"
!call random_seed() ! standard FORTRAN seeding of random numbers

call date_and_time(values=g_c1) ! standard FORTRAN time and date
g_ctim=g_c1 ! save the timestamp
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   call print_time('Attempting to start XTANT at', ind=0) ! prints out the current time, module "Little_subroutines"
endif

! Set all the initial data, read and create files:
! Read input files:
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   if (g_numpar%which_input > 0) then ! it's not the first run
      print*, '# It is run for input files number:', g_numpar%which_input
      call Read_Input_Files(g_matter, g_numpar, g_laser, g_Scell, g_Err, g_numpar%which_input) ! module "Read_input_data"
   else ! it is the first run
      print*, '# It is the first run'
      call Read_Input_Files(g_matter, g_numpar, g_laser, g_Scell, g_Err) ! module "Read_input_data"
   endif
endif ! (g_numpar%MPI_param%process_rank == 0)

!--------------------------------------------------------------
! Master thread shares read info with all the other MPI-processes:
call MPI_share_Read_Input_Files(g_matter, g_numpar, g_laser, g_Scell, g_Err)  ! module "MPI_subroutines"
!--------------------------------------------------------------
if (g_Err%Err) goto 2012   ! if there was an error in the input files, cannot continue, go to the end...
if (g_Err%Stopsignal) goto 2016     ! if the USER does not want to run the calculations, stop

! Printout additional info, if requested:
if (g_numpar%verbose) call print_time_step('Input files read successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! if you set to use OpenMP in compiling: "make"
#ifdef _OPENMP
   call OMP_SET_DYNAMIC(0) ! standard openmp subroutine
   call OMP_SET_NUM_THREADS(g_numpar%NOMP) ! number of threads for openmp defined in INPUT_PARAMETERS.txt
#else ! if you set to use OpenMP in compiling: 'make OMP=no'
!   print*, 'No openmp to deal with...'
#endif

! Starting time, to give enough time for system to thermalize before the pulse:
call set_starting_time(g_laser, g_time, g_numpar%t_start, g_numpar%t_NA, g_numpar%t_Te_Ee) ! module "Little_subroutines"

! And check if user wants to reset it:
call reset_dt(g_numpar, g_matter, g_time)   ! module "Dealing_with_output_files"

! Prepare initial conditions (read supercell and atomic positions from the files):
call set_initial_configuration(g_Scell, g_matter, g_numpar, g_laser, g_MC, g_Err) ! module "Initial_configuration"
if (g_Err%Err) goto 2012   ! if there was an error in preparing the initial configuration, cannot continue, go to the end...
!--------------------------------------------------------------
! Master thread shares read info with all the other MPI-processes:
call MPI_share_initial_configuration(g_Scell, g_matter, g_numpar, g_laser, g_MC, g_Err)   ! module "MPI_subroutines"
!--------------------------------------------------------------
if (g_numpar%verbose) call print_time_step('Initial configuration set successfully:', msec=.true., MPI_param=g_numpar%MPI_param)


! ! Print the title of the program and used parameters on the screen:
! call Print_title(6, g_Scell, g_matter, g_laser, g_numpar, -1) ! module "Dealing_with_output_files"
! if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
!    call print_time('Start at', ind=0) ! prints out the current time, module "Little_subroutines"
! endif


if (g_numpar%verbose) call print_time_step('Getting electron MFPs:', msec=.true., MPI_param=g_numpar%MPI_param)
! Read (or create) electronic mean free paths (both, inelastic and elastic):
call get_MFPs(g_Scell, 1, g_matter, g_laser, g_numpar, g_Scell(1)%TeeV, g_Err) ! module "MC_cross_sections"
!--------------------------------------------------------------
! Master thread shares read info with all the other MPI-processes:
call MPI_share_electron_MFPs(g_matter, g_numpar, g_Err)   ! module "MPI_subroutines"
!--------------------------------------------------------------
if (g_Err%Err) goto 2012   ! if there was an error in the input files, cannot continue, go to the end...
if (g_numpar%verbose) call print_time_step('Electron mean free paths set successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Read (or create) photonic mean free paths:
call get_photon_attenuation(g_matter, g_laser, g_numpar, g_Err) ! module "MC_cross_sections"
!--------------------------------------------------------------
! Master thread shares read info with all the other MPI-processes:
call MPI_share_photon_attenuation(g_matter, g_numpar, g_Err)   ! module "MPI_subroutines"
!--------------------------------------------------------------
if (g_Err%Err) goto 2012   ! if there was an error in the input files, cannot continue, go to the end...
if (g_numpar%verbose) call print_time_step('Photon attenuation lengths set successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

if (.not.g_numpar%do_path_coordinate) then  ! only for real calculations, not for coordinate path
   call save_last_timestep(g_Scell) ! save atomic before making next time-step, module "Atomic_tools"
endif

! Process the laser pulse parameters:
call process_laser_parameters(g_Scell(1), g_matter, g_laser, g_numpar) ! module "Monte_Carlo"
if (g_numpar%verbose) call print_time_step('Laser pulse parameters converted successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Print the title of the program and used parameters on the screen:
call Print_title(6, g_Scell, g_matter, g_laser, g_numpar, -1) ! module "Dealing_with_output_files"
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   call print_time('Start at', ind=0) ! prints out the current time, module "Little_subroutines"
endif

! Create the folder where output files will be stored, and prepare the files:
call prepare_output_files(g_Scell, g_matter, g_laser, g_numpar, g_Scell(1)%TB_Hamil(1,1), g_Scell(1)%TB_Repuls(1,1), g_Err) ! module "Dealing_with_output_files"
if (g_Err%Err) goto 2012   ! if there was an error in preparing the output files, cannot continue, go to the end...
if (g_numpar%verbose) call print_time_step('Output directory prepared successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Create CDF-file with fitted oscillators (Ritchi-Howie), if required:
call printout_CDF_file(g_numpar, g_matter, g_Scell)   ! module "Dealing_with_output_files"

! Printout mean free paths, if required:
call printout_MFP_file(g_numpar, g_matter, g_Scell)   ! module "Dealing_with_output_files"

! Printout pump-laser photon spectrum, if required:
call printout_laser_spectrum(g_laser, g_numpar, g_matter)   ! module "Dealing_with_output_files"

! Collect all ploting scripts into one script to execute all together later:
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   select case (g_numpar%plot_engine)    ! which engine to use for plots
   case default ! Gnuplot
      call collect_gnuplots(trim(adjustl(g_numpar%path_sep)), trim(adjustl(g_numpar%output_path)), skip_execution=.true.) ! module "Gnuplotting"
   case ('py') ! Python
      call collect_python_plots(trim(adjustl(g_numpar%path_sep)), trim(adjustl(g_numpar%output_path)), skip_execution=.true.) ! module "Plots_python"
   case ('no')
      ! Don't create any plots
   end select
endif


!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Project-specific analysis of C60:
! call C60_vdW_vs_Coulomb(g_Scell, g_numpar, g_matter, layers=2) ! Module "TB"
! call C60_crystal_construction(g_Scell, g_matter) ! module "Atomic_tools"
! call Coulomb_beats_vdW(g_Scell, g_numpar) ! see below

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! If user set '-size' option to vary the super-cell size:
if (g_numpar%change_size) then
   call vary_size(g_numpar, g_matter, g_Scell, g_time, Err=g_Err) ! module "XTANT_analyses"
   if (g_Err%Stopsignal .or. g_Err%Err) goto 2012      ! if the USER does not want to run the calculations
endif


!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! If user set to calculate the coordinate path between two phases of material:
if (g_numpar%do_path_coordinate) then
   call coordinate_path(g_numpar, g_matter, g_Scell, g_time, g_Err)  ! module "XTANT_analyses"
   if (g_Err%Err .or. g_Err%Stopsignal) goto 2012      ! if the USER does not want to run the calculations
endif
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! After the initial data are read, and necessary files created,
! now we can proceed with the real calculations

! Contruct TB Hamiltonian, diagonalize to get energy levels, get forces for atoms and supercell:
call get_Hamilonian_and_E(g_Scell, g_numpar, g_matter, 1, g_Err, g_time) ! module "TB"
if (g_numpar%verbose) call print_time_step('Initial Hamiltonian prepared successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Thermalization step for low-energy electrons (used only in relaxation-time approximation):
call Electron_thermalization(g_Scell, g_numpar, skip_thermalization=.true.) ! module "Electron_tools"

! Get global energy of the system at the beginning:
call get_glob_energy(g_Scell, g_matter) ! module "Electron_tools"
! and update the electron distribution:
call update_fe(g_Scell, g_matter, g_numpar, g_time, g_Err) ! module "Electron_tools"
if (g_numpar%verbose) call print_time_step('Initial energy prepared successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Get parameters that use complex Hamiltonian: DOS, CDF, kappa:
if ((abs(g_numpar%optic_model) > 0) .and. (g_numpar%optic_model < 4) ) then ! Trani or similar model, old-style DOS
   call get_optical_parameters(g_numpar, g_matter, g_Scell, g_Err) ! module "Optical_parameters"
   if (g_numpar%verbose) call print_time_step('Optical parameters prepared successfully:', msec=.true., MPI_param=g_numpar%MPI_param)
endif
!call get_DOS(g_numpar, g_matter, g_Scell, g_Err)   ! module "TB"
!if (g_numpar%verbose) call print_time_step('DOS calculated successfully:', msec=.true.)
call use_complex_Hamiltonian(g_numpar, g_matter, g_Scell, 1, g_Err)  ! module "TB_complex"
if (g_numpar%verbose) call print_time_step('Complex-Hamiltonian-dependent parameters (DOS, CDF, k) done:', msec=.true., MPI_param=g_numpar%MPI_param)


! Get current Mulliken charges (average and individual), if required:
call get_Mullikens_all(g_Scell(1), g_matter, g_numpar)
if (g_numpar%verbose) call print_time_step('Mulliken charges calculated successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Get the pressure in the atomic system:
call Get_pressure(g_Scell, g_numpar, g_matter, g_Scell(1)%Pressure, g_Scell(1)%Stress)	! module "TB"
if (g_numpar%verbose) call print_time_step('Pressure calculated successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Calculate the mean square displacement of all atoms:
call get_mean_square_displacement(g_Scell, g_matter, g_Scell(1)%MSD,  g_Scell(1)%MSDP, g_numpar%MSD_power, g_numpar)    ! module "Atomic_tools"
if (g_numpar%verbose) call print_time_step('Mean displacement calculated successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Calculate diffraction peaks:
call get_diffraction_peaks(g_Scell, g_matter, g_numpar)  ! module "Atomic_tools"
if (g_numpar%verbose) call print_time_step('Diffraction peaks calculated successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Calculate electron heat capacity, entropy, and orbital-resolved data:
call get_electronic_thermal_parameters(g_numpar, g_Scell, 1, g_matter, g_Err) ! module "TB"

! And save the (low-energy part of the) distribution on the grid, if required
! (its high-energy part is inside of MC_Propagate subroutine):
call get_low_energy_distribution(g_Scell(1), g_numpar) ! module "Electron_tools"

! Get atomic distribution:
call get_atomic_distribution(g_numpar, g_Scell, 1, g_matter)   ! module "Atomic_thermodynamics"
if (g_numpar%verbose) call print_time_step('Atomic distribution calculated successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

! Update configurational temperature for running average (needed on each timestep):
call update_Ta_config_running_average(g_Scell(1), g_matter, g_numpar)   ! module "Atomic_thermodynamics"
! Update local parameters of fragments (needed for local coupling):
call get_fragments_data(g_Scell, 1, g_numpar, g_matter)  ! module "Atomic_thermodynamics"
call get_fragments_data_for_electrons(g_Scell, 1, g_numpar, g_matter)  ! module "Electron_tools"

! Calculate configurational temperature (implemented only for Pettifor TB):
call Get_configurational_temperature(g_Scell, g_numpar, g_matter)	! module "TB"

! Save initial step in output:
call write_output_files(g_numpar, g_time, g_matter, g_Scell) ! module "Dealing_with_output_files"
if (g_numpar%verbose) call print_time_step('Initial output files set successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Signal possible warning for parameters defined:
! Print on the screen:
call check_all_warnings(6, g_laser, g_Scell, g_numpar)  ! module "Read_input_data"
! Save in the Error file:
call check_all_warnings(g_Err%File_Num, g_laser, g_Scell, g_numpar, g_Err)  ! module "Read_input_data"

!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
! Now we can proceed with time:
! Print out the starting time:
if (.not.g_numpar%nonverbose) call print_time_step('Simulation time:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)   ! module "Little_subroutines"

i_test = 0 !  count number of timesteps
g_dt_save = 0.0d0
do while (g_time .LT. g_numpar%t_total)
   i_test = i_test + 1
   ! If there is a grid for changing time-step, change it:
   call reset_dt(g_numpar, g_matter, g_time)  ! module "Dealing_with_output_files"

   AT_MOVE_1:if (g_numpar%do_atoms) then ! atoms are allowed to be moving:

      ! Do coupling before MD step:
      ! Nonadiabatic electron-ion coupling:
      !call Electron_ion_coupling(g_time, g_matter, g_numpar, g_Scell, g_Err) !  module "TB"
      !if (g_numpar%verbose) call print_time_step('Electron_ion_coupling succesful:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)


      !1111111111111111111111111111111111111111111111111111111
      ! Update atomic data on previous timestep and move further:
      call save_last_timestep(g_Scell) ! module "Atomic_tools"

      ! Make the MD timestep (first part, in the case of Verlet):
      call MD_step(g_Scell, g_matter, g_numpar, g_time, g_Err)  ! module "TB"
      if (g_numpar%verbose) call print_time_step('First step of MD step succesful:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)

      !2222222222222222222222222222222222222222222222222222222
      ! Nonadiabatic electron-ion coupling:
      call Electron_ion_coupling(g_time, g_matter, g_numpar, g_Scell, g_Err) !  module "TB"
      if (g_numpar%verbose) call print_time_step('Electron_ion_coupling succesful:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)

      ! Quenching of atoms (zero-temperature MD):
      call Cooling_atoms(g_numpar, g_matter, g_Scell, g_time, g_numpar%at_cool_dt, g_numpar%at_cool_start, g_numpar%do_cool) ! module "Atomic_tools"
      if (g_numpar%verbose) call print_time_step('Cooling of atoms done successfully:', msec=.true., MPI_param=g_numpar%MPI_param)

      ! Berendsen thermostat (mimicing energy transport; only if included):
      if (g_numpar%Transport_e) then ! for electrons
         ! Include Berendsen thermostat in the electronic system:
         call Electron_transport(1, g_time, g_Scell, g_numpar, g_matter, g_numpar%dt, g_matter%tau_bath_e, g_Err) ! module "Transport"
         if (g_numpar%verbose) call print_time_step('Electron Berendsen thermostat succesful:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)
      endif
      if (g_numpar%Transport) then ! for atoms
         ! Include Berendsen thermostat in the atomic system:
         call Atomic_heat_transport(1, g_Scell, g_matter, g_numpar%dt, g_matter%tau_bath) ! module "Transport"
         ! Include change of the affected layer for calculation of optical constants:
         call Change_affected_layer(1, g_Scell(1)%eps%dd, g_Scell, g_numpar%dt, g_matter%tau_bath)  ! module "Transport"
         if (g_numpar%verbose) call print_time_step('Atomic Berendsen thermostat succesful:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)

      endif
   endif AT_MOVE_1

   !write(*,'(a,f,f,f,f)') 'test 0:', g_time, SUM(g_Scell(1)%fe), SUM(g_Scell(1)%fe_eq), g_Scell(1)%Ne_low

   ! Monte-Carlo for photons, high-energy electrons, and core holes:
   call MC_Propagate(g_MC, g_numpar, g_matter, g_Scell, g_laser, g_time, g_Err) ! module "Monte_Carlo"
   if (g_numpar%verbose) call print_time_step('Monte Carlo model executed successfully:', g_time, msec=.true., MPI_param=g_numpar%MPI_param) ! module "Little_subroutines"

   !write(*,'(a,f,f,f,f)') 'test 1:', g_time, SUM(g_Scell(1)%fe), SUM(g_Scell(1)%fe_eq), g_Scell(1)%Ne_low

   ! Thermalization step for low-energy electrons (used only in nonequilibrium electron kinetics simulation):
   call Electron_thermalization(g_Scell, g_numpar) ! module "Electron_tools"

   !write(*,'(a,f,f,f,f)') 'test 2:', g_time, SUM(g_Scell(1)%fe), SUM(g_Scell(1)%fe_eq), g_Scell(1)%Ne_low

   ! And save the (low-energy part of the) distribution on the grid, if required
   ! (its high-energy part is inside of MC_Propagate subroutine):
   call get_low_energy_distribution(g_Scell(1), g_numpar) ! module "Electron_tools"

   ! Update corresponding energies of the system:
   call update_nrg_after_change(g_Scell, g_matter, g_numpar, g_time, g_Err) ! module "TB"

   !3333333333333333333333333333333333333333333333333333333333
   AT_MOVE_2:if (g_numpar%do_atoms) then ! atoms are allowed to be moving:
      ! Choose which MD propagator to use:
      select case(g_numpar%MD_algo)
      case default  ! velocity Verlet (2d order):
         !velocities update in the Verlet algorithm:
         call save_last_timestep(g_Scell) ! module "Atomic_tools"
         ! Atomic Verlet step:
         call make_time_step_atoms(g_Scell, g_matter, g_numpar, 1)     ! module "Atomic_tools"
         ! Supercell Verlet step:
         call make_time_step_supercell(g_Scell, g_matter, g_numpar, 1) ! supercell Verlet step, module "Atomic_tools"
         ! Update corresponding energies of the system:
         call get_new_energies(g_Scell, g_matter, g_numpar, g_time, g_Err) ! module "TB"

         ! Save numerical acceleration (2d half); (unused!):
         !call numerical_acceleration(g_Scell(1), g_numpar%halfdt, add=.true.)  ! module "Atomic_tools"
      case (1)  ! Yoshida (4th order)
         ! No divided steps, all of them are performed above
      case (2)  ! Martyna (4th order)
         ! No divided steps for atoms, but use Verlet for supercell:
         !velocities update in the Verlet algorithm:
         call save_last_timestep(g_Scell) ! module "Atomic_tools"
         ! Supercell Verlet step:
         call make_time_step_supercell(g_Scell, g_matter, g_numpar, 1) ! supercell Verlet step, module "Atomic_tools"
         ! Update corresponding energies of the system:
         call get_new_energies(g_Scell, g_matter, g_numpar, g_time, g_Err) ! module "TB"
      endselect
   endif AT_MOVE_2
   if (g_numpar%verbose) call print_time_step('Second step of MD step succesful:', g_time, msec=.true., MPI_param=g_numpar%MPI_param) ! module "Little_subroutines"

   !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   g_time = g_time + g_numpar%dt        ! [fs] next time-step
   g_dt_save = g_dt_save + g_numpar%dt  ! [fs] for tracing when to save the output data
   ! Update configurational temperature for running average (needed on each timestep):
   call update_Ta_config_running_average(g_Scell(1), g_matter, g_numpar)   ! module "Atomic_thermodynamics"
   ! Update local parameters of fragments (needed for local coupling):
   call get_fragments_data(g_Scell, 1, g_numpar, g_matter)  ! module "Atomic_thermodynamics"
   call get_fragments_data_for_electrons(g_Scell, 1, g_numpar, g_matter)  ! module "Electron_tools"

   ! Write current data into output files:
   if (g_dt_save .GE. g_numpar%dt_save - 1d-6) then
      ! Print out the curent time-step
      if (.not.g_numpar%nonverbose) call print_time_step('Simulation time:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)   ! module "Little_subroutines"

      ! Get parameters that use complex Hamiltonian: DOS, CDF, kappa:
      if ((abs(g_numpar%optic_model) > 0) .and. (g_numpar%optic_model < 4) ) then ! Trani or similar model, old-style DOS
         call get_optical_parameters(g_numpar, g_matter, g_Scell, g_Err) ! module "Optical_parameters"
      endif
      call use_complex_Hamiltonian(g_numpar, g_matter, g_Scell, 1, g_Err)  ! module "TB_complex"
!      call get_DOS(g_numpar, g_matter, g_Scell, g_Err)   ! module "TB"

      ! Get current Mulliken charges, if required:
      call get_Mullikens_all(g_Scell(1), g_matter, g_numpar)   ! module "TB"

      ! Get current pressure in the system:
      call Get_pressure(g_Scell, g_numpar, g_matter, g_Scell(1)%Pressure, g_Scell(1)%Stress)	! module "TB"

      ! Get atomic distributions and temperatures:
      call get_atomic_distribution(g_numpar, g_Scell, 1, g_matter)   ! module "Atomic_thermodynamics"

      ! Calculate the mean square displacement of all atoms:
      call get_mean_square_displacement(g_Scell, g_matter, g_Scell(1)%MSD, g_Scell(1)%MSDP, g_numpar%MSD_power, g_numpar) ! module "Atomic_tools"

      ! Calculate diffraction peaks:
      call get_diffraction_peaks(g_Scell, g_matter, g_numpar)  ! module "Atomic_tools"

      ! Calculate electron heat capacity, entropy, and orbital-resolved data:
      call get_electronic_thermal_parameters(g_numpar, g_Scell, 1, g_matter, g_Err) ! module "TB"

      ! Calculate configurational temperature:
      call Get_configurational_temperature(g_Scell, g_numpar, g_matter)	! module "TB"

      ! Get testmode additional data (center-of-mass, rotation, total force, etc.):
      call Get_testmode_add_data(g_Scell, 1, g_numpar, g_matter)	! module "Atomic_tools"

      ! Save current output data:
      call write_output_files(g_numpar, g_time, g_matter, g_Scell)    ! module "Dealing_with_output_files"
      ! Communicate with the program (program reads your commands from the communication-file):
      call communicate(g_numpar%FN_communication, g_time, g_numpar, g_matter) ! module "Dealing_with_output_files"
      g_dt_save = 0.0d0

      if (g_numpar%verbose) call print_time_step('Output files written successfully:', g_time, msec=.true., MPI_param=g_numpar%MPI_param)   ! module "Little_subroutines"
   endif
   !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
enddo

!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! Finish execution of the program:
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   call close_file('delete', FN=g_numpar%FN_communication, File_name=g_numpar%Filename_communication) ! module "Dealing_with_files"
endif

2012 continue


! Redo the pump-laser photon spectrum, now including the MC_sampled one:
call printout_laser_spectrum(g_laser, g_numpar, g_matter)   ! module "Dealing_with_output_files"


if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   INQUIRE(UNIT = g_Err%File_Num, opened=file_opened, name=chtest)
   if (file_opened) then
      flush(g_Err%File_Num)
   endif
endif

! Closing the opened files:
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   if (g_Err%Err .or. g_Err%Warn) then ! error or warning was printed in the file
      call close_file('close', FN=g_Err%File_Num) ! module "Dealing_with_files"
   else ! if there was no error or warning, no need to keep the file, delete it
      call close_file('delete', FN=g_Err%File_Num) ! module "Dealing_with_files"
   endif
   call close_save_files()           ! module "Dealing_with_files"
   call close_output_files(g_Scell, g_numpar) ! module "Dealing_with_files"
endif

if (g_numpar%verbose) call print_time_step('Opened files closed successfully', msec=.true., MPI_param=g_numpar%MPI_param)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Convolve output files with finite duration of the probe pulse:
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   if ( (g_Scell(1)%eps%tau > 0.0d0) .and. (.not.g_Err%Err) ) then
      call convolve_output(g_Scell, g_numpar, g_matter)  ! module "Dealing_with_output_files"
      print*, 'Convolution with the probe pulse is performed'
   else
      print*, 'No convolution with the probe pulse was required'
   endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Printing out the duration of the program, starting and ending time and date:
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it

   call parse_time(chtest, c0_in=g_ctim) ! module "Little_subroutines"
   write(*,'(a)') trim(adjustl(m_starline))
   write(*,'(a,a)') 'Duration of execution of program: ', trim(adjustl(chtest))

   call save_duration(g_matter, g_numpar, trim(adjustl(chtest)), ctim=g_ctim) ! module "Dealing_with_output_files"

   call print_time('Started  at', ctim=g_ctim) ! module "Little_subroutines"
   call print_time('Finished at') ! module "Little_subroutines"
   write(*,'(a)') trim(adjustl(m_starline))
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (.not.g_Err%Err) then
   if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      ! Gnuplotting:
      write(*,'(a)')  'Executing scripts to create plots...'
      select case (g_numpar%plot_engine)    ! which engine to use for plots
      case default ! Gnuplot
         call execute_all_gnuplots(g_numpar, trim(adjustl(g_numpar%output_path))//trim(adjustl(g_numpar%path_sep))) ! module "Plots_gnuplot"
         if (g_numpar%verbose) call print_time_step('Gnuplot calles executed', msec=.true., MPI_param=g_numpar%MPI_param)
      case ('py') ! Python
         call execute_all_pyplots(g_numpar, trim(adjustl(g_numpar%output_path))//trim(adjustl(g_numpar%path_sep))) ! module "Plots_python"
         if (g_numpar%verbose) call print_time_step('Python scripts executed', msec=.true., MPI_param=g_numpar%MPI_param)
      case ('no')
         ! No plots to be created
         if (g_numpar%verbose) call print_time_step('No scripts for plotting required', msec=.true., MPI_param=g_numpar%MPI_param)
      end select

   endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Check if there is another set of input files to run next simulation:
 g_numpar%which_input = g_numpar%which_input + 1
 chtest = trim(adjustl(m_INPUT_directory))//g_numpar%path_sep//trim(adjustl(m_INPUT_MATERIAL))
 write(chtest2,'(i3)') g_numpar%which_input
 write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
 if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   inquire(file=trim(adjustl(chtest)),exist=file_exists)    ! check if input file excists
 endif

!--------------------------------------------------------------
! Master thread shares this info with all the other MPI-processes:
#ifdef MPI_USED
    call mpi_bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, g_numpar%MPI_param%ierror)  ! module "mpi"
    call MPI_error_wrapper(g_numpar%MPI_param%process_rank, g_numpar%MPI_param%ierror, 'Next input {file_exists}') ! module "MPI_subroutines"
#endif
    call MPI_barrier_wrapper(g_numpar%MPI_param)  ! module "MPI_subroutines"
!--------------------------------------------------------------

 if (.not.file_exists) then ! check if the short-name file is present
   chtest = trim(adjustl(m_INPUT_directory))//g_numpar%path_sep//trim(adjustl(m_INPUT_ALL))
   write(chtest2,'(i3)') g_numpar%which_input
   write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
   if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      inquire(file=trim(adjustl(chtest)),exist=file_exists)    ! check if input file excists
   endif
!--------------------------------------------------------------
! Master thread shares this info with all the other MPI-processes:
#ifdef MPI_USED
    call mpi_bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, g_numpar%MPI_param%ierror)  ! module "mpi"
    call MPI_error_wrapper(g_numpar%MPI_param%process_rank, g_numpar%MPI_param%ierror, 'Next input {file_exists#2}') ! module "MPI_subroutines"
#endif
    call MPI_barrier_wrapper(g_numpar%MPI_param)  ! module "MPI_subroutines"
!--------------------------------------------------------------
 endif

 if (file_exists) then ! one file exists

    if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      write(*,'(a)') trim(adjustl(m_starline))
      write(*,'(a,a)')  'Another set of input parameter files exists: '
      print*, trim(adjustl(chtest))
    endif

    if (.not.g_numpar%numpar_in_input) then ! a separate file with NumPar is required
      chtest = trim(adjustl(m_INPUT_directory))//g_numpar%path_sep//trim(adjustl(m_NUMERICAL_PARAMETERS))
      write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
      if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         inquire(file=trim(adjustl(chtest)),exist=file_exists)    ! check if input file excists
      endif
!--------------------------------------------------------------
! Master thread shares this info with all the other MPI-processes:
#ifdef MPI_USED
      call mpi_bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, g_numpar%MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(g_numpar%MPI_param%process_rank, g_numpar%MPI_param%ierror, 'Next input {file_exists#3}') ! module "MPI_subroutines"
#endif
      call MPI_barrier_wrapper(g_numpar%MPI_param)  ! module "MPI_subroutines"
!--------------------------------------------------------------

      if (file_exists) then ! second input file exists
         if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'and ', trim(adjustl(chtest))
         endif
      else ! Maybe reuse the original file, if identical parameters are required:
         if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'File ', trim(adjustl(chtest)), ' not found'
            chtest = trim(adjustl(m_INPUT_directory))//g_numpar%path_sep//trim(adjustl(m_NUMERICAL_PARAMETERS))//'.txt'
            inquire(file=trim(adjustl(chtest)),exist=file_exists)    ! check if input file excists
            print*, 'File ', trim(adjustl(chtest)), ' will be reused'
         endif
         !--------------------------------------------------------------
! Master thread shares this info with all the other MPI-processes:
#ifdef MPI_USED
         call mpi_bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, g_numpar%MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(g_numpar%MPI_param%process_rank, g_numpar%MPI_param%ierror, 'Next input {file_exists#4}') ! module "MPI_subroutines"
#endif
         call MPI_barrier_wrapper(g_numpar%MPI_param)  ! module "MPI_subroutines"
!--------------------------------------------------------------
      endif
    endif ! (.not.g_numpar%numpar_in_input)

    if (.not.g_numpar%numpar_in_input .and. .not.file_exists) then ! second input file required but does not exist:
       if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         write(*,'(a)') trim(adjustl(m_starline))
         write(*,'(a,a,a)')  'File ', trim(adjustl(chtest)), ' could not be found.'
         !write(*,'(a)')  'XTANT has done its duty, XTANT can go...'
       endif
    else ! file either not required, or it exists, we can continue
       if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         write(*,'(a)')    'Running XTANT again for these new parameters...'
         write(*,'(a)') trim(adjustl(m_starline))
       endif
       call deallocate_all() ! module "Variables"
       goto 1984 ! go to the beginning and run the program again for the new input files
    endif
 else
    if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      write(*,'(a)') trim(adjustl(m_starline))
      !write(*,'(a)')  'XTANT has done its duty, XTANT can go...'
    endif
 endif

2016 continue
if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
   if (g_Err%Err) then ! add a comforing message if something went wrong :-(
      write(*,'(a)')  'XTANT really tried but could not do it...'
      write(*,'(a)') trim(adjustl(m_starline))
      call print_a_comforting_message(6, g_numpar%path_sep)  ! module "Dealing_with_output_files"
   else ! if finished noramally
      write(*,'(a)')  'XTANT has done its duty, XTANT can go...'
      write(*,'(a)') trim(adjustl(m_starline))
   endif
endif

!-------------------------------------
#ifdef MPI_USED
   if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      print*, 'Finilizing MPI'
   endif
   !call BLACS_EXIT(0)   ! ScaLAPACK library
   call MPI_FINALIZE(g_numpar%MPI_param%ierror) ! MPI
   if (g_numpar%MPI_param%ierror /= 0) then
      write(*, *) 'Error finalizing MPI!'
   endif
#endif
!-------------------------------------


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 contains
 ! nothing
END PROGRAM XTANT
