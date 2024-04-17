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
! This module contains subroutines to read input files:

MODULE Read_input_data
use Objects
use Universal_constants
use Little_subroutines, only : print_time_step, it_is_number, convert_wavelength_to_hw, convert_frequency_to_hw
use Dealing_with_files, only : Path_separator, Count_lines_in_file, close_file, copy_file, read_file, get_file_extension, &
                              ensure_correct_path_separator
use Dealing_with_EADL, only : m_EADL_file, m_EPDL_file, m_EEDL_file, READ_EADL_TYPE_FILE_int, READ_EADL_TYPE_FILE_real, select_imin_imax
use Dealing_with_DFTB, only : m_DFTB_directory, construct_skf_filename, read_skf_file, same_or_different_atom_types, &
                           idnetify_basis_size, m_DFTB_norep_directory, read_skf_file_no_rep
use Dealing_with_BOP, only : m_BOP_directory, m_BOP_file, read_BOP_parameters, idnetify_basis_size_BOP, &
                            read_BOP_repulsive, check_if_repulsion_exists
use Dealing_with_3TB, only : m_3TB_directory, m_3TB_onsite_data, read_3TB_onsite_file , construct_3TB_filenames, &
                            read_3TB_2bdy_file, read_3TB_3bdy_file
use Dealing_with_xTB, only : m_xTB_directory, read_xTB_parameters, identify_basis_size_xTB, identify_AOs_xTB
use Periodic_table, only : Decompose_compound
use Algebra_tools, only : make_cubic_splines, cubic_function
use Dealing_with_CDF, only : read_CDF_file
use Atomic_tools, only : update_atomic_masks_displ

! Open_MP related modules from external libraries:
#ifdef OMP_inside
   USE IFLPORT, only : system
   USE OMP_LIB, only : omp_get_max_threads
#endif

implicit none
PRIVATE


! Modular parameters:
character(25) :: m_INPUT_directory, m_INPUT_MATERIAL, m_NUMERICAL_PARAMETERS, m_INPUT_MINIMUM, m_INPUT_ALL, m_Atomic_parameters, &
                  m_Hubbard_U, m_Communication, m_COPY_INPUT

character(70), parameter :: m_starline = '*************************************************************'
character(70), parameter :: m_dashline = '-------------------------------------------------------------'
character(70), parameter :: m_warnline = ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
character(*), parameter :: m_numbers = '0123456789'
character(*), parameter :: m_LowCase = 'abcdefghijklmnopqrstuvwxyz'
character(*), parameter :: m_UpCase  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

character(15), parameter :: m_INFO_directory = 'INFO'  ! folder with the help-texts
character(15), parameter :: m_INFO_file = 'INFO.txt'  ! file with some XTANT info
character(15), parameter :: m_HELP_file = 'HELP.txt'  ! file with the helpful info
character(15), parameter :: m_QUOTES_file = 'Quotes.txt' ! file with comforting quotes
character(25), parameter :: m_List_ofmaterials = 'List_of_materials.txt'  ! list of materials existing

parameter (m_INPUT_directory = 'INPUT_DATA')    ! directory with all the input data
parameter (m_INPUT_MATERIAL = 'INPUT_MATERIAL') ! old format, material and pulse parameters
parameter (m_NUMERICAL_PARAMETERS = 'NUMERICAL_PARAMETERS') ! old format, numerical parameters
parameter (m_INPUT_MINIMUM = 'INPUT_MINIMUM') ! format with only parameters different from default (to be depricated)
parameter (m_INPUT_ALL = 'INPUT') ! new format with all parameters together
parameter (m_COPY_INPUT = 'Copy_input.txt') ! file with info on preparing multiple input files by copying INPUT
parameter (m_Atomic_parameters = 'Atomic_parameters') ! data-file with atomic parameters
parameter (m_Hubbard_U = 'INPUT_Hubbard_U.dat') ! data-file with Hubbard-U parameters (for SCC calculations)
parameter (m_Communication = 'Communication.txt')  ! file for comunication with the user

character(25), parameter :: m_short_pot = 'TB_short.txt'  ! filename for short-range potential
character(25), parameter :: m_wall_pot = 'TB_wall.txt'  ! obsolete filename for short-range potential

public :: m_INPUT_directory, m_INPUT_MATERIAL, m_NUMERICAL_PARAMETERS, m_Atomic_parameters, m_Hubbard_U
public :: m_INFO_directory, m_INFO_file, m_HELP_file, m_QUOTES_file, m_starline, m_INPUT_MINIMUM, m_INPUT_ALL
public :: Read_Input_Files, get_add_data, m_Communication, m_dashline, printout_warning, check_all_warnings

 contains


! These values will be used, unles changed by the user or in the input file:
subroutine initialize_default_values(matter, numpar, laser, Scell)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   !----------------------------------------------
   integer :: N
   ! Here we set default values, in case some of them are not given by the user:
   matter%Name = '' ! Material name
   numpar%Cell_filename = ''  ! no name given, defaults are hardcoded in nodule "Initial_configuration"
   matter%Chem = '' ! chemical formula of the compound
   if (.not.allocated(Scell)) allocate(Scell(1)) ! So far we only use 1 supercell
   numpar%numpar_in_input = .false.    ! assume separate file with the numerical parameters
   numpar%change_size_min = 0.7d0   ! default starting point for vary_size
   numpar%change_size_max = 2.1d0   ! default ending point for vary_size
   numpar%change_size_step = 300    ! default number of points for vary_size

   numpar%lin_scal = 0   ! do not use linear scaling TB (NOT READY)

   Scell(1)%Te = 300.0d0 ! initial electron temperature [K]
   Scell(1)%TeeV = Scell(1)%Te/g_kb ! [eV] electron temperature
   Scell(1)%Ta = 300.0d0 ! initial atomic temperature [K]
   Scell(1)%TaeV = Scell(1)%Ta/g_kb ! [eV] atomic temperature
   Scell(1)%Ta_var(:) = 0.0d0    ! various definitions of temperatures
   Scell(1)%Ta_r_var(:) = 0.0d0  ! projections of temperatures
   Scell(1)%Ta_conf_run_average(:) = 0.0d0   ! last values of Ta_config

   numpar%t_total = 1000.0d0 ! total duration of simulation [fs]
   call initialize_default_laser(laser, 1) ! initialize 1 pulse by default
   numpar%optic_model = 0 ! no optical calculations by default
   numpar%do_drude = .false.	! excluded optical calculations
   Scell(1)%eps%all_w = .false.	! excluded spectral calculations
   Scell(1)%eps%KK = .false.	! do not use Kramers Kronig relation
   Scell(1)%eps%E_min = 0.05d0 ! starting point of the grid of energy [eV]
   Scell(1)%eps%E_max = 50.0d0 ! ending point of the grid of energy [eV]
   Scell(1)%eps%dE = 0.1d0    ! grid step of energy [eV]
   numpar%drude_ray = 0	! Absorbtion of how many rays (0=exclude, 1=1st ray, (>1)=sum all)
   Scell(1)%eps%l = 800.0d0	! probe-pulse wavelength [nm]
   Scell(1)%eps%tau = -10.0d0	! probe duration FWHM [fs]
   Scell(1)%eps%ReEps0 = 0.0d0	! to start with
   Scell(1)%eps%ImEps0 = 0.0d0	! to start with
   Scell(1)%eps%w = 2.0d0*g_Pi*g_cvel/(Scell(1)%eps%l*1d-9) ! [1/sec] frequency
   Scell(1)%eps%teta = 0.0d0	! Angle of prob-pulse with respect to normal [degrees]
   Scell(1)%eps%teta = Scell(1)%eps%teta*g_Pi/(180.0d0) !c [radians]
   Scell(1)%eps%dd = 100.0d0	! material thickness [nm]
   ! number of unit-cells in X,Y,Z: 
   matter%cell_x = 1
   matter%cell_y = 1
   matter%cell_z = 1
   numpar%At_base = 'EADL' ! where to take atomic data from (EADL, CDF, etc...)
   matter%dens = -1.0d0 ! [g/cm^3] density of the material (negative = use MD supercell to evaluate it)
   numpar%NMC = 30000	! number of iterations in the MC module
#ifdef OMP_inside
   numpar%NOMP = omp_get_max_threads()    ! number of processors available by default
#else ! if you set to use OpenMP in compiling: 'make OMP=no'
   numpar%NOMP = 1   ! unparallelized by default
#endif
   numpar%output_path = ''    ! to start with, no name to output given yet
   numpar%output_extra_name = '' ! no replacement text in the output folder name
   numpar%output_name_add = ''   ! no additional text in the output folder name
   numpar%redo_MFP = .false.     ! no need to recalculate mean free paths by default
   numpar%print_MFP = .false.    ! no need to printout mean free paths by default
   numpar%print_Ta = .false.  ! no need in various atomic temperature definitions
   numpar%ind_starting_V = 2  ! by default, set Maxwellian starting velocities
   numpar%vel_from_file = .false.   ! velosities are not read from file
   numpar%N_basis_size = 0    ! DFTB, BOP or 3TB basis set default (0=s, 1=sp3, 2=sp3d5)
   numpar%do_atoms = .true.   ! Atoms are allowed to move
   matter%W_PR = 25.5d0    ! Parrinello-Rahman super-vell mass coefficient
   numpar%dt = 0.01d0      ! Time step for MD [fs]
   numpar%halfdt = numpar%dt/2.0d0      ! dt/2, often used
   numpar%dtsqare = numpar%dt*numpar%halfdt ! dt*dt/2, often used
   numpar%dt3 = numpar%dt**3/6.0d0            ! dt^3/6, often used
   numpar%dt4 = numpar%dt*numpar%dt3/8.0d0    ! dt^4/48, often used
   numpar%MD_algo = 0       ! 0=Verlet (2d order); 1=Yoshida (4th order)
   numpar%dt_save = 1.0d0	! save data into files every [fs]
   numpar%p_const = .false.	! V=const
   !matter%p_ext = g_P_atm	! External pressure [Pa] (0 = normal atmospheric pressure)
   matter%p_ext = 0.0d0     ! No external pressure by default [Pa]
   numpar%el_ion_scheme = 0	! scheme (0=decoupled electrons; 1=enforced energy conservation; 2=T=const; 3=BO)
   numpar%t_Te_Ee = 1.0d-5	! when to start coupling
   numpar%NA_kind = -1	! -1=Landau; 0=no coupling, 1=dynamical coupling (2=Fermi-golden_rule)
   numpar%Nonadiabat = .true.  ! included
   numpar%tau_fe = 1.0d0   ! Characteristic electron relaxation time [fs]
   numpar%tau_fe_CB = -1.0d0  ! No separate thermalization of CB and VB by default
   numpar%tau_fe_VB = -1.0d0  ! No separate thermalization of CB and VB by default
   numpar%do_partial_thermal = .false. ! no band-resolved thermalization by default
   numpar%scc = .true.  ! included
   numpar%scc_gam_ind = 0  ! Wolf's Coulomb model
   numpar%scc_mix = 1.0d0  ! maximal mixing
   numpar%t_NA = 1.0d-3	! [fs] start of the nonadiabatic
   numpar%acc_window = 5.0d0	! [eV] acceptance window for nonadiabatic coupling:
   numpar%do_DOS = .false.    ! DOS calculation
   numpar%do_kappa = .false.  ! electron heat conductivity calculation
   numpar%do_kappa_dyn = .false.  ! dynamic electron heat conductivity calculation
   numpar%save_CDF = .false.    ! fitted oscillators CDF printout
   numpar%kappa_Te_min = 300.0d0 ! [K]
   numpar%kappa_Te_max = 30000.0d0  ! [K]
   numpar%kappa_dTe = 100.0d0 ! [K]
   numpar%kappa_model = 0  ! default model index
   numpar%do_cool = .false.	! quenching excluded
   numpar%at_cool_start = 2500.0	! starting from when [fs]
   numpar%at_cool_dt = 40.0	! how often [fs]
   numpar%Transport = .false. ! excluded heat transport
   matter%T_bath = 300.0d0	! [K] bath temperature for atoms
   matter%T_bath = matter%T_bath/g_kb	! [eV] thermostat temperature
   matter%T_bath_e = 300.0d0	! [K] bath temperature for electrons
   matter%T_bath_e = matter%T_bath_e/g_kb	! [eV] thermostat temperature
   matter%tau_bath = 300.0d0	! [fs] time constant of cooling for atoms
   matter%tau_bath_e = 300.0d0	! [fs] time constant of cooling for electrons
   numpar%E_cut = 10.0d0 ! [eV] cut-off energy for high
   numpar%E_cut_dynamic = .false. ! do not change E_cut
   numpar%E_work = 1.0d30 ! [eV] work function (exclude electron emission)
   numpar%E_Coulomb = 0.0d0 ! [eV] Coulomb attraction of electrons back to the material
   numpar%save_Ei = .false.	! excluded printout energy levels (band structure)
   numpar%save_DOS = .false.	! excluded calculation and printout of DOS
   numpar%Smear_DOS = 0.05d0	! [eV] default smearing for DOS calculations
   numpar%save_fe = .false.	! excluded printout distribution function
   numpar%save_fe_orb = .false.	! excluded printout orbital-resolved distribution function
   numpar%save_fe_grid = .false.	! excluded printout distribution function on the grid
   numpar%save_PCF = .false.	! excluded printout pair correlation function
   numpar%save_XYZ = .true.	! included printout atomic coordinates in XYZ format
   numpar%save_XYZ_extra = .false.  ! no additional properties to print in XYZ file
   numpar%save_CIF = .true.	! included printout atomic coordinates in CIF format
   numpar%save_raw = .true.	! included printout of raw data for atomic coordinates, relative coordinates, velocities
   numpar%NN_radius = 0.0d0 ! radius of nearest neighbors defined by the user [A]
   numpar%MSD_power = 1     ! by default, print out mean displacement [A^1]
   numpar%save_NN = .false. ! do not print out nearest neighbors numbers
   numpar%do_elastic_MC = .true. ! allow elastic scattering of electrons on atoms within MC module
   numpar%r_periodic(:) = .true. ! use periodic boundaries along each direction of the simulation box
   ! Setting supercell for biomolecules, embedding in water:
   numpar%embed_water = .false.  ! no water added
   numpar%N_water_mol = 100      ! default number of water molecules
   ! number of k-points in each direction (used only for Trani-k!):
   numpar%ixm = 1
   numpar%iym = 1
   numpar%izm = 1
   ! BOP parameters:
   numpar%create_BOP_repulse = .false.
   ! initial n and k of unexcited material (used for DRUDE model only!):
   Scell(1)%eps%n = 1.0d0
   Scell(1)%eps%k = 0.0d0
   ! [me] effective mass of CB electron and VB hole:
   Scell(1)%eps%me_eff = 1.0d0
   Scell(1)%eps%mh_eff = 1.0d0
   Scell(1)%eps%me_eff = Scell(1)%eps%me_eff*g_me	! [kg]
   Scell(1)%eps%mh_eff = Scell(1)%eps%mh_eff*g_me	! [kg]
   ! [fs] mean scattering times of electrons and holes:
   Scell(1)%eps%tau_e = 1.0d0
   Scell(1)%eps%tau_h = 1.0d0
   ! Potential DOS power:
   numpar%power_b = 6.0d0 ! default value assuming U~1/r^5 (for no particular reason, empirically chosen)
   ! Names of the EPICS files:
   numpar%EADL_file = m_EADL_file ! default name, module "Dealing_with_EADL"
   numpar%EPDL_file = m_EPDL_file ! default name, module "Dealing_with_EADL"
   numpar%EEDL_file = m_EEDL_file ! default name (UNUSED), module "Dealing_with_EADL"
end subroutine initialize_default_values


subroutine initialize_default_laser(laser, N)
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   integer, intent(in) :: N	! how many pulses
   integer i
   
   if (.not.allocated(laser)) then
      allocate(laser(N))
   else if (size(laser) /= N) then
      deallocate(laser)
      allocate(laser(N))
   endif
   
   do i = 1,N
      call initialize_default_single_laser(laser, i)
   enddo
end subroutine initialize_default_laser


subroutine initialize_default_single_laser(laser, i) ! Must be already allocated
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   integer, intent(in) :: i	! number of the pulse
   laser(i)%F_in = -1.0d0  ! incoming fluence in [J/cm^2]; negative means unspecified - use absorbed dose instead
   laser(i)%F = 0.0d0  ! ABSORBED DOSE IN [eV/atom]
   laser(i)%hw = 100.0d0  ! PHOTON ENERGY IN [eV]
   laser(i)%FWHM_hw = 0.0d0  ! distribution of photon energies [eV]
   laser(i)%t = 10.0d0	  ! PULSE FWHM-DURATION IN
   laser(i)%KOP = 1  	  ! type of pulse: 0=rectangular, 1=Gaussian, 2=SASE
   !laser(i)%t = laser(i)%t/2.35482	! make a gaussian parameter out of it
   laser(i)%t0 = 0.0d0	  ! POSITION OF THE MAXIMUM OF THE PULSE IN [fs]
end subroutine initialize_default_single_laser


subroutine extend_laser(laser, N_extend)
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   integer, intent(in), optional :: N_extend	! extend the array by this many elements
   !----------------------------
   type(Pulse), dimension(:), allocatable :: temp_laser
   integer :: i, N, add_N
   if (present(N_extend)) then
      add_N = N_extend	! specified extension value
   else
      add_N = 1	! by default, extend the array by 1
   endif
   
   if (.not.allocated(laser)) then 	! it's first time, no parameters given in laser
      allocate(laser(add_N)) 	! just allocate it
      do i = 1, add_N		! and set all default parameters
         call initialize_default_single_laser(laser, i)
      enddo
   else	! there are already parameters for some of the pulses
      N = size(laser)	! number of those pulses were N, they need to be saved
      allocate(temp_laser(N))	! they are saved in a temporary array
      temp_laser = laser	! save them here
      deallocate(laser)		! now we need to extend the dimension: deallocate it first
      allocate(laser(N+add_N))	! and reallocate with the new size
      laser(1:N) = temp_laser(1:N)	! get all the previous data back into their positions
      deallocate(temp_laser)	! free the memory from the temporary array
      do i = N, N+add_N		! set all default parameters for the new elements
         call initialize_default_single_laser(laser, i)
      enddo
   endif
end subroutine extend_laser



!subroutine Read_Input_Files(matter, numpar, laser, TB_Repuls, TB_Hamil, Err)
subroutine Read_Input_Files(matter, numpar, laser, Scell, Err, Numb)
   type(Solid), intent(out) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(out) :: laser	! Laser pulse parameters
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   ! For polymorphic variables:
!    class(TB_repulsive), dimension(:), allocatable, intent(out) :: TB_Repuls  ! parameters of the repulsive part of TB
!    class(TB_Hamiltonian), dimension(:), allocatable, intent(out) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(Error_handling), intent(inout) :: Err	! error save
   integer, intent(in), optional :: Numb ! number of input files to use
   !-----------------------------
   type(User_overwrite_data) :: user_data
   real(8) temp
   integer FN, Reason, count_lines, N, i
   logical :: file_exist, file_opened, read_well, new_format_exists
   character(100) Error_descript, Folder_name, File_name, File_name_NEW
   character(3) chnum
   
   !--------------------------------------------------------------------------
   ! In case the user didn't define something, the default values will be used
   ! Set the default values:
   call initialize_default_values(matter, numpar, laser, Scell)

   !--------------------------------------------------------------------------
   ! Now, read the input file:
   call Path_separator(numpar%path_sep) ! module "Dealing_with_files"
   !Folder_name = 'INPUT_DATA'//numpar%path_sep
   Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep
   numpar%input_path = Folder_name ! save the address with input files

   ! File in minimum-format:
   if (.not.present(Numb)) then ! first run, use default files:
      File_name = trim(adjustl(Folder_name))//trim(adjustl(m_INPUT_MINIMUM))//'.txt'
   else ! it's not the first run, use next set of parameters:
      File_name = trim(adjustl(Folder_name))//trim(adjustl(m_INPUT_MINIMUM))
      write(chnum,'(i3)') Numb
      write(File_name,'(a,a,a,a)') trim(adjustl(File_name)), '_', trim(adjustl(chnum)), '.txt'
   endif
   inquire(file=trim(adjustl(File_name)),exist=new_format_exists)
   
   !-------------------------------
   ! Read input parameters in various formats:
   NEW_FORMAT:if (new_format_exists) then ! minimum format (inconvenient, to be deprecated)
      call read_input_txt(File_name, Scell, matter, numpar, laser, Err) ! see above
      if (Err%Err) goto 3416

   !-------------------------------
   else NEW_FORMAT ! Then use old format of two files
      ! First read material and pulse parameters:
      if (.not.present(Numb)) then ! first run, use default files:
         !File_name = trim(adjustl(Folder_name))//'INPUT_MATERIAL.txt'
         File_name = trim(adjustl(Folder_name))//trim(adjustl(m_INPUT_MATERIAL))//'.txt'
      else ! it's not the first run, use next set of parameters:
         !File_name = trim(adjustl(Folder_name))//'INPUT_MATERIAL'
         File_name = trim(adjustl(Folder_name))//trim(adjustl(m_INPUT_MATERIAL))
         write(chnum,'(i3)') Numb
         write(File_name,'(a,a,a,a)') trim(adjustl(File_name)), '_', trim(adjustl(chnum)), '.txt'
      endif
      inquire(file=trim(adjustl(File_name)),exist=file_exist)

      ! Check the short name of the file, if needed:
      TWO_OR_ONE_FILE:if (.not.file_exist) then ! try the short name:
         if (.not.present(Numb)) then ! first run, use default files:
            File_name = trim(adjustl(Folder_name))//trim(adjustl(m_INPUT_ALL))//'.txt'
         else ! it's not the first run, use next set of parameters:
            File_name = trim(adjustl(Folder_name))//trim(adjustl(m_INPUT_ALL))
            write(chnum,'(i3)') Numb
            write(File_name,'(a,a,a,a)') trim(adjustl(File_name)), '_', trim(adjustl(chnum)), '.txt'
         endif
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
      endif TWO_OR_ONE_FILE


      INPUT_MATERIAL:if (file_exist) then
         call read_input_material(File_name, Scell, matter, numpar, laser, user_data, Err) ! see below
         if (Err%Err) goto 3416

         ! Check multiple input files need to be created:
         if (.not.present(Numb)) then ! first run, use default files:
            call multiply_input_files(trim(adjustl(Folder_name)), trim(adjustl(File_name)), numpar%verbose)   ! below
         endif
      else
         write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, the program terminates'
         call Save_error_details(Err, 1, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3416
      endif INPUT_MATERIAL

      ! Read numerical parameters:
      if (.not.numpar%numpar_in_input) then ! if the parameters were not provided in the INPUT-file, try a separate file:
         if (.not.present(Numb)) then ! first run, use default files:
            File_name = trim(adjustl(Folder_name))//trim(adjustl(m_NUMERICAL_PARAMETERS))//'.txt'
         else ! it's not the first run, use next set of parameters:
            File_name = trim(adjustl(Folder_name))//trim(adjustl(m_NUMERICAL_PARAMETERS))
            write(chnum,'(i3)') Numb
            write(File_name,'(a,a,a,a)') trim(adjustl(File_name)), '_', trim(adjustl(chnum)), '.txt'
         endif
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         ! Maybe reuse the default one, if identical parameters are to be used:
         if (.not.file_exist) then
            File_name = trim(adjustl(Folder_name))//trim(adjustl(m_NUMERICAL_PARAMETERS))//'.txt'
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
         endif

         NUMERICAL_PARAMETERS:if (file_exist) then
            call read_numerical_parameters(File_name, matter, numpar, laser, Scell, user_data, Err) ! see below
            if (Err%Err) goto 3416
         else
            write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, the program terminates'
            call Save_error_details(Err, 1, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3416
         endif NUMERICAL_PARAMETERS
      endif ! (.not.numpar%numpar_in_input)
   endif NEW_FORMAT
   !-------------------------------

   ! Check if the molecule needs to be embedded in water (as requested in input file):
   if (numpar%embed_water) then
      call add_water_to_chemical_formulae(matter, numpar)   ! below
   endif

   if (.not.allocated(Scell)) allocate(Scell(1)) ! for the moment, only one super-cell

   ! Do TB and MD part only if we want (supercell is larger than 0):
   do i = 1, size(Scell)
      ! Read atomic data:
      call read_atomic_parameters(matter, numpar, Err) ! below
      if (Err%Err) goto 3416  ! exit if something went wrong
      if (numpar%user_defined_E_gap > -1.0d-14) then   ! user provided bandgap value, use it:
         Scell(i)%E_gap = numpar%user_defined_E_gap ! [eV]
         ! And redefine the Ip for the valence band:
         matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)) = Scell(i)%E_gap  ![eV]
      else ! assume atomic energy level:
         Scell(i)%E_gap = matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip))  ! [eV] band gap at the beginning
      endif
      Scell(i)%N_Egap = -1 ! just to start with something
      ! Read TB parameters:
      if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then
         call read_TB_parameters(matter, numpar, Scell(i)%TB_Repuls, Scell(i)%TB_Hamil, &
                           Scell(i)%TB_Waals, Scell(i)%TB_Coul, Scell(i)%TB_Expwall, Err) ! below
      else ! do only MC part
          ! Run it like XCASCADE (untested, may not work!)
      endif
   enddo

   ! Check if the user provided atomic data to overwrite the default values:
   call overwrite_atomic_data(user_data, matter)   ! below

   ! Atomic density [1/cm^3] (may be overwritten in the module "Initial_configuration", if negative):
   matter%At_dens = matter%dens/(SUM(matter%Atoms(:)%Ma*matter%Atoms(:)%percentage)/(SUM(matter%Atoms(:)%percentage))*1d3)


   ! Check k-space grid file:
   call read_k_grid(matter, numpar, Err)	! below

   ! Read Hubbard U parameter for self-consistent charge (SCC), if requested:
   if (numpar%scc) then
      ! Get file with Hubbard parameter:
      File_name = trim(adjustl(Folder_name))//trim(adjustl(m_Atomic_parameters))//numpar%path_sep//trim(adjustl(m_Hubbard_U))
      inquire(file=trim(adjustl(File_name)),exist=file_exist)
      if (file_exist) then
         call read_SCC_Hubbard(File_name, matter, numpar%scc)  ! below
      else  ! if there is no Hubbard parameters file, cannot use SCC:
         print*, 'File '//trim(adjustl(File_name))//' not found, SCC cannot be used'
         numpar%scc = .false.
      endif

      ! Consistency check:
      if (size(matter%Atoms) < 2) then ! no need for self-consistent charge calculations
         numpar%scc = .false.
         print*, 'SCC: Self-consistent change calculations for elemental solids are excluded'
      else
         if (numpar%scc_mix < 0.1d0) numpar%scc_mix = 0.1d0
         if (numpar%scc_mix > 1.0d0) numpar%scc_mix = 1.0d0
      endif
   endif

3416 continue !exit in case if input files could not be read
end subroutine Read_Input_Files



subroutine add_water_to_chemical_formulae(matter, numpar)   ! below
   type(Solid), intent(inout) :: matter   ! all material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !-----------------------------
   character(200) :: Folder_name, error_message
   character(100) :: chem_form, temp
   character(3) :: atom_elem
   integer :: N_at, i_H, i_O, INFO, i, N_SC
   integer, dimension(:), allocatable :: at_numbers
   real(8), dimension(:), allocatable :: at_percentage
   character(3), dimension(:), allocatable :: at_short_names ! name of the element
   logical :: H_exists, O_exists

!     print*, 'Chemical formula old: ', trim(adjustl(matter%Chem))

   ! Check if the atoms from the water molecule are already present:
   chem_form = trim(adjustl(matter%Chem)) ! save chemical formula
   Folder_name = trim(adjustl(numpar%input_path))//trim(adjustl(m_Atomic_parameters))   !'Atomic_parameters'

   call Decompose_compound(Folder_name, chem_form, numpar%path_sep, INFO, error_message, matter%N_KAO, at_numbers, at_percentage, at_short_names) ! molude 'Periodic_table'

   if (INFO == 0) then
      H_exists = .false.   ! to start with
      O_exists = .false.   ! to start with

      N_SC = INT(matter%cell_x*matter%cell_y*matter%cell_z) ! number of molecules in the supercell
      matter%Chem = ''  ! to start with
      do i = 1, size(at_short_names)
         if (trim(adjustl(at_short_names(i))) == 'H') then
            ! That's how many hyhdrogens there are in total (molecule + water):
            write(temp,'(i6)') INT(at_percentage(i))*N_SC + numpar%N_water_mol*2
            matter%Chem = trim(adjustl(matter%Chem))//'H'//trim(adjustl(temp))
            H_exists = .true. ! found it, added into the chemical formula
         elseif (trim(adjustl(at_short_names(i))) == 'O') then
            ! That's how many oxygens there are in total (molecule + water):
            write(temp,'(i6)') INT(at_percentage(i))*N_SC + numpar%N_water_mol
            matter%Chem = trim(adjustl(matter%Chem))//'O'//trim(adjustl(temp))
            O_exists = .true. ! found it, added into the chemical formula
         else  ! all other elements just write back as they are (x number of unit cells)
            write(temp,'(i6)') INT(at_percentage(i))*N_SC
            matter%Chem = trim(adjustl(matter%Chem))//trim(adjustl(at_short_names(i)))//trim(adjustl(temp))
         endif
      enddo
      ! Add at the end H or O, if they were not in the formula:
      if (.not.H_exists) then
         write(temp,'(i6)') numpar%N_water_mol*2
         matter%Chem = trim(adjustl(matter%Chem))//'H'//trim(adjustl(temp))
      endif
      if (.not.O_exists) then
         write(temp,'(i6)') numpar%N_water_mol
         matter%Chem = trim(adjustl(matter%Chem))//'O'//trim(adjustl(temp))
      endif
   endif ! (INFO .NE. 0)
!     print*, 'Chemical formula new: ', trim(adjustl(matter%Chem))
!     pause 'Chem done'
end subroutine add_water_to_chemical_formulae



subroutine read_SCC_Hubbard(File_name, matter, SCC)
   character(*), intent(in) :: File_name  ! file with Hubbard parameters
   type(Solid), intent(inout) :: matter   ! all material parameters
   logical, intent(inout) :: SCC ! flag to use SCC or not
   !---------------------------
   real(8) :: U_read
   integer :: TOA, i, j, FN, count_lines, Reason
   logical :: file_opened, read_well, found_el
   character(3) :: El_name

   FN = 103
   open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      print*, 'File '//trim(adjustl(File_name))//' could not be opened, SCC cannot be used'
      SCC = .false.
      goto 3430
   endif

   RF:do i = 1, size(matter%Atoms) ! for all types of atoms
      found_el = .false.
      count_lines = 0   ! to start with
      read(FN,*,IOSTAT=Reason) ! skip first two lines with comments
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         print*, 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         print*, 'SCC cannot be used, proceeding at non-self-consistent level'
         SCC = .false.
         exit RF  ! problem reading file, no need to continue
      endif
      read(FN,*,IOSTAT=Reason) ! skip first two lines with comments
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         print*, 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         print*, 'SCC cannot be used, proceeding at non-self-consistent level'
         SCC = .false.
         exit RF  ! problem reading file, no need to continue
      endif

      SFE:do while (.not.found_el)
         read(FN,*,IOSTAT=Reason) El_name, U_read
         if (.not. read_well) then
            print*, 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            print*, 'SCC cannot be used, proceeding at non-self-consistent level'
            SCC = .false.
            exit RF  ! problem reading file, no need to continue
         endif
         ! Check if this is the element we need:
         if ( trim(adjustl(El_name)) == trim(adjustl(matter%Atoms(i)%Name)) ) then
            found_el = .true.
            !matter%Atoms(i)%Hubbard_U = U_read * 0.5d0   ! alternative definition with 1/2
            matter%Atoms(i)%Hubbard_U = U_read  ! [eV]

            rewind(FN)  ! start the search for the next element
            exit SFE ! found element, go to the next one
         endif
      enddo SFE
   enddo RF

3430 continue
   call close_file('close', FN=FN) ! module "Dealing_with_files"
end subroutine read_SCC_Hubbard



subroutine read_k_grid(matter, numpar, Err)
   type(Solid), intent(in) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------------------------
   integer :: FN, N, count_lines, Reason, i
   character(200) :: error_message, Error_descript
   character(200) :: Folder_name, Path, File_name
   logical :: file_exist, file_opened, read_well
   
   ! Check if we even need the k grid:
   select case (ABS(numpar%optic_model))  ! use multiple k-points, or only gamma
      case (-4, 2, 4:5) ! Models that use multiple k points
         Folder_name = trim(adjustl(numpar%input_path))
         Path = trim(adjustl(Folder_name))//trim(adjustl(matter%Name))
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), 'k_grid.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)

         if (file_exist) then
            FN = 103
            open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then
               print*, 'File '//trim(adjustl(File_name))//' could not be opened, using Monkhorst Pack k-points'
               goto 3426
            endif

            ! Count how many grid points are there:
            call Count_lines_in_file(FN, N)	! module "Dealing_with_files"
            
            ! Knwoing the number of k points, allocate the array:
            allocate(numpar%k_grid(N,3))
            numpar%k_grid = 0.0d0
            
            ! Also adjust the nunbers of k-points correspondingly:
            numpar%ixm = ceiling(dble(N)**(1.0d0/3.0d0))
            numpar%iym = numpar%ixm
            numpar%izm = numpar%ixm
            
            ! Read the k points from the file:
            count_lines = 0
            do i = 1, N
               read(FN,*,IOSTAT=Reason) numpar%k_grid(i,1), numpar%k_grid(i,2), numpar%k_grid(i,3)
               call read_file(Reason, count_lines, read_well)
               if (.not. read_well) then
                  print*, 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
                  print*, 'The remaining ', (N-count_lines), ' k points are set as Gamma point (0,0,0)'
                  goto 3426
               endif
            enddo
            
            close(FN)
         else
            print*, 'k-space grid file not found, using Monkhorst Pack k-points'
            goto 3426
         endif
      case default	! gamma point
         ! no need to even care about k space grid
   end select
3426 continue
end subroutine read_k_grid



subroutine read_atomic_parameters(matter, numpar, Err)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   integer i
   character(200) :: File_name
   logical :: file_exist
   
   select case (trim(adjustl(numpar%At_base)))
   case('CDF', 'cdf', 'CDF_sp') ! read data from corresponding *.cdf file

      ! Check if file with CDF oscillator parameters exists:
      call check_CDF_file_exists(numpar, matter, File_name, file_exist) ! below

      if (file_exist) then
         call get_CDF_data(matter, numpar, Err, File_name) ! see below
      else
         print*, 'File '//trim(adjustl(File_name))//' could not be found,', ' using single-pole CDF approximation'
         numpar%At_base = 'CDF_sp'
         call get_EADL_data(matter, numpar, Err) ! see below
      endif
   case ('XATOM') ! get data from XATOM code
      ! to be integrated with XATOM later...
   case default ! ('EADL', 'BEB'), read data from EADL database
      call get_EADL_data(matter, numpar, Err) ! see below
   end select
   
!    print*, 'TEST:', allocated(matter%Atoms), size(matter%Atoms), numpar%At_base
!    do i = 1, size(matter%Atoms)
!       !print*, trim(adjustl(matter%Chem)), ' '//trim(adjustl(matter%Atoms(i)%Name)), matter%Atoms(i)%N_CDF(:), matter%Atoms(i)%Shl_dsgnr(:), matter%Atoms(i)%Ip(:), matter%Atoms(i)%Ne_shell(:), matter%Atoms(i)%Auger(:)
!       print*, matter%Atoms(i)%Ne_shell(:)
!    enddo
!    call get_table_of_ij_numbers(matter, numpar) ! table of elements number to locate TB parameterization for different combinations of atoms
end subroutine read_atomic_parameters


subroutine get_CDF_data(matter, numpar, Err, File_name)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   character(*), intent(in) :: File_name  ! CDF file
   !===============================================
   integer, dimension(:), allocatable :: at_numbers
   real(8), dimension(:), allocatable :: at_percentage
   character(3), dimension(:), allocatable :: at_short_names ! name of the element
   character(25), dimension(:), allocatable :: at_names ! full name of the element
   real(8), dimension(:), allocatable :: at_masses ! mass of each element [Mp]
   integer, dimension(:), allocatable :: at_NVE    ! number of valence electrons
   integer FN1, INFO, i, j, k, Z
   character(100) :: error_message, Error_descript
   character(100) :: Folder_name, Folder_name2, Path
   logical file_exist, file_opened, read_well
   integer FN, Reason, count_lines, temp
   real(8) retemp

   ! Make sure this file exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exist)

   if (file_exist) then

      !open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), status = 'old')
      FN = 103
      open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      inquire(file=trim(adjustl(File_name)),opened=file_opened)
      if (.not.file_opened) then
         Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened'
         call Save_error_details(Err, 2, Error_descript)
         print*, trim(adjustl(Error_descript))
         return
      endif

      ! Read the file in the CDF format:
      !call read_CDF_file(FN, matter, numpar, Err, trim(adjustl(File_name)), &
      !                  trim(adjustl(m_Atomic_parameters)), trim(adjustl(m_EADL_file)), INFO) ! module "Dealing_with_CDF"
      call read_CDF_file(FN, matter, numpar, Err, trim(adjustl(File_name)), &
                        trim(adjustl(m_Atomic_parameters)), trim(adjustl(numpar%EADL_file)), INFO) ! module "Dealing_with_CDF"

      ! When done, close the CDF file
      close(FN)
   else  ! no file with CDF-oscillator parameters found
      write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, program terminates'
      call Save_error_details(Err, 1, Error_descript)
      print*, trim(adjustl(Error_descript))
      return
   endif
end subroutine get_CDF_data


subroutine check_CDF_file_exists(numpar, matter, File_name, file_exist)
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   type(Solid), intent(in) :: matter   ! all material parameters
   character(*), intent(inout) :: File_name
   logical, intent(inout) :: file_exist
   !------------------
   character(250) :: Folder_name, Path

   Folder_name = trim(adjustl(numpar%input_path))
   Path = trim(adjustl(Folder_name))//trim(adjustl(matter%Name))

   ! Try to find the file with CDF parameters:
   if (LEN(trim(adjustl(numpar%input_CDF_file))) > 0) then  ! there is user-provided name
      ! 1) assuming the full path is given:
      write(File_name, '(a)') trim(adjustl(numpar%input_CDF_file))
      inquire(file=trim(adjustl(File_name)),exist=file_exist)

      if (.not.file_exist) then
         ! 2) assume it is file-name in the input directory:
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(numpar%input_CDF_file))
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
      endif

      if (.not.file_exist) then
         ! 3) assume it is file-name without the extension in the input directory:
         write(File_name, '(a,a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(numpar%input_CDF_file)), '.cdf'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
      endif
   endif

   if (.not.file_exist) then
      ! 4) assume the default name:
      write(File_name, '(a,a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(matter%Name)), '.cdf'
      inquire(file=trim(adjustl(File_name)),exist=file_exist)
   endif

   if (numpar%verbose) call print_time_step('CDF-oscillators reading from:'//trim(adjustl(File_name)), msec=.true.) ! modlue "Little_subroutines"
end subroutine check_CDF_file_exists



subroutine get_EADL_data(matter, numpar, Err)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !===============================================
   integer, dimension(:), allocatable :: at_numbers
   real(8), dimension(:), allocatable :: at_percentage
   character(3), dimension(:), allocatable :: at_short_names ! name of the element
   character(25), dimension(:), allocatable :: at_names ! full name of the element
   real(8), dimension(:), allocatable :: at_masses ! mass of each element [Mp]
   integer, dimension(:), allocatable :: at_NVE    ! number of valence electrons
   integer FN1, INFO, i, j, Z, old_shl_num, mod_shl_num, shl_dsg
   real(8) E_gap
   character(100) :: error_message, Error_descript
   character(100) :: Folder_name, File_name
   logical file_exist, file_opened

   !Folder_name = 'INPUT_DATA'//trim(adjustl(numpar%path_sep))//'Atomic_parameters'
   Folder_name = trim(adjustl(numpar%input_path))//trim(adjustl(m_Atomic_parameters))   !'Atomic_parameters'

   call Decompose_compound(Folder_name, matter%Chem, numpar%path_sep, INFO, error_message, matter%N_KAO, at_numbers, at_percentage, at_short_names, at_names, at_masses, at_NVE) ! molude 'Periodic_table'
   if (INFO .NE. 0) then
      call Save_error_details(Err, INFO, error_message)
      print*, trim(adjustl(error_message))
      goto 3419
   endif
   if (.not.allocated(matter%Atoms)) allocate(matter%Atoms(matter%N_KAO))
   do i = 1, matter%N_KAO ! for all sorts of atoms
      matter%Atoms(i)%Z = at_numbers(i)
      matter%Atoms(i)%Name = at_short_names(i)
      !matter%Atoms(i)%Ma = at_masses(i)*g_Mp ! [kg]
      matter%Atoms(i)%Ma = at_masses(i)*g_amu ! [kg]
      matter%Atoms(i)%percentage = at_percentage(i)
      matter%Atoms(i)%NVB = at_NVE(i)
   enddo

   ! Open eadl.all database:
   !File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//'eadl.all'
   !File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(m_EADL_file))
   File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(numpar%EADL_file))
   inquire(file=trim(adjustl(File_name)),exist=file_exist)
   if (.not.file_exist) then
      Error_descript = 'File '//trim(adjustl(File_name))//' does not exist, the program terminates'
      call Save_error_details(Err, 1, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3419
   endif

!    print*, 'File ', trim(adjustl(File_name)), ' is read'

   !open(NEWUNIT=FN1, FILE = trim(adjustl(File_name)), status = 'old')
   FN1=105
   open(UNIT=FN1, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) then
     do i = 1, matter%N_KAO ! for all atomic kinds of the compound
      Z =  matter%Atoms(i)%Z ! atomic number

      ! First, get the total number of shells of this atom:
      call READ_EADL_TYPE_FILE_int(FN1, File_name, Z, 912, INFO, error_message, matter%Atoms(i)%sh, matter%Atoms(i)%Ne_shell, Shl_num=matter%Atoms(i)%Shl_dsgnr_atomic, Ip=matter%Atoms(i)%Ip)    ! module "Dealing_with_EADL"
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3419
      endif

!       print*, 'READ_EADL_TYPE_FILE_int is done', i

      ! Get values for all parameters:
      call READ_EADL_TYPE_FILE_real(FN1, File_name, Z, 913, matter%Atoms(i)%Ip, INFO=INFO, error_message=Error_descript) ! read binding energies, module "Dealing_with_EADL"
      ! Make sure we start reading it again from the start:
      rewind(FN1)

!       print*, 'READ_EADL_TYPE_FILE_real is done', i


      ! Then correct it to exclude VB:
      call exclude_VB(i, matter%Atoms(i)%sh, matter%Atoms(i)%Ne_shell, matter%Atoms(i)%NVB, mod_shl_num)
      old_shl_num = matter%Atoms(i)%sh ! save old uppermost level to use as initial band gap
      matter%Atoms(i)%sh = mod_shl_num ! new number of shells
!       print*, 'itest', i, matter%Atoms(i)%sh , mod_shl_num 
      E_gap = matter%Atoms(i)%Ip(old_shl_num) ! [eV] save uppermost level
      deallocate(matter%Atoms(i)%Ip) ! to use it next time for new number of shells

      ! Now get the parameters for all the shells (except VB):
      call READ_EADL_TYPE_FILE_int(FN1, File_name, Z, 912, INFO, error_message, matter%Atoms(i)%sh, matter%Atoms(i)%Ne_shell, Shell_name=matter%Atoms(i)%Shell_name, Shl_num=matter%Atoms(i)%Shl_dsgnr, Ip=matter%Atoms(i)%Ip, Ek=matter%Atoms(i)%Ek, Auger=matter%Atoms(i)%Auger, REDO=.false.)    ! module "Dealing_with_EADL"
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3419
      endif

!       print*, i, matter%Atoms(i)%Shl_dsgnr_atomic, matter%Atoms(i)%Shl_dsgnr
!       pause 'Shl_dsgnr_atomic'

!       print*, 'READ_EADL_TYPE_FILE_real is done: Auger', i

      ! Get values for all parameters:
      call READ_EADL_TYPE_FILE_real(FN1, File_name, Z, 913, matter%Atoms(i)%Ip, INFO=INFO, error_message=Error_descript) ! read binding energies, module "Dealing_with_EADL"
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3419
      endif

      call READ_EADL_TYPE_FILE_real(FN1, File_name, Z, 914, matter%Atoms(i)%Ek, INFO=INFO, error_message=Error_descript) ! read kinetic energies, module "Dealing_with_EADL"
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3419
      endif

      if (.not. allocated(matter%Atoms(i)%TOCS)) then
         allocate(matter%Atoms(i)%TOCS(matter%Atoms(i)%sh)) ! allocate type of cross-section to be used for electrons
         matter%Atoms(i)%TOCS = 0 ! do all BEB cross-sections
      endif

      if (.not. allocated(matter%Atoms(i)%TOCSph)) then
         allocate(matter%Atoms(i)%TOCSph(matter%Atoms(i)%sh), source = 0) ! allocate type of cross-section to be used for photons
      endif


      if (.not. allocated(matter%Atoms(i)%El_MFP)) allocate(matter%Atoms(i)%El_MFP(matter%Atoms(i)%sh)) ! allocate electron MFPs
      if (.not. allocated(matter%Atoms(i)%Ph_MFP)) allocate(matter%Atoms(i)%Ph_MFP(matter%Atoms(i)%sh)) ! allocate photon MFPs
      !call READ_EADL_TYPE_FILE_real(FN1, File_name, Z, 921, Target_atoms(i)%Radiat, INFO=INFO, error_message=Error_descript) ! radiative decay, module "Dealing_with_EADL"

      call READ_EADL_TYPE_FILE_real(FN1, File_name, Z, 922, matter%Atoms(i)%Auger, INFO=INFO, error_message=Error_descript) ! read auger-times, module "Dealing_with_EADL"
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3419
      endif

      do j = 1,size(matter%Atoms(i)%Ip)
         if (matter%Atoms(i)%Auger(j) <= 0.0d0) then ! no auger possible/nothing in the database
            matter%Atoms(i)%Auger(j) = 1d30
         else ! there is some value in the database, convert it into our units:
            matter%Atoms(i)%Auger(j)=1d15*g_h/(g_e*matter%Atoms(i)%Auger(j)) ! [fs]
         endif
         ! Check that all deep shell holes decay:
         if (j > 1) then ! all shells except K-shell for now:
            if (matter%Atoms(i)%Auger(j) >= 1d20) matter%Atoms(i)%Auger(j) = matter%Atoms(i)%Auger(j-1)
         endif
!          print*, i, j, matter%Atoms(i)%Auger(j)
      enddo ! j

      ! Now, modify the VB values where possible:
      if (i == 1) then
         matter%Atoms(i)%Shell_name(mod_shl_num) = 'Valence'
         matter%Atoms(i)%Shl_dsgnr(mod_shl_num) = 63
         matter%Atoms(i)%Ip(mod_shl_num) = E_gap ! VB only for 1st kind of atoms
         matter%Atoms(i)%Ek(mod_shl_num) = 0.0d0 ![eV]
         matter%Atoms(i)%Auger(mod_shl_num) = 1d23 ! [fs] no Auger-decays of VB
         matter%Atoms(i)%Ne_shell(mod_shl_num) = matter%Atoms(i)%NVB ! number of valence electrons
      endif
!       print*, 'Total numbert of shells:', matter%Atoms(i)%sh
!       do j = 1, matter%Atoms(i)%sh
!          print*, 'Names:', trim(adjustl(matter%Atoms(i)%Shell_name(j)))
!       enddo
!       print*, 'Shl_dsgnr:', matter%Atoms(i)%Shl_dsgnr
!       print*, 'Ne_shell:', matter%Atoms(i)%Ne_shell
!       print*, 'Ip:', matter%Atoms(i)%Ip
!       print*, 'Ek:', matter%Atoms(i)%Ek
!       print*, 'Auger:', matter%Atoms(i)%Auger
!       pause "get_EADL_data"

        rewind(FN1) ! for the next element, start over from the beginning of the file
     enddo ! i
     close(FN1)
   else 
     Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
     call Save_error_details(Err, 2, Error_descript)
     print*, trim(adjustl(Error_descript))
     goto 3419
   endif

3419 continue
end subroutine get_EADL_data


subroutine exclude_VB(NOA, sh, Ne_shell, NVB, mod_shl_num, N_shl)
   integer, intent(in) :: NOA   ! number of atom in the compound
   integer, intent(inout) :: sh ! number of shells of the element
   real(8), dimension(:), allocatable, intent(inout) :: Ne_shell ! number of electron in each shell
   integer, intent(in) :: NVB   ! number valence electrons according to periodic table
   integer, intent(out) :: mod_shl_num    ! modified number of shells
   integer, intent(in), optional :: N_shl ! if we already counted the number of shells 
   integer i, counter
   real(8) Ne_cur
   real(8), dimension(size(Ne_shell)) :: Ne_temp
   if (.not.present(N_shl)) then
      Ne_cur = 0.0d0 ! to start counting
      counter = 0
      SH_COUNT: do i = sh,1,-1
         Ne_cur = Ne_cur + Ne_shell(i)
         counter = counter + 1
         if (Ne_cur >= NVB) exit SH_COUNT
      enddo SH_COUNT
      if (NOA == 1) then 
         mod_shl_num = sh - counter + 1 ! all deep shells plus one for VB
      else
         mod_shl_num = sh - counter     ! all deep shells, but no VB for other sorts of atoms
      endif
   else
      mod_shl_num = N_shl ! all deep shells plus one for VB
   endif
!    if (mod_shl_num < 1) mod_shl_num = 1 ! for elements with just 1 shell filled (H and He)

!    print*, NVB, mod_shl_num, sh
!    print*, Ne_shell(:)
!    pause 'exclude_VB - 0'

   if (size(Ne_shell) /= mod_shl_num) then ! redefine it:
      Ne_temp = Ne_shell
      if (allocated(Ne_shell)) deallocate(Ne_shell)
      allocate(Ne_shell(mod_shl_num)) ! new size
      RESHAP: do i = 1, mod_shl_num-1 ! deep shells are correct
         Ne_shell(i) = Ne_temp(i)
      enddo RESHAP
      if (mod_shl_num > 0) Ne_shell(mod_shl_num) = NVB ! this is number of electrons in the VB
   endif


!     print*, NVB, mod_shl_num, sh
!     print*, Ne_shell(:)
!     pause 'exclude_VB'
end subroutine exclude_VB



subroutine read_TB_parameters(matter, numpar, TB_Repuls, TB_Hamil, TB_Waals, TB_Coul, TB_Expwall, Err)
   type(Solid), intent(in) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   ! For polymorphic variables:
   class(TB_repulsive), dimension(:,:), allocatable, intent(out) :: TB_Repuls   ! parameters of the repulsive part of TB
   class(TB_Hamiltonian), dimension(:,:), allocatable, intent(out) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   class(TB_vdW),  dimension(:,:), allocatable, intent(out) :: TB_Waals         ! parameters of the van der Waals for TB
   class(TB_Coulomb),  dimension(:,:), allocatable, intent(out) :: TB_Coul	! parameters of the Coulomb together with TB
   class(TB_Exp_wall),  dimension(:,:), allocatable, intent(out) :: TB_Expwall	! parameters of the exponential wall with TB
   type(Error_handling), intent(inout) :: Err	! error save
   !========================================================
   integer FN, count_lines, Reason, INFO, i, j !, N
   character(500) :: Error_descript, Folder_name, File_name, Path, ch_temp
   logical file_exists, file_opened, read_well
   !Folder_name = 'INPUT_DATA'//trim(adjustl(numpar%path_sep))
   Folder_name = trim(adjustl(numpar%input_path))
   Path = trim(adjustl(Folder_name))//trim(adjustl(matter%Name))
   
   do_first:do i = 1, matter%N_KAO
      do_second:do j = 1, matter%N_KAO

         !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         ! First read Hamiltonian (hopping integrals) parametrization:
         write(ch_temp,'(a)') trim(adjustl(matter%Atoms(i)%Name))//'_'//trim(adjustl(matter%Atoms(j)%Name))//'_'
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_Hamiltonian_parameters.txt'
         inquire(file=trim(adjustl(File_name)),exist=file_exists)

         if (.not.file_exists) then ! try general name used for multiple species at once:
            if (numpar%verbose) write(*,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found. '
            write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), 'TB_Hamiltonian_parameters.txt'
            if (numpar%verbose) write(*,'(a)') 'Trying '//trim(adjustl(File_name))//' file instead.'
            inquire(file=trim(adjustl(File_name)),exist=file_exists)
         endif
         
         if (file_exists) then
            FN=106
            open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='READ')
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3421
            endif
            count_lines = 0
            INFO = 0

            ! Read first line to figure out which TB parametrization is used for this material:
            read(FN,*,IOSTAT=Reason) ch_temp
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3421
            endif
            if (j .GT. 1) then
               if (trim(adjustl(ch_temp)) .NE. trim(adjustl(TB_Hamil(1,1)%Param)) ) then
                  write(Error_descript,'(a)') 'Format of TB-Hamiltonian parameters "'//trim(adjustl(ch_temp))//'" does not coinside with parameters in the first file "'//trim(adjustl(TB_Hamil(1,1)%Param))//'". Inconsistent parameterization is not allowed.'
                  call Save_error_details(Err, 5, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3421
               endif
            endif

            ! Make the TB parameters of a selected class, depending on what is read in the file:
            select case (trim(adjustl(ch_temp)))
            case ('Pettifor')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_Pettifor::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Pettifor parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
            case ('Molteni')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_Molteni::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Molteni parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
            case ('Fu')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_Fu::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Fu parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
            case ('Mehl', 'NRL')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_NRL::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Mehl parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))

            case ('DFTB')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_DFTB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
               ! DFTB skf files contain parameters for both Hamiltonian and Repulsive potential, allocate both of them here:
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_DFTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))

            case ('DFTB_no_repulsion', 'DFTB_no_repulsive', 'DFTB_no_rep')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_DFTB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
               ! DFTB skf files does not contain Repulsive potential, allocate special case:
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_DFTB_no::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))

            case ('3TB')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_3TB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for 3TB parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))

               if (matter%N_KAO > 2) then
                  write(Error_descript,'(a,a,$)') '3TB-parametrization does not support more then binary compounds '
                  call Save_error_details(Err, 4, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3421
               endif

            case ('BOP')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_BOP::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for BOP parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
               ! BOP files contain parameters for both Hamiltonian and Repulsive potential, allocate both of them here:
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_BOP::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for BOP parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))

            case ('xTB', 'GFN', 'GFN0')
               if (.not.allocated(TB_Hamil)) then
                  allocate(TB_H_xTB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for xTB parametrization
                  TB_Hamil%Param = ''
               endif
               TB_Hamil(i,j)%Param = trim(adjustl(ch_temp))
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_xTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for xTB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))

            case default
               write(Error_descript,'(a,a,$)') 'Wrong TB-Hamiltonian parametrization class '//trim(adjustl(ch_temp))//' specified in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 4, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3421
            end select

            ! Prior to use TB parameters, we now always have to find out which class they belong to:
            select type (TB_Hamil)
            type is (TB_H_Pettifor)
               Error_descript = ''
               call read_Pettifor_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif

            type is (TB_H_Molteni)
               Error_descript = ''
               call read_Molteni_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif

             type is (TB_H_Fu)
               Error_descript = ''
               !call read_Pettifor_TB_Hamiltonian(FN, numpar%El_num_ij(i,j), TB_Hamil, Error_descript, INFO)
               call read_Fu_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_H_NRL)
               Error_descript = ''
               call read_Mehl_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)   ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_H_DFTB) !in this case, read both Hamiltonian and Repulsive parts together:
               Error_descript = ''
               select type (TB_Repuls)  ! to confirm that repulsive part is consistent with the Hamiltonian
               type is (TB_Rep_DFTB)   ! repulsive parameters provided in the skf-file
                  call read_DFTB_TB_Params(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO) ! below
               type is (TB_Rep_DFTB_no)   ! no repulsive parameters in skf-file
                  call read_DFTB_TB_Params_no_rep(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO) ! below
               endselect
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_H_3TB)
               Error_descript = ''
               call read_3TB_TB_Params(FN, i, j, TB_Hamil, numpar, matter, Error_descript, INFO) ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name))
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_H_BOP) !in this case, read both Hamiltonian and Repulsive parts together:
               Error_descript = ''
               select type (TB_Repuls)  ! to confirm that repulsive part is consistent with the Hamiltonian
               type is (TB_Rep_BOP)
                  call read_BOP_TB_Params(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO) ! below
               endselect
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
             type is (TB_H_xTB) !in this case, read both Hamiltonian and Repulsive parts together:
               Error_descript = ''
               select type (TB_Repuls)  ! to confirm that repulsive part is consistent with the Hamiltonian
               type is (TB_Rep_xTB)
                  call read_xTB_Params(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO) ! below
               endselect
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' called from file '//trim(adjustl(File_name))
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            end select
            close(FN)
         else
            write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, the program terminates'
            call Save_error_details(Err, 1, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3421
         endif


         !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         ! Now read repulsive part parametrization:
         write(ch_temp,'(a)') trim(adjustl(matter%Atoms(i)%Name))//'_'//trim(adjustl(matter%Atoms(j)%Name))//'_'
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_Repulsive_parameters.txt'
         inquire(file=trim(adjustl(File_name)),exist=file_exists)

         if (.not.file_exists) then ! try general name used for multiple species at once:

            if (numpar%verbose) write(*,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found. '
            write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), 'TB_Repulsive_parameters.txt'
            if (numpar%verbose) write(*,'(a)') 'Trying '//trim(adjustl(File_name))//' file instead.'
            inquire(file=trim(adjustl(File_name)),exist=file_exists)
         endif

         if (file_exists) then
            FN=107
            open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='READ')
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3421
            endif
            count_lines = 0
            INFO = 0

            ! Read first line to figure out which TB parametrization is used for this material:
            read(FN,*,IOSTAT=Reason) ch_temp
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3421
            endif

            if (j .GT. 1) then
               if (trim(adjustl(ch_temp)) .NE. trim(adjustl(TB_Repuls(1,1)%Param)) ) then
                  write(Error_descript,'(a,$)') 'Format of TB-repulseive parameters "', trim(adjustl(ch_temp)), '" does not coinside with parameters in the first file "'//trim(adjustl(TB_Repuls(1,1)%Param))//'" .Inconsistent parameterization is not allowed.'
                  call Save_error_details(Err, 5, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3421
               endif
            endif

            ! Make the TB parameters of a selected class, depending on what is read in the file:
            select case (trim(adjustl(ch_temp)))
            case ('Pettifor')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_Pettifor::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for Pettifor parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('Molteni')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_Molteni::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for Molteni parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
             case ('Fu')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_Fu::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for Fu parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('Mehl', 'NRL')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_NRL::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for NRL parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('DFTB')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_DFTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('DFTB_no_repulsion', 'DFTB_no_repulsive', 'DFTB_no_rep')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_DFTB_no::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('3TB')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_3TB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for 3TB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('BOP')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_BOP::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for BOP parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case ('xTB')
               if (.not.allocated(TB_Repuls)) then
                  allocate(TB_Rep_xTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for xTB parametrization
                  TB_Repuls%Param = ''
               endif
               TB_Repuls(i,j)%Param = trim(adjustl(ch_temp))
            case default
               write(Error_descript,'(a,a,a,$)') 'Wrong TB-repulsive parametrization class '//trim(adjustl(ch_temp))//' specified in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 4, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3421
            end select

            ! Prior to use TB parameters, we now always have to find out which class they belong to:
            select type (TB_Repuls)
            type is (TB_Rep_Pettifor)
               Error_descript = ''
               call read_Pettifor_TB_repulsive(FN, i,j, TB_Repuls, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_Rep_Molteni)
               Error_descript = ''
               call read_Molteni_TB_repulsive(FN, i,j, TB_Repuls, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
             type is (TB_Rep_Fu)
               Error_descript = ''
               call read_Fu_TB_repulsive(FN, i,j, TB_Repuls, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_Rep_NRL)
               Error_descript = ''
               ! There is no repulsive part in NRL
            type is (TB_Rep_DFTB)
               Error_descript = ''
               call read_DFTB_TB_repulsive(FN, i,j, TB_Repuls, Error_descript, INFO)    ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3421
               endif
            type is (TB_Rep_DFTB_no)
               Error_descript = ''
               ! Nothing to read, since repulsive potential is not provided in skf-file
            type is (TB_Rep_3TB)
               Error_descript = ''
               ! There is no repulsive part in 3TB
            type is (TB_Rep_BOP)
               Error_descript = ''
               ! Nothing to do with repulsive part in BOP parameterization
            type is (TB_Rep_xTB)
               Error_descript = ''
               !call read_xTB_repulsive(FN, i, j, TB_Repuls, Error_descript, INFO)    ! module "Dealing_with_xTB"
               ! Repulsive part has already been read above together with Hamiltonian parameters
            end select
            close(FN)
            !PAUSE 'READING INPUT'
         else
            write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, the program terminates'
            call Save_error_details(Err, 1, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3421
         endif
         
         
         !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         ! Now read van der Waals part parametrization, if it exists:
         write(ch_temp,'(a)') trim(adjustl(matter%Atoms(i)%Name))//'_'//trim(adjustl(matter%Atoms(j)%Name))//'_'
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_vdW.txt'
         inquire(file=trim(adjustl(File_name)),exist=file_exists)
         
         if (file_exists) then
            FN=108
            open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='READ')
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3422
            endif
            count_lines = 0
            INFO = 0

            ! Read the first line to figure out which TB vdW parametrization is used for this material:
            read(FN,*,IOSTAT=Reason) ch_temp
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3422
            endif
            
            ! Make the TB vdW parameters of a selected class, depending on what is read in the file:
            select case (trim(adjustl(ch_temp)))
            case ('Girifalco')
               if (.not.allocated(TB_Waals)) then
                  allocate(TB_vdW_Girifalco::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for Girifalco parametrization
                  ! Default values:
                  TB_Waals%Param = ''
                  TB_Waals%d0_cut = 0.0d0
                  TB_Waals%dd_cut = 0.0d0
               endif
               TB_Waals(i,j)%Param = trim(adjustl(ch_temp))

            case('LJ', 'Lennar-Jones', 'Lennard_Jones', 'lj', 'Lj', 'LENNARD-JONES', 'LENNARD_JONES', 'lennard-jones', 'lennard_jones')
               if (.not.allocated(TB_Waals)) then
                  allocate(TB_vdW_LJ_cut::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for LJ parametrization
                  ! Default values:
                  TB_Waals%Param = ''
                  TB_Waals%d0_cut = 0.0d0
                  TB_Waals%dd_cut = 0.0d0
                  select type (TB_Waals)
                  type is (TB_vdW_LJ_cut)
                     TB_Waals%eps = 0.0d0
                     TB_Waals%r0 = 0.0d0
                     TB_Waals%n = 1.0d0
                     TB_Waals%d0_short = 0.0d0
                     TB_Waals%dd_short = 0.0d0
                  end select
               endif
               TB_Waals(i,j)%Param = trim(adjustl(ch_temp))

            case('ILJ', 'ilj', 'Improved_LJ', 'Ilj', 'improved_lj', 'I_LJ', 'i_lj')
               if (.not.allocated(TB_Waals)) then
                  allocate(TB_vdW_ILJ_cut::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for LJ parametrization
                  ! Default values:
                  TB_Waals%Param = ''
                  TB_Waals%d0_cut = 0.0d0
                  TB_Waals%dd_cut = 0.0d0
                  select type (TB_Waals)
                  type is (TB_vdW_ILJ_cut)
                     TB_Waals%eps = 0.0d0
                     TB_Waals%r0 = 0.0d0
                     TB_Waals%n = 12.0d0
                     TB_Waals%m = 6.0d0
                     TB_Waals%d0_short = 0.0d0
                     TB_Waals%dd_short = 0.0d0
                  end select
               endif
               TB_Waals(i,j)%Param = trim(adjustl(ch_temp))

            case ('Dumitrica') ! UNFINISHED, DO NOT USE
               if (.not.allocated(TB_Waals)) then
                  allocate(TB_vdW_Dumitrica::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for Dumitrica parametrization
                  ! Default values:
                  TB_Waals%Param = ''
                  TB_Waals%d0_cut = 0.0d0
                  TB_Waals%dd_cut = 0.0d0
               endif
               TB_Waals(i,j)%Param = trim(adjustl(ch_temp))
            case default
               write(Error_descript,'(a,a,a,$)') 'Unknown TB-vdW parametrization class '// &
                     trim(adjustl(ch_temp))//' specified in file '//trim(adjustl(File_name))
               print*, trim(adjustl(Error_descript))
               print*, 'Proceeding without van der Waals forces'
               close(FN) ! close file
               goto 3422
            end select
            
            ! Read the parameters of the dispersion correction (vdW-type additional potential):
            select type (TB_Waals)
            type is (TB_vdW_Girifalco)
               Error_descript = ''
               call read_vdW_Girifalco_TB(FN, i,j, TB_Waals, Error_descript, INFO)  ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3422
               endif
               ! Set default values:
               TB_Waals(i,j)%d0_cut = TB_Waals(i,j)%d_cut   ! long-range cutoff radius [A]
               TB_Waals(i,j)%dd_cut = 0.0d0  ! long-range cutoff width [A]

            type is (TB_vdW_LJ_cut) ! Lennard-Jones (smoothly cut at short and large sitances)
               Error_descript = ''
               call read_vdW_LJ_TB(FN, i,j, TB_Waals, Error_descript, INFO)   ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name))
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3422
               endif
               if (numpar%verbose) write(*,'(a,f,f)') 'LJ: vdW minimum: ', TB_Waals(i,j)%eps, TB_Waals(i,j)%r0

            type is (TB_vdW_ILJ_cut) ! Improved Lennard-Jones (smoothly cut at short and large sitances)
               Error_descript = ''
               call read_vdW_ILJ_TB(FN, i,j, TB_Waals, Error_descript, INFO)   ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name))
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3422
               endif
               if (numpar%verbose) write(*,'(a,f,f)') 'ILJ: vdW minimum: ', TB_Waals(i,j)%eps, TB_Waals(i,j)%r0

            type is (TB_vdW_Dumitrica) ! UNFINISHED, DO NOT USE
               Error_descript = ''
               ! Set default values:
               TB_Waals(i,j)%d0_cut = 0.0d0   ! long-range cutoff radius [A]
               TB_Waals(i,j)%dd_cut = 0.0d0  ! long-range cutoff width [A]
               call read_vdW_Dumitrica_TB(FN, i,j, TB_Waals, Error_descript, INFO)  ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3422
               endif
            end select
            close(FN)
         else
            if (numpar%verbose) print*, 'No van der Waals file found, go on without van der Waals forces'
         endif !(file_exists)
3422     continue


         !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         ! Now read Coulomb parameterization:
         write(ch_temp,'(a)') trim(adjustl(matter%Atoms(i)%Name))//'_'//trim(adjustl(matter%Atoms(j)%Name))//'_'
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_Coulomb.txt'
         inquire(file=trim(adjustl(File_name)),exist=file_exists)
         
         if (file_exists .and. (numpar%E_work < 1d10)) then ! there can be unballanced charge, try to use Coulomb potential
            FN=109
            open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='READ')
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3423
            endif
            count_lines = 0
            INFO = 0
            ! Read the first line to figure out which TB vdW parametrization is used for this material:
            read(FN,*,IOSTAT=Reason) ch_temp
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3423
            endif
            
            ! Make the Coulomb parameters of a selected class, depending on what is read in the file:
            select case (trim(adjustl(ch_temp)))
            case ('Coulomb_cut')
               if (.not.allocated(TB_Coul)) then
                  allocate(TB_Coulomb_cut::TB_Coul(matter%N_KAO,matter%N_KAO)) ! make it for Coulomb parametrization
                  TB_Coul%Param = ''
               endif
               TB_Coul(i,j)%Param = trim(adjustl(ch_temp))
            case ('Cutie') ! testing ONLY
               if (.not.allocated(TB_Coul)) then
                  allocate(Cutie::TB_Coul(matter%N_KAO,matter%N_KAO)) ! make it for Coulomb parametrization
                  TB_Coul%Param = ''                  
               endif
               TB_Coul(i,j)%Param = trim(adjustl(ch_temp))
            case default
               write(Error_descript,'(a,a,a,$)') 'Unknown Coulomb parametrization class '// &
                  trim(adjustl(ch_temp))//' specified in file '//trim(adjustl(File_name))
               print*, trim(adjustl(Error_descript))
               print*, 'Proceeding without Coulomb forces from unballanced charge'
               close(FN) ! close file
               goto 3423
            end select
            
            select type (TB_Coul)
            type is (TB_Coulomb_cut)
               Error_descript = ''
               !call read_Pettifor_TB_repulsive(FN, numpar%El_num_ij(i,j), TB_Repuls, Error_descript, INFO)
               call read_Coulomb_cut_TB(FN, i,j, TB_Coul, Error_descript, INFO)
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3423
               endif
               ! [N*A/e] coupling constant of Coulomb field (e/(4*Pi*e0)):
               TB_Coul(i,j)%k = g_e/(4.0d0*g_Pi*g_e0)*1.0d10
               !print*, TB_Coul(i,j)%Param, TB_Coul(i,j)%k
            end select
            close(FN)
         else
            if (numpar%verbose) then
               print*, 'No Coulomb parameterization file found, or no unbalanced charge possible'
               print*, 'go on without Coulomb forces'
            endif
         endif
3423     continue
         

         !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
         ! Now read additional short-range repulsive (exponential wall) parameterization:
         write(ch_temp,'(a)') trim(adjustl(matter%Atoms(i)%Name))//'_'//trim(adjustl(matter%Atoms(j)%Name))//'_'
         !write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_wall.txt'
         write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), &
                                     trim(adjustl(ch_temp))//trim(adjustl(m_wall_pot))
         inquire(file=trim(adjustl(File_name)),exist=file_exists)
         if (.not.file_exists) then ! try the new name
            !write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_short.txt'
            write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), &
                                        trim(adjustl(ch_temp))//trim(adjustl(m_short_pot))
            inquire(file=trim(adjustl(File_name)),exist=file_exists)
         endif
         ! try the other order of elements:
         if (.not.file_exists) then
            write(ch_temp,'(a)') trim(adjustl(matter%Atoms(j)%Name))//'_'//trim(adjustl(matter%Atoms(i)%Name))//'_'
            !write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_wall.txt'
            write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), &
                                     trim(adjustl(ch_temp))//trim(adjustl(m_wall_pot))
            inquire(file=trim(adjustl(File_name)),exist=file_exists)
         endif
         if (.not.file_exists) then ! try the new name
            !write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), trim(adjustl(ch_temp))//'TB_short.txt'
            write(File_name, '(a,a,a)') trim(adjustl(Path)), trim(adjustl(numpar%path_sep)), &
                                        trim(adjustl(ch_temp))//trim(adjustl(m_short_pot))
            inquire(file=trim(adjustl(File_name)),exist=file_exists)
         endif

         if (file_exists) then
            FN=110
            open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='READ')
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3425
            endif
            count_lines = 0
            INFO = 0
            ! Read the first line to figure out which exponential wall parametrization is used for this material:
            read(FN,*,IOSTAT=Reason) ch_temp
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3425
            endif
            !print*, 'test 0:', trim(adjustl(ch_temp))

            ! Check if there is a path to another file to be read:
            select case (trim(adjustl(ch_temp)))
            case ('PATH', 'Path', 'path')
               read(FN,'(a)',IOSTAT=Reason) ch_temp
               call ensure_correct_path_separator(ch_temp, numpar%path_sep)  ! module "Dealing_with_files"
               inquire(file=trim(adjustl(ch_temp)),exist=file_exists)
               !print*, 'test 1:', trim(adjustl(ch_temp)), file_exists

               if (.not.file_exists) then
                  write(Error_descript,'(a,a,a,$)') 'Path in short-range (exponential wall) file '// &
                     trim(adjustl(File_name))//' not found: '//trim(adjustl(ch_temp))
                  print*, trim(adjustl(Error_descript))
                  print*, 'Proceeding without additional short-range (exponential wall) forces'
                  close(FN) ! close file
                  goto 3425
               else
                  close(FN) ! close file, to open the other one
                  write(File_name, '(a,a,a)') trim(adjustl(ch_temp))
                  open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='READ')
                  inquire(file=trim(adjustl(File_name)),opened=file_opened)

                  if (.not. file_opened) then
                     write(Error_descript,'(a,a,a,$)') 'File in short-range (exponential wall) file '// &
                     trim(adjustl(ch_temp))//' could not be opened: '//trim(adjustl(File_name))
                     print*, trim(adjustl(Error_descript))
                     print*, 'Proceeding without additional short-range (exponential wall) forces'
                     close(FN) ! close file
                     goto 3425
                  else
                     read(FN,*,IOSTAT=Reason) ch_temp
                     call read_file(Reason, count_lines, read_well)
                     if (.not. read_well) then
                        write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
                        call Save_error_details(Err, 3, Error_descript)
                        print*, trim(adjustl(Error_descript))
                        goto 3425
                     endif ! (.not. read_well)
                  endif ! (.not. file_opened)
               endif ! (.not.file_exists)
            end select ! (trim(adjustl(ch_temp)))

            ! Make the exponential wall  parameters of a selected class, depending on what is read in the file:
            select case (trim(adjustl(ch_temp)))
            case ('Simple_wall', 'SIMPLE_WALL', 'simple_wall')
               if (.not.allocated(TB_Expwall)) then
                  allocate(TB_Exp_wall_simple::TB_Expwall(matter%N_KAO,matter%N_KAO)) ! make it for exponential wall  parametrization
                  TB_Expwall%Param = ''
               endif
               TB_Expwall(i,j)%Param = trim(adjustl(ch_temp))
               !print*, i, j, TB_Expwall(i,j)%Param, trim(adjustl(ch_temp))
            case ('General', 'general', 'GENERAL')
               if (.not.allocated(TB_Expwall)) then
                  allocate(TB_Short_Rep::TB_Expwall(matter%N_KAO,matter%N_KAO)) ! make it for exponential wall  parametrization
                  TB_Expwall%Param = ''
               endif
               TB_Expwall(i,j)%Param = trim(adjustl(ch_temp))
               !print*, i, j, TB_Expwall(i,j)%Param, trim(adjustl(ch_temp))
            case default
               write(Error_descript,'(a,a,a,$)') 'Unknown short-range (exponential wall) parametrization class '//trim(adjustl(ch_temp))//' specified in file '//trim(adjustl(File_name))
               print*, trim(adjustl(Error_descript))
               print*, 'Proceeding without additional short-range (exponential wall) forces'
               close(FN) ! close file
               goto 3425
            end select
            
            ! Prior to use Exponential wall parameters, we now always have to find out which class the belong to:
            select type (TB_Expwall)
            type is (TB_Exp_wall_simple)
               Error_descript = ''
               call read_Exponential_wall_TB(FN, i,j, TB_Expwall, Error_descript, INFO)   ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name)) 
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3425
               endif
            type is (TB_Short_Rep)
               Error_descript = ''
               call read_Short_Rep_TB(FN, i,j, TB_Expwall, Error_descript, INFO)   ! below
               if (INFO .NE. 0) then
                  Err%Err_descript = trim(adjustl(Error_descript))//' in file '//trim(adjustl(File_name))
                  call Save_error_details(Err, INFO, Err%Err_descript)
                  print*, trim(adjustl(Err%Err_descript))
                  goto 3425
               endif
            end select
            close(FN)
         else
            if (numpar%verbose) then
               print*, 'No exponential wall parameterization file found'
               print*, 'go on without an exponential wall at short distances'
            endif
         endif
3425     continue

      enddo do_second
   enddo do_first

   ! For 3TB parameterization, compounds requires rewriting on-site terms:
   select type (TB_Hamil)
   type is (TB_H_3TB)   ! UNFINISHED
      if (matter%N_KAO > 1) then
         ! [OS 0] :
!          TB_Hamil(1,1)%Hhavg = TB_Hamil(1,2)%Hhavg
!          TB_Hamil(1,1)%Hhcf = TB_Hamil(1,2)%Hhcf
!          TB_Hamil(2,2)%Hhavg = TB_Hamil(2,1)%Hhavg
!          TB_Hamil(2,2)%Hhcf = TB_Hamil(2,1)%Hhcf
         ! [OS 1] :
         TB_Hamil(1,1)%Hhavg = TB_Hamil(2,1)%Hhavg
         TB_Hamil(1,1)%Hhcf = TB_Hamil(2,1)%Hhcf
         TB_Hamil(2,2)%Hhavg = TB_Hamil(1,2)%Hhavg
         TB_Hamil(2,2)%Hhcf = TB_Hamil(1,2)%Hhcf
      endif
   endselect


3421 continue
end subroutine read_TB_parameters



subroutine read_Short_Rep_TB(FN, i,j, TB_Expwall, Error_descript, INFO)   ! below
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Short_Rep), dimension(:,:), intent(inout) ::  TB_Expwall ! parameters of the exponential wall potential
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------
   integer :: count_lines, Reason
   logical :: read_well
   character(30) :: text
   count_lines = 1
   INFO = 0

   ! Set default values:
   TB_Expwall(i,j)%f_exp%use_it = .false.       ! no exp by default
   TB_Expwall(i,j)%f_inv_exp%use_it = .false.   ! no inverse exp by default
   TB_Expwall(i,j)%f_ZBL%use_it = .false.       ! no ZBL potential by default
   TB_Expwall(i,j)%f_cut_inv%use_it = .false.   ! no short-range cutoff by default
   TB_Expwall(i,j)%f_cut%d0 = 0.0d0 ! cut off at zero, no repulsion by default
   TB_Expwall(i,j)%f_cut%dd = 0.01d0 ! short cut-off by default
   TB_Expwall(i,j)%f_tab%use_it = .false.       ! no tabulated potential by default

   read_well = .true.   ! to start with
   RD: do while (read_well)
      read(FN,*,IOSTAT=Reason) text
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) exit RD  ! end of file, stop reading

      call interpret_short_range_data(FN, count_lines, read_well, text, TB_Expwall(i,j), INFO, Error_descript) ! below
      if (.not. read_well) exit RD  ! end of file, stop reading
!       print*, 'o:', i, j, text, TB_Expwall(i,j)%f_inv_exp%use_it
   enddo RD
!
!    print*, i, j, TB_Expwall(:,:)%f_inv_exp%use_it
!    pause 'read_Short_Rep_TB'
end subroutine read_Short_Rep_TB


subroutine interpret_short_range_data(FN, count_lines, read_well, text, TB_Expwall, INFO, Error_descript)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(inout) :: count_lines
   logical, intent(inout) :: read_well
   character(*), intent(in) :: text
   type(TB_Short_Rep), intent(inout) ::  TB_Expwall ! parameters of the exponential wall potential
   integer, intent(inout) :: INFO   ! error description
   character(*), intent(inout) :: Error_descript   ! error save
   !-------------------------------------
   integer :: N_pow, i, Reason

   select case (trim(adjustl(text)))
   case ('CUTOFF', 'Cutoff', 'cutoff', 'CUT_OFF', 'Cut_off', 'cut_off', 'CUT-OFF', 'Cut-off', 'cut-off', 'FERMI', 'Fermi', 'fermi')
      read(FN,*,IOSTAT=Reason) TB_Expwall%f_cut%d0, TB_Expwall%f_cut%dd
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         return   ! exit the function if there is nothing else to do
      endif

   case ('CUTOFF_SHORT', 'Cutoff_short', 'cutoff_short', 'CUT_OFF_SHORT', 'Cut_off_short', 'cut_off_short', &
         'CUT-OFF_short', 'Cut-off_short', 'cut-off_short', 'FERMI_INV', 'Fermi_inv', 'fermi_inv')
      read(FN,*,IOSTAT=Reason) TB_Expwall%f_cut_inv%d0, TB_Expwall%f_cut_inv%dd
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         return   ! exit the function if there is nothing else to do
      endif

   case ('EXP', 'Exp', 'exp', 'Exponential', 'exponential', 'EXPONENTIAL')
      TB_Expwall%f_exp%use_it = .true.
      read(FN,*,IOSTAT=Reason) TB_Expwall%f_exp%Phi, TB_Expwall%f_exp%r0, TB_Expwall%f_exp%a
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         TB_Expwall%f_exp%use_it = .false. ! not to use
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         return   ! exit the function if there is nothing else to do
      endif

   case ('INVEXP', 'InvExp', 'invexp', 'Invexp', 'INV_EXP', 'Inv_Exp', 'Inv_exp', 'inv_exp')
      TB_Expwall%f_inv_exp%use_it = .true.
      read(FN,*,IOSTAT=Reason) TB_Expwall%f_inv_exp%C, TB_Expwall%f_inv_exp%r0
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         TB_Expwall%f_inv_exp%use_it = .false.   ! not to use
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         return   ! exit the function if there is nothing else to do
      endif

   case ('POWER', 'Power', 'power', 'POW', 'Pow', 'pow')
      read(FN,*,IOSTAT=Reason) N_pow   ! number of power-functions
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         allocate(TB_Expwall%f_pow(N_pow))
         TB_Expwall%f_pow(:)%use_it = .true.
         do i = 1, N_pow   ! read for all functions
            read(FN,*,IOSTAT=Reason) TB_Expwall%f_pow(i)%Phi, TB_Expwall%f_pow(i)%r0, TB_Expwall%f_pow(i)%m
!           print*, TB_Expwall%f_pow%Phi, TB_Expwall%f_pow%r0, TB_Expwall%f_pow%m
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               deallocate(TB_Expwall%f_pow)   ! not to use
               write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
               INFO = 3
               return   ! exit the function if there is nothing else to do
            endif
         enddo
      else
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         return   ! exit the function if there is nothing else to do
      endif ! (read_well)

   case ('ZBL', 'zbl')
      TB_Expwall%f_ZBL%use_it = .true.

   case ('TAB', 'Tab', 'tab', 'TABLE', 'Table', 'table', 'TABULATED', 'Tabulated', 'tabulated')
      call Process_tabulated_potential(FN, TB_Expwall, count_lines, read_well, INFO, Error_descript)   ! below
      if (INFO /= 0) return

   end select
end subroutine interpret_short_range_data



subroutine Process_tabulated_potential(FN, TB_Expwall, count_lines, read_well, INFO, Error_descript)
   integer, intent(in) :: FN  ! file number to read from
   type(TB_Short_Rep), intent(inout) ::  TB_Expwall ! parameters of the exponential wall potential
   integer, intent(inout) :: count_lines
   logical, intent(inout) :: read_well
   integer, intent(inout) :: INFO   ! error description
   character(*), intent(inout) :: Error_descript   ! error save
   !-------------------------------------
   integer :: N_pow, i, Reason
   character(100) :: text_line, ch_var

   INFO = 0 ! no error to start with
   read_well = .true.   ! to start with
   TB_Expwall%f_tab%use_spline = .false.   ! by default, use finite difference

   read(FN, '(a)', IOSTAT=Reason) text_line
   call read_file(Reason, count_lines, read_well)
   read(text_line,*,IOSTAT=Reason) N_pow, ch_var    ! number of power-functions; type of
   if (Reason == 0) then   ! read the marker of potential
      select case (trim(adjustl(ch_var)))
      case ('diff', 'difference', 'findif', 'fin_dif')
         TB_Expwall%f_tab%use_spline = .false.   ! use finite difference
      case ('spline', 'SPLINE')
         TB_Expwall%f_tab%use_spline = .true.    ! use spline
      case default   ! use spline
         TB_Expwall%f_tab%use_spline = .false.   ! be default, use finite difference
      end select
   else  ! try to read just one number
      read(text_line,*,IOSTAT=Reason) N_pow  ! number of power-functions
      count_lines = count_lines - 1 ! interpreting the same line again
      call read_file(Reason, count_lines, read_well)
   endif

   if (read_well) then
      TB_Expwall%f_tab%use_it = .true.
      allocate(TB_Expwall%f_tab%R(N_pow), source = 0.0d0)
      allocate(TB_Expwall%f_tab%E(N_pow), source = 0.0d0)
      do i = 1, N_pow   ! read for all functions
         read(FN,*,IOSTAT=Reason) TB_Expwall%f_tab%R(i), TB_Expwall%f_tab%E(i)
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            TB_Expwall%f_tab%use_it = .false.
            deallocate(TB_Expwall%f_tab%R, TB_Expwall%f_tab%E)   ! not to use
            write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
            INFO = 3
            return   ! exit the function if there is nothing else to do
         endif
      enddo

      ! Create cubic splines:
      !call make_natural_cubic_splines(TB_Expwall%f_tab%R, TB_Expwall%f_tab%E, &
      !      TB_Expwall%f_tab%a, TB_Expwall%f_tab%b, TB_Expwall%f_tab%c, TB_Expwall%f_tab%d)  ! module "Algebra_tools"
      call make_cubic_splines(TB_Expwall%f_tab%R, TB_Expwall%f_tab%E, &
            TB_Expwall%f_tab%a, TB_Expwall%f_tab%b, TB_Expwall%f_tab%c, TB_Expwall%f_tab%d)  ! module "Algebra_tools"

      ! Test spline:
!       do i = 1, N_pow-1
!          write(*,'(a, i4, f, es, es, es, es)') 'Spline:', i, TB_Expwall%f_tab%R(i), TB_Expwall%f_tab%E(i), cubic_function(0.0d0,TB_Expwall%f_tab%a(i), TB_Expwall%f_tab%b(i), TB_Expwall%f_tab%c(i), TB_Expwall%f_tab%d(i)), &
!          cubic_function((TB_Expwall%f_tab%R(i+1)-TB_Expwall%f_tab%R(i))*0.5d0,TB_Expwall%f_tab%a(i), TB_Expwall%f_tab%b(i), TB_Expwall%f_tab%c(i), TB_Expwall%f_tab%d(i)), &
!          cubic_function((TB_Expwall%f_tab%R(i+1)-TB_Expwall%f_tab%R(i)),TB_Expwall%f_tab%a(i), TB_Expwall%f_tab%b(i), TB_Expwall%f_tab%c(i), TB_Expwall%f_tab%d(i))
!       enddo
!       pause 'Process_tabulated_potential'

   else
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      return   ! exit the function if there is nothing else to do
   endif
end subroutine Process_tabulated_potential




subroutine read_Exponential_wall_TB(FN, i,j, TB_Expwall, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Exp_wall_simple), dimension(:,:), intent(inout) ::  TB_Expwall ! parameters of the exponential wall potential
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   read(FN,*,IOSTAT=Reason) TB_Expwall(i,j)%C
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3426
   endif

   read(FN,*,IOSTAT=Reason) TB_Expwall(i,j)%r0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3426
   endif

   read(FN,*,IOSTAT=Reason) TB_Expwall(i,j)%d0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3426
   endif

   read(FN,*,IOSTAT=Reason) TB_Expwall(i,j)%dd
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3426
   endif
   
3426 continue
end subroutine read_Exponential_wall_TB




subroutine read_Coulomb_cut_TB(FN, i,j, TB_Coul, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Coulomb_cut), dimension(:,:), intent(inout) ::  TB_Coul ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   read(FN,*,IOSTAT=Reason) TB_Coul(i,j)%dm
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3424
   endif

   read(FN,*,IOSTAT=Reason) TB_Coul(i,j)%dd
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3424
   endif

3424 continue
end subroutine read_Coulomb_cut_TB



subroutine read_vdW_Girifalco_TB(FN, i,j, TB_Waals, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_vdW_Girifalco), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%C12
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%C6
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
!    read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%r_L
!    call read_file(Reason, count_lines, read_well)
!    if (.not. read_well) then
!       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
!       INFO = 3
!       goto 3423
!    endif
!    
!    read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d_L
!    call read_file(Reason, count_lines, read_well)
!    if (.not. read_well) then
!       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
!       INFO = 3
!       goto 3423
!    endif
! 
!    read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%r_S
!    call read_file(Reason, count_lines, read_well)
!    if (.not. read_well) then
!       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
!       INFO = 3
!       goto 3423
!    endif
!    
!    read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d_S
!    call read_file(Reason, count_lines, read_well)
!    if (.not. read_well) then
!       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
!       INFO = 3
!       goto 3423
!    endif
! 
!    read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%r_LJ
!    call read_file(Reason, count_lines, read_well)
!    if (.not. read_well) then
!       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
!       INFO = 3
!       goto 3423
!    endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%dm
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d_cut
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%a
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%b
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%c
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%dsm
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%ds_cut
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%as
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%bs
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%cs
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%ds
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%es
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%fs
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
!    print*, TB_Waals(i,j)
!    pause

3423 continue
end subroutine read_vdW_Girifalco_TB


subroutine read_vdW_LJ_TB(FN, i,j, TB_Waals, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_vdW_LJ_cut), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO  ! error description
   !---------------------------
   character(20) :: LJ_type
   character(200) :: read_line
   real(8) :: A, B, n
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Form of LJ:
   read(FN,*,IOSTAT=Reason) LJ_type
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

   ! LJ coefficients:
   read(FN,'(a)',IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) A, B, n   ! try to read it from the text line
   if (Reason /= 0) then   ! try to read two coefficients
      read(read_line,*,IOSTAT=Reason) A, B   ! try to read it into only 2 variables
      n = 6.0d0 ! default value
   endif
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

   ! Various forms of LJ potential are supported
   ! https://en.wikipedia.org/wiki/Lennard-Jones_potential
   select case (trim(adjustl(LJ_type)))
   case ('SE', 'Se', 'se', 'sE', 'Segma-E', 'sigma-e', 'SIGMA-E') ! Sigma-epsylon form
      TB_Waals(i,j)%eps = A/8.0d0         ! [eV] prefactor
      TB_Waals(i,j)%r0  = sqrt(2.0d0)*B   ! [A] radius
   case ('AB', 'Ab', 'ab', 'aB') ! AB form
      TB_Waals(i,j)%eps = B**2/(4.0d0*A)        ! [eV] prefactor
      TB_Waals(i,j)%r0  = (A/B)**(1.0d0/6.0d0)  ! [A] radius
   case ('n-exp', 'N-exp', 'N-Exp', 'N-EXP') ! n-exp form
      TB_Waals(i,j)%eps = A   ! [eV] prefactor
      TB_Waals(i,j)%r0  = B   ! [A] radius
   end select
   TB_Waals(i,j)%n   = n   ! power

   ! Short-range cutoff parameters:
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d0_short, TB_Waals(i,j)%dd_short ! [A] cutoff radiues, [A] cutoff width
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

   ! Long-range cutoff parameters:
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d0_cut, TB_Waals(i,j)%dd_cut ! [A] cutoff radiues, [A] cutoff width
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

3430 continue
end subroutine read_vdW_LJ_TB



subroutine read_vdW_ILJ_TB(FN, i,j, TB_Waals, Error_descript, INFO)
   ! Improved Lennard-Jones: https://www.mdpi.com/1420-3049/26/13/3906
   ! V=eps*( m/(n-m)*(r0/r)^(n) - n/(n-m)*(r0/r)^m )
   ! reducing to LJ for n=12, m=6
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_vdW_ILJ_cut), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO  ! error description
   !---------------------------
   character(20) :: LJ_type
   character(200) :: read_line
   real(8) :: A, B, n, m
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! ILJ coefficients:
   read(FN,'(a)',IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) A, B, n, m   ! try to read it from the text line
   if (Reason /= 0) then   ! try to read two coefficients
      read(read_line,*,IOSTAT=Reason) A, B, n   ! try to read it into only 3 variables
      if (Reason /= 0) then   ! try to read two coefficients
         read(read_line,*,IOSTAT=Reason) A, B   ! try to read it into only 2 variables
         n = 12.0d0  ! default value to reduce to standard LJ
         m = 6.0d0   ! default value to reduce to standard LJ
      else
         m = 6.0d0   ! default value to reduce to standard LJ
      endif
   endif
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

   TB_Waals(i,j)%eps = A   ! [eV] prefactor
   TB_Waals(i,j)%r0  = B   ! [A] radius
   TB_Waals(i,j)%n   = n   ! power of first term
   TB_Waals(i,j)%m   = m   ! power of second term

   ! Short-range cutoff parameters:
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d0_short, TB_Waals(i,j)%dd_short ! [A] cutoff radiues, [A] cutoff width
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

   ! Long-range cutoff parameters:
   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%d0_cut, TB_Waals(i,j)%dd_cut ! [A] cutoff radiues, [A] cutoff width
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3430
   endif

3430 continue
end subroutine read_vdW_ILJ_TB


subroutine read_vdW_Dumitrica_TB(FN, i,j, TB_Waals, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_vdW_Dumitrica), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%C6
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Waals(i,j)%alpha
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   
3423 continue
end subroutine read_vdW_Dumitrica_TB



subroutine read_Molteni_TB_repulsive(FN, i,j, TB_Repuls, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Rep_Molteni), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Read the repulsive parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%NP
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   select case (TB_Repuls(i,j)%NP)
   case (1:2) ! Molteni
      call read_Molteni_repuls(FN, TB_Repuls, i, j, read_well, count_lines, Reason, Error_descript, INFO)
   case default ! Allen
      call read_Allen_repuls(FN, TB_Repuls, i, j, read_well, count_lines, Reason, Error_descript, INFO)
   end select


3423 continue
end subroutine read_Molteni_TB_repulsive


subroutine read_Allen_repuls(FN, TB_Repuls, i, j, read_well, count_lines, Reason, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Rep_Molteni), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   logical, intent(inout) :: read_well
   integer, intent(inout) :: count_lines, Reason
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%b
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

  read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%r0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%rcut
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   TB_Repuls(i,j)%rcut = TB_Repuls(i,j)%rcut*TB_Repuls(i,j)%r0

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%d
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   3423 continue
end subroutine read_Allen_repuls



subroutine read_Molteni_repuls(FN, TB_Repuls, i, j, read_well, count_lines, Reason, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Rep_Molteni), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   logical, intent(inout) :: read_well
   integer, intent(inout) :: count_lines, Reason
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%phi1
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%phi2
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%r0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%rcut
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
   TB_Repuls(i,j)%rcut = TB_Repuls(i,j)%rcut*TB_Repuls(i,j)%r0

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%d
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   select case (TB_Repuls(i,j)%NP)
   case (2) ! rational:
      read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%m
   case default ! exp:
      read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%alpha
   end select
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   3423 continue
end subroutine read_Molteni_repuls


subroutine read_Pettifor_TB_repulsive(FN, i, j, TB_Repuls, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Rep_Pettifor), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Read the repulsive parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%E0_TB
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%phi0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%m
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%mc
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%d0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%d1
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%dm
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%dc
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(5)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

3423 continue
end subroutine read_Pettifor_TB_repulsive



subroutine read_Fu_TB_repulsive(FN, i, j, TB_Repuls, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Rep_Fu), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Read the repulsive parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%E0_TB
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%phi0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%m
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%mc
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%d0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%d1
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%dm
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%dc
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%c0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%a0(5)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%C_a
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3423
   endif
3423 continue
end subroutine read_Fu_TB_repulsive


subroutine read_DFTB_TB_repulsive(FN, i,j, TB_Repuls, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j ! numbers of pair of elements for which we read the data
   type(TB_Rep_DFTB), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Skip first line as it are already defined within the Hamiltonian file
   read(FN,*,IOSTAT=Reason)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3415
   endif

   read(FN,*,IOSTAT=Reason) TB_Repuls(i,j)%ToP   ! type of parameterization: 0=polinomial, 1=spline
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3415
   endif
   
3415 continue
end subroutine read_DFTB_TB_repulsive


!------------------------------------------------------------
subroutine read_Molteni_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_Molteni), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Read the Hamiltonian parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Es
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Ep
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Esa
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(5)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%r0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%n
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rcut
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif
   TB_Hamil(i,j)%rcut = TB_Hamil(i,j)%rcut*TB_Hamil(i,j)%r0

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%d
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif
3422 continue
end subroutine read_Molteni_TB_Hamiltonian


subroutine read_Pettifor_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_Pettifor), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Read the Hamiltonian parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Es
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Ep
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%r0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%n
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%r1
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rm 
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

3422 continue
end subroutine read_Pettifor_TB_Hamiltonian



subroutine read_Fu_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_Fu), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason
   logical read_well
   count_lines = 1

   ! Read the Hamiltonian parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Es
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Ep
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%V0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%r0
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%n
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%r1
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rm 
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nc(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(1)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(2)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(3)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c0(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c1(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c2(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%c3(4)
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%C_a
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3422
   endif

3422 continue

end subroutine read_Fu_TB_Hamiltonian




subroutine read_Mehl_TB_Hamiltonian(FN, i,j, TB_Hamil, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_NRL), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   integer count_lines, Reason, i_cur, ind
   logical read_well
   count_lines = 1

   ! Read the Hamiltonian parameters of TB:
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%ind_split, TB_Hamil(i,j)%ind_overlap, ind
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3425
   endif
   
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%Rc, TB_Hamil(i,j)%lden
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3425
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%lambd
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
      INFO = 3
      goto 3425
   endif
   
   
   do i_cur = 1, 3	! s, p or d states
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%al(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%bl(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%cl(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%dl(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   enddo
   
   TB_Hamil(i,j)%al(4) = 0.0d0
   TB_Hamil(i,j)%bl(4) = 0.0d0
   TB_Hamil(i,j)%cl(4) = 0.0d0
   TB_Hamil(i,j)%dl(4) = 0.0d0
   select case (TB_Hamil(i,j)%ind_split)	! in case ther eis splitting of d states into t2g and e2:
   case (1)	! there is splitting
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%al(4)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%bl(4)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%cl(4)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%dl(4)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   case default
      ! There is no splitting
   end select
   
   
   do i_cur = 1, 10	! all states (ss sigma) thru (dd delta)
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%ellm(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%fllm(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      if (ind <= 0) then
         TB_Hamil(i,j)%gllm(i_cur) = 0.0d0
      else 	! the variable is not zero
         read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%gllm(i_cur)
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
            INFO = 3
            goto 3425
         endif
      endif
      
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%hllm(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   enddo
   
   
   do i_cur = 1, 10	! all states (ss sigma) thru (dd delta)
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%pllm(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%qllm(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   
      if (ind <= 0) then
         TB_Hamil(i,j)%rllm(i_cur) = 0.0d0
      else 	! the variable is not zero
         read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rllm(i_cur)
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
            INFO = 3
            goto 3425
         endif
      endif
      
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%sllm(i_cur)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3425
      endif
   enddo
   
!    
!    print*, TB_Hamil(i,j)
!    
!    PAUSE 'read_Mehl_TB_Hamiltonian'
   
   3425 continue
end subroutine read_Mehl_TB_Hamiltonian



subroutine read_DFTB_TB_Params(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_DFTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_DFTB), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter	! all material parameters
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------------------------------------
   character(200) :: Folder_name, File_name
   character(200) :: path_to_skf
   integer count_lines, Reason, i_cur, ind, FN_skf, ToA, N_basis_siz
   logical file_exist, file_opened, read_well
   INFO = 0
   count_lines = 2
   
   read(FN,*,IOSTAT=Reason) path_to_skf
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3426
   endif

   ! Define the path to the skf file:
   select case (trim(adjustl(path_to_skf)))  ! how to set it
   case ('PATH', 'Path', 'path') ! then specify exactly the path to the skf file
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%param_name  ! read the full path
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3426
      endif
      Folder_name = trim(adjustl(TB_Hamil(i,j)%param_name))    ! folder with chosen parameters sets
   case default   ! it is a parameterization within the predefined directory 'DFTB'
      TB_Hamil(i,j)%param_name = trim(adjustl(path_to_skf)) ! name of the directory with skf files
      Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_DFTB_directory))//numpar%path_sep ! default folder
      Folder_name = trim(adjustl(Folder_name))//trim(adjustl(TB_Hamil(i,j)%param_name))    ! folder with chosen parameters sets
   endselect

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rcut, TB_Hamil(i,j)%d  ! [A] cut off, and width of cut-off region [A]
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3426
   endif

   ! Make sure the slash is correct:
   call ensure_correct_path_separator(Folder_name, numpar%path_sep)  ! module "Dealing_with_files"

   ! Assume it is a file name:
   inquire(file=trim(adjustl(Folder_name)),exist=file_exist)
   if (file_exist) then ! such a file exists, use it
      File_name = trim(adjustl(Folder_name))
   else  ! no file => assume it is a directory name, and use default file name:
      ! Construct name of the skf file:
      call construct_skf_filename( trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                                 File_name)    ! module "Dealing_with_DFTB"
      File_name = trim(adjustl(Folder_name))//numpar%path_sep//trim(adjustl(File_name))
   endif

   ! Check if such DFTB parameterization exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exist)
   if (.not.file_exist) then
      Error_descript = 'File '//trim(adjustl(File_name))//' not found, the program terminates'
      INFO = 1
      goto 3426
   endif
   FN_skf=111
   open(UNIT=FN_skf, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      INFO = 2
      goto 3426
   endif
   
   ToA = same_or_different_atom_types(trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name))) ! module "Dealing_with_DFTB"
   call read_skf_file(FN_skf, TB_Hamil(i,j), TB_Repuls(i,j), ToA, Error_descript)    ! module "Dealing_with_DFTB"
   if (LEN(trim(adjustl(Error_descript))) > 0) then
      INFO = 5
      goto 3426
   endif
   
   ! Check which basis set is used: 0=s, 1=sp3, 2=sp3d5:
   if ((i == matter%N_KAO) .and. (j == matter%N_KAO)) then  ! only when all parameters for all elements are read from files:
      call idnetify_basis_size(TB_Hamil, N_basis_siz)  ! module "Dealing_with_DFTB"'
      numpar%N_basis_size = max(numpar%N_basis_size,N_basis_siz)
   endif
   
3426 continue 
   ! Close files that have been read through:
   call close_file('close', FN=FN_skf) ! module "Dealing_with_files"
   call close_file('close', FN=FN) ! module "Dealing_with_files"
end subroutine read_DFTB_TB_Params



subroutine read_DFTB_TB_Params_no_rep(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_DFTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_DFTB_no), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter	! all material parameters
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------------------------------------
   character(200) :: Folder_name, File_name, Inner_folder_name, path_to_skf
   integer count_lines, Reason, i_cur, ind, FN_skf, ToA, N_basis_siz
   logical file_exist, file_opened, read_well
   INFO = 0
   count_lines = 2

   read(FN,*,IOSTAT=Reason) path_to_skf    ! name of the directory with skf files
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3426
   endif


   ! Define the path to the skf file:
   select case (trim(adjustl(path_to_skf)))  ! how to set it
   case ('PATH', 'Path', 'path') ! then specify exactly the path to the skf file
      read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%param_name  ! read the full path
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3426
      endif
      Folder_name = trim(adjustl(TB_Hamil(i,j)%param_name))    ! folder with chosen parameters sets
   case default   ! it is a parameterization within the predefined directory 'DFTB'
      TB_Hamil(i,j)%param_name = trim(adjustl(path_to_skf)) ! name of the directory with skf files
      ! folder with all DFTB data:
      Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_DFTB_norep_directory))//numpar%path_sep
      Folder_name = trim(adjustl(Folder_name))//trim(adjustl(TB_Hamil(i,j)%param_name))
   endselect

   ! Make sure the slash is correct:
   call ensure_correct_path_separator(Folder_name, numpar%path_sep)  ! module "Dealing_with_files"

   ! folder with chosen parameters sets:
   select case (trim(adjustl(TB_Hamil(i,j)%param_name)))
   case default ! e.g. '1element'
      ! no inner folders, just skf-files
   case ('2elements')
      read(FN,*,IOSTAT=Reason) Inner_folder_name   ! read the inner folder name
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
         INFO = 3
         goto 3426
      endif
      Folder_name = trim(adjustl(Folder_name))//numpar%path_sep//trim(adjustl(Inner_folder_name))
      TB_Hamil(i,j)%param_name = trim(adjustl(TB_Hamil(i,j)%param_name))//numpar%path_sep//trim(adjustl(Inner_folder_name))
   endselect

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rcut, TB_Hamil(i,j)%d  ! [A] cut off, and width of cut-off region [A]
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3426
   endif

   ! Assume it is a file name:
   inquire(file=trim(adjustl(Folder_name)),exist=file_exist)
   if (file_exist) then ! such a file exists, use it
      File_name = trim(adjustl(Folder_name))
   else  ! no file => assume it is a directory name, and use default file name:
      ! Construct name of the skf file:
      call construct_skf_filename( trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                                 File_name, '_no_repulsion')   ! module "Dealing_with_DFTB"
      File_name = trim(adjustl(Folder_name))//numpar%path_sep//trim(adjustl(File_name))
   endif

   ! Check if such DFTB parameterization exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exist)
   if (.not.file_exist) then
      Error_descript = 'File '//trim(adjustl(File_name))//' not found, the program terminates'
      INFO = 1
      goto 3426
   endif
   FN_skf=111
   open(UNIT=FN_skf, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      INFO = 2
      goto 3426
   endif

   ToA = same_or_different_atom_types(trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name))) ! module "Dealing_with_DFTB"
   call read_skf_file_no_rep(FN_skf, TB_Hamil(i,j), TB_Repuls(i,j), ToA, Error_descript)    ! module "Dealing_with_DFTB"
   if (LEN(trim(adjustl(Error_descript))) > 0) then
      INFO = 5
      goto 3426
   endif

   ! Check which basis set is used: 0=s, 1=sp3, 2=sp3d5:
   if ((i == matter%N_KAO) .and. (j == matter%N_KAO)) then  ! only when all parameters for all elements are read from files:
      call idnetify_basis_size(TB_Hamil, N_basis_siz)  ! module "Dealing_with_DFTB"'
      numpar%N_basis_size = max(numpar%N_basis_size,N_basis_siz)
   endif

3426 continue
   ! Close files that have been read through:
   call close_file('close', FN=FN_skf) ! module "Dealing_with_files"
   call close_file('close', FN=FN) ! module "Dealing_with_files"
end subroutine read_DFTB_TB_Params_no_rep



subroutine read_3TB_TB_Params(FN, i, j, TB_Hamil, numpar, matter, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_3TB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter	! all material parameters
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------------------------------------
   character(100) :: Folder_name, File_name
   character(200) :: Filename_onsite, Filename_2body, Filename_3body
   integer count_lines, Reason, FN_onsite, FN_2bdy, FN_3bdy, N_basis_siz
   logical file_exist, file_opened, read_well, file_exists

   ! To start with:
   INFO = 0
   count_lines = 1
   FN_onsite=110
   FN_2bdy=111
   FN_3bdy=112

   ! Read the second line of the input file:
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rcut, TB_Hamil(i,j)%d  ! [A] cut off, and width of cut-off region [A]
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3500
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rc  ! distance rescaling coefficient
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3500
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%include_3body  ! flag, include 3-body terms or not
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      TB_Hamil(i,j)%include_3body = .false.  ! by default, exclude 3-body terms
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%nullify_diag_cf  ! flag, exclude diagonal part of crystal field
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      TB_Hamil(i,j)%nullify_diag_cf = .true.  ! by default, exclude diagonal part of crystal field
   endif


   ! Folder with all 3TB data:
   Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_3TB_directory))//numpar%path_sep

   ! Read onsite parameters:
   Filename_onsite = trim(adjustl(Folder_name))//trim(adjustl(m_3TB_onsite_data))
   inquire(file=trim(adjustl(Filename_onsite)),exist=file_exists)
   if (.not.file_exists) then
      INFO = 1 ! file with onsite parameters not found
      Error_descript = 'File '//trim(adjustl(Filename_onsite))//' not found, the program terminates'
      goto 3500
   endif
   open(UNIT=FN_onsite, FILE = trim(adjustl(Filename_onsite)), status = 'old', action='read')
   inquire(file=trim(adjustl(Filename_onsite)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(Filename_onsite))//' could not be opened, the program terminates'
      INFO = 2
      goto 3500
   endif
   ! Having constructed the name of the file with onsite parameters, read it:
   call read_3TB_onsite_file(FN_onsite, trim(adjustl(matter%Atoms(i)%Name)), TB_Hamil(i,j), &
                              N_basis_siz, Error_descript)    ! module "Dealing_with_3TB"
   if (LEN(trim(adjustl(Error_descript))) > 0) then
      INFO = 5
      goto 3500
   endif

   ! Save the basis size as the maximal among all elements:
   numpar%N_basis_size = max(numpar%N_basis_size,N_basis_siz)



   ! Construct name of the 3TB files with 2-body and 3-body parameters:
   call construct_3TB_filenames(trim(adjustl(Folder_name)), trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                     numpar%path_sep, Filename_2body, Filename_3body, INFO)  ! module "Dealing_with_3TB"
   ! Check if such 3TB parameterization exists:
   if (INFO /= 0) then
      select case(INFO)
      case (1)
         Error_descript = 'File '//trim(adjustl(Filename_2body))//' not found, the program terminates'
         INFO = 1
      case (2)
         Error_descript = 'File '//trim(adjustl(Filename_3body))//' not found, the program terminates'
         INFO = 1
      case (3)
         Error_descript = 'File '//trim(adjustl(Filename_2body))//' not found, the program terminates'
         INFO = 1
      case (4)
         Error_descript = 'File '//trim(adjustl(Filename_3body))//' not found, the program terminates'
         INFO = 1
      end select
      goto 3500
   endif

   ! Open and read from the parameters files:
   ! 2-body parameters:
   open(UNIT=FN_2bdy, FILE = trim(adjustl(Filename_2body)), status = 'old', action='read')
   inquire(file=trim(adjustl(Filename_2body)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(Filename_2body))//' could not be opened, the program terminates'
      INFO = 2
      goto 3500
   endif
   ! File exists and opened, read from it:
   call read_3TB_2bdy_file(FN_2bdy, trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                           TB_Hamil(i,j), Error_descript)    ! module "Dealing_with_3TB"
   if (LEN(trim(adjustl(Error_descript))) > 0) then
      INFO = 5
      goto 3500
   endif


   ! 3-body parameters:
   if (TB_Hamil(i,j)%include_3body) then  ! only if user defined it to include

      open(UNIT=FN_3bdy, FILE = trim(adjustl(Filename_3body)), status = 'old', action='read')
      inquire(file=trim(adjustl(Filename_3body)),opened=file_opened)
      if (.not.file_opened) then
         Error_descript = 'File '//trim(adjustl(Filename_3body))//' could not be opened, the program terminates'
         INFO = 2
         goto 3500
      endif
      ! File exists and opened, read from it:
      call read_3TB_3bdy_file(FN_3bdy, trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                           TB_Hamil(i,j), Error_descript)    ! module "Dealing_with_3TB"
      if (LEN(trim(adjustl(Error_descript))) > 0) then
         INFO = 5
         goto 3500
      endif
   endif ! (TB_Hamil(i,j)%include_3body)


3500 continue
   ! Close files that have been read through:
   call close_file('close', FN=FN_onsite) ! module "Dealing_with_files"
   call close_file('close', FN=FN_2bdy) ! module "Dealing_with_files"
   call close_file('close', FN=FN_3bdy) ! module "Dealing_with_files"
   call close_file('close', FN=FN) ! module "Dealing_with_files"


   ! While testing :
   if (INFO .NE. 0) then
      print*, trim(adjustl(Error_descript))
   endif
   !PAUSE 'pause read_3TB_TB_Params'

end subroutine read_3TB_TB_Params



subroutine read_BOP_TB_Params(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_BOP), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_BOP), dimension(:,:), intent(inout) ::  TB_Repuls    ! parameters of the repulsive potential
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter	! all material parameters
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------------------------------------
   character(200) :: Folder_name, File_name
   real(8) :: bond_length   ! to construct repulsive part of BOP potential, we need to know dimer bond length [A]
   integer count_lines, Reason, FN_BOP, N_basis_siz
   logical :: file_exist, file_opened, read_well, file_exists, data_exists
   INFO = 0
   count_lines = 2
   
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rcut, TB_Hamil(i,j)%dcut  ! [A] cut off, and width of cut-off region [A]
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3427
   endif

   Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_BOP_directory))//numpar%path_sep  ! folder with all BOP data
   File_name = trim(adjustl(Folder_name))//numpar%path_sep//trim(adjustl(m_BOP_file))   ! file with BOP parameters, module "Dealing_with_BOP"

   ! Check if such BOP parameterization exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exist)

   if (.not.file_exist) then
      Error_descript = 'File '//trim(adjustl(File_name))//' not found, the program terminates'
      INFO = 1
      goto 3427
   endif
   FN_BOP=1111
   open(UNIT=FN_BOP, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      INFO = 2
      goto 3427
   endif

   if (j >= i) then ! reading lower triangle together wit the upper one:
      call read_BOP_parameters(FN_BOP, trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                                            TB_Hamil, i, j, Error_descript)    ! module "Dealing_with_BOP"
   endif

   if (LEN(trim(adjustl(Error_descript))) > 0) then
      INFO = 5
      goto 3427
   endif
   
   ! Check which basis set is used: 0=s, 1=sp3, 2=sp3d5:
   if (i == j) then
      call idnetify_basis_size_BOP(TB_Hamil(i,j), N_basis_siz)  ! module "Dealing_with_BOP"'
      numpar%N_basis_size = max(numpar%N_basis_size, N_basis_siz)
   endif

   ! Now, deal with the repulsive part of BOP:
   if (j >= i) then ! reading lower triangle together wit the upper one:
      ! Make sure either file with repulsive potential exists, or we can create it:
      call check_if_repulsion_exists(INT(matter%Atoms(i)%Z), trim(adjustl(matter%Atoms(i)%Name)), &
                INT(matter%Atoms(j)%Z), trim(adjustl(matter%Atoms(j)%Name)), &
                Folder_name, numpar%path_sep, file_exists, data_exists, bond_length, Error_descript) ! module "Dealing_with_BOP"
      if (LEN(trim(adjustl(Error_descript))) > 0) then
         INFO = 4
         goto 3427
      endif

      ! If we have file with repulsive potential, read it from it:
      if (file_exists) then ! read from it:
         call read_BOP_repulsive(TB_Repuls, i, j, Folder_name, numpar%path_sep, &
            trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), Error_descript) ! module "Dealing_with_BOP"

      elseif (data_exists) then ! construct new repulsive potential:
         numpar%BOP_bond_length = bond_length   ! [A] save to reuse later
         numpar%create_BOP_repulse = .true. ! marker to construct BOP repulsive potential
         numpar%BOP_Folder_name = Folder_name   ! directory where to find it

      else ! no way to access repulsive part
         INFO = 5
         goto 3427
      endif
   endif
   
3427 continue 
   ! Close files that have been read through:
   call close_file('close', FN=FN_BOP) ! module "Dealing_with_files"
   call close_file('close', FN=FN) ! module "Dealing_with_files"
end subroutine read_BOP_TB_Params




subroutine read_xTB_Params(FN, i,j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO)
   integer, intent(in) :: FN ! file number where to read from
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_xTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_xTB), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter	! all material parameters
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------------------------------------
   character(100) :: Folder_name, File_name
   integer count_lines, Reason, i_cur, ind, FN_skf, ToA, N_basis_siz
   logical file_exist, file_opened, read_well
   INFO = 0 ! to start with no error
   count_lines = 2

   ! name of the xTB parameterization (currently, only GFN0 is supported!); number of GTO primitives
   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%param_name, TB_Hamil(i,j)%Nprim
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3428
   endif

   read(FN,*,IOSTAT=Reason) TB_Hamil(i,j)%rcut, TB_Hamil(i,j)%d  ! [A] cut off, and width of cut-off region [A]
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
       write(Error_descript,'(a,i3)') 'Could not read line ', count_lines
       INFO = 3
       goto 3428
   endif

   Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_xTB_directory))//numpar%path_sep ! folder with xTB data

   ! Read xTB parameters:
   call read_xTB_parameters(Folder_name, i, j, TB_Hamil, TB_Repuls, numpar, matter, Error_descript, INFO)    ! module "Dealing_with_xTB"

   ! Check which basis set is used: 0=s; 1=ss*; 2=sp3; 3=sp3s*; 4=sp3d5; 5=sp3d5s*
   if ((i == matter%N_KAO) .and. (j == matter%N_KAO)) then  ! only when all parameters for all elements are read from files:
      call identify_basis_size_xTB(TB_Hamil, N_basis_siz)  ! module "Dealing_with_DFTB"'
      ! Save the index of the basis set:
      numpar%N_basis_size = max(numpar%N_basis_size,N_basis_siz)

      ! Identify the parameters of the AO:
      call identify_AOs_xTB(TB_Hamil)   ! module "Dealing_with_xTB"


   endif

3428 continue
   ! Close files that have been read through:
   call close_file('close', FN=FN) ! module "Dealing_with_files"
end subroutine read_xTB_Params



subroutine read_numerical_parameters(File_name, matter, numpar, laser, Scell, user_data, Err, add_data, count_lines_in)
   character(*), intent(in) :: File_name
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(User_overwrite_data), intent(inout) :: user_data   ! atomic data provided by the user
   type(Error_handling), intent(inout) :: Err	! error save
   logical, intent(in), optional :: add_data   ! read stricktly numpar data, no additional data
   integer, intent(inout), optional :: count_lines_in ! what line we already arrived to in the file (if not starting)
   !---------------------------------
   integer FN, N, Reason, count_lines, i, NSC, temp1, temp2, temp3
   logical file_opened, read_well, old_file, add_data_present
   character(200) Error_descript, temp_ch, temp_ch2, read_line

   NSC = 1 ! for now, we only have 1 supercell...
   if (present(add_data)) then
      add_data_present = add_data ! user defines whether to read additional data
   else
      add_data_present = .true. ! by default, read possible additional data
   endif

   if (present(count_lines_in)) then
      count_lines = count_lines_in   ! mark that we start from this line in the file
   else
      count_lines = 0   ! mark that we start from top of the file
   endif

   !inquire(file=trim(adjustl(File_name)),opened=file_opened)
   inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN) ! if file is opened, use its number
   if (file_opened) then ! file is already opened, continue reading from it
      old_file = .true. ! mark that this file was laready opened
   else ! file is not open, open and read from it:
      old_file = .false. ! mark that this file was not opened yet
      FN=108
      open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      inquire(file=trim(adjustl(File_name)),opened=file_opened)
      if (.not.file_opened) then
         Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
         call Save_error_details(Err, 2, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3418
      endif
   endif

   ! number of unit-cells in X,Y,Z:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) matter%cell_x, matter%cell_y, matter%cell_z
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! periodicity along X,Y,Z directions:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) temp1, temp2, temp3
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   numpar%r_periodic(:) = .true.	! periodic by default
   if (temp1 == 0) numpar%r_periodic(1) = .false.	! along X
   if (temp2 == 0) numpar%r_periodic(2) = .false.	! along Y   
   if (temp3 == 0) numpar%r_periodic(3) = .false.	! along Z

   ! where to take atomic data from (EADL, CDF):
   numpar%user_defined_E_gap = -1.0d0 ! default
   !read(FN,*,IOSTAT=Reason) numpar%At_base
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%At_base, numpar%input_CDF_file, numpar%user_defined_E_gap
   if (Reason /= 0) then ! try different order of variables
      read(read_line,*,IOSTAT=Reason) numpar%At_base, numpar%user_defined_E_gap, numpar%input_CDF_file
   endif
   if (Reason /= 0) then ! try 2 variables: including E_gap:
      numpar%user_defined_E_gap = -1.0d0  ! default
      numpar%input_CDF_file = ''  ! nullify it
      read(read_line,*,IOSTAT=Reason) numpar%At_base, numpar%user_defined_E_gap
   endif
   if (Reason /= 0) then ! try to read 2 variables: including path to CDF-file:
      numpar%user_defined_E_gap = -1.0d0 ! default
      numpar%input_CDF_file = ''  ! nullify it
      read(read_line,*,IOSTAT=Reason) numpar%At_base, numpar%input_CDF_file
   endif
   if (Reason /= 0) then ! try to read just single variable:
      numpar%input_CDF_file = ''  ! nullify it
      read(read_line,*,IOSTAT=Reason) numpar%At_base
      if (numpar%verbose) write(*,'(a)') 'No valid filename with CDF oscillators provided, assuming default'
   else
      ! Make sure that, if there is a path, the separator is correct:
      call ensure_correct_path_separator(numpar%input_CDF_file, numpar%path_sep) ! module "Dealing_with_files"
   endif
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! [g/cm^3] density of the material (used in MC in case of EADL parameters):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) matter%dens
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! number of iterations in the MC module:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%NMC
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (numpar%NMC < 0) numpar%NMC = 0  ! by default, no MC

   ! number of threads for OPENMP:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%NOMP
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (numpar%NOMP < 1) then ! use default: maximum number of available threads
#ifdef OMP_inside
      numpar%NOMP = omp_get_max_threads() ! number of processors available by default
#else ! if you set to use OpenMP in compiling: 'make OMP=no'
      numpar%NOMP = 1
#endif
   endif
   
   ! MD algorithm (0=Verlet, 2d order; 1=Yoshida, 4th order)
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%MD_algo
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! Include (1) or exclude (0) atopmic motion:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%do_atoms = .true.	! Atoms move
   else
      numpar%do_atoms = .false.	! Frozen atoms
   endif

   ! Parinello-Rahman super-vell mass coefficient
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) matter%W_PR
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (matter%W_PR < 0.0d0) matter%W_PR = 1.0d0 ! use default value

   ! Time step for MD [fs]:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%MD_step_grid_file ! file with time grid, or timestep for md [fs]
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   ! If read well, interpret it and set timestep or time-grid:
   call set_MD_step_grid(numpar%MD_step_grid_file, numpar, read_well, Error_descript)    ! below

   ! save data into files every 'dt_save_time' [fs]
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%dt_save  ! save data into files every 'dt_save_time' [fs]
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! it's = 1 if P=const, or = 0 if V=const:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N	! It's = 0 if P=const, or = 1 if V=const
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%p_const = .true.	! P=const
   else
      numpar%p_const = .false.	! V=const
   endif

   ! external pressure [Pa] (0 = normal atmospheric pressure):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) matter%p_ext  ! External pressure [Pa]
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
!    if (matter%p_ext < 0.0d0) matter%p_ext = g_P_atm	! atmospheric pressure

   ! SCC:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%scc, numpar%scc_gam_ind, numpar%scc_mix
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! 4 numbers:
   ! 1) Scheme of propagation of electronic ensemble:
   ! (0=decoupled electrons; 1=enforced energy conservation; 2=T=const; 3=BO; 4=relaxation time appeoximation);
   ! 2) Global characteristic relaxation time of ALL electrons [fs];
   ! 3) Characteristic relaxation time of electrons in CB (above Fermi-level) [fs];
   ! 3) Characteristic relaxation time of electrons in VB (below Fermi-level) [fs];
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%el_ion_scheme, numpar%tau_fe, numpar%tau_fe_CB, numpar%tau_fe_VB
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   numpar%do_partial_thermal = ( (numpar%tau_fe_CB > -1.0d-8) .and. (numpar%tau_fe_VB > -1.0d-8) )
   if (numpar%tau_fe_CB < 0.0d0) numpar%tau_fe_CB = 0.0d0   ! eliminate nigative values (even within precision)
   if (numpar%tau_fe_VB < 0.0d0) numpar%tau_fe_VB = 0.0d0   ! eliminate nigative values (even within precision)

   ! -1=nonperturbative (default), 0=no coupling, 1=dynamical coupling, 2=Fermi golden rule (DO NOT USE!):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%NA_kind
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (numpar%NA_kind .EQ. 0) then
      numpar%Nonadiabat = .false. ! excluded
   else
      numpar%Nonadiabat = .true.  ! included
   endif

   ! [fs] when to switch on the nonadiabatic coupling:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%t_NA, numpar%M2_scaling ! [fs] start of the nonadiabatic coupling; scaling factor
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! [eV] acceptance window and quasidegeneracy window for nonadiabatic coupling:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%acc_window, numpar%degeneracy_eV
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! atoms quenching (0=no, 1=yes); starting from when [fs]; how often [fs]:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N, numpar%at_cool_start, numpar%at_cool_dt ! include atomic cooling? When to start? How often?
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%do_cool = .true.	! included
   else
      numpar%do_cool = .false.	! excluded
   endif

   ! 0=no heat transport, 1=include heat transport; thermostat temperature for ATOMS [K]:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N, matter%T_bath, matter%tau_bath
   if (Reason == 0) then
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
      if (Err%Err) goto 3418
      if (N .EQ. 1) then
         numpar%Transport = .true.	 ! included
      else
         numpar%Transport = .false. ! excluded
      endif
      matter%T_bath = matter%T_bath/g_kb	! [eV] thermostat temperature for atoms
   else ! maybe there is a filename given to read from instead of numbers:
      read(read_line,*,IOSTAT=Reason) numpar%At_bath_step_grid_file  ! name of file with parameters
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
      if (Err%Err) goto 3418
      ! If read well, try to interpret it and set parameters on time-grid:
      call set_Bath_grid_atoms(numpar, read_well, Error_descript)    ! below
      if (.not. read_well) then
         goto 3418
      endif
   endif

   ! 0=no heat transport, 1=include heat transport; thermostat temperature for ELECTRONS [K]:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N, matter%T_bath_e, matter%tau_bath_e
   if (Reason == 0) then
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
      if (Err%Err) goto 3418
      if (N .EQ. 1) then
         numpar%Transport_e = .true.	 ! included
      else
         numpar%Transport_e = .false. ! excluded
      endif
      matter%T_bath_e = matter%T_bath_e/g_kb	! [eV] thermostat temperature for electrons
   else
      ! maybe there is a filename given to read from instead of numbers:
      read(read_line,*,IOSTAT=Reason) numpar%El_bath_step_grid_file  ! name of file with parameters
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
      if (Err%Err) goto 3418
      ! If read well, try to interpret it and set parameters on time-grid:
      call set_Bath_grid_electrons(numpar, read_well, Error_descript)    ! below
      if (.not. read_well) then
         goto 3418
      endif
   endif

   ! [eV] cut-off energy, separating low-energy-electrons from high-energy-electrons:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%E_cut  ! [eV] cut-off energy for high-energy-electrons
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (numpar%E_cut <= 0.0d0) then ! use dynamical evolution of E_cut adjusting it to top-most CB level
      numpar%E_cut_dynamic = .true.  ! change E_cut
   else
      numpar%E_cut_dynamic = .false. ! do not change E_cut
   endif

   ! [eV] work function, for electron emission:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%E_work  ! [eV] work function for electron emission
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (numpar%E_work <= 0.0d0) then ! it is a condition on the counter of collisions:
      ! don't forget to exclude electrons that made more collisions than allowed
   else if (numpar%E_work <= numpar%E_cut) then ! exclude it from the calculations:
      numpar%E_work = 1.0d30
   endif
   
   ! save electron energy levels (1) or not (0):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N  ! save electron energy levels (1) or not (0)
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%save_Ei = .true.	! included
   else
      numpar%save_Ei = .false.	! excluded
   endif
   
   ! save DOS (1=Gamma; 2=k-points) or not (0):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N, numpar%Smear_DOS, numpar%DOS_splitting ! save DOS (1) or not (0), smearing width, do partial DOS or no
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (abs(N) == 2) then
      ! if the user wanted complex k-points, leave it be, otherwise, overwrite it with the option of complex k for DOS
      ! * Note that it SWITCHES OFF the probe-pulse calculations in this case
      ! which means it is currently not possible to get gamma-point probe pulse, and multiple-k-points for DOS.
      ! They are connected: either both are calculated for multiple (and the same) k-points,
      ! or DOS calculation takes precedence, and probe is switched off:
      if ((numpar%optic_model /= 2) .and. (abs(numpar%optic_model) /= 4) .and. (numpar%optic_model /= 5)) then
         numpar%optic_model = -5    ! complex, for given number of k-points for DOS calculations
      endif
      numpar%save_DOS = .true.   ! included
   elseif (N == 1) then ! calculate with the same model used for the probe pulse (gamma-point or multiple k-points)
      numpar%save_DOS = .true.   ! included
   else ! no DOS output requested
      numpar%save_DOS = .false.  ! excluded
   endif
!    print*, 'DOS_in:', numpar%save_DOS, numpar%optic_model

   ! save Mulliken or not, and within which model: (0) no; (1) for atom types; 
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%Mulliken_model
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! save electron electron distribution (1) or not (0):
   !read(FN,*,IOSTAT=Reason) N  ! save electron distribution function (1) or not (0)
   read(FN, '(a)', IOSTAT=Reason) read_line  ! read parameters to interpret them in a subroutine
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   call interprete_distribution_input(read_line, numpar, Scell(1), read_well) ! below
   if (.not. read_well) then
      write(Error_descript,'(a,i5,a,$)') 'Could not interprete line ', count_lines, ' in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3418
   endif

   ! save atomic pair correlation function (1) or not (0):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N  ! save atomic pair correlation function (1) or not (0)
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%save_PCF = .true.	! included
   else
      numpar%save_PCF = .false.	! excluded
   endif

   ! save atomic positions in XYZ (1) or not (0):
   !read(FN,*,IOSTAT=Reason) N  ! save atomic positions in XYZ (1) or not (0)
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,'(a)',IOSTAT=Reason) temp_ch   ! read parameters to interpret them below
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   read(temp_ch,*,IOSTAT=Reason) N, temp_ch2
   if (Reason /= 0) then ! something wrong with input
      ! try reading just the first flag, no grid parameters
      read(temp_ch,*,IOSTAT=Reason) N
      count_lines = count_lines - 1 ! still reading the same line
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
      if (Err%Err) goto 3418
      ! And no additional input:
      temp_ch2 = ''
   endif
   if (N .EQ. 1) then
      numpar%save_XYZ = .true.	! included
   else
      numpar%save_XYZ = .false.	! excluded
   endif
   !print*, 'temp_ch2', temp_ch2

   call interpret_additional_XYZ_input(temp_ch2, numpar%save_XYZ_extra) ! below
   !pause 'interpret_additional_XYZ_input'

   ! save atomic positions in CIF (1) or not (0):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N  ! save atomic positions in CIF (1) or not (0)
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%save_CIF = .true.	! included printout atomic coordinates in CIF format
   else
      numpar%save_CIF = .false.	! excluded printout atomic coordinates in CIF format
   endif

   ! save raw data for atomic positions and velocities (1) or not (0):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (N .EQ. 1) then
      numpar%save_raw = .true.	! included printout raw data
   else
      numpar%save_raw = .false.	! excluded printout raw data
   endif

   ! read information about displacement:
   ! command, filename, or power of mean displacement to print out (set integer N: <u^N>-<u0^N>):
   read(FN, '(a)', IOSTAT=Reason) read_line
   !read(read_line,*,IOSTAT=Reason) numpar%MSD_power
   call interprete_displacement_command(read_line, Scell(1), numpar, Reason) ! below

   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   ! save number of nearest neighbors within the digen radius (>0) or not (<=0):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%NN_radius
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   if (numpar%NN_radius > 1.0e-2) then
      numpar%save_NN = .true.	! included printout nearest neighbors
   else
      numpar%save_NN = .false.	! excluded printout nearest neighbors
   endif

   !  which format to use to plot figures: eps, jpeg, gif, png, pdf
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%fig_extention
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418
   select case ( trim(adjustl(numpar%fig_extention)) )
   case ('JPEG', 'JPEg', 'JPeg', 'Jpeg', 'jpeg', 'JPG', 'jpg')
      numpar%fig_extention = 'jpeg'
      numpar%ind_fig_extention = 2
   case ('GIF', 'GIf', 'Gif', 'gif')
      numpar%fig_extention = 'gif'
      numpar%ind_fig_extention = 3
   case ('PNG', 'PNg', 'Png', 'png')
      numpar%fig_extention = 'png'
      numpar%ind_fig_extention = 4
   case ('PDF', 'PDf', 'Pdf', 'pdf')
      numpar%fig_extention = 'pdf'
      numpar%ind_fig_extention = 5
   case default ! eps
      numpar%fig_extention = 'eps'
      numpar%ind_fig_extention = 1
   end select

   ! number of k-points in each direction (used only for Trani-k!):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%ixm, numpar%iym, numpar%izm
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
   if (Err%Err) goto 3418

   !----------------------------------------
   ! Read optional data provided by the user (e.g., to overwrite default atomic data):
   if (add_data_present) then
      call read_user_additional_data(FN, count_lines, matter, Scell, NSC, user_data)   ! below
   endif

   ! Close this file only if it was opened within this subroutine:
3418 continue
   if (.not.old_file .and. file_opened) close(FN)
end subroutine read_numerical_parameters



subroutine interprete_displacement_command(read_line, Scell, numpar, Reason)
   character(*), intent(in) :: read_line
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   integer, intent(inout) :: Reason
   !-------------------------------
   character(200) :: Folder_name, File_name, sting_temp
   character(10) :: ch_temp
   logical :: file_exists
   integer :: FN1, ch_int, i

   ! Try to interprete it as a number (power of mean displacement)
   ! to make back-compatible with the legacy format:
   read(read_line,*,IOSTAT=Reason) numpar%MSD_power
   if (Reason == 0) then
      if (numpar%verbose) print*, 'Atomic displacement analysis is set in legacy format'
      return ! it was a number, nothing more to do
   endif

   !---------------
   ! If it was not a command, check if it was a file
   Folder_name = trim(adjustl(m_INPUT_directory))//numpar%path_sep   ! directory with input files
   read(read_line,*,IOSTAT=Reason) sting_temp   ! read first block
   File_name = trim(adjustl(Folder_name))//trim(adjustl(sting_temp))  ! file name with masks
   inquire(file=trim(adjustl(File_name)),exist=file_exists)
   if (file_exists) then
      FN1 = 4001
      open(UNIT=FN1, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      read(FN1,*,IOSTAT=Reason) ch_temp, ch_int
      if (Reason == 0) then   ! if read number of masks well
         ! Allocate displacements data objects:
         allocate(Scell%Displ(ch_int))
         ! Read parameters of the masks:
         do i = 1, ch_int  ! for all masks:
            read(FN1,'(a)',IOSTAT=Reason) sting_temp
            if (Reason == 0) then
               call read_displacement_command(trim(adjustl(sting_temp)), Scell, numpar, Reason, i, FN1) ! below

            else
               print*, 'Could not read atomic masks from file: ', trim(adjustl(File_name))
               print*, 'Using default mean displacement only'
               Reason = -1
               exit
            endif
         enddo ! i
      endif ! (Reason == 0)
      close (FN1)
      if (Reason == 0) then
         if (numpar%verbose) print*, 'Atomic masks for displacements are read from file: ', trim(adjustl(File_name))
         return ! read all, nothing more to do
      endif
   else ! check if it is a command:
      ! if it is not a number, check if it is a command:
      call read_displacement_command(trim(adjustl(read_line)), Scell, numpar, Reason, 1) ! below
      if (Reason == 0) then
         if (numpar%verbose) print*, 'Atomic masks for displacements are set by command'
         return ! it was a number, nothing more to do
      else
         print*, 'Could not interprete command for atomic masks'
         print*, 'Using default mean displacement only'
      endif
   endif

   !---------------
   ! 4) if it was none of the above, just use default: mean displacement N=1:
   numpar%MSD_power = 1
   Reason = 0
   if (numpar%verbose) print*, 'Default atomic displacement analysis is used'
end subroutine interprete_displacement_command


subroutine read_displacement_command(read_line, Scell, numpar, Reason, mask_num, FN1)
   character(*), intent(in) :: read_line
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   integer, intent(inout) :: Reason
   integer, intent(in) :: mask_num  ! number of mask
   integer, intent(in), optional :: FN1 ! file number with masks definition
   !-------------------------------
   integer :: Nsiz, N_char, N_len, i, j
   character(100) :: ch_comm, ch_val
   character(200) :: ch_temp

   Reason = -1 ! to start with

   !---------------
   ! Make sure the array is allocated
   if (.not.allocated(Scell%Displ)) then
      allocate(Scell%Displ(1))
   endif

   !---------------
   ! Set defaults:
   Scell%Displ(mask_num)%MSD_power = 1.0d0   ! default: linear mean displacement
   write(ch_temp, '(i6)') mask_num
   Scell%Displ(mask_num)%mask_name = 'mask_'//trim(adjustl(ch_temp)) ! default namae
   Scell%Displ(mask_num)%print_r = .true.   ! all axis-resolved displacement

   !---------------
   ! Interprete the input line:
   N_len = LEN(read_line)  ! length of the string
   N_char = 1  ! to start with
   do while (N_char < N_len)  ! read the entire line
      read(read_line(N_char:N_len),*,IOSTAT=Reason) ch_temp ! read this block
      N_char = N_char + LEN(trim(ch_temp)) + 1  ! mark begining of the next block to read

      ! Split the string into command and its value:
      call split_command_separator(trim(adjustl(ch_temp)), ':', ch_comm, ch_val)  ! below
      ! If the command is incorrect, nothing more to do
      if (LEN(trim(adjustl(ch_comm))) == 0) cycle
      if (LEN(trim(adjustl(ch_val))) == 0) cycle

      select case(trim(adjustl(ch_comm)))
      !==========
      case ('name', 'Name', 'NAME')   ! mask name
         ! Make sure the name does not repeat:
         call number_atomic_mask(ch_val, Scell, mask_num) ! below
         Scell%Displ(mask_num)%mask_name = trim(adjustl(ch_val))
         Reason = 0

         ! Identify the mask format, and read extra parameters if any:
         select case (trim(adjustl(Scell%Displ(mask_num)%mask_name(1:3)))) ! define section
         case ('sec', 'Sec', 'SEC') ! Spatial section
            ! Read one more line with the definition of the section:
            if (present(FN1)) then ! read next line from this file
               call read_and_define_section(FN1, Scell, mask_num) ! below
            endif

         case ('all', 'All', 'ALL') ! all atoms
            ! Nothing to do, all atoms are included
         end select

      !==========
      case ('power', 'Power', 'POWER')   ! power of mean displacement
         read(ch_val,*,IOSTAT=Reason) Scell%Displ(mask_num)%MSD_power
         Reason = 0

      !==========
      case ('axis', 'Axis', 'AXIS', 'axes', 'Axes', 'AXES')   ! axis-resolved data
      ! (specification unused, all axis are printed out)
         do i = 1, LEN(trim(adjustl(ch_val)))   ! interprete all letters
            select case(trim(adjustl(ch_val(i:i))))
            case ('X', 'x')
               Scell%Displ(mask_num)%print_r(1) = .true.
            case ('Y', 'y')
               Scell%Displ(mask_num)%print_r(2) = .true.
            case ('Z', 'z')
               Scell%Displ(mask_num)%print_r(3) = .true.
            endselect
         enddo
         Reason = 0
      end select

      if (LEN(trim(adjustl(ch_temp))) < 1) exit ! nothing more to read
   enddo

   !pause 'read_displacement_command'
end subroutine read_displacement_command


subroutine number_atomic_mask(mask_name, Scell, mask_num)
   character(*), intent(inout) :: mask_name
   type(Super_cell), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: mask_num
   !---------------
   integer :: j, rep_num, Reason
   character(100) :: string_part1, string_part2, string_temp, ch_mask_num

   do j = 1, mask_num
      string_temp = trim(adjustl(Scell%Displ(j)%mask_name))  ! to work with
      ! If the name repeats, add a number to it at the end:
      if ( trim(adjustl(mask_name)) == trim(adjustl(string_temp)) ) then
         ! Find if there is already a number assigned:
         call split_command_separator(trim(adjustl(string_temp)), '_', string_part1, string_part2, back=.true.)  ! below

         if (LEN(trim(adjustl(string_part2))) == 0) then ! no numnber in the name
            mask_name = trim(adjustl(mask_name))//'_1'   ! make it the first
         else  ! there is a number, add to it
            read(string_part2, *, iostat=Reason) rep_num
            if (Reason == 0) then ! read well
               write(ch_mask_num, '(i5)') rep_num+1   ! next number
               mask_name = trim(adjustl(string_part1))//'_'//trim(adjustl(ch_mask_num))   ! add this number to the name
            else  ! not a number, just underscore in the name
               mask_name = trim(adjustl(string_temp))//'_1' ! add the first number
            endif
         endif
      endif
   enddo
end subroutine number_atomic_mask



subroutine read_and_define_section(FN1, Scell, mask_num)
   integer, intent(in) :: FN1 ! file number
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: mask_num  ! mask number
   !---------------------------
   character(200) :: sting_temp, sting_temp2
   character(200), dimension(3) :: string_part
   character(100) :: ch_comm, ch_val, block1, block2
   integer :: Reason, i

   read(FN1,'(a)',IOSTAT=Reason) sting_temp
   if (Reason /= 0) then   ! couldn't read, use default
      ! No mask parameters, cannot do section, use all instead
      Scell%Displ(mask_num)%mask_name = 'all'
      return
   else  ! read something, attempt to interprete it:
      ! Define default parameters of the section:
      Scell%Displ(mask_num)%logical_and = .false.  ! to start with
      Scell%Displ(mask_num)%logical_or = .false.   ! to start with
      Scell%Displ(mask_num)%r_start = -1.0d10   ! start at -infinity
      Scell%Displ(mask_num)%r_end = 1.0d10      ! end at +infinity

      !---------------
      ! Check if there is a separator:
      string_part = ''  ! to start with

      call split_command_separator(trim(adjustl(sting_temp)), ';', string_part(1), string_part(2))  ! below
      if (LEN(trim(adjustl(string_part(2)))) /= 0) then
         sting_temp2 = string_part(2)
         call split_command_separator(trim(adjustl(sting_temp2)), ';', string_part(2), string_part(3))  ! below
      endif

      ! Read and interprete all 3 parts:
      do i = 1, 3
         block1 = '' ! to start with
         block2 = '' ! to start with
         ch_comm = '' ! to start with
         ch_val = '' ! to start with

         if (LEN(trim(adjustl(string_part(i)))) > 0) then
            sting_temp = trim(adjustl(string_part(i)))
         else  ! no text to interprete here
            if (i > 1) cycle  ! nothing else to do, skip it
         endif

         !---------------
         ! Check if there is an 'and':
         call split_command_separator(trim(adjustl(sting_temp)), 'and', ch_comm, ch_val)  ! below
         ! Check other possible ways:
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'AND', ch_comm, ch_val)  ! below
         endif
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'And', ch_comm, ch_val)  ! below
         endif
         ! If there was 'and', save two blocks:
         if (LEN(trim(adjustl(ch_comm))) /= 0) then
            Scell%Displ(mask_num)%logical_and = .true.
            block1 = ch_comm
            block2 = ch_val
         endif

         !---------------
         ! Check if there is an 'or':
         call split_command_separator(trim(adjustl(sting_temp)), 'or', ch_comm, ch_val)  ! below
         ! Check other possible ways:
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'OR', ch_comm, ch_val)  ! below
         endif
         if (LEN(trim(adjustl(ch_comm))) == 0) then
            call split_command_separator(trim(adjustl(sting_temp)), 'Or', ch_comm, ch_val)  ! below
         endif
         ! If there was 'or', save two blocks:
         if (LEN(trim(adjustl(ch_comm))) /= 0) then
            Scell%Displ(mask_num)%logical_or = .true.
            block1 = ch_comm
            block2 = ch_val
         endif

         !---------------
         ! interprete the line:
         !if (Scell%Displ(mask_num)%logical_and .or. Scell%Displ(mask_num)%logical_or) then ! interprete 2 blocks
         if ((LEN(trim(adjustl(block1))) > 0) .and. (LEN(trim(adjustl(block2))) > 0) ) then  ! interprete 2 blocks
            ! block 1:
            call identify_section_axis(Scell, block1, mask_num, 1) ! below

            ! block 2:
            call identify_section_axis(Scell, block2, mask_num, 2) ! below

         else ! interprete single line
            ! whole line:
            call identify_section_axis(Scell, sting_temp, mask_num, 1) ! below
         endif
      enddo ! i = 1,3
      !print*, 'read_and_define_section: ', Scell%Displ(mask_num)%logical_and, Scell%Displ(mask_num)%logical_or, &
      !trim(adjustl(block1))//' ', trim(adjustl(block2))
!       print*, 'sta', mask_num, Scell%Displ(mask_num)%r_start(1, :)
!       print*, 'end', mask_num, Scell%Displ(mask_num)%r_end(1, :)
!       print*, 'sta', mask_num, Scell%Displ(mask_num)%r_start(2, :)
!       print*, 'end', mask_num, Scell%Displ(mask_num)%r_end(2, :)

   endif
!    pause 'PAUSE read_and_define_section'
end subroutine read_and_define_section


subroutine identify_section_axis(Scell, read_line, mask_num, axis_ind)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   character(*), intent(in) :: read_line  ! line with two command separated by a given symbol
   integer, intent(in) :: mask_num, axis_ind
   !-----------------
   character(100) :: command1, command2, ch_num
   integer :: Reason, len1, len2, axis_num, ind_larger, ind_smaller
   real(8) :: r_temp

   ! section along X:
   call split_command_separator(trim(adjustl(read_line)), 'x', command1, command2)  ! below
   if (LEN(trim(adjustl(command1))) == 0) then   ! another possible way of writing it
      call split_command_separator(trim(adjustl(read_line)), 'X', command1, command2)  ! below
   endif
   command1 = trim(adjustl(command1))
   command2 = trim(adjustl(command2))
   len1 = LEN(trim(adjustl(command1)))
   len2 = LEN(trim(adjustl(command2)))

   if ( (len1 /= 0) .or. (len2 /=0) ) then   ! section along X:
      axis_num = 1
      Scell%Displ(mask_num)%axis_ind(axis_ind) = axis_num  ! X

      ! Block #1:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command1)), '>')
      ind_smaller = INDEX(trim(adjustl(command1)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command1, Scell, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command1, Scell, mask_num, axis_ind, axis_num, 2)   ! below

      ! Block #2:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command2)), '>')
      ind_smaller = INDEX(trim(adjustl(command2)), '<')

      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command2, Scell, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command2, Scell, mask_num, axis_ind, axis_num, 2)   ! below
   endif

   ! section along Y:
   call split_command_separator(trim(adjustl(read_line)), 'y', command1, command2)  ! below
   if (LEN(trim(adjustl(command1))) == 0) then   ! another possible way of writing it
      call split_command_separator(trim(adjustl(read_line)), 'Y', command1, command2)  ! below
   endif
   command1 = trim(adjustl(command1))
   command2 = trim(adjustl(command2))
   len1 = LEN(trim(adjustl(command1)))
   len2 = LEN(trim(adjustl(command2)))
   if ( (len1 /= 0) .or. (len2 /=0) ) then ! section along Y
      axis_num = 2
      Scell%Displ(mask_num)%axis_ind(axis_ind) = axis_num  ! Y

      ! Block #1:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command1)), '>')
      ind_smaller = INDEX(trim(adjustl(command1)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command1, Scell, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command1, Scell, mask_num, axis_ind, axis_num, 2)   ! below

      ! Block #2:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command2)), '>')
      ind_smaller = INDEX(trim(adjustl(command2)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command2, Scell, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command2, Scell, mask_num, axis_ind, axis_num, 2)   ! below
   endif

   ! section along Z:
   call split_command_separator(trim(adjustl(read_line)), 'z', command1, command2)  ! below
   if (LEN(trim(adjustl(command1))) == 0) then   ! another possible way of writing it
      call split_command_separator(trim(adjustl(read_line)), 'Z', command1, command2)  ! below
   endif
   command1 = trim(adjustl(command1))
   command2 = trim(adjustl(command2))
   len1 = LEN(trim(adjustl(command1)))
   len2 = LEN(trim(adjustl(command2)))
   if ( (len1 /= 0) .or. (len2 /=0) ) then   ! section along Z
      axis_num = 3
      Scell%Displ(mask_num)%axis_ind(axis_ind) = axis_num  ! Z

      ! Block #1:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command1)), '>')
      ind_smaller = INDEX(trim(adjustl(command1)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command1, Scell, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command1, Scell, mask_num, axis_ind, axis_num, 2)   ! below

      ! Block #2:
      ! Check if there is an angebraic sign:
      ind_larger = INDEX(trim(adjustl(command2)), '>')
      ind_smaller = INDEX(trim(adjustl(command2)), '<')
      ! Assign the start of the axis, if it is:
      call assign_atomic_section(ind_larger, command2, Scell, mask_num, axis_ind, axis_num, 1)   ! below
      ! Assign the end of the axis, if it is:
      call assign_atomic_section(ind_smaller, command2, Scell, mask_num, axis_ind, axis_num, 2)   ! below
   endif
end subroutine identify_section_axis



subroutine assign_atomic_section(ind_given, command, Scell, mask_num, ind_sec, ind_axis, start_or_end)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: ind_given, mask_num, ind_sec, ind_axis, start_or_end
   character(*), intent(in) :: command
   !-----------------
   integer :: Reason, len1
   real(8) :: r_temp
   logical :: left_boundary

   left_boundary = .false. ! to start with

   if (ind_given > 0) then ! there is a '>' or '<'
      Reason = -1 ! to start with
      len1 = LEN(trim(adjustl(command)))

      if (ind_given == 1) then  ! it is on the left of the number
         read(command(2:len1), *, IOSTAT=Reason) r_temp
         if (start_or_end == 1) then   ! '>' it is the lower boundary
            left_boundary = .true.
         endif
      elseif (ind_given == len1) then ! it is on the right of the number
         read(command(1:len1-1), *, IOSTAT=Reason) r_temp
         if (start_or_end == 2) then   ! '<' it is the lower boundary
            left_boundary = .true.
         endif
      endif ! (ind_given == len1)
      !print*, 'r_temp', r_temp

      if (Reason == 0) then   ! read it well, use it
         if (left_boundary) then ! '>'
            Scell%Displ(mask_num)%r_start(ind_sec, ind_axis) = max ( r_temp, Scell%Displ(mask_num)%r_start(ind_sec, ind_axis) )
         else  ! '<'
            Scell%Displ(mask_num)%r_end(ind_sec, ind_axis) = min ( r_temp, Scell%Displ(mask_num)%r_end(ind_sec, ind_axis) )
         endif
      endif ! (Reason == 0)
      !print*, Scell%Displ(mask_num)%r_start(ind_sec, ind_axis), Scell%Displ(mask_num)%r_end(ind_sec, ind_axis)
   endif ! (ind_given > 0)
end subroutine assign_atomic_section


subroutine split_command_separator(read_line, separator, command1, command2, back)
   character(*), intent(in) :: read_line  ! line with two command separated by a given symbol
   character(*), intent(in) :: separator  ! separator symbol
   character(*), intent(out) :: command1  ! command #1, before separator
   character(*), intent(out) :: command2  ! command #2, after separator
   logical, intent(in), optional :: back  ! to search from the backend
   !-----------------------
   integer :: N_sep, sep_len
   logical :: back_search

   if (present(back)) then ! read what user set
      back_search = back
   else  ! by default, search from the start
      back_search = .false.
   endif

   ! Find the position of the separator:
   N_sep = 0   ! to start with
   N_sep = INDEX(read_line, separator, back=back_search) ! intrinsic
   sep_len = LEN(trim(adjustl(separator)))

   ! If separator is there, read commands:
   if (N_sep > 0) then
      command1 = read_line(1:N_sep-1)
      command2 = read_line(N_sep+sep_len:)
   else ! no separator in the string
      command1 = ''  ! undefined
      command2 = ''  ! undefined
   endif
end subroutine split_command_separator



subroutine interpret_additional_XYZ_input(temp_ch2, save_XYZ_extra)
   character(*), intent(in) :: temp_ch2
   logical, dimension(:), intent(inout) :: save_XYZ_extra
   !-------------------
   integer :: i, Nsiz
   Nsiz = LEN(trim(adjustl(temp_ch2)))
   save_XYZ_extra = .false. ! to start with, no additional input
   RAPXYZ:do i = 1, Nsiz
      !print*, 'interpret_additional_XYZ_input: ', trim(adjustl(temp_ch2(i:i)))
      select case (trim(adjustl(temp_ch2(i:i))))
      case ('!')  ! comment line starts, don't read from here on
         exit RAPXYZ
      case ('M', 'm')   ! to printout atomic mass
         save_XYZ_extra(1) = .true.
      case ('Q', 'q')   ! to printout atomic charge
         save_XYZ_extra(2) = .true.
      case ('E', 'e')   ! to printout atomic kinetic energy
         save_XYZ_extra(3) = .true.
      end select
   enddo RAPXYZ
   !print*, 'save_XYZ_extra=', save_XYZ_extra
end subroutine interpret_additional_XYZ_input



subroutine interprete_distribution_input(temp_ch, numpar, Scell, read_well)
   character(*), intent(in) :: temp_ch ! line read from the input file
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   logical, intent(inout) :: read_well
   !--------------------
   integer :: count_lines, Reason, N, Nsiz, i, Na
   real(8) :: dE, Emax, Emin, dE_min, dEa, Ea_max, dEa_out, Ea_max_out
   character(100) :: ch_temp1, ch_temp2
   logical :: read_well_at

   count_lines = 0
   read_well = .true.   ! to start with
   ! Default values:
   numpar%save_fe = .false.
   numpar%save_fe_orb = .false.
   numpar%save_fe_grid = .false.
   dE_min = 1.0d-4   ! [eV] minimal allowed grid step
   dE = 0.1d0  ! [eV] energy grid step
   Emin = -30.0d0 ! [eV] energy grid start
   Emax = 100.0d0 ! [eV] energy grid end
   numpar%save_fa = .false.

   !-----------------------------------
   ! Read electronic distribution part:
   read(temp_ch,*,IOSTAT=Reason) N, dE, Emax
   call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"

   if (.not. read_well) then ! something wrong with the user-defined grid
      ! try reading just the first flag, no grid parameters
      read(temp_ch,*,IOSTAT=Reason) N
      call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
      if (.not. read_well) then ! something wrong with the user-defined grid
         print*, 'Trouble reading the line with electronic distribution parameters, fe will not be printed out!'
         N = 0 ! by default, no printing out
      endif
   else
      if ((dE < dE_min) .and. (abs(N) == 2)) then
         dE = dE_min
         print*, 'Energy grid step for fe is too small, using the default value:', dE_min
      endif
   endif

   select case (N)
   case (1) ! printout distribution on TB energy levels
      numpar%save_fe = .true.
   case (-1)   ! printout orbital-resolved distribution
      numpar%save_fe = .true.
      numpar%save_fe_orb = .true.
   case (2) ! printout distribution on TB energy levels and on the user-defined grid
      numpar%save_fe = .true.
      numpar%save_fe_grid = .true.
   case (-2) ! printout distribution on the user-defined grid, but not on TB energy levels
      numpar%save_fe_grid = .true.
   case (3) ! printout distribution on TB energy levels, orbital-resolved, and on the user-defined grid
      numpar%save_fe = .true.
      numpar%save_fe_orb = .true.
      numpar%save_fe_grid = .true.
   case (-3) ! printout orbital-resolved distribution, and on the user-defined grid
      !numpar%save_fe = .true.   ! but not the total one, for some reason
      numpar%save_fe_orb = .true.
      numpar%save_fe_grid = .true.
   end select

   if (numpar%save_fe_grid) then
      ! Now we know the grid parameters:
      Nsiz = INT((Emax-Emin)/dE)+1
      ! Set the default grids:
      allocate(Scell%E_fe_grid(Nsiz), source=0.0d0)
      allocate(Scell%fe_on_grid(Nsiz), source=0.0d0)
      allocate(Scell%fe_high_on_grid(Nsiz), source=0.0d0)
      allocate(numpar%high_DOS(Nsiz), source=0.0d0)
      allocate(Scell%fe_norm_on_grid(Nsiz), source=0.0d0)
      allocate(Scell%fe_norm_high_on_grid(Nsiz), source=0.0d0)
      ! Create the grid:
      Scell%E_fe_grid(1) = Emin
      do i = 2, Nsiz+1
         Scell%E_fe_grid(i) = Scell%E_fe_grid(i-1) + dE
      enddo
   endif


   !-----------------------------------
   ! Read atomic distribution part, if any:
   ! check if there is a command for atomic distribution provided:
   ch_temp1 = ''  ! to start with
   ch_temp2 = ''  ! to start with
   call split_command_separator(trim(adjustl(temp_ch)), 'fa ', ch_temp1, ch_temp2)  ! below
   if (LEN(trim(adjustl(ch_temp1))) > 0) then   ! printout atomic distribution
      numpar%save_fa = .true.
      !print*, 'Atomic didstribution will be printed out', trim(adjustl(ch_temp1))//':', trim(adjustl(ch_temp2))
   else
      !print*, 'Atomic didstribution will not be printed out', trim(adjustl(ch_temp1))//':', trim(adjustl(ch_temp2))
   endif

   read_well_at = .false.  ! to start with
   if (numpar%save_fa) then
      !read(ch_temp2,*,IOSTAT=Reason) dEa, Ea_max
      read(ch_temp2,*,IOSTAT=Reason) dEa_out, Ea_max_out
      call read_file(Reason, count_lines, read_well_at)    ! module "Dealing_with_files"
   endif
   if (.not.read_well_at) then   ! use defaults
      Ea_max_out = 10.0d0  ! [eV] default value
      dEa_out = 0.01d0    ! [eV] default value
   endif

   ! Assume equidistrant grid:
   ! For internal use:
   Ea_max = 1.0d0 ! to srtart with
   Nsiz = 200
   dEa = Ea_max/dble(Nsiz)
   allocate(Scell%fa(Nsiz), source = 0.0d0)
   allocate(Scell%fa_eq(Nsiz), source = 0.0d0)
   allocate(Scell%fa_pot(Nsiz), source = 0.0d0)
   allocate(Scell%fa_tot(Nsiz), source = 0.0d0)
   allocate(Scell%fa_eq_pot(Nsiz), source = 0.0d0)
   allocate(Scell%Ea_grid(Nsiz))
   allocate(Scell%Ea_pot_grid(Nsiz))
   allocate(Scell%Ea_tot_grid(Nsiz))
   ! Set the grid:
   Scell%Ea_grid(1) = 0.0d0 ! starting point
   do i = 2, Nsiz
      Scell%Ea_grid(i) = Scell%Ea_grid(i-1) + dEa
      !print*, i, Scell%Ea_grid(i)
   enddo ! i
   Scell%Ea_pot_grid(:) = Scell%Ea_grid(:) - 10.0d0   ! to start with
   Scell%Ea_tot_grid(:) = Scell%Ea_grid(:) - 10.0d0   ! to start with

   ! For printout:
   Nsiz = INT(Ea_max_out/dEa_out)+1
   allocate(Scell%fa_out(Nsiz), source = 0.0d0)
   allocate(Scell%fa_eq_out(Nsiz), source = 0.0d0)
   allocate(Scell%fa_pot_out(Nsiz), source = 0.0d0)
   allocate(Scell%fa_tot_out(Nsiz), source = 0.0d0)
   allocate(Scell%fa_eq_pot_out(Nsiz), source = 0.0d0)
   allocate(Scell%Ea_grid_out(Nsiz))
   allocate(Scell%Ea_pot_grid_out(Nsiz))
   allocate(Scell%Ea_tot_grid_out(Nsiz))
   ! Set the grid:
   Scell%Ea_grid_out(1) = 0.0d0 ! starting point
   do i = 2, Nsiz
      Scell%Ea_grid_out(i) = Scell%Ea_grid_out(i-1) + dEa_out
      !print*, i, Scell%Ea_grid_out(i)
   enddo ! i
   Scell%Ea_pot_grid_out = Scell%Ea_grid_out - 10.0d0  ! to start with
   Scell%Ea_tot_grid_out = Scell%Ea_grid_out - 10.0d0  ! to start with
   !pause 'interprete_distribution_input'
end subroutine interprete_distribution_input



subroutine set_Bath_grid_electrons(numpar, read_well_out, Error_descript)
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   logical, intent(inout) :: read_well_out
   character(*), intent(inout) :: Error_descript
   !-----------------------------------------
   character(200) :: Path, full_file_name
   logical :: file_exist, read_well
   integer :: FN, Nsiz, count_lines, Reason, i

   read_well_out = .false.   ! to start with
   numpar%Transport_e = .false. ! to start with

   Path = trim(adjustl(m_INPUT_directory))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(numpar%El_bath_step_grid_file))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(numpar%El_bath_reset_grid)) deallocate(numpar%El_bath_reset_grid) ! make sure it's possible to allocate
      allocate(numpar%El_bath_reset_grid(Nsiz)) ! allocate it
      if (allocated(numpar%El_bath_grid_Ta)) deallocate(numpar%El_bath_grid_Ta) ! make sure it's possible to allocate
      allocate(numpar%El_bath_grid_Ta(Nsiz)) ! allocate it
      if (allocated(numpar%El_bath_grid_tau)) deallocate(numpar%El_bath_grid_tau) ! make sure it's possible to allocate
      allocate(numpar%El_bath_grid_tau(Nsiz)) ! allocate it

      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
         read(FN,*,IOSTAT=Reason) numpar%El_bath_reset_grid(i), numpar%El_bath_grid_Ta(i), numpar%El_bath_grid_tau(i)
         call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
         if (.not. read_well) then ! something wrong with the user-defined grid
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(numpar%El_bath_step_grid_file))//' could not read line ', count_lines
            goto 9997  ! couldn't read the data, exit the cycle
         endif
      enddo
      numpar%i_El_bath_dt = 1 ! to start from
      read_well_out = .true.     ! we read the grid from the file well

9997 call close_file('close', FN=FN) ! module "Dealing_with_files"
   endif ! file_exist
end subroutine set_Bath_grid_electrons


subroutine set_Bath_grid_atoms(numpar, read_well_out, Error_descript)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   logical, intent(inout) :: read_well_out
   character(*), intent(inout) :: Error_descript
   !-----------------------------------------
   character(200) :: Path, full_file_name
   logical :: file_exist, read_well
   integer :: FN, Nsiz, count_lines, Reason, i

   read_well_out = .false.   ! to start with
   numpar%Transport = .false. ! to start with

   Path = trim(adjustl(m_INPUT_directory))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(numpar%At_bath_step_grid_file))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(numpar%At_bath_reset_grid)) deallocate(numpar%At_bath_reset_grid) ! make sure it's possible to allocate
      allocate(numpar%At_bath_reset_grid(Nsiz)) ! allocate it
      if (allocated(numpar%At_bath_grid_Ta)) deallocate(numpar%At_bath_grid_Ta) ! make sure it's possible to allocate
      allocate(numpar%At_bath_grid_Ta(Nsiz)) ! allocate it
      if (allocated(numpar%At_bath_grid_tau)) deallocate(numpar%At_bath_grid_tau) ! make sure it's possible to allocate
      allocate(numpar%At_bath_grid_tau(Nsiz)) ! allocate it

      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
         read(FN,*,IOSTAT=Reason) numpar%At_bath_reset_grid(i), numpar%At_bath_grid_Ta(i), numpar%At_bath_grid_tau(i)
         call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
         if (.not. read_well) then ! something wrong with the user-defined grid
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(numpar%At_bath_step_grid_file))// &
                              ' could not read line ', count_lines
            goto 9998  ! couldn't read the data, exit the cycle
         endif
      enddo
      numpar%i_At_bath_dt = 1 ! to start from
      read_well_out = .true.     ! we read the grid from the file well

9998 call close_file('close', FN=FN) ! module "Dealing_with_files"
   endif ! file_exist
end subroutine set_Bath_grid_atoms



subroutine read_user_additional_data(FN, count_lines, matter, Scell, NSC, user_data)
   integer, intent(in) :: FN  ! file to read from
   integer, intent(inout) :: count_lines  ! line we are reading
   type(Solid), intent(in) :: matter   ! all material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(User_overwrite_data), intent(inout) :: user_data ! data to be read
   !--------------------------
   logical :: read_well
   integer :: Reason
   character(100) :: text

   ! Default values:
   user_data%do_overwrite = .false. ! to start with, no data to use

   read_well = .true.   ! to start with
   RD: do while (read_well)
      read(FN,*,IOSTAT=Reason) text
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) exit RD  ! end of file, stop reading

      call interpret_user_data(FN, count_lines, read_well, text, matter, Scell, NSC, user_data) ! below
      if (.not. read_well) exit RD  ! end of file, stop reading
   enddo RD

end subroutine read_user_additional_data


subroutine interpret_user_data(FN, count_lines, read_well, text, matter, Scell, NSC, user_data)
   integer, intent(in) :: FN  ! file to read from
   integer, intent(inout) :: count_lines  ! line we are reading
   logical, intent(inout) :: read_well ! marker if the line read well
   character(*), intent(in) :: text ! what was read in the previous line
   type(Solid), intent(in) :: matter	! all material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(User_overwrite_data), intent(inout) :: user_data ! data to be read
   !-------------------------
   real(8) :: temp
   integer :: Nat, Nsh, i, temp_int(2), Reason
   character(3) :: temp_ch
   character(150) :: Error_descript

   ! Nat = size(matter%Atoms) ! number of atoms is unknown yet, use just some large value here:
   Nat = 100
   ! Find the maximal number of shells:
   Nsh = 0  ! to start with
   do i = 1, Nat
      !Nsh = max(Nsh, size(matter%Atoms(i)%Ip)) ! numer of shells is undefined yet!
      Nsh = 50
   enddo

   select case (trim(adjustl(text)))
   !--------------------------------
   case ('DRUDE', 'Drude', 'drude')
      ! initial n and k of unexcited material (used for DRUDE model only!):
      read(FN,*,IOSTAT=Reason) Scell(NSC)%eps%n, Scell(NSC)%eps%k	! initial n and k coeffs
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         print*, 'Could not read line (#1)', count_lines, ' under "DRUDE" option'
         print*, 'Default values are used: n,k=', Scell(1)%eps%n, Scell(1)%eps%k
      endif

      ! [me] effective mass of CB electron and VB hole:
      read(FN,*,IOSTAT=Reason) Scell(NSC)%eps%me_eff, Scell(NSC)%eps%mh_eff
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         print*, 'Could not read line (#2)', count_lines, ' under "DRUDE" option'
         print*, 'Default values are used: me,mh=', Scell(1)%eps%me_eff, Scell(1)%eps%mh_eff
      else
         Scell(NSC)%eps%me_eff = Scell(NSC)%eps%me_eff*g_me	! [kg]
         Scell(NSC)%eps%mh_eff = Scell(NSC)%eps%mh_eff*g_me	! [kg]
      endif

      ! [fs] mean scattering times of electrons and holes:
      read(FN,*,IOSTAT=Reason) Scell(NSC)%eps%tau_e, Scell(NSC)%eps%tau_h
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         print*, 'Could not read line (#3)', count_lines, ' under "DRUDE" option'
         print*, 'Default values are used: te,th=', Scell(1)%eps%tau_e, Scell(1)%eps%tau_h
      endif

   !--------------------------------
   case ('NAME', 'Name', 'name')
      read(FN,*,IOSTAT=Reason) temp_int(1), temp_ch
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         if (.not.allocated(user_data%name)) then
            allocate(user_data%name(Nat))
            user_data%name = ''  ! to start with
         endif
         user_data%name(temp_int(1)) = temp_ch
         user_data%do_overwrite = .true.
      endif

   !--------------------------------
   case ('MASS', 'Mass', 'mass')
      read(FN,*,IOSTAT=Reason) temp_int(1), temp
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         if (.not.allocated(user_data%mass)) allocate(user_data%mass(Nat), source = -1.0d0)
         user_data%mass(temp_int(1)) = temp
         user_data%do_overwrite = .true.
      endif

   !--------------------------------
   case ('NO_AUGER', 'No_Auger', 'no_auger', 'Exclude_Auger', 'exclude_auger', 'EXCLUDE_AUGER')
      if (.not.allocated(user_data%auger)) allocate(user_data%auger(Nat,Nsh), source = 1.0d30)
      user_data%auger = 1.0d30   ! exclude Auger by setting time to infinity
      user_data%do_overwrite = .true.

   !--------------------------------
   case ('AUGER', 'Auger', 'auger')
      read(FN,*,IOSTAT=Reason) temp_int(:), temp
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         if (.not.allocated(user_data%auger)) allocate(user_data%auger(Nat,Nsh), source = -1.0d30)
         user_data%auger(temp_int(1),temp_int(2)) = temp
         user_data%do_overwrite = .true.
      endif
!       print*, 'AUGER', temp_int(1), temp_int(2), user_data%auger(temp_int(1),temp_int(2))
!       print*, 'SIZE:', size(user_data%auger,1), size(user_data%auger,2), Nat, Nsh
!       print*, 'ALL:', allocated(user_data%auger), user_data%auger(:,:)
!       pause 'AUGER'

   !--------------------------------
   case ('IP', 'Ip', 'ip', 'ionization_potential', 'Ionization_potential', 'IONIZATION_POTENTIAL')
      read(FN,*,IOSTAT=Reason) temp_int(:), temp
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         if (.not.allocated(user_data%Ip)) allocate(user_data%Ip(Nat,Nsh), source = -1.0d30)
         user_data%Ip(temp_int(1),temp_int(2)) = temp
         user_data%do_overwrite = .true.
      endif

   !--------------------------------
   case ('Ne', 'ne', 'NE', 'NUMBER_OF_ELECTRONS', 'Number_of_electrons', 'Population')
      read(FN,*,IOSTAT=Reason) temp_int(:), temp
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         if (.not.allocated(user_data%Ne_shell)) allocate(user_data%Ne_shell(Nat,Nsh), source = -1.0d0)
         user_data%Ne_shell(temp_int(1),temp_int(2)) = temp
         user_data%do_overwrite = .true.
      endif

   !--------------------------------
   case ('Ek', 'ek', 'EK', 'Kinetic_energy', 'KINETIC_ENERGY', 'kinetic_energy')
      read(FN,*,IOSTAT=Reason) temp_int(:), temp
      call read_file(Reason, count_lines, read_well)
      if (read_well) then
         if (.not.allocated(user_data%Ek)) allocate(user_data%Ek(Nat,Nsh), source = -1.0d30)
         user_data%Ek(temp_int(1),temp_int(2)) = temp
         user_data%do_overwrite = .true.
      endif
   end select

end subroutine interpret_user_data


subroutine overwrite_atomic_data(user_data, matter)
   type(User_overwrite_data), intent(in) :: user_data ! data to be read
   type(Solid), intent(inout) :: matter	! all material parameters
   !----------------------------
   integer :: i, Nat, j, Nsh

   ! Check if the user overwrite atomic data:
   if (user_data%do_overwrite) then
      Nat = size(matter%Atoms) ! number of atoms

      ! Element name:
      if (allocated(user_data%Name)) then
         do i = 1, Nat   ! check for all elements
            if (LEN(trim(adjustl(user_data%Name(i)))) > 0) then ! overwrite the name
               matter%Atoms(i)%Name = user_data%Name(i)
            endif
         enddo
      endif

      ! Mass:
      if (allocated(user_data%mass)) then
         do i = 1, Nat   ! check for all elements
            if (user_data%mass(i) > 0.0d0) then ! overwrite the mass
               matter%Atoms(i)%Ma = user_data%mass(i) * g_amu ! [kg]
            endif
         enddo
      endif

      ! Auger decay times:
      if (allocated(user_data%auger)) then
         do i = 1, Nat
            Nsh = size(matter%Atoms(i)%Ip)   ! number of shells
            do j = 1, Nsh
               if (user_data%Auger(i,j) > 0.0d0) then ! overwrite the auger decay time
                  matter%Atoms(i)%Auger(j) = user_data%auger(i,j) ! [fs]
               endif
            enddo ! j
         enddo ! i
      endif

      ! Ionization potentials:
      if (allocated(user_data%Ip)) then
         do i = 1, Nat
            Nsh = size(matter%Atoms(i)%Ip)   ! number of shells
            do j = 1, Nsh
               if (user_data%Ip(i,j) > 0.0d0) then ! overwrite the ionization potential
                  matter%Atoms(i)%Ip(j) = user_data%Ip(i,j) ! [eV]
               endif
            enddo ! j
         enddo ! i
      endif

      ! Kinetic energies of electronic in atomic shells:
      if (allocated(user_data%Ip)) then
         do i = 1, Nat
            Nsh = size(matter%Atoms(i)%Ek)   ! number of shells
            do j = 1, Nsh
               if (user_data%Ek(i,j) > 0.0d0) then ! overwrite the kinetic energy
                  matter%Atoms(i)%Ek(j) = user_data%Ek(i,j) ! [eV]
               endif
            enddo ! j
         enddo ! i
      endif

      ! Number of electrons in the atomic shell:
      if (allocated(user_data%Ip)) then
         do i = 1, Nat
            Nsh = size(matter%Atoms(i)%Ek)   ! number of shells
            do j = 1, Nsh
               if (user_data%Ne_shell(i,j) >= 0.0d0) then ! overwrite the electron population
                  matter%Atoms(i)%Ne_shell(j) = user_data%Ne_shell(i,j)
               endif
            enddo ! j
         enddo ! i
      endif

   endif ! (user_data%do_overwrite)
end subroutine overwrite_atomic_data


subroutine set_MD_step_grid(File_name, numpar, read_well_out, Error_descript)
   character(*), intent(in) :: File_name    ! file name with input data (or timestep of MD [fs])
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   logical, intent(inout) :: read_well_out
   character(*), intent(inout) :: Error_descript
   !-----------------------------------------
   character(200) :: Path, full_file_name
   logical :: file_exist, read_well
   integer :: FN, Nsiz, count_lines, Reason, i

   read_well_out = .false.   ! to start with
   Path = trim(adjustl(m_INPUT_directory))//numpar%path_sep    ! where to find the file with the data
   full_file_name = trim(adjustl(Path))//trim(adjustl(File_name))   ! to read the file from the INPUT_DATA directory
   inquire(file=trim(adjustl(full_file_name)),exist=file_exist) ! check if input file is there
   !print*, File_name, file_exist

   if (file_exist) then ! try to read it, if there is a grid provided
      open(newunit = FN, FILE = trim(adjustl(full_file_name)), status = 'old', action='read')
      ! Find the grid size from the file:
      call Count_lines_in_file(FN, Nsiz) ! module "Dealing_with_files"
      ! Knowing the size, create the grid-array and read from the file:
      if (allocated(numpar%dt_MD_reset_grid)) deallocate(numpar%dt_MD_reset_grid) ! make sure it's possible to allocate
      allocate(numpar%dt_MD_reset_grid(Nsiz)) ! allocate it
      if (allocated(numpar%dt_MD_grid)) deallocate(numpar%dt_MD_grid) ! make sure it's possible to allocate
      allocate(numpar%dt_MD_grid(Nsiz)) ! allocate it

      ! Read data on the grid from the file:
      count_lines = 0   ! just to start counting lines in the file
      do i = 1, Nsiz    ! read grid line by line from the file
         read(FN,*,IOSTAT=Reason) numpar%dt_MD_reset_grid(i), numpar%dt_MD_grid(i)    ! grid data from the file
         call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
         if (.not. read_well) then ! something wrong with the user-defined grid
            write(Error_descript,'(a,i3)') 'In the file '//trim(adjustl(File_name))//' could not read line ', count_lines
            goto 9993  ! couldn't read the data, exit the cycle
         endif
         !print*, i, numpar%dt_MD_reset_grid(i), numpar%dt_MD_grid(i)
      enddo
      read_well_out = .true.    ! we read the grid from the file well
      numpar%i_dt = 0   ! to start with
      numpar%dt = numpar%dt_MD_grid(1)   ! to start with

!       print*, 'set_MD_step_grid:'
!       print*, numpar%dt_MD_reset_grid(:)
!       print*, 'dt=', numpar%dt_MD_grid(:)

9993 call close_file('close', FN=FN) ! module "Dealing_with_files"
   else ! If there is no input file, check if the teimstep is provided instead:
      count_lines = 0   ! just to start counting lines in the file
      read(File_name,*,IOSTAT=Reason) numpar%dt  ! Time step for MD [fs]
      call read_file(Reason, count_lines, read_well)    ! module "Dealing_with_files"
      if (read_well) read_well_out = .true.    ! we read the grid from the file well
      numpar%i_dt = -1   ! to mark that the reset option is unused
      ! Set often-used values:
   endif ! file_exist
   numpar%halfdt = numpar%dt/2.0d0           ! dt/2, often used
   numpar%dtsqare = numpar%dt*numpar%halfdt  ! dt*dt/2, often used
   numpar%dt3 = numpar%dt**3/6.0d0           ! dt^3/6, often used
   numpar%dt4 = numpar%dt*numpar%dt3/8.0d0   ! dt^4/48, often used

   !pause 'set_MD_step_grid'
end subroutine set_MD_step_grid



subroutine read_input_material(File_name, Scell, matter, numpar, laser, user_data, Err)
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   character(*), intent(in) :: File_name
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(User_overwrite_data), intent(inout) :: user_data
   type(Error_handling), intent(inout) :: Err	! error save
   real(8) read_var(3) ! just to read variables from file
   integer FN, N, Reason, count_lines, i
   logical file_opened, read_well
   character(100) Error_descript, text
   character(500) read_line

   !open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), status = 'old')
   FN=109
   open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      call Save_error_details(Err, 2, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3417
   endif

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Material parameters:
   count_lines = 0   ! to start with

   ! Material name (and possibly the file name with coordinates):
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) matter%Name, numpar%Cell_filename
   if (Reason /= 0) then ! try to read just single variable:
      numpar%Cell_filename = ''  ! nullify it
      !read(FN,*,IOSTAT=Reason) matter%Name
      read(read_line,*,IOSTAT=Reason) matter%Name
      if (numpar%verbose) write(*,'(a)') 'No valid filename with coordinates provided, assuming default'
   else
      ! Make sure that, if there is a path, the separator is correct:
      call ensure_correct_path_separator(numpar%Cell_filename, numpar%path_sep) ! module "Dealing_with_files"
      ! Check if such a file exists:
      call check_coordinates_filename(numpar%Cell_filename, numpar%verbose) ! see below
   endif
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
   if (Err%Err) goto 3417
   !print*, trim(adjustl(matter%Name)), ' : ', trim(adjustl(numpar%Cell_filename))

   ! chemical formula of the compound (used in MC in case of EADL parameters):
   read(FN,*,IOSTAT=Reason) matter%Chem
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
   if (Err%Err) goto 3417

   if (.not.allocated(Scell)) allocate(Scell(1)) ! just for start, 1 supercell
   do i = 1, size(Scell) ! for all supercells

      ! initial electron temperature [K] or filename with the distribution function:
      read(FN, '(a)', IOSTAT=Reason) read_line
      read(read_line,*,IOSTAT=Reason) Scell(i)%Te
      if (Reason == 0) then ! there was a number, interpret it as electronic temperature
         call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                                 add_error_info='Line: '//trim(adjustl(read_line))) !below
         if (Err%Err) goto 3417
         numpar%fe_filename = '' ! user provided no filename for distribution (use Fermi with given tempreature instead)
      else ! maybe there was a name of the file with electronic distribution:
         Scell(i)%Te = -1.0d0 ! Just to indicate nonequilibrium distribution
         read(read_line,*,IOSTAT=Reason) numpar%fe_filename  ! read filename
         call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                                 add_error_info='Line: '//trim(adjustl(read_line))) !below
         if (Err%Err) goto 3417
      endif
      Scell(i)%TeeV = Scell(i)%Te/g_kb ! [eV] electron temperature
      ! Printout warning if electron temperature is too high:
      if (Scell(i)%TeeV >= 5.0) then
         write(text,'(f16.3)',IOSTAT=Reason) Scell(i)%TeeV
         call printout_warning(6, 2, text_to_add=trim(adjustl(text)) ) ! below
      endif

      ! initial atomic temperature [K]:
      read(FN, '(a)', IOSTAT=Reason) read_line
      read(read_line,*,IOSTAT=Reason) Scell(i)%Ta
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
      if (Err%Err) goto 3417
      Scell(i)%TaeV = Scell(i)%Ta/g_kb ! [eV] atomic temperature
      ! Printout warning if atomic temperature is too high:
      if (Scell(i)%TaeV >= 1.0) then
         write(text,'(f16.3)',IOSTAT=Reason) Scell(i)%TaeV
         call printout_warning(6, 3, text_to_add=trim(adjustl(text)) ) ! below
      endif

   enddo !Scell


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Basic parameters of the simulation:

   ! Start of the simulation [fs]:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%t_start
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
   if (Err%Err) goto 3417

   ! total duration of simulation [fs]:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) numpar%t_total
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
   if (Err%Err) goto 3417

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Laser parameters:

   ! Number of laser pulses to model:
   read(FN, '(a)', IOSTAT=Reason) read_line
   read(read_line,*,IOSTAT=Reason) N          ! How many pulses
   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
   if (Err%Err) goto 3417

   ! For multiple pulses, parameters for each:
   PULS:if (N >= 0) then ! If there is at least one pulse:
    if (allocated(laser)) deallocate(laser)
    allocate(laser(N))  ! that's how many pulses
    do i = 1, N         ! read parameters for all pulses
      text = ''   ! to start with
      read_var = -1.0d0 ! to start with

      ! Read the laser pulse fluence or dose:
      read(FN, '(a)', IOSTAT=Reason) read_line  ! to interprete the text below
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
      if (Err%Err) goto 3417

      ! Check if there are specifications provided, or just the number (assuming default):
      if (it_is_number(read_line)) then ! function from module "Little_subroutines"
         read(read_line,*,IOSTAT=Reason) read_var(:)
         if (numpar%verbose) write(*,'(a, f, f, f)') &
               'No specification of fluence provided, assuming absorbed dose [eV/atom]', read_var(:)
      else  ! read the specifications of the input (fluence vs dose etc.):
         read(read_line,*,IOSTAT=Reason) text, read_var(:)
         if (numpar%verbose) write(*,'(a, f, f, f)') &
               'Specification of fluence provided: '//trim(adjustl(text)), read_var(:)
      endif
      ! Check if at least one number was read well:
      if (read_var(1) < 0.0d0) then ! could not read the numberes, something went wrong
         write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript)), read_var(:)
         print*, 'Line: ', trim(adjustl(read_line))
         goto 3417
      endif
      ! Check if there was single fluence, or a grid:
      if (Reason == 0) then ! three numbers are given, set few inputs with varying dose:
         ! For this simulation run, use the first value given:
         laser(i)%F = read_var(1)   ! ABSORBED DOSE IN [eV/atom]
         ! For other simulation runs, create the input files accordingly:
         call prepare_multiple_inputs(numpar, File_name, read_var)  ! below

      else ! probably, there weren't three parameters in the line, so read just the fluence
         laser(i)%F = read_var(1)   ! ABSORBED DOSE IN [eV/atom]
      endif

      ! Check if there are additional pulse specifications:
      call check_pulse_specifications(trim(adjustl(text)), laser(i), read_var(1))   ! below
      ! Printout warning if absorbed dose is too high:
      if (laser(i)%F >= 10.0) then
         write(text,'(f16.3)',IOSTAT=Reason) laser(i)%F
         call printout_warning(6, 4, text_to_add=trim(adjustl(text)) ) ! below
      endif


      read(FN, '(a)', IOSTAT=Reason) read_line
      ! Few options for input of the photon energy or wavelength:
      ! May contain photon energy or wavelength of the laser, may contain units, may contain FWHM of spectral distribuition and units:
      call get_photon_parameters(read_line, laser, i, count_lines, File_name, Err) ! below
      if (Err%Err) goto 3417
      ! Printout warning if photon energy is too high:
      if (laser(i)%hw >= 1.0d5) then
         write(text,'(f16.3)',IOSTAT=Reason) laser(i)%hw
         call printout_warning(6, 1, text_to_add=trim(adjustl(text)) ) ! below
      endif


      read(FN, '(a)', IOSTAT=Reason) read_line
      read(read_line,*,IOSTAT=Reason) laser(i)%t	  ! PULSE FWHM-DURATION IN [fs]
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
      if (Err%Err) goto 3417

      read(FN, '(a)', IOSTAT=Reason) read_line
      read(read_line,*,IOSTAT=Reason) laser(i)%KOP ! type of pulse: 0=rectangular, 1=Gaussian, 2=SASE
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, add_error_info='Line: '//trim(adjustl(read_line))) !below
      if (Err%Err) goto 3417

      read(FN, '(a)', IOSTAT=Reason) read_line
      read(read_line,*,IOSTAT=Reason) laser(i)%t0  ! POSITION OF THE MAXIMUM OF THE PULSE IN [fs]
      call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                              add_error_info='Line: '//trim(adjustl(read_line)))  ! below
      if (Err%Err) goto 3417

      if (laser(i)%KOP .EQ. 1) laser(i)%t = laser(i)%t/2.35482	! make a gaussian parameter out of it
     enddo ! have read parameters for all pulses
   endif PULS
   

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Additional parameters:

   ! Check if there are additional options provided:
   read_well = .true.   ! to start with
   RDID: do while (read_well)
      !read(FN,*,IOSTAT=Reason) text
      read_line = ''
      ! First, read it as unformatted:
      read(FN, *, IOSTAT=Reason) read_line
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) exit RDID  ! end of file, stop reading

      ! If there is anyting in this line, try reading it as formatted:
      if ( LEN(trim(adjustl(read_line))) > 0 ) then
         backspace(FN)  ! reread the same line
         count_lines = count_lines - 1
         read(FN, '(a)', IOSTAT=Reason) read_line
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) exit RDID  ! end of file, stop reading
      endif

      read(read_line,*,IOSTAT=Reason) text
      count_lines = count_lines-1   ! reread the same line
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) exit RDID  ! end of file, stop reading

      ! Check if additional INPUT options are provided:
      call interpret_user_data_INPUT(FN, trim(adjustl(File_name)), count_lines, read_line, Scell, numpar, Err) ! below
      if (Err%Err) exit RDID  ! end of file, stop reading

      ! Check if optional numerical or model data provided by the user (e.g., to overwrite default atomic data):
      call interpret_user_data(FN, count_lines, read_well, text, matter, Scell, 1, user_data) ! below
      if (.not. read_well) exit RDID  ! end of file, stop reading

      ! Check if the numerical parameters are provided in this file:
      select case (trim(adjustl(text)))
      case ('NUMPAR', 'NumPar', 'numpar', &
            'NUMERICAL_PARAMETERS', 'Numerical_Parameters', 'numerical_parameters', &
            'NUMERICS', 'Numerics', 'numerics')
         numpar%numpar_in_input = .true.  ! mark that the num.parameters were read from here
         call read_numerical_parameters(trim(adjustl(File_name)), matter, numpar, laser, Scell, &
                                       user_data, Err, add_data=.false., count_lines_in=count_lines) ! below
         if (Err%Err) exit RDID  ! end of file, stop reading
      end select
   enddo RDID

   if (numpar%verbose) call print_time_step('Verbose option is on, XTANT is going to be a chatterbox', msec=.true.) ! modlue "Little_subroutines"
   ! Close this file, it has been read through:
3417  if (file_opened) close(FN)
end subroutine read_input_material



subroutine get_photon_parameters(read_line, laser, i, count_lines, File_name, Err)
   character(*), intent(in) :: read_line  ! line read from the input file
   type(Pulse), dimension(:), intent(inout) :: laser ! Laser pulse parameters
   integer, intent(in) :: i ! laser pulse index
   integer, intent(inout) :: count_lines  ! conter in which line the error appeared
   character(*), intent(in) :: File_name  ! which file we are reading now
   type(Error_handling), intent(inout) :: Err   ! error save
   !--------------------
   integer :: Reason
   real(8) :: temp1, temp2, temp2_rel
   character(12) :: ch_temp1, ch_temp2

   ! A few options to set input:
   !--------------------
   ! 1) lets try to read 4 variables: photon energy/wavelength, its units, FWHM spectral width, its units:
   read(read_line,*,IOSTAT=Reason) temp1, ch_temp1, temp2, ch_temp2
   if (Reason == 0) then ! 4 variables in order, interprete them:
      ! Interprete units of the mean:
      call photon_units(ch_temp1, temp1, laser(i)%hw, .true.) ! below
      ! Construct the spread (FWHM):
      if (temp2 <= 0.0d0) then
         temp2_rel = 0.0d0
      elseif (trim(adjustl(ch_temp2(1:1))) == '%') then  ! given in relative units
         temp2_rel = abs(temp2) * 1.0d-2
      else  ! given in absolute units -> get it in relative units
         temp2_rel = abs(temp2/temp1)
      endif
      laser(i)%FWHM_hw = temp2_rel * laser(i)%hw

      return ! done, can exit the subroutine
   endif

   !--------------------
   ! 2) lets try to read 3 variable: photon energy/wavelength, FWHM, units (same for both):
   read(read_line,*,IOSTAT=Reason) temp1, temp2, ch_temp1
   !print*, 'get_photon_parameters:', temp1, temp2, ch_temp1
   if (Reason == 0) then ! 3 variables in order, interprete them:
      ! Construct the spread (FWHM):
      if (temp2 <= 0.0d0) then
         temp2_rel = 0.0d0
         laser(i)%hw = temp1  ! assume photon energy is given in [eV]
      elseif (trim(adjustl(ch_temp1(1:1))) == '%') then  ! given in relative units
         laser(i)%hw = temp1  ! assume photon energy is given in [eV]
         temp2_rel = abs(temp2) * 1.0d-2
      else  ! given in absolute units -> get it in relative units
         ! Interprete units of the mean:
         call photon_units(ch_temp1, temp1, laser(i)%hw, .true.) ! below
         temp2_rel = abs(temp2/temp1)
      endif
      laser(i)%FWHM_hw = temp2_rel * laser(i)%hw

      return ! done, can exit the subroutine
   endif

   !--------------------
   ! 3) lets try to read 3 variable: photon energy/wavelength, units, FWHM (assuming the same units as in photon energy/wavelength):
   read(read_line,*,IOSTAT=Reason) temp1, ch_temp1, temp2
   if (Reason == 0) then  ! 3 variables in order, interprete them:
      ! Interprete units of the mean:
      call photon_units(ch_temp1, temp1, laser(i)%hw, .true.) ! below
      ! Construct the spread (FWHM):
      if (temp2 <= 0.0d0) then
         temp2_rel = 0.0d0
      elseif (trim(adjustl(ch_temp1(1:1))) == '%') then  ! given in relative units
         temp2_rel = abs(temp2) * 1.0d-2
      else  ! given in absolute units -> get it in relative units
         temp2_rel = abs(temp2/temp1)
      endif
      laser(i)%FWHM_hw = temp2_rel * laser(i)%hw

      return ! done, can exit the subroutine
   endif

   !--------------------
   ! 4) lets try to read 2 variables: photon energy/wavelength, FWHM (both assuming eV)
   read(read_line,*,IOSTAT=Reason) temp1, temp2
   if (Reason == 0) then  ! 3 variables in order, interprete them:
      ! Interprete units of the mean:
      laser(i)%hw = temp1
      ! The spread (FWHM):
      laser(i)%FWHM_hw = abs(temp2)

      return ! done, can exit the subroutine
   endif

   !--------------------
   ! 5) lets try to read 2 variables: photon energy/wavelength and its units (assuming monochromatic pulse, no FWHM)
   read(read_line,*,IOSTAT=Reason) temp1, ch_temp1
   if (Reason == 0) then  ! 3 variables in order, interprete them:
      ! Interprete units of the mean:
      call photon_units(ch_temp1, temp1, laser(i)%hw, .true.) ! below
      ! No spread (FWHM):
      laser(i)%FWHM_hw = 0.0d0
      return ! done, can exit the subroutine
   endif

   !--------------------
   ! 6) lets try to read 1 variable: photon energy (assuming [eV], monochromatic pulse, no FWHM)
   read(read_line,*,IOSTAT=Reason) laser(i)%hw
   if (Reason == 0) then  ! 3 variables in order, interprete them:
      ! No spread (FWHM):
      laser(i)%FWHM_hw = 0.0d0
      return ! done, can exit the subroutine
   endif

   !--------------------
   ! After all options, check if it read well:
2121   call check_if_read_well(Reason, count_lines, trim(adjustl(File_name)), Err, &
                           add_error_info='Line: '//trim(adjustl(read_line)))  ! below

end subroutine get_photon_parameters



subroutine photon_units(ch_temp, var_in, var_out, print_message, ch_temp_alt)
   character(*), intent(in) :: ch_temp   ! units to be interpred
   real(8), intent(in) :: var_in    ! variable to be renormalized according to the units
   real(8), intent(out) :: var_out  ! variable renormalized according to the units
   logical, intent(in) :: print_message
   character(*), intent(in), optional :: ch_temp_alt   ! alternative units that may be used if the main ones fail
   !--------------------

   call interprete_photon_units(ch_temp, var_in, var_out, print_message) ! below

   ! If relative units are used '%' then read the absolute units:
   if ((var_out < 0.0d0) .and. present(ch_temp_alt)) then
      call interprete_photon_units(ch_temp_alt, var_in, var_out, print_message)   ! below
   endif

   var_out = abs(var_out)
end subroutine photon_units


subroutine interprete_photon_units(ch_temp, var_in, var_out, print_message)
   character(*), intent(in) :: ch_temp    ! units to be interpred
   real(8), intent(in) :: var_in
   real(8), intent(out) :: var_out
   logical, intent(in) :: print_message
   !--------------------
   real(8) :: nm, coef_out, f

   !print*, 'interprete_photon_units:', ch_temp, var_in, var_out
   if (trim(adjustl(ch_temp(1:1))) == '!') then ! it is a comment line, not units
      var_out = -var_in
      return   ! no need to continue
   endif


   select case (trim(adjustl(ch_temp)))
   case ('eV', 'EV', 'ev') ! [eV]
      var_out = var_in

   case ('keV', 'KEV', 'kev') ! -> [eV]
      coef_out = 1.0d3
      var_out = coef_out * var_in

   case ('RY', 'Ry', 'ry') ! -> [eV]
      coef_out = g_Ry
      var_out = coef_out * var_in

   case ('au', 'a.u.', 'AU', 'A.U.') ! -> [eV]
      coef_out = g_au2ev
      var_out = coef_out * var_in

   case ('nm', 'Nm', 'NM') ! wavelength -> [eV]
      ! Convertion coefficient from 1 nm to eV:
      nm = 1.0d0
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('A', 'a', 'angstrom', 'Angstrom') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 0.1d0
      ! Convertion coefficient from [nm] to [eV]:
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('cm', 'CM', 'Cm', 'cantimeter') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 1.0d7
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('m', 'M', 'meter') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 1.0d9
      ! Convertion coefficient from [nm] to [eV]:
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('mm', 'Mm', 'MM', 'millimeter') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 1.0d6
      ! Convertion coefficient from [nm] to [eV]:
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('mkm', 'Mkm', 'MKM', 'micron', 'Micron', 'MICRON') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 1.0d3
      ! Convertion coefficient from [nm] to [eV]:
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('pm', 'PM', 'Pm', 'picometer') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 1.0d-3
      ! Convertion coefficient from [nm] to [eV]:
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"
   case ('a0', 'A0', 'Bohr', 'BOHR', 'bohr') ! wavelength -> [eV]
      ! Convert wavelength into [nm]:
      nm = 0.1d0 * g_au2A  ! a0 -> nm
      ! Convertion coefficient from [nm] to [eV]:
      var_out = convert_wavelength_to_hw(nm*var_in)   ! module "Little_subroutines"

   case ('hertz', 'Hz', 'hz', 'HZ', 'Hertz', 'HERTZ') ! frequency -> [eV]
      ! Convert frequency into [Hz]:
      f = 1.0d0
      ! Convertion coefficient from [Hz] to [eV]:
      var_out = convert_frequency_to_hw(f*var_in)   ! module "Little_subroutines"

   case ('GHz', 'Ghz', 'GHZ', 'ghz') ! frequency -> [eV]
      ! Convert frequency into [Hz]:
      f = 1.0d9
      ! Convertion coefficient from [Hz] to [eV]:
      var_out = convert_frequency_to_hw(f*var_in)   ! module "Little_subroutines"

   case ('THz', 'Thz', 'THZ', 'thz') ! frequency -> [eV]
      ! Convert frequency into [Hz]:
      f = 1.0d12
      ! Convertion coefficient from [Hz] to [eV]:
      var_out = convert_frequency_to_hw(f*var_in)   ! module "Little_subroutines"

   case ('PHz', 'Phz', 'PHZ', 'phz') ! frequency -> [eV]
      ! Convert frequency into [Hz]:
      f = 1.0d15
      ! Convertion coefficient from [Hz] to [eV]:
      var_out = convert_frequency_to_hw(f*var_in)   ! module "Little_subroutines"

   case ('%', 'percent') ! wavelength -> [eV]
      ! Take care to find the absolute values externally:
      coef_out = -1.0d-2
      var_out = coef_out*var_in

   case default   ! no renormalization
      var_out = var_in

      if (print_message) then
         write(6, '(a)') 'Could not interprete units of photon pulse "'//trim(adjustl(ch_temp))//'"'
         write(6, '(a)') 'Assuming defult: photon energy in [eV]'
      endif
   end select
end subroutine interprete_photon_units




subroutine check_if_read_well(Reason, count_lines, File_name, Err, add_error_info)
   integer, intent(in) :: Reason ! IO flag from intrinsic fortran: read(FN,*,IOSTAT=Reason)
   integer, intent(inout) :: count_lines  ! conter in which line the error appeared
   character(*), intent(in) :: File_name  ! which file we are reading now
   type(Error_handling), intent(inout) :: Err   ! error description construct
   character(*), intent(in), optional :: add_error_info  ! additional info about the error to print
   !-------------------
   logical :: read_well
   character(200) :: Error_descript

   ! Check the Reason, if it read well:
   call read_file(Reason, count_lines, read_well)  ! module "Dealing_with_files"

   if (.not. read_well) then
      write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
      call Save_error_details(Err, 3, Error_descript)
      ! Print main error message on the screen
      print*, trim(adjustl(Error_descript))
      ! Print additional error info:
      if (present(add_error_info)) then
         print*, trim(adjustl(add_error_info))
      endif
   endif
end subroutine check_if_read_well


subroutine check_pulse_specifications(text, laser, read_var)
   character(*), intent(in) :: text ! specification text
   type(Pulse), intent(inout) :: laser ! Laser pulse parameters
   real(8), intent(in) :: read_var  ! value provided (incoming fluence or absorbed dose)
   !------------------------
   integer :: i, Nsiz

   Nsiz = LEN(text)
   if (Nsiz > 0) then ! there was an additional text with specifications:
      do i = 1, Nsiz
         select case (text(i:i))
         case ('d', 'D')
            laser%F = read_var   ! absorbed dose [eV/atom]
            laser%F_in = -1.0d0  ! default, undefined
         case ('f', 'F')
            laser%F_in = read_var   ! incoming fluence [J/cm^2]
            laser%F = -1.0d0  ! default, undefined
         end select
      enddo
   endif
   !print*, 'check_pulse_specifications', Nsiz, text, laser%F_in, laser%F
end subroutine check_pulse_specifications



subroutine check_coordinates_filename(Cell_filename, verbose)  ! check if the name has correct extension
   character(*), intent(inout) :: Cell_filename
   logical, intent(in) :: verbose
   !------------------
   character(200) :: filename_extension

   call get_file_extension(trim(adjustl(Cell_filename)), filename_extension)  ! module "Dealing_with_files"

   if (LEN(trim(adjustl(filename_extension))) <= 0) then ! not a valid filename
      Cell_filename = ''
      if (verbose) write(*,'(a)') 'No valid filename with coordinates provided, using default'
   endif
end subroutine check_coordinates_filename



subroutine interpret_user_data_INPUT(FN, File_name, count_lines, string_in, Scell, numpar, Err)
   integer, intent(in) :: FN  ! file to read from
   character(*), intent(in) :: File_name
   integer, intent(inout) :: count_lines  ! line we are reading
   character(*), intent(in) :: string_in ! what was read in the previous line
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !--------------------------
   real(8), dimension(3) :: read_var   ! fluence: starting, ending, step; all in [eV/atom]
   integer :: Reason, num_phon, i, N
   real(8) :: i_min, i_max
   logical :: read_well
   character(200) :: Error_descript, temp_ch, string
   character(20) :: temp_ch1, temp_ch2, temp_ch3


   !print*, trim(adjustl(string_in))

   read_well = .true.   ! to start with
   read_var = 0.0d0     ! unused variable in this case

   ! Read first argument in the line:
   read(string_in,*,IOSTAT=Reason) string

   select case (trim(adjustl(string)))
   !----------------------------------
   case ('output_name', 'Output_name', 'Outout_Name', 'OUTPUT_NAME', 'outname', 'Outname', 'OUTNAME')
      ! User-defined additional text in the output folder name INSTEAD of the default name:
      ! Read second argument in the line:
      read(string_in,*,IOSTAT=Reason) string, numpar%output_extra_name
      if (Reason /= 0) then ! did not read well, use default:
         numpar%output_extra_name = ''
         write(*,'(a)') 'No valid output directory name provided, using default'
      else ! read it well
         ! Make sure the slash is correct, and no forbidden characters are present:
         call ensure_correct_path_separator(numpar%output_extra_name, numpar%path_sep, no_slash=.true.)  ! module "Dealing_with_files"
         !write(*,'(a)') 'Output directory name provided: '//trim(adjustl(numpar%output_extra_name))
      endif

   case ('output_add', 'Output_add', 'Outout_Add', 'OUTPUT_ADD', 'outadd', 'Outdate', 'OUTADD')
      ! User-defined additional text in the output folder AFTER the default name:
      ! Read second argument in the line:
      read(string_in,*,IOSTAT=Reason) string, numpar%output_name_add
      if (Reason /= 0) then ! did not read well, use default:
         numpar%output_name_add = ''
         write(*,'(a)') 'No valid add-text for output name provided, default unchanged'
      else ! read it well
         ! Make sure the slash is correct, and no forbidden characters are present:
         call ensure_correct_path_separator(numpar%output_name_add, numpar%path_sep, no_slash=.true.)  ! module "Dealing_with_files"
      endif

   !----------------------------------
   case ('EADLname', 'EADL_name', 'EADL_Name', 'EADL_NAME')
      ! User-defined name of the file with EADL database:
      read(string_in,*,IOSTAT=Reason) string, temp_ch
      if (Reason /= 0) then ! did not read well, use default:
         write(*,'(a)') 'No valid EADL-database filename provided, using default: ', trim(adjustl(numpar%EADL_file))
      else ! read it well
         ! Make sure the slash is correct:
         call ensure_correct_path_separator(temp_ch, numpar%path_sep)  ! module "Dealing_with_files"
         numpar%EADL_file = trim(adjustl(temp_ch))
      endif

   case ('EPDLname', 'EPDL_name', 'EPDL_Name', 'EPDL_NAME')
      ! User-defined name of the file with EADL database:
      read(string_in,*,IOSTAT=Reason) string, temp_ch
      if (Reason /= 0) then ! did not read well, use default:
         write(*,'(a)') 'No valid EPDL-database filename provided, using default: ', trim(adjustl(numpar%EPDL_file))
      else ! read it well
         ! Make sure the slash is correct:
         call ensure_correct_path_separator(temp_ch, numpar%path_sep)  ! module "Dealing_with_files"
         numpar%EPDL_file = trim(adjustl(temp_ch))
      endif

   !----------------------------------
   case ('Set_V2', 'set_V2', 'set_v2')
      ! Choice of initial atomic velocity distribution:
      numpar%ind_starting_V = 2  ! Maxwell
   case ('Set_V1', 'set_V1', 'set_v1')
      ! Choice of initial atomic velocity distribution:
      numpar%ind_starting_V = 1  ! linear
   case ('Set_V0', 'set_V0', 'set_v0')
      ! Choice of initial atomic velocity distribution:
      numpar%ind_starting_V = 0  ! equal

   case ('Set_power_b', 'SET_POWER_B', 'Set_Power_b', 'Set_Power_B')
      backspace(FN)
      read(FN,*,IOSTAT=Reason) temp_ch1, N
      if (Reason ==0) then ! read well, use the number:
         numpar%power_b = N
      endif
      !print*, numpar%power_b

   !----------------------------------
   case ('print_Ta', 'Print_Ta', 'PRINT_TA', 'PRINT_Ta')
      ! Printout various definitions of atomic temperature:
      numpar%print_Ta = .true.

   !----------------------------------
   case ('print_MFP', 'Print_MFP', 'PRINT_MFP')
      ! Printout mean free paths:
      numpar%print_MFP = .true.

   !----------------------------------
   case ('redo_MFP', 'MFP', 'REDO_MFP')
      ! recalcaulte mean free paths:
      numpar%redo_MFP = .true.

   !----------------------------------
   case ('el-ph', 'EL-PH', 'El-Ph', 'Coupling', 'COUPLING', 'coupling')
      read(FN,*,IOSTAT=Reason) num_phon   ! number of simulations for average electron-phonon coupling parameter
      call read_file(Reason, count_lines, read_well)
      if (read_well) then ! pass the number of simulations
         call prepare_multiple_inputs(numpar, File_name, read_var, .true., num_phon, string=trim(adjustl(string)))   ! below
      else ! use default nuber of iterations
         call prepare_multiple_inputs(numpar, File_name, read_var, .true., string=trim(adjustl(string))) ! below
         backspace(FN)
      endif

   !----------------------------------
   case ('WATER', 'EMBED_WATER', 'EMBED_IN_WATER', 'Water', 'water', 'embed_water', 'Embed_Water', 'Embed_in_water')
      numpar%embed_water = .true.   ! save the flag for water embedding
      read(FN,*,IOSTAT=Reason) numpar%N_water_mol   ! number of water molecules to use
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         numpar%N_water_mol = 100 ! use default
         print*, 'Using default number of water molecules to embed (100)'
         backspace(FN)
         return
      elseif (numpar%N_water_mol < 0) then   ! don't use water
         numpar%embed_water = .false.
         numpar%N_water_mol = 0
         print*, 'Cannot embed with negative numbner of water molecules!'
      endif

   !----------------------------------
   case ('DOS', 'Dos', 'dos', 'do_DOS', 'get_DOS')
      numpar%do_DOS = .true.  ! calculate DOS

   !----------------------------------
   case ('print_CDF', 'save_CDF', 'get_CDF', 'Print_CDF', 'Save_CDF', 'Get_CDF')
      numpar%save_CDF = .true.   ! printout CDF file with oscillators

   !----------------------------------
   case ('KAPPA', 'Kappa', 'kappa', 'conductivity', 'do_kappa', 'Do_kappa', 'Get_kappa', 'get_kappa', &
         'KAPPA_DYN', 'Kappa_dyn', 'kappa_dyn', 'Kappa_Dyn', 'Kappa_dynamical', 'kappa_dynamical')

      select case (trim(adjustl(string)))
      case ('KAPPA', 'Kappa', 'kappa', 'conductivity', 'do_kappa', 'Do_kappa', 'Get_kappa', 'get_kappa')
         print*, 'Electronic heat conductivity will be calculated (static)'
         numpar%do_kappa = .true.   ! statically calculate K (electron heat conductivity vs Te)
      case ('KAPPA_DYN', 'Kappa_dyn', 'kappa_dyn', 'Kappa_Dyn', 'Kappa_dynamical', 'kappa_dynamical')
         print*, 'Electronic heat conductivity will be calculated (dynamic)'
         numpar%do_kappa_dyn = .true.   ! dynamically calculate K (for transient Te)
      end select

      read(FN,'(a)',IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well)
      read(temp_ch,*,IOSTAT=Reason) numpar%kappa_Te_min, numpar%kappa_Te_max, numpar%kappa_dTe
      if (Reason /= 0) then   ! use default values
         numpar%kappa_Te_min = 300.0d0
         numpar%kappa_Te_max = 30000.0d0
         numpar%kappa_dTe = 100.0d0
         numpar%kappa_model = 0
         if (numpar%do_kappa) then
            print*, 'With default parameters of the model (Kubo-Greenwood)'
         elseif (numpar%do_kappa_dyn) then
            print*, 'With default parameters of the model (dynamical coupling)'
         endif
         backspace(FN)
         return
      else  ! check if the model index provided
         read(FN,*,IOSTAT=Reason) numpar%kappa_model
         call read_file(Reason, count_lines, read_well)
         if (Reason /= 0) then   ! model index not provided
            if (numpar%do_kappa) then
               print*, 'With default model (numerical Onsager coefficients)'
            elseif (numpar%do_kappa_dyn) then
               print*, 'With default model (dynamical coupling)'
            endif
            backspace(FN)
            return
         endif
      endif

   !----------------------------------
   case ('PROBE', 'Probe', 'probe')
      ! Calculate optical parameters (and electronic heat conductivity, ir requested), and with which model:
      read(FN,*,IOSTAT=Reason) numpar%optic_model, N, read_var
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         return
      endif

      SCL:do i = 1, size(Scell) ! for all supercells
         if (numpar%optic_model /= 0) then ! yes, calculate optical coefficients:
            numpar%do_drude = .true.   ! included
            Scell(i)%eps%KK = .false.  ! no K-K relations
            select case(N) ! which model to use
            case (2)    ! use Kramers Kronig relations for spectrum
               Scell(i)%eps%KK = .true.
               Scell(i)%eps%all_w = .true.
            case (1)    ! calculate spectrum, but directly, without using Kramers Kronig relations
               Scell(i)%eps%all_w = .true.
            case default
               Scell(i)%eps%all_w = .false.
            end select
         else
            numpar%do_drude = .false.  ! not included
         endif

         Scell(i)%eps%E_min = read_var(1) ! starting point of the grid of energy [eV]
         Scell(i)%eps%E_max = read_var(2) ! ending point of the grid of energy [eV]
         Scell(i)%eps%dE = read_var(3)    ! grid step of energy [eV]

         ! Absorbtion of how many rays (0=exclude, 1=1st ray, (>1)=sum all); probe-pulse wavelength [nm]; probe duration FWHM [fs]
         read(FN,*,IOSTAT=Reason) numpar%drude_ray, Scell(i)%eps%l, Scell(i)%eps%tau
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            return
         endif
         !if (.not.numpar%do_drude) Scell(i)%eps%tau = -0.0d0 ! to exclude convolution if there is no probe pulse
         Scell(i)%eps%ReEps0 = 0.0d0	! to start with
         Scell(i)%eps%ImEps0 = 0.0d0	! to start with
         Scell(i)%eps%w = 2.0d0*g_Pi*g_cvel/(Scell(i)%eps%l*1d-9) ! [1/sec] frequency

         ! Angle of prob-pulse with respect to normal [degrees]; material thickness [nm]:
         read(FN,*,IOSTAT=Reason) Scell(i)%eps%teta, Scell(i)%eps%dd
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i5,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            return
         endif
         Scell(i)%eps%teta = Scell(i)%eps%teta*g_Pi/(180.0d0) !c [radians]
      enddo SCL

   !----------------------------------
   case ('size', 'Size', 'SIZE')
      print*, 'Supercell size variation will be performed to plot potential energy curve'
      numpar%change_size = .true. ! do changing size
      ! Optional number of points for vary_size:
      read(FN,*,IOSTAT=Reason) i_min, i_max, N
      call read_file(Reason, count_lines, read_well)
      if (read_well) then ! save and use it later
         numpar%change_size_min = i_min
         numpar%change_size_max = i_max
         numpar%change_size_step = N
         write(temp_ch1, '(f12.2)') numpar%change_size_min
         write(temp_ch2, '(f12.2)') numpar%change_size_max
         write(temp_ch3,'(i8)') numpar%change_size_step
         write(temp_ch, '(a)') trim(adjustl(temp_ch1))//' : '//trim(adjustl(temp_ch2))//' : '//trim(adjustl(temp_ch3))
         print*, 'With parameters of the supercell (min, max, grid): '//trim(adjustl(temp_ch))
      else ! reread the line next time, if it was not an integer number
         BACKSPACE(FN)
         print*, 'With default parameters of the supercell min, max, and step'
      endif
      !write(*,'(a)') trim(adjustl(m_starline))

   !----------------------------------
   case default
      ! Check if the user needs any additional info (by setting the flags):
      call interprete_additional_data(string, numpar%path_sep, change_size=numpar%change_size, contin=Err%Stopsignal, &
                  allow_rotate=numpar%allow_rotate, verbose=numpar%verbose, nonverbose=numpar%nonverbose) ! module "Read_input_data"

   endselect
end subroutine interpret_user_data_INPUT



subroutine prepare_multiple_inputs(numpar, File_name, read_var, do_eph, num_phon, string)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   character(*), intent(in) :: File_name  ! given filename
   real(8), dimension(3), intent(in) :: read_var   ! fluence: starting, ending, step; all in [eV/atom]
   logical, intent(in), optional :: do_eph  ! for electron-phonon coupling, it requires special multiple input
   integer, intent(in), optional :: num_phon ! how many simulation for average electron-phonon coupling
   character(*), intent(in), optional :: string
   !-----------------------------------
   real(8) :: dose_cur, RN, x_cur, x_0, pulse
   integer :: N, i, FN, FN2, file_len, count_lines, j, N_lines, Reason, sz
   character(250) :: Save_file, Cur_file, Num_par_file
   character(250), allocatable, dimension(:) :: File_content
   character(5) :: chtest2
   character(25) :: chtest, ch_temp
   logical :: read_well, do_el_phon, file_exists

   ! Check, if user requested average electron-phonon coupling calculations:
   if (present(do_eph)) then
      do_el_phon = do_eph
   else  ! not average electron-phonon calculations
      do_el_phon = .false.
   endif

   ! Check, if it is fluence grid, or electron-phonon grid:
   if (do_el_phon) then
      if (present(num_phon)) then   ! user provided number of points
         N = num_phon-1
         if (N < 1) N = 0
      else  ! use default value:
         N = 10
      endif
   else ! then must be the fluence points:
      ! Get the total number of fluence points:
      N = ceiling(abs(read_var(2) - read_var(1))/abs(read_var(3)))
   endif

   ! Set the name of the numerical_parameters file:
   Num_par_file = trim(adjustl(numpar%input_path))//trim(adjustl(m_NUMERICAL_PARAMETERS))//'.txt'

   ! Set temporary file name:
   file_len = LEN(trim(adjustl(File_name)))  ! length of the string
   Save_file = trim(adjustl(File_name(1:file_len-4)))//'_SAVE.txt'

   ! Create the temporary save file
   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      call copy_file(File_name, Save_file, 1, add_com='* /YQ') ! module "Dealing_with_files"
   else
      call copy_file(File_name, Save_file) ! module "Dealing_with_files"
   endif

   ! Open the file:
   FN=4000
   open(UNIT=FN, FILE = trim(adjustl(Save_file)), status = 'old', action='read')
   call Count_lines_in_file(FN, N_lines)  ! module "Dealing_with_files"
   ! Knowing the size, allocate the array:
   allocate(File_content(N_lines))
   ! Read the saved file:
   count_lines = 0   ! to start with
   do j = 1, N_lines
      read(FN, '(a)', IOSTAT=Reason) File_content(j)   ! read the current line
      call read_file(Reason, count_lines, read_well)   ! modlue "Dealing_with_files"
      if ( (.not.read_well) .and. (numpar%path_sep == '/') ) then ! if it is Linux
         backspace(FN)  ! to reread the line
         count_lines = count_lines - 1 ! reread the same line, don't count it as the next one
         read(FN, '(a)', IOSTAT=Reason) File_content(j)(1:sz) ! read it again, now knowing the size
         call read_file(Reason, count_lines, read_well) ! modlue "Dealing_with_files"
      endif
      if (.not.read_well) then
         print*, 'Problem in prepare_multiple_inputs: cannot read line ', count_lines, ' in file '//trim(adjustl(Cur_file))
         goto 3440
      endif
   enddo
   ! Close files:
   call close_file('delete', FN=FN) ! module "Dealing_with_files"

   ! Make a new copy of the input file with the changed fluence only:
   do i = 1, N
      ! Set the new dose:
      !dose_cur = min(read_var(1),read_var(2)) + dble(i)*abs(read_var(3))
      dose_cur = read_var(1) + dble(i) * sign( read_var(3), (read_var(2)-read_var(1)) )

      ! Set the new file name:
      FN2=FN+i
      write(chtest2,'(i5)') i
      Cur_file = trim(adjustl(File_name(1:file_len-4)))//'_'//trim(adjustl(chtest2))//'.txt'
      open(UNIT=FN2, FILE = trim(adjustl(Cur_file)))

      ! Copy data into this new file, accept for the fluence line:
      do j = 1, N_lines

         if (do_el_phon) then ! it is electron-phonon calculations, requires changing a few parameters:
            select case (j)
            case default
               if (present(string)) then
                  if (trim(adjustl(File_content(j))) /= string) then ! don't make extra copies of the file
                     write(FN2, '(a)') trim(adjustl(File_content(j)))   ! just copy this line
                  endif
               else
                  write(FN2, '(a)') trim(adjustl(File_content(j)))   ! just copy this line
               endif
            case (5) ! start of simulation [fs]
               ! Get the user-provided simulation time:
               read( File_content(j), *, IOSTAT=Reason) x_0
               ! Sample randomly simulation starting time:
               call random_number(RN)  ! [0:1]
               x_cur = x_0 * (1.0d0 + 0.2d0*(RN-0.5d0))
               write(chtest, '(f16.5)') x_cur
               write(FN2, '(a,a)') trim(adjustl(chtest)), '          ! start of simulation [fs]'
            case (6) ! end of simulation [fs]
               ! Get the user-provided simulation time:
               read( File_content(j), *, IOSTAT=Reason) x_0
               ! Sample randomly simulation starting time:
               call random_number(RN)  ! [0:1]
               x_cur = x_0 * (1.0d0 + 0.5d0*(RN-0.5d0))
               pulse = x_cur  ! save it for later, to set pulse duration equal to this value
               write(chtest, '(f16.5)') x_cur
               write(FN2, '(a,a)') trim(adjustl(chtest)), '          ! end of simulation time [fs]'
            case (8) ! absorbed dose per this pulse [eV/atom] (min, max, step)
               ! Get the user-provided simulation time:
               ch_temp = '' ! to start with
               if ( .not. it_is_number(File_content(j)) ) then ! there is pulse specificatin parameter
                  read( File_content(j), *, IOSTAT=Reason) ch_temp, x_0
               else
                  read( File_content(j), *, IOSTAT=Reason) x_0
               endif
               ! Sample randomly simulation starting time:
               call random_number(RN)  ! [0:1]
               x_cur = x_0 * (1.0d0 + 0.2d0*(RN-0.5d0))
               write(chtest, '(f16.5)') x_cur
               chtest = trim(adjustl(ch_temp))//'  '//trim(adjustl(chtest))   ! amend the specification parameter, if any
               write(FN2, '(a,a)') trim(adjustl(chtest)), '          ! absorbed dose per this pulse [eV/atom]'
            case (10)   ! pulse FWHM-duration [fs]
               write(chtest, '(f16.5)') pulse
               write(FN2, '(a,a)') trim(adjustl(chtest)), '          ! pulse FWHM-duration [fs]'
            end select

         else ! it is fluence grid:
            if (j == 8) then  ! dose is set in the line #8
               ch_temp = '' ! to start with
               if ( .not. it_is_number(File_content(j)) ) then ! there is pulse specificatin parameter
                  read( File_content(j), *, IOSTAT=Reason) ch_temp
               endif
               write(chtest, '(f16.5)') dose_cur
               chtest = trim(adjustl(ch_temp))//'  '//trim(adjustl(chtest))   ! amend the specification parameter, if any
               write(FN2, '(a,a)') trim(adjustl(chtest)), '          ! absorbed dose [eV/atom]'
            else  ! not the dose
               if (present(string)) then
                  if (trim(adjustl(File_content(j))) /= string) then ! don't make extra copies of the file
                     write(FN2, '(a)') trim(adjustl(File_content(j)))   ! just copy this line
                  endif
               else
                  write(FN2, '(a)') trim(adjustl(File_content(j)))   ! just copy this line
               endif
            endif
         endif ! (do_el_phon)
      enddo ! j = 1, N_lines

      call close_file('close', FN=FN2) ! module "Dealing_with_files"

      ! Also, make a copy of the numerical_parameters file:
!       Cur_file = trim(adjustl(numpar%input_path))//trim(adjustl(m_NUMERICAL_PARAMETERS))
!       write(Cur_file,'(a,a,a,a)') trim(adjustl(Cur_file)), '_', trim(adjustl(chtest2)), '.txt'
!       inquire(file=trim(adjustl(Cur_file)),exist=file_exist) ! check if input file is there
!       if (file_exist) then ! Create the second file:
!          if (numpar%path_sep .EQ. '\') then	! if it is Windows
!             call copy_file(Num_par_file, Cur_file, 1, add_com='* /YQ') ! module "Dealing_with_files"
!          else
!             call copy_file(Num_par_file, Cur_file) ! module "Dealing_with_files"
!          endif
!       endif
   enddo

3440 continue
end subroutine prepare_multiple_inputs



subroutine multiply_input_files(Folder_name, File_name_in, verbose)
   character(*), intent(in) :: Folder_name, File_name_in  ! input directory and file
   logical, intent(in) :: verbose
   !-----------------------
   character(300), dimension(:), allocatable :: File_content
   character(200) :: File_name, Copy_file_name
   character(300) :: read_line, replace_line(100,50)
   character(10) :: temp_ch
   integer :: FN, FN1, FN2, N, i, j, k, Nsiz, count_lines, Reason, sz, line_num(100,50)
   integer :: N_lines, i_block, i_line
   logical :: file_exists, file_opened, read_well, was_closed

   ! Check if the copy-instructions file is present at all:
   Copy_file_name = trim(adjustl(Folder_name))//trim(adjustl(m_COPY_INPUT))
   inquire(file=trim(adjustl(Copy_file_name)),exist=file_exists)
   if (.not. file_exists) then
      return ! nothing else to do here
   else
      if (verbose) then
         write(*,'(a)') ' Multiple input files will be created automatically from '//trim(adjustl(m_COPY_INPUT))
      endif
   endif

   !-----------------------
   ! 1) Read input file:
   inquire(file=trim(adjustl(File_name_in)), opened=file_opened, number=FN)
   if (.not.file_opened) then
      was_closed = .true.
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name_in)), action = 'read')
   else
      was_closed = .false.
      rewind(FN)
   endif
   ! Get how many lines are in the file
   call Count_lines_in_file(FN, N_lines)  ! module "Dealing_with_files"
   allocate(File_content(N_lines))
   ! Read the file:
   count_lines = 0   ! to start with
   do j = 1, N_lines
      read(FN, '(a)', IOSTAT=Reason) File_content(j)   ! read the current line
      call read_file(Reason, count_lines, read_well)   ! modlue "Dealing_with_files"
      if ( (.not.read_well) .and. (Folder_name(LEN(Folder_name):LEN(Folder_name)) == '/') ) then ! if it is Linux
         backspace(FN)  ! to reread the line
         count_lines = count_lines - 1 ! reread the same line, don't count it as the next one
         read(FN, '(a)', IOSTAT=Reason) File_content(j)(1:sz) ! read it again, not knowing the size
         call read_file(Reason, count_lines, read_well) ! modlue "Dealing_with_files"
      endif
      if (.not.read_well) then
         print*, 'Problem in multiply_input_files: cannot read line ', count_lines, ' in file '//trim(adjustl(File_name_in))
         return
      endif
   enddo
   if (was_closed) call close_file('close', FN=FN) ! module "Dealing_with_files"

   !-----------------------
   ! 2) Read copy-file data:
   open(NEWUNIT=FN1, FILE = trim(adjustl(Copy_file_name)), action = 'read')
   ! Get how many lines are in the file:
   call Count_lines_in_file(FN1, N_lines)  ! module "Dealing_with_files"
   ! Count, how many copies are required:
   replace_line = '' ! to start with
   line_num = 0   ! to start with
   i_block = 0 ! to start with
   i_line = 0 ! to start with
   Nsiz = 0 ! to start with
   count_lines = 0 ! to start with
   RDCL:do i = 1, N_lines
      read(FN1, '(a)', IOSTAT=Reason) read_line   ! read the current line
      call read_file(Reason, count_lines, read_well)   ! modlue "Dealing_with_files"
      if (Reason .LT. 0) then ! end of file reached
         exit RDCL
      elseif (.not.read_well) then
         print*, 'Problem in multiply_input_files: cannot read line ', count_lines, ' in file '//trim(adjustl(m_COPY_INPUT))
         return
      endif

      if (trim(adjustl(read_line(1:1))) /= '!') then ! it is not a commen line, try to interprete it
         select case (trim(adjustl(read_line)))
         case ('NEW', 'New', 'new', 'COPY', 'Copy', 'copy') ! count as new copy
            i_block = i_block + 1   ! count blocks
            i_line = 0  ! restart line counter
         case('', '!') ! skipline
         case default
            if (i_block > 0) then
               i_line = i_line + 1
               read(read_line, *, IOSTAT=Reason) line_num(i_block,i_line)
               call read_file(Reason, count_lines, read_well)   ! modlue "Dealing_with_files"
               if (read_well) then
                  if (line_num(i_block,i_line) < 10) then
                     replace_line(i_block,i_line) = trim(adjustl(read_line(3:)))
                  elseif (line_num(i_block,i_line) < 100) then
                     replace_line(i_block,i_line) = trim(adjustl(read_line(4:)))
                  endif
               else ! nullify the wrong reading
                  line_num(i_block,i_line) = 0
               endif ! (read_well)
               !print*, i_block, i_line, line_num(i_block,i_line), trim(adjustl(replace_line(i_block,i_line)))
            endif ! (i_block > 0)
         end select
      endif ! (trim(adjustl(read_line(1:1))) /= '!')
   enddo RDCL ! i = 1, N_lines
   call close_file('close', FN=FN1) ! module "Dealing_with_files"
   Nsiz = i_block

   !-----------------------
   ! 3) Prepare copies of the input file (with requested modifications):
   do i = 1, Nsiz
      write(temp_ch,'(i0)') i
      File_name = File_name_in(1:LEN(File_name_in)-4)//'_'//trim(adjustl(temp_ch))//'.txt'
      open(NEWUNIT=FN2, FILE = trim(adjustl(File_name)))

      ! Copy into this file everything from the INPUT, replacing only what's required:
      k = 1 ! to start with
      do j = 1, size(File_content)
         if (j == line_num(i,k)) then ! replace this line
            write(FN2, '(a)') trim(adjustl(replace_line(i,k)))
            k = k + 1   ! counter
         else ! copy the line unchanged
            write(FN2, '(a)') trim(adjustl(File_content(j)))
         endif
      enddo

      call close_file('close', FN=FN2) ! module "Dealing_with_files"
   enddo ! i = 1, Nsiz


end subroutine multiply_input_files



! Reads additional data from the command line passed along with the XTANT:
subroutine get_add_data(path_sep, change_size, contin, allow_rotate, verbose, nonverbose)
   character(1), intent(inout) :: path_sep
   logical, intent(inout) :: change_size
   logical, intent(out) :: contin
   logical, intent(out) :: allow_rotate
   logical, intent(out) :: verbose, nonverbose
   !---------------
   character(1000) :: string
   integer :: i_arg, count_args, N_arg
   logical :: read_well

   ! Default values:
   change_size = .false. ! don't do changing size
   verbose = .false.    ! don't print a lot of stuff
   nonverbose = .false. ! print as normal

   ! Identify the OS by the system-used path separator:
   call Path_separator(path_sep) ! module "Dealing_with_files"

   ! Count how many arguments the user provided:
   N_arg = COMMAND_ARGUMENT_COUNT() ! Fortran intrinsic function

   read_well = .true.   ! to start with
   count_args = 0 ! to start with

   ALLARG:do i_arg = 1, N_arg ! read all the arguments passed
      ! Read the argument provided:
      call GET_COMMAND_ARGUMENT(i_arg,string)  ! intrinsic

      ! Act on the command passed:
      call interprete_additional_data(string, path_sep, change_size, contin, allow_rotate, verbose, nonverbose)  ! below

   enddo ALLARG
end subroutine get_add_data


subroutine interprete_additional_data(string, path_sep, change_size, contin, allow_rotate, verbose, nonverbose)
   character(*), intent(in) :: string
   character(1), intent(inout) :: path_sep
   logical, intent(inout), optional :: change_size
   logical, intent(out), optional :: contin
   logical, intent(out), optional :: allow_rotate
   logical, intent(out), optional :: verbose, nonverbose
   !---------------
   character(1000) :: read_string, printline, ch_temp, string_read
   character(200) :: file_name
   integer :: FN, Reason, count_lines
   logical :: file_opened, read_text_well, read_well, file_exists

   string_read = trim(adjustl(string))
   ! Trim minus signs if present (for legacy reasons):
   do while (string_read(1:1) == '-')
      string_read = string_read(2:)
   enddo

   ! Interpret the command:
   select case (trim(adjustl(string_read)))
   case ('verbose', 'VERBOSE', 'Verbose')
      print*, 'Verbose on: XTANT will print markers for testing and debugging'
      if (present(verbose)) verbose = .true.
      write(*,'(a)') trim(adjustl(m_starline))

   case ('nonverbose', 'NONVERBOSE', 'Nonverbose', 'noverbose', 'NOVERBOSE')
      print*, 'Nonverbose on: XTANT will print almost nothing'
      if (present(nonverbose)) nonverbose = .true.
      write(*,'(a)') trim(adjustl(m_starline))

   case ('allow_rotation', 'allow_rotate', 'no_ang_removal')
      print*, 'The angular momenta of the sample will not be removed'
      if (present(allow_rotate)) allow_rotate = .true. ! don't remove angular momentum from initial conditions
      write(*,'(a)') trim(adjustl(m_starline))

   case ('size', 'Size', 'SIZE')
      print*, 'Supercell size variation will be performed to plot potential energy curve'
      if (present(change_size)) change_size = .true. ! do changing size
      write(*,'(a)') trim(adjustl(m_starline))

   case ('test', 'TEST', 'Test')
      print*, 'Wow, it really works!'
      if (present(contin)) contin = .true.  ! stop calculations, user only wanted some info
      write(*,'(a)') trim(adjustl(m_starline))

   case ('Matter', 'matter', 'MATTER', 'Materials', 'materials', 'list', 'LIST', 'List')
      write(*,'(a)') trim(adjustl(m_starline))
      ! Create file with list of available materials:
      call Get_list_of_materials(path_sep)  ! below

      write(*,'(a)') trim(adjustl(m_starline))
      if (present(contin)) contin = .true.  ! stop calculations, user only wanted some info

   case ('help', 'HELP', 'Help')
      ! Filename with help:
      file_name = trim(adjustl(m_INPUT_directory))//path_sep//trim(adjustl(m_INFO_directory))//path_sep//trim(adjustl(m_HELP_file))

      inquire(file=trim(adjustl(file_name)),exist=file_exists)
      if (.not.file_exists) then ! no file, cannot print help
         write(*,'(a)') 'Could not find file ', trim(adjustl(file_name))
         write(*,'(a)') 'Cannot help, sorry. Read the manual.'
      else ! (.not.file_exists)
         FN=200
         open(UNIT=FN, FILE = trim(adjustl(file_name)), status = 'old', action='READ')
         inquire(file=trim(adjustl(file_name)),opened=file_opened)
         if (.not.file_opened) then
            write(*,'(a)') 'Could not open file ', trim(adjustl(file_name))
            write(*,'(a)') 'Cannot help, sorry. Read the manual.'
         else ! (.not.file_opened)
            read_text_well = .true. ! to start with
            count_lines = 0   ! to start with
            do while (read_text_well)
               read(FN,'(a)',IOSTAT=Reason) printline
               call read_file(Reason, count_lines, read_text_well)   ! module "Dealing_with_files"
               if (Reason > 0) then   ! something wrong in the line
                  write(*,'(a)') 'Problem reading file '//trim(adjustl(file_name))
                  write(ch_temp, '(i)') count_lines
                  write(*,'(a)') 'in line '//trim(adjustl(ch_temp))
                  read_well = .false.
               elseif (Reason < 0) then ! end of file reached ...
                  close(FN)
               else
                  write(*,'(A)') trim(adjustl(printline))
               endif
            enddo
         endif ! (.not.file_opened)
      endif ! (.not.file_exists)

      write(*,'(a)') trim(adjustl(m_starline))
      if (present(contin)) contin = .true.  ! stop calculations, user only wanted some help

   case ('info', 'INFO', 'Info')
      ! Filename with help:
      file_name = trim(adjustl(m_INPUT_directory))//path_sep//trim(adjustl(m_INFO_directory))//path_sep//trim(adjustl(m_INFO_file))

      inquire(file=trim(adjustl(file_name)),exist=file_exists)
      if (.not.file_exists) then ! no file, cannot print help
         write(*,'(a)') 'Could not find file ', trim(adjustl(file_name))
         write(*,'(a)') 'Cannot help, sorry. Read the manual.'
      else ! (.not.file_exists)
         FN=201
         open(UNIT=FN, FILE = trim(adjustl(file_name)), status = 'old', action='READ')
         inquire(file=trim(adjustl(file_name)),opened=file_opened)
         if (.not.file_opened) then
            write(*,'(a)') 'Could not open file ', trim(adjustl(file_name))
            write(*,'(a)') 'Cannot help, sorry. Read the manual.'
         else ! (.not.file_opened)
            read_text_well = .true. ! to start with
            count_lines = 0   ! to start with
            do while (read_text_well)
               read(FN,'(a)',IOSTAT=Reason) printline
               call read_file(Reason, count_lines, read_text_well)   ! module "Dealing_with_files"
               if (Reason > 0) then   ! something wrong in the line
                  write(*,'(a)') 'Problem reading file '//trim(adjustl(file_name))
                  write(ch_temp, '(i)') count_lines
                  write(*,'(a)') 'in line '//trim(adjustl(ch_temp))
                  read_well = .false.
               elseif (Reason < 0) then ! end of file reached ...
                  close(FN)
               else
                  write(*,'(A)') trim(adjustl(printline))
               endif
            enddo
         endif ! (.not.file_opened)
      endif ! (.not.file_exists)

      write(*,'(a)') trim(adjustl(m_starline))
      if (present(contin)) contin = .true.  ! stop calculations, user only wanted some info
   case default
   end select
end subroutine interprete_additional_data


subroutine Get_list_of_materials(path_sep)
   character(*), intent(in) :: path_sep ! file name
   !--------------------
   integer :: FN_temp, FN, count_lines, Reason, i, open_status, iret
   character(200) :: Error_descript, File_name, read_line, command, File_scratch
   logical :: file_opened, file_exist

   File_scratch = trim(adjustl(m_INPUT_directory))//path_sep//'Scratch.txt'
   File_name = trim(adjustl(m_INPUT_directory))//path_sep//trim(adjustl(m_List_ofmaterials))

   ! Make a list of the materials available:
   if (path_sep == '\') then	! if it is Windows
      command = 'dir '//trim(adjustl(m_INPUT_directory))//' /b > '//trim(adjustl(File_scratch))
   else ! linux:
      command = "ls -t "//trim(adjustl(m_INPUT_directory))//" > "//trim(adjustl(File_scratch))
   endif
#ifdef OMP_inside
   iret = system(trim(adjustl(command)))   ! execute the command to save file names in the temp file
#else
   call system(trim(adjustl(command))) ! execute the command to save file names in the temp file
#endif

   FN=200
   open(UNIT=FN, FILE = trim(adjustl(File_name)))
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      write(*,'(a)') 'Could not open file ', trim(adjustl(File_name))
      write(*,'(a)') 'Cannot help, sorry. Check materials manually.'
      goto 9992
   endif

   ! Open the files with directory names:
   open(NEWUNIT=FN_temp, file=trim(adjustl(File_scratch)), iostat=open_status, action='read')
   if ( open_status /= 0 ) then
      print *, 'Could not open ',trim(adjustl(File_scratch)),' for listing.', ' Unit = ', FN_temp
      write(*,'(a)') 'Cannot help, sorry. Check materials manually.'
      goto 9992
   endif

   ! Printout the list of materials:
   write(*,'(a)') ' Materials already available in XTANT:'
   count_lines = 0   ! to start with
   RDLST:do
      read(FN_temp,'(a)',IOSTAT=Reason) read_line
      count_lines = count_lines + 1
      if (Reason < 0) then ! end of file
         exit RDLST
      else  ! line in the fil
         ! Check if the line is commented out:
         if (trim(adjustl(read_line(1:1))) == '!') then ! it is a comment
            ! ignor this folder / file
         else ! don't ignore this name
            ! Check if it is a directory:
            inquire(DIRECTORY=trim(adjustl(m_INPUT_directory))//path_sep//trim(adjustl(read_line)), exist=file_exist)
            if (file_exist) then ! if it is a directory, it means it can be a material:
               select case ( trim(adjustl(read_line)) )  ! check if it is material or just forled
               case (m_Atomic_parameters, m_INFO_directory, m_HELP_file, m_DFTB_directory, m_DFTB_norep_directory, &
                  m_3TB_directory, m_BOP_directory, m_xTB_directory) ! work directory, not a material
                  ! skip this line
               case default   ! material name
                  write(FN,'(a)') trim(adjustl(read_line))  ! save in the file
                  write(*,'(a)') trim(adjustl(read_line))   ! print on the screen
               endselect
            endif ! (file_exist)
         endif ! (trim(adjustl(read_line)) == '!')
      endif ! (Reason < 0)
   enddo RDLST

9992 continue
   call close_file('delete', trim(adjustl(File_scratch)))   ! module "Dealing_with_files"
   call close_file('close', FN=FN)  ! module "Dealing_with_files"
end subroutine Get_list_of_materials


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine check_all_warnings(print_to, laser, Scell, Err)
   integer, intent(in) :: print_to  ! file number to print to
   type(Super_cell), dimension(:), intent(in) :: Scell   ! suoer-cell with all the atoms inside
   type(Pulse), dimension(:), intent(in) :: laser        ! Laser pulse parameters
   type(Error_handling), intent(inout), optional :: Err  ! errors and warnings save
   !---------------------
   character(100) :: text
   integer :: i, Reason

   do i = 1, size(laser)   ! for all pulses
      ! Printout warning if absorbed dose is too high:
      if (laser(i)%F >= 10.0) then
         write(text,'(f16.3)',IOSTAT=Reason) laser(i)%F
         call printout_warning(print_to, 4, text_to_add=trim(adjustl(text)) ) ! below
         if (present(Err)) call Save_error_details(Err, 0, '', empty=.true., Warn=.true.) ! module "Objects"
      endif

      ! Printout warning if photon energy is too high:
      if (laser(i)%hw >= 1.0d5) then
         write(text,'(f16.3)',IOSTAT=Reason) laser(i)%hw
         call printout_warning(print_to, 1, text_to_add=trim(adjustl(text)) ) ! below
         if (present(Err)) call Save_error_details(Err, 0, '', empty=.true., Warn=.true.) ! module "Objects"
      endif
   enddo

   do i = 1, size(Scell)   ! for all supercells
      ! Printout warning if electron temperature is too high:
      if (Scell(i)%TeeV >= 5.0) then
         write(text,'(f16.3)',IOSTAT=Reason) Scell(i)%TeeV
         call printout_warning(print_to, 2, text_to_add=trim(adjustl(text)) ) ! below
         if (present(Err)) call Save_error_details(Err, 0, '', empty=.true., Warn=.true.) ! module "Objects"
      endif

      ! Printout warning if atomic temperature is too high:
      if (Scell(i)%TaeV >= 1.0) then
         write(text,'(f16.3)',IOSTAT=Reason) Scell(i)%TaeV
         call printout_warning(print_to, 3, text_to_add=trim(adjustl(text)) ) ! below
         if (present(Err)) call Save_error_details(Err, 0, '', empty=.true., Warn=.true.) ! module "Objects"
      endif
   enddo

end subroutine check_all_warnings


subroutine printout_warning(print_to, ind, text_to_print, text_to_add) ! standardized format for warning printouts
   integer, intent(in) :: print_to  ! file number to print to
   integer, intent(in) :: ind    ! index of the warning to be printed out
   character(*), intent(in), optional :: text_to_print   ! optional text to printout TOGETHER with the warning
   character(*), intent(in), optional :: text_to_add     ! optional text to printout INSIDE the warning
   !-------------------------
   character(500) :: ch_lng
   character(5) :: ch_sh

   ! Warning title:
   write(print_to, '(a)') ''
   write(print_to, '(a)') trim(adjustl(m_warnline))
   write(ch_sh, '(i5)') ind
   !write(print_to, '(a)') '   WARNING #'//trim(adjustl(ch_sh))//':'
   write(print_to, '(a)') '                     >>> WARNING <<<'
   if (present(text_to_print)) then
      write(print_to, '(a)') trim(adjustl(text_to_print))
   endif

   !---------------------
   select case (ind)
   !---------------------
   case (1) ! too high photon energy
      if (present(text_to_add)) then
         write(ch_lng,'(a)') '('//trim(adjustl(text_to_add))//' eV)'
      else
         write(ch_lng,'(a)') ''
      endif
      write(print_to, '(a)') 'Photon energy is too high '//trim(adjustl(ch_lng))
      write(print_to, '(a)') 'Relativistic effects are not included in XTANT-3'
      write(print_to, '(a)') 'Thus, high-energy electron kinetics is unreliable!'
      write(print_to, '(a)') 'Proceed with caution, or reduce the photon energy below ~100 keV'
   !---------------------
   case (2) ! too high Te
      if (present(text_to_add)) then
         write(ch_lng,'(a)') '('//trim(adjustl(text_to_add))//' eV)'
      else
         write(ch_lng,'(a)') ''
      endif
      write(print_to, '(a)') 'Electronic temperature is too high '//trim(adjustl(ch_lng))
      write(print_to, '(a)') 'Tight binding approximation may not handle it well.'
      write(print_to, '(a)') 'Band structure and interatomic forces may not be reliable!'
      write(print_to, '(a)') 'Proceed with caution, or reduce the electron tempreature below ~5-10 eV'
   !---------------------
   case (3) ! too high Ta
      if (present(text_to_add)) then
         write(ch_lng,'(a)') '('//trim(adjustl(text_to_add))//' eV)'
      else
         write(ch_lng,'(a)') ''
      endif
      write(print_to, '(a)') 'Atomic temperature is too high '//trim(adjustl(ch_lng))
      write(print_to, '(a)') 'Tight binding molecular dynamics may not handle it well.'
      write(print_to, '(a)') 'Atoma coming too close may pose trouble, if short-range repulsion is not taken care of.'
      write(print_to, '(a)') 'Proceed with caution, or reduce the atomic tempreature below ~1 eV'
   !---------------------
   case (4) ! too high dose / fluence
      if (present(text_to_add)) then
         write(ch_lng,'(a)') '('//trim(adjustl(text_to_add))//' eV/atom)'
      else
         write(ch_lng,'(a)') ''
      endif
      write(print_to, '(a)') 'Absorbed dose is too high '//trim(adjustl(ch_lng))
      write(print_to, '(a)') 'Tight binding may not handle it well.'
      write(print_to, '(a)') 'Resulting in too high excitation may pose trouble!'
      write(print_to, '(a)') 'Proceed with caution, or reduce the dose (fluence) below ~10 eV/atom'

   !---------------------
   case (5) ! too low photon energy
      if (present(text_to_add)) then
         write(ch_lng,'(a)') '('//trim(adjustl(text_to_add))//' eV)'
      else
         write(ch_lng,'(a)') ''
      endif
      write(print_to, '(a)') 'Photon energy is too low (<30 eV) '//trim(adjustl(ch_lng))
      write(print_to, '(a)') 'Conversion from incoming fluence to dose may not work well!'
      write(print_to, '(a)') 'Used atmoic photon cross section is unreliable at such energies.'
      write(print_to, '(a)') 'Proceed with caution, or increase photon energy above ~30 eV,'
      write(print_to, '(a)') 'or change the cross-section from atomic to CDF-based.'
   !---------------------
   end select

   ! Warning ending:
   write(print_to, '(a)') trim(adjustl(m_warnline))
   write(print_to, '(a)') ''
end subroutine printout_warning



!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! Obsolete format, to be deprecated
! Alternative format of input file:

subroutine read_input_txt(File_name, Scell, matter, numpar, laser, Err)
   character(*), intent(in) :: File_name
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Error_handling), intent(inout) :: Err	! error save
   !--------------------------------------------------------
   integer :: FN, count_lines, Reason, i
   character(200) :: Error_descript, read_line
   logical :: file_opened

   Error_descript = 'Obsolete option to read '//trim(adjustl(File_name))//' is omitted, please use INPUT.txt format'
   call Save_error_details(Err, 0, Error_descript)
   print*, trim(adjustl(Error_descript))

end subroutine read_input_txt


end MODULE Read_input_data
