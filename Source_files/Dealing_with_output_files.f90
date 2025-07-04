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
! This module contains subroutines to read input files:

MODULE Dealing_with_output_files
! Open_MP related modules from external libraries:
#ifdef _OPENMP
   USE OMP_LIB, only : OMP_GET_MAX_THREADS
#endif
#ifndef __GFORTRAN__
   USE IFLPORT, only : system, chdir
#endif

use Universal_constants
use Objects
use Atomic_tools, only : pair_correlation_function
use Variables, only : g_numpar, g_matter
use Little_subroutines, only : number_of_types_of_orbitals, name_of_orbitals, set_starting_time, order_of_time, convolution, &
                              convert_hw_to_wavelength, convert_wavelength_to_hw, find_order_of_number, print_time
use Dealing_with_files, only : get_file_stat, copy_file, read_file, close_file, Count_lines_in_file
use Dealing_with_EADL, only : define_PQN
use Gnuplotting
use Read_input_data, only : m_INPUT_directory, m_INFO_directory, m_INFO_file, m_HELP_file, m_starline, m_dashline, &
                           m_INPUT_MINIMUM, m_INPUT_MATERIAL, m_NUMERICAL_PARAMETERS, m_INPUT_ALL, m_Communication, &
                           m_QUOTES_file, printout_warning
use Dealing_with_CDF, only : write_CDF_file

#ifdef MPI_USED
use MPI_subroutines, only : MPI_barrier_wrapper, broadcast_variable
#endif

implicit none
PRIVATE

character(30), parameter :: m_XTANT_version = 'XTANT-3 (version 09.06.2025)'
character(30), parameter :: m_Error_log_file = 'OUTPUT_Error_log.txt'

public :: write_output_files, convolve_output, reset_dt, print_title, prepare_output_files, communicate
public :: close_save_files, close_output_files, save_duration, execute_all_gnuplots, write_energies
public :: XTANT_label, m_Error_log_file, printout_CDF_file, print_a_comforting_message, printout_MFP_file, printout_laser_spectrum

 contains


subroutine write_output_files(numpar, time, matter, Scell)
   type(Super_cell), dimension(:), intent(inout):: Scell ! super-cell with all the atoms inside
   real(8), intent(in) :: time ! time instance [fs]
   type(Solid), intent(inout) :: matter ! parameters of the material
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   !------------------------------------------------------
   type(Energies) :: nrg   ! [eV] energies in the super-cell
   real(8) :: Pressure
   real(8), dimension(3,3) :: Stress
   integer NSC


   ! Prepare output:
   if (numpar%save_PCF) call pair_correlation_function(Scell(1)%MDatoms, matter, Scell, 1, numpar%MPI_param) ! module "Atomic_tools"


   !--------------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return
   endif
   !--------------------------------------------------------------------------


   do NSC = 1, size(Scell)
      ! All subroutines for saving output data into files are within this file below:
      call update_save_files(time, Scell(NSC)%MDatoms, matter, numpar, Scell(NSC))
      call write_temperatures_n_displacements(numpar%FN_temperatures, time, Scell(NSC)%Te, Scell(NSC)%Ta,  &
                                                      Scell(NSC)%Ta_sub, Scell(NSC)%MSD, Scell(NSC)%MSDP)
      ! Renormalization to printing units:
      Pressure = Scell(NSC)%Pressure * 1.0d-9
      Stress = Scell(NSC)%Stress * 1.0d-9
      nrg = Scell(NSC)%nrg
      nrg%E_coul_scc = nrg%E_coul_scc/dble(Scell(NSC)%Na)  ! -> per atom

      ! All writing subroutines are in this file below:
      call write_pressure(numpar%FN_pressure, time, Pressure, Stress)   ! pressure tensore
      call write_energies(numpar%FN_energies, time, nrg) ! energies in various subsystems
      call write_numbers(numpar%FN_numbers, time, Scell(NSC))  ! total numbers
      call write_orb_resolved(numpar%FN_orb_resolved, time, Scell(NSC), matter) ! orbital-resolved electronic data
      call write_holes(numpar%FN_deep_holes, time, matter, Scell(NSC))  ! core holes
      if (numpar%save_raw) call write_atomic_relatives(numpar%FN_atoms_S, Scell(NSC)%MDatoms)   ! atomic coords and velocities
      call write_super_cell(numpar%FN_supercell, time, Scell(NSC))   ! supercell parameters
      call write_electron_properties(numpar%FN_electron_properties, time, Scell, NSC, Scell(NSC)%Ei, matter, numpar, &
               numpar%FN_Ce, numpar%FN_kappa, numpar%FN_kappa_dyn, numpar%FN_Se, numpar%FN_Te, numpar%FN_mu) ! TB electron parameters

      call write_atomic_properties(time, Scell, NSC, matter, numpar) ! atomic parameters

      ! Section of atoms according to masks, if any:
      if (allocated(Scell(1)%Displ)) then
         call write_sectional_displacements(numpar%FN_displacements, time, Scell(NSC), matter) ! atomic displaecements
      endif

      if (numpar%save_XYZ) call write_atomic_xyz(numpar%FN_atoms_R, Scell(1)%MDatoms, matter, Scell(1)%supce(:,:), &
               print_mass=numpar%save_XYZ_extra(1), print_charge=numpar%save_XYZ_extra(2), print_Ekin=numpar%save_XYZ_extra(3))   ! below

      if (numpar%save_CIF) call write_atomic_cif(numpar%FN_cif, Scell(1)%supce(:,:), Scell(1)%MDatoms, matter, time) ! CIF format

      if (numpar%save_Ei) then
         if (numpar%scc) then ! Energy levels include SCC term:
            call save_energy_levels(numpar%FN_Ei, time, Scell(1)%Ei_scc_part)
         else  ! non-SCC (uncorrected energy levels):
            call save_energy_levels(numpar%FN_Ei, time, Scell(1)%Ei)
         endif
      endif

      if (numpar%save_DOS) then  ! Material DOS
         select case (numpar%DOS_splitting)
         case (1) ! with partial DOS
            call save_DOS(numpar%FN_DOS, time, Scell(1)%DOS, Scell(1)%partial_DOS)
         case default   ! no partial dos
            call save_DOS(numpar%FN_DOS, time, Scell(1)%DOS)
         end select
      endif

      select case (numpar%DOS_splitting)  ! orbital-resolved data
         case (1) ! with partial DOS
         call write_coulping(numpar%FN_coupling, time, Scell, NSC, numpar) ! electron-ion coupling
      end select

      if (numpar%save_fe) then
         if (numpar%do_partial_thermal) then ! Electron distribution functions on TB energy levels
            call save_distribution(numpar%FN_fe, numpar, Scell(1), time, Scell(1)%Ei, Scell(1)%fe, Scell(1)%fe_eq, &
                                    Scell(1)%fe_eq_VB, Scell(1)%fe_eq_CB)
         else
            call save_distribution(numpar%FN_fe, numpar, Scell(1), time, Scell(1)%Ei, Scell(1)%fe, Scell(1)%fe_eq)
         endif
      endif

      if (numpar%save_fa) then
         call save_atomic_distribution(numpar%FN_fa, numpar, Scell(1), time, Scell(1)%Ea_grid_out, Scell(1)%fa_out, Scell(1)%fa_eq_out)
         call save_atomic_distribution(numpar%FN_fa_pot, numpar, Scell(1), time, Scell(1)%Ea_pot_grid_out, &
                                          Scell(1)%fa_pot_out, Scell(1)%fa_eq_pot_out)
         call save_atomic_distribution(numpar%FN_fa_tot, numpar, Scell(1), time, Scell(1)%Ea_tot_grid_out, Scell(1)%fa_tot_out)
      endif


      if (numpar%save_fe_grid) call electronic_distribution_on_grid(Scell(1), numpar, time)  ! distribution on grid
      if (numpar%save_PCF) call write_PCF(numpar%FN_PCF, Scell(1)%MDatoms, matter, Scell, 1) ! pair correlation function
      if (numpar%do_drude) call write_optical_coefs(numpar%FN_optics, time, Scell(1)%eps)    ! optical coeffs
      if (Scell(1)%eps%all_w) call write_optical_all_hw(numpar%FN_all_w, time, Scell(1)%eps) ! CDF spectrum
      if (numpar%save_NN) call save_nearest_neighbors(numpar%FN_neighbors, Scell, 1, time)   ! atomic nearest neighbors

      if (allocated(numpar%NN_radii)) then
         call save_nearest_neighbors_element(numpar%FN_element_NN, numpar, Scell, 1, time) ! element-specific nearest neighbors
      endif

      if (numpar%save_diff_peaks) then ! selected diffraction peaks and powder spectrum
         call save_diffraction_peaks(numpar%FN_diff_peaks, time, Scell(1))    ! below
         call save_diffraction_powder(numpar%FN_diff_powder, time, Scell(1))    ! below
      endif


      if (numpar%save_testmode) then   ! testmode additional data (center of mass, rotation, total force, etc.)
         call save_testmode_data(numpar%FN_testmode, time, Scell(1))  ! below
      endif

   enddo ! NSC
end subroutine write_output_files



subroutine printout_CDF_file(numpar, matter, Scell)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Super_cell), dimension(:), intent(in):: Scell ! super-cell with all the atoms inside
   !----------------------------
   character(200) :: file_name
   integer :: NSC, FN
   parameter (NSC = 1)  ! one supercell

   !--------------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything wrong here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return
   endif
   !--------------------------------------------------------------------------


   if (numpar%save_CDF) then ! printout CDF file
      file_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//'OUTPUT_Ritchie_CDF_'// &
                  trim(adjustl(matter%Name))//'.cdf'
      FN = 9998
      open(UNIT=FN, FILE = trim(adjustl(file_name)), status = 'new')
      ! Printout CDF-oscillators coefficients:
      call write_CDF_file(FN, trim(adjustl(matter%Name)), trim(adjustl(matter%Chem)), matter%dens, &
               (Scell(NSC)%E_VB_top-Scell(NSC)%E_VB_bottom), 0.0d0, matter%Atoms)   ! module "Dealing_with_CDF"

      call close_file('close', FN=FN)  ! module "Dealing_with_files"
   endif
end subroutine printout_CDF_file



subroutine printout_laser_spectrum(laser, numpar, matter)
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   !type(Super_cell), dimension(:), intent(in):: Scell ! super-cell with all the atoms inside
   !----------------------------
   real(8) :: coef, t0, t_last, x_tics
   integer :: i_pulse, N_pulse, FN, i, N_grid
   character(300) :: file_name, text_var, file_spectrum, gnu_photon_spectrum
   character(18) :: temp, ch_temp
   character(11) :: sh_cmd, call_slash

   !--------------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything wrong here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return
   endif
   !--------------------------------------------------------------------------


   ! Number of pulses:
   N_pulse = size(laser)

   ! Do for each pulse:
   do i_pulse = 1, N_pulse
      ! Sort the number of absorbed photons, if spectrum is used:
      if (allocated(laser(i_pulse)%Spectrum)) then ! photon spectrum given

         if (N_pulse > 1) then
            write(temp,'(i0)') i_pulse ! number of the pulse
            temp = '_pulse_'//trim(adjustl(temp))
         else
            temp = ''   ! empty
         endif

         text_var = 'OUTPUT_'
         write(text_var,'(a)') trim(adjustl(text_var))//'photon_spectrum'//trim(adjustl(temp))

         file_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(text_var))//'.dat'
         file_spectrum = trim(adjustl(text_var))//'.dat'

         FN = 9989
         open(UNIT=FN, FILE = trim(adjustl(file_name)))

         ! Printout Photon spectra:
         ! 1) Create the comment lines:
         write(FN,'(a)') '#Energy   Incoming   Absorbed    Absorbed_sampled'
         write(FN,'(a)') '#[eV]     [arb.units] [1/eV]   [1/eV]'

         ! 2) Save spectra:
         coef = maxval(laser(i_pulse)%Spectrum_abs(:)) / maxval(laser(i_pulse)%Spectrum(2,:))
         N_grid = size(laser(i_pulse)%Spectrum,2)

         ! Renormalize MC-sampled photon spectrumm:
         laser(i_pulse)%Spectrum_MC = laser(i_pulse)%Spectrum_MC/dble(numpar%NMC)

         do i = 1, N_grid
            write(FN,'(f16.4, es24.8, es24.8, es24.8)') laser(i_pulse)%Spectrum(1,i), &
                                    laser(i_pulse)%Spectrum(2,i) * coef, &
                                    laser(i_pulse)%Spectrum_abs(i), laser(i_pulse)%Spectrum_MC(i)/dble(size(matter%Atoms))
         enddo ! k
         call close_file('close', FN=FN)  ! module "Dealing_with_files"


         !=======================================================
         ! Gnuplot the data:
         if (numpar%path_sep .EQ. '\') then	! if it is Windows
            call_slash = 'call '
            sh_cmd = '.cmd'
         else ! It is linux
            call_slash = './'
            sh_cmd = '.sh'
         endif

         ! Grid size:
         t0 = laser(i_pulse)%Spectrum(1,1)
         t_last = laser(i_pulse)%Spectrum(1,N_grid)

         ! Photon spectrum gnuplotting:
         gnu_photon_spectrum = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))// &
                  'OUTPUT_photon_spectrum'//trim(adjustl(temp))//trim(adjustl(sh_cmd))
         open(NEWUNIT=FN, FILE = trim(adjustl(gnu_photon_spectrum)), action="write", status="replace")

         ! Find order of the number, and set number of tics as tenth of it:
         call order_of_time((t_last - t0), ch_temp, x_tics=x_tics)	! module "Little_subroutines"

         call write_gnuplot_script_header_new(FN, numpar%ind_fig_extention, 3.0d0, x_tics,  'MFPs', &
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
               '\" u 1:2 w l lw \"$LW\" dashtype \"_\" title \"Incomming" ,\'

            write(FN, '(a)') '\"'//trim(adjustl(file_spectrum))// &
                        '\" u 1:3 w l lw \"$LW\" title \"Absorbed\" ,\'

            write(FN, '(a)') '\"'//trim(adjustl(file_spectrum))// &
                        '\" u 1:4 w l lt rgb \"#000000\" lw \"$LW\" dashtype \"_.\" title \"MC Sampled\" '
         endif

         call write_gnuplot_script_ending_new(FN, gnu_photon_spectrum, numpar%path_sep) ! module "Gnuplotting"
         close(FN)

      endif ! (allocated(laser(i_pulse)%Spectrum))
   enddo ! i_pulse

end subroutine printout_laser_spectrum



subroutine printout_MFP_file(numpar, matter, Scell)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Super_cell), dimension(:), intent(in):: Scell ! super-cell with all the atoms inside
   !----------------------------
   real(8) :: t0, t_last, x_tics
   character(300) :: file_name, text_var, file_electron_IMFP, file_electron_EMFP, file_photon_MFP, gnu_electron_MFP, gnu_photon_MFP
   character(11) :: chtemp11, sh_cmd, call_slash
   character(8) :: temp, time_order, col, col_VB
   integer :: NSC, FN, Nsiz, i, Nshl, j, k, N_grid, j_start, count_col
   parameter (NSC = 1)  ! one supercell


   !--------------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything wrong here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return
   endif
   !--------------------------------------------------------------------------


   if (numpar%print_MFP) then ! printout MFP file
      !-----------------------
      ! Inelastic electron MFP:
      text_var = 'OUTPUT_'
      write(text_var,'(a)') trim(adjustl(text_var))//trim(adjustl(matter%Name))//'_Electron_IMFP'

      file_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(text_var))//'.dat'
      file_electron_IMFP = trim(adjustl(text_var))//'.dat'

      FN = 9998
      open(UNIT=FN, FILE = trim(adjustl(file_name)), status = 'new')

      ! Printout MFPs:
      ! 1) Create the comment lines:
      write(FN,'(a)',advance='no') '#Energy   '
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         do j = 1, Nshl ! for all shells of this atom
            if ((i == 1) .and. (j == Nshl)) then ! valence
               ! skip and print at the end
               !write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))
            else ! core shell
               write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(i)%Shell_name(j)))
            endif
         enddo ! j
      enddo ! i
      Nshl = size(matter%Atoms(1)%Ip)
      write(FN,'(a)') '   '//trim(adjustl(matter%Atoms(1)%Shell_name(Nshl)))//'  Total'

      ! 2) Save MFPs for all shells and total:
      N_grid = size(matter%Atoms(1)%El_MFP(1)%E)
      do k = 1, N_grid
         write(FN,'(f16.4)',advance='no') matter%Atoms(1)%El_MFP(1)%E(k)
         ATOMS:do i = 1, size(matter%Atoms) ! for all atoms
            Nshl = size(matter%Atoms(i)%Ip)
            SHELLS:do j = 1, Nshl ! for all shells of this atom
               if ((i == 1) .and. (j == Nshl)) then ! valence
                  ! skip and print at the end
                  !write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))
               else ! core shell
                  write(FN,'(es24.8)',advance='no') 1.0d0/matter%Atoms(i)%El_MFP(j)%L(k)  ! MFP for ionizing this shell
               endif
            enddo SHELLS
         enddo ATOMS
         Nshl = size(matter%Atoms(1)%Ip)
         write(FN,'(es24.8, es24.8)') 1.0d0/matter%Atoms(1)%El_MFP(Nshl)%L(k), 1.0d0/matter%El_MFP_tot%L(k)  ! total
      enddo ! k
      call close_file('close', FN=FN)  ! module "Dealing_with_files"


      !-----------------------
      ! Elastic electron MFP:
      write(text_var,'(a,a,a)') 'OUTPUT_', trim(adjustl(matter%Name))//'_Electron_EMFP'
      file_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(text_var))//'.dat'
      file_electron_EMFP = trim(adjustl(text_var))//'.dat'

      FN = 9695
      open(UNIT=FN, FILE = trim(adjustl(file_name)))
      ! Printout MFPs:
      ! 1) Create the comment lines:
      write(FN,'(a)',advance='no') '#Energy   '
      do i = 1, size(matter%Atoms) ! for all atoms
         write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Name))
      enddo ! i
      write(FN,'(a)') '  Total'

      ! 2) Save MFPs for all shells and total:
      N_grid = size(matter%El_EMFP_tot%E)
      do k = 1, N_grid
         write(FN,'(f16.4)',advance='no') matter%El_EMFP_tot%E(k)
         do i = 1, size(matter%Atoms) ! for all atoms
            if (matter%Atoms(i)%El_EMFP%L(k) > 0.0d0) then
               write(FN,'(es24.8)',advance='no') 1.0d0/matter%Atoms(i)%El_EMFP%L(k) ! elastic MFP
            else  ! zero CS
               write(FN,'(es24.8)',advance='no') 1.0d20  ! elastic MFP
            endif
         enddo
         write(FN,'(es24.8)') 1.0d0/matter%El_EMFP_tot%L(k)  ! total
      enddo
      call close_file('close', FN=FN)  ! module "Dealing_with_files"


      !-----------------------
      ! Photon attenuation length:
      text_var = 'OUTPUT_'
      write(text_var,'(a)') trim(adjustl(text_var))//trim(adjustl(matter%Name))//'_Photon_IMFP'

      file_name = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(text_var))//'.dat'
      file_photon_MFP = trim(adjustl(text_var))//'.dat'

      FN = 9997
      open(UNIT=FN, FILE = trim(adjustl(file_name)), status = 'new')

      ! Printout MFPs:
      ! 1) Create the comment lines:
      write(FN,'(a)',advance='no') '#Energy   '
      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         do j = 1, Nshl ! for all shells of this atom
            if ((i == 1) .and. (j == Nshl)) then ! valence
               ! skip and print at the end
               !write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))
            else ! core shell
               write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Name))//'-'//trim(adjustl(matter%Atoms(i)%Shell_name(j)))
            endif
         enddo ! j
      enddo ! i
      Nshl = size(matter%Atoms(1)%Ip)
      write(FN,'(a)') '   '//trim(adjustl(matter%Atoms(1)%Shell_name(Nshl)))//'  Total'

      ! 2) Save MFPs for all shells and total:
      N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
      do k = 1, N_grid
         write(FN,'(f16.4)',advance='no') matter%Atoms(1)%Ph_MFP(1)%E(k)
         ATOMS2:do i = 1, size(matter%Atoms) ! for all atoms
            Nshl = size(matter%Atoms(i)%Ip)
            SHELLS2:do j = 1, Nshl ! for all shells of this atom
               if ((i == 1) .and. (j == Nshl)) then ! valence
                  ! skip and print at the end
                  !write(FN,'(a)',advance='no') '   '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))
               else ! core shell
                  write(FN,'(es24.8)',advance='no') 1.0d0/matter%Atoms(i)%Ph_MFP(j)%L(k)  ! MFP for ionizing this shell
               endif
            enddo SHELLS2
         enddo ATOMS2
         Nshl = size(matter%Atoms(1)%Ip)
         write(FN,'(es24.8, es24.8)') 1.0d0/matter%Atoms(1)%Ph_MFP(Nshl)%L(k), 1.0d0/matter%Ph_MFP_tot%L(k)  ! total
      enddo ! k
      call close_file('close', FN=FN)  ! module "Dealing_with_files"


      !=======================================================
      ! Gnuplot the data:
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         call_slash = 'call '
         sh_cmd = '.cmd'
      else ! It is linux
         call_slash = './'
         sh_cmd = '.sh'
      endif

      N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
      !t0 = matter%Atoms(1)%Ph_MFP(1)%E(1)
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
      !call write_gnuplot_script_ending(FN, gnu_electron_MFP, 1)   ! below
      call write_gnuplot_script_ending_new(FN, gnu_electron_MFP, numpar%path_sep) ! module "Gnuplotting"
      close(FN)

      !--------------------
      ! 2) Photon attenuation lengths gnuplotting:
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
      !call write_gnuplot_script_ending(FN, gnu_photon_MFP, 1)  ! below
      call write_gnuplot_script_ending_new(FN, gnu_photon_MFP, numpar%path_sep) ! module "Gnuplotting"
      close(FN)

   endif
end subroutine printout_MFP_file



subroutine electronic_distribution_on_grid(Scell, numpar, tim)
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: tim ! [fs] timestep
   !----------------
   integer :: i, j, Nsiz, Nei
   real(8) :: N_steps

   ! Add the high-energy part of the distribution (obtained from MC):
   Scell%fe_on_grid = Scell%fe_on_grid + Scell%fe_high_on_grid
   Scell%fe_norm_on_grid = Scell%fe_norm_on_grid + Scell%fe_norm_high_on_grid

   ! Average over the number of time-steps it was collected over:
   N_steps = max( 1.0d0, dble(numpar%fe_aver_num) )   ! at least one step, to not change anything if nothing happened
   Scell%fe_on_grid = Scell%fe_on_grid / N_steps
   Scell%fe_norm_on_grid = Scell%fe_norm_on_grid / N_steps

   ! Now save the distribution in the file:
   call save_distribution_on_grid(numpar%FN_fe_on_grid, tim, Scell%E_fe_grid, Scell%fe_on_grid, Scell%fe_norm_on_grid)  ! below

   ! Reset the high-energy electron part for the next step:
   Scell%fe_on_grid = 0.0d0
   Scell%fe_high_on_grid = 0.0d0
   Scell%fe_norm_on_grid = 0.0d0
   Scell%fe_norm_high_on_grid = 0.0d0
   numpar%fe_aver_num = 0   ! to restart counting time-steps
end subroutine electronic_distribution_on_grid


subroutine save_distribution_on_grid(FN, tim, wr, fe, fe_norm)
   integer, intent(in) :: FN
   real(8), intent(in) :: tim
   real(8), dimension(:), intent(in) :: wr
   real(8), dimension(:), intent(in) :: fe, fe_norm
   integer i
   write(FN,'(a,f25.16)') '#', tim
   do i = 1, size(fe)
      write(FN,'(f25.16,es25.16E4,es25.16E4)') wr(i), fe(i), fe_norm(i)
   enddo
   write(FN,*) ''
   write(FN,*) ''
end subroutine save_distribution_on_grid




subroutine convolve_output(Scell, numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   !------------------
   integer i, FN, j, Nsiz
   logical file_exist, file_opened, file_named
   character(200) :: File_name, file_path
      
   SCL:do i = 1, size(Scell) ! for all supercells
      !if (Scell(i)%eps%tau > 0.0d0) then ! convolve output files:
         ! Subroutine "convolution" is from module "Little_subroutines"
         ! SIDENOTE: for some reason, using NEWUNIT creates 'unnamed' files
         ! at least in intel-fortran under Windows (had no problem under Linux!),
         ! therefore I have to use fixed number for the file unit,
         ! to be able to refer to its name within the "convolution" subroutine.
         ! The file names given here must exactly coinside with the names given 
         ! below in the subroutine "create_output_files".
         FN = 9999
         file_path = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))

         File_name = trim(adjustl(file_path))//'OUTPUT_optical_coefficients.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))   
            call convolution(FN, Scell(i)%eps%tau)       ! optical coefficients
            close(FN)
         endif

         File_name = trim(adjustl(file_path))//'OUTPUT_electron_properties.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))  
            call convolution(FN, Scell(i)%eps%tau) ! electron properties
            close(FN)
         endif

         File_name = trim(adjustl(file_path))//'OUTPUT_electron_entropy.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))
            call convolution(FN, Scell(i)%eps%tau) ! electron entropy
            close(FN)
         endif

         if (numpar%print_Ta) then
            File_name = trim(adjustl(file_path))//'OUTPUT_atomic_entropy.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
            if (file_exist) then
               open(UNIT=FN, FILE = trim(adjustl(File_name)))
               call convolution(FN, Scell(i)%eps%tau) ! atomic entropy
               close(FN)
            endif

            File_name = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
            if (file_exist) then
               open(UNIT=FN, FILE = trim(adjustl(File_name)))
               call convolution(FN, Scell(i)%eps%tau) ! atomic temperature
               close(FN)
            endif
         endif

         if (numpar%do_partial_thermal) then
            File_name = trim(adjustl(file_path))//'OUTPUT_electron_temperatures.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
            if (file_exist) then
               open(UNIT=FN, FILE = trim(adjustl(File_name)))
               call convolution(FN, Scell(i)%eps%tau) ! electron temperature
               close(FN)
            endif

            File_name = trim(adjustl(file_path))//'OUTPUT_electron_chempotentials.dat'
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
            if (file_exist) then
               open(UNIT=FN, FILE = trim(adjustl(File_name)))
               call convolution(FN, Scell(i)%eps%tau) ! electron chemical poential
               close(FN)
            endif
         endif

         File_name = trim(adjustl(file_path))//'OUTPUT_electron_hole_numbers.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))   
            call convolution(FN, Scell(i)%eps%tau)      ! numbers of particles
            close(FN)
         endif
         
         File_name = trim(adjustl(file_path))//'OUTPUT_energies.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))   
            call convolution(FN, Scell(i)%eps%tau)     ! energies
            close(FN)
         endif

         File_name = trim(adjustl(file_path))//'OUTPUT_temperatures.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))   
            call convolution(FN, Scell(i)%eps%tau) ! temperatures
            close(FN)
         endif

         ! Atomic masks for sectional displacements:
         if (allocated(Scell(i)%Displ)) then
            Nsiz = size(Scell(i)%Displ)   ! how many masks
            do j = 1, Nsiz    ! for all masks
               File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(i)%Displ(j)%mask_name))//'.dat'
               inquire(file=trim(adjustl(File_name)),exist=file_exist)
               if (file_exist) then
                  open(UNIT=FN, FILE = trim(adjustl(File_name)))
                  call convolution(FN, Scell(i)%eps%tau) ! displacements
                  close(FN)
               endif
            enddo ! i
         endif
         
         File_name = trim(adjustl(file_path))//'OUTPUT_pressure_and_stress.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))   
            call convolution(FN, Scell(i)%eps%tau) ! pressure and stress tensor
            close(FN)
         endif

         File_name = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))    
            call convolution(FN, Scell(i)%eps%tau)   ! core holes
            close(FN)
         endif

         File_name = trim(adjustl(file_path))//'OUTPUT_supercell.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         if (file_exist) then
            open(UNIT=FN, FILE = trim(adjustl(File_name)))    
            call convolution(FN, Scell(i)%eps%tau)    ! supercell
            close(FN)
         endif
      !endif ! tau>0
   enddo SCL
end subroutine convolve_output



subroutine write_holes(FN, time, matter, Scell)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Solid), intent(in) :: matter	! Material parameters
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   integer i, j, Nshl, Na, temp1, temp2
   write(FN, '(f25.16)', advance='no') time
   Na = size(matter%Atoms)
   ATOMS:do i = 1, Na ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      SHELLS:do j = 1, Nshl ! for all shells of this atom
         if (Na == 1) then
            if (j == matter%Atoms(Na)%sh) then ! last shell:
               write(FN, '(es25.16)') Scell%MChole(i)%Noh(j)
            else ! not last shell:
               write(FN, '(es25.16)', advance='no') Scell%MChole(i)%Noh(j)
            endif
         else
!             if ((i == 1) .and. (j == Nshl)) then ! VB:
!                ! skip it
!             else ! atomic shell:
               temp1 = Na    ! default
               temp2 = matter%Atoms(Na)%sh
               if (temp2 == 0) then
                  temp1 = Na - 1
                  temp2 = matter%Atoms(temp1)%sh
               endif
               if ((i == temp1) .and. (j == temp2)) then ! last shell:
                  write(FN, '(es25.16)') Scell%MChole(i)%Noh(j)
               else ! not last shell:
                  write(FN, '(es25.16)', advance='no') Scell%MChole(i)%Noh(j)
               endif
!             endif
         endif
      enddo SHELLS
   enddo ATOMS
   !write(FN, '(a)', advance='yes') ' '
end subroutine write_holes



subroutine write_optical_all_hw(FN, tim, eps)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: tim
   type(Drude), intent(in) :: eps	! epsylon, Drude dielectric function and its parameters
   integer i
   do i = 1,size(eps%Eps_hw,2)
      write(FN, '(f25.16, es25.16, es25.16, es25.16, f25.16, f25.16, f25.16, f25.16, f25.16, es25.16, es25.16, &
                  es25.16, es25.16, es25.16, es25.16, es25.16)') &
                  eps%Eps_hw(1,i), eps%Eps_hw(2,i), eps%Eps_hw(3,i), &
                  eps%Eps_hw(4,i), eps%Eps_hw(5,i), eps%Eps_hw(6,i), &
                  eps%Eps_hw(7,i), eps%Eps_hw(8,i), eps%Eps_hw(9,i), &
                  eps%Eps_hw(10,i), eps%Eps_hw(11:,i)
   enddo
   write(FN, '(a)')
   write(FN, '(a)')
   ! Reminder:
!    Eps_hw(1,i) = hw     ! energy
!    Eps_hw(2,i) = Re_eps ! real part of CDF
!    Eps_hw(3,i) = Im_eps ! imaginary part of CDF
!    Eps_hw(4,i) = LF  ! loss function
!    Eps_hw(5,i) = R   ! reflectivity
!    Eps_hw(6,i) = T   ! transmission
!    Eps_hw(7,i) = A   ! absorption
!    Eps_hw(8,i) = n   ! optical n
!    Eps_hw(9,i) = k   ! optical k
!    Eps_hw(10,i) = k   ! dc-conductivity
!    Eps_hw(11,i) = Re_E_xx
!    Eps_hw(12,i) = Im_E_xx
!    Eps_hw(13,i) = Re_E_yy
!    Eps_hw(14,i) = Im_E_yy
!    Eps_hw(15,i) = Re_E_zz
!    Eps_hw(16,i) = Im_E_zz
end subroutine write_optical_all_hw


subroutine write_optical_coefs(FN, tim, eps)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: tim
   type(Drude), intent(in) :: eps	! epsylon, Drude dielectric function and its parameters
   write(FN, '(f25.16,es25.16,es25.16,es25.16,f25.16,f25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') &
      tim, eps%R, eps%T, eps%A, eps%n, eps%k, eps%ReEps, eps%ImEps, eps%dc_cond, eps%Eps_xx, eps%Eps_yy, eps%Eps_zz
end subroutine write_optical_coefs


subroutine write_PCF(FN, atoms, matter, Scell, NSC)
   integer, intent(in) :: FN	! file number
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   type(Solid), intent(inout) :: matter	! Material parameters
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   integer i

   do i = 1,size(matter%PCF,2)
      write(FN, '(f25.16,es25.16)') matter%PCF(1,i), matter%PCF(2,i)
   enddo
   write(FN, '(a)') ''
   write(FN, '(a)') ''
end subroutine write_PCF


subroutine save_distribution(FN, numpar, Scell, tim, wr, fe, fe_eq, fe_eq_VB, fe_eq_CB)
   integer, intent(in) :: FN
   type(Numerics_param), intent(in) :: numpar  ! numerical parameters
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   real(8), intent(in) :: tim
   real(8), dimension(:), intent(in) :: wr
   real(8), dimension(:), intent(in) :: fe
   real(8), dimension(:), allocatable, intent(in) :: fe_eq
   real(8), dimension(:), intent(in), optional :: fe_eq_VB, fe_eq_CB ! equivalent distr. in VB and CB
   integer i, j, k
   write(FN,'(a,f25.16)') '#', tim
   if (allocated(fe_eq)) then ! there is equivalent-temperature Fermi distribution
      if (present(fe_eq_VB) .and. present(fe_eq_CB)) then
         if (numpar%save_fe_orb) then  ! orbital-resolved distributions
            do i = 1, size(fe)
               write(FN,'(f25.16,f25.16,f25.16,f25.16,f25.16)', advance='no') wr(i), fe(i), fe_eq(i), fe_eq_VB(i), fe_eq_CB(i)
               ! Add orbital-resolved distributions:
               do j = 1, size(Scell%Orb_data)  ! all kinds of atoms
                  do k = 1, size(Scell%Orb_data(j)%Ne)  ! all orbital types
                     write(FN,'(f25.16)', advance='no') Scell%Orb_data(j)%fe(k,i)
                  enddo ! k
               enddo ! j
               write(FN,'(a)') ''
            enddo ! i
         else  ! no orbital-resolved distributions
            do i = 1, size(fe)
               write(FN,'(f25.16,f25.16,f25.16,f25.16,f25.16)') wr(i), fe(i), fe_eq(i), fe_eq_VB(i), fe_eq_CB(i)
            enddo
         endif ! numpar%save_fe_orb
      else  ! without band-resolved part
         if (numpar%save_fe_orb) then  ! orbital-resolved distributions
            do i = 1, size(fe)
               write(FN,'(f25.16,f25.16,f25.16)', advance='no') wr(i), fe(i), fe_eq(i)
               ! Add orbital-resolved distributions:
               do j = 1, size(Scell%Orb_data)  ! all kinds of atoms
                  do k = 1, size(Scell%Orb_data(j)%Ne)  ! all orbital types
                     write(FN,'(f25.16)', advance='no') Scell%Orb_data(j)%fe(k,i)
                  enddo ! k
               enddo ! j
               write(FN,'(a)') ''
            enddo ! i
         else  ! no orbital-resolved distributions
            do i = 1, size(fe)
               write(FN,'(f25.16,f25.16,f25.16)') wr(i), fe(i), fe_eq(i)
            enddo
         endif ! numpar%save_fe_orb
      endif ! (present(fe_eq_VB) .and. present(fe_eq_CB))
   else  ! fe is Fermi, no equivalent distribution needed
      if (numpar%save_fe_orb) then  ! orbital-resolved distributions
         do i = 1, size(fe)
            write(FN,'(f25.16,f25.16)', advance='no') wr(i), fe(i)
            ! Add orbital-resolved distributions:
            do j = 1, size(Scell%Orb_data)  ! all kinds of atoms
               do k = 1, size(Scell%Orb_data(j)%Ne)  ! all orbital types
                  write(FN,'(f25.16)', advance='no') Scell%Orb_data(j)%fe(k,i)
               enddo ! k
            enddo ! j
            write(FN,'(a)') ''
         enddo ! i
      else  ! no orbital-resolved distributions
         do i = 1, size(fe)
            write(FN,'(f25.16,f25.16)') wr(i), fe(i)
         enddo
      endif ! numpar%save_fe_orb
   endif ! (allocated(fe_eq))
   write(FN,*) ''
   write(FN,*) ''
end subroutine save_distribution



subroutine save_atomic_distribution(FN, numpar, Scell, tim, Ea_grid, fa, fa_eq)
   integer, intent(in) :: FN
   type(Numerics_param), intent(in) :: numpar  ! numerical parameters
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   real(8), intent(in) :: tim
   real(8), dimension(:), intent(in) :: Ea_grid
   real(8), dimension(:), intent(in) :: fa
   real(8), dimension(:), allocatable, intent(in), optional :: fa_eq
   integer i, j, k
   logical :: there_is_eq

   if (present(fa_eq)) then
      if (allocated(fa_eq)) then ! there is equivalent-temperature Maxwell distribution
         there_is_eq = .true.
      else
         there_is_eq = .false.
      endif
   else
      there_is_eq = .false.
   endif

   write(FN,'(a,f25.16)') '#', tim
   if (there_is_eq) then ! there is equivalent-temperature Maxwell distribution
      do i = 1, size(fa)
         write(FN,'(f25.16,f25.16,f25.16)') Ea_grid(i), fa(i), fa_eq(i)
      enddo
   else  ! fe is Maxwell, no equivalent distribution needed
      do i = 1, size(fa)
         write(FN,'(f25.16,f25.16)') Ea_grid(i), fa(i)
      enddo
   endif ! (allocated(fe_eq))
   write(FN,*) ''
   write(FN,*) ''
end subroutine save_atomic_distribution



subroutine save_energy_levels(FN, tim, wr)
   integer, intent(in) :: FN
   real(8), intent(in) :: tim
   real(8), dimension(:), intent(in) :: wr
   integer i
   write(FN,'(f25.16)', advance='no') tim
   do i = 1, size(wr)
      write(FN,'(f25.16)', advance='no') wr(i)
   enddo
   write(FN,*) '' 
end subroutine save_energy_levels



subroutine save_DOS(FN, tim, DOS, partial_DOS)
   integer, intent(in) :: FN
   real(8), intent(in) :: tim
   real(8), dimension(:,:), intent(in) :: DOS
   real(8), dimension(:,:,:), intent(in), optional :: partial_DOS
   integer i, Nsiz, j, Nat, k, Ntype
   Nsiz = size(DOS,2)
   do i = 1, Nsiz
      if (present(partial_DOS)) then
      Nat = size(partial_DOS,1)
      Ntype = size(partial_DOS,2)
         write(FN,'(f25.16,f25.16)',advance='no') DOS(1,i), DOS(2,i)
         do j = 1, Nat
            do k = 1, Ntype
               write(FN,'(f25.16)', advance='no') partial_DOS(j,k,i)
            enddo
         enddo
         write(FN,'(a)') ''
      else
         write(FN,'(f25.16,f25.16)') DOS(1,i), DOS(2,i)
      endif
   enddo
   write(FN,*) ''
   write(FN,*) ''
end subroutine save_DOS


subroutine save_nearest_neighbors(FN, Scell, NSC, tim)
   integer, intent(in) :: FN    ! file to write into
   real(8), intent(in) :: tim   ! current simulation time
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   !--------------------------
   integer :: i
   real(8) :: Nat, NofN, NoNN(7)
   
   Nat = dble(size(Scell(NSC)%MDatoms)) ! number of atoms
   NofN = dble(SUM(Scell(NSC)%Near_neighbors_user(:)))/Nat    ! Average number of nearest neighbors
   do i = 1, 7
      NoNN(i) = dble( COUNT(Scell(NSC)%Near_neighbors_user == i-1) ) / Nat
   enddo
   
   write(FN, '(es25.16,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6)') tim, NofN, NoNN(:)
end subroutine save_nearest_neighbors


subroutine save_nearest_neighbors_element(FN, numpar, Scell, NSC, tim) ! element-specific nearest neighbors
   integer, dimension(:), intent(in) :: FN    ! file to write into
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   real(8), intent(in) :: tim   ! current simulation time
   !--------------------------
   integer :: i, j
   real(8) :: Nat, NofN, NoNN(7)

   do i = 1, size(numpar%NN_radii) ! for all requested elements
      write(FN(i), '(es25.16,f10.6)', advance='no') tim, Scell(NSC)%NN_numbers(i)%total

      do j = 1, size(Scell(NSC)%NN_numbers(i)%NNN)     ! for all elements in the compound:
         write(FN(i), '(f10.6)', advance='no') Scell(NSC)%NN_numbers(i)%NNN(j)
      enddo ! j
      write(FN(i), '(a)') '' ! end line
   enddo ! i
end subroutine save_nearest_neighbors_element


subroutine write_atomic_xyz(FN, atoms, matter, Supce, print_mass, print_charge, print_Ekin)
   integer, intent(in) :: FN	! file number
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   type(Solid), intent(in) :: matter	! Material parameters
   real(8), dimension(3,3), intent(in) :: Supce	! [A]  supercell vectors [a(x,y,z),b(x,y,z),c(x,y,z)]
   logical, intent(in), optional :: print_mass, print_charge, print_Ekin ! flags what to printout
   !-------------------------------------
   integer i
   character(10) :: Numb_out
   logical :: do_mass, do_charge, do_Ekin

   ! Check for additional properties provided:
   if (present(print_mass)) then
      do_mass = print_mass
   else
      do_mass = .false.
   endif

   if (present(print_charge)) then
      do_charge = print_charge
   else
      do_charge = .false.
   endif

   if (present(print_Ekin)) then
      do_Ekin = print_Ekin
   else
      do_Ekin = .false.
   endif


   ! Write out the data block:
   write(Numb_out, '(i10)') size(atoms)
   write(FN, '(a)') trim(adjustl(Numb_out))
   write(FN, '(a,f,f,f,f,f,f,f,f,f,a)', advance='no') 'Lattice="', Supce(1,1), Supce(1,2), Supce(1,3), &
                                                     Supce(2,1), Supce(2,2), Supce(2,3), &
                                                     Supce(3,1), Supce(3,2), Supce(3,3), '" Properties=species:S:1:pos:R:3'

   ! optional additional data:
   if (do_charge) then
      write(FN, '(a)', advance='no') ':charge:R:1'
   endif
   if (do_Ekin) then
      write(FN, '(a)', advance='no') ':kinetic_energy:R:1'
   endif
   if (do_mass) then
      write(FN, '(a)', advance='no') ':mass:R:1'
   endif
   write(FN, '(a)') ''  ! to end the line


   ! Atomic data block:
   do i = 1, size(atoms)
      write(FN, '(a,es25.16,es25.16,es25.16)' , advance='no') trim(adjustl(matter%Atoms(atoms(i)%KOA)%Name)), &
                                                               atoms(i)%R(1), atoms(i)%R(2), atoms(i)%R(3)
      if (do_mass) write(FN, '(es25.16)', advance='no') matter%Atoms(atoms(i)%KOA)%Ma
      !if (do_charge) write(FN, '(es25.16)', advance='no') matter%Atoms(atoms(i)%KOA)%mulliken_q
      if (do_charge) write(FN, '(es25.16)', advance='no') atoms(i)%q
      if (do_Ekin) write(FN, '(es25.16)', advance='no') atoms(i)%Ekin
      write(FN, '(a)') ''  ! to end the line
   enddo
end subroutine write_atomic_xyz


subroutine write_atomic_cif(FN_out, Supce, atoms, matter, tim)
   integer, intent(in) :: FN_out	! file number to write to
   real(8), dimension(3,3), intent(in) :: Supce	! [A]  supercell vectors [a(x,y,z),b(x,y,z),c(x,y,z)]
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   type(Solid), intent(in) :: matter	! Material parameters
   real(8), intent(in) :: tim	! [fs] timestep
   !----------------------------------
   real(8) :: a, b, c, alpha, beta, gamma
   integer :: i, Nat
   character(12) :: i_char
   Nat = size(atoms)
!    a = DSQRT(Supce(1,1)*Supce(1,1) + Supce(1,2)*Supce(1,2) + Supce(1,3)*Supce(1,3))
!    b = DSQRT(Supce(2,1)*Supce(2,1) + Supce(2,2)*Supce(2,2) + Supce(2,3)*Supce(2,3))
!    c = DSQRT(Supce(3,1)*Supce(3,1) + Supce(3,2)*Supce(3,2) + Supce(3,3)*Supce(3,3))
   a = DSQRT(Supce(1,1)*Supce(1,1) + Supce(2,1)*Supce(2,1) + Supce(3,1)*Supce(3,1))
   b = DSQRT(Supce(1,2)*Supce(1,2) + Supce(2,2)*Supce(2,2) + Supce(3,2)*Supce(3,2))
   c = DSQRT(Supce(1,3)*Supce(1,3) + Supce(2,3)*Supce(2,3) + Supce(3,3)*Supce(3,3))
!    alpha = (Supce(3,1)*Supce(2,1) + Supce(3,2)*Supce(2,2) + Supce(3,3)*Supce(2,3))/(c*b)
!    beta = (Supce(1,1)*Supce(3,1) + Supce(1,2)*Supce(3,2) + Supce(1,3)*Supce(3,3))/(a*c)
!    gamma = (Supce(1,1)*Supce(2,1) + Supce(1,2)*Supce(2,2) + Supce(1,3)*Supce(2,3))/(a*b)
   alpha = (Supce(1,3)*Supce(1,2) + Supce(2,3)*Supce(2,2) + Supce(3,3)*Supce(3,2))/(c*b)
   beta = (Supce(1,1)*Supce(1,3) + Supce(2,1)*Supce(2,3) + Supce(3,1)*Supce(3,3))/(a*c)
   gamma = (Supce(1,1)*Supce(1,2) + Supce(2,1)*Supce(2,2) + Supce(3,1)*Supce(3,2))/(a*b)

   write(i_char,'(f12.3)') tim
   write(FN_out, '(a)') 'Data_for_x_ray_diffraction_'//trim(adjustl(i_char))
   write(FN_out, '(a,f)') '_cell_length_a', a
   write(FN_out, '(a,f)') '_cell_length_b', b
   write(FN_out, '(a,f)') '_cell_length_c', c
   write(FN_out, '(a,f)') '_cell_angle_alpha', ACOS(alpha)*180.0d0/g_Pi
   write(FN_out, '(a,f)') '_cell_angle_beta', ACOS(beta)*180.0d0/g_Pi
   write(FN_out, '(a,f)') '_cell_angle_gamma', ACOS(gamma)*180.0d0/g_Pi
   write(FN_out, '(a)') "_symmetry_space_group_name_H-M 'P 1' "
   write(FN_out, '(a)') "_symmetry_Int_Tables_number    '1' "
   write(FN_out, '(a)') "loop_ "
   write(FN_out, '(a)') '_symmetry_equiv_pos_site_id'
   write(FN_out, '(a)') '_symmetry_equiv_pos_as_xyz'
   write(FN_out, '(a)') "1     'x, y, z' "
   write(FN_out, '(a)') "loop_ "
   write(FN_out, '(a)') "_atom_site_label"
   write(FN_out, '(a)') "_atom_site_type_symbol"
   write(FN_out, '(a)') "_atom_site_fract_x"
   write(FN_out, '(a)') "_atom_site_fract_y"
   write(FN_out, '(a)') "_atom_site_fract_z"
   do i = 1, Nat
      write(i_char,'(i12)') i
      write(FN_out, '(a,es25.16,es25.16,es25.16)')  trim(adjustl(i_char))//'	'//trim(adjustl(matter%Atoms(atoms(i)%KOA)%Name)), &
         atoms(i)%S(1), atoms(i)%S(2), atoms(i)%S(3)
   enddo
end subroutine write_atomic_cif



subroutine write_electron_properties(FN, time, Scell, NSC, Ei, matter, numpar, FN_Ce, FN_kappa, FN_kappa_dyn, FN_Se, FN_Te, FN_mu)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: Ei	! energy levels
   type(Solid), intent(in) :: matter	! Material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   ! File numbers for band-resolved Ce and kappa, electron entropy, electron temperatures and chem.potentials:
   integer, intent(in) :: FN_Ce, FN_kappa, FN_kappa_dyn, FN_Se, FN_Te, FN_mu
   !------------------------
   integer i, Nat, n_at, Nsiz, norb, N_types, i_at, i_types, i_G1

   ! Write electron properties:
   write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)', advance='no') time, &
      !Scell(NSC)%Ne_low/dble(Scell(NSC)%Ne)*100.0d0, Scell(NSC)%mu, Scell(NSC)%E_gap, Scell(NSC)%Ce, Scell(NSC)%G_ei, &
      Scell(NSC)%Ne_low/Scell(NSC)%Na, Scell(NSC)%mu, Scell(NSC)%E_gap, Scell(NSC)%Ce, Scell(NSC)%G_ei, &
      Scell(NSC)%E_VB_bottom, Scell(NSC)%E_VB_top, Scell(NSC)%E_bottom, Scell(NSC)%E_top
   Nat = size(matter%Atoms(:)) ! number of elements
   do i = 1, Nat    ! index starting from 11
      !write(FN,'(es25.16)', advance='no') (matter%Atoms(i)%NVB - matter%Atoms(i)%mulliken_Ne)
      write(FN,'(es25.16)', advance='no') matter%Atoms(i)%mulliken_q
   enddo
   write(FN,'(a)') ''

   ! Write band-resolved electron heat capacity:
   n_at = size(Scell(NSC)%MDatoms) ! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
   norb =  Nsiz/n_at ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"
   ! Total Ce:
   write(FN_Ce, '(es25.16,es25.16)', advance='no') time, Scell(NSC)%Ce
   ! All shells resolved:
   do i_at = 1, Nat
      do i_types = 1, N_types
         i_G1 = (i_at-1) * N_types + i_types
         write(FN_Ce,'(es25.16)',advance='no') Scell(NSC)%Ce_part(i_G1)
      enddo   ! i_types
   enddo ! i_at
   write(FN_Ce,'(a)') ''

   ! Write electron heat conductivity if requesed:
   if (numpar%do_kappa) then
      do i = 1, size(Scell(NSC)%kappa_e_vs_Te)  ! electron temperature dependence
         write(FN_kappa, '(es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') Scell(NSC)%kappa_Te_grid(i), &
            1.0d0/ ( 1.0d0/Scell(NSC)%kappa_e_vs_Te(i) + 1.0d0/Scell(NSC)%kappa_ee_vs_Te(i) ), & ! total conductivity
            Scell(NSC)%kappa_e_vs_Te(i), Scell(NSC)%kappa_ee_vs_Te(i), & ! electron-phonon, electron-electron contributions
            Scell(NSC)%kappa_mu_grid(i), Scell(NSC)%kappa_Ce_grid(i) ! chem.potential, electron heat capacity
      enddo
      write(FN_kappa, '(a)')
      write(FN_kappa, '(a)')
   endif
   ! The dynamic version of kappa:
   if (numpar%do_kappa_dyn) then
      write(FN_kappa_dyn, '(es25.16, es25.16, es25.16, es25.16)') time, Scell(NSC)%kappa_e, & ! total conductivity
            Scell(NSC)%kappa_e_vs_Te(1), Scell(NSC)%kappa_ee_vs_Te(1) ! electron-phonon, electron-electron contributions
   endif


   if (numpar%do_partial_thermal) then
      ! Write electron entropy:
      write(FN_Se, '(es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') time, Scell(NSC)%Se, Scell(NSC)%Se_eq, &
               Scell(NSC)%Se_VB, Scell(NSC)%Se_eq_VB, Scell(NSC)%Se_CB, Scell(NSC)%Se_eq_CB
      ! Write electron temperatures and chemical potentials:
      write(FN_Te, '(es25.16, es25.16, es25.16, es25.16)') time, Scell(NSC)%Te, Scell(NSC)%Te_VB, Scell(NSC)%Te_CB
      write(FN_mu, '(es25.16, es25.16, es25.16, es25.16)') time, Scell(NSC)%mu, Scell(NSC)%mu_VB, Scell(NSC)%mu_CB
   else
      ! Write electron entropy:
      write(FN_Se, '(es25.16, es25.16, es25.16)') time, Scell(NSC)%Se, Scell(NSC)%Se_eq
!       ! Write electron temperatures and chemical potentials:
!       write(FN_Te, '(es25.16, es25.16)') time, Scell(NSC)%Te
!       write(FN_mu, '(es25.16, es25.16)') time, Scell(NSC)%mu
   endif

end subroutine write_electron_properties


subroutine write_atomic_properties(time, Scell, NSC, matter, numpar) ! atomic parameters
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in) :: matter	! Material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !------------------------------

   ! Atomic temperatures (various definitions):
   if (numpar%print_Ta) then

      ! Atomic entropy:
      write(numpar%FN_Sa, '(es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') time, Scell(NSC)%Sa, Scell(NSC)%Sa_eq, &
      Scell(NSC)%Sa_eq_num, Scell(NSC)%Sa_conf, Scell(NSC)%Sa_tot

      ! kinetic; entropic; distributional; fluctuational; "potential"; configurational etc.:
      write(numpar%FN_Ta, '(es25.16, &
            es25.16, es25.16, es25.16, es25.16, &
            es25.16, es25.16, es25.16, es25.16, &
            es25.16, es25.16, es25.16)') &
      time, &
      Scell(NSC)%Ta_var(1), Scell(NSC)%Ta_var(2), Scell(NSC)%Ta_var(3), Scell(NSC)%Ta_var(4), &
      Scell(NSC)%Ta_var(5), Scell(NSC)%Ta_var(6), Scell(NSC)%Ta_var(7), Scell(NSC)%Ta_var(8), &
      Scell(NSC)%Fv, Scell(NSC)%Tconf, Scell(NSC)%Tconf2

      ! partial temperatures along X,Y,Z:
      write(numpar%FN_Ta_part, '(es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') time, &
      Scell(NSC)%Ta_r_var(1:6)
   endif
end subroutine write_atomic_properties


subroutine write_coulping_header(FN, Scell, NSC, matter, numpar)
   integer, intent(in) :: FN	! file number
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in) :: matter	! Material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   !--------------------------------
   integer :: i, j, N_at, N_types, Nsiz, i_at, i_at2, i_types, i_types2, i_G1, i_G2, nat, norb
   real(8) :: G_part
   character(2) :: chtemp1, chtemp2
   
   N_at = matter%N_KAO    ! number of kinds of atoms
   ! Find number of orbitals per atom:
   nat = size(Scell(1)%MDatoms) ! number of atoms
   Nsiz = size(Scell(1)%Ha,1) ! total number of orbitals
   norb =  Nsiz/nat ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

   ! Total coupling:
   write(FN, '(a)', advance='no') ' #Time   Total   '
   ! Atoms resolved:
   do i = 1, N_at
      do j = 1, N_at
         write(FN,'(a,a,a)', advance='no') trim(adjustl(matter%Atoms(i)%Name)), '-', trim(adjustl(matter%Atoms(j)%Name))//' '
      enddo
   enddo
   ! All shells resolved:
   do i_at = 1, N_at
      do i_types = 1, N_types
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"  
         do i_at2 = 1, N_at
            do i_types2 = 1, N_types
               chtemp2 = name_of_orbitals(norb, i_types2) ! module "Little_subroutines"
               write(FN,'(a,a,a)',advance='no') trim(adjustl(matter%Atoms(i_at)%Name))//'_'//trim(adjustl(chtemp1)), '--', &
                                        trim(adjustl(matter%Atoms(i_at2)%Name))//'_'//trim(adjustl(chtemp2))//'  '
            enddo   ! i_types2
         enddo ! i_at2
      enddo   ! i_types
   enddo ! i_at
   write(FN,'(a)') ''
end subroutine write_coulping_header


subroutine write_Ce_header(FN, Scell, NSC, matter)
   integer, intent(in) :: FN	! file number
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in) :: matter	! Material parameters
   !--------------------------------
   integer :: N_at, N_types, Nsiz, i_at, i_types, nat, norb
   character(2) :: chtemp1

   N_at = matter%N_KAO    ! number of kinds of atoms
   ! Find number of orbitals per atom:
   nat = size(Scell(NSC)%MDatoms) ! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
   norb =  Nsiz/nat ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

   ! Total Ce:
   write(FN, '(a)', advance='no') ' #Time   Total   '
   ! All shells resolved:
   do i_at = 1, N_at
      do i_types = 1, N_types
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
         write(FN,'(a)',advance='no') trim(adjustl(matter%Atoms(i_at)%Name))//'_'//trim(adjustl(chtemp1))//'   '
      enddo   ! i_types
   enddo ! i_at
   write(FN,'(a)') ''
end subroutine write_Ce_header


subroutine write_orb_resolved_header(FN, Scell, matter)
   integer, intent(in) :: FN	! file number
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter	! Material parameters
   !--------------------------------
   integer :: N_at, N_types, Nsiz, i_at, i_types, nat, norb
   character(2) :: chtemp1

   N_at = matter%N_KAO    ! number of kinds of atoms
   ! Find number of orbitals per atom:
   nat = size(Scell%MDatoms) ! number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   norb =  Nsiz/nat ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"

   ! Ne:
   write(FN, '(a)', advance='no') ' #Time   Total_Ne   '
   ! All shells resolved:
   do i_at = 1, N_at
      do i_types = 1, N_types
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
         write(FN,'(a)',advance='no') trim(adjustl(matter%Atoms(i_at)%Name))//'_'//trim(adjustl(chtemp1))//'   '
      enddo   ! i_types
   enddo ! i_at

   ! Ee:
   write(FN, '(a)', advance='no') 'Total_Ee   '
   ! All shells resolved:
   do i_at = 1, N_at
      do i_types = 1, N_types
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"
         write(FN,'(a)',advance='no') 'Ee:'//trim(adjustl(matter%Atoms(i_at)%Name))//'_'//trim(adjustl(chtemp1))//'   '
      enddo   ! i_types
   enddo ! i_at

   write(FN,'(a)') ''
end subroutine write_orb_resolved_header


subroutine write_orb_resolved(FN, time, Scell, matter) ! orbital-resolved electronic data
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter	! Material parameters
   !--------------------------------
   real(8) :: Ne_tot, Ee_tot
   integer :: Nat, n_at, N_types, Nsiz, i_at, i_types, norb
   character(2) :: chtemp1

   ! Write band-resolved electron heat capacity:
   Nat = matter%N_KAO      ! number of kinds of atoms
   n_at = size(Scell%MDatoms) ! number of atoms
   Nsiz = size(Scell%Ha,1) ! total number of orbitals
   norb =  Nsiz/n_at ! orbitals per atom
   ! Find number of different orbital types:
   N_types = number_of_types_of_orbitals(norb)  ! module "Little_subroutines"
   ! Total Ne and Ee:
   Ne_tot = 0.0d0 ! to start with
   Ee_tot = 0.0d0 ! to start with
   do i_at = 1, Nat
      do i_types = 1, N_types
         Ne_tot = Ne_tot  + Scell%Orb_data(i_at)%Ne(i_types)   ! [1/atom]
         Ee_tot = Ee_tot  + Scell%Orb_data(i_at)%Ee(i_types)   ! [eV/atom]
      enddo ! i_types
   enddo ! i_at

   ! Electron densities:
   write(FN, '(f16.8,f16.8)', advance='no') time, Ne_tot
   ! All shells resolved:
   do i_at = 1, Nat
      do i_types = 1, N_types
         write(FN,'(f16.8)',advance='no') Scell%Orb_data(i_at)%Ne(i_types)
      enddo   ! i_types
   enddo ! i_at
   ! Electron energies:
   write(FN, '(f16.8)', advance='no') Ee_tot
   ! All shells resolved:
   do i_at = 1, Nat
      do i_types = 1, N_types
         write(FN,'(f16.8)',advance='no') Scell%Orb_data(i_at)%Ee(i_types)
      enddo   ! i_types
   enddo ! i_at
   write(FN,'(a)') ''
end subroutine write_orb_resolved


subroutine write_coulping(FN, time, Scell, NSC, numpar)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
!    type(Solid), intent(in) :: matter	! Material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   !--------------------------------
   integer :: i, j, N_at, N_types, Nsiz
   real(8) :: G_part
   
   N_at = size(numpar%DOS_weights,1)    ! number of kinds of atoms
   N_types = size(numpar%DOS_weights,2) ! number of atomic shells (basis set size)
   Nsiz = size(Scell(NSC)%G_ei_partial, 1)
   ! Total coupling:
   write(FN, '(es25.16, es25.16)', advance='no') time, Scell(NSC)%G_ei
   ! Atoms resolved:
   do i = 1, N_at
      do j = 1, N_at
         write(FN,'(es25.16)', advance='no') SUM(Scell(NSC)%G_ei_partial((i-1)*N_types+1:(i-1)*N_types+N_types, &
                                                                         (j-1)*N_types+1:(j-1)*N_types+N_types) )
      enddo
   enddo
   ! Partial:
   do i = 1, Nsiz
      do j = 1, Nsiz
         write(FN,'(es25.16)', advance='no') Scell(NSC)%G_ei_partial(i,j)
      enddo
   enddo
   write(FN,'(a)') ''
end subroutine write_coulping



subroutine write_atomic_relatives(FN, atoms)
   integer, intent(in) :: FN	! file number
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   integer i
   do i = 1, size(atoms)
      write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') &
         atoms(i)%R(1), atoms(i)%R(2), atoms(i)%R(3), atoms(i)%V(1), atoms(i)%V(2), atoms(i)%V(3), &
         atoms(i)%S(1), atoms(i)%S(2), atoms(i)%S(3), atoms(i)%SV(1), atoms(i)%SV(2), atoms(i)%SV(3)
   enddo
   write(FN, '(a)') ''
   write(FN, '(a)') ''
end subroutine write_atomic_relatives



subroutine write_numbers(FN, time, Scell)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   !write(FN,'(f25.16,f25.16,es,es25.16,es25.16,es25.16,es25.16)') time, Scell%Ne_low/dble(Scell%Na), Scell%Ne_CB/dble(Scell%Na), &
   write(FN,'(f25.16,f25.16,es,es25.16,es25.16,es25.16,es25.16)') time, &
      (Scell%Ne_low-Scell%Ne_CB)/dble(Scell%Na), &    ! valence-band electrons (below E_fermi)
      Scell%Ne_CB/dble(Scell%Na), &                   ! conduction-band electrons (above E_fermi)
      Scell%Ne_high/dble(Scell%Na), &                 ! high-energy electrons (in MC)
      Scell%Nh/dble(Scell%Na), &                      ! all core holes
      (dble(Scell%Ne)-(Scell%Ne_low+Scell%Ne_high-Scell%Nh))/dble(Scell%Na), &   ! error in particle conservation
      Scell%Nph/dble(Scell%Na)                        ! photons
end subroutine write_numbers



subroutine write_pressure(FN, time, Pressure, Stress)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   real(8), intent(in) :: Pressure
   real(8), dimension(3,3), intent(in) :: Stress
  write(FN,'(es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') time, Pressure, Stress(:,:)
end subroutine write_pressure


subroutine write_super_cell(FN, time, Scell)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   write(FN,'(es25.16, f18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, &
        es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6, es18.6)') &
        time, Scell%V, Scell%supce(:,:), Scell%Vsupce(:,:)
end subroutine write_super_cell


subroutine write_energies(FN, time, nrg)
   integer, intent(in) :: FN	! file number to write to
   real(8), intent(in) :: time	! [fs]
   type(Energies), intent(in) :: nrg

   write(FN, '(es25.16, es25.16, es25.16, es26.16E4, es25.16, es25.16, es25.16, es25.16, es26.16E4, 1X, es25.16E4)') time, &
   nrg%E_tot, &   ! Band energy
   nrg%Eh_tot, &  ! Energy of holes
   nrg%At_pot + nrg%E_vdW + nrg%E_coul_scc + nrg%E_coul + nrg%E_expwall, & ! potential energy (incl. short-range, vdW, Coulomb, etc.)
   nrg%At_kin, &  ! Kinetic energy of atoms
   nrg%Total, &   ! Total atomic energy
   nrg%Total + nrg%E_supce + nrg%El_high, &   ! Energy of atoms (incl. supercell) and electrons
   nrg%Total + nrg%E_supce + nrg%El_high + nrg%Eh_tot, & ! Total energy (incl. holes)
   nrg%E_vdW, &   ! van der Waals
   nrg%E_expwall  ! Short-range repulsive
   !print*, nrg%Total, nrg%E_supce, nrg%El_high
end subroutine write_energies


subroutine write_output_file_one(FN, par1, par2, par3, par4, par5, par6, par7)
   integer, intent(in) :: FN	! file number to write to
   real(8), intent(in) :: par1, par2	! at least two parameters must be written into the file
   real(8), intent(in), optional ::  par3, par4, par5, par6, par7	! parameters to be written into the file
   if (present(par3) .AND. present(par4) .AND. present(par5) .AND. present(par6) .AND. present(par7)) then
      write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') par1, par2, par3, par4, par5, par6, par7
   elseif (present(par3) .AND. present(par4) .AND. present(par5) .AND. present(par6)) then
      write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') par1, par2, par3, par4, par5, par6
   elseif (present(par3) .AND. present(par4) .AND. present(par5)) then
      write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16)') par1, par2, par3, par4, par5
   elseif (present(par3) .AND. present(par4)) then
      write(FN, '(es25.16,es25.16,es25.16,es25.16)') par1, par2, par3, par4
   elseif (present(par3)) then
      write(FN, '(es25.16,es25.16,es25.16)') par1, par2, par3
   else
      write(FN, '(es25.16,es25.16)') par1, par2
   endif
end subroutine write_output_file_one


subroutine write_temperatures_n_displacements(FN, time, Te, Ta, Ta_sub, MSD, MSDP)
   integer, intent(in) :: FN	! file number to write to
   real(8), intent(in) :: time	! [fs]
   real(8), intent(in) :: Te, Ta, MSD   ! electron temperature, atomic temperature, atomic mean displacement
   real(8), dimension(:), intent(in), allocatable :: Ta_sub, MSDP  ! temperature and atomic mean displacement for each element
   if (.not.allocated(MSDP)) then
      write(FN, '(es25.16,$)') time, Te, Ta, Ta_sub(:), MSD
   else
      write(FN, '(es25.16,$)') time, Te, Ta, Ta_sub(:), MSD, MSDP(:)
   endif
   write(FN,'(a)')  ! make a new line
end subroutine write_temperatures_n_displacements



subroutine save_testmode_data(FN, time, Scell)  ! center-of-mass, rotation, total forces, etc.
   integer, intent(in) :: FN      ! file number to write to
   real(8), intent(in) :: time   ! [fs]
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   !-------------------------------
   write(FN, '(es25.16, $)') time, &   ! time
                             Scell%V_CoM, & ! center of mass velosity (3)
                             Scell%I_tot, & ! moment of inertia tensor (3x3)
                             Scell%F_tot    ! total force of all atoms (3)
   write(FN,'(a)')  ! make a new line
end subroutine save_testmode_data



subroutine write_sectional_displacements(FN_displacements, time, Scell, matter) ! atomic displaecements
   integer, dimension(:), intent(in) :: FN_displacements   ! file numbers to write to
   real(8), intent(in) :: time   ! [fs]
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! parameters of the material
   !-------------------------------
   integer :: Nsiz, i, j, N_at

   ! Atomic masks for sectional displacements:
   Nsiz = size(Scell%Displ)   ! how many masks
   do i = 1, Nsiz    ! for all masks
      write(FN_displacements(i), '(es25.16,$)') time, Scell%Displ(i)%mean_disp, Scell%Displ(i)%mean_disp_r(:)
      ! Now for kinds of atoms:
      N_at = matter%N_KAO    ! number of kinds of atoms
      if (allocated(Scell%Displ(i)%mean_disp_sort)) then
         do j = 1, N_at
            write(FN_displacements(i), '(es25.16,$)') Scell%Displ(i)%mean_disp_sort(j), &
            Scell%Displ(i)%mean_disp_r_sort(j,1), Scell%Displ(i)%mean_disp_r_sort(j,2), Scell%Displ(i)%mean_disp_r_sort(j,3)
         enddo
      else
         print*, 'mean_disp_sort is not allocated, cannot printout sectional_displacements'
      endif
      write(FN_displacements(i),'(a)') ! make a new line
   enddo ! i
end subroutine write_sectional_displacements




subroutine save_diffraction_peaks(FN, time, Scell)
   integer, intent(in) :: FN  ! file number to save to
   real(8), intent(in) :: time   ! [fs]
   type(Super_cell), intent(in):: Scell ! super-cell with all the atoms inside
   !--------------
   integer :: i

   write(FN, '(es)', advance = 'no') time
   do i = 1, size(Scell%diff_peaks%I_diff_peak)
      write(FN, '(es)', advance = 'no') Scell%diff_peaks%I_diff_peak(i)
   enddo
   write(FN, '(a)') ''  ! next line
end subroutine save_diffraction_peaks



subroutine save_diffraction_powder(FN, time, Scell)
   integer, intent(in) :: FN  ! file number to save to
   real(8), intent(in) :: time   ! [fs]
   type(Super_cell), intent(in):: Scell ! super-cell with all the atoms inside
   !--------------
   integer :: i

   write(FN,'(a,f25.16)') '#', time
   do i = 1, size(Scell%diff_peaks%two_theta)
      write(FN, '(f25.16,es)') Scell%diff_peaks%two_theta(i)*g_rad2deg, Scell%diff_peaks%I_powder(i)
   enddo
   write(FN, '(a)') ''
   write(FN, '(a)') ''
end subroutine save_diffraction_powder



subroutine prepare_output_files(Scell, matter, laser, numpar, TB_Hamil, TB_Repuls, Err)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(TB_repulsive), intent(in) :: TB_Repuls  ! parameters of the repulsive part of TB
   type(TB_Hamiltonian), intent(in) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(Error_handling), intent(inout) :: Err	! error save
   !========================================================
   character(200) :: File_name, Error_descript, chtest, chtest1, chtest_MDgrid, chtest_AtBathgrid, chtest_ElBathgrid
   character(3) :: chtest2
   integer INFO
   integer :: MOD_TIM ! time when the communication.txt file was last modified
   logical :: file_opened, file_exist, NP_file_exists, IM_file_exists


   !--------------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything wrong here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return
   endif
   !--------------------------------------------------------------------------


   ! Create directory where the output files will be saved:
   call create_output_folder(Scell, matter, laser, numpar)	! module "Dealing_with_output_files"

   ! Save input files, so that repeating the same calculations would be easy:
   if (numpar%which_input >= 1) then
      
      ! Check if new format of input file exists:
      !chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT'
      chtest = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_INPUT_MINIMUM))

      write(chtest2,'(i3)') numpar%which_input
      write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
      inquire(file=trim(adjustl(chtest)),exist=file_exist)
      
      if (.not.file_exist) then ! use old format of input files:
         !chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT_MATERIAL'
         chtest = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_INPUT_MATERIAL))
         write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
         inquire(file=trim(adjustl(chtest)),exist=IM_file_exists)
         if (.not.IM_file_exists) then ! check if short-named file exists
            chtest = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_INPUT_ALL))
            write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
         endif

         !chtest1 = 'INPUT_DATA'//numpar%path_sep//'NUMERICAL_PARAMETERS'
         chtest1 = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_NUMERICAL_PARAMETERS))
         write(chtest1,'(a,a,a,a)') trim(adjustl(chtest1)), '_', trim(adjustl(chtest2)), '.txt'
         inquire(file=trim(adjustl(chtest1)),exist=NP_file_exists)
         if (.not.NP_file_exists) then ! check if unnumbered file with NumPars:
            chtest1 = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_NUMERICAL_PARAMETERS))//'.txt'
            inquire(file=trim(adjustl(chtest1)),exist=NP_file_exists)
         endif

         ! And a file with MD_grid if user provided:
         if (allocated(numpar%dt_MD_reset_grid)) then
            chtest_MDgrid = 'INPUT_DATA'//numpar%path_sep//trim(adjustl(numpar%MD_step_grid_file))
         endif
         ! File with electronic thermostat parameters:
         if (allocated(numpar%El_bath_reset_grid)) then
            chtest_ElBathgrid = 'INPUT_DATA'//numpar%path_sep//trim(adjustl(numpar%El_bath_step_grid_file))
         endif
         ! File with atomic thermostat parameters:
         if (allocated(numpar%At_bath_reset_grid)) then
            chtest_AtBathgrid = 'INPUT_DATA'//numpar%path_sep//trim(adjustl(numpar%At_bath_step_grid_file))
         endif
      endif
   else
      ! Check if new format of input file exists:
      !chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT.txt'
      chtest = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_INPUT_MINIMUM))//'.txt'
      inquire(file=trim(adjustl(chtest)),exist=file_exist)
      
      if (.not.file_exist) then ! use old format of input files:
         !chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT_MATERIAL.txt'
         chtest = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_INPUT_MATERIAL))//'.txt'
         inquire(file=trim(adjustl(chtest)),exist=IM_file_exists)
         if (.not.IM_file_exists) then ! check if short-named file exists
            chtest = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_INPUT_ALL))//'.txt'
         endif

         !chtest1 = 'INPUT_DATA'//numpar%path_sep//'NUMERICAL_PARAMETERS.txt'
         chtest1 = trim(adjustl(m_INPUT_directory))//numpar%path_sep//trim(adjustl(m_NUMERICAL_PARAMETERS))//'.txt'
         inquire(file=trim(adjustl(chtest1)),exist=NP_file_exists)

         ! And a file with MD_grid if user provided:
         if (allocated(numpar%dt_MD_reset_grid)) then
            chtest_MDgrid = 'INPUT_DATA'//numpar%path_sep//trim(adjustl(numpar%MD_step_grid_file))
         endif
         ! File with electronic thermostat parameters:
         if (allocated(numpar%El_bath_reset_grid)) then
            chtest_ElBathgrid = 'INPUT_DATA'//numpar%path_sep//trim(adjustl(numpar%El_bath_step_grid_file))
         endif
         ! File with atomic thermostat parameters:
         if (allocated(numpar%At_bath_reset_grid)) then
            chtest_AtBathgrid = 'INPUT_DATA'//numpar%path_sep//trim(adjustl(numpar%At_bath_step_grid_file))
         endif
      endif
   endif

   if (numpar%path_sep .EQ. '\') then	! if it is Windows

      call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_files"

      if (.not.file_exist .and. NP_file_exists) call copy_file(trim(adjustl(chtest1)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      ! And file with MD grid, if user provided:
      if (allocated(numpar%dt_MD_reset_grid)) then
         call copy_file(trim(adjustl(chtest_MDgrid)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      endif
      if (allocated(numpar%El_bath_reset_grid)) then
         call copy_file(trim(adjustl(chtest_ElBathgrid)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      endif
      if (allocated(numpar%At_bath_reset_grid)) then
         call copy_file(trim(adjustl(chtest_AtBathgrid)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      endif
   else ! it is linux
      call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"

      if (.not.file_exist .and. NP_file_exists) call copy_file(trim(adjustl(chtest1)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      ! And file with MD grid, if user provided:
      if (allocated(numpar%dt_MD_reset_grid)) then
         call copy_file(trim(adjustl(chtest_MDgrid)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      endif
      if (allocated(numpar%El_bath_reset_grid)) then
         call copy_file(trim(adjustl(chtest_ElBathgrid)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      endif
      if (allocated(numpar%At_bath_reset_grid)) then
         call copy_file(trim(adjustl(chtest_AtBathgrid)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      endif
   endif

   ! Create the file with parameters (the same that's printed on the screen):
   call output_parameters_file(Scell,matter,laser,numpar,TB_Hamil,TB_Repuls,Err)	! and save the input data for output, module "Dealing_with_output_files"

   ! Prepare a file for communication with the user:
   !File_name = trim(adjustl(numpar%output_path))//numpar%path_sep//'Comunication.txt'
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep//trim(adjustl(m_Communication))
   numpar%Filename_communication = File_name ! save it to reuse later
   numpar%FN_communication = 110
   open(UNIT=numpar%FN_communication, FILE = trim(adjustl(File_name)), status = 'replace')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      INFO = 2
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      call Save_error_details(Err, INFO, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 9999
   else
      call get_file_stat(trim(adjustl(File_name)), Last_modification_time=numpar%MOD_TIME) ! get the time when it was last modified
   endif
   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      if (file_opened) close(numpar%FN_communication) ! for windows, close the file to let the user write into it
      ! Check if the file was modified since the last time:
      call get_file_stat(trim(adjustl(File_name)), Last_modification_time=MOD_TIM) ! module 'Dealing_with_files'

      if (MOD_TIM /= numpar%MOD_TIME) then ! open file again only if it was modified by the user
         numpar%MOD_TIME = MOD_TIM ! save new time of the last modification
!          print*, MOD_TIM, numpar%MOD_TIME
      endif
   else ! it is linux
   endif

   ! Prepare all the output files (create and write titles:)
   call create_output_files(Scell,matter,laser,numpar)      ! below

9999 continue
end subroutine prepare_output_files


subroutine make_save_files(path)
   character(len=*), intent(in) :: path
   character(200) file_name
   integer FN, FN2, FN3
   FN = 700	! number is fixed
   file_name = trim(adjustl(path))//'SAVE_atoms.dat'
   open(UNIT=FN, FILE = trim(adjustl(file_name)))
   FN2 = 701	! number is fixed
   file_name = trim(adjustl(path))//'SAVE_supercell.dat'
   open(UNIT=FN2, FILE = trim(adjustl(file_name)))
   FN3 = 702
   !file_name = trim(adjustl(path))//'SAVE_parameters.dat'
   file_name = trim(adjustl(path))//'SAVE_el_distribution.dat'
   open(UNIT=FN3, FILE = trim(adjustl(file_name)))
end subroutine make_save_files



subroutine update_save_files(time, atoms, matter, numpar, Scell)
   real(8), intent(in) :: time	! [fs]
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   type(Solid), intent(in) :: matter	! Material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   integer i

   ! SAVE_atoms.dat :
   rewind(700)	! overwrite the old state
   do i = 1,size(atoms)
      write(700, '(i3,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') atoms(i)%KOA, atoms(i)%S(:), atoms(i)%S0(:), atoms(i)%SV(:), atoms(i)%SV0(:)
   enddo

   ! SAVE_supercell.dat :
   rewind(701)	! overwrite the old state
   write(701,*) Scell%supce(:,:)
   write(701,'(a)') ''
   write(701,*) Scell%supce0(:,:)
   write(701,'(a)') ''
   write(701,*) Scell%Vsupce(:,:)
   write(701,'(a)') ''
   write(701,*) Scell%Vsupce0(:,:)
   write(701,'(a)') ''

   ! SAVE_el_distribution.dat :
   rewind(702)	! overwrite the old state
   !write(702,'(es25.16,es25.16,es25.16,es25.16,es25.16)') time, Scell%Te, Scell%mu, Scell%Ne_low, Scell%Ta
   write(702,'(a)') '# Electron distribution'
   do i = 1, size(Scell%fe)
      write(702,'(f25.16, f25.16)') Scell%Ei(i), Scell%fe(i)
   enddo
end subroutine update_save_files


subroutine close_save_files()
   close(700)
   close(701)
   close(702)
end subroutine close_save_files


subroutine close_output_files(Scell, numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !-------------------
   logical :: file_opened
   integer :: Nsiz, i

   close(numpar%FN_temperatures)
   close(numpar%FN_pressure)
   close(numpar%FN_electron_properties)
   close(numpar%FN_Ce)
   close(numpar%FN_energies)
   close(numpar%FN_supercell)
   close(numpar%FN_numbers)
   close(numpar%FN_orb_resolved)
   close(numpar%FN_deep_holes)
   if (numpar%do_kappa) close(numpar%FN_kappa)
   if (numpar%do_kappa_dyn) close(numpar%FN_kappa_dyn)
   if (numpar%save_raw) close(numpar%FN_atoms_S)
   if (numpar%do_drude) close(numpar%FN_optics)
   if (numpar%save_XYZ) close(numpar%FN_atoms_R)
   if (numpar%save_CIF) close(numpar%FN_cif)
   if (numpar%save_Ei)  close(numpar%FN_Ei)
   if (numpar%save_DOS)  close(numpar%FN_DOS)
   if (numpar%DOS_splitting == 1) close(numpar%FN_coupling)
   if (numpar%save_fa) then
      close(numpar%FN_fa)
      close(numpar%FN_fa_pot)
      close(numpar%FN_fa_tot)
   endif
   if (numpar%save_fe)  close(numpar%FN_fe)
   if (numpar%save_fe_grid)  close(numpar%FN_fe_on_grid)
   if (numpar%save_PCF) close(numpar%FN_PCF)
   if (Scell(1)%eps%all_w) close(numpar%FN_all_w)
   if (numpar%save_NN) close(numpar%FN_neighbors)
   if (allocated(numpar%NN_radii)) then
      Nsiz = size(numpar%FN_element_NN)   ! how many masks
      do i = 1, Nsiz
         close(numpar%FN_element_NN(i))
      enddo
   endif

   close(numpar%FN_Se)
   if (numpar%do_partial_thermal) then
      close(numpar%FN_Te)
      close(numpar%FN_mu)
   endif
   close(numpar%FN_Sa)
   if (numpar%print_Ta) then
      close(numpar%FN_Ta)
      close(numpar%FN_Ta_part)
   endif

   if ( (allocated(Scell(1)%Displ)) .and. (allocated(numpar%FN_displacements)) ) then
      Nsiz = size(Scell(1)%Displ)   ! how many masks
      do i = 1, Nsiz ! for all atomic masks
         inquire(unit=numpar%FN_displacements(i),opened=file_opened)
         if (file_opened) close(numpar%FN_displacements(i))
      enddo
   endif

   if (numpar%save_diff_peaks) then
      close(numpar%FN_diff_peaks)
      close(numpar%FN_diff_powder)
   endif
end subroutine close_output_files


subroutine create_output_files(Scell, matter, laser, numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   character(200) :: file_path, file_name
   integer :: FN, i, j, Nshl, Nsiz, N_at
   ! OUTPUT files, name and address:
   character(100) :: file_temperatures	! time [fs], Te [K], Ta [K]
   character(100) :: file_pressure	! time [fs], stress_tensore(3,3) [GPa], Pressure [GPa]
   character(100) :: file_energies	! energies [eV]
   character(100) :: file_atoms_R	! atomic coordinates and velocities
   character(100) :: file_atoms_S	! atomic coordinates and velocities
   character(100) :: file_atoms_cif	! atomic coordinates in cif-format (standard for constructing diffraction patterns)
   character(100) :: file_atomic_entropy	! atomic entropy
   character(100) :: file_atomic_temperatures ! atomic temperatures (varioous definitions)
   character(100) :: file_supercell	! supercell vectors
   character(100) :: file_electron_properties	! electron properties
   character(200) :: file_electron_heat_capacity	! band-resolved electron heat capacity
   character(200) :: file_electron_heat_conductivity  ! electron heat conductivity
   character(100) :: file_electron_entropy	! electron entropy
   character(100) :: file_electron_temperatures ! electron temperatures (for band-resolved calculations)
   character(100) :: file_electron_chempot ! electron chemical potentials (for band-resolved calculations)
   character(100) :: file_numbers	! total numbers of electrons and holes
   character(100) :: file_orb_resolved ! orbital-resolved electronic data
   character(100) :: file_deep_holes	! number of deep-shell holes in each shell
   character(100) :: file_Ei		! energy levels
   character(100) :: file_DOS	! DOS
   character(100) :: file_coupling  ! partial coupling parameter
   character(100) :: file_fe		! electron distribution (low-energy part)
   !character(100) :: file_fe_partial   ! band-resolved electron distributions (low-energy part)
   character(100) :: file_fe_on_grid   ! electron distribution (full: low- + high-energy)
   character(100) :: file_fa		! atomic distribution
   character(100) :: file_PCF		! pair correlation function
   character(100) :: file_optics	! optical coefficients
   character(100) :: file_all_w		! optical coeffs for all hw
   character(100) :: file_NN		! nearest neighbors
   character(100), dimension(:), allocatable :: file_sect_displ, file_sect_displ_short  ! sectional displacements
   character(100), dimension(:), allocatable :: file_element_NN, file_element_NN_short      ! element-specific nearest neighbors
   character(100) :: file_diff_peaks, file_diff_powder      ! selected diffraction peaks, powder diffraction
   character(100) :: file_testmode		! testmode file
   character(100) :: chtemp
   character(200) :: chtemp2
   character(11) :: chtemp11, text1, text2, text3

   call make_save_files(trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep)))

   file_path = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))

   file_temperatures = trim(adjustl(file_path))//'OUTPUT_temperatures.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_temperatures)))
   numpar%FN_temperatures = FN
   call create_file_header(numpar%FN_temperatures, '#Time	Te	Ta(:)	Displacement(:)')
   if (numpar%MSD_power > 1) then
      write(chtemp11,'(i2)') numpar%MSD_power
      call create_file_header(numpar%FN_temperatures, '#[fs]	[K]	[K](:)	[A^'//trim(adjustl(chtemp11))//'](:)')
   else
      call create_file_header(numpar%FN_temperatures, '#[fs]	[K]	[K](:)	[A](:)')
   endif
!    call create_file_header(numpar%FN_temperatures, '#Time	Te	Ta(kin)	Ta(conf)	Displacement')
!    call create_file_header(numpar%FN_temperatures, '#[fs]	[K]	[K]	[K]	[A]')

   ! Atomic masks for sectional displacements:
   if (allocated(Scell(1)%Displ)) then
      Nsiz = size(Scell(1)%Displ)   ! how many masks
      if (.not. allocated(numpar%FN_displacements)) then  ! allocate and open files
         allocate(numpar%FN_displacements(Nsiz))
      endif
      if (.not. allocated(file_sect_displ)) then  ! allocate and open files
         allocate(file_sect_displ(Nsiz))
         allocate(file_sect_displ_short(Nsiz))
      endif
      do i = 1, Nsiz ! for all masks
         file_sect_displ_short(i) = 'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(i)%mask_name))//'.dat'
         file_sect_displ(i) = trim(adjustl(file_path))//trim(adjustl(file_sect_displ_short(i)))
         open(NEWUNIT=FN, FILE = trim(adjustl(file_sect_displ(i))))
         numpar%FN_displacements(i) = FN

         N_at = matter%N_KAO    ! number of kinds of atoms
         chtemp = ''  ! to start with
         do j = 1, N_at
            chtemp = trim(adjustl(chtemp))//'   '//trim(adjustl(matter%Atoms(j)%Name))//':total   X  Y  Z'
         enddo
         call create_file_header(numpar%FN_displacements(i), '#Time  Total X  Y  Z  '//trim(adjustl(chtemp)) )
         if (INT(Scell(1)%Displ(i)%MSD_power) > 1) then
            write(chtemp11,'(i2)') INT(Scell(1)%Displ(i)%MSD_power)
            call create_file_header(numpar%FN_displacements(i), '#[fs]   [A^'//trim(adjustl(chtemp11))//'](:)')
         else
            call create_file_header(numpar%FN_displacements(i), '#[fs]   [A](:)')
         endif
      enddo
   endif


   if (numpar%save_diff_peaks) then ! Diffraction:
      ! selected diffraction peaks:
      file_diff_peaks = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_diff_peaks)))
      numpar%FN_diff_peaks = FN
      ! Create the header, containing all the peaks:
      chtemp2 = ''      ! to start with
      do i = 1, size(Scell(1)%diff_peaks%I_diff_peak)
         write(text1, '(i0)') Scell(1)%diff_peaks%ijk_diff_peak(1,i)
         write(text2, '(i0)') Scell(1)%diff_peaks%ijk_diff_peak(2,i)
         write(text3, '(i0)') Scell(1)%diff_peaks%ijk_diff_peak(3,i)
         chtemp2 = trim(adjustl(chtemp2))//'    ('//trim(adjustl(text1))//trim(adjustl(text2))//trim(adjustl(text3))//')'//' '
      enddo
      !print*, trim(adjustl(chtemp2))
      call create_file_header(numpar%FN_diff_peaks, '#Time  '//trim(adjustl(chtemp2)) )    ! below
      call create_file_header(numpar%FN_diff_peaks, '#[fs]    [arb.units]')

      ! Powder diffraction:
      file_diff_powder = trim(adjustl(file_path))//'OUTPUT_diffraction_powder.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_diff_powder)))
      numpar%FN_diff_powder = FN
   endif


   if (numpar%save_testmode) then   ! testmode additional data (center of mass, rotation, total force, etc.)
      file_testmode = trim(adjustl(file_path))//'OUTPUT_testmode_data.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_testmode)))
      numpar%FN_testmode = FN
      call create_file_header(numpar%FN_testmode, '#Time	V_CoM I_tot F_tot')
      call create_file_header(numpar%FN_testmode, '#[fs]	[A/fs:3]	[3x3]	[N:3]')
   endif

   
   file_pressure = trim(adjustl(file_path))//'OUTPUT_pressure_and_stress.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_pressure)))
   numpar%FN_pressure = FN
   call create_file_header(numpar%FN_pressure, '#Time	Pressure	Pressure(x,x)	Pressure(x,y)	Pressure(x,z)	Pressure(y,x)	Pressure(y,y)	Pressure(y,z)	Pressure(z,x)	Pressure(z,y)	Pressure(z,z)')
   call create_file_header(numpar%FN_pressure, '#[fs]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]')

   file_electron_properties = trim(adjustl(file_path))//'OUTPUT_electron_properties.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_properties)))
   numpar%FN_electron_properties = FN
   call create_file_header(numpar%FN_electron_properties, &
      '#Time	Ne	mu	band_gap	Ce	Coupling_parameter	VB_bottom	VB_top	CB_bottom	CB_top Mullikens(:)')
   call create_file_header(numpar%FN_electron_properties, &
      '#[fs]	[1/atom]	[eV]	[eV]	[J/(m^3K)]	[W/(m^3K)]	[eV]	[eV]	[eV]	[eV]  [e](:)')

   file_electron_heat_capacity = trim(adjustl(file_path))//'OUTPUT_electron_Ce.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_heat_capacity)))
   numpar%FN_Ce = FN
   call write_Ce_header(numpar%FN_Ce, Scell, 1, matter) ! below

   file_electron_entropy = trim(adjustl(file_path))//'OUTPUT_electron_entropy.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_entropy)))
   numpar%FN_Se = FN
   if (numpar%do_partial_thermal) then
      call create_file_header(numpar%FN_Se, '#Time Se  Se_eq   Se_VB Se_eq_VB Se_CB Se_eq_CB')
      call create_file_header(numpar%FN_Se, '#[fs]  [eV/K]   [eV/K]   [eV/K]   [eV/K]   [eV/K]   [eV/K]')
   else
      call create_file_header(numpar%FN_Se, '#Time Se  Se_eq')
      call create_file_header(numpar%FN_Se, '#[fs]  [eV/K]   [eV/K]')
   endif

   if (numpar%print_Ta) then
      file_atomic_entropy = trim(adjustl(file_path))//'OUTPUT_atomic_entropy.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_atomic_entropy)))
      numpar%FN_Sa = FN
      call create_file_header(numpar%FN_Sa, '#Time Sa  Sa_eq   Sa_eq_num   Sa_conf  Sa_tot')
      call create_file_header(numpar%FN_Sa, '#[fs]  [eV/K]   [eV/K]  [eV/K]   [eV/K] [eV/K]')

      file_atomic_temperatures = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_atomic_temperatures)))
      numpar%FN_Ta = FN
      call create_file_header(numpar%FN_Ta, '#Time kin   entropic distr fluct  pot   virial  sin^2(1)  config  F*v   conf  hyperconf')
      call create_file_header(numpar%FN_Ta, '#[fs] [K]   [K]  [K]   [K]   [K]   [K]   [K]   [K] [-]   [K]   [K]')

      file_atomic_temperatures = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_partial.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_atomic_temperatures)))
      numpar%FN_Ta_part = FN
      call create_file_header(numpar%FN_Ta_part, '#Time kin:X   kin:Y  kin:Z vir:X   vir:Y   vir:Z')
      call create_file_header(numpar%FN_Ta_part, '#[fs]  [K]   [K]  [K]   [K]   [K]   [K]')
   endif

   if (numpar%do_partial_thermal) then
      file_electron_temperatures = trim(adjustl(file_path))//'OUTPUT_electron_temperatures.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_temperatures)))
      numpar%FN_Te = FN
      call create_file_header(numpar%FN_Te, '#Time Te  Te_VB Te_CB')
      call create_file_header(numpar%FN_Te, '#[fs]  [K]   [K]   [K]')

      file_electron_chempot = trim(adjustl(file_path))//'OUTPUT_electron_chempotentials.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_chempot)))
      numpar%FN_mu = FN
      call create_file_header(numpar%FN_mu, '#Time mu  mu_VB mu_CB')
      call create_file_header(numpar%FN_mu, '#[fs]  [eV]   [eV]   [eV]')
!    else
!       call create_file_header(numpar%FN_Te, '#Time Te')
!       call create_file_header(numpar%FN_Te, '#[fs]  [K]')
   endif

   if (numpar%do_kappa) then
      file_electron_heat_conductivity = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_heat_conductivity)))
      numpar%FN_kappa = FN
      !call create_file_header(numpar%FN_kappa, '#Time kappa')
      !call create_file_header(numpar%FN_kappa, '#[fs]  [W/(K*m)]')
      call create_file_header(numpar%FN_kappa, '#Te   kappa_tot   kappa_e_ph  kappa_e_e mu Ce')
      call create_file_header(numpar%FN_kappa, '#[K]  [W/(K*m)]   [W/(K*m)]   [W/(K*m)]   [eV]  [J/(m^3*K)]')
   endif

   if (numpar%do_kappa_dyn) then
      file_electron_heat_conductivity = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_dyn.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_heat_conductivity)))
      numpar%FN_kappa_dyn = FN
      call create_file_header(numpar%FN_kappa_dyn, '#time   kappa_tot   kappa_e_ph  kappa_e_e')
      call create_file_header(numpar%FN_kappa_dyn, '#[fs]  [W/(K*m)]   [W/(K*m)]   [W/(K*m)]')
   endif

   file_energies = trim(adjustl(file_path))//'OUTPUT_energies.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_energies)))
   numpar%FN_energies = FN
   call create_file_header(numpar%FN_energies, &
   '#Time	Electrons	Holes	Potential	Kinetic	Atoms	Atoms_n_electrons	Atom_all_electrons	Total	van_der_Waals   Short-range')
   call create_file_header(numpar%FN_energies, &
   '#[fs]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]  [eV/atom]')

   file_numbers = trim(adjustl(file_path))//'OUTPUT_electron_hole_numbers.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_numbers)))
   numpar%FN_numbers = FN
   write(numpar%FN_numbers, '(a)') '#Time	VB_electrons	CB_electrons	High_energy_electrons	Deep_holes	Error	Photons'
   write(numpar%FN_numbers, '(a)') '#[fs]	[1/atom]	[1/atom]	[1/atom]	[1/atom]	[1/atom]	[1/atom]'


   file_orb_resolved = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_data.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_orb_resolved)))
   numpar%FN_orb_resolved = FN
   call write_orb_resolved_header(numpar%FN_orb_resolved, Scell(1), matter) ! below


   file_deep_holes = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_deep_holes)))
   numpar%FN_deep_holes = FN
   write(numpar%FN_deep_holes, '(a)', advance='no') '#Time	'
   ATOMS:do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      SHELLS:do j = 1, Nshl ! for all shells of this atom
            if ((i .NE. 1) .or. (j .NE. Nshl)) then ! atomic shell:
               call define_PQN(matter%Atoms(i)%Shl_dsgnr(j), chtemp11) ! module "Dealing_with_EADL"
               write(chtemp,'(a)') trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(chtemp11))
               write(numpar%FN_deep_holes, '(a)', advance='no') trim(adjustl(chtemp))//'	'
            else ! VB:
               ! skip it
            endif
      enddo SHELLS
   enddo ATOMS
   write(numpar%FN_deep_holes, '(a)', advance='yes') ' '


   if (numpar%DOS_splitting == 1) then
      file_coupling = trim(adjustl(file_path))//'OUTPUT_coupling.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_coupling)))
      numpar%FN_coupling = FN
      call write_coulping_header(numpar%FN_coupling, Scell, 1, matter, numpar) ! below
   endif


   if (numpar%save_raw) then
      file_atoms_S = trim(adjustl(file_path))//'OUTPUT_coordinates_and_velocities.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_atoms_S)))
      numpar%FN_atoms_S = FN
!       write(numpar%FN_atoms_S, '(a)') 'X	Y	Z	Vx	Vy	Vz	Sx	Sy	Sz	SVx	SVy	SVz'
!       write(numpar%FN_atoms_S, '(a)') ''
   endif

   file_supercell = trim(adjustl(file_path))//'OUTPUT_supercell.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_supercell)))
   numpar%FN_supercell = FN
   write(numpar%FN_supercell,'(a)') '#Time	Volume	h11	h12	h13	h21	h22	h23	h31	h32	h33	Vh11	Vh12	Vh13	Vh21	Vh22	Vh23	Vh31	Vh32	Vh33'
   write(numpar%FN_supercell,'(a)') '#[fs]	[A^3]	[A]	[A]	[A]	[A]	[A]	[A]	[A]	[A]	[A]	[A/fs]	[A/fs]	[A/fs]	[A/fs]	[A/fs]	[A/fs]	[A/fs]	[A/fs]	[A/fs]'


   if (numpar%do_drude) then
      file_optics = trim(adjustl(file_path))//'OUTPUT_optical_coefficients.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_optics)))
      numpar%FN_optics = FN
      write(numpar%FN_optics, '(a)') '#Time	Reflectivity	Transmission	Absorption	n	k	Re(e)	Im(e)	dc-conductivity	Re_Exx	ImE_xx	ReE_yy	ImE_yy	ReE_zz	ImE_zz'
   endif

   if (numpar%save_XYZ) then
      file_atoms_R = trim(adjustl(file_path))//'OUTPUT_atomic_coordinates.xyz'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_atoms_R)))
      numpar%FN_atoms_R = FN
   endif
   
   if (numpar%save_CIF) then
      file_atoms_cif = trim(adjustl(file_path))//'OUTPUT_atomic_coordinates.cif'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_atoms_cif)))
      numpar%FN_cif = FN
   endif

   if (numpar%save_Ei) then
      file_Ei = trim(adjustl(file_path))//'OUTPUT_energy_levels.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_Ei)))
      numpar%FN_Ei = FN
   endif
   
   if (numpar%save_DOS) then
      file_DOS = trim(adjustl(file_path))//'OUTPUT_DOS.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_DOS)))
      numpar%FN_DOS = FN
   endif
   
   if (numpar%save_fe) then
      file_fe = trim(adjustl(file_path))//'OUTPUT_electron_distribution.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_fe)))
      numpar%FN_fe = FN
   endif

   if (numpar%save_fe_grid) then
      file_fe_on_grid = trim(adjustl(file_path))//'OUTPUT_electron_distribution_on_grid.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_fe_on_grid)))
      numpar%FN_fe_on_grid = FN
   endif

   if (numpar%save_fa) then
      ! kinetic energy distribution:
      file_fa = trim(adjustl(file_path))//'OUTPUT_atomic_distribution.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_fa)))
      numpar%FN_fa = FN
      ! potential energy distribution:
      file_fa = trim(adjustl(file_path))//'OUTPUT_atomic_distribution_pot.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_fa)))
      numpar%FN_fa_pot = FN
      ! and total energy distribution:
      file_fa = trim(adjustl(file_path))//'OUTPUT_atomic_distribution_tot.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_fa)))
      numpar%FN_fa_tot = FN
   endif

   if (numpar%save_PCF) then
      file_PCF = trim(adjustl(file_path))//'OUTPUT_pair_correlation_function.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_PCF)))
      numpar%FN_PCF = FN
   endif
   
   if (numpar%save_NN) then
      file_NN = trim(adjustl(file_path))//'OUTPUT_nearest_neighbors.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_NN)))
      numpar%FN_neighbors = FN
      write(numpar%FN_neighbors, '(a)') '#Time Average N(0)    N(1)    N(2)    N(3)    N(4)    N(5)    N(6)'
   endif

   if (allocated(numpar%NN_radii)) then
      if (.not.allocated(file_element_NN)) allocate(file_element_NN(size(numpar%NN_radii)))
      if (.not.allocated(file_element_NN_short)) allocate(file_element_NN_short(size(numpar%NN_radii)))
      if (.not.allocated(numpar%FN_element_NN)) allocate(numpar%FN_element_NN(size(numpar%NN_radii)))

      do i = 1, size(numpar%NN_radii) ! for all requested elements

         file_element_NN(i) = 'OUTPUT_nearest_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'.dat'
         file_element_NN_short(i) = file_element_NN(i)      ! save short name for gnuplotting
         file_element_NN(i) = trim(adjustl(file_path))//trim(adjustl(file_element_NN(i)))

         open(NEWUNIT=FN, FILE = trim(adjustl(file_element_NN(i))))
         numpar%FN_element_NN(i) = FN
         write(numpar%FN_element_NN(i), '(a)', advance='no') '#Time     Total '
         do j = 1, size(numpar%NN_radii)     ! for all elements in the compound:
            write(numpar%FN_element_NN(i), '(a)', advance='no') trim(adjustl(matter%Atoms(j)%Name))//'  '
         enddo ! j
         write(numpar%FN_element_NN(i), '(a)') '' ! end line
      enddo ! i
   endif

   do i = 1, size(Scell)
      if (Scell(i)%eps%all_w) then
         file_all_w = trim(adjustl(file_path))//'OUTPUT_dielectric_function.dat'
         open(NEWUNIT=FN, FILE = trim(adjustl(file_all_w)))
         numpar%FN_all_w = FN
      endif
   enddo

   ! Prepare gnuplot scripts to plot all the output data:
   call create_gnuplot_scripts(Scell, matter, numpar, laser, file_path, &
   'OUTPUT_temperatures.dat', &
   'OUTPUT_pressure_and_stress.dat', &
   'OUTPUT_energies.dat', &
   file_atoms_R, file_atoms_S, &
   'OUTPUT_supercell.dat', &
   'OUTPUT_electron_properties.dat', &
   'OUTPUT_electron_heat_conductivity.dat', &
   'OUTPUT_electron_heat_conductivity_dyn.dat', &
   'OUTPUT_electron_hole_numbers.dat', &
   'OUTPUT_orbital_resolved_data.dat', &
   'OUTPUT_deep_shell_holes.dat', &
   'OUTPUT_optical_coefficients.dat', &
   file_Ei, file_PCF, &
   'OUTPUT_nearest_neighbors.dat', &
   file_element_NN_short, &
   'OUTPUT_electron_entropy.dat', &
   'OUTPUT_electron_temperatures.dat', &
   'OUTPUT_electron_chempotentials.dat', &
   'OUTPUT_atomic_entropy.dat', &
   'OUTPUT_atomic_temperatures.dat', &
   'OUTPUT_atomic_temperatures_partial.dat', &
   file_sect_displ_short, &
   'OUTPUT_diffraction_peaks.dat', &
   'OUTPUT_diffraction_powder.dat', &
   'OUTPUT_testmode_data.dat')  ! below

   ! clean up:
   if (allocated(file_sect_displ)) deallocate(file_sect_displ, file_sect_displ_short)
   if (allocated(file_element_NN)) deallocate(file_element_NN)
end subroutine create_output_files


subroutine create_file_header(FN, text)
   integer, intent(in) :: FN		! file number
   character(*), intent(in) :: text	! what to write in this file
   write(FN,'(a)') trim(adjustl(text))
end subroutine create_file_header


subroutine create_gnuplot_scripts(Scell,matter,numpar,laser, file_path, file_temperatures, file_pressure, file_energies, &
file_atoms_R, file_atoms_S, file_supercell, file_electron_properties, file_heat_capacity, file_heat_capacity_dyn, &
file_numbers, file_orb, file_deep_holes, file_optics, file_Ei, file_PCF, file_NN, file_element_NN, file_electron_entropy, file_Te, file_mu, &
file_atomic_entropy, file_atomic_temperatures, file_atomic_temperatures_part, file_sect_displ, &
file_diffraction_peaks, file_diffraction_powder, file_testmode)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   character(*), intent(in) :: file_path
   character(*), intent(in) :: file_temperatures	! time [fs], Te [K], Ta [K]
   character(*), intent(in) :: file_pressure	! pressure and stress tensore
   character(*), intent(in) :: file_energies	! energies [eV]
   character(*), intent(in) :: file_atoms_R	! atomic coordinates and velocities
   character(*), intent(in) :: file_atoms_S	! atomic coordinates and velocities
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
   character(*), dimension(:), intent(in) :: file_element_NN      ! element-specific nearest neighbors
   character(*), intent(in) :: file_electron_entropy  ! electron entropy
   character(*), intent(in) :: file_Te ! electron temperatures
   character(*), intent(in) :: file_mu ! electron chem.potentials
   character(*), intent(in) :: file_atomic_entropy ! atomic entropy
   character(*), intent(in) :: file_atomic_temperatures ! atomic temperatures (various definitions)
   character(*), intent(in) :: file_atomic_temperatures_part  ! partial atomic temperatures (X, Y, Z)
   character(*), dimension(:), intent(in) :: file_sect_displ
   character(*), intent(in) :: file_diffraction_peaks, file_diffraction_powder  ! selected diffraction peaks; powder diffraction
   character(*), intent(in) :: file_testmode	! testmode data
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

      call write_energy_levels_gnuplot(FN, Scell, 'OUTPUT_energy_levels.dat')
      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)
   endif

   
   ! Energies:
   File_name  = trim(adjustl(file_path))//'OUTPUT_energies_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_energies(numpar, File_name, file_energies, t0, t_last, 'OUTPUT_energies.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Temepratures:
   File_name  = trim(adjustl(file_path))//'OUTPUT_temperatures_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_temperatures(numpar, File_name, file_temperatures, t0, t_last, 'OUTPUT_temepratures.'//trim(adjustl(numpar%fig_extention))) ! below
   
   ! Mean square displacement:
   File_name  = trim(adjustl(file_path))//'OUTPUT_mean_displacement_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_MSD(File_name, file_temperatures, t0, t_last, 'OUTPUT_mean_displacement.'//trim(adjustl(numpar%fig_extention)), &
                numpar%MSD_power) ! below


   ! Atomic masks for sectional displacements:
   if (allocated(Scell(1)%Displ)) then
      Nsiz = size(Scell(1)%Displ)   ! how many masks
      do j = 1, Nsiz    ! for all masks
         File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(j)%mask_name))// &
                  '_Gnuplot'//trim(adjustl(sh_cmd))
         call gnu_displacements(File_name, file_sect_displ(j), t0, t_last, 'OUTPUT_mean_displacement_'// &
               trim(adjustl(Scell(1)%Displ(j)%mask_name))//'.'//trim(adjustl(numpar%fig_extention)), &
               Scell(1)%Displ(j)%MSD_power) ! below
         ! Partial by elements, if there is more than one:
         File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(j)%mask_name))// &
                  '_partial_Gnuplot'//trim(adjustl(sh_cmd))
         if (matter%N_KAO > 1) then
            call gnu_displacements_partial(File_name, file_sect_displ(j), t0, t_last, 'OUTPUT_mean_displacement_'// &
               trim(adjustl(Scell(1)%Displ(j)%mask_name))//'_partial.'//trim(adjustl(numpar%fig_extention)), &
               Scell(1)%Displ(j)%MSD_power, matter) ! below
         endif
      enddo ! j
   endif


   ! Diffraction:
   if (numpar%save_diff_peaks) then
      ! Diffraction peaks:
      File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_diffraction_peaks(Scell(1), File_name, file_diffraction_peaks, t0, t_last, &
                                    'OUTPUT_diffraction_peaks.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Powder diffraction:
      File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_powder_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 10.0d0, 'Powder', '2theta (deg)', 'Peak Intensity (a.u.)', 'OUTPUT_diffraction_powder.gif', numpar%path_sep, setkey=0)

      call write_diffraction_powder_gnuplot(FN, Scell, numpar, trim(adjustl(file_diffraction_powder)))      ! below

      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)
   endif
   
   ! Pressure:
   File_name  = trim(adjustl(file_path))//'OUTPUT_pressure_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_pressure(File_name, file_pressure, t0, t_last, 'OUTPUT_pressure.'//trim(adjustl(numpar%fig_extention))) ! below
   
   ! Stress tensor:
   File_name  = trim(adjustl(file_path))//'OUTPUT_stress_tensor_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_stress(File_name, file_pressure, t0, t_last, 'OUTPUT_pressure_tensor.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Numbers of particles:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_and_holes_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_numbers(File_name, file_numbers, t0, t_last, 'OUTPUT_electrons_holes.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Orbital-resolved electron parameters:
   File_name  = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_orbital_resolved(Scell(1), matter, numpar, File_name, file_orb, t0, t_last, 'OUTPUT_orbital_resolved_Ne.'// &
               trim(adjustl(numpar%fig_extention))) ! below

   ! Numbers of CB electrons:
   File_name  = trim(adjustl(file_path))//'OUTPUT_CB_electrons_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_CB_electrons(File_name, file_numbers, t0, t_last, 'OUTPUT_CB_electrons.'//trim(adjustl(numpar%fig_extention)))

   ! Numbers of deep-shell holes:
   File_name  = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_holes(File_name, file_deep_holes, t0, t_last, matter, 'OUTPUT_deep_shell_holes.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Band gap:
   File_name  = trim(adjustl(file_path))//'OUTPUT_Egap_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_Egap(File_name, file_electron_properties, t0, t_last, 'OUTPUT_Egap.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Chemical potential and Ne:
   File_name  = trim(adjustl(file_path))//'OUTPUT_mu_and_Ne_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_mu(File_name, file_electron_properties, t0, t_last, 'OUTPUT_mu_and_Ne.'//trim(adjustl(numpar%fig_extention))) ! below
   
   ! Boundaries of the bands:
   File_name  = trim(adjustl(file_path))//'OUTPUT_bands_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_Ebands(File_name, file_electron_properties, t0, t_last, 'OUTPUT_bands.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron heat capacity:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electron_Ce'//trim(adjustl(sh_cmd))
   call gnu_capacity(File_name, file_electron_properties, t0, t_last, 'OUTPUT_electron_Ce.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron heat conductivity:
   if (numpar%do_kappa) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity'//trim(adjustl(sh_cmd))
      call gnu_heat_conductivity(File_name, file_heat_capacity, &
      numpar%kappa_Te_min, numpar%kappa_Te_max, &
      'OUTPUT_electron_heat_conductivity.'//trim(adjustl(numpar%fig_extention))) ! below
   endif
   if (numpar%do_kappa_dyn) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_dyn'//trim(adjustl(sh_cmd))
      call gnu_heat_conductivity_dyn(File_name, file_heat_capacity_dyn, &
      numpar%kappa_Te_min, numpar%kappa_Te_max, &
      'OUTPUT_electron_heat_conductivity_dyn.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Electron entropy:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electron_entropy'//trim(adjustl(sh_cmd))
   call gnu_entropy(File_name, file_electron_entropy, t0, t_last, 'OUTPUT_electron_entropy.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron temperatures and chemical potential (for band-resolved calculations):
   if (numpar%do_partial_thermal) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_temperatures'//trim(adjustl(sh_cmd))
      call gnu_el_temperatures(File_name, file_Te, t0, t_last, &
               'OUTPUT_electron_temperatures.'//trim(adjustl(numpar%fig_extention))) ! below

      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_chempotentials'//trim(adjustl(sh_cmd))
      call gnu_chempots(File_name, file_mu, t0, t_last, &
               'OUTPUT_electron_chempotentials.'//trim(adjustl(numpar%fig_extention))) ! below
   endif


   ! Atomic temperatures (various definitions):
   if (numpar%print_Ta) then
      ! Atomic entropy:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_entropy'//trim(adjustl(sh_cmd))
      call gnu_entropy_atomic(File_name, file_atomic_entropy, t0, t_last, 'OUTPUT_atomic_entropy.'//trim(adjustl(numpar%fig_extention))) ! below

      File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures'//trim(adjustl(sh_cmd))
      call gnu_at_temperatures(File_name, file_atomic_temperatures, t0, t_last, &
               'OUTPUT_atomic_temperatures.'//trim(adjustl(numpar%fig_extention))) ! below

      File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_partial'//trim(adjustl(sh_cmd))
      call gnu_at_temperatures_part(File_name, file_atomic_temperatures_part, t0, t_last, &
               'OUTPUT_atomic_temperatures_partial.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Testmode additional data:
   if (numpar%save_testmode) then   ! if we want to gnuplot testmode data...
      !File_name  = trim(adjustl(file_path))//'OUTPUT_center_of_mass'//trim(adjustl(sh_cmd))
      !call gnu_center_of_mass(File_name, file_testmode, t0, t_last, 'OUTPUT_center_of_mass.'//trim(adjustl(numpar%fig_extention))) ! below
   endif


   ! Electron-ion coupling parameter:
   File_name  = trim(adjustl(file_path))//'OUTPUT_coupling_parameter'//trim(adjustl(sh_cmd))
   call gnu_coupling(File_name, file_electron_properties, t0, t_last, 'OUTPUT_coupling.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Volume:
   File_name  = trim(adjustl(file_path))//'OUTPUT_volume_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_volume(File_name, file_supercell, t0, t_last, 'OUTPUT_volume.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Mulliken charges:
   if (numpar%Mulliken_model >= 1) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_Mulliken_charges_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_Mulliken_charges(File_name, file_electron_properties, t0, t_last, 'OUTPUT_Mulliken_charges.'//trim(adjustl(numpar%fig_extention))) ! below
   endif
   
   ! Nearest neighbors:
   if (numpar%save_NN) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_neighbors_Gnuplot'//trim(adjustl(sh_cmd))
      call gnu_nearest_neighbors(File_name, file_NN, t0, t_last, 'OUTPUT_nearest_neighbors.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   ! Element-specific nearest neighbors:
   if (allocated(numpar%NN_radii)) then
      do i = 1, size(numpar%NN_radii) ! for all requested elements
         File_name  = trim(adjustl(file_path))//'OUTPUT_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'_Gnuplot'//trim(adjustl(sh_cmd))
         call gnu_nearest_neighbors_elements(File_name, file_element_NN(i), trim(adjustl(numpar%NN_radii(i)%Name)), matter, t0, t_last, &
              'OUTPUT_nearest_neighbors_'//trim(adjustl(numpar%NN_radii(i)%Name))//'.'//trim(adjustl(numpar%fig_extention))) ! below
      enddo ! i
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
      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)
   endif

   ! Orbital-resoloved distribution function of electrons:
   if (numpar%save_fe_orb) then
      ! Distribution function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_fe_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 5.0d0, 'Distribution', 'Energy (eV)', 'Electron distribution (a.u.)', 'OUTPUT_orbital_resolved_fe.gif', numpar%path_sep, setkey=0)
      call write_orb_distribution_gnuplot(FN, Scell, numpar, matter, 'OUTPUT_electron_distribution.dat')   ! below
      call write_gnuplot_script_ending(FN, File_name, 1)
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
      call write_gnuplot_script_ending(FN, File_name, 1)
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
      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)

      ! 2) Distribution of potential energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atoms_distribution_pot_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, x_tics, 'Distribution', 'Energy (eV)', 'Atomic distribution (a.u.)', 'OUTPUT_atomic_distribution_pot.gif', numpar%path_sep, setkey=0)
      call write_atomic_distribution_gnuplot(FN, Scell, numpar, &
            'OUTPUT_atomic_distribution_pot.dat', its_pot=.true. )   ! below
      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)

      ! 3) Distribution of total energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_atoms_distribution_tot_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, x_tics, 'Distribution', 'Energy (eV)', 'Atomic distribution (a.u.)', 'OUTPUT_atomic_distribution_tot.gif', numpar%path_sep, setkey=0)
      call write_atomic_distribution_gnuplot(FN, Scell, numpar, &
            'OUTPUT_atomic_distribution_tot.dat', its_pot=.true., no_maxwell=.true.)   ! below
      call write_gnuplot_script_ending(FN, File_name, 1)
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
         call write_gnuplot_script_ending(FN, File_name, 1)
         close(FN)
      endselect
   endif

   ! Optical coefficients
   if (numpar%do_drude) then
      File_name  = trim(adjustl(file_path))//'OUTPUT_optical_coefficients'//trim(adjustl(sh_cmd))
      call gnu_optical_coefficients(File_name, file_optics, t0, t_last, 'OUTPUT_optical_coefficients.'//trim(adjustl(numpar%fig_extention))) ! below
      ! also n and k:
      File_name  = trim(adjustl(file_path))//'OUTPUT_optical_n_and_k'//trim(adjustl(sh_cmd))
      call gnu_n_and_k(File_name, file_optics, t0, t_last, 'OUTPUT_optical_n_and_k.'//trim(adjustl(numpar%fig_extention))) ! below
   endif

   !ccccccccccccccccccccccccccccc
   ! Create also convolved plots:
   CONV:if (Scell(1)%eps%tau > 0.0d0) then ! convolved files too:
      ! Energies:
      File_name  = trim(adjustl(file_path))//'OUTPUT_energies_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_energies(numpar, File_name, trim(adjustl(file_energies(1:len(trim(adjustl(file_energies)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_energies_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Temepratures:
      File_name  = trim(adjustl(file_path))//'OUTPUT_temperatures_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_temperatures(numpar, File_name, trim(adjustl(file_temperatures(1:len(trim(adjustl(file_temperatures)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_temepratures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      
      ! Mean displacement:
      File_name  = trim(adjustl(file_path))//'OUTPUT_mean_displacement_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_MSD(File_name, trim(adjustl(file_temperatures(1:len(trim(adjustl(file_temperatures)))-4)))//'_CONVOLVED.dat', t0, t_last, &
            'OUTPUT_mean_displacement_CONVOLVED.'//trim(adjustl(numpar%fig_extention)), numpar%MSD_power) ! below
      
      ! Atomic masks for sectional displacements:
      if (allocated(Scell(1)%Displ)) then
         Nsiz = size(Scell(1)%Displ)   ! how many masks
         do i = 1, Nsiz    ! for all masks
            File_name = trim(adjustl(file_path))//'OUTPUT_displacements_'//trim(adjustl(Scell(1)%Displ(i)%mask_name))// &
                  '_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
            call gnu_displacements(File_name, trim(adjustl(file_sect_displ(i)(1:len(trim(adjustl(file_sect_displ(i))))-4) )), &
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
                  //trim(adjustl(numpar%fig_extention)), Scell(1)%Displ(i)%MSD_power, matter) ! below
            endif
         enddo ! i
      endif

      ! Diffraction:
      if (numpar%save_diff_peaks) then
         ! Diffraction peaks:
         File_name  = trim(adjustl(file_path))//'OUTPUT_diffraction_peaks_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_diffraction_peaks(Scell(1), File_name, &
            trim(adjustl(file_diffraction_peaks(1:len(trim(adjustl(file_diffraction_peaks)))-4)))//'_CONVOLVED.dat' , &
            t0, t_last, 'OUTPUT_diffraction_peaks_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif



      ! Pressure:
      File_name  = trim(adjustl(file_path))//'OUTPUT_pressure_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_pressure(File_name, trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_pressure_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Stress tensor:
      File_name  = trim(adjustl(file_path))//'OUTPUT_stress_tensor_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_stress(File_name, trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_pressure_tensor_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      
      ! Numbers of particles:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_and_holes_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_numbers(File_name, trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electrons_holes_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Orbital-resolved electron parameters:
      File_name  = trim(adjustl(file_path))//'OUTPUT_orbital_resolved_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_orbital_resolved(Scell(1), matter, numpar, File_name, &
               trim(adjustl(file_orb(1:len(trim(adjustl(file_orb)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_orbital_resolved_Ne_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Numbers of CB electrons:
      File_name  = trim(adjustl(file_path))//'OUTPUT_CB_electrons_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_CB_electrons(File_name, trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_CB_electrons_CONVOLVED.'//trim(adjustl(numpar%fig_extention)))

      ! Numbers of deep-shell holes:
      File_name  = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_holes(File_name, trim(adjustl(file_deep_holes(1:len(trim(adjustl(file_deep_holes)))-4)))//'_CONVOLVED.dat', t0, t_last, matter, 'OUTPUT_deep_shell_holes_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Band gap:
      File_name  = trim(adjustl(file_path))//'OUTPUT_Egap_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_Egap(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_Egap_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Chemical potential and Ne:
      File_name  = trim(adjustl(file_path))//'OUTPUT_mu_and_Ne_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_mu(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_mu_and_Ne_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      
      ! Boundaries of the bands:
      File_name  = trim(adjustl(file_path))//'OUTPUT_bands_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_Ebands(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_bands_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron heat capacity:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_Ce_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_capacity(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electron_Ce_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron heat conductivity:
      if (numpar%do_kappa) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_heat_conductivity(File_name, &
         trim(adjustl(file_heat_capacity(1:len(trim(adjustl(file_heat_capacity)))-4)))//'_CONVOLVED.dat', &
         t0, t_last, 'OUTPUT_electron_heat_conductivity.'//trim(adjustl(numpar%fig_extention))) ! below
      endif
      if (numpar%do_kappa_dyn) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity_dyn_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_heat_conductivity_dyn(File_name, &
         trim(adjustl(file_heat_capacity_dyn(1:len(trim(adjustl(file_heat_capacity_dyn)))-4)))//'_CONVOLVED.dat', &
         t0, t_last, 'OUTPUT_electron_heat_conductivity_dyn.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Electron entropy:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_entropy(File_name, trim(adjustl(file_electron_entropy(1:len(trim(adjustl(file_electron_entropy)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electron_entropy_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron temperatures and chemical potential (for band-resolved calculations):
      if (numpar%do_partial_thermal) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_temperatures_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_el_temperatures(File_name, trim(adjustl(file_Te(1:len(trim(adjustl(file_Te)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_electron_temperatures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_electron_chempotentials_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_chempots(File_name, trim(adjustl(file_mu(1:len(trim(adjustl(file_mu)))-4)))//'_CONVOLVED.dat', &
               t0, t_last, 'OUTPUT_electron_chempotentials_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Atomic temperatures (various definitions):
      if (numpar%print_Ta) then
          ! Atomic entropy:
         File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_entropy_atomic(File_name, trim(adjustl(file_atomic_entropy(1:len(trim(adjustl(file_atomic_entropy)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_atomic_entropy_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_at_temperatures(File_name, &
            trim(adjustl(file_atomic_temperatures(1:len(trim(adjustl(file_atomic_temperatures)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_atomic_temperatures_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

         File_name  = trim(adjustl(file_path))//'OUTPUT_atomic_temperatures_partial_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_at_temperatures_part(File_name, &
            trim(adjustl(file_atomic_temperatures_part(1:len(trim(adjustl(file_atomic_temperatures_part)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_atomic_temperatures_partial.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Electron-ion coupling parameter:
      File_name  = trim(adjustl(file_path))//'OUTPUT_coupling_parameter_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_coupling(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_coupling_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Volume:
      File_name  = trim(adjustl(file_path))//'OUTPUT_volume_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_volume(File_name, trim(adjustl(file_supercell(1:len(trim(adjustl(file_supercell)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_volume_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      
      ! Mulliken charges:
      if (numpar%Mulliken_model >= 1) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_Mulliken_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_Mulliken_charges(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_Mulliken_charges_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif

      ! Nearest neighbors:
      if (numpar%save_NN) then
         File_name  = trim(adjustl(file_path))//'OUTPUT_neighbors_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_nearest_neighbors(File_name, trim(adjustl(file_NN(1:len(trim(adjustl(file_NN)))-4)))//'_CONVOLVED.dat', &
            t0, t_last, 'OUTPUT_nearest_neighbors_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      endif
      
      if (numpar%do_drude) then
         ! optical coefficients:
         File_name  = trim(adjustl(file_path))//'OUTPUT_optical_coefficients_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_optical_coefficients(File_name, trim(adjustl(file_optics(1:len(trim(adjustl(file_optics)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_optical_coefficients_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
         ! also n and k:
         File_name  = trim(adjustl(file_path))//'OUTPUT_optical_n_and_k_CONVOLVED'//trim(adjustl(sh_cmd))
         call gnu_n_and_k(File_name, trim(adjustl(file_optics(1:len(trim(adjustl(file_optics)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_optical_n_and_k_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
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
   !call write_gnuplot_script_header(FN, 1, 3, 'Energies','Time (fs)', 'Energies (eV)', trim(adjustl(file_path))//'OUTPUT_energies.'//trim(adjustl(g_numpar%fig_extention)), 1)
   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Energies','Time (fs)', 'Energies (eV/atom)', trim(adjustl(eps_name)), 1)
   
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_energies


subroutine gnu_temperatures(numpar, File_name, file_temperatures, t0, t_last, eps_name)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
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
   
   if (g_matter%N_KAO == 1) then
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
         do i = 1, g_matter%N_KAO - 1
            write(chtemp,'(a)') g_matter%Atoms(i)%Name//' atoms'
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', 3+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         enddo
         write(chtemp,'(a)') g_matter%Atoms(g_matter%N_KAO)%Name//' atoms'
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', 3+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_temperatures)), '\"u 1:2 w l lw \"$LW\" title \"Electrons\" ,\'
         write(FN, '(a,a,a)') '\"', trim(adjustl(file_temperatures)), '\"u 1:3 w l lw \"$LW\" title \"Atoms average\" ,\'
         do i = 1, g_matter%N_KAO - 1
            write(chtemp,'(a)') g_matter%Atoms(i)%Name//' atoms'
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', 3+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) ,'\" ,\'
         enddo
         write(chtemp,'(a)') g_matter%Atoms(g_matter%N_KAO)%Name//' atoms'
         write(FN, '(a,i3,a,a)') '\"\" u 1:', 3+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//'\"'
      endif
   endif
   
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_temperatures


subroutine gnu_MSD(File_name, file_MSD, t0, t_last, eps_name, MSD_power)
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
      call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement','Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement','Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"
   endif
   
   if (g_matter%N_KAO == 1) then
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:5 w l lw LW title "Displacement" '
!          write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:5 w l lw LW title "Displacement" '	! if included Tconf
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:5 w l lw \"$LW\" title \"Displacement\" '
      endif
   else ! more than one element:
      i_start = 4 + g_matter%N_KAO
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] "' , trim(adjustl(file_MSD)), ' "u 1:', i_start ,' w l lw LW title "Average" ,\'
         do i = 1, g_matter%N_KAO - 1
            write(chtemp,'(a)') g_matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         enddo
         write(chtemp,'(a)') g_matter%Atoms(g_matter%N_KAO)%Name
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " '
      else
         write(FN, '(a,es25.16,a,a,a,i3,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_MSD)), '\"u 1:', i_start, ' w l lw \"$LW\" title \"Average\" ,\'
         do i = 1, g_matter%N_KAO - 1
            write(chtemp,'(a)') g_matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) ,'\" ,\'
         enddo
         write(chtemp,'(a)') g_matter%Atoms(g_matter%N_KAO)%Name
         write(FN, '(a,i3,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//'\"'
      endif
   endif ! (g_matter%N_KAO == 1)
   
   
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_MSD



subroutine gnu_displacements(File_name, file_MSD, t0, t_last, eps_name, MSD_power)
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

      call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), g_numpar%path_sep, 2)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), g_numpar%path_sep, 2)	! module "Gnuplotting"
   endif

   i_start = 2
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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

   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_displacements



subroutine gnu_displacements_partial(File_name, file_MSD, t0, t_last, eps_name, MSD_power, matter)
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_MSD	! input file
   real(8), intent(in) :: t0, t_last ! time instance [fs]
   character(*), intent(in) :: eps_name ! name of the figure
   real(8), intent(in) :: MSD_power ! power of MSD
   type(Solid), intent(in) :: matter     ! material parameters
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

      call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A^'//trim(adjustl(MSD_text))//')', trim(adjustl(eps_name)), g_numpar%path_sep, 2)	! module "Gnuplotting"
   else
      call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mean displacement', 'Time (fs)', &
        'Mean displacement (A)', trim(adjustl(eps_name)), g_numpar%path_sep, 2)	! module "Gnuplotting"
   endif

   i_start = 6
   do i = 1, matter%N_KAO
      chtemp = trim(adjustl(matter%Atoms(i)%Name))
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_displacements_partial



subroutine gnu_diffraction_peaks(Scell, File_name, file_diffraction_peaks, t0, t_last, fig_name)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   character(*), intent(in) :: File_name   ! file to create
   character(*), intent(in) :: file_diffraction_peaks ! input file
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Diffraction peak', &
            'Time (fs)', 'Peak intensity (arb. units)', trim(adjustl(fig_name)), g_numpar%path_sep, 0)      ! module "Gnuplotting"


   if (size(Scell%diff_peaks%I_diff_peak) == 1) then  ! only one peak to plot
      peak_name = make_diff_peak_name(Scell, 1) ! below
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_diffraction_peaks)), ' "u 1:2 w l lw LW title "'//trim(adjustl(peak_name))//'" '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_diffraction_peaks)), '\"u 1:2 w l lw \"$LW\" title \"'//trim(adjustl(peak_name))//'\" '
      endif
   else ! more than one peak:

      peak_name = make_diff_peak_name(Scell, 1) ! below
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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

   call write_gnuplot_script_ending(FN, File_name, 1)
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



subroutine gnu_Mulliken_charges(File_name, file_electron_properties, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Mulliken charge','Time (fs)', 'Mulliken charge', trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"
   
   if (g_matter%N_KAO == 1) then
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:11 w l lw LW title "Charge" '
      else
         write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:11 w l lw \"$LW\" title \"Charge\" '
      endif
   else ! more than one element:
      i_start = 10
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(chtemp,'(a)') g_matter%Atoms(1)%Name ! first element
         write(FN, '(a,es25.16,a,a,a,i3,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), '" u 1:', i_start+1 ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         do i = 2, g_matter%N_KAO - 1   ! intermediate elements
            write(chtemp,'(a)') g_matter%Atoms(i)%Name
            write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " ,\'
         enddo
         write(chtemp,'(a)') g_matter%Atoms(g_matter%N_KAO)%Name    ! last element
         write(FN, '(a,i3,a,a,a)') ' "" u 1:', i_start+i ,' w l lw LW title " ', trim(adjustl(chtemp))  ,' " '
      else
         write(chtemp,'(a)') g_matter%Atoms(1)%Name ! first element
         write(FN, '(a,es25.16,a,a,a,i3,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:', i_start+1, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) ,'\" ,\'
         do i = 2, g_matter%N_KAO - 1   ! intermediate elements
            write(chtemp,'(a)') g_matter%Atoms(i)%Name 
            write(FN, '(a,i3,a,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp)) , '\" ,\'
         enddo
         write(chtemp,'(a)') g_matter%Atoms(g_matter%N_KAO)%Name  ! last element
         write(FN, '(a,i3,a,a)') '\"\" u 1:', i_start+i, ' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))//'\"'
      endif
   endif ! (g_matter%N_KAO == 1)
   
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_Mulliken_charges


subroutine gnu_nearest_neighbors(File_name, file_NN, t0, t_last, fig_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Nearest neighbors','Time (fs)', 'Nearest neighbors fraction', trim(adjustl(fig_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"
   
   if (g_numpar%path_sep == '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_nearest_neighbors




subroutine gnu_nearest_neighbors_elements(File_name, file_NN, Name, matter, t0, t_last, fig_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Nearest neighbors', 'Time (fs)', &
      'Nearest neighbors of '//trim(adjustl(Name)), &
      trim(adjustl(fig_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"

   if (g_numpar%path_sep == '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_nearest_neighbors_elements


subroutine gnu_pressure(File_name, file_pressure, t0, t_last, eps_name)
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
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Pressure','Time (fs)', 'Pressure (GPa)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_pressure)), ' "u 1:2 w l lw LW title "Pressure" '
   else
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_pressure)), '\"u 1:2 w l lw \"$LW\" title \"Pressure\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_pressure


subroutine gnu_stress(File_name, file_pressure, t0, t_last, eps_name)
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
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Pressure tensor','Time (fs)', 'Pressure tensor (GPa)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_stress



subroutine gnu_numbers(File_name, file_numbers, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Numbers', 'Time (fs)', &
                        'Particles (1/atom)', trim(adjustl(eps_name)), g_numpar%path_sep, 1) ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_numbers



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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Numbers', 'Time (fs)', &
                        'Electrons (1/atom)', trim(adjustl(eps_name)), g_numpar%path_sep, 1) ! module "Gnuplotting"

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

   call write_gnuplot_script_ending(FN, trim(adjustl(File_name)), 1) ! below
   close(FN)
end subroutine gnu_orbital_resolved




subroutine gnu_CB_electrons(File_name, file_numbers, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'CB_electrons','Time (fs)', 'Electrons (per atom)', trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_numbers)), ' "u 1:($3/4*100) w l lw LW title "CB electrons" ,\'
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_numbers)), ' "u 1:($3) w l lw LW title "CB electrons" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_numbers)), ' "u 1:($7) w l lw LW title "Photons"'
   else
      !write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_numbers)), '\"u 1:(\$3/4*100) w l lw \"$LW\" title \"CB electrons\" ,\'
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_numbers)), '\"u 1:(\$3) w l lw \"$LW\" title \"CB electrons\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_numbers)), '\"u 1:(\$7) w l lw \"$LW\" title \"Photons\"'
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_CB_electrons


subroutine gnu_holes(File_name, file_deep_holes, t0, t_last, matter, eps_name)
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
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Holes','Time (fs)', 'Particles (total)', &
         trim(adjustl(eps_name)), g_numpar%path_sep, setkey=0, fontsize=font_siz)  ! module "Gnuplotting"
   
   counter = 0 ! to start with
   first_line = .true.  ! to start from the first line
   last_atom_no_shells = .false.   ! to start with
   W_vs_L:if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_holes


subroutine gnu_Egap(File_name, file_electron_properties, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Egap','Time (fs)', 'Bandgap (eV)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), '"u 1:4 w l lw LW title "Band gap" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:4 w l lw \"$LW\" title \"Band gap \" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_Egap



subroutine gnu_mu(File_name, file_electron_properties, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'mu and Ne', &
            'Time (fs)', 'Energy (eV)', trim(adjustl(eps_name)), g_numpar%path_sep, 0) ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_mu


subroutine gnu_Ebands(File_name, file_electron_properties, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'bands boundaries','Time (fs)', 'Energy (eV)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_Ebands



subroutine gnu_capacity(File_name, file_electron_properties, t0, t_last, eps_name)
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
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron Ce','Time (fs)', 'Heat capacity (J/(m^3 K))', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:5 w l lw LW title "Electron heat capacity"  '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:5 w l lw \"$LW\" title \"Electron heat capacity\"  '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_capacity


subroutine gnu_heat_conductivity_dyn(File_name, file_heat_capacity, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron Kappa','Time (fs)', 'Heat conductivity (W/(m K))', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_heat_capacity)), ' "u 1:2 w l lw LW title "kappa total",\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_heat_capacity)), '" u 1:3 w l lw LW title "kappa e-ph" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_heat_capacity)), '" u 1:4 w l lw LW title "kappa e-e" '

   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_heat_capacity)), '\"u 1:2 w l lw \"$LW\" title \"kappa total\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_heat_capacity)), '\" u 1:3 w l lw \"$LW\" title \"kappa e-ph\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_heat_capacity)), '\" u 1:4 w l lw \"$LW\" title \"kappa e-e\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_heat_conductivity_dyn


subroutine gnu_heat_conductivity(File_name, file_heat_capacity, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron K','Electron temperature (K)', &
            'Heat conductivity (W/(m^3 K))', trim(adjustl(eps_name)), g_numpar%path_sep, 0)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_heat_capacity)), &
               ' "u 1:2 w l lw LW title "Electron heat conductivity"  '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_heat_capacity)), &
               '\"u 1:2 w l lw \"$LW\" title \"Electron heat conductivity\"  '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_heat_conductivity


subroutine gnu_entropy(File_name, file_electron_entropy, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron entropy','Time (fs)', 'Electron entropy (eV/K)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_entropy)), '" u 1:3 w l lw LW title "Equilibrium" ,\'

      if (g_numpar%do_partial_thermal) then ! partial entropies too
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

      if (g_numpar%do_partial_thermal) then ! partial entropies too
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:2 w l lw \"$LW\" title \"Nonequilibrium\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:5 w l lw \"$LW\" title \"Equilibrium VB\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:4 w l lw \"$LW\" title \"Nonequilibrium VB\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:7 w l lw \"$LW\" title \"Equilibrium CB\" ,\'
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:6 w l lw \"$LW\" title \"Nonequilibrium CB\" '
      else  ! only total, no partial
         write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:2 w l lw \"$LW\" title \"Nonequilibrium\" '
      endif
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_entropy



subroutine gnu_entropy_atomic(File_name, file_atomic_entropy, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Atomic entropy', 'Time (fs)', &
         'Atomic entropy (eV/K)', trim(adjustl(eps_name)), g_numpar%path_sep, 1)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_entropy_atomic


subroutine gnu_el_temperatures(File_name, file_Te, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Electron tempereature', 'Time (fs)', 'Electron temperature (K)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_Te)), '" u 1:2 w l lw LW title "Total (kinetic)" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Te)), '" u 1:3 w l lw LW title "Valence" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_Te)), '" u 1:4 w l lw LW title "Conduction" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_Te)), '\" u 1:2 w l lw \"$LW\" title \"Total (kinetic)\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Te)), '\" u 1:3 w l lw \"$LW\" title \"Valence\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_Te)), '\" u 1:4 w l lw \"$LW\" title \"Conduction\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_el_temperatures



subroutine gnu_at_temperatures(File_name, file_Ta, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Atomic tempereature', &
         'Time (fs)', 'Atomic temperature (K)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_at_temperatures



subroutine gnu_at_temperatures_part(File_name, file_Ta, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Atomic tempereature', &
         'Time (fs)', 'Atomic temperature (K)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_at_temperatures_part



subroutine gnu_chempots(File_name, file_mu, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Chemical potential', 'Time (fs)', 'Electron chemical potential (eV)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)   ! module "Gnuplotting"

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_mu)), '" u 1:2 w l lw LW title "Total" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_mu)), '" u 1:3 w l lw LW title "Valence" ,\'
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_mu)), '" u 1:4 w l lw LW title "Conduction" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_mu)), '\" u 1:2 w l lw \"$LW\" title \"Total\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_mu)), '\" u 1:3 w l lw \"$LW\" title \"Valence\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_mu)), '\" u 1:4 w l lw \"$LW\" title \"Conduction\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_chempots



subroutine gnu_coupling(File_name, file_electron_properties, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Coupling parameter','Time (fs)', &
      'Coupling parameter (W/(m^3 K))', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), &
         ' "u 1:6 w l lw LW title "Electron-ion coupling" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), &
         '\"u 1:6 w l lw \"$LW\" title \"Electron-ion coupling\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_coupling


subroutine gnu_volume(File_name, file_supercell, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Volume','Time (fs)', 'Volume (A^3)', &
      trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_supercell)), ' "u 1:2 w l lw LW title "Supercell volume" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"', trim(adjustl(file_supercell)), '\"u 1:2 w l lw \"$LW\" title \"Supercell volume\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_volume


subroutine gnu_optical_coefficients(File_name, file_optics, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Optical coefficients','Time (fs)', &
      'Optical coefficients', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
      
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_optics)), ' "u 1:2 w l lw LW title "Reflectivity" ,\'
      write(FN, '(a,a,a,a,a)') ' "', trim(adjustl(file_optics)), ' "u 1:3 w l lw LW title "Transmission" ,\'
      write(FN, '(a,a,a,a,a)') ' "', trim(adjustl(file_optics)), ' "u 1:4 w l lw LW title "Absorption" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_optics)), '\"u 1:2 w l lw \"$LW\" title \"Reflectivity\" ,\'
      write(FN, '(a,a,a,a,a)') '\"', trim(adjustl(file_optics)), '\"u 1:3 w l lw \"$LW\" title \"Transmission\" ,\'
      write(FN, '(a,a,a,a,a)') '\"', trim(adjustl(file_optics)), '\"u 1:4 w l lw \"$LW\" title \"Absorption \" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_optical_coefficients


subroutine gnu_n_and_k(File_name, file_optics, t0, t_last, eps_name)
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

   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Optical n and k','Time (fs)', &
      'Optical parameters', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_optics)), ' "u 1:5 w l lw LW title "n" ,\'
      write(FN, '(a,a,a)') ' "', trim(adjustl(file_optics)), ' "u 1:6 w l lw LW title "k" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_optics)), '\"u 1:5 w l lw \"$LW\" title \"n\" ,\'
      write(FN, '(a,a,a)') '\"', trim(adjustl(file_optics)), '\"u 1:6 w l lw \"$LW\" title \"k \" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_n_and_k





subroutine write_gnuplot_script_ending(FN, File_name, ind)
   integer, intent(in) :: FN, ind
   character(*), intent(in) :: File_name
   character(300) :: command
   integer :: iret

   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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


subroutine execute_all_gnuplots(file_path)
   character(*), intent(in) :: file_path
   character(300) :: command
   integer :: iret
   
   !call chdir(trim(adjustl(file_path)))
   command = trim(adjustl(file_path))
   iret = chdir(command)
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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


subroutine write_energy_levels_gnuplot(FN, Scell, file_Ei)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   character(*), intent(in) :: file_Ei  ! file with energy levels
   integer i, M, NSC
   character(30) :: ch_temp

   do NSC = 1, size(Scell)
      M = size(Scell(NSC)%Ei)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') 25.0d0      ! Scell(NSC)%E_top

      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a,a,a,i5,a)') 'p [][:'//trim(adjustl(ch_temp))//'] "', trim(adjustl(file_Ei)), '"u 1:', 2, ' pt 7 ps 0.2 ,\'
         do i = 3, M
            write(FN, '(a,a,a,i5,a)') ' "', trim(adjustl(file_Ei)), '"u 1:', i, ' pt 7 ps 0.2 ,\'
         enddo
         write(FN, '(a,a,a,i5,a)') ' "', trim(adjustl(file_Ei)), '"u 1:', M+1, ' pt 7 ps 0.2'
      else
         write(FN, '(a,a,a,i5,a)') 'p [][:'//trim(adjustl(ch_temp))//'] \"', trim(adjustl(file_Ei)), &
                                       '\"u 1:', 2, ' pt 7 ps \"$LW\" ,\'
         do i = 3, M
            write(FN, '(a,a,a,i5,a)') '\"', trim(adjustl(file_Ei)), '\"u 1:', i, ' pt 7 ps \"$LW\" ,\'
         enddo
         write(FN, '(a,a,a,i5,a)') '\"', trim(adjustl(file_Ei)), '\"u 1:', M+1, ' pt 7 ps \"$LW\"'
      endif
   enddo
end subroutine write_energy_levels_gnuplot




subroutine write_diffraction_powder_gnuplot(FN, Scell, numpar, file_powder, min_x, max_x)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_powder  ! file with powder diffraction spectrum
   real(8), intent(in), optional :: min_x, max_x      ! start and end of x-grid
   !-----------------------
   real(8) :: x_start, x_end
   integer :: i, M, NSC
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4
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

      if (g_numpar%path_sep .EQ. '\') then      ! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_powder))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][:] "'//trim(adjustl(file_powder))// &
                  '" index (i-1) u 1:2 w l lw 2 lt rgb "black" title sprintf("%i fs",(i*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_powder))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][:] \"'//trim(adjustl(file_powder))// &
                  '\" index (i-1) u 1:2 w l lw 2 lt rgb \"black\" title sprintf(\"%i fs\",(i*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
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
         case (3:4)
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

      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         !write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:5] "'//trim(adjustl(file_fe))// &

         if (no_maxw) then ! just distribution, without equivalent Maxwell:
            write(FN, '(a)') 'p [:][0:5] "'//trim(adjustl(file_fe))// &
                  '" index (i) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         else ! with equivalent Maxwell
            write(FN, '(a)') 'p [:][0:5] "'//trim(adjustl(file_fe))// &
                  '" index (i) u 1:3 w l lw 2 lt rgb "grey" title "'// trim(adjustl(ch_titel)) //'" ,\'

            write(FN, '(a)') ' "'//trim(adjustl(file_fe))// &
                  '" index (i) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif ! no_maxw
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'

         !write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:5] \"'//trim(adjustl(file_fe))// &
         if (no_maxw) then ! just distribution, without equivalent Maxwell:
            write(FN, '(a)') 'p [:][0:5] \"'//trim(adjustl(file_fe))// &
            '\" index (i) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         else
            write(FN, '(a)') 'p [:][0:5] \"'//trim(adjustl(file_fe))// &
                  '\" index (i) u 1:3 w l lw 2 lt rgb \"grey\" title \"'// trim(adjustl(ch_titel)) //'\" ,\'

            write(FN, '(a)') ' \"'//trim(adjustl(file_fe))// &
                  '\" index (i) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i*' // trim(adjustl(ch_temp3)) // '+'// &
                  trim(adjustl(ch_temp2))// ')) '
         endif ! no_maxw
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_atomic_distribution_gnuplot


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

      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
                  '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i*' // trim(adjustl(ch_temp3)) // '+'// &
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
                  '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i*' // trim(adjustl(ch_temp3)) // '+'// &
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
            case (3:4)
               col = 3  ! column number after which orbital-resolved data start
            case default
               col = 2  ! column number after which orbital-resolved data start
         endselect
      endif

      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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
                      '",(i*' // trim(adjustl(ch_temp3)) //'+'// trim(adjustl(ch_temp2))// ')) ,\'
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
                      '\",(i*'// trim(adjustl(ch_temp3)) //'+'// trim(adjustl(ch_temp2))// ')) ,\'
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
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'set logscale y'
         write(FN, '(a)') 'set format y "10^{%L}"'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:200] "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i-1)'// &
                  '*' // trim(adjustl(ch_temp3)) // trim(adjustl(ch_temp2))// ')'
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'set logscale y'
         write(FN, '(a)') 'set format y \"10^{%L}\"'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:200] \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i-1)'// &
                  '*' // trim(adjustl(ch_temp3)) // trim(adjustl(ch_temp2))// ')'
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

      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
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




subroutine output_parameters_file(Scell,matter,laser,numpar,TB_Hamil,TB_Repuls,Err)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Pulse), dimension(:), intent(in) :: laser	! Laser pulse parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(TB_repulsive), intent(in) :: TB_Repuls  ! parameters of the repulsive part of TB
   type(TB_Hamiltonian), intent(in) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(Error_handling), intent(inout) :: Err	! error save
   !===============================================
   character :: path_sep
   integer INFO, FN, i, j, Nshl
   character(10) :: text(3) ! temporary variable
   character(100) :: chtemp(6)
   character(200) :: File_name, Error_descript
   logical :: file_opened

!    path_sep = numpar%path_sep
!    File_name = trim(adjustl(numpar%output_path))//path_sep
!    File_name = trim(adjustl(File_name))//'!OUTPUT_'//trim(adjustl(matter%Name))//'_Parameters.txt'
!    FN = 111
!    open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
!    inquire(file=trim(adjustl(File_name)),opened=file_opened)
   !if (.not.file_opened) then
!    print*, trim(adjustl(File_name)), file_opened, FN

   ! Check if Parameters file is opened, and try to open if it wasn't:
   call open_parameters_file(numpar, matter, FN, INFO)   ! below

   if (INFO < 0) then
      INFO = 2
      Error_descript = 'Error in output_parameters_file: file could not be opened, the program terminates'
      call Save_error_details(Err, INFO, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 9999
   endif

   numpar%FN_parameters = FN ! save this file number with parameters
#ifdef _OPENMP
   call Print_title(FN, Scell, matter, laser, numpar, 1) ! below
#else
#ifdef MPI_USED
   call Print_title(FN, Scell, matter, laser, numpar, 10) ! below
#endif
   call Print_title(FN, Scell, matter, laser, numpar, 4) ! below
#endif

   if (INFO >= 0) then
      write(chtemp(1), '(f12.2)') dble(matter%hw_plasma)
      write(FN, '(a,a,a)') ' Estimated plasmon energy: ', trim(adjustl(chtemp(1))), ' [eV]'

      write(FN, '(a)') ' Atomic data used for '//trim(adjustl(matter%Name))//' are:'

      do i = 1, matter%N_KAO
         write(FN,'(a)') trim(adjustl(m_dashline))
         write(chtemp(1), '(i12)') i
         write(chtemp(2), '(f6.2)') matter%Atoms(i)%percentage
         write(FN, '(a,$)') 'Element #'//trim(adjustl(chtemp(1)))//' is '//trim(adjustl(matter%Atoms(i)%Name))// &
            ' contributing to the compound with '//trim(adjustl(chtemp(2)))
         write(FN, '(a)') ''
         write(chtemp(1), '(i4)') INT(matter%Atoms(i)%Z)
         write(chtemp(2), '(es25.5)') matter%Atoms(i)%Ma
         write(FN, '(a,a,a,a,a)') 'Atomic number: ', trim(adjustl(chtemp(1))), ', mass: ', trim(adjustl(chtemp(2))), ' [kg]'

         write(chtemp(1), '(f12.2)') dble(matter%Atoms(i)%NVB)
         write(chtemp(2), '(f12.2)') dble(Scell(1)%Ne/Scell(1)%Na)
         write(FN, '(a,a,a,a,a)') 'Number of valence electrons: ', trim(adjustl(chtemp(1))), ' (band: ', &
                                    trim(adjustl(chtemp(2))), '/atom)'

         write(FN, '(a)') 'Shell#	Designator	Ne	Ip [eV]	Ek [eV]	Auger [fs]	'
         Nshl = size(matter%Atoms(i)%Ip)
         do j = 1, Nshl
            write(chtemp(1), '(i12)') j
            write(chtemp(2), '(i12)') matter%Atoms(i)%Shl_dsgnr(j)
            write(chtemp(3), '(f9.1)') matter%Atoms(i)%Ne_shell(j)
            write(chtemp(4), '(f9.1)') matter%Atoms(i)%Ip(j)
            write(chtemp(5), '(f9.1)') matter%Atoms(i)%Ek(j)
            if (matter%Atoms(i)%Auger(j) .LT. 1d10) then
               write(chtemp(6), '(es14.3)') matter%Atoms(i)%Auger(j)
            else
               write(chtemp(6), '(es14.5)') matter%Atoms(i)%Auger(j)
            endif
            write(FN, '(a,$)') trim(adjustl(chtemp(1)))//'	', trim(adjustl(matter%Atoms(i)%Shell_name(j)))//' 	', &
                                 trim(adjustl(chtemp(2)))//'	', trim(adjustl(chtemp(3)))//'	', &
                                 trim(adjustl(chtemp(4)))//'	', trim(adjustl(chtemp(5)))//'	', trim(adjustl(chtemp(6)))
            write(FN, '(a)') ''
         enddo !j
      enddo !i
      !write(FN,'(a)') '*************************************************************'
      write(FN,'(a)') trim(adjustl(m_starline))
   endif

   ! Save pulse parameters conversion, if required:
   call printout_fluence_dose_conversion(Scell(1), laser, numpar, matter, INFO, Err%File_Num)   ! below
   if (INFO < 0) then
      INFO = 9
      Error_descript = 'Specified photon energy does not allow for conversion of fluence into dose, the program terminates'
      call Save_error_details(Err, INFO, Error_descript)
      write(6, '(a)') trim(adjustl(m_dashline))
      write(6, '(a)') trim(adjustl(Error_descript))
      write(6, '(a)') 'The photon energy is too small, linear-absorption approximation does not work:'
      write(6, '(a)') 'Either specify the DOSE in the input file or increase photon energy'
      write(6, '(a)') trim(adjustl(m_dashline))
      goto 9999
   endif
9999 continue
end subroutine output_parameters_file


! Create the folder where the results will be storred:
subroutine create_output_folder(Scell, matter, laser, numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   integer i, iret
   character(200) :: File_name, File_name2, command, matter_name
   character(100) :: ch1, ch2, ch3, ch4
   logical :: file_exist

   ! Check embedding in water:
   if (numpar%embed_water) then ! if material embedded in water, add it to the name
      matter_name = trim(adjustl(matter%Name))//'_in_water'
   else  ! just material name
      matter_name = trim(adjustl(matter%Name))
   endif

   !------------------------------
   ! 1) If user defined the output name:
   UDN:if (LEN(trim(adjustl(numpar%output_extra_name))) > 0) then ! user-defined name
      write(File_name,'(a,a,a)') 'OUTPUT_', trim(adjustl(matter_name))//'_', trim(adjustl(numpar%output_extra_name))

   !------------------------------
   else UDN ! 2) If user did not defined the output name, construct the default name:

    LAS:if (maxval(laser(:)%F) .GT. 0.0d0) then

      if (allocated(laser(1)%Spectrum)) then    ! spectrum vs. single photon energy
         ch1 = 'spectrum'
      else ! single photon energy
         write(ch1,'(f7.1)') (laser(1)%hw)   ! photon energy
      endif

      if (laser(1)%KOP .EQ. 1) then
         write(ch2,'(f6.1)') (laser(1)%t*2.35482d0)	! pulse duration
      else
         write(ch2,'(f6.1)') laser(1)%t		! pulse duration
      endif

      if (laser(1)%F < 1.0e-2) then ! low dose
         write(ch3,'(f7.6)') laser(1)%F   ! dose [eV/atom]
      elseif (laser(1)%F < 1.0e4) then ! mideum dose
         write(ch3,'(f7.2)') laser(1)%F   ! dose [eV/atom]
      else ! high dose
         write(ch3,'(es12.2)') laser(1)%F ! dose [eV/atom]
      endif

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         if (size(laser) .GT. 1) then
            write(ch4,'(i2)') size(laser)
            write(File_name,'(a,a,a,a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw_', trim(adjustl(ch1)), '_t_', &
                  trim(adjustl(ch2)), '_F_', trim(adjustl(ch3)), '_', trim(adjustl(ch4)), '_pulses'
         else ! singe pulse
            write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw_', trim(adjustl(ch1)), '_t_', &
                  trim(adjustl(ch2)), '_F_', trim(adjustl(ch3))
         endif
      else ! it is linux
         if (size(laser) .GT. 1) then
            write(ch4,'(i2)') size(laser)
            write(File_name,'(a,a,a,a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw=', trim(adjustl(ch1)), '_t=', &
               trim(adjustl(ch2)), '_F=', trim(adjustl(ch3)), '_', trim(adjustl(ch4)), '_pulses'
         else ! singe pulse
            write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw=', trim(adjustl(ch1)), '_t=', &
               trim(adjustl(ch2)), '_F=', trim(adjustl(ch3))
         endif
      endif 
    else LAS ! no pulse
      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         do i = 1,size(Scell)
            write(ch1,'(f8.1)') Scell(i)%Te ! electron temperature [K]
            write(ch2,'(f8.1)') Scell(i)%Ta ! atomic temperature [K]
         enddo
         if (numpar%Nonadiabat) then
            write(ch3,'(a)') 'with_coupling'
         else
            write(ch3,'(a)') 'no_coupling'
         endif
         write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_Te_', trim(adjustl(ch1)), '_Ta_', &
            trim(adjustl(ch2)), '_', trim(adjustl(ch3))
      else ! it is linux
         do i = 1,size(Scell)
            write(ch1,'(f8.1)') Scell(i)%Te ! electron temperature [K]
            write(ch2,'(f8.1)') Scell(i)%Ta ! atomic temperature [K]
         enddo
         if (numpar%Nonadiabat) then
            write(ch3,'(a)') 'with_coupling'
         else
            write(ch3,'(a)') 'no_coupling'
         endif
         write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_Te=', trim(adjustl(ch1)), '_Ta=', &
            trim(adjustl(ch2)), '_', trim(adjustl(ch3))
      endif
    endif LAS

    ! Do TB and MD part only if we want (supercell is larger than 0):
    if (matter%cell_x*matter%cell_y*matter%cell_z .LE. 0) then
      write(File_name,'(a,a)') trim(adjustl(File_name)), '_MC_only'
    endif
   endif UDN


   ! User-defined addition to the output folder name, if any:
   if (LEN(trim(adjustl(numpar%output_name_add))) > 0) then
      write(File_name,'(a)') trim(adjustl(File_name))//'_'//trim(adjustl(numpar%output_name_add))
   endif


   File_name2 = File_name
   i = 0
#ifndef __GFORTRAN__
   ! for intel fortran compiler:
   inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
#else
   ! for gfortran compiler:
   inquire(FILE=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
#endif

   do while (file_exist)
      i = i + 1
      write(ch1,'(i6)') i
      write(File_name2,'(a,a,a)') trim(adjustl(File_name)), '_v', trim(adjustl(ch1))
#ifndef __GFORTRAN__
      ! for intel fortran compiler:
      inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
#else
      ! for gfortran compiler:
      inquire(FILE=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
#endif
   enddo
   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      command='md "'//trim(adjustl(File_name2))//'"' ! to create a folder use this command
      !CALL system(command)  ! create the folder
   else
      command='mkdir '//trim(adjustl(File_name2)) ! to create a folder use this command
      !CALL system(command)  ! create the folder
   endif
   iret = system(command)
   numpar%output_path = File_name2
end subroutine create_output_folder


! Coomunication with the used via file:
subroutine communicate(FN, time, numpar, matter)
   integer, intent(in) :: FN ! file number to read from
   real(8), intent(in) :: time ! current time [fs]
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   type(Solid), intent(inout) :: matter ! parameters of the material
   integer :: Reason, i, MOD_TIM, sz
   character(200) :: readline, given_line, File_name, error_part
   real(8) given_num
   logical :: read_well, read_well_2, file_opened, smth_read_master
   
   !--------------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything wrong here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      goto 7779
   endif

   smth_read_master = .false.

   File_name = numpar%Filename_communication
   inquire(UNIT=FN,opened=file_opened)
   if (file_opened) close(FN) ! for windows, we have to close the file to let the user write into it
   ! Check if the file was modified since the last time:
   call get_file_stat(trim(adjustl(File_name)), Last_modification_time=MOD_TIM) ! module 'Dealing_with_files'
   
   if (MOD_TIM /= numpar%MOD_TIME) then ! open file again only if it was modified by the user
      numpar%MOD_TIME = MOD_TIM ! save new time of the last modification
      open(UNIT=FN,FILE=trim(adjustl(File_name)),ERR=7777)
7777     continue ! in case if the program could not open the file
   endif
   
   inquire(UNIT=FN,opened=file_opened)

   COM_OPEN:if (file_opened) then ! read it
      rewind(FN,IOSTAT=Reason)  ! to start reading from the start
      i = 1 ! to start with
      read_well = .true.   ! to start with
      Reason = 1  ! to start with
      do while (Reason >= 0) ! read all lines if there is more than one
         call pars_comunications_file(FN, i, given_line, given_num, Reason) ! below
         if (Reason == 0) then
            call act_on_comunication(given_line, given_num, numpar, matter, time)   ! below
            smth_read_master = .true.
         endif
      enddo
      rewind(FN,IOSTAT=Reason)
      write(FN,'(a)',IOSTAT=Reason) ''
      rewind(FN,IOSTAT=Reason)

      call get_file_stat(trim(adjustl(File_name)), Last_modification_time=MOD_TIM) ! module 'Dealing_with_files'
      if (MOD_TIM /= numpar%MOD_TIME) then ! if it was modified by the user, then
         numpar%MOD_TIME = MOD_TIM         ! save new time of the last modification
      endif

      close(FN,ERR=7778) ! we have to close the file to let the user write into it
7778  continue ! in case if the program could not close the file
   endif COM_OPEN

7779 continue
!----------------------------
! If anything changed in the master-process, tell it to all the others:
#ifdef MPI_USED
   error_part = 'Error in "communicate":'
   call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {smth_read_master}', smth_read_master) ! module "MPI_subroutines"
   if (smth_read_master) then
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%MOD_TIME}', numpar%MOD_TIME) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%verbose}', numpar%verbose) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%t_total}', numpar%t_total) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%dt}', numpar%dt) ! module "MPI_subroutines"
      call reset_support_times(numpar)   ! above
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%dt}', numpar%dt) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%dt_save}', numpar%dt_save) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%NOMP}', numpar%NOMP) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%Transport}', numpar%Transport) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {matter%tau_bath}', matter%tau_bath) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {matter%T_bath}', matter%T_bath) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {numpar%Transport_e}', numpar%Transport_e) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {matter%tau_bath_e}', matter%tau_bath_e) ! module "MPI_subroutines"
      call broadcast_variable(numpar%MPI_param, trim(adjustl(error_part))//' {matter%T_bath_e}', matter%T_bath_e) ! module "MPI_subroutines"
   endif
#endif
end subroutine communicate



subroutine save_duration(matter, numpar, chtext, ctim)
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   character(*), intent(in) :: chtext ! time duration to print out
   integer, dimension(8), intent(in), optional :: ctim
   !---------------------------------
   integer FN, INFO
   logical file_opened, file_exists
   character(200) :: File_name
   character(1) path_sep

   call open_parameters_file(numpar, matter, FN, INFO)   ! below

   if (INFO < 0) return

   if (present(ctim)) then ! print out the absolute time
      call print_time('Started  at', ctim=ctim, FN_in=FN) ! module "Little_subroutines"
      call print_time('Finished at', FN_in=FN) ! module "Little_subroutines"
   endif

   write(FN,'(a,a)') 'Duration of execution of program: ', trim(adjustl(chtext))
   write(FN,'(a)') trim(adjustl(m_starline))
end subroutine save_duration


subroutine open_parameters_file(numpar, matter, FN, INFO)
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   type(Solid), intent(in) :: matter ! parameters of the material
   integer, intent(out) :: FN    ! index of the file with parameters
   integer, intent(out) :: INFO  ! -1 = no output directory; 0=file opened, 1=file created and opened
   !-------------------------
   character(500) :: File_name
   character(1) :: path_sep
   logical :: file_opened, file_exists

   FN = 0   ! no file
   INFO = 0 ! to start with

   ! Check if output directory exists:
   if (LEN(trim(adjustl(numpar%output_path))) == 0) then
      INFO = -1   ! no output directory
      return   ! nothing else to do
   else

#ifndef __GFORTRAN__
      ! for intel fortran compiler:
      inquire(DIRECTORY=trim(adjustl(numpar%output_path)),exist=file_exists)    ! check if input directory excists
#else
      ! for gfortran compiler:
      inquire(FILE=trim(adjustl(numpar%output_path)),exist=file_exists)    ! check if input directory excists
#endif

      if (.not.file_exists) then
         INFO = -1   ! no output directory
         return   ! nothing else to do
      endif
   endif

   ! If directory exists, check the existence of the Parameters file:
   path_sep = trim(adjustl(numpar%path_sep))
   File_name = trim(adjustl(numpar%output_path))//path_sep
   File_name = trim(adjustl(File_name))//'!OUTPUT_'//trim(adjustl(matter%Name))//'_Parameters.txt'
   inquire(file=trim(adjustl(File_name)),exist=file_exists)
   if (.not.file_exists) then ! no such file exists, create it:
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
      INFO = 1
      return
   endif

   ! Check if file is already opened, and if yes, get its number:
   inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN)

   ! If file is not opened, open it to a new number:
   if (.not.file_opened) then
      !open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), status="old", position="append", action="write")
   endif
end subroutine open_parameters_file



subroutine printout_fluence_dose_conversion(Scell, laser, numpar, matter, INFO, FN_err)
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   type(Solid), intent(in) :: matter ! parameters of the material
   integer, intent(inout) :: INFO   ! flag
   integer, intent(in) :: FN_err ! error-file number
   !---------------------------
   character(30) :: F_in, Dose, PN, hw
   integer :: i, Nsiz, FN, i_CS

   INFO = 0 ! to start with
   ! How many laser pulses:
   Nsiz = size(laser)

   ! Which inelastic cross section to use (BEB vs CDF):
   do i = 1, size(matter%Atoms)
      i_CS = minval(matter%Atoms(i)%TOCS(:))
   enddo

   ! To write into Parameters file, find or open it:
   call open_parameters_file(numpar, matter, FN, INFO)   ! above

   ! For each laser pulse:
   do i = 1, Nsiz
      !print*, 'printout_fluence_dose_conversion:', i, laser(i)%F_in, laser(i)%F

      if (laser(i)%F_in > 0.0d0) then
         write(hw, '(f24.3)') laser(i)%hw ! photon energy as a text
         write(PN, '(i0)') i  ! pulse index as a text

         if ( (laser(i)%F_in > 1.0d6) .or. (laser(i)%F_in < 1.0d-6) ) then
            write(F_in, '(es26.10)') laser(i)%F_in
         else
            write(F_in, '(f24.6)') laser(i)%F_in
         endif

         if ( (laser(i)%F > 1.0d6) .or. (laser(i)%F < 1.0d-6) ) then
            write(Dose, '(es26.10)') laser(i)%F
         else
            write(Dose, '(f24.6)') laser(i)%F
         endif

         !---------------------
         ! Print on the screen:
         write(6, '(a)') ' Pulse #'//trim(adjustl(PN))//': incoming fluence '//trim(adjustl(F_in))//' [J/cm^2]'
         write(6, '(a)') ' corresponds to absorbed dose '//trim(adjustl(Dose))//' [eV/atom]'
         write(6, '(a)') trim(adjustl(m_starline))
         ! Only for EADL atomic cross section, make a worning (but not for CDF):

         if (laser(i)%hw < Scell%E_gap) then
            INFO = -1
            write(6, '(a)') 'ERROR: Photon energy is too low (<E_gap): '//trim(adjustl(hw))//' [eV]'
            write(6, '(a)') 'Conversion from incoming fluence to dose CANNOT be done!'
         elseif ( (laser(i)%hw < 30.0d0) .and. (i_CS < 1) ) then ! Print warning for too low photon energy:
            INFO = 1
            !write(6, '(a)') 'WARNING: Photon energy is too low (<30 eV): '//trim(adjustl(hw))//' [eV]'
            !write(6, '(a)') 'Conversion from incoming fluence to dose may not work well!'
            call printout_warning(6, 5, text_to_add=trim(adjustl(hw)) ) ! module "Read_input_data"
         endif

         ! Print in the file too:
         if (INFO >= 0) then
            write(FN, '(a)') ' Pulse #'//trim(adjustl(PN))//': incoming fluence '//trim(adjustl(F_in))//' [J/cm^2]'
            write(FN, '(a)') ' corresponds to absorbed dose '//trim(adjustl(Dose))//' [eV/atom]'
            write(FN, '(a)') trim(adjustl(m_starline))
            ! Print warning for too low photon energy:
            ! Only for EADL atomic cross section, make a worning (but not for CDF):
            if (laser(i)%hw < Scell%E_gap) then
               write(FN, '(a)') 'ERROR: Photon energy is too low (<E_gap): '//trim(adjustl(hw))//' [eV]'
               write(FN, '(a)') 'Conversion from incoming fluence to dose CANNOT be done!'
            elseif ( (laser(i)%hw < 30.0d0) .and. (i_CS < 1) ) then ! Print warning for too low photon energy:
               !write(FN, '(a)') 'WARNING: Photon energy is too low (<30 eV): '//trim(adjustl(hw))//' [eV]'
               !write(FN, '(a)') 'Conversion from incoming fluence to dose may not work well!'
               call printout_warning(FN, 5, text_to_add=trim(adjustl(hw)) ) ! module "Read_input_data"
               ! Annd in the error-file too:
               call printout_warning(FN_err, 5, text_to_add=trim(adjustl(hw)) ) ! module "Read_input_data"
            endif
         endif

      endif
   enddo
end subroutine printout_fluence_dose_conversion



subroutine reset_dt(numpar, matter, tim_cur)
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   type(Solid), intent(inout) :: matter ! parameters of the material
   real(8), intent(in) :: tim_cur        ! current time step of the simulation
   real(8) :: est
   est = 1.0d-6 ! precision

   ! Simulation time step:
   if ((numpar%i_dt > 0) .and. (numpar%i_dt <= size(numpar%dt_MD_grid))) then   ! only if there is an option to change dt
      if (tim_cur >= numpar%dt_MD_reset_grid(numpar%i_dt)-est) then ! time to change dt
         numpar%dt = numpar%dt_MD_grid(numpar%i_dt)              ! to this value
         numpar%i_dt = numpar%i_dt + 1 ! next step to read from
         call reset_support_times(numpar)   ! below
         if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'Time-step of MD simulation is changed to', numpar%dt
         endif
      endif
   elseif (numpar%i_dt == 0) then   ! its before the simulation start, reset the starting time
      numpar%i_dt = numpar%i_dt + 1 ! next step to read from
      !numpar%t_start = numpar%dt_MD_reset_grid(1)   ! to start from
      numpar%dt = numpar%dt_MD_grid(1)           ! to start from
      call reset_support_times(numpar)   ! below
   endif

   ! Atomic thermostat parameters:
   if (allocated(numpar%At_bath_reset_grid) .and. (numpar%i_At_bath_dt <= size(numpar%At_bath_grid_Ta)) ) then
      if (tim_cur >= numpar%At_bath_reset_grid(numpar%i_At_bath_dt)-est) then ! time to change dt

         matter%T_bath = numpar%At_bath_grid_Ta(numpar%i_At_bath_dt) ! new bath temperature [K]
         matter%T_bath = matter%T_bath/g_kb  ! [eV] thermostat temperature for atoms

         matter%tau_bath = numpar%At_bath_grid_tau(numpar%i_At_bath_dt) ! new characteristic time [fs]
         if ((matter%tau_bath > 1.0d14) .or. (matter%T_bath < -1.0d-6)) then  ! there is no bath, too slow to couple
            numpar%Transport = .false. ! excluded
            if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
               print*, 'Atomic thermostat is off'
            endif
         else
            numpar%Transport = .true.	 ! included
            if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
               print*, 'Atomic thermostat parameters are changed to', &
                  matter%T_bath*g_kb, matter%tau_bath
            endif
         endif

         numpar%i_At_bath_dt = numpar%i_At_bath_dt + 1 ! next step to read from
      endif
   endif

   ! Electronic thermostat  parameters:
   if (allocated(numpar%El_bath_reset_grid) .and. (numpar%i_El_bath_dt <= size(numpar%El_bath_grid_Ta)) ) then
      if (tim_cur >= numpar%El_bath_reset_grid(numpar%i_El_bath_dt)-est) then ! time to change dt

         matter%T_bath_e = numpar%El_bath_grid_Ta(numpar%i_El_bath_dt) ! new bath temperature [K]
         matter%T_bath_e = matter%T_bath_e/g_kb  ! [eV] thermostat temperature for atoms

         matter%tau_bath_e = numpar%El_bath_grid_tau(numpar%i_El_bath_dt) ! new characteristic time [fs]
         if ((matter%tau_bath_e > 1.0d14) .or. (matter%T_bath_e < -1.0d-6)) then  ! there is no bath, too slow to couple
            numpar%Transport_e = .false. ! excluded
            if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
               print*, 'Electronic thermostat is off'
            endif
         else
            numpar%Transport_e = .true.	 ! included
            if (g_numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
               print*, 'Electronic thermostat parameters are changed to', &
                  matter%T_bath_e*g_kb, matter%tau_bath_e
            endif
         endif

         numpar%i_El_bath_dt = numpar%i_El_bath_dt + 1 ! next step to read from
      endif
   endif
end subroutine reset_dt


pure subroutine reset_support_times(numpar)
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   numpar%halfdt = numpar%dt/2.0d0      ! dt/2, often used
   numpar%dtsqare = numpar%dt*numpar%halfdt ! dt*dt/2, often used
   numpar%dt3 = numpar%dt**3/6.0d0            ! dt^3/6, often used
   numpar%dt4 = numpar%dt*numpar%dt3/8.0d0    ! dt^4/48, often used
end subroutine reset_support_times


subroutine act_on_comunication(given_line, given_num, numpar, matter, time)
   !logical, intent(in) :: read_well ! did we read something meaningful from the comunication file?
   character(*), intent(in) :: given_line ! line read from the file
   real(8), intent(in) :: given_num  ! number read from the file
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   type(Solid), intent(inout) :: matter ! parameters of the material
   real(8), intent(in) :: time ! current time [fs]
   integer FN, noth, lngt, INFO
   logical file_opened, read_well
   character(200) :: File_name, temp1, temp2, given_line_processed
   character(1) path_sep

   read_well = .true.   ! by default, we could read everything well; change later if problem
   path_sep = trim(adjustl(numpar%path_sep))
   lngt = LEN(trim(adjustl(given_line)))    ! length of the line

   ! Check if the last character is TAB:
   if (given_line(lngt:lngt) == char(9) ) then  ! remove the TAB from the line
      given_line_processed = given_line(1:lngt-1)
   else
      given_line_processed = given_line
   endif

   if (read_well) then
!       File_name = trim(adjustl(numpar%output_path))//path_sep
!       File_name = trim(adjustl(File_name))//'!OUTPUT_'//trim(adjustl(matter%Name))//'_Parameters.txt'
!       inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN)
!       if (.not.file_opened) then
!          FN = 300
!          !open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
!          open(UNIT=FN, FILE = trim(adjustl(File_name)), status="old", position="append", action="write")
!       endif
      call open_parameters_file(numpar, matter, FN, INFO)   ! below
      if (INFO < 0) then
         write(6,'(a)') 'act_on_comunication Error: Could not open file with parameters to write in to'
         write(6,'(a)') 'Cannot change the parameters, proceeding as before...'
         return
      endif

      select case(trim(adjustl(given_line_processed)))
      case ('verbose', 'VERBOSE', 'Verbose')
         numpar%verbose = .true.
         write(6,'(a)') 'Verbose option on: XTANT will print a lot of markers for testing and debugging'
         write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' verbose option was switched on'
      !-------------------
      case ('verbose_off', 'VERBOSE_OFF', 'Verbose_off', 'no_verbose', 'noverbose', 'no-verbose')
         numpar%verbose = .false.
         write(6,'(a)') 'Verbose option off: XTANT will not print markers'
         write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' verbose option was switched off'

      !-------------------
      case ('time', 'TIME', 'Time', 'TIme', 'TIMe', 'tIme', 'emit', 'Vremya')
         numpar%t_total = given_num ! total duration of simulation [fs]
         write(6,'(a,f10.3)') 'Duration of simulation is changed to', given_num
         write(FN,'(a,f10.3,a,f10.3)') 'At time instance of ', time, ' duration of simulation is changed to ', given_num

      !-------------------
      case ('MDdt', 'dtMD', 'mddt', 'dtmd', 'MDDT', 'DTMD', 'MD_dt', 'md_dt')
         numpar%dt = given_num ! Time step for MD [fs]
         call reset_support_times(numpar)   ! above
         !numpar%halfdt = numpar%dt/2.0d0      ! dt/2, often used
         !numpar%dtsqare = numpar%dt*numpar%halfdt ! dt*dt/2, often used
         !numpar%dt3 = numpar%dt**3/6.0d0            ! dt^3/6, often used
         !numpar%dt4 = numpar%dt*numpar%dt3/8.0d0    ! dt^4/48, often used
         write(6,'(a,f9.3)') 'Time-step of MD simulation is changed to', given_num
         write(FN,'(a,f10.3,a,f9.3)') 'At time instance of ', time, ' time-step of MD simulation is changed to ', given_num 

      !-------------------
      case ('SAVEdt', 'savedt', 'dtsave', 'dtSAVE', 'Savedt', 'SaveDT', 'SaveDt', 'SAVE_dt', 'Save_dt', 'save_dt')
         numpar%dt_save = given_num ! save data into files every 'dt_save_time' [fs]
         write(6,'(a,f9.3)') 'Time-step of saving output files is changed to', given_num
         write(FN,'(a,f10.3,a,f9.3)') 'At time instance of ', time, ' time-step of saving output files is changed to ', given_num
      
      !-------------------
      case ('OMP', 'omp', 'NOMP', 'nomp', 'Nomp', 'N_OMP', 'n_omp')
         ! Reset the OpenMP parallelization options:
         numpar%NOMP = given_num
         write(temp1,'(f12.3)') time
         write(temp2,'(i10)') INT(given_num)
         
#ifdef _OPENMP
         noth = OMP_GET_MAX_THREADS()   ! to chech if the function worked
         call set_OMP_number( numpar%NOMP, .true., 6, 'Reset number of threads in OpenMP to '//trim(adjustl(temp2)) )    ! below
         if ( noth /= OMP_GET_MAX_THREADS() ) then
            write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs], number of threads in OpenMP is changed to ', &
               trim(adjustl(temp2))
         else
            write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), &
               '[fs]: unsuccessful attempt to change number of threads in OpenMP to ',  trim(adjustl(temp2))
            write(6,'(a)') 'Number of threads in OpenMP is unchanged: ',  trim(adjustl(temp2))
         endif
#else
         write(FN,'(a)') ' The code compiled without OpenMP, cannot set parallelization'
         write(6,'(a)') 'The code compiled without OpenMP, cannot set parallelization'
#endif

      !-------------------
      case ('Thermostat_dt_a', 'THERMOSTAT_DT_A', 'thermostat_dt_a')
         if (given_num < 0.0d0) then
            numpar%Transport = .false. ! excluded atomic thermostat
            write(6,'(a)') 'Atomic thermostat is switched off'
            write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' atomic thermostat is switched off'
         else
            numpar%Transport = .true. ! included atomic thermostat
            matter%tau_bath = given_num   ! [fs] time constant of cooling for atoms
            write(6,'(a,f12.3)') 'Atomic thermostat time is changed to', given_num
            write(FN,'(a,f10.3,a,f12.3)') 'At time instance of ', time, ' atomic thermostat time is changed to ', given_num
         endif
      !-------------------
      case ('Thermostat_Ta', 'THERMOSTAT_Ta', 'thermostat_Ta')
         if (given_num < 0.0d0) then
            numpar%Transport = .false. ! excluded atomic thermostat
            write(6,'(a)') 'Atomic thermostat is switched off'
            write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' atomic thermostat is switched off'
         else
            numpar%Transport = .true. ! included atomic thermostat
            matter%T_bath = given_num   ! [K] bath temperature for atoms
            matter%T_bath = matter%T_bath/g_kb  ! [eV] thermostat temperature
            write(6,'(a,f12.3)') 'Atomic thermostat temperature is changed to', given_num
            write(FN,'(a,f10.3,a,f12.3)') 'At time instance of ', time, ' atomic thermostat temperature is changed to ', given_num
         endif

      !-------------------
      case ('Thermostat_dt_e', 'THERMOSTAT_DT_E', 'thermostat_dt_e')
         if (given_num < 0.0d0) then
            numpar%Transport_e = .false. ! excluded atomic thermostat
            write(6,'(a)') 'Electronic thermostat is switched off'
            write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' electronic thermostat is switched off'
         else
            numpar%Transport_e = .true. ! included atomic thermostat
            matter%tau_bath_e = given_num   ! [fs] time constant of cooling for atoms
            write(6,'(a,f10.3)') 'Electronic thermostat time is changed to', given_num
            write(FN,'(a,f10.3,a,f12.3)') 'At time instance of ', time, ' electronic thermostat time is changed to ', given_num
         endif
      !-------------------
      case ('Thermostat_Te', 'THERMOSTAT_Te', 'thermostat_Te')
         if (given_num < 0.0d0) then
            numpar%Transport_e = .false. ! excluded atomic thermostat
            write(6,'(a)') 'Electronic thermostat is switched off'
            write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' electronic thermostat is switched off'
         else
            numpar%Transport_e = .true. ! included atomic thermostat
            matter%T_bath_e = given_num   ! [K] bath temperature for atoms
            matter%T_bath_e = matter%T_bath_e/g_kb  ! [eV] thermostat temperature
            write(6,'(a,f12.3)') 'Electronic thermostat temperature is changed to', given_num
            write(FN,'(a,f10.3,a,f12.3)') 'At time instance of ', time, ' electronic thermostat temperature is changed to ', given_num
         endif

      !-------------------
      case default
         print*, 'Could not interpret what is read from the file: ', trim(adjustl(given_line)), given_num
      end select
   else
      print*, 'Could not read well from the file: ', trim(adjustl(given_line)), given_num
   endif
end subroutine act_on_comunication



subroutine set_OMP_number(NOMP, prnt, FN, lin)
   integer, intent(inout) :: NOMP  ! number of threads to be set; negative means = maximal threads available
   logical, intent(in) :: prnt  ! do we want to print out anything?
   integer, intent(in) :: FN    ! file number to print into
   character(*), intent(in), optional :: lin    ! a line to print out
   !------------------------------------
   character(10) :: temp2
   
#ifdef _OPENMP
   call OMP_SET_DYNAMIC(0) ! standard openmp subroutine
   if (NOMP <= 0) then ! use all available processors / threads:
      NOMP = OMP_GET_MAX_THREADS() ! number of threads for openmp defined in INPUT_PARAMETERS.txt
   endif
   call OMP_SET_NUM_THREADS(NOMP) ! number of threads for openmp defined in INPUT_PARAMETERS.txt
   if (prnt) then
      if (present(lin)) then
         write(FN,'(a)') trim(adjustl(lin))
      else  ! printout default message
         write(temp2,'(i10)') INT(NOMP)
         write(FN,'(a,a)') ' The code was compiled with OpenMP parallelization, THREADS: ', trim(adjustl(temp2))
      endif
   endif
#else
   if (prnt) then
      if (present(lin)) then
         write(FN,'(a)') trim(adjustl(lin))
      else  ! printout default message
         write(FN,'(a)') ' The code was compiled without using OpenMP'
      endif
   endif
#endif
end subroutine set_OMP_number



subroutine pars_comunications_file(FN, i, out_line, out_num, Reason)
   integer, intent(in) :: FN
   integer, intent(inout) :: i
   character(*), intent(out) :: out_line
   real(8), intent(out) :: out_num
   integer, intent(out) :: Reason
   !---------------------------------
   character(200) :: read_line
   logical :: read_well
   read_well = .false.
   out_line = ''
   out_num = 0.0d0

   read(FN, '(a)', IOSTAT=Reason) read_line
   call pars_comunications(read_line, out_line, out_num, read_well)  ! below
end subroutine pars_comunications_file


subroutine pars_comunications(readline, out_line, out_num, read_well)
   character(*), intent(in) :: readline
   character(*), intent(out) :: out_line
   real(8), intent(out) :: out_num
   logical, intent(out) :: read_well
   !---------------------------------
   integer :: Reason, i
   read_well = .false.
   out_line = ''
   out_num = 0.0d0

   i = 1    ! to start with
   read(readline, *, IOSTAT=Reason) out_line, out_num
   !print*, 'Reason', Reason, out_line, out_num
   if (Reason /= 0) then ! try again but with one variable
      read(readline, *, IOSTAT=Reason) out_line
      !print*, 'Reason2', Reason, out_line
   endif
   call read_file(Reason, i, read_well)  ! module "Dealing_with_files"
   if (Reason .LT. 0) then
      print*, 'No descriptor or value found in the communication file'
   else if (Reason .GT. 0) then
      print*, 'Given number interpreted as', out_num, ', it does not match the variable type'
   endif
   if (.not.read_well) then
      print*, 'Wrong format of input in Comunication, could not interpret.'
   endif
end subroutine pars_comunications



subroutine Print_title(print_to, Scell, matter, laser, numpar, label_ind)
   integer, intent(in) :: print_to ! the screen, or file
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! material parameters
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   integer, intent(in) :: label_ind ! which label to print
   !type(TB_repulsive), dimension(:), intent(in) :: TB_Repuls  ! parameters of the repulsive part of TB
   !type(TB_Hamiltonian), dimension(:), intent(in) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !--------------
   integer :: i, j
   character(100) :: text, text1, text2, text3
   logical :: optional_output
   real(8) :: lambda, temp

   !--------------------------------------------------------------------
   ! Make sure non-master MPI processes aren't doing anything wrong here
   if (numpar%MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return
   endif
   !--------------------------------------------------------------------

   write(print_to,'(a)') trim(adjustl(m_starline))
   call XTANT_label(print_to, label_ind) ! below
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '*  XTANT: X-ray-induced Thermal And Nonthermal Transitions  '
   write(print_to,'(a,a)') '* Current version of the code: ', trim(adjustl(m_XTANT_version))
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '  A hybrid approach consisting of: '
   write(print_to,'(a)') ' (1) Monte Carlo '
   write(print_to,'(a)') ' (2) Boltzmann collision integrals '
   write(print_to,'(a)') ' (3) Transferable Tight Binding '
   write(print_to,'(a)') ' (4) Molecular Dynamics '
   write(print_to,'(a)') ' DOI: https://doi.org/10.5281/zenodo.8392569'
   write(print_to,'(a)') trim(adjustl(m_starline))

   !ooooooooooooooooooooooooooooooooooooooooooooo
   write(print_to,'(a)') '  Calculations performed for the following parameters:'
   write(print_to,'(a,a)') ' Target material: ', trim(adjustl(matter%Name))
   write(print_to,'(a)') ' Chemical formula interpreted as: '
   do i = 1, size(matter%Atoms)
      write(text,'(f12.6)') matter%Atoms(i)%percentage
      write(text1,'(i3)') INT(matter%Atoms(i)%Z)
      write(print_to,'(a,a,a,a,a)') ' '//trim(adjustl(text)), ' of ', trim(adjustl(matter%Atoms(i)%Name)), &
         ' (element #', trim(adjustl(text1))//')'
   enddo

   if (numpar%embed_water) then
      write(text1,'(i6)') numpar%N_water_mol
      write(print_to,'(a)') ' (Note that the material was embedded in water with # of molecules: '//trim(adjustl(text1))//')'
   endif

   do i = 1, size(Scell)
      if (numpar%fe_input_exists) then
         write(print_to,'(a,a)') ' Initial electron distribution read from file: ', trim(adjustl(numpar%fe_filename))
      else
         write(print_to,'(a,f12.3,a)') ' Initial electron temperature	' , Scell(i)%Te, ' [K]'
      endif
      write(print_to,'(a,f12.3,a)') ' Initial atomic temperature	' , Scell(i)%Ta, ' [K]'
   enddo
   if ((size(laser) == 0) .or. (maxval(laser(:)%t) <= 0.0d0) .or. (maxval(laser(:)%hw) <= 0.0d0) .or. (maxval(laser(:)%F) <= 0.0d0)) then
      write(print_to,'(a)') ' No FEL-pulse is calculated'
   else
      write(text, '(i10)') size(laser)
      write(print_to,'(a,a)') ' Number of FEL-pulses included  ', trim(adjustl(text))
      do i = 1, size(laser)
         if (size(laser) .GT. 1) then
            write(print_to,'(a,i2,a)') ' Parameters of the pulse number #', i, ' are:' 
         endif
         write(print_to,'(a,f12.3,a)') ' Photon energy ' , laser(i)%hw, ' [eV]'

         ! Get the wavelength corresponding to the given photon energy:
         lambda = convert_hw_to_wavelength(laser(i)%hw)  ! module "Little_subroutines"
         write(print_to,'(a,f12.3,a)') ' Corresponding to the wavelength ' , lambda, ' [nm]'
         !write(print_to,'(a,f12.3,a)') ' Corresponding to photon energy  ' , convert_wavelength_to_hw(lambda), ' [eV]'

         if (laser(i)%FWHM_hw > 0.0d0) then
            write(print_to,'(a,f12.3,a)') ' with the spectral width ' , laser(i)%FWHM_hw, ' [eV] FWHM'
         else
            write(print_to,'(a)') ' Monochromatic pulse is used'
         endif

         select case (laser(i)%KOP)
         case (0)
            write(print_to,'(a)') ' Flat top pulse is used with'
            write(print_to,'(a,f12.3,a)') ' Pulse duration ' , laser(i)%t, ' [fs]'
         case (2)
            write(print_to,'(a)') ' SASE pulse is used with'
            write(print_to,'(a,f12.3,a)') ' Pulse duration ' , laser(i)%t, ' [fs]'
         case default
            write(print_to,'(a)') ' Gaussian pulse is used with'
            write(print_to,'(a,f12.3,a)') ' Pulse duration ' , laser(i)%t*2.35482, ' [fs]'
         end select
         write(print_to,'(a,f12.3,a)') ' Pulse maximum at ' , laser(i)%t0, ' [fs]'
         write(print_to,'(a,f12.5,a)') ' Absorbed dose    ' , laser(i)%F, ' [eV/atom]'
      enddo
   endif ! FEL included or not?
   
   !ooooooooooooooooooooooooooooooooooooooooooooo
   write(print_to,'(a)') trim(adjustl(m_starline))
   SCL:do i = 1, size(Scell)
      select case (abs(numpar%optic_model))
      case (1) ! within the Drude model
         write(print_to,'(a)') '  Probe-pulse is calculated within Drude model'
         write(print_to,'(a)') ' with the following parameters of the probe:'
         write(print_to,'(a, f7.1, a, f5.1, a)') ' Wavelength: ', Scell(i)%eps%l, '[nm]; Angle:', &
            Scell(i)%eps%teta/g_Pi*(180.0d0), '[degrees]'
         write(print_to,'(a, f7.1, a)') ' Thickness of the sample: ', Scell(i)%eps%dd, ' [nm]'
         write(print_to,'(a, es12.3, es12.3)') ' Effective mass of electron and hole: ', Scell(i)%eps%me_eff, Scell(i)%eps%mh_eff
         write(print_to,'(a, es12.3, es12.3)') ' Effective scattering time of electron and of hole: ', Scell(i)%eps%tau_e, Scell(i)%eps%tau_h
      case (2:3)  ! Trani model
         write(print_to,'(a)') '  Probe-pulse is calculated with RPA (Trani et al.) approach ' ! [PRB 72, 075423 (2005)]'
         write(print_to,'(a)') ' with the following parameters of the probe:'
         write(print_to,'(a, f7.1, a, f5.1, a)') ' Wavelength: ', Scell(i)%eps%l, ' [nm]; Angle:', &
            Scell(i)%eps%teta/g_Pi*(180.0d0), '    [degrees]'
         write(print_to,'(a, f7.1, a)') ' Thickness of the sample: ', Scell(i)%eps%dd, ' [nm]'
         if (numpar%optic_model .EQ. 2) then
            write(text1, '(i10)') numpar%ixm
            write(text2, '(i10)') numpar%iym
            write(text3, '(i10)') numpar%izm
            if (allocated(numpar%k_grid)) then
               write(print_to,'(a,a,a,a,a,a)') ' Number of k-points (on user-defined grid): ', &
                  trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
            else
               write(print_to,'(a,a,a,a,a,a)') ' Number of k-points (on Monkhorst-Pack grid): ', &
                  trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
            endif
         else
            write(print_to,'(a)') ' Calculations are performed for Gamma-point'
         endif
      case (4:5) ! Graf-Vogl model
         if (numpar%optic_model < 0) then
            write(print_to,'(a)') '  Probe-pulse is calculated with Graf-Vogl approach'
         else if (numpar%optic_model == 5) then
            write(print_to,'(a)') '  Probe-pulse is calculated with Kubo-Greenwood non-orthogonal'
         else
            write(print_to,'(a)') '  Probe-pulse is calculated with Kubo-Greenwood orthogonalized'
         endif
         write(print_to,'(a)') ' with the following parameters of the probe:'
         write(print_to,'(a, f7.1, a, f5.1, a)') ' Wavelength: ', Scell(i)%eps%l, ' [nm]; Angle:', &
            Scell(i)%eps%teta/g_Pi*(180.0d0), '    [degrees]'
         write(print_to,'(a, f7.1, a)') ' Thickness of the sample: ', Scell(i)%eps%dd, ' [nm]'
         if (numpar%ixm*numpar%iym*numpar%izm .EQ. 1) then
            write(print_to,'(a)') ' Calculations are performed for Gamma-point'
         else
            write(text1, '(i10)') numpar%ixm
            write(text2, '(i10)') numpar%iym
            write(text3, '(i10)') numpar%izm
            if (allocated(numpar%k_grid)) then
               write(print_to,'(a,a,a,a,a,a)') ' Number of k-points (on user-defined grid): ', &
                  trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
            else
               write(print_to,'(a,a,a,a,a,a)') ' Number of k-points (on Monkhorst-Pack grid): ', &
                  trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
            endif
         endif
      case default ! no optical coefficients needed
         write(print_to,'(a)') '  No probe-pulse is calculated'
      end select

      if (Scell(i)%eps%all_w) then
         if (Scell(i)%eps%KK) then
            write(print_to,'(a)') ' The spectrum is calculated via Im(CDF) using Kramers Kronig relations'
         else
            write(print_to,'(a)') ' The spectrum is calculated directly for both, Re(CDF) and Im(CDF)'
         endif
      endif
   enddo SCL

   !ooooooooooooooooooooooooooooooooooooooooooooo
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '  The model parameters used are:'
   write(text, '(f15.5)') numpar%t_total
   write(print_to,'(a,a,a)') ' Duration of modelling ' , trim(adjustl(text)), ' [fs]'

   if (matter%cell_x*matter%cell_y*matter%cell_z .LE. 0) then
      write(print_to,'(a)') ' TBMD part is switched off, only MC modelling is performed'
   else
      write(print_to,'(a)') ' Tight Binding parametrization schemes used are:'         
      call print_Hamiltonian_info(print_to, Scell(1), matter)  ! below
      !write(print_to,'(a,a)') ' Hamiltonian:      ', trim(adjustl(Scell(1)%TB_Hamil(1,1)%Param))
      !write(print_to,'(a,a)') ' Repulsive energy: ', trim(adjustl(Scell(1)%TB_Repuls(1,1)%Param))
      
      ASSOCIATE (ARRAY => Scell(1)%TB_Hamil(1,1)) ! this is the sintax we have to use to check the class of defined types
         select type(ARRAY)
         type is (TB_H_DFTB) ! TB parametrization
            write(print_to,'(a,a)') ' With skf-parametrization: (', trim(adjustl(ARRAY%param_name))//')'
            select case (numpar%basis_size_ind)
            case (0)
               text = 's'
            case (1)    ! sp3
               text = 'sp3'
            case default    ! sp3d5
               text = 'sp3d5'
            endselect
            write(print_to,'(a,a)') ' With the basis set: ', trim(adjustl(text))
          type is (TB_H_3TB) ! TB parametrization
            select case (numpar%basis_size_ind)
            case (0)
               text = 's'
            case (1)    ! sp3
               text = 'sp3'
            case default    ! sp3d5
               text = 'sp3d5'
            endselect
            write(print_to,'(a,a)') ' With the basis set: ', trim(adjustl(text))
            if (ARRAY%include_3body) then
               write(print_to,'(a,a)') ' With 3-body terms included'
            else
               write(print_to,'(a,a)') ' Only 2-body terms included (no 3-body terms)'
            endif
          type is (TB_H_BOP) ! TB parametrization
            select case (numpar%basis_size_ind)
            case (0)
               text = 's'
            case (1)    ! sp3
               text = 'sp3'
            case default    ! sp3d5
               text = 'sp3d5'
            endselect
            write(print_to,'(a,a)') ' With the basis set: ', trim(adjustl(text))
          type is (TB_H_xTB) ! TB parametrization
            select case (numpar%basis_size_ind)
            case (0)    ! s
               text = 'Cartesian s'
            case (1)    ! s s*
               text = 'Cartesian ss*'
            case (2)    ! sp3
               text = 'Cartesian sp3'
            case (3)    ! sp3s*
               text = 'Cartesian sp3s*'
            case (4)    ! sp3d5
               text = 'Cartesian sp3d6'
            case (5)    ! sp3d5s*
               text = 'Cartesian sp3d6s*'
            endselect
            write(print_to,'(a,a)') ' With the basis set: ', trim(adjustl(text))
            write(print_to,'(a,i1,a)') ' Using ', ARRAY%Nprim, ' GTO for STO'
            write(print_to,'(a)') ' The following orbitals are considerred for the elements: '
            do i = 1, size(matter%Atoms)
               ASSOCIATE (ARRAY2 => Scell(1)%TB_Hamil(i,i))
                  select type(ARRAY2)
                  type is (TB_H_xTB) ! TB parametrization
                     write(print_to,'(a)') ' For '//trim(adjustl(matter%Atoms(i)%Name))//': '//trim(adjustl(ARRAY2%AO_names))
                  endselect
               END ASSOCIATE
            enddo
         endselect
      END ASSOCIATE

      if (numpar%scc) then
         write(print_to,'(a)') ' Second-order TB: including self-consistent charge calculations:'
         select case (numpar%scc_gam_ind)
         case (-1)
            text1 = ' Garrity-Choudhary'
         case (1)
            text1 = ' Klopman-Ohno'
         case (2)
            text1 = ' Mataga-Nishimoto'
         case default
            text1 = " Garrity-Choudhary with Wolf's Coulomb "
         end select
         write(print_to,'(a,a)') ' Model for gamma-function used: ', trim(adjustl(text1))
         write(text1, '(f6.2)') numpar%scc_mix
         write(print_to,'(a,a)') ' Mixing factor for SCC : ', trim(adjustl(text1))
      else
         write(print_to,'(a)') ' Zero-order TB: non-self-consistent-charge calculations'
      endif

      
      if (allocated(Scell(1)%TB_Waals)) then ! if we have vdW potential defined
         ! Find (first awailable) name of vdW parameterization:
         VDWP:do i = 1, size(Scell(1)%TB_Waals,1)
            do j = 1, size(Scell(1)%TB_Waals,2)
               text1 = trim(adjustl(Scell(1)%TB_Waals(i,j)%Param))
               if (LEN(trim(adjustl(text1))) > 0) exit VDWP ! found some name
            enddo
         enddo VDWP
         write(print_to,'(a,a)') ' van der Waals energy: ', trim(adjustl(text1))
      else !For this material vdW class is undefined
         write(print_to,'(a,a)') ' No van der Waals potential was defined'
      endif
      if (allocated(Scell(1)%TB_Coul)) then ! if we have Coulomb potential defined
         write(print_to,'(a,a)') ' Coulomb energy: ', trim(adjustl(Scell(1)%TB_Coul(1,1)%Param))
      else !For this material vdW class is undefined
         write(print_to,'(a,a)') ' No Coulomb potential was defined or unbalanced charge allowed'
      endif
      if (allocated(Scell(1)%TB_Expwall)) then ! if we have exponential wall potential defined
         FPN:do i = 1, size(Scell(1)%TB_Expwall,1)
            do j = 1, size(Scell(1)%TB_Expwall,2)
               if (LEN(trim(adjustl(Scell(1)%TB_Expwall(i,j)%Param))) > 0) then
                  write(text1, '(a)') trim(adjustl(Scell(1)%TB_Expwall(i,j)%Param))
                  exit FPN
               endif
            enddo
         enddo FPN
         write(print_to,'(a,a)') ' Short-range repulsion (exponential wall): ', trim(adjustl(text1))
      else !For this material exponential wall class is undefined
         write(print_to,'(a,a)') ' No additional short-range repulsion was defined for close interatomic distances'
      endif

      ! What kind of supercell is used:
      select case (numpar%save_files_used)
      case default ! constructed from unit cells
         write(text1, '(i10)') matter%cell_x
         write(text2, '(i10)') matter%cell_y
         write(text3, '(i10)') matter%cell_z
         write(print_to,'(a,a,a,a,a,a)') ' Super-cell size in unit-cells: ', &
            trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))

         ! Which initial velocity distribution:
         select case (numpar%ind_starting_V)
         case default
            write(print_to,'(a)') ' Maxwellian distribution of initial velocities is used'
         case (1)
            write(print_to,'(a)') ' Uniform distribution of initial velocities is used'
         case (0)
            write(print_to,'(a)') ' Delta distribution of initial velocities is used'
         end select
      case (1)  ! save files are used
         write(print_to,'(a)') ' Super-cell parameters are set in SAVE files'
         if (numpar%vel_from_file) then
            write(print_to,'(a)') ' Atomic velocities read from the file'
         else
            ! Which initial velocity distribution:
            select case (numpar%ind_starting_V)
            case default
               write(print_to,'(a)') ' Maxwellian distribution of initial velocities is used'
            case (1)
               write(print_to,'(a)') ' Uniform distribution of initial velocities is used'
            case (0)
               write(print_to,'(a)') ' Delta distribution of initial velocities is used'
            end select
         endif
      case (2)  ! path coordinate
         write(print_to,'(a)') ' Coordinate path calculations are performed, defined by PATH files'
         if (numpar%vel_from_file) then
            write(print_to,'(a)') ' Atomic velocities read from the file'
         else
            ! Which initial velocity distribution:
            select case (numpar%ind_starting_V)
            case default
               write(print_to,'(a)') ' Maxwellian distribution of initial velocities is used'
            case (1)
               write(print_to,'(a)') ' Uniform distribution of initial velocities is used'
            case (0)
               write(print_to,'(a)') ' Delta distribution of initial velocities is used'
            end select
         endif
      endselect

      write(text1, '(i10)') Scell(1)%Na
      write(print_to,'(a,a)') ' Number of atoms in the supercell: ', trim(adjustl(text1))
      if (numpar%r_periodic(1)) then		! periodic (not free surface) along X
         write(print_to,'(a,a)') ' Boundary condition along X-axis: ', 'periodic'
      else
         write(print_to,'(a,a)') ' Boundary condition along X-axis: ', 'free'
      endif
      if (numpar%r_periodic(2)) then		! periodic (not free surface) along Y
         write(print_to,'(a,a)') ' Boundary condition along Y-axis: ', 'periodic'
      else
         write(print_to,'(a,a)') ' Boundary condition along Y-axis: ', 'free'
      endif
      if (numpar%r_periodic(3)) then		! periodic (not free surface) along Z
         write(print_to,'(a,a)') ' Boundary condition along Z-axis: ', 'periodic'
      else
         write(print_to,'(a,a)') ' Boundary condition along Z-axis: ', 'free'
      endif
   endif
   write(print_to,'(a,a)') ' Electron cross sections used are from: ', trim(adjustl(numpar%At_base))
   if (matter%dens < 1e6) then   ! real density:
      write(print_to,'(a,f10.3,a)') ' Density of the material: ', matter%dens,' [g/cm^3]'
      write(print_to,'(a,es12.3,a)') ' The used atomic density (used in MC cross sections): ', matter%At_dens, ' [1/cm^3]'
   else  ! in artificial cases, the density may be wild:
      write(print_to,'(a,es25.3,a)') ' Density of the material: ', matter%dens,' [g/cm^3]'
      write(print_to,'(a,es12.3,a)') ' Atomic density (used in MC cross sections): ', matter%At_dens, ' [1/cm^3]'
   endif

   write(print_to,'(a,a)') ' EADL database used: ', trim(adjustl(numpar%EADL_file))
   write(print_to,'(a,a)') ' EPDL database used: ', trim(adjustl(numpar%EPDL_file))

   !ooooooooooooooooooooooooooooooooooooooooooooo
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '  The following numerical parameters are used:'
   write(print_to,'(a,i6)') ' Number of iterations in the MC module: ', numpar%NMC
   if (numpar%do_elastic_MC) then ! allow elastic scattering of electrons on atoms within MC module
      write(print_to,'(a)') ' Elastic high-energy-electron scattering is included in MC via Motts cross section'
   else
      write(print_to,'(a)') ' Elastic high-energy-electron scattering is excluded in MC'
   endif

#ifdef _OPENMP
   write(print_to,'(a,i6)') ' Number of threads in OPENMP: ', numpar%NOMP
#else
#ifdef MPI_USED
      write(print_to,'(a,i6)') ' Number of processes in MPI: ', numpar%MPI_param%size_of_cluster
#else
      write(print_to,'(a)') ' The code is compiled without pparallelization'
#endif
#endif

   AT_MOVE:if (numpar%do_atoms) then ! atoms are moving:
      select case (numpar%MD_algo)
      case (1)
         write(print_to,'(a)') ' MD algorithm used: Yoshida (4th order)'
      case (2)
         write(print_to,'(a)') ' MD algorithm used: Martyna-Tuckerman (4th order)'
      case default
         write(print_to,'(a)') ' MD algorithm used: velocity Verlet (2th order)'
      endselect
      if (allocated(numpar%dt_MD_reset_grid)) then
         write(print_to,'(a)') ' Time step in MD simulations is set from file: '//trim(adjustl(numpar%MD_step_grid_file))
      else
         write(print_to,'(a,f9.3,a)') ' Time step in MD simulations: ', numpar%dt,' [fs]'
      endif
      write(print_to,'(a,f9.3,a)') ' Output data are saved every: ', max(numpar%dt_save,numpar%dt),' [fs]'
      if (numpar%p_const) then	! P=const
         write(print_to,'(a)') ' Constant pressure simulation (Parrinello-Rahman scheme, NPH) '
         ! Get the supercell mass in units of total atoms mass
         temp = 0.0d0
         do j = 1, Scell(1)%Na
            temp = temp + matter%Atoms(Scell(1)%MDatoms(j)%KOA)%Ma ! total mass of all atoms in supercell
         enddo
         write(text1, '(f12.3)') matter%W_PR/temp
         write(print_to,'(a)') ' With the mass coefficient (M_box/W_PR) W_PR: '//trim(adjustl(text1))
         write(print_to,'(a,f12.3,a)') ' External pressure: ', matter%p_ext/1.0d9,' [GPa]'
      else ! V=const
         write(print_to,'(a)') ' Constant volume simulation (NVE)'
      endif
   else AT_MOVE
      write(print_to,'(a)') ' Atoms were FROZEN instead of moving in MD!'
   endif AT_MOVE

   !ooooooooooooooooooooooooooooooooooooooooooooo
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '  The schemes for electron populations used are:'

   write(print_to,'(a)') ' Scheme used for low-energy electrons relaxation: '
   select case (numpar%el_ion_scheme)
   case (0)
      write(print_to,'(a)') ' Decoupled electrons and atoms (instant electron thermalization)'
   case (1)
      write(print_to,'(a)') ' Enforced total energy conservation'
   case (2)
      write(print_to,'(a)') ' Enforced constant temperature of electrons'
   case (3)
      write(print_to,'(a)') ' True Born-Oppenheimer (constant electron populations)'
   case (4)
      if (numpar%tau_fe < numpar%dt/30.0d0) then ! it's basically instantaneous
         write(text1, '(f13.6)') 0.0e0
      elseif (numpar%tau_fe < 1e6) then
         write(text1, '(f13.6)') numpar%tau_fe
      else
         write(text1, '(es16.6)') numpar%tau_fe
      endif
      write(print_to,'(a)') ' Relaxation-time approximation for electron thermalization'
      write(print_to,'(a)') ' with the total characteristic time '//trim(adjustl(text1))//' [fs]'

      !if ((numpar%tau_fe_CB > -1.0e-7) .and. (numpar%tau_fe_VB > -1.0e-7)) then ! Partial thermalization is on:
      if (numpar%do_partial_thermal) then ! Partial thermalization is on:
         write(print_to,'(a)') ' Band-resolved relaxation is applied with characteristic times:'

         if (numpar%tau_fe_VB < numpar%dt/30.0d0) then ! it's basically instantaneous
            write(text1, '(f13.6)') 0.0e0
         elseif (numpar%tau_fe_VB < 1e6) then
            write(text1, '(f13.6)') numpar%tau_fe_VB
         else
            write(text1, '(es16.6)') numpar%tau_fe_VB
         endif
         write(print_to,'(a)') ' VB: Valence band relaxation time: '//trim(adjustl(text1))//' [fs]'

         if (numpar%tau_fe_CB < numpar%dt/30.0d0) then ! it's basically instantaneous
            write(text1, '(f13.6)') 0.0e0
         elseif (numpar%tau_fe_CB < 1e6) then
            write(text1, '(f13.6)') numpar%tau_fe_CB
         else
            write(text1, '(es16.6)') numpar%tau_fe_CB
         endif
         write(print_to,'(a)') ' CB: Conduction band relaxation time: '//trim(adjustl(text1))//' [fs]'
      else
         write(print_to,'(a)') ' No band-resolved relaxation is used'
      endif

   end select

   write(print_to,'(a)') ' Scheme used for electron-ion (electron-phonon) coupling: '
   if (numpar%NA_kind == 0) then
      write(print_to,'(a)') ' No nonadiabatic coupling'
   else
      select case (numpar%NA_kind)
      case (-1)
         write(print_to,'(a)') ' Landau full probability'
      case (2)
         write(print_to,'(a)') ' Fermi golden rule coupling'
      case (3)
         write(print_to,'(a)') ' Incomplete FGR with sin-function'
      case default
         write(print_to,'(a)') ' Dynamical nonadiabatic coupling'
      end select
      write(print_to,'(a, f10.1, a)') ' switched on at: ', numpar%t_NA, ' [fs]'
      write(print_to,'(a, f7.1, a)') ' with the acceptance window: ', numpar%acc_window, ' [eV]'
      write(print_to,'(a, f8.5, a)') ' degeneracy tolerance: ', numpar%degeneracy_eV, ' [eV]'
      write(text,'(f8.5)') numpar%M2_scaling
      if (numpar%M2_scaling == 4.0d0) text = trim(adjustl(text))//' (default)'
      write(print_to,'(a,a)') ' and scaling factor of: ', trim(adjustl(text))
      select case (numpar%ind_at_distr)
      case (1)
         write(print_to,'(a)') ' using transient nonequilibrium atomic distribution'
      case default
         write(print_to,'(a)') ' using equivalent Maxwellian atomic distribution'
      endselect
   endif

   if (numpar%do_cool) then
      write(print_to,'(a)') ' Quenching of atoms for resolidification is included'
      write(print_to,'(a, f7.1, a, f7.1, a)') ' Starting at: ', numpar%at_cool_start, ' [fs], with the step of: ', numpar%at_cool_dt, ' [fs]'
   else
      write(print_to,'(a)') ' No quenching of atoms for resolidification'
   endif

   if (allocated(numpar%El_bath_reset_grid)) then
      write(print_to,'(a)') ' Berendsen thermostat is used for electrons'
      write(print_to,'(a)') ' with parameters set in the file: '//trim(adjustl(numpar%El_bath_step_grid_file))
   elseif (g_numpar%Transport_e) then ! for electrons
      write(text,'(f10.1)') matter%T_bath_e*g_kb
      write(text1,'(f10.1)') matter%tau_bath_e
      write(print_to,'(a)') ' Berendsen thermostat is used for electrons'
      write(print_to,'(a)') ' Electronic bath temperature: '//trim(adjustl(text))//' [K], time constant: '//trim(adjustl(text1))//' [fs]'
   else
      write(print_to,'(a)') ' No electronic thermostat is used'
   endif

   if (allocated(numpar%At_bath_reset_grid)) then
      write(print_to,'(a)') ' Berendsen thermostat is used for atoms'
      write(print_to,'(a)') ' with parameters set in the file: '//trim(adjustl(numpar%At_bath_step_grid_file))
   elseif (g_numpar%Transport) then ! for atoms
      write(text,'(f10.1)') matter%T_bath*g_kb
      write(text1,'(f10.1)') matter%tau_bath
      write(print_to,'(a)') ' Berendsen thermostat is used for atoms'
      write(print_to,'(a)') ' Atomic bath temperature: '//trim(adjustl(text))//' [K], time constant: '//trim(adjustl(text1))//' [fs]'
   else
      write(print_to,'(a)') ' No atomic thermostat is used'
   endif

   write(text1, '(f7.1)') numpar%E_cut
   write(print_to,'(a, a, a)') ' Electron energy cut-off, separating high-energy- from low-energy-electrons: ', trim(adjustl(text1)), ' [eV]'
   select case (numpar%el_ion_scheme)
   case (3:4)
      write(print_to,'(a)') ' But it maybe dynamically adjusted to the top of CB (nonequilibrium simulation)'
   endselect

   if (numpar%E_work >= 1.0d25) then
      write(print_to,'(a)') ' No electron emission is allowed in the calculation'
   else if (numpar%E_work >= 0.0d0) then
      write(print_to,'(a, f7.1, a)') ' Electron is considerred to be emitted if its energy is above: ', numpar%E_work, ' [eV]'
   else ! < 0, => number of collisions is the conduction, instead of work function
      write(print_to,'(a, f2.0, a)') ' Electron is considerred to be emitted after ', ABS(numpar%E_work), ' collisions'
   endif


   !ooooooooooooooooooooooooooooooooooooooooooooo
   ! Optional output:
   write(print_to,'(a)') trim(adjustl(m_starline)) ! Output stuff:
   optional_output = .false.  ! to start with
   write(print_to,'(a)') '  Optional output:'

   if (numpar%save_testmode) then   ! testmode is on, some data for testing are printed out
      write(print_to,'(a)') ' Testmode is on, center-of-mass, rotation, total forces, etc.'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%print_MFP) then ! printout MFP file
      write(print_to,'(a)') ' Electron and photon mean free paths'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_CDF) then ! printout CDF file
      write(print_to,'(a)') ' CDF (complex dielectric function) parameters used in MC'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_Ei) then
      write(print_to,'(a)') ' Electron energy levels (molecular orbitals)'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_DOS) then
      write(print_to,'(a, f7.5, a)') ' Density of states (DOS); smearing used: ', numpar%Smear_DOS, ' [eV]'
      select case (ABS(numpar%optic_model))	! use multiple k-points, or only gamma
      case (2,4:5)   ! multiple k points
         write(text1, '(i10)') numpar%ixm
         write(text2, '(i10)') numpar%iym
         write(text3, '(i10)') numpar%izm
         if (allocated(numpar%k_grid)) then
            write(print_to,'(a,a,a,a,a,a)') ' calculated on the user-defined grid for points: ', &
               trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
         else
            write(print_to,'(a,a,a,a,a,a)') ' calculated on Monkhorst-Pack grid for points: ', &
               trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
         endif
      case default   ! gamma point
         write(print_to,'(a)') ' calculated at the Gamma point'
      end select
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%Mulliken_model >= 1) then
      write(print_to,'(a)') ' Average Mulliken charges on various elements'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_fe) then
      write(print_to,'(a)') ' Electron distribution on energy levels'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_fe_orb) then
      write(print_to,'(a)') ' Orbital-resolved electron distribution on energy levels'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_fe_grid) then
      write(print_to,'(a)') ' Electron distribution on the grid'
      optional_output = .true.   ! there is at least some optional output
   endif

    if (numpar%save_fa) then
      write(print_to,'(a)') ' Atomic distribution'
      optional_output = .true.   ! there is at least some optional output
   endif

    if (numpar%print_Ta) then
      write(print_to,'(a)') ' Various definitions of atomic temperatures'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_PCF) then
      write(print_to,'(a)') ' Atomic pair correlation function'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_XYZ) then
      write(print_to,'(a)') ' Atomic coordinates in XYZ-format'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_CIF) then
      write(print_to,'(a)') ' Atomic coordinates in CIF-format'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_raw) then
      write(print_to,'(a)') ' Atomic coordinates and velocities (raw data)'
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%save_NN) then
      write(text1, '(f6.2)') numpar%NN_radius
      write(print_to,'(a,a,a)') ' Nearest neighbors numbers within the radius of ', trim(adjustl(text1)), ' [A]'
      optional_output = .true.   ! there is at least some optional output
   endif
   if (allocated(numpar%NN_radii)) then
      write(print_to,'(a)') ' with element-specific radii of'
      do i = 1, size(numpar%NN_radii)
         write(text1, '(f6.2)') numpar%NN_radii(i)%r_cut
         write(print_to,'(a,a)') ' '//numpar%NN_radii(i)%Name//' : ', trim(adjustl(text1))//' [A]'
      enddo
      optional_output = .true.   ! there is at least some optional output
   endif

   if (numpar%do_kappa) then
      write(print_to,'(a)') ' Electronic heat conductivity (and Ce, mu) vs Te'
      write(text1, '(i10)') numpar%ixm
      write(text2, '(i10)') numpar%iym
      write(text3, '(i10)') numpar%izm
      if (allocated(numpar%k_grid)) then
         write(print_to,'(a,a,a,a,a,a)') ' Averaged over k-points (on user-defined grid): ', &
               trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
      else
         write(print_to,'(a,a,a,a,a,a)') ' Averaged over k-points (on Monkhorst-Pack grid): ', &
               trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
      endif
      optional_output = .true.   ! there is at least some optional output
   !else
   !   write(print_to,'(a)') ' No calculation of electronic heat conductivity'
   endif

   if (numpar%do_kappa_dyn) then
      write(print_to,'(a)') ' Electronic heat conductivity (dynamical)'
      optional_output = .true.   ! there is at least some optional output
   !else
   !   write(print_to,'(a)') ' No calculation of electronic heat conductivity (внтфьшсфд)'
   endif

   if (numpar%save_diff_peaks) then
      write(text1, '(f16.2)') Scell(1)%diff_peaks%hw
      write(text2, '(f16.4)') Scell(1)%diff_peaks%l*1.0d10
      write(print_to,'(a)') ' Diffraction peaks for X-ray ('//trim(adjustl(text1))//' [eV], '//trim(adjustl(text2))//' [A]):'
      do i = 1, size(Scell(1)%diff_peaks%I_diff_peak)
         write(text1, '(i0)') Scell(1)%diff_peaks%ijk_diff_peak(1,i)
         write(text2, '(i0)') Scell(1)%diff_peaks%ijk_diff_peak(2,i)
         write(text3, '(i0)') Scell(1)%diff_peaks%ijk_diff_peak(3,i)
         write(print_to,'(a)', advance='no') '('//trim(adjustl(text1))//trim(adjustl(text2))//trim(adjustl(text3))//')'//' '
      enddo
      write(print_to,'(a)') ''
      optional_output = .true.   ! there is at least some optional output
   endif

   if (.not.optional_output) then ! there ws no optional output, report it
      write(print_to,'(a)') ' none requested by the user'
   endif

   if (LEN(trim(adjustl(numpar%output_path))) > 0) then  ! output directory name already defined, print it out:
      write(print_to,'(a, a)') ' Output saved in the directory: ', trim(adjustl(numpar%output_path))
   endif

9999   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine Print_title



subroutine print_Hamiltonian_info(print_to, Scell, matter)
   integer, intent(in) :: print_to
   type(Super_cell), intent(in) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! materil parameters
   !--------------------
   integer :: i, j, Nsiz
   logical :: different_TB_param

   different_TB_param = .false. ! to start with
   Nsiz = size(Scell%TB_Hamil,1)

   if (Nsiz > 1) then  ! more than one parameterization is possible:
      ! Check if ther eis more than one parameterization, or are they the same:
      CHKTB:do i = 1, Nsiz
         do j = 1, Nsiz
            if (trim(adjustl(Scell%TB_Hamil(i,j)%Param)) /= trim(adjustl(Scell%TB_Hamil(1,1)%Param))) then
               different_TB_param = .true.
               exit CHKTB
            endif
         enddo
      enddo CHKTB
   endif

   ! If there is more then one TB parameterization:
   if (different_TB_param) then
      write(print_to,'(a)') ' Hamiltonians:      '
      do i = 1, Nsiz
         do j = i, Nsiz
            write(print_to,'(a,a)') trim(adjustl((matter%Atoms(i)%Name)))//'-'//trim(adjustl((matter%Atoms(j)%Name)))//':  ', &
                                    trim(adjustl(Scell%TB_Hamil(i,j)%Param))
         enddo
      enddo

      write(print_to,'(a)') ' Repulsive parts:   '
      do i = 1, Nsiz
         do j = i, Nsiz
            write(print_to,'(a,a)') trim(adjustl((matter%Atoms(i)%Name)))//'-'//trim(adjustl((matter%Atoms(j)%Name)))//':  ', &
                                    trim(adjustl(Scell%TB_Repuls(i,j)%Param))
         enddo
      enddo
   else ! Only one parameterization:
      write(print_to,'(a,a)') ' Hamiltonian:      ', trim(adjustl(Scell%TB_Hamil(1,1)%Param))
      write(print_to,'(a,a)') ' Repulsive part:   ', trim(adjustl(Scell%TB_Repuls(1,1)%Param))
   endif
end subroutine print_Hamiltonian_info


subroutine XTANT_label(print_to, ind)
   integer, intent(in) :: print_to ! the screen, or file number to print to
   integer, intent(in) :: ind    ! which label to print
   select case (ind)
   case default   ! regular
      call XTANT_label_3(print_to, ind)  ! below
   case (-2)   ! old regular
      call XTANT_label_old(print_to)  ! below
   case (-1)       ! none
      ! Print nothing
   case (0)       ! small
      call XTANT_label_small(print_to) ! below
   case (2)       ! shaded
      call XTANT_label_shade(print_to) ! below
   case (3)       ! filled
      call XTANT_label_filled(print_to) ! below
   case (4)       ! starred
      call XTANT_label_starred(print_to) ! below
   endselect
end subroutine XTANT_label



subroutine XTANT_label_3(print_to, ind)
   integer, intent(in) :: print_to ! the screen, or file
   integer, intent(in) :: ind    ! which label to print
   !---------------------
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '*   __    __  _______     __       __    _   _______        *'
   write(print_to,'(a)') '*   \ \  / / |__   __|   /  \     |  \  | | |__   __|       *'
   write(print_to,'(a)') '*    \ \/ /     | |     / /\ \    |   \ | |    | |  ____    *'
   write(print_to,'(a)') '*     )  (      | |    / /__\ \   | |\ \| |    | | |__  )   *'
   write(print_to,'(a)') '*    / /\ \     | |   / ______ \  | | \   |    | |   / /    *'
   write(print_to,'(a)') '*   /_/  \_\    |_|  /_/      \_\ |_|  \__|    |_|   \ \    *'
   write(print_to,'(a)') '*                                                  ___) )   *'
   write(print_to,'(a)') '*                                                 /____/    *'
   write(print_to,'(a)') '*                                                           *'
   write(print_to,'(a)') trim(adjustl(m_starline))
   select case (ind)
   case (1) ! OpenMP
   write(print_to,'(a)') '*      (Version compiled with OpenMP parallelization)       *'
   case default   ! MPI
   write(print_to,'(a)') '*       (Version compiled with MPI parallelization)         *'
   end select
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine XTANT_label_3


subroutine XTANT_label_small(print_to)
   integer, intent(in) :: print_to ! the screen, or file
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '       __    __  _______    __      __   _   _______ '
   write(print_to,'(a)') '       \ \  / / |__   __|  /  \    |  \ | | |__   __|'
   write(print_to,'(a)') '        \ \/ /     | |    / /\ \   |   \| |    | |   '
   write(print_to,'(a)') '        / /\ \     | |   / ___  \  | |\   |    | |   '
   write(print_to,'(a)') '       /_/  \_\    |_|  /_/    \_\ |_| \__|    |_| 3 '
   write(print_to,'(a)') ' '
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine XTANT_label_small

subroutine XTANT_label_old(print_to)
   integer, intent(in) :: print_to ! the screen, or file
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '     __    __  _______     __       __    _   _______ '
   write(print_to,'(a)') '     \ \  / / |__   __|   /  \     |  \  | | |__   __|'
   write(print_to,'(a)') '      \ \/ /     | |     / /\ \    |   \ | |    | |   '
   write(print_to,'(a)') '       )  (      | |    / /__\ \   | |\ \| |    | |   '
   write(print_to,'(a)') '      / /\ \     | |   / ______ \  | | \   |    | |   '
   write(print_to,'(a)') '     /_/  \_\    |_|  /_/      \_\ |_|  \__|    |_| 3 '
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine XTANT_label_old

subroutine XTANT_label_shade(print_to)
   integer, intent(in) :: print_to ! the screen, or file
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '     __    __  _______     __       __    _   _______ '
   write(print_to,'(a)') '     \ \  /// |__   _||   / \\     |  \  ||| |__   _||'
   write(print_to,'(a)') '      \ \///     | |     / /\\\    |   \ |||    | |   '
   write(print_to,'(a)') '       ) |(      | |    / /__\\\   | |\ \|||    | |   '
   write(print_to,'(a)') '      / /\\\     | |   / ______\\  | | \  ||    | |   '
   write(print_to,'(a)') '     /_/  \_\    |_|  /_/      \\\ |_|  \__|    |_| 3 '
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine XTANT_label_shade

subroutine XTANT_label_filled(print_to)
   integer, intent(in) :: print_to ! the screen, or file
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '     __    __  _______     __       __    _   _______  '
   write(print_to,'(a)') '     \*\  /*/ |__***__|   /**\     |**\  |*| |__***__| '
   write(print_to,'(a)') '      \*\/*/     |*|     /*/\*\    |***\ |*|    |*|    '
   write(print_to,'(a)') '       )**(      |*|    /*/__\*\   |*|\*\|*|    |*|    '
   write(print_to,'(a)') '      /*/\*\     |*|   /*______*\  |*| \***|    |*|    '
   write(print_to,'(a)') '     /_/  \_\    |_|  /_/      \_\ |_|  \__|    |_| *3 '
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine XTANT_label_filled

subroutine XTANT_label_starred(print_to)
   integer, intent(in) :: print_to ! the screen, or file
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') '     **   ** ********   **     ***   ** ******** '
   write(print_to,'(a)') '      ** **     **     ****    ****  **    **    '
   write(print_to,'(a)') '       ***      **    **  **   ** ** **    **    '
   write(print_to,'(a)') '      ** **     **   ********  **  ****    **    '
   write(print_to,'(a)') '     **   **    **  **      ** **   ***    **  3 '
   write(print_to,'(a)') '     (No OpenMP version, single-thread compiled) '
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine XTANT_label_starred


subroutine print_a_comforting_message(print_to, path_sep)
   integer, intent(in) :: print_to ! the screen, or file
   character(*), intent(in) :: path_sep ! file name
   !----------------------
   character(500), dimension(:), allocatable :: quote
   character(50), dimension(:), allocatable :: author
   character(500) :: file_name
   character(50) :: ch_temp
   integer :: FN, Nsiz, Reason, count_lines, i
   logical :: file_exists, file_opened, read_well, read_text_well
   real(8) :: RN
   !----------------------

   file_name = trim(adjustl(m_INPUT_directory))//path_sep//trim(adjustl(m_INFO_directory))//path_sep//trim(adjustl(m_QUOTES_file))
   write(print_to,'(a)') trim(adjustl(m_starline))
   write(print_to,'(a)') 'Sorry that something went wrong.'
   write(print_to,'(a)') 'To lighten the mood a little, here is a quote for you:'

   inquire(file=trim(adjustl(file_name)),exist=file_exists)
   if (.not.file_exists) then ! no file, cannot print help
      write(print_to,'(a)') 'Could not find file ', trim(adjustl(file_name))
      write(print_to,'(a)') 'Cannot print a quote, sorry.'
   else ! (.not.file_exists)
      FN=201
      open(UNIT=FN, FILE = trim(adjustl(file_name)), status = 'old', action='READ')
      inquire(file=trim(adjustl(file_name)),opened=file_opened)
      if (.not.file_opened) then
         write(print_to,'(a)') 'Could not open file ', trim(adjustl(file_name))
         write(print_to,'(a)') 'Cannot print a quote, sorry.'
      else ! (.not.file_opened)
         ! Allocate arrays with quotes:
         call Count_lines_in_file(FN, Nsiz)  ! module "Dealing_with_files"
         Nsiz = Nsiz/3
         allocate(quote(Nsiz))
         allocate(author(Nsiz))
         ! Now, read the quotes:
         read_text_well = .true. ! to start with
         count_lines = 0   ! to start with
         do i = 1, Nsiz
            read(FN,'(a)',IOSTAT=Reason) quote(i)
            read(FN,'(a)',IOSTAT=Reason) author(i)
            read(FN,'(a)',IOSTAT=Reason)  ! skip an empty line
            !print*, i, quote(i), author(i)
            call read_file(Reason, count_lines, read_text_well)   ! module "Dealing_with_files"
            if (Reason > 0) then   ! something wrong in the line
               write(print_to,'(a)') 'Problem reading file '//trim(adjustl(file_name))
               write(ch_temp, '(i)') count_lines
               write(print_to,'(a)') 'in line '//trim(adjustl(ch_temp))
               read_well = .false.
            elseif (Reason < 0) then ! end of file reached ...
               close(FN)
            endif
         enddo

         ! Once read all quates, now pick one at random to printout:
         call random_number(RN)
         i = ceiling(RN*Nsiz)
         write(print_to,'(A)') trim(adjustl(quote(i)))
         write(print_to,'(A)') '-- '//trim(adjustl(author(i)))

      endif ! (.not.file_opened)
   endif ! (.not.file_exists)
   write(print_to,'(a)') trim(adjustl(m_starline))
end subroutine print_a_comforting_message



END MODULE Dealing_with_output_files
