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


character(50), parameter :: m_Python_plot_all = 'OUTPUT_Python_plot_all.py'


public :: collect_python_plots, execute_all_pyplots, create_python_plot_scripts, &
            Plot_electron_MFP_python, Plot_photon_MFP_python, Plot_laser_spectrum_python


 contains




subroutine create_python_plot_scripts(Scell,matter,numpar,laser, file_path, file_temperatures, file_pressure, file_energies, &
file_atoms_R, file_atoms_S, file_supercell, file_electron_properties, file_heat_capacity, file_heat_capacity_dyn, &
file_numbers, file_orb, file_deep_holes, file_optics, file_Ei, file_PCF, file_NN, file_element_NN, file_electron_entropy, file_Te, file_mu, &
file_atomic_entropy, file_atomic_temperatures, file_atomic_temperatures_part, file_sect_displ, &
file_diffraction_peaks, file_diffraction_peaks_part, file_diffraction_powder, file_diffraction_peaks_DW, file_Debye_temperature, &
file_testmode)
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
   character(*), intent(in) :: file_diffraction_peaks, file_diffraction_powder, file_diffraction_peaks_DW  ! diffraction peaks
   character(*), dimension(:), intent(in) :: file_diffraction_peaks_part
   character(*), intent(in) :: file_Debye_temperature ! Debye temperatures
   character(*), intent(in) :: file_testmode    ! testmode data
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



   ! Energy levels:
   if (numpar%save_Ei) then
      call Python_plot_energy_levels(numpar, t0, t_last, Scell, 'OUTPUT_energy_levels.dat', 'OUTPUT_energy_levels.py')  ! below
   endif


   ! Energies:
   call Python_plot_energies(numpar, file_energies, t0, t_last, 'OUTPUT_energies.py') ! below

end subroutine create_python_plot_scripts



subroutine Python_plot_energies(numpar, file_energies, t0, t_last, script_name) ! below
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   character(*), intent(in) :: file_Energies, script_name  ! file with energy levels, script
   !----------------
   integer :: FN
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_lables
   character(300) :: File_name

   ! Py script file:
   File_name  = trim(adjustl(numpar%output_path))//trim(adjustl(numpar%path_sep))//trim(adjustl(script_name))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")

   ! Prepare column indices and names:
   allocate(col_nums(4), source = 0)
   allocate(col_lables(4))

   col_nums(1) = 3
   col_lables(1) = '"Potential energy"'
   col_nums(2) = 5
   col_lables(2) = '"Atomic energy"'
   col_nums(3) = 6
   col_lables(3) = '"Atoms and electrons"'
   col_nums(4) = 7
   col_lables(4) = '"Total energy"'

   call Create_python_plot(FN, file_energies, col_nums, col_lables, &
      'Time (fs)', 'Energy (eV/atom)', 'Energies', &
      "best", 'OUTPUT_energies', trim(adjustl(numpar%fig_extention)), &
      x_min=t0, x_max=t_last)     ! below

   close(FN)
end subroutine Python_plot_energies




subroutine Python_plot_energy_levels(numpar, t0, t_last, Scell, file_Ei, script_name)
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: t0, t_last      ! starting and ending time
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   character(*), intent(in) :: file_Ei, script_name  ! file with energy levels, script
   !----------------
   integer i, M, NSC, FN
   integer, dimension(:), allocatable :: col_nums
   character(30), dimension(:), allocatable :: col_lables
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
      call Create_python_plot(FN, file_Ei, col_nums, col_lables, &
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
   character(30), dimension(:), allocatable :: col_lables, col_labels2
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
   allocate(col_lables(N_tot_col+1))
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
            col_lables(count_col) = '"'//trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))//'"'
         endif
      enddo
   enddo
   ! Add the valence band MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_lables(count_col) = '"Valence"'
   ! And the total MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_lables(count_col) = '"Total inelastic"'

   ! Add the elastic part:
   allocate(col_nums2(1), source = 0)
   allocate(col_labels2(1))
   col_nums2(1) = 1+size(matter%Atoms)
   col_labels2(1) = '"Elastic"'

   ! Use them to plot the data:
   call Create_python_plot(FN, file_electron_IMFP, col_nums, col_lables, &
      'Electron energy (eV)', 'Mean free path (A)', 'Electron mean free paths', &
      "best", 'OUTPUT_MFP_electron', trim(adjustl(numpar%fig_extention)), &
      x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, &
      set_x_log=.true., set_y_log=.true., &
      Data_file2=file_electron_EMFP, col_nums2=col_nums2, col_labels2=col_labels2, &
      colors_inverted=.true. )     ! below

   close(FN)
   deallocate(col_nums, col_lables, col_nums2, col_labels2)
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
   character(30), dimension(:), allocatable :: col_lables
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
   allocate(col_lables(N_tot_col+1))
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
            col_lables(count_col) = '"'//trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(matter%Atoms(i)%Shell_name(j)))//'"'
         endif
      enddo
   enddo
   ! Add the valence band MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_lables(count_col) = '"Valence"'
   ! And the total MFP:
   count_col = count_col + 1  ! number of columns
   col_nums(count_col) = count_col
   col_lables(count_col) = '"Total"'

   ! Use them to plot the data:
   call Create_python_plot(FN, file_photon_MFP, col_nums, col_lables, &
      'Photon energy (eV)', 'Attenuation length (A)', 'Photon attenuation length', &
      "best", 'OUTPUT_MFP_photon', trim(adjustl(numpar%fig_extention)), &
      x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, &
      set_x_log=.true., set_y_log=.true., colors_inverted=.true. )     ! below

   close(FN)
   deallocate(col_nums, col_lables)
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
   character(30), dimension(:), allocatable :: col_lables, linestyle
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
   allocate(col_lables(3))
   allocate(linestyle(3))
   ! Fill the column numbers and labels:
   col_nums(1) = 1
   col_nums(2) = 2
   col_nums(3) = 3
   ! Names of the columns (skipping the "Energy"):
   col_lables(1) = '"Incoming"'
   col_lables(2) = '"Absorbed"'
   col_lables(3) = '"MC Sampled"'
   ! Style of curves:
   linestyle(1)='"--"'
   linestyle(2)='"-"'
   linestyle(3)='"-."'

   ! Use them to plot the data:
   call Create_python_plot(FN, file_spectrum, col_nums, col_lables, &
      'Photon energy (eV)', 'Photon spectrum (arb. units)', 'Photon spectrum', &
      "best", 'OUTPUT_photon_spectrum', trim(adjustl(numpar%fig_extention)), &
      x_min=x_min, x_max=x_max, y_min=y_min, &
      l_style=linestyle, &
      colors_inverted=.false. )     ! below

   close(FN)
   deallocate(col_nums, col_lables, linestyle)
end subroutine Plot_laser_spectrum_python







!===================================================================
! General routines to use python for plotting:
subroutine Create_python_plot(FN, Data_file, col_nums, col_lables, &
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
   character(*), dimension(:), allocatable, intent(in) :: col_lables    ! array of column labels
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
   character(10000) :: col_nums_txt, col_lables_txt, col_nums_txt2, col_lables_txt2, linestyle_txt
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
   write(FN,'(a)') 'df = pd.read_csv(r"'//trim(adjustl(Data_file))//'", '//"sep=r'\s+', "//'header=None, comment="#")'

   ! Log-scale, if required:
   if (present(set_x_log)) then
      if (set_x_log) write(FN,'(a)') 'plt.xscale("log")'
   endif
   if (present(set_y_log)) then
      if (set_y_log) write(FN,'(a)')'plt.yscale("log")'
   endif

   ! Set a list of which columns to plot:
   col_nums_txt = '['
   do i = 1, size(col_nums)
      ! Column numbers:
      if (i > 1) then   ! add come in between
         col_nums_txt = trim(adjustl(col_nums_txt))//','
      endif
      write(temp_txt,'(i)') col_nums(i)
      col_nums_txt = trim(adjustl(col_nums_txt))//' '//trim(adjustl(temp_txt))
   enddo
   col_nums_txt = trim(adjustl(col_nums_txt))//']'

   ! Column lables, if required:
   if (allocated(col_lables)) then ! the legend is required:
      col_lables_txt = '['
      do i = 1, size(col_lables)
         ! Column titles:
         if (i > 1) then   ! add come in between
            col_lables_txt = trim(adjustl(col_lables_txt))//','
         endif
         col_lables_txt = trim(adjustl(col_lables_txt))//' '//trim(adjustl(col_lables(i)))
      enddo
      col_lables_txt = trim(adjustl(col_lables_txt))//']'
   endif

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

   ! Prepare the plot:
   write(FN,'(a)') '# Prepare the plot:'
   write(FN,'(a)') 'columns_to_plot = '//trim(adjustl(col_nums_txt))
   if (allocated(col_lables)) then ! the legend is required
      write(FN,'(a)') 'labels = '//trim(adjustl(col_lables_txt))
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
   write(FN,'(a)') 'linestyles = '//trim(adjustl(linestyle_txt))

   ! Make the plot:
   write(FN,'(a)') 'for i, col in enumerate(columns_to_plot):'          ! cycle for all curves
   write(FN,'(a)') '    color = colors[i % len(colors)]'                ! repeat colors
   write(FN,'(a)') '    ls    = linestyles[i % len(linestyles)]'        ! repeat line styles
   if (allocated(col_lables)) then
      write(FN,'(a)') '    label = labels[i % len(labels)]'                ! repeat colors
   endif
   write(FN,'(a)') '    plt.plot(df.iloc[:, 0], df.iloc[:, col],'       ! plot columns
   write(FN,'(a)') '    color=color,'                                   ! set color
   if (allocated(col_lables)) then                                      ! the legend is required
      write(FN,'(a)') '    label=label,'
   endif
   write(FN,'(a)') '    linestyle=ls)'                                  ! line style


   !-----------------------
   ! If we want to add data from another file on the same plot:
   if (present(Data_file2) .and. present(col_nums2) .and. present(col_labels2)) then !# Add a curve from the second file
      write(FN,'(a)') '# Add a curve from the second file'
      write(FN,'(a)') 'df2 = pd.read_csv(r"'//trim(adjustl(Data_file2))//'", '//"sep=r'\s+', "//'header=None, comment="#")'

      ! Set a list of which columns to plot:
      col_nums_txt2 = '['
      do i = 1, size(col_nums2)
         ! Column numbers:
         if (i > 1) then   ! add come in between
            col_nums_txt2 = trim(adjustl(col_nums_txt2))//','
         endif
         write(temp_txt,'(i)') col_nums2(i)
         col_nums_txt2 = trim(adjustl(col_nums_txt2))//' '//trim(adjustl(temp_txt))
      enddo
      col_nums_txt2 = trim(adjustl(col_nums_txt2))//']'
      ! Column lables, if required:
      if (allocated(col_labels2)) then ! the legend is required:
         col_lables_txt2 = '['
         do i = 1, size(col_labels2)
            ! Column numbers:
            if (i > 1) then   ! add come in between
               col_lables_txt2 = trim(adjustl(col_lables_txt2))//','
            endif
            col_lables_txt2 = trim(adjustl(col_lables_txt2))//' '//col_labels2(i)
         enddo
         col_lables_txt2 = trim(adjustl(col_lables_txt2))//']'
      endif

      write(FN,'(a)') 'columns_to_plot2 = '//trim(adjustl(col_nums_txt2))
      write(FN,'(a)') 'labels2 = '//trim(adjustl(col_lables_txt2))
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
   if (allocated(col_lables)) then ! the legend is required
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
      print *, 'Could not open ',trim(adjustl(File_name)),' for py-plotting.', ' Unit = ', FN
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

   ! Collect all the py-names into this file:
   do i = 1,N_f
      if (trim(adjustl(All_files(i))) /= trim(adjustl(Py_plot_all_files))) then ! exclude the file itself
         if (trim(adjustl(path_sep)) == '\') then	! if it is Windows
            write(FN,'(a)') 'subprocess.run(["python", "'//trim(adjustl(All_files(i)))//'"])'
         else ! it is linux
            write(FN,'(a)') 'subprocess.run(["python3", "'//trim(adjustl(All_files(i)))//'"])'
         endif
      endif
   enddo
   close (FN)

   !--------------
   if (skip_exec) return   ! If execution of py-plot scripts is not requested, we are done;
   ! otherwise, execute all the gnuplot scripts, if requested:
   idir = chdir(trim(adjustl(output_path))) ! go into the directory with output files

   if (trim(adjustl(path_sep)) == '\') then	! if it is Windows
      iret = system('python '//trim(adjustl(Py_plot_all_files)))   ! create the folder
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
