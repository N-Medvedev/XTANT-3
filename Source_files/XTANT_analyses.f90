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
! This module contains XTANT phase-path and cohesive analyses

MODULE XTANT_analyses

use Universal_constants
use Objects
use Little_subroutines, only : print_time_step, parse_yes_no
use Algebra_tools, only : DET_3X3
use Atomic_tools, only : Coordinates_rel_to_abs, velocities_abs_to_rel, shortest_distance, get_diffraction_peaks, &
                         get_mean_square_displacement, Coordinates_abs_to_rel
use TB, only : get_Hamilonian_and_E, vdW_interplane, get_electronic_thermal_parameters, get_DOS, get_Mullikens_all
use Electron_tools, only : Electron_thermalization, get_glob_energy
use Dealing_with_output_files, only : write_output_files, write_energies
use Dealing_with_files, only : close_file
use Optical_parameters, only : get_optical_parameters
use ZBL_potential, only : get_total_ZBL
use Atomic_thermodynamics, only : get_atomic_distribution, update_Ta_config_running_average

implicit none
PRIVATE

! Modular parameters:
character(50), parameter :: m_Coordinate_path_energy = 'OUTPUT_coordinate_path.dat'
character(50), parameter :: m_Cohesive_energy = 'OUTPUT_Energy.dat'
character(50), parameter :: m_Cohesive_split_target_energy = 'OUTPUT_split_target_energy.dat'


public :: coordinate_path, vary_size

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 contains



!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! Use this for obtaining coordinate path between two phases:
subroutine coordinate_path(numpar, matter, Scell, cur_time, Err) ! THIS SUBROUTINE USES GLOBAL VARIABLES
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter ! material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err	! error save
   !------------------------------------------
!    integer, intent(in) :: ind ! 0=NVE, 1=NPH
   integer :: i, i_step, i_at, Nat, N_steps, SCN
   type(Atom), dimension(:), allocatable :: MDAtoms ! if more then one supercell
   type(Atom), dimension(:), allocatable :: MDAtoms0 ! if more then one supercell
   real(8), dimension(3,3) :: supce, supce0 	! [A] length of super-cell
   real(8), dimension(3,3) :: Vsupce, Vsupce0    ! Derivatives of Super-cell vectors (velosities)
   real(8) :: sc_fact, cur_time
   
   if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      write(6, '(a)') 'Starting subroutine coordinate_path ...'
      open(UNIT = 100, FILE = trim(adjustl(m_Coordinate_path_energy))) !<-
   endif
   Nat = size(Scell(1)%MDatoms)   ! number of atoms in the supercell
   allocate(MDAtoms(Nat))
   allocate(MDAtoms0(Nat))
   SCN = 1
   do i = 1, Nat ! to use below
      MDAtoms(i)%S0(:) = Scell(SCN)%MDAtoms(i)%S0(:)
      MDAtoms(i)%SV0(:) = Scell(SCN)%MDAtoms(i)%SV0(:)
      MDAtoms(i)%S(:) = Scell(SCN)%MDAtoms(i)%S(:)
      MDAtoms(i)%SV(:) = Scell(SCN)%MDAtoms(i)%SV(:)
      ! Take care of boundary crossing:
      if ( abs(MDAtoms(i)%S(1) - MDAtoms(i)%S0(1)) > 0.5 ) then
         if (MDAtoms(i)%S(1) > MDAtoms(i)%S0(1)) then
            MDAtoms(i)%S0(1) = MDAtoms(i)%S0(1) + 1.0d0
         else
            MDAtoms(i)%S(1) = MDAtoms(i)%S(1) + 1.0d0
         endif
      endif
      
      if ( abs(MDAtoms(i)%S(2) - MDAtoms(i)%S0(2)) > 0.5 ) then
         if (MDAtoms(i)%S(2) > MDAtoms(i)%S0(2)) then
            MDAtoms(i)%S0(2) = MDAtoms(i)%S0(2) + 1.0d0
         else
            MDAtoms(i)%S(2) = MDAtoms(i)%S(2) + 1.0d0
         endif
      endif
      
      if ( abs(MDAtoms(i)%S(3) - MDAtoms(i)%S0(3)) > 0.5 ) then
         if (MDAtoms(i)%S(3) > MDAtoms(i)%S0(3)) then
            MDAtoms(i)%S0(3) = MDAtoms(i)%S0(3) + 1.0d0
         else
            MDAtoms(i)%S(3) = MDAtoms(i)%S(3) + 1.0d0
         endif
      endif
      
!       write(6,'(i3,f,f,f,f,f,f)') i, Scell(SCN)%MDAtoms(i)%S0(:), Scell(SCN)%MDAtoms(i)%S(:)
      !write(6,'(i3,f,f,f,f,f,f)') i, MDAtoms(i)%S0(:), MDAtoms(i)%S(:) 
   enddo
   supce0 = Scell(1)%supce0
   supce = Scell(1)%supce
   Vsupce0 = Scell(1)%Vsupce0
   Vsupce = Scell(1)%Vsupce
   
   N_steps = 100
   
   do i_step = 1, N_steps+1
      i = i_step
      sc_fact = dble(i_step-1)/dble(N_steps)
      cur_time = sc_fact
!       write(6,'(a,f)') 'Step:', cur_time
      
      ! set coordinates and supercell:
      Scell(1)%supce = supce0 - (supce0 - supce) * sc_fact
      Scell(SCN)%Vsupce = Vsupce0 - (Vsupce0 - Vsupce) * sc_fact
      do i_at = 1, Nat
            Scell(SCN)%MDAtoms(i_at)%S0(:) = MDAtoms(i_at)%S0(:) + (MDAtoms(i_at)%S(:) - MDAtoms(i_at)%S0(:)) * sc_fact
            Scell(SCN)%MDAtoms(i_at)%SV0(:) = MDAtoms(i_at)%SV0(:) + (MDAtoms(i_at)%SV0(:) - MDAtoms(i_at)%SV0(:)) * sc_fact
            Scell(SCN)%MDAtoms(i_at)%S(:) = MDAtoms(i_at)%S0(:) + (MDAtoms(i_at)%S(:) - MDAtoms(i_at)%S0(:)) * sc_fact
            Scell(SCN)%MDAtoms(i_at)%SV(:) = MDAtoms(i_at)%SV0(:) + (MDAtoms(i_at)%SV(:) - MDAtoms(i_at)%SV0(:)) * sc_fact
      enddo
      call Coordinates_rel_to_abs(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"
      call velocities_abs_to_rel(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"
      
      ! Contruct TB Hamiltonian, diagonalize to get energy levels, get forces for atoms and supercell:
      call get_Hamilonian_and_E(Scell, numpar, matter, 1, Err, cur_time) ! module "TB"
      ! Thermalization step for low-energy electrons (used only in relaxation-time approximation):
      call Electron_thermalization(Scell, numpar, skip_thermalization=.true.) ! module "Electron_tools"

      ! Get global energy of the system at the beginning:
      call get_glob_energy(Scell, matter) ! module "Electron_tools"

!        write(100,'(es25.16,es25.16,es25.16,es25.16)') cur_time, Scell(1)%nrg%Total+Scell(1)%nrg%E_supce+Scell(1)%nrg%El_high+Scell(1)%nrg%Eh_tot+Scell(1)%nrg%E_vdW, Scell(1)%nrg%E_rep, Scell(1)%nrg%El_low
       
       if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         call write_energies(6, cur_time, Scell(1)%nrg)   ! module "Dealing_with_output_files"
         call write_energies(100, cur_time, Scell(1)%nrg)   ! module "Dealing_with_output_files"
       endif
       call get_electronic_thermal_parameters(numpar, Scell, 1, matter, Err) ! module "TB"

       ! Save initial step in output:
       call write_output_files(numpar, cur_time, matter, Scell) ! module "Dealing_with_output_files"
       
       call print_time_step('Coordinate path point:', cur_time, msec=.true., MPI_param=numpar%MPI_param)   ! module "Little_subroutines"
       
   enddo
   
   if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      !close(100)
      call close_file("Close", FN=100)     ! module "Dealing_with_files"
      write(6, '(a)') 'Subroutine coordinate_path completed, file '// trim(adjustl(m_Coordinate_path_energy)) //' is created'
      write(6, '(a)') 'XTANT is terminating now...'
   endif
   !Err%Err = .true.   ! not to continue with the real calculations
   Err%Stopsignal = .true.  ! not to continue with the real calculations
end subroutine coordinate_path



!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! Use this for testing and finding potential energy minimum as a function of supercell size:
subroutine vary_size(numpar, matter, Scell, cur_time, do_forces, Err)   !  THIS SUBROUTINE USES GLOBAL VARIABLES
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter ! material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   real(8), intent(inout) :: cur_time  ! [fs] current time
   integer, intent(in), optional  :: do_forces
   type(Error_handling), intent(inout), optional :: Err	! error save
   !-----------------------------------------
   real(8) :: r_sh, x, y, z, E_vdW_interplane, cur_time_save, z_sh, z_sh0, temp, E_ZBL
   real(8) :: d_i, i_min, i_max, rescale_factor
   integer i, j, at1, at2, N_points, i_test
   character(13) :: char1
   logical yesno

   if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      open(UNIT = 100, FILE = trim(adjustl(m_Cohesive_energy)))
      if (present(do_forces)) then
         write(100,'(a)') '#Distance   E_total  E_rep El_low   F_rep F_att'
      else
         write(100,'(a)') '#Distance   E_total  E_rep El_low   E_vdW E_ZBL Z_size'
      endif
   endif
   cur_time_save = cur_time
   z_sh0 = 0.0d0
   
   !----------------------------------------------
!    ! Project-specific (Graphene on SiC), removing graphene from substrate:
!    z_sh = -0.05d0 + dble(i_test)/5000.0d0
!    temp = maxval(Scell(1)%MDatoms(:)%S(3))
!    do i = 1,size(Scell(1)%MDatoms) ! find the nearest neighbour
!       if (Scell(1)%MDatoms(i)%S(3) == temp) then
!          call Shift_all_atoms(matter, Scell, 1, shz=z_sh0-z_sh, N_start=i, N_end=i) ! above
!          print*, 'ATOM #', i, Scell(1)%MDatoms(i)%S(:)
!       endif
!    enddo
!    !----------------------------------------------
   
   ! Set grid points:
   i_min = numpar%change_size_min
   i_max = numpar%change_size_max
   N_points = numpar%change_size_step-1
   N_points = max(1,N_points) ! to make sure it is not smaller than 1
   d_i = (i_max - i_min)/dble(N_points)   ! step size in units of Supce

   !print*, 'vary_size:', i_min, i_max, N_points, d_i

   ! Calculate the mean square displacement of all atoms:
   call get_mean_square_displacement(Scell, matter, Scell(1)%MSD,  Scell(1)%MSDP, numpar%MSD_power, numpar) ! module "Atomic_tools"
   if (numpar%verbose) call print_time_step('Mean displacement calculated succesfully:', msec=.true., MPI_param=numpar%MPI_param)

   ! Calculate diffraction peaks:
   call get_diffraction_peaks(Scell, matter, numpar)  ! module "Atomic_tools"

   !do i_test = 1,300 !<-
   do i_test = 1, numpar%change_size_step+1
      !----------------------------------------------
      ! General feature, changing size:
      
       !Scell(1)%supce = Scell(1)%supce0*(0.7d0 + dble(i_test)/200.0d0) !<-
       rescale_factor = i_min  +  d_i*dble(i_test-1) !<-
       Scell(1)%supce = Scell(1)%supce0 * rescale_factor   !<-
       !print*, Scell(1)%supce0(1,1)*(0.7d0 + dble(i_test)/200.0d0), Scell(1)%supce(1,1)
      
!       print*, 'Scell0', Scell(1)%supce0
!       print*, 'Scell', Scell(1)%supce
!       pause 'CELL'
      
      !----------------------------------------------
      ! Project-specific (C60 crystal), shifting one C60 ball relative to the other:
!        z_sh = -0.05d0 + dble(i_test)/5000.0d0
!        call Shift_all_atoms(matter, Scell, 1, shz=z_sh0-z_sh, N_start=61, N_end=120) ! above
!        print*, 'Z=', z_sh, Scell(1)%MDatoms(1)%S(3), Scell(1)%MDatoms(31)%S(3), Scell(1)%MDatoms(91)%S(3)
!        z_sh0 = z_sh
      !----------------------------------------------
      ! Project-specific (Graphene on SiC), removing graphene from substrate:
!       z_sh = -0.02d0 + dble(i_test)/5000.0d0
!       temp = maxval(Scell(1)%MDatoms(:)%S(3))
!       do i = 1,size(Scell(1)%MDatoms) ! find the nearest neighbour
!          if (Scell(1)%MDatoms(i)%S(3) == temp) then
!             call Shift_all_atoms(matter, Scell, 1, shz=z_sh0-z_sh, N_start=i, N_end=i) ! above
!             print*, 'ATOM #', i, Scell(1)%MDatoms(i)%S(:)
!          endif
!       enddo
!       z_sh0 = z_sh
!       !----------------------------------------------

      call Det_3x3(Scell(1)%supce,Scell(1)%V) !<- modlue "Algebra_tools"

      call Coordinates_rel_to_abs(Scell, 1, if_old=.true.)	! from the module "Atomic_tools"!<-
      
      cur_time = 1d9   ! to start with
      r_sh = 1d10    ! to start with
      at1 = 1  ! to start with
      at2 = 2  ! to start with
      do j = 1,size(Scell(1)%MDatoms)-1 ! find the nearest neighbour
         do i = j+1,size(Scell(1)%MDatoms) ! find the nearest neighbour
            call shortest_distance(Scell, 1, Scell(1)%MDatoms, j, i, r_sh) ! module 'Atomic_tools'
            if (cur_time > r_sh) then
               cur_time = r_sh ! [A] nearest neighbor distance
               at1 = j
               at2 = i
            endif
         enddo
      enddo
      !call change_r_cut_TB_Hamiltonian(1.70d0*(Scell(1)%supce(3,3)*0.25d0)/1.3d0, TB_Waals=Scell(1)%TB_Waals) !<-

      ! Contruct TB Hamiltonian, diagonalize to get energy levels, get forces for atoms and supercell:
      call get_Hamilonian_and_E(Scell, numpar, matter, 1, Err, cur_time) ! module "TB"
      if (numpar%verbose) call print_time_step('Hamiltonian constructed and diagonalized', msec=.true., MPI_param=numpar%MPI_param)

      ! Thermalization step for low-energy electrons (used only in relaxation-time approximation):
      call Electron_thermalization(Scell, numpar, skip_thermalization=.true.) ! module "Electron_tools"

      ! Get global energy of the system at the beginning:
      call get_glob_energy(Scell, matter) ! module "Electron_tools"

      ! Get initial optical coefficients:
      call get_optical_parameters(numpar, matter, Scell, Err) ! module "Optical_parameters"
      
      ! Get initial DOS:
      call get_DOS(numpar, matter, Scell, Err)	! module "TB"

      call get_Mullikens_all(Scell(1), matter, numpar)   ! module "TB"
      call get_electronic_thermal_parameters(numpar, Scell, 1, matter, Err) ! module "TB"

      ! Get atomic distribution:
      ! Update configurational temperature for running average (needed on each timestep):
      call update_Ta_config_running_average(Scell(1), matter, numpar)   ! module "Atomic_thermodynamics"
      call get_atomic_distribution(numpar, Scell, 1, matter)   ! module "Atomic_thermodynamics"

      ! Save initial step in output:
      call write_output_files(numpar, cur_time, matter, Scell) ! module "Dealing_with_output_files"

      ! Get interplane energy for vdW potential:
      E_vdW_interplane = vdW_interplane(Scell(1)%TB_Waals, Scell, 1, numpar, matter)/dble(Scell(1)%Na) !module "TB"

      ! Get ZBL potential is requested:
      call get_total_ZBL(Scell, 1, matter, numpar, E_ZBL) ! module "ZBL_potential"
      E_ZBL = E_ZBL/dble(Scell(1)%Na)   ! [eV] => [eV/atom]

      if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
         if (present(do_forces)) then
            write(*,'(a,X1i0,a,X1i0,f14.6,f14.6,f14.6)') 'Supercell size:', i_test-1, &
            ' '//trim(adjustl(matter%Atoms(Scell(1)%MDAtoms(at1)%KOA)%Name))//'-'// &
            trim(adjustl(matter%Atoms(Scell(1)%MDAtoms(at2)%KOA)%Name)) , rescale_factor, cur_time, &
            Scell(1)%nrg%Total+Scell(1)%nrg%E_supce+Scell(1)%nrg%El_high+Scell(1)%nrg%Eh_tot+Scell(1)%nrg%E_vdW
            write(100,'(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') &
               cur_time, Scell(1)%nrg%Total+Scell(1)%nrg%E_supce+Scell(1)%nrg%El_high+Scell(1)%nrg%Eh_tot+Scell(1)%nrg%E_vdW, &
               Scell(1)%nrg%E_rep, Scell(1)%nrg%El_low, Scell(1)%MDatoms(do_forces)%forces%rep(:), &
               Scell(1)%MDatoms(do_forces)%forces%att(:)
         else
            write(*,'(a,X1i0,a,X1i0,X1i0,f14.6,f14.6,f14.6)') 'Supercell size:', i_test-1, &
            ' '//trim(adjustl(matter%Atoms(Scell(1)%MDAtoms(at1)%KOA)%Name))//'-'// &
            trim(adjustl(matter%Atoms(Scell(1)%MDAtoms(at2)%KOA)%Name)), &
            at1, at2, rescale_factor, cur_time, &
            Scell(1)%nrg%Total+Scell(1)%nrg%E_supce+Scell(1)%nrg%El_high+Scell(1)%nrg%Eh_tot+Scell(1)%nrg%E_vdW
            write(100,'(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') cur_time, &
               Scell(1)%nrg%Total+Scell(1)%nrg%E_supce+Scell(1)%nrg%El_high+Scell(1)%nrg%Eh_tot+Scell(1)%nrg%E_vdW, &
               Scell(1)%nrg%E_rep, Scell(1)%nrg%El_low, E_vdW_interplane, E_ZBL, Scell(1)%supce(3,3)
         endif
      endif ! (numpar%MPI_param%process_rank == 0)
   enddo
   cur_time = cur_time_save
   Scell(1)%supce = Scell(1)%supce0

   ! Uncomment here if you want to be able to proceed with regular calculations after "size",
   ! this option has never been used, so now by default it is depricated.
!    write(*,'(a)') '*************************************************************'
!    print*, ' Would you like to proceed with XTANT calculation? (y/n)',char(13)
!    read(*,*) char1
   if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      call close_file("Close", FN=100)     ! module "Dealing_with_files"
      write(*,'(a)') '*************************************************************'
   endif
   char1 = 'n' ! by default, stop calculations here
   call parse_yes_no(trim(adjustl(char1)), yesno) ! module "Little_subroutines"
   Err%Stopsignal = .not.yesno
end subroutine vary_size


END MODULE XTANT_analyses
