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
! This module contains subroutines to read input files:

MODULE Dealing_with_output_files
! Open_MP related modules from external libraries:
#ifdef OMP_inside
   USE OMP_LIB, only : OMP_GET_MAX_THREADS
#endif
USE IFLPORT, only : system, chdir

use Universal_constants
use Objects
use Atomic_tools, only : pair_correlation_function
use Variables, only : g_numpar, g_matter
use Little_subroutines, only : number_of_types_of_orbitals, name_of_orbitals, set_starting_time, order_of_time, convolution
use Dealing_with_files, only : get_file_stat, copy_file, read_file
!use Electron_tools
use Dealing_with_EADL, only : define_PQN
use Gnuplotting
use Read_input_data, only : m_INPUT_directory, m_INFO_directory, m_INFO_file, m_HELP_file

implicit none
PRIVATE

public :: write_output_files, convolve_output, reset_dt, print_title, prepare_output_files, communicate
public :: close_save_files, close_output_files, save_duration, execute_all_gnuplots, write_energies

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
      call write_pressure(numpar%FN_pressure, time, Pressure, Stress)
      call write_energies(numpar%FN_energies, time, nrg)
      call write_numbers(numpar%FN_numbers, time, Scell(NSC))
      call write_holes(numpar%FN_deep_holes, time, matter, Scell(NSC))
      if (numpar%save_raw) call write_atomic_relatives(numpar%FN_atoms_S, Scell(NSC)%MDatoms)
      call write_super_cell(numpar%FN_supercell, time, Scell(NSC))
      call write_electron_properties(numpar%FN_electron_properties, time, Scell, NSC, Scell(NSC)%Ei, matter, numpar, &
               numpar%FN_Ce, numpar%FN_kappa, numpar%FN_Se)
      if (numpar%save_XYZ) call write_atomic_xyz(numpar%FN_atoms_R, Scell(1)%MDatoms, matter, Scell(1)%supce(:,:))
      if (numpar%save_CIF) call write_atomic_cif(numpar%FN_cif, Scell(1)%supce(:,:), Scell(1)%MDatoms, matter, time)
      if (numpar%save_Ei) then
         if (numpar%scc) then ! Energy levels include SCC term:
            call save_energy_levels(numpar%FN_Ei, time, Scell(1)%Ei_scc_part)
         else  ! non-SCC (uncorrected energy levels):
            call save_energy_levels(numpar%FN_Ei, time, Scell(1)%Ei)
         endif
      endif
      if (numpar%save_DOS) then
         select case (numpar%DOS_splitting)
         case (1) ! with partial DOS
            call save_DOS(numpar%FN_DOS, time, Scell(1)%DOS, Scell(1)%partial_DOS)
         case default   ! no partial dos
            call save_DOS(numpar%FN_DOS, time, Scell(1)%DOS)
         end select
      endif
      select case (numpar%DOS_splitting)
         case (1) ! with partial DOS
         call write_coulping(numpar%FN_coupling, time, Scell, NSC, numpar)
      end select
      if (numpar%save_fe) call save_distribution(numpar%FN_fe, time, Scell(1)%Ei, Scell(1)%fe, Scell(1)%fe_eq)
      if (numpar%save_fe_grid) call electronic_distribution_on_grid(Scell(1), numpar, time)  ! below
      if (numpar%save_PCF) call write_PCF(numpar%FN_PCF, Scell(1)%MDatoms, matter, Scell, 1)
      if (numpar%do_drude) call write_optical_coefs(numpar%FN_optics, time, Scell(1)%eps)
      if (Scell(1)%eps%all_w) call write_optical_all_hw(numpar%FN_all_w, time, Scell(1)%eps)
      if (numpar%save_NN) call save_nearest_neighbors(numpar%FN_neighbors, Scell, 1, time)
      
   enddo
end subroutine write_output_files



subroutine electronic_distribution_on_grid(Scell, numpar, tim)
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar ! numerical parameters, including MC energy cut-off
   real(8), intent(in) :: tim ! [fs] timestep
   !----------------
   integer :: i, j, Nsiz, Nei
   real(8) :: N_steps

   ! Add the high-energy part of the distribution (obtained from MC):
   Scell%fe_on_grid = Scell%fe_on_grid + Scell%fe_high_on_grid

   ! Average over the number of time-steps it was collected over:
   N_steps = max( 1.0d0, dble(numpar%fe_aver_num) )   ! at least one step, to not change anything if nothing happened
   Scell%fe_on_grid = Scell%fe_on_grid / N_steps

   ! Now save the distribution in the file:
   call save_distribution_on_grid(numpar%FN_fe_on_grid, tim, Scell%E_fe_grid, Scell%fe_on_grid)  ! below

   ! Reset the high-energy electron part for the next step:
   Scell%fe_on_grid = 0.0d0
   Scell%fe_high_on_grid = 0.0d0
   numpar%fe_aver_num = 0   ! to restart counting time-steps
end subroutine electronic_distribution_on_grid


subroutine save_distribution_on_grid(FN, tim, wr, fe)
   integer, intent(in) :: FN
   real(8), intent(in) :: tim
   real(8), dimension(:), intent(in) :: wr
   real(8), dimension(:), intent(in) :: fe
   integer i
   write(FN,'(a,f25.16)') '#', tim
   do i = 1, size(fe)
      write(FN,'(f25.16,es25.16)') wr(i), fe(i)
   enddo
   write(FN,*) ''
   write(FN,*) ''
end subroutine save_distribution_on_grid




subroutine convolve_output(Scell, numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   integer i, FN
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
      write(FN, '(f25.16, es25.16, es25.16, es25.16, f25.16, f25.16, f25.16, f25.16, f25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16, es25.16)') eps%Eps_hw(1,i), eps%Eps_hw(2,i), eps%Eps_hw(3,i), eps%Eps_hw(4,i), eps%Eps_hw(5,i), eps%Eps_hw(6,i), eps%Eps_hw(7,i), eps%Eps_hw(8,i), eps%Eps_hw(9,i), eps%Eps_hw(10,i), eps%Eps_hw(11:,i)
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
   write(FN, '(f25.16,es25.16,es25.16,es25.16,f25.16,f25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') tim, eps%R, eps%T, eps%A, eps%n, eps%k, eps%ReEps, eps%ImEps, eps%dc_cond, eps%Eps_xx, eps%Eps_yy, eps%Eps_zz
end subroutine write_optical_coefs


subroutine write_PCF(FN, atoms, matter, Scell, NSC)
   integer, intent(in) :: FN	! file number
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   type(Solid), intent(inout) :: matter	! Material parameters
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   integer i
   call pair_correlation_function(atoms, matter, Scell, NSC) ! module "Atomic_tools"
   do i = 1,size(matter%PCF,2)
      write(FN, '(f25.16,es25.16)') matter%PCF(1,i), matter%PCF(2,i)
   enddo
   write(FN, '(a)') ''
   write(FN, '(a)') ''
end subroutine write_PCF


subroutine save_distribution(FN, tim, wr, fe, fe_eq)
   integer, intent(in) :: FN
   real(8), intent(in) :: tim
   real(8), dimension(:), intent(in) :: wr
   real(8), dimension(:), intent(in) :: fe
   real(8), dimension(:), allocatable, intent(in) :: fe_eq
   integer i
   write(FN,'(a,f25.16)') '#', tim
   if (allocated(fe_eq)) then ! there is equivalent-temperature Fermi distribution
      do i = 1, size(fe)
         write(FN,'(f25.16,f25.16,f25.16)') wr(i), fe(i), fe_eq(i)
      enddo
   else  ! fe is Fermi, no equivalent distribution needed
      do i = 1, size(fe)
         write(FN,'(f25.16,f25.16)') wr(i), fe(i)
      enddo
   endif
   write(FN,*) ''
   write(FN,*) ''
end subroutine save_distribution



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


subroutine write_atomic_xyz(FN, atoms, matter, Supce)
   integer, intent(in) :: FN	! file number
   type(Atom), dimension(:), intent(in) :: atoms	! atomic parameters
   type(Solid), intent(in) :: matter	! Material parameters
   real(8), dimension(3,3), intent(in) :: Supce	! [A]  supercell vectors [a(x,y,z),b(x,y,z),c(x,y,z)]
   !-------------------------------------
   integer i
   character(10) :: Numb_out
   write(Numb_out, '(i10)') size(atoms)
   write(FN, '(a)') trim(adjustl(Numb_out))
   !write(FN, '(a)') trim(adjustl(matter%Name)) ! Material name, not needed
   write(FN, '(a,f,f,f,f,f,f,f,f,f,a)') 'Lattice="', Supce(:,:), '" Properties=species:S:1:pos:R:3'
   do i = 1, size(atoms)
      write(FN, '(a,es25.16,es25.16,es25.16)') trim(adjustl(matter%Atoms(atoms(i)%KOA)%Name)), atoms(i)%R(1), atoms(i)%R(2), atoms(i)%R(3)
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
      write(FN_out, '(a,es25.16,es25.16,es25.16)')  trim(adjustl(i_char))//'	'//trim(adjustl(matter%Atoms(atoms(i)%KOA)%Name)), atoms(i)%S(1), atoms(i)%S(2), atoms(i)%S(3)
   enddo
end subroutine write_atomic_cif



subroutine write_electron_properties(FN, time, Scell, NSC, Ei, matter, numpar, FN_Ce, FN_kappa, FN_Se)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: Ei	! energy levels
   type(Solid), intent(in) :: matter	! Material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   integer, intent(in) :: FN_Ce, FN_kappa, FN_Se  ! file number for band-resolved Ce and kappa, and electron entropy
   !------------------------
   integer i, Nat, n_at, Nsiz, norb, N_types, i_at, i_types, i_G1

   ! Write electron properties:
   write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)', advance='no') time, Scell(NSC)%Ne_low/real(Scell(NSC)%Ne)*100.0d0, Scell(NSC)%mu, Scell(NSC)%E_gap, Scell(NSC)%Ce, Scell(NSC)%G_ei, Scell(NSC)%E_VB_bottom, Scell(NSC)%E_VB_top, Scell(NSC)%E_bottom, Scell(NSC)%E_top
   Nat = size(matter%Atoms(:)) ! number of elements
   do i = 1, Nat    ! index starting from 11
      write(FN,'(es25.16)', advance='no') (matter%Atoms(i)%NVB - matter%Atoms(i)%mulliken_Ne)
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
      ! Total kappa:
      write(FN_kappa, '(es25.16,es25.16)', advance='no') time, Scell(NSC)%kappa_e
      ! All shells resolved:
      do i_at = 1, Nat
         do i_types = 1, N_types
            i_G1 = (i_at-1) * N_types + i_types
            write(FN_kappa,'(es25.16)',advance='no') Scell(NSC)%kappa_e_part(i_G1)
         enddo   ! i_types
      enddo ! i_at
      write(FN_kappa,'(a)') ''
   endif

   ! Write electron entropy:
   write(FN_Se, '(es25.16, es25.16, es25.16)') time, Scell(NSC)%Se, Scell(NSC)%Se_eq

end subroutine write_electron_properties




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
         !write(*,'(a,a,a)') trim(adjustl(matter%Atoms(i)%Name)), '-', trim(adjustl(matter%Atoms(j)%Name))//' '
      enddo
   enddo
   ! All shells resolved:
   do i_at = 1, N_at
      do i_types = 1, N_types
!          i_G1 = (i_at-1) * N_types + i_types
         chtemp1 = name_of_orbitals(norb, i_types) ! module "Little_subroutines"  
         do i_at2 = 1, N_at
            do i_types2 = 1, N_types
!                i_G2 = (i_at2-1) * N_types + i_types2
               chtemp2 = name_of_orbitals(norb, i_types2) ! module "Little_subroutines"
               write(FN,'(a,a,a)',advance='no') trim(adjustl(matter%Atoms(i_at)%Name))//'_'//trim(adjustl(chtemp1)), '--', &
                                        trim(adjustl(matter%Atoms(i_at2)%Name))//'_'//trim(adjustl(chtemp2))//'  '
                !write(*,'(a,a,a)') trim(adjustl(matter%Atoms(i_at)%Name))//'_'//trim(adjustl(chtemp1)), '--', &
                !                        trim(adjustl(matter%Atoms(i_at2)%Name))//'_'//trim(adjustl(chtemp2))//'  '
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
   ! All shells resolved:
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
      write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') atoms(i)%R(1), atoms(i)%R(2), atoms(i)%R(3), atoms(i)%V(1), atoms(i)%V(2), atoms(i)%V(3), atoms(i)%S(1), atoms(i)%S(2), atoms(i)%S(3), atoms(i)%SV(1), atoms(i)%SV(2), atoms(i)%SV(3)
   enddo
   write(FN, '(a)') ''
   write(FN, '(a)') ''
end subroutine write_atomic_relatives


subroutine write_numbers(FN, time, Scell)
   integer, intent(in) :: FN	! file number
   real(8), intent(in) :: time	! [fs]
   type(Super_cell), intent(in) :: Scell ! suoer-cell with all the atoms inside
   write(FN,'(f25.16,f25.16,es,es25.16,es25.16,es25.16,es25.16)') time, Scell%Ne_low/real(Scell%Na), Scell%Ne_CB/real(Scell%Na), Scell%Ne_high/real(Scell%Na), Scell%Nh/real(Scell%Na), (real(Scell%Ne)-(Scell%Ne_low+Scell%Ne_high-Scell%Nh))/real(Scell%Na), Scell%Nph/real(Scell%Na) 
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

   write(FN, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') time, nrg%E_tot, nrg%Eh_tot, nrg%At_pot+nrg%E_vdW+nrg%E_coul_scc, nrg%At_kin, nrg%Total, & !+nrg%E_supce,
   nrg%Total+nrg%E_supce+nrg%El_high, &
   nrg%Total+nrg%E_supce+nrg%El_high+nrg%Eh_tot, nrg%E_vdW
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
   write(FN,'(a)')  ! make new line
end subroutine write_temperatures_n_displacements


subroutine prepare_output_files(Scell,matter,laser,numpar,TB_Hamil,TB_Repuls,Err)
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
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
   logical :: file_opened, file_exist

   ! Create directory where the output files will be saved:
   call create_output_folder(Scell,matter,laser,numpar)	! module "Dealing_with_output_files"

   ! Save input files, so that repeating the same calculations would be easy:
   if (numpar%which_input >= 1) then
      
      ! Check if new format of input file exists:
      chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT'
      write(chtest2,'(i3)') numpar%which_input
      write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
      inquire(file=trim(adjustl(chtest)),exist=file_exist)
      
      if (.not.file_exist) then ! use old format of input files:
         chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT_MATERIAL'
         write(chtest2,'(i3)') numpar%which_input
         write(chtest,'(a,a,a,a)') trim(adjustl(chtest)), '_', trim(adjustl(chtest2)), '.txt'
         chtest1 = 'INPUT_DATA'//numpar%path_sep//'NUMERICAL_PARAMETERS'
         write(chtest1,'(a,a,a,a)') trim(adjustl(chtest1)), '_', trim(adjustl(chtest2)), '.txt'
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
      chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT.txt'
      inquire(file=trim(adjustl(chtest)),exist=file_exist)
      
      if (.not.file_exist) then ! use old format of input files:
         chtest = 'INPUT_DATA'//numpar%path_sep//'INPUT_MATERIAL.txt'
         chtest1 = 'INPUT_DATA'//numpar%path_sep//'NUMERICAL_PARAMETERS.txt'
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
      !call copy_file('INPUT_DATA'//numpar%path_sep//'INPUT_MATERIAL.txt',trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      !call copy_file('INPUT_DATA'//numpar%path_sep//'NUMERICAL_PARAMETERS.txt',trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
      if (.not.file_exist) call copy_file(trim(adjustl(chtest1)),trim(adjustl(numpar%output_path)),1) ! module "Dealing_with_output_files"
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
      !call copy_file('INPUT_DATA'//numpar%path_sep//'INPUT_MATERIAL.txt',trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      call copy_file(trim(adjustl(chtest)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      !call copy_file('INPUT_DATA'//numpar%path_sep//'NUMERICAL_PARAMETERS.txt',trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
      if (.not.file_exist) call copy_file(trim(adjustl(chtest1)),trim(adjustl(numpar%output_path))) ! module "Dealing_with_output_files"
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
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep//'Comunication.txt'
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
   call create_output_files(Scell,matter,laser,numpar)

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
   type(Super_cell), intent(in) :: Scell ! suoer-cell with all the atoms inside
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
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   !logical file_opened
   close(numpar%FN_temperatures)
   close(numpar%FN_pressure)
   close(numpar%FN_electron_properties)
   close(numpar%FN_Ce)
   close(numpar%FN_energies)
   close(numpar%FN_supercell)
   close(numpar%FN_numbers)
   close(numpar%FN_deep_holes)
   if (numpar%do_kappa) close(numpar%FN_kappa)
   if (numpar%save_raw) close(numpar%FN_atoms_S)
   if (numpar%do_drude) close(numpar%FN_optics)
   if (numpar%save_XYZ) close(numpar%FN_atoms_R)
   if (numpar%save_CIF) close(numpar%FN_cif)
   if (numpar%save_Ei)  close(numpar%FN_Ei)
   if (numpar%save_DOS)  close(numpar%FN_DOS)
   if (numpar%DOS_splitting == 1) close(numpar%FN_coupling)
   if (numpar%save_fe)  close(numpar%FN_fe)
   if (numpar%save_fe_grid)  close(numpar%FN_fe_on_grid)
   if (numpar%save_PCF) close(numpar%FN_PCF)
   if (Scell(1)%eps%all_w) close(numpar%FN_all_w)
   if (numpar%save_NN) close(numpar%FN_neighbors)
end subroutine close_output_files


subroutine create_output_files(Scell,matter,laser,numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   character(200) :: file_path, file_name
   integer :: FN, i, j, Nshl
   ! OUTPUT files, name and address:
   character(100) :: file_temperatures	! time [fs], Te [K], Ta [K]
   character(100) :: file_pressure	! time [fs], stress_tensore(3,3) [GPa], Pressure [GPa]
   character(100) :: file_energies	! energies [eV]
   character(100) :: file_atoms_R	! atomic coordinates and velocities
   character(100) :: file_atoms_S	! atomic coordinates and velocities
   character(100) :: file_atoms_cif	! atomic coordinates in cif-format (standard for constructing diffraction patterns)
   character(100) :: file_supercell	! supercell vectors
   character(100) :: file_electron_properties	! electron properties
   character(200) :: file_electron_heat_capacity	! band-resolved electron heat capacity
   character(200) :: file_electron_heat_conductivity  ! electron heat conductivity
   character(100) :: file_electron_entropy	! electron entropy
   character(100) :: file_numbers	! total numbers of electrons and holes
   character(100) :: file_deep_holes	! number of deep-shell holes in each shell
   character(100) :: file_Ei		! energy levels
   character(100) :: file_DOS	! DOS
   character(100) :: file_coupling  ! partial coupling parameter
   character(100) :: file_fe		! electron distribution (low-energy part)
   character(100) :: file_fe_on_grid   ! electron distribution (full: low- + high-energy)
   character(100) :: file_PCF		! pair correlation function
   character(100) :: file_optics	! optical coefficients
   character(100) :: file_all_w		! optical coeffs for all hw
   character(100) :: file_NN		! nearest neighbors
   character(100) :: chtemp
   character(11) :: chtemp11

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
   
   file_pressure = trim(adjustl(file_path))//'OUTPUT_pressure_and_stress.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_pressure)))
   numpar%FN_pressure = FN
   call create_file_header(numpar%FN_pressure, '#Time	Pressure	Pressure(x,x)	Pressure(x,y)	Pressure(x,z)	Pressure(y,x)	Pressure(y,y)	Pressure(y,z)	Pressure(z,x)	Pressure(z,y)	Pressure(z,z)')
   call create_file_header(numpar%FN_pressure, '#[fs]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]	[GPa]')

   file_electron_properties = trim(adjustl(file_path))//'OUTPUT_electron_properties.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_properties)))
   numpar%FN_electron_properties = FN
   call create_file_header(numpar%FN_electron_properties, '#Time	Ne	mu	band_gap	Ce	Coupling_parameter	VB_bottom	VB_top	CB_bottom	CB_top Mullikens(:)')
   call create_file_header(numpar%FN_electron_properties, '#[fs]	[%]	[eV]	[eV]	[J/(m^3K)]	[W/(m^3K)]	[eV]	[eV]	[eV]	[eV]  [e](:)')

   file_electron_heat_capacity = trim(adjustl(file_path))//'OUTPUT_electron_Ce.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_heat_capacity)))
   numpar%FN_Ce = FN
   call write_Ce_header(numpar%FN_Ce, Scell, 1, matter) ! below

   file_electron_entropy = trim(adjustl(file_path))//'OUTPUT_electron_entropy.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_entropy)))
   numpar%FN_Se = FN
   call create_file_header(numpar%FN_Se, '#Time Se  Se_eq')
   call create_file_header(numpar%FN_electron_properties, '#[fs]  [eV/K]   [eV/K]')

   if (numpar%do_kappa) then
      file_electron_heat_conductivity = trim(adjustl(file_path))//'OUTPUT_electron_heat_conductivity.dat'
      open(NEWUNIT=FN, FILE = trim(adjustl(file_electron_heat_conductivity)))
      numpar%FN_kappa = FN
      ! We can use the same header here as for Ce:
      call write_Ce_header(numpar%FN_kappa, Scell, 1, matter) ! below
   endif

   file_energies = trim(adjustl(file_path))//'OUTPUT_energies.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_energies)))
   numpar%FN_energies = FN
   call create_file_header(numpar%FN_energies, '#Time	Electrons	Holes	Potential	Kinetic	Atoms	Atoms_n_electrons	Atom_all_electrons	Total	van_der_Waals')
   call create_file_header(numpar%FN_energies, '#[fs]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]	[eV/atom]')

   file_numbers = trim(adjustl(file_path))//'OUTPUT_electron_hole_numbers.dat'
   open(NEWUNIT=FN, FILE = trim(adjustl(file_numbers)))
   numpar%FN_numbers = FN
   write(numpar%FN_numbers, '(a)') '#Time	VB_electrons	CB_electrons	High_energy_electrons	Deep_holes	Error	Photons'
   write(numpar%FN_numbers, '(a)') '#[fs]	[1/atom]	[1/atom]	[1/atom]	[1/atom]	[1/atom]	[1/atom]'

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

   do i = 1, size(Scell)
      if (Scell(i)%eps%all_w) then
         file_all_w = trim(adjustl(file_path))//'OUTPUT_dielectric_function.dat'
         open(NEWUNIT=FN, FILE = trim(adjustl(file_all_w)))
         numpar%FN_all_w = FN
      endif
   enddo

   call create_gnuplot_scripts(Scell,matter,numpar,laser, file_path, 'OUTPUT_temperatures.dat',  'OUTPUT_pressure_and_stress.dat', 'OUTPUT_energies.dat', file_atoms_R, file_atoms_S, 'OUTPUT_supercell.dat', 'OUTPUT_electron_properties.dat', 'OUTPUT_electron_hole_numbers.dat', 'OUTPUT_deep_shell_holes.dat', 'OUTPUT_optical_coefficients.dat', file_Ei, file_PCF, 'OUTPUT_nearest_neighbors.dat', 'OUTPUT_electron_entropy.dat')  ! below
   !call create_gnuplot_scripts(matter,numpar,laser, file_path, file_temperatures, file_energies, file_atoms_R, file_atoms_S, file_supercell, file_electron_properties, file_numbers, file_deep_holes, file_Ei, file_PCF)
end subroutine create_output_files


subroutine create_file_header(FN, text)
   integer, intent(in) :: FN		! file number
   character(*), intent(in) :: text	! what to write in this file
   write(FN,'(a)') trim(adjustl(text))
end subroutine create_file_header


subroutine create_gnuplot_scripts(Scell,matter,numpar,laser, file_path, file_temperatures, file_pressure, file_energies, file_atoms_R, file_atoms_S, file_supercell, file_electron_properties, file_numbers, file_deep_holes, file_optics, file_Ei, file_PCF, file_NN, file_electron_entropy)
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(in) :: matter
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), intent(in) :: laser		! Laser pulse parameters
   character(*) :: file_path
   character(*) :: file_temperatures	! time [fs], Te [K], Ta [K]
   character(*) :: file_pressure	! pressure and stress tensore
   character(*) :: file_energies	! energies [eV]
   character(*) :: file_atoms_R	! atomic coordinates and velocities
   character(*) :: file_atoms_S	! atomic coordinates and velocities
   character(*) :: file_supercell	! supercell vectors
   character(*) :: file_electron_properties	! electron properties
   character(*) :: file_numbers	! total numbers of electrons and holes
   character(*) :: file_deep_holes	! deep-shell holes
   character(*) :: file_optics		! optical coefficients
   character(*) :: file_Ei		! energy levels
   character(*) :: file_PCF		! pair correlation function
   character(*) :: file_NN      ! nearest neighbors
   character(*) :: file_electron_entropy  ! electron netropy
   !----------------
   character(200) :: File_name, File_name2
   real(8) :: t0, t_last, x_tics
   integer FN, i, j, Nshl, counter, iret
   character(200) :: chtemp, command
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
   
   ! Shell-script for execution of all gnuplot files:
   File_name  = trim(adjustl(file_path))//'OUTPUT_Gnuplot_all'//trim(adjustl(sh_cmd))
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
   
   if (numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a)') '@echo off'

      write(FN, '(a)') 'echo Executing OUTPUT_energies_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_energies_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_temperatures_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_temperatures_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_mean_displacement_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_mean_displacement_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_pressure_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_pressure_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_stress_tensor_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_stress_tensor_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_electrons_and_holes_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_electrons_and_holes_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_CB_electrons_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_CB_electrons_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_deep_shell_holes_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_deep_shell_holes_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_volume_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_volume_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_mu_and_Egap_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_mu_and_Egap_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_bands_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_bands_Gnuplot'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_electron_Ce'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_electron_Ce'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_coupling_parameter'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_coupling_parameter'//trim(adjustl(sh_cmd))

      write(FN, '(a)') 'echo Executing OUTPUT_electron_entropy'//trim(adjustl(sh_cmd))
      write(FN, '(a)') 'call OUTPUT_electron_entropy'//trim(adjustl(sh_cmd))

      if (numpar%do_drude) then 
         write(FN, '(a)') 'echo Executing OUTPUT_optical_coefficients'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_optical_coefficients'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'echo Executing OUTPUT_optical_n_and_k'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_optical_n_and_k'//trim(adjustl(sh_cmd))
      endif
      if (numpar%save_Ei) then
         write(FN, '(a)') 'echo Executing OUTPUT_energy_levels_Gnuplot'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_energy_levels_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (numpar%save_fe) then
         write(FN, '(a)') 'echo Executing OUTPUT_electron_distribution_Gnuplot'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_electron_distribution_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (numpar%save_fe_grid) then
         write(FN, '(a)') 'echo Executing OUTPUT_electron_distribution_on_grid_Gnuplot'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_electron_distribution_on_grid_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (numpar%DOS_splitting >= 1) then   ! Mulliken charges
         if (numpar%Mulliken_model >= 1) then
            write(FN, '(a)') 'echo Executing OUTPUT_Mulliken_charges_Gnuplot'//trim(adjustl(sh_cmd))
            write(FN, '(a)') 'call OUTPUT_Mulliken_charges_Gnuplot'//trim(adjustl(sh_cmd))
         endif
      endif
      if (numpar%save_NN) then
         write(FN, '(a)') 'echo Executing OUTPUT_neighbors_Gnuplot'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_neighbors_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (Scell(1)%eps%tau > 0.0d0) then ! convolved files too:
         write(FN, '(a)') 'call Executing convolved files...'
         write(FN, '(a)') 'call OUTPUT_optical_coefficients_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_optical_n_and_k_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_energies_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_temperatures_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_mean_displacement_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))         
         write(FN, '(a)') 'call OUTPUT_pressure_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_stress_tensor_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_electrons_and_holes_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_CB_electrons_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_deep_shell_holes_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_volume_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_mu_and_Egap_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_bands_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_electron_Ce_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_coupling_parameter_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') 'call OUTPUT_electron_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
         if (numpar%Mulliken_model >= 1) then
            write(FN, '(a)') 'call OUTPUT_Mulliken_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         endif
      endif
   else ! It is linux
      write(FN, '(a)') '#!/bin/bash'
      !write(FN, '(a)') 'cd ../'
      write(FN, '(a)') './OUTPUT_energies_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_temperatures_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_mean_displacement_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_pressure_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_stress_tensor_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_electrons_and_holes_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_CB_electrons_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_deep_shell_holes_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_volume_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_mu_and_Egap_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_bands_Gnuplot'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_electron_Ce'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_coupling_parameter'//trim(adjustl(sh_cmd))
      write(FN, '(a)') './OUTPUT_electron_entropy'//trim(adjustl(sh_cmd))
      if (numpar%do_drude) then 
         write(FN, '(a)') './OUTPUT_optical_coefficients'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_optical_n_and_k'//trim(adjustl(sh_cmd))
      endif
      if (numpar%save_Ei) then
         write(FN, '(a)') './OUTPUT_energy_levels_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (numpar%save_fe) then
         write(FN, '(a)') './OUTPUT_electron_distribution_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (numpar%save_fe_grid) then
         write(FN, '(a)') './OUTPUT_electron_distribution_on_grid_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (numpar%DOS_splitting >= 1) then   ! Mulliken charges
         if (numpar%Mulliken_model >= 1) then
            write(FN, '(a)') './OUTPUT_Mulliken_charges_Gnuplot'//trim(adjustl(sh_cmd))
         endif
      endif
      if (numpar%save_NN) then
         write(FN, '(a)') './OUTPUT_neighbors_Gnuplot'//trim(adjustl(sh_cmd))
      endif
      if (Scell(1)%eps%tau > 0.0d0) then ! convolved files too:
         write(FN, '(a)') './OUTPUT_optical_coefficients_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_optical_n_and_k_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_energies_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_temperatures_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_mean_displacement_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_pressure_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_stress_tensor_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_electrons_and_holes_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_CB_electrons_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_deep_shell_holes_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_volume_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_mu_and_Egap_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_bands_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_electron_Ce_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_coupling_parameter_CONVOLVED'//trim(adjustl(sh_cmd))
         write(FN, '(a)') './OUTPUT_electron_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
         if (numpar%Mulliken_model >= 1) then
            write(FN, '(a)') './OUTPUT_Mulliken_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
         endif
      endif
   endif
   write(FN, '(a)') 'echo Gnuplot scripts are executed, see created plots.'
   close(FN)
   if (numpar%path_sep .EQ. '\') then	! if it is Windows
   else ! It is linux
      !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
      command = 'chmod +x '//trim(adjustl(File_name))
      iret = system(command)
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
   
   ! Pressure:
   File_name  = trim(adjustl(file_path))//'OUTPUT_pressure_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_pressure(File_name, file_pressure, t0, t_last, 'OUTPUT_pressure.'//trim(adjustl(numpar%fig_extention))) ! below
   
   ! Stress tensor:
   File_name  = trim(adjustl(file_path))//'OUTPUT_stress_tensor_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_stress(File_name, file_pressure, t0, t_last, 'OUTPUT_pressure_tensor.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Numbers of particles:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_and_holes_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_numbers(File_name, file_numbers, t0, t_last, 'OUTPUT_electrons_holes.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Numbers of CB electrons:
   File_name  = trim(adjustl(file_path))//'OUTPUT_CB_electrons_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_CB_electrons(File_name, file_numbers, t0, t_last, 'OUTPUT_CB_electrons.'//trim(adjustl(numpar%fig_extention)))

   ! Numbers of deep-shell holes:
   File_name  = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_holes(File_name, file_deep_holes, t0, t_last, matter, 'OUTPUT_deep_shell_holes.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Chemical potential and band gap:
   File_name  = trim(adjustl(file_path))//'OUTPUT_mu_and_Egap_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_Egap(File_name, file_electron_properties, t0, t_last, 'OUTPUT_mu_and_Egap.'//trim(adjustl(numpar%fig_extention))) ! below
   
   ! Boundaries of the bands:
   File_name  = trim(adjustl(file_path))//'OUTPUT_bands_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_Ebands(File_name, file_electron_properties, t0, t_last, 'OUTPUT_bands.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron heat capacity:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electron_Ce'//trim(adjustl(sh_cmd))
   call gnu_capacity(File_name, file_electron_properties, t0, t_last, 'OUTPUT_electron_Ce.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron entropy:
   File_name  = trim(adjustl(file_path))//'OUTPUT_electron_entropy'//trim(adjustl(sh_cmd))
   call gnu_entropy(File_name, file_electron_entropy, t0, t_last, 'OUTPUT_electron_entropy.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Electron-ion coupling parameter:
   File_name  = trim(adjustl(file_path))//'OUTPUT_coupling_parameter'//trim(adjustl(sh_cmd))
   call gnu_coupling(File_name, file_electron_properties, t0, t_last, 'OUTPUT_coupling.'//trim(adjustl(numpar%fig_extention))) ! below

   ! Volume:
   File_name  = trim(adjustl(file_path))//'OUTPUT_volume_Gnuplot'//trim(adjustl(sh_cmd))
   call gnu_volume(File_name, file_supercell, t0, t_last, 'OUTPUT_volume.'//trim(adjustl(numpar%fig_extention))) ! below

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
   
   ! Distribution function of electrons:
   if (numpar%save_fe) then
      !! Find order of the number, and set number of tics as tenth of it:
      !call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

      ! Distribution function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_distribution_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 5.0d0, 'Distribution', 'Energy (eV)', 'Electron distribution (a.u.)', 'OUTPUT_electron_distribution.gif', numpar%path_sep, setkey=0)
      !call write_energy_levels_gnuplot(FN, Scell, 'OUTPUT_electron_distribution.dat')
      call write_distribution_gnuplot(FN, Scell, numpar, 'OUTPUT_electron_distribution.dat')   ! below
      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)
   endif

   ! Distribution function of all electrons on the grid:
   if (numpar%save_fe_grid) then
      ! Distribution function can only be plotted as animated gif:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_distribution_on_grid_Gnuplot'//trim(adjustl(sh_cmd))
      open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
      call write_gnuplot_script_header_new(FN, 6, 1.0d0, 10.0d0, 'Distribution', 'Energy (eV)', 'Electron density (1/(V*E))', 'OUTPUT_electron_distribution_on_grid.gif', numpar%path_sep, setkey=0)
      !call write_energy_levels_gnuplot(FN, Scell, 'OUTPUT_electron_distribution.dat')
      call write_distribution_on_grid_gnuplot(FN, Scell, numpar, 'OUTPUT_electron_distribution_on_grid.dat')   ! below
      call write_gnuplot_script_ending(FN, File_name, 1)
      close(FN)
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
      
      ! Pressure:
      File_name  = trim(adjustl(file_path))//'OUTPUT_pressure_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_pressure(File_name, trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_pressure_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Stress tensor:
      File_name  = trim(adjustl(file_path))//'OUTPUT_stress_tensor_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_stress(File_name, trim(adjustl(file_pressure(1:len(trim(adjustl(file_pressure)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_pressure_tensor_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      
      ! Numbers of particles:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electrons_and_holes_Gnu_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_numbers(File_name, trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electrons_holes_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Numbers of CB electrons:
      File_name  = trim(adjustl(file_path))//'OUTPUT_CB_electrons_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_CB_electrons(File_name, trim(adjustl(file_numbers(1:len(trim(adjustl(file_numbers)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_CB_electrons_CONVOLVED.'//trim(adjustl(numpar%fig_extention)))

      ! Numbers of deep-shell holes:
      File_name  = trim(adjustl(file_path))//'OUTPUT_deep_shell_holes_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_holes(File_name, trim(adjustl(file_deep_holes(1:len(trim(adjustl(file_deep_holes)))-4)))//'_CONVOLVED.dat', t0, t_last, matter, 'OUTPUT_deep_shell_holes_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Chemical potential and band gap:
      File_name  = trim(adjustl(file_path))//'OUTPUT_mu_and_Egap_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_Egap(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_mu_and_Egap_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below
      
      ! Boundaries of the bands:
      File_name  = trim(adjustl(file_path))//'OUTPUT_bands_Gnuplot_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_Ebands(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_bands_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron heat capacity:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_Ce_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_capacity(File_name, trim(adjustl(file_electron_properties(1:len(trim(adjustl(file_electron_properties)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electron_Ce_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

      ! Electron entropy:
      File_name  = trim(adjustl(file_path))//'OUTPUT_electron_entropy_CONVOLVED'//trim(adjustl(sh_cmd))
      call gnu_entropy(File_name,      trim(adjustl(file_electron_entropy(1:len(trim(adjustl(file_electron_entropy)))-4)))//'_CONVOLVED.dat', t0, t_last, 'OUTPUT_electron_entropy_CONVOLVED.'//trim(adjustl(numpar%fig_extention))) ! below

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

   !call write_gnuplot_script_header(FN, 1, 3, 'Numbers','Time (fs)', 'Particles per atoms (arb.units)', trim(adjustl(file_path))//'OUTPUT_electrons_holes.'//trim(adjustl(g_numpar%fig_extention)))
   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Numbers','Time (fs)', 'Particles per atoms (arb.units)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Numbers','Time (fs)', 'Particles (per atoms)', trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"

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
   integer :: FN, counter, Nshl, i, j, Na
   character(100) :: chtemp
   character(11) :: chtemp11
   real(8) :: x_tics
   character(8) :: temp, time_order
   logical :: first_line
   
   open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), action="write", status="replace")
   
    ! Find order of the number, and set number of tics as tenth of it:
   call order_of_time((t_last - t0), time_order, temp, x_tics)	! module "Little_subroutines"

   !call write_gnuplot_script_header(FN, 1, 3, 'Holes','Time (fs)', 'Particles per atoms (arb.units)', trim(adjustl(file_path))//'OUTPUT_deep_shell_holes.'//trim(adjustl(g_numpar%fig_extention)))
   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Holes','Time (fs)', 'Particles per atoms (arb.units)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Holes','Time (fs)', 'Particles (total)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   counter = 0 ! to start with
   first_line = .true.  ! to start from the first line
   W_vs_L:if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
     ATOMS0:do i = 1, size(matter%Atoms) ! for all atoms
         Na = size(matter%Atoms)
         Nshl = size(matter%Atoms(i)%Ip)
         SHELLS0:do j = 1, Nshl ! for all shells of this atom
            if ((i .NE. 1) .or. (j .NE. Nshl)) then ! atomic shell:
               counter = counter + 1
               call define_PQN(matter%Atoms(i)%Shl_dsgnr(j), chtemp11) ! module "Dealing_with_EADL"
               write(chtemp,'(a,a,a)') trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(chtemp11))
               select case(size(matter%Atoms))
               case (1)
                  !if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] "' , trim(adjustl(file_deep_holes)), '" u 1:', 1+j ,' w l lw LW title "', trim(adjustl(chtemp))  ,'"'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") ' "', trim(adjustl(file_deep_holes)), '" u 1:', 1+j ,' w l lw LW title "', trim(adjustl(chtemp))  ,'"'
                  endif
                  if ((i .NE. size(matter%Atoms)) .OR. (j .LT. Nshl-1)) then
                     write(FN, '(a)') ',\'
                  else
                     write(FN, '(a)') ''
                  endif
               case default
                  !if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] "' , trim(adjustl(file_deep_holes)), ' "u 1:', 1+counter,' w l lw LW title " ', trim(adjustl(chtemp))  ,' "'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") ' "', trim(adjustl(file_deep_holes)), '" u 1:', 1+counter,' w l lw LW title "', trim(adjustl(chtemp))  ,'"'
                  endif
                  if ((i .NE. size(matter%Atoms)) .OR. (j .LT. Nshl)) then
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
      ATOMS:do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         SHELLS:do j = 1, Nshl ! for all shells of this atom
            if ((i .NE. 1) .or. (j .NE. Nshl)) then ! atomic shell:
               counter = counter + 1
               call define_PQN(matter%Atoms(i)%Shl_dsgnr(j), chtemp11) ! module "Dealing_with_EADL"
               write(chtemp,'(a,a,a)') trim(adjustl(matter%Atoms(i)%Name))//' '//trim(adjustl(chtemp11))
               select case(size(matter%Atoms))
               case (1)
!                   if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] \"' , trim(adjustl(file_deep_holes)), '\"u 1:', 1+j ,' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))  ,'\"'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") '\"', trim(adjustl(file_deep_holes)), '\"u 1:', 1+j ,' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))  ,'\"'
                  endif
                  if ((i .NE. size(matter%Atoms)) .OR. (j .LT. Nshl-1)) then
                     write(FN, '(a)') ',\'
                  else
                     write(FN, '(a)') ''
                  endif
               case default
!                   if ((i == 1) .and. (j == 1)) then
                  if (first_line) then
                     first_line = .false. ! first line is done, don't repeat it
                     write(FN, '(a,es25.16,a,a,a,i3,a,a,a)', ADVANCE = "NO") 'p [', t0, ':][] \"' , trim(adjustl(file_deep_holes)), '\"u 1:', 1+counter,' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))  ,'\"'
                  else
                     write(FN, '(a,a,a,i3,a,a,a)', ADVANCE = "NO") '\"', trim(adjustl(file_deep_holes)), '\"u 1:', 1+counter,' w l lw \"$LW\" title \" ', trim(adjustl(chtemp))  ,'\"'
                  endif
                  if ((i .NE. size(matter%Atoms)) .OR. (j .LT. Nshl)) then
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

   !call write_gnuplot_script_header(FN, 1, 3, 'mu and Egap','Time (fs)', 'Energy (eV)', trim(adjustl(file_path))//'OUTPUT_mu_and_Egap.'//trim(adjustl(g_numpar%fig_extention)))
   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'mu and Egap','Time (fs)', 'Energy (eV)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'mu and Egap','Time (fs)', 'Energy (eV)', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:3 w l lw LW title "Chemical potnential" ,\'
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), ' "u 1:4 w l lw LW title "Band gap" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:3 w l lw \"$LW\" title \"Chemical potnential\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:4 w l lw \"$LW\" title \"Band gap \" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_Egap


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
      write(FN, '(a,a,a,i12,a)') ' "', trim(adjustl(file_electron_properties)), ' "u 1:3 w l lt rgb "#000000" lw LW title "Chemical potnential" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:10 w l lw \"$LW\" title \"Top of CB \" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:9 w l lw \"$LW\" title \"Bottom of CB\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:8 w l lw \"$LW\" title \"Top of VB\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:7 w l lw \"$LW\" title \"Bottom of VB\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_properties)), '\"u 1:3 w l lt rgb \"#000000\"  lw \"$LW\" title \"Chemical potnential\" '
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
      write(FN, '(a,a,a,i12,a)') '"', trim(adjustl(file_electron_entropy)), '" u 1:2 w l lw LW title "Nonequilibrium" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_entropy)), '\" u 1:3 w l lw \"$LW\" title \"Equilibrium\" ,\'
      write(FN, '(a,a,a,i12,a)') '\"', trim(adjustl(file_electron_entropy)), '\" u 1:2 w l lw \"$LW\" title \"Nonequilibrium\" '
   endif
   call write_gnuplot_script_ending(FN, File_name, 1)
   close(FN)
end subroutine gnu_entropy


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

   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Coupling parameter','Time (fs)', 'Coupling parameter (W/(m^3 K))', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Coupling parameter','Time (fs)', 'Coupling parameter (W/(m^3 K))', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_electron_properties)), ' "u 1:6 w l lw LW title "Electron-ion coupling" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_electron_properties)), '\"u 1:6 w l lw \"$LW\" title \"Electron-ion coupling\" '
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

   !call write_gnuplot_script_header(FN, 1, 3, 'Volume','Time (fs)', 'Volume (A^3)', trim(adjustl(file_path))//'OUTPUT_volume.'//trim(adjustl(g_numpar%fig_extention)))
   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Volume','Time (fs)', 'Volume (A^3)', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Volume','Time (fs)', 'Volume (A^3)', trim(adjustl(eps_name)), g_numpar%path_sep, 1)	! module "Gnuplotting"
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] "' , trim(adjustl(file_supercell)), ' "u 1:2 w l lw LW title "Supercell volume" '
   else ! It is linux
      write(FN, '(a,es25.16,a,a,a)') 'p [', t0, ':][] \"' , trim(adjustl(file_supercell)), '\"u 1:2 w l lw \"$LW\" title \"Supercell volume\" '
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

   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Optical coefficients','Time (fs)', 'Optical coefficients', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Optical coefficients','Time (fs)', 'Optical coefficients', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
      
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

   !call write_gnuplot_script_header(FN, 1, 3.0d0, 'Optical n and k','Time (fs)', 'Optical parameters', trim(adjustl(eps_name)))
   call write_gnuplot_script_header_new(FN, g_numpar%ind_fig_extention, 3.0d0, x_tics, 'Optical n and k','Time (fs)', 'Optical parameters', trim(adjustl(eps_name)), g_numpar%path_sep, 0)	! module "Gnuplotting"
   
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
   character(100) :: command
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
   character(100) :: command
   integer :: iret
   
   !call chdir(trim(adjustl(file_path)))
   command = trim(adjustl(file_path))
   iret = chdir(command)
   
   if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
      !call system("OUTPUT_Gnuplot_all.cmd")
      command = "OUTPUT_Gnuplot_all.cmd"
      iret = system(command)
   else ! linux:
      !call system("./OUTPUT_Gnuplot_all.sh")
      command = "./OUTPUT_Gnuplot_all.sh"
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
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
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


subroutine write_distribution_gnuplot(FN, Scell, numpar, file_fe)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_fe  ! file with electronic distribution function
   !-----------------------
   integer :: i, M, NSC
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4
   logical :: do_fe_eq

   do NSC = 1, size(Scell)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') 25.0d0      ! Scell(NSC)%E_top
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f)') numpar%dt_save

      select case (numpar%el_ion_scheme)
         case (3:4)
            do_fe_eq = .true.
         case default
            do_fe_eq = .false.
      endselect
      ! minimal energy grid:
      write(ch_temp4,'(f)') -25.0d0  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)

      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         !if (do_fe_eq) then  ! plot also equivalent Fermi distribution

            write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][0:2] "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:3 w l lw 2 lt rgb "grey" title "Equivalent Fermi" ,\'
            write(FN, '(a)') ' "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i-1+'// &
                  trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
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
            write(FN, '(a)') ' \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i-1+'// &
                  trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
         !else
         !   write(FN, '(a)') 'p [:'//trim(adjustl(ch_temp))//'][0:2] \"'//trim(adjustl(file_fe))// &
         !         '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i-1+'// &
         !         trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
         !endif
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_distribution_gnuplot


subroutine write_distribution_on_grid_gnuplot(FN, Scell, numpar, file_fe)
   integer, intent(in) :: FN            ! file to write into
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
   type(Numerics_param), intent(in) :: numpar   ! all numerical parameters
   character(*), intent(in) :: file_fe  ! file with electronic distribution function
   !-----------------------
   integer :: i, M, NSC
   character(30) :: ch_temp, ch_temp2, ch_temp3, ch_temp4
   logical :: do_fe_eq

   do NSC = 1, size(Scell)
      ! Choose the maximal energy, up to what energy levels should be plotted [eV]:
      write(ch_temp,'(f)') 100.0d0      ! Scell(NSC)%E_top
      write(ch_temp2,'(f)') numpar%t_start
      write(ch_temp3,'(f)') numpar%dt_save

      select case (numpar%el_ion_scheme)
         case (3:4)
            do_fe_eq = .true.
         case default
            do_fe_eq = .false.
      endselect
      ! minimal energy grid:
      write(ch_temp4,'(f)') -25.0d0  ! (FLOOR(Scell(NSC)%E_bottom/10.0d0)*10.0)
      if (g_numpar%path_sep .EQ. '\') then	! if it is Windows
         write(FN, '(a)') 'stats "'//trim(adjustl(file_fe))//'" nooutput'
         write(FN, '(a)') 'set logscale y'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] "'//trim(adjustl(file_fe))// &
                  '" index (i-1) u 1:2 pt 7 ps 1 title sprintf("%i fs",(i-1+'// &
                  trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
      else  ! Linux
         write(FN, '(a)') 'stats \"'//trim(adjustl(file_fe))//'\" nooutput'
         write(FN, '(a)') 'set logscale y'
         write(FN, '(a)') 'do for [i=1:int(STATS_blocks)] {'
         write(FN, '(a)') 'p ['//trim(adjustl(ch_temp4))//':'//trim(adjustl(ch_temp))//'][1e-6:] \"'//trim(adjustl(file_fe))// &
                  '\" index (i-1) u 1:2 pt 7 ps 1 title sprintf(\"%i fs\",(i-1+'// &
                  trim(adjustl(ch_temp2))// ')/' // trim(adjustl(ch_temp3)) //') '
      endif
      write(FN, '(a)') '}'
   enddo
end subroutine write_distribution_on_grid_gnuplot


subroutine output_parameters_file(Scell,matter,laser,numpar,TB_Hamil,TB_Repuls,Err)
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
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
   path_sep = numpar%path_sep
   File_name = trim(adjustl(numpar%output_path))//path_sep
   File_name = trim(adjustl(File_name))//'!OUTPUT_'//trim(adjustl(matter%Name))//'_Parameters.txt'

   !open(NEWUNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
   FN = 111
   open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      INFO = 2
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      call Save_error_details(Err, INFO, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 9999
   endif
   numpar%FN_parameters = FN ! save this file number with parameters
   call Print_title(FN,Scell, matter,laser,numpar)
   !close(FN)
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) then
      write(FN, '(a)') 'Atomic data used for '//trim(adjustl(matter%Name))//' are:'
      do i = 1, matter%N_KAO
         write(chtemp(1), '(i12)') i
         write(chtemp(2), '(f6.2)') matter%Atoms(i)%percentage
         write(FN, '(a,$)') 'Element #'//trim(adjustl(chtemp(1)))//' is '//trim(adjustl(matter%Atoms(i)%Name))//' contributing to the compound with '//trim(adjustl(chtemp(2)))
         write(FN, '(a)') ''
         write(chtemp(1), '(i4)') INT(matter%Atoms(i)%Z)
         write(chtemp(2), '(es25.5)') matter%Atoms(i)%Ma
         write(FN, '(a,a,a,a,a)') 'Atomic number: ', trim(adjustl(chtemp(1))), ', mass: ', trim(adjustl(chtemp(2))), ' [kg]'

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
      write(FN,'(a)') '*************************************************************'
   endif
9999 continue
end subroutine output_parameters_file


! Create the folder where the results will be storred:
subroutine create_output_folder(Scell,matter,laser,numpar)
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
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

   LAS:if (maxval(laser(:)%F) .GT. 0.0d0) then
      write(ch1,'(f7.1)') (laser(1)%hw)	! photon energy
      if (laser(1)%KOP .EQ. 1) then
         write(ch2,'(f6.1)') (laser(1)%t*2.35482d0)	! pulse duration
      else
         write(ch2,'(f6.1)') laser(1)%t		! pulse duration
      endif
      write(ch3,'(f6.2)') laser(1)%F	! dose [eV/atom]

      if (numpar%path_sep .EQ. '\') then	! if it is Windows
         if (size(laser) .GT. 1) then
            write(ch4,'(i2)') size(laser)
            write(File_name,'(a,a,a,a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw_', trim(adjustl(ch1)), '_t_', trim(adjustl(ch2)), '_F_', trim(adjustl(ch3)), '_', trim(adjustl(ch4)), '_pulses'
         else ! singe pulse
            write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw_', trim(adjustl(ch1)), '_t_', trim(adjustl(ch2)), '_F_', trim(adjustl(ch3))
         endif
      else ! it is linux
         if (size(laser) .GT. 1) then
            write(ch4,'(i2)') size(laser)
            write(File_name,'(a,a,a,a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw=', trim(adjustl(ch1)), '_t=', trim(adjustl(ch2)), '_F=', trim(adjustl(ch3)), '_', trim(adjustl(ch4)), '_pulses'
         else ! singe pulse
            write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_hw=', trim(adjustl(ch1)), '_t=', trim(adjustl(ch2)), '_F=', trim(adjustl(ch3))
         endif
      endif 
   else LAS
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
         write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_Te_', trim(adjustl(ch1)), '_Ta_', trim(adjustl(ch2)), '_', trim(adjustl(ch3))
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
         write(File_name,'(a,a,a,a,a,a,a,a)') 'OUTPUT_', trim(adjustl(matter_name)), '_Te=', trim(adjustl(ch1)), '_Ta=', trim(adjustl(ch2)), '_', trim(adjustl(ch3))
      endif
   endif LAS

   ! Do TB and MD part only if we want (supercell is larger than 0):
   if (matter%cell_x*matter%cell_y*matter%cell_z .LE. 0) then
      write(File_name,'(a,a)') trim(adjustl(File_name)), '_MC_only'
   endif 
   File_name2 = File_name
   i = 0
   inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
   do while (file_exist)
      i = i + 1
      write(ch1,'(i6)') i
      write(File_name2,'(a,a,a)') trim(adjustl(File_name)), '_v', trim(adjustl(ch1))
      inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
   enddo
   command='mkdir '//trim(adjustl(File_name2)) ! to create a folder use this command
   !CALL system(command)  ! create the folder
   iret = system(command)
   numpar%output_path = File_name2
end subroutine create_output_folder


! Coomunication with the used via file:
subroutine communicate(FN, time, numpar, matter)
   integer, intent(in) :: FN ! file number to read from
   real(8), intent(in) :: time ! current time [fs]
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   type(Solid), intent(in) :: matter ! parameters of the material
   integer :: Reason, i, MOD_TIM, sz
   character(200) :: readline, given_line, File_name
   real(8) given_num
   logical :: read_well, read_well_2, file_opened
   
   File_name = trim(adjustl(numpar%output_path))//numpar%path_sep//'Comunication.txt'
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

!    prinT*, 'communicate', MOD_TIM, numpar%MOD_TIME, file_opened

   COM_OPEN:if (file_opened) then ! read it
      rewind(FN)  ! to start reading from the start
      i = 1 ! to start with
      read_well = .true.   ! to start with
      Reason = 1  ! to start with
      do while (Reason >= 0) ! read all lines if there is more than one
         call pars_comunications_file(FN, i, given_line, given_num, Reason) ! below
         if (Reason == 0) call act_on_comunication(read_well_2, given_line, given_num, numpar, matter, time)   ! below
      enddo
      rewind(FN)
      write(FN,'(a)') ''
      rewind(FN)

      call get_file_stat(trim(adjustl(File_name)), Last_modification_time=MOD_TIM) ! module 'Dealing_with_files'
      if (MOD_TIM /= numpar%MOD_TIME) then ! if it was modified by the user, then
         numpar%MOD_TIME = MOD_TIM         ! save new time of the last modification
      endif

      close(FN) ! we have to close the file to let the user write into it
   endif COM_OPEN
end subroutine communicate



subroutine save_duration(matter, numpar, chtext)
   type(Solid), intent(in) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   character(*), intent(in) :: chtext ! time duration to print out
   integer FN
   logical file_opened, file_exists
   character(200) :: File_name
   character(1) path_sep
   path_sep = trim(adjustl(numpar%path_sep))
   File_name = trim(adjustl(numpar%output_path))//path_sep
   File_name = trim(adjustl(File_name))//'!OUTPUT_'//trim(adjustl(matter%Name))//'_Parameters.txt'
   inquire(file=trim(adjustl(File_name)),exist=file_exists)
   if (.not.file_exists) return  ! no such file exists, nowhere to print
   inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN)
   if (.not.file_opened) then
      FN = 300
      open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
   endif
   !write(FN,'(a)') '-----------------------------------------------------------'
   write(FN,'(a,a)') 'Duration of execution of program: ', trim(adjustl(chtext))
   write(FN,'(a)') '*************************************************************'
end subroutine save_duration



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
         print*, 'Time-step of MD simulation is changed to', numpar%dt
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
         if (matter%tau_bath > 1.0d14) then  ! there is no bath, too slow to couple
            numpar%Transport = .false. ! excluded
            print*, 'Atomic thermostat is off'
         else
            numpar%Transport = .true.	 ! included
            print*, 'Atomic thermostat parameters are changed to', &
                  matter%T_bath*g_kb, matter%tau_bath
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
         if (matter%tau_bath_e > 1.0d14) then  ! there is no bath, too slow to couple
            numpar%Transport_e = .false. ! excluded
            print*, 'Electronic thermostat is off'
         else
            numpar%Transport_e = .true.	 ! included
            print*, 'Electronic thermostat parameters are changed to', &
                  matter%T_bath_e*g_kb, matter%tau_bath_e
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


subroutine act_on_comunication(read_well, given_line, given_num, numpar, matter, time)
   logical, intent(in) :: read_well ! did we read something meaningful from the comunication file?
   character(*), intent(in) :: given_line ! line read from the file
   real(8), intent(in) :: given_num  ! number read from the file
   type(Numerics_param), intent(inout) :: numpar ! all numerical parameters
   type(Solid), intent(in) :: matter ! parameters of the material
   real(8), intent(in) :: time ! current time [fs]
   integer FN, noth, lngt
   logical file_opened
   character(200) :: File_name, temp1, temp2, given_line_processed
   character(1) path_sep

   path_sep = trim(adjustl(numpar%path_sep))
   lngt = LEN(trim(adjustl(given_line)))    ! length of the line

   ! Check if the last character is TAB:
   if (given_line(lngt:lngt) == char(9) ) then  ! remove the TAB from the line
      given_line_processed = given_line(1:lngt-1)
   else
      given_line_processed = given_line
   endif

   if (read_well) then
      File_name = trim(adjustl(numpar%output_path))//path_sep
      File_name = trim(adjustl(File_name))//'!OUTPUT_'//trim(adjustl(matter%Name))//'_Parameters.txt'
      inquire(file=trim(adjustl(File_name)),opened=file_opened, number=FN)
      if (.not.file_opened) then
         FN = 300
         open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'new')
      endif

      select case(trim(adjustl(given_line_processed)))
      case ('verbose', 'VERBOSE', 'Verbose')
         print*, 'Verbose option on: XTANT will print a lot of markers for testing and debugging'
         numpar%verbose = .true.
         write(FN,'(a,f10.3,a)') 'At time instance of ', time, ' verbose option was switched on'

      case ('time', 'TIME', 'Time', 'TIme', 'TIMe', 'tIme', 'emit', 'Vremya')
         numpar%t_total = given_num ! total duration of simulation [fs]
         print*, 'Duration of simulation is changed to', given_num
         write(FN,'(a,f10.3,a,f10.3)') 'At time instance of ', time, ' duration of simulation is changed to ', given_num

      case ('MDdt', 'dtMD', 'mddt', 'dtmd', 'MDDT', 'DTMD')
         numpar%dt = given_num ! Time step for MD [fs]
         call reset_support_times(numpar)   ! above
         !numpar%halfdt = numpar%dt/2.0d0      ! dt/2, often used
         !numpar%dtsqare = numpar%dt*numpar%halfdt ! dt*dt/2, often used
         !numpar%dt3 = numpar%dt**3/6.0d0            ! dt^3/6, often used
         !numpar%dt4 = numpar%dt*numpar%dt3/8.0d0    ! dt^4/48, often used
         print*, 'Time-step of MD simulation is changed to', given_num
         write(FN,'(a,f10.3,a,f9.3)') 'At time instance of ', time, ' time-step of MD simulation is changed to ', given_num 

      case ('SAVEdt', 'savedt', 'dtsave', 'dtSAVE', 'Savedt', 'SaveDT', 'SaveDt')
         numpar%dt_save = given_num ! save data into files every 'dt_save_time' [fs]
         print*, 'Time-step of saving output files is changed to', given_num
         write(FN,'(a,f10.3,a,f9.3)') 'At time instance of ', time, ' time-step of saving output files is changed to ', given_num
      
      case ('OMP', 'omp', 'NOMP', 'nomp', 'Nomp', 'N_OMP', 'n_omp')
         ! Reset the OpenMP parallelization options:
         numpar%NOMP = given_num
         write(temp1,'(f12.3)') time
         write(temp2,'(i10)') INT(given_num)
         
#ifdef OMP_inside
         noth = OMP_GET_MAX_THREADS()   ! to chech if the function worked
         call set_OMP_number( numpar%NOMP, .true., 6, 'Reset number of threads in OpenMP to '//trim(adjustl(temp2)) )    ! below
         if ( noth /= OMP_GET_MAX_THREADS() ) then
            write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs], number of threads in OpenMP is changed to ',  trim(adjustl(temp2))
         else
            write(FN,'(a,a,a,a)') 'At time instant of ',  trim(adjustl(temp1)), '[fs]: unsuccessful attempt to change number of threads in OpenMP to ',  trim(adjustl(temp2))
            write(6,'(a)') 'Number of threads in OpenMP is unchanged: ',  trim(adjustl(temp2))
         endif
#else
         write(FN,'(a)') ' The code compiled without OpenMP, cannot set parallelization'
         write(6,'(a)') 'The code compiled without OpenMP, cannot set parallelization'
#endif

      case ('Te', 'te', 'TE') ! DO NOT USE: this option is not finished yet!
!          print*, 'Time-step of saving output files is changed to', given_num
!          write(FN,'(a,f10.3,a,f9.3)') 'At time instance of ', time, ' electronic temperature is changed to ', given_num
      case ('pulse', 'PULSE', 'Pulse') ! DO NOT USE: this option is not finished yet!
!          print*, 'Parameters of the pulse number', int(given_num), 'are changed'
!          write(FN,'(a,f10.3,a,i4,a)') 'At time instance of ', time, ' parameters of the pulse number ', int(given_num), ' are changed'
      case default
         print*, 'Could not interpret what is read from the file: ', trim(adjustl(given_line)), given_num
      end select
   endif
end subroutine act_on_comunication



subroutine set_OMP_number(NOMP, prnt, FN, lin)
   integer, intent(inout) :: NOMP  ! number of threads to be set; negative means = maximal threads available
   logical, intent(in) :: prnt  ! do we want to print out anything?
   integer, intent(in) :: FN    ! file number to print into
   character(*), intent(in), optional :: lin    ! a line to print out
   !------------------------------------
   character(10) :: temp2
   
#ifdef OMP_inside
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
   call read_file(Reason, i, read_well)  ! module "Dealing_with_files"
   if (Reason .LT. 0) then
      print*, 'No descriptor or value found in the communication file'
   else if (Reason .GT. 0) then
      print*, 'Given number interpreted as', out_num, ', it does not match the variable type'
   endif
   if (.not.read_well) then
      print*, 'Comunication format must be as follows:'
      print*, 'Two columns: 1) descriptor; 2) value'
      print*, 'Allowed descriptors: Time; dt; Save_dt; OMP'
   endif
end subroutine pars_comunications


subroutine pars_comunications_file(FN, i, out_line, out_num, Reason)
   integer, intent(in) :: FN
   integer, intent(inout) :: i
   character(*), intent(out) :: out_line
   real(8), intent(out) :: out_num
   integer, intent(out) :: Reason
   !---------------------------------
   logical :: read_well
   read_well = .false.
   out_line = ''
   out_num = 0.0d0

   read(FN, *, IOSTAT=Reason) out_line, out_num
   call read_file(Reason, i, read_well)  ! module "Dealing_with_files"
   if (Reason .LT. 0) then
      !print*, 'No descriptor or value found in the communication file'
   else if (Reason .GT. 0) then
      !print*, 'Given number interpreted as', out_num, ', it does not match the variable type'
      print*, 'Wrong format of input, could not interpret.'
      print*, 'Comunication format must be as follows:'
      print*, 'Two columns: 1) descriptor; 2) value'
      print*, 'Allowed descriptors: Time; dt; Save_dt; OMP'
   endif
end subroutine pars_comunications_file




subroutine Print_title(print_to, Scell, matter, laser, numpar)
   integer, intent(in) :: print_to ! the screen, or file
   type(Super_cell), dimension(:), intent(in) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(in) :: matter ! material parameters
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar ! all numerical parameters
   !type(TB_repulsive), dimension(:), intent(in) :: TB_Repuls  ! parameters of the repulsive part of TB
   !type(TB_Hamiltonian), dimension(:), intent(in) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   integer i 
   character(100) :: text, text1, text2, text3, starline

   starline = '*************************************************************'

   write(print_to,'(a)') trim(adjustl(starline))
   write(print_to,'(a)') '*  XTANT: X-ray-induced Thermal And Nonthermal Transitions  *'
   write(print_to,'(a)') trim(adjustl(starline))
   write(print_to,'(a)') '  A hybrid approach consisting of: '
   write(print_to,'(a)') ' (1) Monte Carlo '
   write(print_to,'(a)') ' (2) Transferable Tight Binding '
   write(print_to,'(a)') ' (3) Molecular Dynamics '
   write(print_to,'(a)') ' (4) Boltzmann collision integrals '
   write(print_to,'(a,a)') ' Applied for ', trim(adjustl(matter%Name))
   write(print_to,'(a)') trim(adjustl(starline))
   write(print_to,'(a)') ' Chemical formula of target material interpreted as: '
   do i = 1, size(matter%Atoms)
      write(text,'(f12.6)') matter%Atoms(i)%percentage
      write(text1,'(i3)') INT(matter%Atoms(i)%Z)
      write(print_to,'(a,a,a,a,a)') ' '//trim(adjustl(text)), ' of ', trim(adjustl(matter%Atoms(i)%Name)), ' (element #', trim(adjustl(text1))//')'
   enddo

   if (numpar%embed_water) then
      write(text1,'(i6)') numpar%N_water_mol
      write(print_to,'(a)') ' (Note that the material was embedded in water with # of molecules: '//trim(adjustl(text1))//')'
   endif

   write(print_to,'(a,a,a)') ' Calculations performed for the following parameters:'
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
         write(print_to,'(a,f12.5,a)') ' Absorbed fluence ' , laser(i)%F, ' [eV/atom]'
      enddo
   endif ! FEL included or not?
   
   SCL:do i = 1, size(Scell)
      select case (numpar%optic_model)
      case (1)	! within the Drude model
         write(print_to,'(a)') ' Probe-pulse is calculated within Drude model'
         write(print_to,'(a)') ' with the following parameters of the probe:'
         write(print_to,'(a, f7.1, a, f5.1, a)') ' Wavelength: ', Scell(i)%eps%l, '[nm]; Angle:', Scell(i)%eps%teta/g_Pi*(180.0d0), '[degrees]'
         write(print_to,'(a, f7.1, a)') ' Thickness of the sample: ', Scell(i)%eps%dd, ' [nm]'
         write(print_to,'(a, es12.3, es12.3)') ' Effective mass of electron and hole: ', Scell(i)%eps%me_eff, Scell(i)%eps%mh_eff
         write(print_to,'(a, es12.3, es12.3)') ' Effective scattering time of electron and of hole: ', Scell(i)%eps%tau_e, Scell(i)%eps%tau_h
      case (2:3)	! Trani model
         write(print_to,'(a)') ' Probe-pulse is calculated within Trani approach  [PRB 72, 075423 (2005)]'
         write(print_to,'(a)') ' with the following parameters of the probe:'
         write(print_to,'(a, f7.1, a, f5.1, a)') ' Wavelength: ', Scell(i)%eps%l, ' [nm]; Angle:', Scell(i)%eps%teta/g_Pi*(180.0d0), '    [degrees]'
         write(print_to,'(a, f7.1, a)') ' Thickness of the sample: ', Scell(i)%eps%dd, ' [nm]'
         if (numpar%optic_model .EQ. 2) then
            write(text1, '(i10)') numpar%ixm
            write(text2, '(i10)') numpar%iym
            write(text3, '(i10)') numpar%izm
             if (allocated(numpar%k_grid)) then
               write(print_to,'(a,a,a,a,a,a)') ' Number of k-points (on user-defined grid): ', trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
            else
               write(print_to,'(a,a,a,a,a,a)') ' Number of k-points (on Monkhorst Pack grid): ', trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
            endif
         else
            write(print_to,'(a)') ' Calculations are performed for Gamma-point'
         endif
      case default ! no optical coefficients needed
         write(print_to,'(a)') ' No probe-pulse is calculated'
      end select
      if (Scell(i)%eps%all_w) then
         if (Scell(i)%eps%KK) then
            write(print_to,'(a)') ' The spectrum is calculated via Im(CDF) using Kramers Kronig relations'
         else
            write(print_to,'(a)') ' The spectrum is calculated directly for both, Re(CDF) and Im(CDF)'
         endif
      endif
   enddo SCL
   write(print_to,'(a)') trim(adjustl(starline))
   write(text, '(f15.5)') numpar%t_total
   write(print_to,'(a,a,a)') ' Duration of modelling ' , trim(adjustl(text)), ' [fs]'

   if (matter%cell_x*matter%cell_y*matter%cell_z .LE. 0) then
      write(print_to,'(a)') ' TBMD part is switched off, only MC modelling is performed'
   else
      write(print_to,'(a)') ' Tight Binding parametrization schemes used are:'         
      write(print_to,'(a,a)') ' Hamiltonian:      ', trim(adjustl(Scell(1)%TB_Hamil(1,1)%Param))
      write(print_to,'(a,a)') ' Repulsive energy: ', trim(adjustl(Scell(1)%TB_Repuls(1,1)%Param))
      
      ASSOCIATE (ARRAY => Scell(1)%TB_Hamil(1,1)) ! this is the sintax we have to use to check the class of defined types
         select type(ARRAY)
         type is (TB_H_DFTB) ! TB parametrization
            write(print_to,'(a,a)') ' With skf-parametrization: (', trim(adjustl(ARRAY%param_name))//')'
            select case (numpar%N_basis_size)
            case (0)
               text = 's'
            case (1)    ! sp3
               text = 'sp3'
            case default    ! sp3d5
               text = 'sp3d5'
            endselect
            write(print_to,'(a,a)') ' With the basis set: ', trim(adjustl(text))
          type is (TB_H_3TB) ! TB parametrization
            select case (numpar%N_basis_size)
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
            select case (numpar%N_basis_size)
            case (0)
               text = 's'
            case (1)    ! sp3
               text = 'sp3'
            case default    ! sp3d5
               text = 'sp3d5'
            endselect
            write(print_to,'(a,a)') ' With the basis set: ', trim(adjustl(text))
          type is (TB_H_xTB) ! TB parametrization
            select case (numpar%N_basis_size)
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
         write(print_to,'(a,a)') ' van der Waals energy: ', trim(adjustl(Scell(1)%TB_Waals(1,1)%Param))
      else !For this material vdW class is undefined
         write(print_to,'(a,a)') ' No van der Waals potential was defined'
      endif
      if (allocated(Scell(1)%TB_Coul)) then ! if we have Coulomb potential defined
         write(print_to,'(a,a)') ' Coulomb energy: ', trim(adjustl(Scell(1)%TB_Coul(1,1)%Param))
      else !For this material vdW class is undefined
         write(print_to,'(a,a)') ' No Coulomb potential was defined or unballanced charge allowed'
      endif
      if (allocated(Scell(1)%TB_Expwall)) then ! if we have exponential wall potential defined
         write(print_to,'(a,a)') ' Exponential wall energy: ', trim(adjustl(Scell(1)%TB_Expwall(1,1)%Param))
      else !For this material exponential wall class is undefined
         write(print_to,'(a,a)') ' No exponential wall potential was defined for close interatomic distances'
      endif

      ! What kind of supercell is used:
      select case (numpar%save_files_used)
      case default ! constructed from unit cells
         write(text1, '(i10)') matter%cell_x
         write(text2, '(i10)') matter%cell_y
         write(text3, '(i10)') matter%cell_z
         write(print_to,'(a,a,a,a,a,a)') ' Super-cell size in unit-cells: ', trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))
      case (1)  ! save files are used
         write(print_to,'(a)') ' Super-cell parameters are set in SAVE files'
      case (2)  ! path coordinate
         write(print_to,'(a)') ' Coordinate path calculations are performed, defined by PATH files'
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
   else  ! in artificial cases, the dnsity may be wild:
      write(print_to,'(a,es25.3,a)') ' Density of the material: ', matter%dens,' [g/cm^3]'
      write(print_to,'(a,es12.3,a)') ' The used atomic density (used in MC cross sections): ', matter%At_dens, ' [1/cm^3]'
   endif
   write(print_to,'(a)') ' The following numerical parameters are used:'
   write(print_to,'(a,i6)') ' Number of iterations in the MC module: ', numpar%NMC
   if (numpar%do_elastic_MC) then ! allow elastic scattering of electrons on atoms within MC module
      write(print_to,'(a)') ' Elastic high-energy-electron scattering is included in MC via Motts cross section'
   else
      write(print_to,'(a)') ' Elastic high-energy-electron scattering is excluded in MC'
   endif

#ifdef OMP_inside
   write(print_to,'(a,i6)') ' Number of threads for OPENMP: ', numpar%NOMP
#else ! if you set to use OpenMP in compiling: 'make OMP=no'
   write(print_to,'(a)') ' The code is compiled without OPENMP'
#endif

   AT_MOVE:if (numpar%do_atoms) then ! atoms are moving:
      select case (numpar%MD_algo)
      case (1)
         write(print_to,'(a)') ' MD algorithm used: Yoshida (4th order)'
      case (2)
         write(print_to,'(a)') ' MD algorithm used: Martyna (4th order)'
      case default
         write(print_to,'(a)') ' MD algorithm used: velocity Verlet (2th order)'
      endselect
      if (allocated(numpar%dt_MD_reset_grid)) then
         write(print_to,'(a)') ' Time step in MD simulations is set from file: '//trim(adjustl(numpar%MD_step_grid_file))
      else
         write(print_to,'(a,f9.3,a)') ' Time step in MD simulations: ', numpar%dt,' [fs]'
      endif
      write(print_to,'(a,f9.3,a)') ' Output data are saved every: ', numpar%dt_save,' [fs]'
      if (numpar%p_const) then	! P=const
         write(print_to,'(a)') ' Constant pressure simulation (Parrinello-Rahman scheme, NPH) '
         write(print_to,'(a,f9.1,a)') ' External pressure: ', matter%p_ext,' [Pa]'
      else ! V=const
         write(print_to,'(a)') ' Constant volume simulation (NVE)'
      endif
      write(print_to,'(a)') ' Scheme used for low-energy electrons modelling: '
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
         if (numpar%tau_fe < 1e6) then
            write(text1, '(f13.6)') numpar%tau_fe
         else
            write(text1, '(es16.6)') numpar%tau_fe
         endif
         write(print_to,'(a)') ' Relaxation-time approximation for electron thermalization'
         write(print_to,'(a)') ' with the characteristic time '//trim(adjustl(text1))//' [fs]'
      end select
   else AT_MOVE
      write(print_to,'(a)') ' Atoms were FROZEN instead of moving in MD!'
   endif AT_MOVE

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
      write(print_to,'(a, f7.1, a)') ' switched on at: ', numpar%t_NA, ' [fs]'
      write(print_to,'(a, f7.1, a)') ' with the acceptance window: ', numpar%acc_window, ' [eV]'
      write(print_to,'(a, f8.5, a)') ' degeneracy tollerance: ', numpar%degeneracy_eV, ' [eV]'
      write(print_to,'(a, f8.5)') ' and scaling factor of: ', numpar%M2_scaling
   endif

   if (numpar%do_kappa) then
      write(print_to,'(a)') ' Calculation of electronic heat conductivity is included'
   else
      write(print_to,'(a)') ' No calculation of electronic heat conductivity'
   endif

   if (numpar%do_cool) then
      write(print_to,'(a)') ' Quenching of atoms for resolidification is included'
      write(print_to,'(a, f7.1, a, f7.1, a)') ' Starting at: ', numpar%at_cool_start, ' [fs], with the step of: ', numpar%at_cool_dt, ' [fs]'
   else
      write(print_to,'(a)') ' No quenching of atoms for resolidification'
   endif


   if (allocated(numpar%El_bath_reset_grid)) then
      write(print_to,'(a)') ' Berendsen thermostat is used for electnros'
      write(print_to,'(a)') ' with parameters set in the file: '//trim(adjustl(numpar%El_bath_step_grid_file))
   elseif (g_numpar%Transport_e) then ! for electrons
      write(text,'(f10.1)') matter%T_bath_e*g_kb
      write(text1,'(f10.1)') matter%tau_bath_e
      write(print_to,'(a)') ' Berendsen thermostat is used for electnros'
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

   write(print_to,'(a, f7.1, a)') ' Electron energy cut-off, separating high-energy- from low-energy-electrons: ', numpar%E_cut, ' [eV]'
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
   if (numpar%save_Ei) then
      write(print_to,'(a)') ' Calculated energy levels are saved in output file'
   endif
   if (numpar%save_DOS) then
      write(print_to,'(a, f7.5, a)') ' Calculated DOS is saved in output file; smearing used: ', numpar%Smear_DOS, ' [eV]'
      select case (ABS(numpar%optic_model))	! use multiple k-points, or only gamma
         case (2)	! multiple k points
            write(text1, '(i10)') numpar%ixm
            write(text2, '(i10)') numpar%iym
            write(text3, '(i10)') numpar%izm
            if (allocated(numpar%k_grid)) then
               write(print_to,'(a,a,a,a,a,a)') ' It is calculated on the user-defined grid for points: ', trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))   
            else
               write(print_to,'(a,a,a,a,a,a)') ' It is calculated on Monkhorst Pack grid for points: ', trim(adjustl(text1)),'x',trim(adjustl(text2)),'x',trim(adjustl(text3))   
            endif
         case default	! gamma point
            write(print_to,'(a)') ' It is calculated at the Gamma point'
         end select
   endif
   if (numpar%save_fe) then
      write(print_to,'(a)') ' Calculated electron electron distributions are saved in output file'
   endif
   if (numpar%save_fe_grid) then
      write(print_to,'(a)') ' Calculated electron electron distribution on grid is saved in output file'
   endif
   if (numpar%save_PCF) then
      write(print_to,'(a)') ' Calculated atomic pair correlation functions are saved in output file'
   endif
   if (numpar%save_XYZ) then
      write(print_to,'(a)') ' Calculated atomic positions in XYZ-format are saved in output file'
   endif
   if (numpar%save_CIF) then
      write(print_to,'(a)') ' Calculated atomic positions in CIF-format are saved in output file'
   endif
   if (numpar%save_NN) then
      write(text1, '(f6.2)') numpar%NN_radius
      write(print_to,'(a,a,a)') ' Nearest neighbors numbers within the radius of ', trim(adjustl(text1)), ' [A] are saved'
   endif

9999   write(print_to,'(a)') trim(adjustl(starline))
end subroutine Print_title

END MODULE Dealing_with_output_files
