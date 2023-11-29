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
! This module contains subroutines to set initial conditions:

MODULE Initial_configuration
use Universal_constants
use Objects
use Algebra_tools, only : Det_3x3, Reciproc
use Dealing_with_files, only : Count_lines_in_file, Count_columns_in_file, close_file, read_file, get_file_extension, number_of_columns
use Atomic_tools, only : Coordinates_abs_to_rel, remove_angular_momentum, get_fragments_indices, remove_momentum, &
                  Get_random_velocity, check_periodic_boundaries, Make_free_surfaces, Coordinates_abs_to_rel_single, &
                  shortest_distance, velocities_rel_to_abs, velocities_abs_to_rel, Coordinates_rel_to_abs, &
                  check_periodic_boundaries_single, Coordinates_rel_to_abs_single, deflect_velosity
use TB, only : get_DOS_masks, get_Hamilonian_and_E
use Electron_tools, only : get_glob_energy
use Dealing_with_BOP, only : m_repulsive, m_N_BOP_rep_grid
use ZBL_potential, only : ZBL_pot
use TB_xTB, only : identify_xTB_orbitals_per_atom
use Little_subroutines, only : linear_interpolation, Find_in_array_monoton, deallocate_array
use Dealing_with_eXYZ, only : interpret_XYZ_comment_line
use Periodic_table, only : Decompose_compound
use Read_input_data, only : m_Atomic_parameters, m_dashline
use Dealing_with_POSCAR, only : read_POSCAR, get_KOA_from_element
use Dealing_with_mol2, only : read_mol2

implicit none
PRIVATE


real(8) :: m_H2O_dist, m_H2O_theta, m_one_third
parameter (m_H2O_dist = 0.943d0)    ! [A] distance between H and O atoms in H2O molecule
parameter (m_H2O_theta = 106.0d0 * g_Pi/180.0d0)   ! [deg] H-O-H angle in H2O molecule
parameter (m_one_third = 1.0d0/3.0d0)



public :: create_BOP_repulsive, set_initial_configuration


 contains



subroutine create_BOP_repulsive(Scell, matter, numpar, TB_Repuls, i, j, Folder_name, path_sep, Name1, Name2, bond_length_in, Elem1, Elem2, Err)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(TB_Rep_BOP), dimension(:,:), intent(inout) ::  TB_Repuls    ! parameters of the repulsive potential
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   character(*), intent(in) :: Folder_name    ! directory where to find BOP parameters
   character(1), intent(in) :: path_sep
   character(*), intent(in) :: Name1, Name2 ! element names
   real(8), intent(in) :: bond_length_in   ! [A] bond length for dimer
   real(8), intent(in) :: Elem1, Elem2  ! atomic numbers of the two elements we need the parameters for
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------
   character(300) :: File_name, Error_descript, chtemp(2)
   integer :: FN_BL, k, NSC, n1, icur
   real(8) :: r_start, r_stop, dr, supcesize, Pot_shift, d_bond, ZBL_length, TB_d, ZBL_d, bond_length
   real(8), dimension(m_N_BOP_rep_grid) :: Ref_Pot, V_rep
   integer, dimension(:), allocatable :: int_vec

   ! Region around the bond length where we smoothen the potentials:
   d_bond = 0.3d0   ! [A]
   ! Shifted bond length for fitting:
   bond_length = bond_length_in + d_bond    ! [A]

   ! Create repulsive potential:
   if (.not.allocated(TB_Repuls(i,j)%R)) then
      allocate(TB_Repuls(i,j)%R(m_N_BOP_rep_grid))
      ! Create grid:
      r_start = max(0.25d0, min(bond_length-0.5d0, 1.0d0) ) ! start of the grid [A]
      r_stop = bond_length + d_bond ! end of the grid [A]
      dr = (r_stop - r_start)/dble(m_N_BOP_rep_grid-1)   ! step set to have fixed number of points equal to m_N_BOP_rep_grid
      ! Save the grid:
      TB_Repuls(i,j)%R(1) = r_start ! [A]
      do k = 2, m_N_BOP_rep_grid
         TB_Repuls(i,j)%R(k) = TB_Repuls(i,j)%R(k-1) + dr   ! [A]
      enddo
   endif
   if (.not.allocated(TB_Repuls(i,j)%V_rep)) allocate(TB_Repuls(i,j)%V_rep(m_N_BOP_rep_grid))

   ! Get the reference potential (ZBL):
   do k = 1, m_N_BOP_rep_grid
      Ref_Pot(k) =  ZBL_pot(Elem1, Elem2, TB_Repuls(i,j)%R(k))  ! module "ZBL_potential"
!       print*, k, TB_Repuls(i,j)%R(k), Ref_Pot(k)
   enddo

   ! Set the equilibrium distance between the dimer atoms to get the correct bond length:
   supcesize = 10.0d0 * TB_Repuls(i,j)%R(m_N_BOP_rep_grid)  ! large supercell size to exclude periodicity
   Scell(1)%supce = RESHAPE( (/ supcesize, 0.0d0, 0.0d0,  &
                                0.0d0, supcesize, 0.0d0,  &
                                0.0d0, 0.0d0, supcesize /), (/3,3/) )
   Scell(1)%supce0 = Scell(1)%supce

   ! Place dimer along X axis at the distance of bond length:
   allocate(Scell(1)%MDatoms(2))    ! dimer
   Scell(1)%Na = 2
   Scell(1)%Ne = SUM(matter%Atoms(:)%NVB*matter%Atoms(:)%percentage)/SUM(matter%Atoms(:)%percentage)*Scell(1)%Na

   if (numpar%verbose) then
      write(*, '(a)', advance='no') 'Number of valence electrons: '
      allocate(int_vec(size(matter%Atoms(:)%NVB)), source=matter%Atoms(:)%NVB)
      write(*,*) int_vec
      write(*,*) ' (total: ', Scell(1)%Ne, ')'
   endif

   Scell(1)%Ne_low = Scell(1)%Ne ! at the start, all electrons are low-energy
   Scell(1)%Ne_high = 0.0d0 ! no high-energy electrons at the start
   Scell(1)%Ne_emit = 0.0d0 ! no emitted electrons at the start
   ! Allocate arrays for calculation of hamiltonian and related stuff:
   allocate(Scell(1)%Near_neighbor_list(Scell(1)%Na,Scell(1)%Na))  ! nearest neighbors
   allocate(Scell(1)%Near_neighbor_dist(Scell(1)%Na,Scell(1)%Na,4))  ! [A] distances
   allocate(Scell(1)%Near_neighbor_dist_s(Scell(1)%Na,Scell(1)%Na,3)) ! relative dist.
   allocate(Scell(1)%Near_neighbor_size(Scell(1)%Na)) ! how many nearest neighbours
   allocate(Scell(1)%Near_neighbors_user(Scell(1)%Na))
   ASSOCIATE (ARRAY => Scell(i)%TB_Hamil(:,:))
   select type(ARRAY)
   type is (TB_H_BOP)   ! it can be various basis sets:
     select case (numpar%N_basis_size)    ! find which one is used now:
     case (0)    ! s
        n1 = 1.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
     case (1)    ! sp3
        n1 = 4.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
     case default    ! sp3d5
        n1 = 9.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
     endselect
   endselect
   END ASSOCIATE
   allocate(Scell(1)%Ha(n1,n1))   ! hamiltonian size
   allocate(Scell(1)%Ha0(n1,n1)) ! hamiltonian0 size
   allocate(Scell(1)%H_non(n1,n1))	! nondiagonalized Hamiltonian
   allocate(Scell(1)%H_non0(n1,n1))	! nondiagonalized Hamiltonian
   allocate(Scell(1)%Ei(n1))  ! energy levels, eigenvalues of the hamiltonian matrix
   allocate(Scell(1)%Ei0(n1))  ! energy levels0, eigenvalues of the hamiltonian matrix
   allocate(Scell(1)%Aij(n1,n1))	! coefficients used for forces in TB
   allocate(Scell(1)%fe(size(Scell(1)%Ei))) ! electron distribution function (Fermi-function)
   allocate(Scell(1)%fe_eq(size(Scell(1)%Ei))) ! equivalent electron distribution function (Fermi-function)
   if (numpar%do_kappa) then
      allocate(Scell(1)%I_ij(size(Scell(1)%Ei))) ! electron-ion collision integral
      allocate(Scell(1)%Ce_i(size(Scell(1)%Ei))) ! electron-energy resolved heat capacity
   endif
   Scell(1)%MDatoms(1)%KOA = i
   Scell(1)%MDatoms(2)%KOA = j
   Scell(1)%MDatoms(1)%S(:) = 0.1d0
   !Scell(1)%MDatoms(2)%S(1) = 0.2d0
   !Scell(1)%MDatoms(2)%S(1) = Scell(1)%MDatoms(1)%S(1) + TB_Repuls(i,j)%R(m_N_BOP_rep_grid) / supcesize
   Scell(1)%MDatoms(2)%S(1) = Scell(1)%MDatoms(1)%S(1) + bond_length / supcesize
   Scell(1)%MDatoms(2)%S(2:3) = 0.1d0
   call Det_3x3(Scell(1)%supce, Scell(1)%V)
   call Coordinates_rel_to_abs(Scell, 1, if_old=.true.)    ! from the module "Atomic_tools"
   call get_DOS_masks(Scell, matter, numpar)  ! module "TB"

   ! Get the energy shift:
   ! Contruct TB Hamiltonian, diagonalize to get energy levels, get forces for atoms and supercell:
   call get_Hamilonian_and_E(Scell, numpar, matter, 1, Err, 0.0d0) ! module "TB"
   ! Get global energy of the system at the beginning:
   call get_glob_energy(Scell, matter) ! module "Electron_tools"

   ! Energy shift:
   !Pot_shift = -Ref_Pot(m_N_BOP_rep_grid) + Scell(1)%nrg%El_low   ! [eV]
   Pot_shift = -ZBL_pot(Elem1, Elem2, bond_length) + Scell(1)%nrg%El_low   ! [eV]
   ! Shift potential accordingly, to produce correct minimum at bond length:
   Ref_Pot = Ref_Pot + Pot_shift

   ! Also define TB potential at the point of bond length + d:
   Scell(1)%MDatoms(2)%S(1) = Scell(1)%MDatoms(1)%S(1) + (bond_length+d_bond) / supcesize
   call Det_3x3(Scell(1)%supce, Scell(1)%V)
   call Coordinates_rel_to_abs(Scell, 1, if_old=.true.)    ! from the module "Atomic_tools"
   call get_DOS_masks(Scell, matter, numpar)  ! module "TB"
   call get_Hamilonian_and_E(Scell, numpar, matter, 1, Err, 0.0d0) ! module "TB"
   call get_glob_energy(Scell, matter) ! module "Electron_tools"
   TB_d = Scell(1)%nrg%El_low   ! [eV]

   ! Also define ZBL potential at the point of bond length - d:
   call Find_in_array_monoton(abs(Ref_Pot), abs(TB_d), icur) ! module "Little_subroutines"
   call linear_interpolation(Ref_Pot, TB_Repuls(i,j)%R, TB_d, ZBL_length, icur) ! module "Little_subroutines"

!    do k = 1, m_N_BOP_rep_grid
!       print*, k, TB_Repuls(i,j)%R(k), Ref_Pot(k)
!    enddo
!     print*, 'ZBL_length', ZBL_length, TB_d, icur, Ref_Pot(icur), TB_Repuls(i,j)%R(icur)


   ! Calculate the repulsive term, such that total potential equals the referenced one:
   do k = 1, m_N_BOP_rep_grid
      ! Set the distance according to grid point for repulsive potential:
      Scell(1)%MDatoms(2)%S(1) = Scell(1)%MDatoms(1)%S(1) + TB_Repuls(i,j)%R(k) / supcesize
      call Coordinates_rel_to_abs(Scell, 1, if_old=.true.)    ! from the module "Atomic_tools"
      call get_Hamilonian_and_E(Scell, numpar, matter, 1, Err, 0.0d0) ! module "TB"
      ! Get global energy of the system at the beginning:
      call get_glob_energy(Scell, matter) ! module "Electron_tools"
      ! Set the repulsive potential:
      if ( TB_Repuls(i,j)%R(k) <= ZBL_length ) then
         V_rep(k) = Ref_Pot(k) - Scell(1)%nrg%El_low    ! [eV]
         ZBL_d = V_rep(k)
      else
         ZBL_d = 1.0d0 / (1.0d0/(Ref_Pot(k)-TB_d) + 1.0d0/(Scell(1)%nrg%El_low-TB_d))
         V_rep(k) = (TB_d + ZBL_d) - Scell(1)%nrg%El_low   ! [eV]
      endif
!       write(*,'(i3,f,f,es,es,es,f,f)') k, TB_Repuls(i,j)%R(k), V_rep(k),  (1.0d0/(Ref_Pot(k)-TB_d) + 1.0d0/(Scell(1)%nrg%El_low-TB_d)), ZBL_d, TB_d - Scell(1)%nrg%El_low, TB_Repuls(i,j)%R(k), ZBL_length
   enddo
   TB_Repuls(i,j)%V_rep = V_rep ! save it
!    print*, 'Pot_shift', Pot_shift, supcesize, sqrt(SUM((Scell(1)%MDatoms(1)%R(:)-Scell(1)%MDatoms(2)%R(:))**2))

   ! Restore the parameters:
   deallocate(Scell(1)%MDatoms)
   deallocate(Scell(1)%Near_neighbor_list)  ! nearest neighbors
   deallocate(Scell(1)%Near_neighbor_dist)  ! [A] distances
   deallocate(Scell(1)%Near_neighbor_dist_s) ! relative dist.
   deallocate(Scell(1)%Near_neighbor_size) ! how many nearest neighbours
   deallocate(Scell(1)%Near_neighbors_user)
   deallocate(Scell(1)%Ha)
   deallocate(Scell(1)%Ha0)
   deallocate(Scell(1)%H_non)
   deallocate(Scell(1)%H_non0)
   deallocate(Scell(1)%Ei)
   deallocate(Scell(1)%Ei0)
   deallocate(Scell(1)%Aij)
   deallocate(Scell(1)%fe)
   deallocate(Scell(1)%fe_eq)
   deallocate(Scell(1)%G_ei_partial)
   deallocate(Scell(1)%Ce_part)
   call deallocate_array(Scell(1)%I_ij)      ! module "Little_subroutines"
   call deallocate_array(Scell(1)%Norm_WF)   ! module "Little_subroutines"
   call deallocate_array(Scell(1)%Ce_i)      ! module "Little_subroutines"
   call deallocate_array(Scell(1)%kappa_e_part)   ! module "Little_subroutines"
   deallocate(numpar%mask_DOS)

!    pause 'create_BOP_repulsive'

   if (j /= i) then ! and the lower triangle
      if (.not.allocated( TB_Repuls(j,i)%R)) allocate(TB_Repuls(j,i)%R(m_N_BOP_rep_grid))
      if (.not.allocated( TB_Repuls(j,i)%V_rep)) allocate(TB_Repuls(j,i)%V_rep(m_N_BOP_rep_grid))
      TB_Repuls(j,i)%R(:) = TB_Repuls(i,j)%R(:)
      TB_Repuls(j,i)%V_rep(:) = TB_Repuls(i,j)%V_rep(:)
   endif

   ! File with repulsive BOP potential to be created:
   File_name = trim(adjustl(Folder_name))//path_sep// &
      trim(adjustl(Name1))//'_'//trim(adjustl(Name2))//trim(adjustl(m_repulsive))   ! file with repulsive BOP parameters
   FN_BL = 113
   open(UNIT=FN_BL, FILE = trim(adjustl(File_name)))
   ! Write into the file:
   do k = 1, m_N_BOP_rep_grid    ! for all grid points
      write(FN_BL,'(f24.16, es24.16)') TB_Repuls(i,j)%R(k), TB_Repuls(i,j)%V_rep(k)
   enddo
   call close_file('close', FN=FN_BL) ! module "Dealing_with_files"

3415 continue
end subroutine create_BOP_repulsive



subroutine print_message_about_input_files(Error_descript)
   character(200), intent(in), optional :: Error_descript
   write(*,'(a)') trim(adjustl(m_dashline))
   if (present(Error_descript)) write(*,'(a)') trim(adjustl(Error_descript))
   write(*,'(a)') 'No file with supercell was found. The file(s) must be set in one of the formats:'
   write(*,'(a)') '1) Path-coordinates (in internal XTANT SAVE-files format: PHASE_1_supercell.dat and PHASE_2_supercell.dat)'
   write(*,'(a)') '2) SAVE-files  (internal XTANT SAVE-files format: SAVE_atoms.dat and SAVE_supercell.dat)'
   write(*,'(a)') '3) XYZ format (default name: Cell.xyz)'
   write(*,'(a)') '4) POSCAR file (default name: Cell.poscar)'
   write(*,'(a)') '5) mol2 file (default name: Cell.mol2)'
   write(*,'(a)') '6) Unit-cell coordinates (files Unit_cell_equilibrium.txt and Unit_cell_atom_relative_coordinates.txt)'
   write(*,'(a)') trim(adjustl(m_dashline))
end subroutine print_message_about_input_files



subroutine set_initial_configuration(Scell, matter, numpar, laser, MC, Err)
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(MC_data), dimension(:), allocatable, intent(inout) :: MC ! all MC parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !========================================================
   integer i, Nsc, Natoms, FN, FN2, Reason, count_lines, N, j, k, n1, FN3, FN4, FN_XYZ, FN_POSCAR, FN_mol2
   character(200) :: File_name, File_name2, Error_descript, File_name_S1, File_name_S2, File_name_XYZ, File_name_POSCAR, &
                     File_name_mol2, Cell_filename
   character(10) :: file_extension
   logical :: file_exist, file_opened, read_well, file_exist_1, file_exist_2, XYZ_file_exists, POSCAR_file_exists, mol2_file_exists
   real(8) RN, temp, Mass, V2, Ta

   Nsc = 1 !in the present version of the code, there is always only one super-cell

   ! If file with BOP repulsive potential does not exist, create it:
   if (numpar%create_BOP_repulse) then
      do i = 1, size(Scell(NSC)%TB_Repuls,1)
         do j = 1, size(Scell(NSC)%TB_Repuls,2)
            ASSOCIATE (TB_Repuls => Scell(NSC)%TB_Repuls)
               select type(TB_Repuls)
               type is (TB_Rep_BOP)
                  call create_BOP_repulsive(Scell, matter, numpar, TB_Repuls, i, j, numpar%BOP_Folder_name, numpar%path_sep, &
                     trim(adjustl(matter%Atoms(i)%Name)), trim(adjustl(matter%Atoms(j)%Name)), &
                     numpar%BOP_bond_length, matter%Atoms(i)%Z, matter%Atoms(j)%Z, Err) ! above
               endselect
            END ASSOCIATE
         enddo ! j
      enddo ! i
   endif

   MD:if (matter%cell_x*matter%cell_y*matter%cell_z .GT. 0) then  ! only then there is something to do with atoms:
      if (.not.allocated(Scell)) allocate(Scell(Nsc))
      ALL_SC: do i = 1, Nsc
         
         numpar%do_path_coordinate = .false. ! to check files with phase 1 and 2
         
         ! Supercell vectors:

         ! Check if the user provided the filename with coordinates and supercell:
         call get_file_extension(trim(adjustl(numpar%Cell_filename)), file_extension)  ! module "Dealing_with_files"
         XYZ_file_exists = .false.  ! to start with
         POSCAR_file_exists = .false.  ! to start with
         mol2_file_exists = .false.  ! to start with
         Cell_filename = ''   ! default
         if (LEN(trim(adjustl(file_extension))) > 0) then ! no filename was provided by the user, use defaults
            Cell_filename = trim(adjustl(numpar%Cell_filename))
            select case(trim(adjustl(file_extension)))
            case ('XYZ', 'XYz', 'Xyz', 'xyz')
               XYZ_file_exists = .true.
            case ('POSCAR', 'Poscar', 'poscar', 'PosCar')
               POSCAR_file_exists = .true.
            case ('mol2', 'MOL2', 'Mol2')
               mol2_file_exists = .true.
            case default
               write(*,'(a)') 'Extension of provided file '//trim(adjustl(numpar%Cell_filename))//' not supported; using default instead'
               XYZ_file_exists = .false.  ! to start with
               POSCAR_file_exists = .false.  ! to start with
               mol2_file_exists = .false.  ! to start with
               Cell_filename = ''   ! default
            end select
         endif

         ! Check if there is extended XYZ-format with the unit/super-cell:
         if (.not.XYZ_file_exists) then ! there is no name given, use default
            Cell_filename = 'Cell.xyz'   ! default name
         else
            Cell_filename = trim(adjustl(numpar%Cell_filename))
         endif
         FN_XYZ = 9004
         write(File_name_XYZ, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
                                         trim(adjustl(Cell_filename))
         inquire(file=trim(adjustl(File_name_XYZ)),exist=XYZ_file_exists)

         ! Check if there is POSCAR-format with the unit/super-cell:
         if (.not.POSCAR_file_exists) then ! there is no name given, use default
            Cell_filename = 'Cell.poscar'   ! default name
         else
            Cell_filename = trim(adjustl(numpar%Cell_filename))
         endif
         FN_POSCAR = 9005
         write(File_name_POSCAR, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
                                            trim(adjustl(Cell_filename))
         inquire(file=trim(adjustl(File_name_POSCAR)),exist=POSCAR_file_exists)

         ! Check if there is mol2-format with the unit/super-cell:
         if (.not.mol2_file_exists) then ! there is no name given, use default
            Cell_filename = 'Cell.mol2'   ! default name
         else
            Cell_filename = trim(adjustl(numpar%Cell_filename))
         endif
         FN_mol2 = 9006
         write(File_name_mol2, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
                                            trim(adjustl(Cell_filename))
         inquire(file=trim(adjustl(File_name_mol2)),exist=mol2_file_exists)


         ! Check if user set to calculate along path coordinate:
         FN3 = 9002
         FN4 = 9003
         write(File_name_S1, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'PHASE_1_supercell.dat'
         inquire(file=trim(adjustl(File_name_S1)),exist=file_exist_1)
         write(File_name_S2, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'PHASE_2_supercell.dat'
         inquire(file=trim(adjustl(File_name_S2)),exist=file_exist_2)
         numpar%do_path_coordinate = (file_exist_1 .and. file_exist_2)
         
         ! Check if there is file with Supercell vectors:
         FN = 9000
         write(File_name, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'SAVE_supercell.dat'
         inquire(file=trim(adjustl(File_name)),exist=file_exist)
         
         ! Select among different possibilities to set the super-cell:
         ! In the following priorities:
         ! 1) Path-coordinates  (in internal XTANT SAVE-files format)
         ! 2) SAVE-files  (internal XTANT SAVE-files format)
         ! 3) Cell-file  (extended XYZ format)
         ! 4) POSCAR-file (vasp format)
         ! 5) mol2-file  (SYBYL molecules format; severely restricted here to the bare minimum!)
         ! 6) unit-cell coordinates  (old internal XTANT format)

         SAVED_SUPCELL:if (numpar%do_path_coordinate) then ! read phase 1 and 2 supercells:
            
            ! Read phase 1 parameters:
            inquire(file=trim(adjustl(File_name_S1)),exist=file_exist)
            INPUT_PHASE_1:if (file_exist) then
               open(UNIT=FN3, FILE = trim(adjustl(File_name_S1)), status = 'old', action='read')
               inquire(file=trim(adjustl(File_name_S1)),opened=file_opened)
               if (.not.file_opened) then
                  Error_descript = 'File '//trim(adjustl(File_name_S1))//' could not be opened, the program terminates'
                  call Save_error_details(Err, 2, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3416
               endif

               ! Read the supercell parameters of the initial phase:
               call get_supercell_vectors(FN3, File_name_S1, Scell, i, 1, matter, Err, ind=0) ! see below

            else INPUT_PHASE_1
               write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name_S1))//' could not be found, the program terminates'
               call Save_error_details(Err, 1, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3416
            endif INPUT_PHASE_1
            
            ! Read phase 2 parameters:
            inquire(file=trim(adjustl(File_name_S2)),exist=file_exist)
            INPUT_PHASE_2:if (file_exist) then
               open(UNIT=FN4, FILE = trim(adjustl(File_name_S2)), status = 'old', action='read')
               inquire(file=trim(adjustl(File_name_S2)),opened=file_opened)
               if (.not.file_opened) then
                  Error_descript = 'File '//trim(adjustl(File_name_S2))//' could not be opened, the program terminates'
                  call Save_error_details(Err, 2, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3416
               endif
               
               ! Read the supercell parameters of the final phase:
               call get_supercell_vectors(FN4, File_name_S2, Scell, i, 1, matter, Err, ind=1) ! see below
               
            else INPUT_PHASE_2
               write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name_S2))//' could not be found, the program terminates'
               call Save_error_details(Err, 1, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3416
            endif INPUT_PHASE_2
            
            numpar%do_path_coordinate = .false. ! to check again files with atomic coordinates below
            inquire(file=trim(adjustl(File_name_S1)),opened=file_opened)
            if (file_opened) close (FN3)
            inquire(file=trim(adjustl(File_name_S2)),opened=file_opened)
            if (file_opened) close (FN4)
         
         !----------------------------
         elseif (file_exist) then SAVED_SUPCELL  ! read from this file with transient Super cell:
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
            INPUT_SUPCELL:if (file_exist) then
               open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
               inquire(file=trim(adjustl(File_name)),opened=file_opened)
               if (.not.file_opened) then
                  Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
                  call Save_error_details(Err, 2, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3416
               endif

               call get_supercell_vectors(FN, File_name, Scell, i, 1, matter, Err) ! see below

            else INPUT_SUPCELL
               write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, the program terminates'
               call Save_error_details(Err, 1, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3416
            endif INPUT_SUPCELL

         !----------------------------
         elseif (XYZ_file_exists) then SAVED_SUPCELL  ! read from XYZ file:
            ! Will read it together with atomic coordinates from the same file below

         !----------------------------
         elseif (POSCAR_file_exists) then SAVED_SUPCELL  ! read from POSCAR file:
            ! Will read it together with atomic coordinates from the same file below

         !----------------------------
         elseif (mol2_file_exists) then SAVED_SUPCELL  ! read from mol2 file:
            ! Will read it together with atomic coordinates from the same file below

         !----------------------------
         else SAVED_SUPCELL   ! no supercell, create from unit cell
            write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), &
                                trim(adjustl(matter%Name))//trim(adjustl(numpar%path_sep)), 'Unit_cell_equilibrium.txt'
            inquire(file=trim(adjustl(File_name)),exist=file_exist)
            INPUT_SUPCELL2:if (file_exist) then
               open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'old', action='read')
               inquire(file=trim(adjustl(File_name)),opened=file_opened)
               if (.not.file_opened) then
                  Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
                  call Save_error_details(Err, 2, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3416
               endif

               call get_supercell_vectors(FN, File_name, Scell, i, 2, matter, Err) ! see below

            else INPUT_SUPCELL2
               write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name))//' could not be found, the program terminates'
               call Save_error_details(Err, 1, Error_descript)
               call print_message_about_input_files(Error_descript)  ! above
               goto 3416
            endif INPUT_SUPCELL2
         endif SAVED_SUPCELL
         inquire(file=trim(adjustl(File_name)),opened=file_opened)
         if (file_opened) close (FN)

         
         ! Check how to set the atomic coordinates:
         ! a) If user wants path coordinates:
         write(File_name_S1, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'PHASE_1_atoms.dat'
         inquire(file=trim(adjustl(File_name_S1)),exist=file_exist_1)
         write(File_name_S2, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'PHASE_2_atoms.dat'
         inquire(file=trim(adjustl(File_name_S2)),exist=file_exist_2)
         ! Check if user set to calculate along path coordinate:
         numpar%do_path_coordinate = (file_exist_1 .and. file_exist_2)

         ! b) if user set atomic positions in relative units within the supercell:
         FN2 = 9001
         write(File_name2, '(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'SAVE_atoms.dat'
         inquire(file=trim(adjustl(File_name2)),exist=file_exist)
            
         ! Select among different possibilities to set the atomic cell:
         ! In the following priorities:
         ! 1) Path-coordinates  (in internal XTANT SAVE-files format)
         ! 2) SAVE-files  (internal XTANT SAVE-files format)
         ! 3) XYZ file  (extended XYZ format)
         ! 4) POSCAR file
         ! 5) mol2 file
         ! 6) unit-cell coordinates  (old internal XTANT format)

         SAVED_ATOMS:if (numpar%do_path_coordinate) then ! read from the files with initial and final configurations to do the path coordinate plots

            if (numpar%verbose) print*, 'Atomic coordinates from files: ', trim(adjustl(File_name_S1)), trim(adjustl(File_name_S2))

            ! Save the flag for output:
            numpar%save_files_used = 2  ! path coordinate

            inquire(file=trim(adjustl(File_name_S1)),opened=file_opened)
            if (file_opened) close (FN3)
            inquire(file=trim(adjustl(File_name_S2)),opened=file_opened)
            if (file_opened) close (FN4)
            
            ! Get the phase 1 coordinates:
            open(UNIT=FN3, FILE = trim(adjustl(File_name_S1)), status = 'old', action='read')
            call get_initial_atomic_coord(FN3, File_name_S1, Scell, i, 1, matter, numpar, Err, ind = 0) ! below
            ! Get the phase 2 coordinates:
            open(UNIT=FN4, FILE = trim(adjustl(File_name_S2)), status = 'old', action='read')
            call get_initial_atomic_coord(FN4, File_name_S2, Scell, i, 1, matter, numpar, Err, ind = 1) ! below
            
            ! Get atomic temperature set by the velocities given in the SAVE file:
            Natoms = size(Scell(i)%MDatoms)	! number of atoms
            Ta = 0.0d0 ! atomic temperature
            do j = 1,Natoms	! all atoms:
               V2 = SUM(Scell(i)%MDatoms(j)%V(:)*Scell(i)%MDatoms(j)%V(:))*1d10 ! abs value of velocity [A/fs]^2 -> [m/s]^2
               Mass = matter%Atoms(Scell(i)%MDatoms(j)%KOA)%Ma ! atomic mass
               Ta = Ta + Mass*V2/2.0d0/g_e ! Temperature [eV], Eq.(2.62) from H.Jeschke PhD thesis, p.49
            enddo
            Ta = Ta*2.0d0/(3.0d0*real(Natoms) - 6.0d0) ! [eV] proper normalization
            Ta = Ta*g_kb	! [eV] -> [K]

            if (max(Ta,Scell(i)%Ta)/min(Ta+1d-6,Scell(i)%Ta+1d-6) > 1.5d0) then ! if given temperature is too different from the initial one
               ! Set initial velocities according to the given input temperature:
               call set_initial_velocities(matter,Scell,i,Scell(i)%MDatoms,numpar,numpar%allow_rotate) ! below
            endif


         !----------------------------
         elseif (file_exist) then SAVED_ATOMS    ! read from SAVE file with atomic coordinates:

            if (numpar%verbose) print*, 'Atomic coordinates from file: ', trim(adjustl(File_name2))

            ! Save the flag for output:
            numpar%save_files_used = 1  ! Save files read

            open(UNIT=FN2, FILE = trim(adjustl(File_name2)), status = 'old', action='read')
            call get_initial_atomic_coord(FN2, File_name2, Scell, i, 1, matter, numpar, Err) ! below

            ! Get atomic temperature set by the velocities given in the SAVE file:
            Natoms = size(Scell(i)%MDatoms)	! number of atoms
            Ta = 0.0d0 ! atomic temperature
            do j = 1,Natoms	! all atoms:
               V2 = SUM(Scell(i)%MDatoms(j)%V(:)*Scell(i)%MDatoms(j)%V(:))*1d10 ! abs value of velocity [A/fs]^2 -> [m/s]^2
               Mass = matter%Atoms(Scell(i)%MDatoms(j)%KOA)%Ma ! atomic mass
               Ta = Ta + Mass*V2/2.0d0/g_e ! Temperature [eV], Eq.(2.62) from H.Jeschke PhD thesis, p.49
            enddo
            Ta = Ta*2.0d0/(3.0d0*dble(Natoms) - 6.0d0) ! [eV] proper normalization
            Ta = Ta*g_kb	! [eV] -> [K]

            if (max(Ta,Scell(i)%Ta)/min(Ta+1d-6,Scell(i)%Ta+1d-6) > 1.5d0) then ! if given temperature is too different from the initial one
               ! Set initial velocities according to the given input temperature:
               call set_initial_velocities(matter,Scell,i,Scell(i)%MDatoms,numpar,numpar%allow_rotate)  ! below
            endif


         !----------------------------
         elseif (XYZ_file_exists) then SAVED_ATOMS ! XYZ file contains atomic coordinates

            if (numpar%verbose) print*, 'Atomic coordinates from file: ', trim(adjustl(File_name_XYZ))

            open(UNIT=FN_XYZ, FILE = trim(adjustl(File_name_XYZ)), status = 'old', action='read')
            inquire(file=trim(adjustl(File_name_XYZ)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name_XYZ))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3416
            endif

            ! 1) Read the unit cell:
            call read_XYZ(FN_XYZ, File_name_XYZ, Scell, i, matter, numpar, Err) ! see below
            if ( trim(adjustl(Err%Err_descript)) /= '' ) then
               goto 3416
            endif

            ! 2) Make the supercell, if required:
            call get_initial_atomic_coord(FN2, File_name2, Scell, i, 3, matter, numpar, Err) ! below
            if ( trim(adjustl(Err%Err_descript)) /= '' ) then
               goto 3416
            endif

            ! 3) Set initial velocities:
            call set_initial_velocities(matter,Scell,i,Scell(i)%MDatoms,numpar,numpar%allow_rotate) ! below

            inquire(file=trim(adjustl(File_name_XYZ)),opened=file_opened)
            if (file_opened) close (FN_XYZ)


         !----------------------------
         elseif (POSCAR_file_exists) then SAVED_ATOMS ! POSCAR file contains atomic coordinates

            if (numpar%verbose) print*, 'Atomic coordinates from file: ', trim(adjustl(File_name_POSCAR))

            open(UNIT=FN_POSCAR, FILE = trim(adjustl(File_name_POSCAR)), status = 'old', action='read')
            inquire(file=trim(adjustl(File_name_POSCAR)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name_POSCAR))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3416
            endif

            ! 1) Read the unit cell:
            call read_POSCAR(FN_POSCAR, File_name_POSCAR, Scell, i, matter, numpar, Err) ! module "Dealing_with_POSCAR"
            if ( trim(adjustl(Err%Err_descript)) /= '' ) then
               goto 3416
            endif

            ! 2) Cinstruct the supercell from given cells, if required:
            call get_initial_atomic_coord(FN2, File_name2, Scell, i, 3, matter, numpar, Err) ! below
            if ( trim(adjustl(Err%Err_descript)) /= '' ) then
               goto 3416
            endif

            ! 3) Set initial velocities:
            call set_initial_velocities(matter,Scell,i,Scell(i)%MDatoms,numpar,numpar%allow_rotate) ! below

            inquire(file=trim(adjustl(File_name_POSCAR)),opened=file_opened)
            if (file_opened) close (FN_POSCAR)

          !----------------------------
          elseif (mol2_file_exists) then SAVED_ATOMS ! mol2 file contains atomic coordinates

            if (numpar%verbose) print*, 'Atomic coordinates from file: ', trim(adjustl(File_name_mol2))

            open(UNIT=FN_mol2, FILE = trim(adjustl(File_name_mol2)), status = 'old', action='read')
            inquire(file=trim(adjustl(File_name_mol2)),opened=file_opened)
            if (.not.file_opened) then
               Error_descript = 'File '//trim(adjustl(File_name_mol2))//' could not be opened, the program terminates'
               call Save_error_details(Err, 2, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3416
            endif

            ! 1) Read the unit cell:
            call read_mol2(FN_mol2, File_name_mol2, Scell, i, matter, numpar, Err) ! module "Dealing_with_mol2"
            if ( trim(adjustl(Err%Err_descript)) /= '' ) then
               goto 3416
            endif

            ! 2) Cinstruct the supercell from given cells, if required:
            call get_initial_atomic_coord(FN2, File_name2, Scell, i, 3, matter, numpar, Err) ! below
            if ( trim(adjustl(Err%Err_descript)) /= '' ) then
               goto 3416
            endif

            ! 3) Set initial velocities:
            call set_initial_velocities(matter, Scell, i, Scell(i)%MDatoms, numpar, numpar%allow_rotate) ! below

            inquire(file=trim(adjustl(File_name_mol2)), opened=file_opened)
            if (file_opened) close (FN_mol2)

         !----------------------------
         else SAVED_ATOMS  ! unit-cell file
            ! Save the flag for output:
            numpar%save_files_used = 0  ! Save files read

            ! c) if user set to construct supercell from unit cells:
            write(File_name2,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
                                       'Unit_cell_atom_relative_coordinates.txt'

            if (numpar%verbose) print*, 'Atomic coordinates from file: ', trim(adjustl(File_name2))

            inquire(file=trim(adjustl(File_name2)),exist=file_exist)
            INPUT_ATOMS:if (file_exist) then
               open(UNIT=FN2, FILE = trim(adjustl(File_name2)), status = 'old', action='read')
               inquire(file=trim(adjustl(File_name2)),opened=file_opened)
               if (.not.file_opened) then
                  Error_descript = 'File '//trim(adjustl(File_name2))//' could not be opened, the program terminates'
                  call Save_error_details(Err, 2, Error_descript)
                  print*, trim(adjustl(Error_descript))
                  goto 3416
               endif

               call get_initial_atomic_coord(FN2, File_name2, Scell, i, 2, matter, numpar, Err) ! below
               ! Set initial velocities:
               call set_initial_velocities(matter,Scell,i,Scell(i)%MDatoms,numpar,numpar%allow_rotate) ! below

            else INPUT_ATOMS
               write(Error_descript,'(a,$)') 'File '//trim(adjustl(File_name2))//' could not be found, the program terminates'
               call Save_error_details(Err, 1, Error_descript)
               call print_message_about_input_files(Error_descript)  ! above
               goto 3416
            endif INPUT_ATOMS
         endif SAVED_ATOMS
         inquire(file=trim(adjustl(File_name2)),opened=file_opened)
         if (file_opened) close (FN2)
         

         !---------------------------
         ! Check embedding in water:
         if (numpar%embed_water) then
            call embed_molecule_in_water(Scell, matter, numpar)  ! below
         endif


         !---------------------------
         ! Check periodicity:
         call Make_free_surfaces(Scell, numpar, matter)	! module "Atomic_tools"
         !call Coordinates_rel_to_abs(Scell, i, if_old=.true.)	! from the module "Atomic_tools"
         ! Save atomic coordinates at their equilibrium positions:
         do j = 1, Scell(i)%Na
            Scell(i)%MDatoms(j)%R_eq(:) = Scell(i)%MDatoms(j)%R(:)   ! save coords at equilibrium positions to ger mean square displacements
            Scell(i)%MDatoms(j)%S_eq(:) = Scell(i)%MDatoms(j)%S(:)   ! save coords at equilibrium positions to ger mean square displacements
            !print*, 'After Make_free_surfaces', Scell(i)%MDAtoms(j)%R(:), Scell(i)%MDAtoms(j)%S(:), Scell(i)%supce
         enddo ! j
         !---------------------------
         

         ! Allocate nearest neighbor lists:
         if (.not.allocated(Scell(i)%Near_neighbor_list)) allocate(Scell(i)%Near_neighbor_list(Scell(i)%Na,Scell(i)%Na))  ! nearest neighbors
         if (.not.allocated(Scell(i)%Near_neighbor_dist)) allocate(Scell(i)%Near_neighbor_dist(Scell(i)%Na,Scell(i)%Na,4))  ! [A] distances
         if (.not.allocated(Scell(i)%Near_neighbor_dist_s)) allocate(Scell(i)%Near_neighbor_dist_s(Scell(i)%Na,Scell(i)%Na,3)) ! relative dist.
         if (.not.allocated(Scell(i)%Near_neighbor_size)) allocate(Scell(i)%Near_neighbor_size(Scell(i)%Na)) ! how many nearest neighbours
         if (numpar%save_NN) then   ! if user wants to study number of nearest neighbors within defined radius
            if (.not.allocated(Scell(i)%Near_neighbors_user)) allocate(Scell(i)%Near_neighbors_user(Scell(i)%Na)) ! user-defined nearest neighbours
         endif

         ! Allocate Hamiltonan matrices:
         !if (.not.allocated(Scell(i)%Ha)) allocate(Scell(i)%Ha(Scell(i)%Ne,Scell(i)%Ne))   ! hamiltonian size
         ASSOCIATE (ARRAY => Scell(i)%TB_Hamil(:,:))
            select type(ARRAY)
            type is (TB_H_Pettifor)
               n1 = size(ARRAY(i,i)%V0)*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
            type is (TB_H_Molteni)
               n1 = size(ARRAY(i,i)%V0)*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
            type is (TB_H_Fu)
               n1 = size(ARRAY(i,i)%V0)*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
            type is (TB_H_NRL)	! it is always sp3d5 basis set, so 9 orbitals per atom:
               n1 = 9.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
!                if (.not.allocated(Scell(i)%Vij)) allocate(Scell(i)%Vij(n1,n1))	! radial part of hamiltonian size
               if (.not.allocated(Scell(i)%Sij)) allocate(Scell(i)%Sij(n1,n1))	! Overlap matrix for non-orthogonal TB
               if (.not.allocated(Scell(i)%Hij)) allocate(Scell(i)%Hij(n1,n1))	! Non-orthogonal TB Hamiltonian
               if (.not.allocated(Scell(i)%Hij_sol)) allocate(Scell(i)%Hij_sol(n1,n1))	! eigenvectors of nondiagonalized Hamiltonian
            type is (TB_H_DFTB)   ! it can be various basis sets:
               select case (numpar%N_basis_size)    ! find which one is used now:
               case (0)    ! s
                  n1 = 1.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               case (1)    ! sp3
                  n1 = 4.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               case default    ! sp3d5
                  n1 = 9.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               endselect
               if (.not.allocated(Scell(i)%Sij)) allocate(Scell(i)%Sij(n1,n1))	! Overlap matrix for non-orthogonal TB
               if (.not.allocated(Scell(i)%Hij)) allocate(Scell(i)%Hij(n1,n1))	! Non-orthogonal TB Hamiltonian
               if (.not.allocated(Scell(i)%Hij_sol)) allocate(Scell(i)%Hij_sol(n1,n1))	! eigenvectors of nondiagonalized Hamiltonian
            type is (TB_H_3TB)   ! it can be various basis sets:
               select case (numpar%N_basis_size)    ! find which one is used now:
               case (0)    ! s
                  n1 = 1.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               case (1)    ! sp3
                  n1 = 4.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               case default    ! sp3d5
                  n1 = 9.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               endselect
               if (.not.allocated(Scell(i)%Sij)) allocate(Scell(i)%Sij(n1,n1))	! Overlap matrix for non-orthogonal TB
               if (.not.allocated(Scell(i)%Hij)) allocate(Scell(i)%Hij(n1,n1))	! Non-orthogonal TB Hamiltonian
               if (.not.allocated(Scell(i)%Hij_sol)) allocate(Scell(i)%Hij_sol(n1,n1))	! eigenvectors of nondiagonalized Hamiltonian
            type is (TB_H_BOP)   ! it can be various basis sets:
               select case (numpar%N_basis_size)    ! find which one is used now:
               case (0)    ! s
                  n1 = 1.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               case (1)    ! sp3
                  n1 = 4.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               case default    ! sp3d5
                  n1 = 9.0d0*Scell(i)%Na ! number of energy levels is defined by the number of TB parameters included
               endselect
               if (.not.allocated(Scell(i)%Sij)) allocate(Scell(i)%Sij(n1,n1))	! Overlap matrix for non-orthogonal TB
               if (.not.allocated(Scell(i)%Hij)) allocate(Scell(i)%Hij(n1,n1))	! Non-orthogonal TB Hamiltonian
               if (.not.allocated(Scell(i)%Hij_sol)) allocate(Scell(i)%Hij_sol(n1,n1))	! eigenvectors of nondiagonalized Hamiltonian
            type is (TB_H_xTB)
               ! number of energy levels is defined by the number of TB parameters included:
               n1 = identify_xTB_orbitals_per_atom(numpar%N_basis_size) * Scell(i)%Na ! module "TB_xTB"
               if (.not.allocated(Scell(i)%Sij)) allocate(Scell(i)%Sij(n1,n1))	! Overlap matrix for non-orthogonal TB
               if (.not.allocated(Scell(i)%Hij)) allocate(Scell(i)%Hij(n1,n1))	! Non-orthogonal TB Hamiltonian
               if (.not.allocated(Scell(i)%Hij_sol)) allocate(Scell(i)%Hij_sol(n1,n1))	! eigenvectors of nondiagonalized Hamiltonian
            end select
         END ASSOCIATE

         if (.not.allocated(Scell(i)%Ha)) allocate(Scell(i)%Ha(n1,n1))   ! hamiltonian size
         if (.not.allocated(Scell(i)%Ha0)) allocate(Scell(i)%Ha0(n1,n1)) ! hamiltonian0 size
         if (.not.allocated(Scell(i)%H_non)) allocate(Scell(i)%H_non(n1,n1))	! nondiagonalized Hamiltonian
         if (.not.allocated(Scell(i)%H_non0)) allocate(Scell(i)%H_non0(n1,n1))	! nondiagonalized Hamiltonian
         if (.not.allocated(Scell(i)%Ei)) allocate(Scell(i)%Ei(n1))  ! energy levels, eigenvalues of the hamiltonian matrix
         if (.not.allocated(Scell(i)%Ei0)) allocate(Scell(i)%Ei0(n1))  ! energy levels0, eigenvalues of the hamiltonian matrix
         if (.not.allocated(Scell(i)%Mij)) allocate(Scell(i)%Mij(n1,n1), source = 0.0d0)   ! dynamical coupling matrix
         if (.not.allocated(Scell(i)%Aij)) allocate(Scell(i)%Aij(n1,n1))	! coefficients used for forces in TB
         if ((numpar%scc) .and. (.not.allocated(Scell(i)%Ei_scc_part)) ) then
            allocate(Scell(i)%Ei_scc_part(n1))  ! energy levels of non-SCC part of the hamiltonian
         endif
         if (allocated(Scell(i)%Sij) .and. .not.allocated(Scell(i)%eigen_S)) allocate(Scell(i)%eigen_S(n1)) ! eigenvalues of Sij
         
         ! Electron distribution function:
         if (.not. allocated(Scell(i)%fe)) allocate(Scell(i)%fe(size(Scell(i)%Ei))) ! electron distribution function
         if (.not. allocated(Scell(i)%fe_eq)) allocate(Scell(i)%fe_eq(size(Scell(i)%Ei))) ! equivalent distribution function (Fermi-function)
         ! Check if there is a file with the initial distribution:
         write(File_name2,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, 'SAVE_el_distribution.dat' ! default filename
         inquire(file=trim(adjustl(File_name2)),exist=file_exist_1)  ! check if default file exists

         if ((Scell(i)%Te < 0.0d0) .or. file_exist_1) then ! distribution must be provided in the file
            if (.not.file_exist_1) then ! no default file, must be provided by the user then
               write(File_name2,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(numpar%fe_filename)) ! user-provided filename
               inquire(file=trim(adjustl(File_name2)),exist=file_exist)
            else ! default file
               numpar%fe_filename = 'SAVE_el_distribution.dat'
               file_exist = file_exist_1
               Scell(i)%Te = -1.0e0
            endif

            INPUT_DISTR:if (file_exist) then
               numpar%fe_input_exists = .true.  ! distribution was provided
               open(UNIT=FN2, FILE = trim(adjustl(File_name2)), status = 'old', action='read')
               inquire(file=trim(adjustl(File_name2)),opened=file_opened)
               if (.not.file_opened) then
                  numpar%fe_input_exists = .false.  ! no distribution given, use Fermi from the start
               endif
               ! Read distribution from the file:
               call read_electron_distribution(FN2, numpar%fe_input, numpar%fe_input_exists)   ! below
            else INPUT_DISTR
               numpar%fe_input_exists = .false.  ! no distribution given, use Fermi from the start
            endif INPUT_DISTR

            ! If for any reason the distribution could not be read from the file, use Fermi distribution:
            if (.not.numpar%fe_input_exists) then
               ! Assume electronic temperature equal to the atomic one:
               Scell(i)%Te = Scell(i)%Ta  ! [K]
               Scell(i)%TeeV = Scell(i)%Te/g_kb ! [eV] electron temperature
               print*, 'File '//trim(adjustl(File_name2))//' could not be opened, use electron temperature: ', Scell(i)%Te
            endif
         else
            ! No distribution file provided by the user, use Fermi with given temperature instead
         endif

         ! If partial band-specific distributions are required:
         if ((numpar%tau_fe_CB > -1.0d-8) .and. (numpar%tau_fe_VB > -1.0d-8)) then ! Partial thermalization is on
            ! Define equivalent distribution functions (Fermi-functions) in VB and CB:
            if (.not. allocated(Scell(i)%fe_eq_VB)) allocate(Scell(i)%fe_eq_VB(size(Scell(i)%Ei)), source = 0.0d0)
            if (.not. allocated(Scell(i)%fe_eq_CB)) allocate(Scell(i)%fe_eq_CB(size(Scell(i)%Ei)), source = 0.0d0)
         endif

         if (numpar%do_kappa) then
            if (.not. allocated(Scell(i)%I_ij)) allocate(Scell(i)%I_ij(size(Scell(i)%Ei))) ! scattering integral
            if (.not. allocated(Scell(i)%Ce_i)) allocate(Scell(i)%Ce_i(size(Scell(i)%Ei))) ! electron heat capacity
         endif
!          if (.not. allocated(Scell(i)%Norm_WF)) allocate(Scell(i)%Norm_WF(size(Scell(i)%Ei))) ! normalization coefficient of the wave function

         ! DOS masks:
         call get_DOS_masks(Scell, matter, numpar)  ! module "TB"
         
         do j = 1,size(laser) ! for each pulse:
            laser(j)%Fabs = laser(j)%F*dble(Scell(i)%Na) ! total absorbed energy by supercell [eV]
            laser(j)%Nph = laser(j)%Fabs/laser(i)%hw     ! number of photons absorbed in supercell
         enddo

         temp = 0.0d0
         do j = 1, Scell(i)%Na
            temp = temp + matter%Atoms(Scell(i)%MDatoms(j)%KOA)%Ma ! total mass of all atoms in supercell
         enddo
         matter%W_PR = temp/matter%W_PR ! Mass of unit cell in the Parrinello-Rahman method [kg]

      enddo ALL_SC
   endif MD

   ! Initialize MC data:
   do Nsc = 1, size(Scell)
      if (.not. allocated(MC)) then
         Scell(Nsc)%Nph = 0.0d0		! number of absorbed photons
         Scell(Nsc)%Ne_high = 0.0d0	! no electrons at the beginning
         Scell(Nsc)%Ne_emit = 0.0d0	! no emitted electrons at the beginning
         Scell(Nsc)%Nh = 0.0d0		! no holes at the beginning
         allocate(Scell(Nsc)%MChole(size(matter%Atoms)))
         do i = 1, size(matter%Atoms) ! for each kind of atoms:
            allocate(Scell(Nsc)%MChole(i)%Noh(matter%Atoms(i)%sh))
            Scell(Nsc)%MChole(i)%Noh(:) = 0.0d0 ! no holes in any shell
         enddo
         if (numpar%NMC > 0) then
            allocate(MC(numpar%NMC))	! all MC arrays for photons, electrons and holes
         !if (size(MC) > 0) then
            do i = 1, size(MC)
               MC(i)%noe = 0.0d0
               MC(i)%noe_emit = 0.0d0
               MC(i)%noh_tot = 0.0d0
               allocate(MC(i)%electrons(Scell(Nsc)%Ne))
               MC(i)%electrons(:)%E = 0.0d0
               MC(i)%electrons(:)%ti = 1d25
               MC(i)%electrons(:)%colls = 0
               allocate(MC(i)%holes(Scell(Nsc)%Ne))
               MC(i)%holes(:)%E = 0.0d0
               MC(i)%holes(:)%ti = 1d26
            enddo
         endif !(size(MC) > 0)
      endif
   enddo

   ! If we didn't set the density for MC, use it from the MD part:
   if (matter%dens <= 0.0d0) then
      matter%At_dens = dble(Scell(1)%Na)/(Scell(1)%V*1.0d-24)
      matter%dens = matter%At_dens*(SUM(matter%Atoms(:)%Ma*matter%Atoms(:)%percentage)/(SUM(matter%Atoms(:)%percentage))*1d3) ! just in case there was no better given density (no cdf file was used)
   else
      matter%At_dens = matter%dens/(SUM(matter%Atoms(:)%Ma*matter%Atoms(:)%percentage)/(SUM(matter%Atoms(:)%percentage))*1d3)   ! atomic density [1/cm^3]
   endif

3416 continue
!     do i = 1, Scell(1)%Na
!        write(6,'(i4,f,f,f,f,f,f)') i, Scell(1)%MDAtoms(i)%S0(:), Scell(1)%MDAtoms(i)%S(:)
!     enddo ! j
!    pause 'set_initial_configuration'
end subroutine set_initial_configuration



subroutine read_XYZ(FN_XYZ, File_name_XYZ, Scell, SCN, matter, numpar, Err) ! extended XYZ format
   integer, intent(in) :: FN_XYZ ! extended XYZ file number (must be already open)
   character(*), intent(in) :: File_name_XYZ ! extended XYZ file name
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: SCN ! number of the supercell (always =1)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------

   integer :: count_lines, Reason
   logical :: read_well
   character(1000) :: line_2, Error_descript

   count_lines = 0   ! to start with

   ! First line in XYZ, number of atoms:
   read(FN_XYZ,*,IOSTAT=Reason) Scell(SCN)%Na
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_XYZ))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3418
   endif

   ! Second line, in extended XYZ format contains important information:
   line_2 = ''
   read(FN_XYZ,'(a)',IOSTAT=Reason) line_2  ! read the full line, then interpret it
   call read_file(Reason, count_lines, read_well)
   if (.not. read_well) then
      write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_XYZ))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3418
   endif

   ! Having read this line, act accordingly:
   call interpret_XYZ_comment(FN_XYZ, File_name_XYZ, count_lines, line_2, Scell, SCN, matter, numpar, Err)   ! below

   ! Get the volume of the now-defined supercell:
   call Det_3x3(Scell(SCN)%supce,Scell(SCN)%V) ! module "Algebra_tools"

   3418 continue
end subroutine read_XYZ


subroutine interpret_XYZ_comment(FN_XYZ, File_name_XYZ, count_lines, line_2, Scell, SCN, matter, numpar, Err)
   integer, intent(in) :: FN_XYZ, SCN
   integer, intent(inout) :: count_lines
   character(*), intent(in) :: File_name_XYZ ! extended XYZ file name
   character(*), intent(in) :: line_2  ! line #2 from extended XYZ file
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !--------------------
   integer :: ind_S, ind_R, ind_V, ind_atoms
   real(8) :: SC_X, SC_Y
   character(200) :: Error_descript

   call interpret_XYZ_comment_line(line_2, Scell(SCN)%Supce, ind_S, ind_R, ind_V, ind_atoms, SC_X, SC_Y, Error_descript)     ! module "Dealing_with_eXYZ"
   if (trim(adjustl(Error_descript)) /= '') then
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3420
   endif

   if (numpar%verbose) print*, 'In interpret_XYZ_comment, the ind_atoms=', ind_atoms

   ! If atomic coordinates are provided:
   select case (ind_atoms)
   case (1) ! there are data for atomic species and coordinates
      Scell(SCN)%Supce0 = Scell(SCN)%Supce   ! initial
      if (numpar%verbose) print*, 'Reading defined atomic coordinates from xyz-file'
      call read_XYZ_coords(FN_XYZ, File_name_XYZ, count_lines, Scell, SCN, matter, ind_S, ind_R, ind_V, Err) ! below
   case (0) ! to be set randomly
      if (numpar%verbose) print*, 'Setting random atomic coordinates defined in xyz-file'
      call read_XYZ_random(FN_XYZ, File_name_XYZ, count_lines, Scell, SCN, matter, SC_X, SC_Y, numpar, Err) ! below
   case default
      Error_descript = 'Could not interpret the data in file '//trim(adjustl(File_name_XYZ))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
   endselect

3420 continue
end subroutine interpret_XYZ_comment


subroutine read_XYZ_random(FN_XYZ, File_name_XYZ, count_lines, Scell, SCN, matter, SC_X_in, SC_Y_in, numpar, Err)
   integer, intent(in) :: FN_XYZ, SCN
   integer, intent(inout) :: count_lines
   character(*), intent(in) :: File_name_XYZ ! extended XYZ file name
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   real(8), intent(in) :: SC_X_in, SC_Y_in   ! [A] supercell size along X and Y
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !------------------------
   real(8) :: SC_x, SC_y, eps, V_at, RN(3), lower_Z, upper_Z, a_r, dR_min, iter_max
   real(8), dimension(:), allocatable :: SC_z, rho
   integer, dimension(:), allocatable :: NOA, KOA
   character(3), dimension(:), allocatable :: El_name
   real(8), dimension(:), allocatable :: at_r_cov, at_mass
   real(8), dimension(:), allocatable :: at_r_cov_read, at_mass_read
   integer :: i, j, N_at, INFO, N_temp, Reason, i_counter, at_counter, iter
   character(300) :: Error_descript, Folder_name
   logical :: read_well, redo_placement

   ! Folder with the periodic table:
   Folder_name = trim(adjustl(numpar%input_path))//trim(adjustl(m_Atomic_parameters))

   ! Get the number of various elements in the target:
   call Count_lines_in_file(FN_XYZ, N_at)  ! module "Dealing_with_files"
   ! Get back to the same line in the file:
   rewind(FN_XYZ)
   read(FN_XYZ,*)
   read(FN_XYZ,*)

   ! Knowing number of different elements, allocate arrays:
   allocate(SC_z(N_at))
   allocate(rho(N_at))
   allocate(NOA(N_at))
   allocate(KOA(N_at))
   allocate(El_name(N_at))
   allocate(at_r_cov(N_at))
   allocate(at_mass(N_at))

   ! Read the data for each element:
   ELS:do i = 1, N_at
      read(FN_XYZ,*,IOSTAT=Reason) El_name(i), rho(i), NOA(i)
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_XYZ))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3421
      endif

      call Decompose_compound(Folder_name, El_name(i), numpar%path_sep, INFO, Error_descript, N_temp, &
                              at_masses=at_mass_read, at_r_cov=at_r_cov_read) ! molude 'Periodic_table'
      if (INFO /= 0) then ! problem in Periodic table reading
         Error_descript = trim(adjustl(Error_descript))//' Called from the file '//trim(adjustl(File_name_XYZ))
         call Save_error_details(Err, 4, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3421
      endif

      ! Save the data for reuse below:
      at_r_cov(i) = at_r_cov_read(1)
      at_mass(i) = at_mass_read(1)

      ! Find the index from the element name:
      call get_KOA_from_element(El_name(i), matter, KOA(i)) ! module "Dealing_with_POSCAR"
      if (KOA(i) <= 0) then
         write(Error_descript,'(a,a,a)') 'Inconsistency: in the target, there is no element ', &
                                          trim(adjustl(El_name(i))), ' from file '//trim(adjustl(File_name_XYZ))
         call Save_error_details(Err, 3, Error_descript)
         print*, '---------------------------'
         print*, trim(adjustl(Error_descript))
         goto 3421
      endif
      !Scell(SCN)%MDAtoms(i)%KOA = KOA
   enddo ELS

   dR_min = minval(at_r_cov(:))*0.99d0  ! minimum allowed distance between atoms

   ! Knowing the data we want to set, we can set it:
   Scell(SCN)%Na = SUM(NOA)
   allocate(Scell(SCN)%MDAtoms(Scell(SCN)%Na))
   ! Supercell:
   eps = 1.0d-10
   if ((SC_X_in > eps) .and. (SC_Y_in > eps)) then
      SC_x = SC_X_in
      SC_y = SC_Y_in
   endif
   ! The supercell size is defined from the density and number of atoms:
   do i = 1, N_at ! for each layer
      V_at = NOA(i) * g_amu * at_mass(i) / rho(i) * 1.0d27   ! Volume of the part of the supercell [A^3]
      if ((SC_X_in > eps) .and. (SC_Y_in > eps)) then ! X and Y are fixed, find Z-dimension of the supercell:
         SC_z(i) = V_at / (SC_x*SC_y)  ! [A]
      else  ! cubic supercell:
         SC_z(i) = (V_at)**m_one_third ! [A]
         SC_x = SC_z(1)
         SC_y = SC_z(1)
      endif
   enddo
   ! Now we have the total supercell size:
   Scell(SCN)%supce(:,:) = 0.0d0 ! off-diagonals
   Scell(SCN)%supce(1,1) = SC_x
   Scell(SCN)%supce(2,2) = SC_y
   Scell(SCN)%supce(3,3) = SUM(SC_z)
   Scell(SCN)%supce0 = Scell(SCN)%supce

   if (numpar%verbose) print*, 'Supercell sizes defined:', Scell(SCN)%supce(1,1), Scell(SCN)%supce(2,2), Scell(SCN)%supce(3,3)


   ! Now, place atoms randomly, according to the conditions specified:
   iter_max = 10000
   at_counter = 1
   i_counter = NOA(at_counter)
   lower_Z = 0.0d0   ! to start with
   upper_Z = SC_z(1) ! to start with
   do i = 1, Scell(SCN)%Na ! define all atoms
      ! Which type of atom is this one:
      if (numpar%verbose) print*, 'Random setting of atom #', i, Scell(SCN)%Na, i_counter
      if (i > i_counter) then
         at_counter = at_counter + 1   ! next type of atoms
         i_counter = i_counter + NOA(at_counter)   ! next number of atoms to start the new set
         lower_Z = upper_Z   ! next layer
         upper_Z = upper_Z + SC_z(at_counter) ! next layer
      endif
      ! Specify atomic parameters:
      Scell(SCN)%MDAtoms(i)%KOA = KOA(at_counter)
      !print*, 'Atom', i, Scell(SCN)%MDAtoms(i)%KOA, at_counter

      redo_placement = .true.  ! to start with
      iter = 0 ! count iteration of attempted placement
      do while (redo_placement)
         redo_placement = .false. ! assume we place it well, no need to redo it
         ! Set random coordinates of the new water molecule:
         call random_number(RN)  ! random numbers for relative coordinates along X,Y,Z
         Scell(SCN)%MDAtoms(i)%R(1) = Scell(SCN)%supce(1,1)*RN(1)
         Scell(SCN)%MDAtoms(i)%R(2) = Scell(SCN)%supce(2,2)*RN(2)
         Scell(SCN)%MDAtoms(i)%R(3) = lower_Z + (upper_Z - lower_Z)*RN(3)
         Scell(SCN)%MDAtoms(i)%R0(:) = Scell(SCN)%MDAtoms(i)%R(:)
         ! Get the relative coordinates from the absolute ones provided:
         !call Coordinates_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
         call Coordinates_abs_to_rel_single(Scell, SCN, i, if_old = .true.) ! module "Atomic_tools"
         ! Check that it is not overlapping with existing atoms:
         CHKI:do j = 1, i-1
            ! Get the relative distance to this atom:
            call shortest_distance(Scell(SCN), i, j, a_r) ! module "Atomic_tools"
            ! Check if it is not too short:
            if (a_r < dR_min) then ! overlapping atoms, place a molecule in a new place:
               redo_placement = .true. ! assume we place it well, no need to redo it
               iter = iter + 1   ! next iteration
               exit CHKI
            endif
         enddo CHKI
         ! If there is no place to place the atoms, maybe the supercell must be a little larger:
         if (iter > iter_max) then
            print*, 'Increasing supercell size:', Scell(SCN)%supce(1,1), '->', Scell(SCN)%supce(1,1)*1.01d0, Scell(SCN)%supce(3,3)
            iter = 0 ! restart
            Scell(SCN)%supce(1,1) = Scell(SCN)%supce(1,1)*1.01d0
            Scell(SCN)%supce(2,2) = Scell(SCN)%supce(2,2)*1.01d0
            Scell(SCN)%supce0 = Scell(SCN)%supce
            ! Absolute coordinate change since Supercell changed:
            call Coordinates_rel_to_abs(Scell, SCN, if_old=.true.)   ! from the module "Atomic_tools"
         endif
      enddo ! do while (redo_placement)
   enddo ! i = 1, Scell(SCN)%Na

   ! Clean up:
   deallocate(SC_z, rho, NOA, KOA, El_name, at_r_cov, at_mass)
3421 continue
end subroutine read_XYZ_random




subroutine read_XYZ_coords(FN, File_name_XYZ, count_lines, Scell, SCN, matter, ind_S, ind_R, ind_V, Err)
   integer, intent(in) :: FN, SCN, ind_S, ind_R, ind_V
   integer, intent(inout) :: count_lines
   character(*), intent(in) :: File_name_XYZ ! extended XYZ file name
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !------------------------
   integer :: i, Reason, KOA
   real(8), dimension(3) :: coord, vel
   logical :: read_well
   character(3) :: El_name
   character(200) :: Error_descript

   allocate(Scell(SCN)%MDAtoms(Scell(SCN)%Na))

   if (ind_S < 0) then
      write(Error_descript,'(a,i3,a,$)') 'Could not interprete line #', count_lines, ' in file '//trim(adjustl(File_name_XYZ))
      call Save_error_details(Err, 3, Error_descript)
      print*, trim(adjustl(Error_descript))
      goto 3419
   endif

   ! Next lines contain at least the atomic type and coordinates:
   do i = 1, Scell(SCN)%Na
      if (ind_S == 0) then ! KOA is set in the file
         if (ind_V == 1) then ! with velocities:
            read(FN,*,IOSTAT=Reason) Scell(SCN)%MDAtoms(i)%KOA, coord(:), vel(:)
         else  ! no velocities, only coordinates:
            read(FN,*,IOSTAT=Reason) Scell(SCN)%MDAtoms(i)%KOA, coord(:)
         endif
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_XYZ))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3419
         endif

      else  ! element name is set
         if (ind_V == 1) then ! with velocities:
            read(FN,*,IOSTAT=Reason) El_name, coord(:), vel(:)
         else  ! no velocities, only coordinates:
            read(FN,*,IOSTAT=Reason) El_name, coord(:)
         endif
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name_XYZ))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3419
         endif

         ! Find the index from the element name:
         call get_KOA_from_element(El_name, matter, KOA) ! module "Dealing_with_POSCAR"
         if (KOA <= 0) then
            write(Error_descript,'(a,i3,a,$)') 'In the target, there is no element ', trim(adjustl(El_name)), ' from file '//trim(adjustl(File_name_XYZ))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3419
         endif
         Scell(SCN)%MDAtoms(i)%KOA = KOA
      endif ! (ind_S == 0)

      ! Sort the atomic coordinates that were read:
      if (ind_R == 1) then
         Scell(SCN)%MDAtoms(i)%R(:) = coord(:)
         Scell(SCN)%MDAtoms(i)%R0(:) = Scell(SCN)%MDAtoms(i)%R(:)
         ! Get the relative coordinates from the absolute ones provided:
         call Coordinates_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
         !print*, 'read_XYZ_coords', Scell(SCN)%MDAtoms(i)%R(:), Scell(SCN)%MDAtoms(i)%S(:), Scell(SCN)%supce
      else
         Scell(SCN)%MDAtoms(i)%S(:) = coord(:)
         Scell(SCN)%MDAtoms(i)%S0(:) = Scell(SCN)%MDAtoms(i)%S(:)
      endif

      ! Sort the atomic velocities that were read:
      if (ind_V == 1) then
         Scell(SCN)%MDAtoms(i)%V(:) = vel(:)
         Scell(SCN)%MDAtoms(i)%V0(:) = Scell(SCN)%MDAtoms(i)%V(:)
         ! Get the relative velocities from the absolute ones provided:
         call velocities_abs_to_rel(Scell, SCN, if_old=.true.) ! module "Atomic_tools"
      elseif (ind_V == 2) then
         Scell(SCN)%MDAtoms(i)%SV(:) = vel(:)
         Scell(SCN)%MDAtoms(i)%SV0(:) = Scell(SCN)%MDAtoms(i)%SV(:)
      endif
   enddo ! i = 1, Scell(SCN)%Na
   3419 continue
end subroutine read_XYZ_coords



subroutine read_electron_distribution(FN2, fe_input, fe_input_exists)
   integer, intent(in) :: FN2 ! file with distribution
   real(8), dimension(:), allocatable, intent(inout) :: fe_input  ! initial distribution function
   logical, intent(inout) :: fe_input_exists ! flag to use the distribution from a file
   !--------------------
   integer :: N_cols, N_lines, i, Reason, count_lines
   real(8), dimension(:), allocatable :: temp_fe
   logical :: read_well

   ! Count how many grid points are there:
   call Count_columns_in_file(FN2, N_cols, skip_lines=1) ! module "Dealing_with_files"
   ! Count how many columns are in the file (assume the last one is the distriution):
   call Count_lines_in_file(FN2, N_lines, skip_lines=1)   ! module "Dealing_with_files"

   ! Knowing the size, allocate the arrays:
   allocate(fe_input(N_lines), source=0.0d0)
   allocate(temp_fe(N_cols))

   ! Read the data from the file:
   read(FN2,*,IOSTAT=Reason) ! skip the first comment line
   count_lines = 1
   ! Read the rest as distribution:
   READ_DISTR:do i = 1, N_lines
      read(FN2,*,IOSTAT=Reason) temp_fe
      call read_file(Reason, count_lines, read_well) ! module "Dealing_with_files"
      if (.not. read_well) then
         fe_input_exists = .false.  ! could not read distribution
         exit READ_DISTR
      else  ! save distribution function:
         fe_input(i) = temp_fe(N_cols)
         !print*, 'read_electron_distribution:', i, fe_input(i)
      endif
   enddo READ_DISTR
!    pause 'read_electron_distribution DONE'

   ! Clean up:
   deallocate(temp_fe)
   call close_file('close', FN=FN2) ! module "Dealing_with_files"
end subroutine read_electron_distribution



subroutine embed_molecule_in_water(Scell, matter, numpar)  ! below
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   !---------------------------
   integer :: i, j, i_mol, i_h2o, N_at, N_h2o, SCN, N_tot, j_h2o, iter, i_H, i_O, iter_max
   real(8) :: V_mol, V_tot, dV, s_center(3), RN(3), dS_min(3), dS(3), dS_min_abs, dV_resc, dR_min
   real(8) :: theta, phi, cos_phi, u0, v0, w0, u, v, w, a_r
   logical :: redo_placement
   type(Atom), dimension(:), allocatable :: MDAtoms ! all atoms in MD

   !-----------------------
   ! 0) Define initial parameters:
   iter_max = 10000
   dR_min = 3.0d0*g_a0  ! [A] atoms no closer than this distance
   dV_resc = (1.1d0)**(1.0d0/3.0d0)  ! rescaling factor of increase of volume in case if needed (molecule overlap)
   SCN = 1  ! so far, only one supercell
   iter = 1 ! to start with
   N_at = Scell(SCN)%Na ! number of atoms in bio/molecule
   N_h2o = 3*numpar%N_water_mol  ! number of atoms from water molecules (2*H+1*O)
   N_tot = N_at + N_h2o ! total number of atoms after adding water
   allocate(MDAtoms(N_at))
   MDAtoms = Scell(SCN)%MDAtoms  ! save initial for reuse
   ! Increase the size of MDAtoms array:
   deallocate(Scell(SCN)%MDAtoms)
   allocate(Scell(SCN)%MDAtoms(N_tot))
   ! Now find the indices of H and O atomic types:
   i = 1
   i_H = i
   do while ( trim(adjustl(matter%Atoms(i)%Name)) /= 'H')
      i = i + 1
      i_H = i
   enddo
   i = 1
   i_O = i
   do while ( trim(adjustl(matter%Atoms(i)%Name)) /= 'O')
      i = i + 1
      i_O = i
   enddo

!    print*, i_H, i_O, matter%Atoms(i_H)%Name, matter%Atoms(i_O)%Name
!    pause

   !-----------------------
   ! 1) Estimate size of the new box that will include water:
   ! 1.a) estimated volume of the molecule (assuming the radius of each atom ~ 0.53 A):
   V_mol = N_at * 4.0d0/3.0d0*g_Pi*g_a0**3   ! [A^3]
   ! 1.b) estimated volume of water molecules (assuming water density of 1 g/cc):
   !V_tot = V_mol + N_h2o*(1d3 / g_amu / 1d30 / 3.0d0) ! total volume [A^3]
   ! 1.c) Volume coefficient to expand the supercell:
   call Det_3x3(Scell(SCN)%supce,Scell(SCN)%V) ! module "Algebra_tools"
   ! estimated volume of water molecules (assuming water density of 1 g/cc):
   V_tot = Scell(SCN)%V + N_h2o*( ((2.0d0*1.0d0 + 16.0d0)*g_amu)/3.0d0 / 1d3 * 1d30) ! total volume [A^3]
   dV = (V_tot/Scell(SCN)%V)**(1.0d0/3.0d0)
   ! 1.d) Rescale the supercell:
2023 continue

   if (iter >= iter_max) then ! were unable to fit all molecules into the box, increasing the box side (lowering density)
      print*, 'All water molecules did not fit in the box, incrasing box size:'
      print*, 'dV=', dV, Scell(SCN)%supce(3,3)
   endif

   Scell(SCN)%supce(:,:) = Scell(SCN)%supce(:,:) * dV
   Scell(SCN)%supce0(:,:) = Scell(SCN)%supce(:,:)

   ! Renew the volume of the supercell:
   call Det_3x3(Scell(SCN)%supce,Scell(SCN)%V) ! module "Algebra_tools"

   ! Update relative coordinates:
   Scell(SCN)%MDAtoms(1:N_at) = MDAtoms(1:N_at)
   call Coordinates_abs_to_rel(Scell, SCN, .true.) ! module "Algebra_tools"

   !-----------------------
   ! 2) Place material (bio/molecule) centerred inside the box:
   do i = 1,3  ! x,y,z
      s_center(i) = 0.5d0 * ( maxval(Scell(SCN)%MDAtoms(:)%S(i)) + minval(Scell(SCN)%MDAtoms(:)%S(i)) )
   enddo
   ! Shift coordinates of all atoms:
   do i = 1, Scell(SCN)%Na
      Scell(SCN)%MDAtoms(i)%S(:) = 0.5d0 + (Scell(SCN)%MDAtoms(i)%S(:) - s_center(:))
      Scell(SCN)%MDAtoms(i)%S0 = Scell(SCN)%MDAtoms(i)%S ! previous timestep reset
   enddo
   ! Update absolute coordinates:
   call Coordinates_rel_to_abs(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"

   !-----------------------
   ! 3) Place remaining water molecules around the bio/molecule:
   N_tot = N_at   ! to start with
   j_h2o = N_at   ! to start with
   do i_h2o = 1, numpar%N_water_mol ! water molecules
      ! first, place O atom:
      j_h2o = j_h2o + 1 ! absolute number in the array of MDAtoms
      redo_placement = .true.  ! to start with
      iter = 1 ! count iteration of attempted placement
      do while (redo_placement)
         redo_placement = .false. ! assume we place it well, no need to redo it
         ! Set random coordinates of the new water molecule:
         call random_number(RN)  ! random numbers for relative coordinates along X,Y,Z
         Scell(SCN)%MDAtoms(j_h2o)%S(:) = RN(:)
         Scell(SCN)%MDAtoms(j_h2o)%S0(:) = Scell(SCN)%MDAtoms(j_h2o)%S(:)  ! previous timestep set
         call Coordinates_rel_to_abs_single(Scell, SCN, j_h2o, .true.)     ! module "Atomic_tools"

         ! Check that it is not overlapping with existing atoms:
         CHKI:do i = 1, N_tot
            ! Get the relative distance to this atom:
            if (j_h2o /= i) then
               !dS(:) = abs(Scell(SCN)%MDAtoms(j_h2o)%S(:) - Scell(SCN)%MDAtoms(i)%S(:))
               call shortest_distance(Scell(SCN), j_h2o, i, a_r) ! module "Atomic_tools"

               ! Check if it is not too short:
               if (a_r < dR_min) then ! overlapping atoms, place a molecule in a new place:
                  redo_placement = .true. ! assume we place it well, no need to redo it
                  iter = iter + 1   ! next iteration
!                   print*, 'CHKI', j_h2o, i, sqrt(SUM(dS(:)*dS(:)))
                  exit CHKI
               endif
            endif
         enddo CHKI
         ! If it's not possible to place a molecule in such a small volume, increase the volume:
!          print*, 'O  :', i_h2o, iter
         if (iter >= iter_max) then
            dV = dV_resc  ! increase the supercell volume
            goto 2023   ! try again placement in larger volume
         endif
      enddo
      ! Set nesessary atomic parameters:
      Scell(SCN)%MDatoms(j_h2o)%KOA = i_O

      ! then, place first H atom:
      j_h2o = j_h2o + 1
      redo_placement = .true.  ! to start with
      iter = 1 ! count iteration of attempted placement
      do while (redo_placement)
         redo_placement = .false. ! assume we place it well, no need to redo it
         ! Set random coordinates of the new water molecule:
         call random_number(RN)  ! random numbers for relative coordinates along X,Y,Z
         theta = 2.0d0*g_Pi*RN(1) ! angle
         phi = -g_half_Pi + g_Pi*RN(2)  ! second angle
         cos_phi = cos(phi)
         Scell(SCN)%MDAtoms(j_h2o)%R(1) = Scell(SCN)%MDAtoms(j_h2o-1)%R(1) + m_H2O_dist*cos_phi*cos(theta) ! [A]
         Scell(SCN)%MDAtoms(j_h2o)%R(2) = Scell(SCN)%MDAtoms(j_h2o-1)%R(2) + m_H2O_dist*cos_phi*sin(theta) ! [A]
         Scell(SCN)%MDAtoms(j_h2o)%R(3) = Scell(SCN)%MDAtoms(j_h2o-1)%R(3) + m_H2O_dist*sin(phi) ! [A]
         Scell(SCN)%MDAtoms(j_h2o)%R0(:) = Scell(SCN)%MDAtoms(j_h2o)%R(:)  ! previous timestep set

         call Coordinates_abs_to_rel_single(Scell, SCN, j_h2o, .true.)   ! module "Atomic_tools"

         ! Save cosine directions of the first H-O for the O-(second H) below:
         u0 = (Scell(SCN)%MDAtoms(j_h2o)%R(1) - Scell(SCN)%MDAtoms(j_h2o-1)%R(1))/m_H2O_dist
         v0 = (Scell(SCN)%MDAtoms(j_h2o)%R(2) - Scell(SCN)%MDAtoms(j_h2o-1)%R(2))/m_H2O_dist
         w0 = (Scell(SCN)%MDAtoms(j_h2o)%R(3) - Scell(SCN)%MDAtoms(j_h2o-1)%R(3))/m_H2O_dist

         if (abs(w0) > 1.0d0) then
            print*, 'ERROR in embed_molecule_in_water:'
            print*, 'BIG:', w0, j_h2o, SUM(Scell(SCN)%supce(1,:)*Scell(SCN)%supce(1,:))
            print*, 'N-1:', Scell(SCN)%MDAtoms(j_h2o-1)%R(:)
            print*, 'N  :', Scell(SCN)%MDAtoms(j_h2o)%R(:)
            print*, 'N1S:', Scell(SCN)%MDAtoms(j_h2o-1)%S(:)
            print*, 'N S:', Scell(SCN)%MDAtoms(j_h2o)%S(:)
         endif

         ! put atoms back into the supercell, if needed:
         call check_periodic_boundaries_single(matter, Scell, SCN, j_h2o)  ! module "Atomic_tools"

         ! Check that it is not overlapping with existing atoms:
         CHKI2:do i = 1, N_tot
            ! Get the relative distance to this atom:
            if (j_h2o /= i) then
!                dS(:) = abs(Scell(SCN)%MDAtoms(j_h2o)%S(:) - Scell(SCN)%MDAtoms(i)%S(:))
!                ! Check if it is not too short:
!                if ( sqrt( SUM(dS(:)*dS(:)) ) < dS_min_abs ) then ! overlapping atoms, place a molecule in a new place:

               call shortest_distance(Scell(SCN), j_h2o, i, a_r) ! module "Atomic_tools"

               ! Check if it is not too short:
               if (a_r < dR_min) then ! overlapping atoms, place a molecule in a new place:

                  redo_placement = .true. ! assume we place it well, no need to redo it
                  iter = iter + 1   ! next iteration
!                   print*, 'CHKI2', j_h2o, i, sqrt(SUM(dS(:)*dS(:))), dS_min_abs
!                   pause
                  exit CHKI2
               endif
            endif
         enddo CHKI2
         ! If it's not possible to place a molecule in such a small volume, increase the volume:
!          print*, 'H1 :', i_h2o, iter
         if (iter >= iter_max) then
            dV = dV_resc  ! increase the supercell volume
            goto 2023   ! try again placement in larger volume
         endif
      enddo
      ! Set nesessary atomic parameters:
      Scell(SCN)%MDatoms(j_h2o)%KOA = i_H

      ! next, place second H atom:
      j_h2o = j_h2o + 1
      redo_placement = .true.  ! to start with
      iter = 1 ! count iteration of attempted placement
      do while (redo_placement)
         redo_placement = .false. ! assume we place it well, no need to redo it
         ! Set coordinates of the new water molecule:
         call random_number(RN)  ! random numbers for relative coordinates along X,Y,Z
         theta = m_H2O_theta     ! angle H-O-H is fixed
         phi = -g_half_Pi + g_Pi*RN(2) ! second angle is random
         ! Cosine directions to the O-(second H), randomly oriented but with fixed H-O-H angle:
!          print*, 'v', u0, v0, w0
         call deflect_velosity(u0, v0, w0, theta, phi, u, v, w)   ! module "Atomic_tools"
         ! Coordinates of the second H atom:
         Scell(SCN)%MDAtoms(j_h2o)%R(1) = Scell(SCN)%MDAtoms(j_h2o-2)%R(1) + m_H2O_dist*u ! [A]
         Scell(SCN)%MDAtoms(j_h2o)%R(2) = Scell(SCN)%MDAtoms(j_h2o-2)%R(2) + m_H2O_dist*v ! [A]
         Scell(SCN)%MDAtoms(j_h2o)%R(3) = Scell(SCN)%MDAtoms(j_h2o-2)%R(3) + m_H2O_dist*w ! [A]
         Scell(SCN)%MDAtoms(j_h2o)%R0(:) = Scell(SCN)%MDAtoms(j_h2o)%R(:)  ! previous timestep set
         call Coordinates_abs_to_rel_single(Scell, SCN, j_h2o, .true.)   ! module "Atomic_tools"
         ! put atoms back into the supercell, if needed:
         call check_periodic_boundaries_single(matter, Scell, SCN, j_h2o)  ! module "Atomic_tools"

         ! Check that it is not overlapping with existing atoms:
         CHKI3:do i = 1, N_tot
            ! Get the relative distance to this atom:
            if (j_h2o /= i) then
!                dS(:) = abs(Scell(SCN)%MDAtoms(j_h2o)%S(:) - Scell(SCN)%MDAtoms(i)%S(:))
               call shortest_distance(Scell(SCN), j_h2o, i, a_r) ! module "Atomic_tools"
               ! Check if it is not too short:
               !if ( sqrt( SUM(dS(:)*dS(:)) ) < dS_min_abs ) then ! overlapping atoms, place a molecule in a new place:
               if (a_r < dR_min) then ! overlapping atoms, place a molecule in a new place:
                  redo_placement = .true. ! assume we place it well, no need to redo it
                  iter = iter + 1   ! next iteration
!                   print*, 'CHKI3', j_h2o, i, sqrt(SUM(dS(:)*dS(:)))
                  exit CHKI3
               endif
            endif
         enddo CHKI3
         ! If it's not possible to place a molecule in such a small volume, increase the volume:
!          print*, 'H2 :', i_h2o, iter
         if (iter >= iter_max) then
            dV = dV_resc  ! increase the supercell volume
            goto 2023   ! try again placement in larger volume
         endif
      enddo
      ! Set nesessary atomic parameters:
      Scell(SCN)%MDatoms(j_h2o)%KOA = i_H

      N_tot = N_tot + 3 ! we added one water molecule (3 new atoms)
   enddo

   !-----------------------
   ! 4) Redefine parameters of the supercell:
   Scell(SCN)%Na = N_at + N_h2o  ! new number of atoms (bio/molecule + water)
   Scell(SCN)%Ne = SUM(matter%Atoms(:)%NVB*matter%Atoms(:)%percentage)/SUM(matter%Atoms(:)%percentage)*Scell(SCN)%Na
   Scell(SCN)%Ne_low = Scell(SCN)%Ne ! at the start, all electrons are low-energy
   Scell(SCN)%Ne_high = 0.0d0 ! no high-energy electrons at the start
   Scell(SCN)%Ne_emit = 0.0d0 ! no emitted electrons at the start

   if (numpar%verbose) then
      write(*, '(a)', advance='no') 'Number of valence electrons: '
      write(*,*) matter%Atoms(:)%NVB
      write(*,*) ' (total: ', Scell(1)%Ne, ')'
   endif

!     print*, 'Ne ', matter%Atoms(:)%NVB, Scell(SCN)%Na, Scell(SCN)%Ne
!     print*, 'Per', matter%Atoms(:)%percentage, SUM(matter%Atoms(:)%percentage)
!      pause 'embed_molecule_in_water'

   ! Redefine velosities of all atoms:
   call set_initial_velocities(matter, Scell, 1, Scell(1)%MDatoms, numpar, numpar%allow_rotate) ! below

   ! save coords at equilibrium positions to ger mean square displacements later:
   do j = 1, Scell(SCN)%Na
      Scell(SCN)%MDatoms(j)%R_eq(:) = Scell(SCN)%MDatoms(j)%R0(:)
      Scell(SCN)%MDatoms(j)%S_eq(:) = Scell(SCN)%MDatoms(j)%S0(:)
   enddo ! j
   ! For Martyna algorithm (only start from zeros for now...):
   do i = 1, Scell(SCN)%Na
      Scell(SCN)%MDAtoms(i)%A = 0.0d0
      Scell(SCN)%MDAtoms(i)%A_tild(:) = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_F = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_J = 0.0d0
      Scell(SCN)%MDAtoms(i)%A0 = 0.0d0
      Scell(SCN)%MDAtoms(i)%A_tild0 = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_F0 = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_J0 = 0.0d0
   enddo

   deallocate(MDAtoms)
end subroutine embed_molecule_in_water



subroutine get_initial_atomic_coord(FN, File_name, Scell, SCN, which_one, matter, numpar, Err, ind)
   integer, intent(in) :: FN, which_one, SCN ! file number; type of file to read from (2=unit-cell, 1=super-cell); number of supercell
   character(*), intent(in) :: File_name ! file with the super-cell parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   integer, intent(in), optional :: ind ! read files for phase path tracing
   !=====================================
   integer :: INFO, Na, i, j
   integer Reason, count_lines
   character(200) Error_descript
   logical read_well
   type(Atom), dimension(:), allocatable :: MDAtoms ! if more then one supercell
   real(8), dimension(size(matter%Atoms)) :: perc
   real(8), dimension(3,3) :: unit_cell   ! temporary storage of supercell vectors

   select case (which_one)
   case (1) ! saved all atomic coordinates
      count_lines = 0
      call Count_lines_in_file(FN, Na) ! that's how many atoms we have
      
      if (present(ind)) then    ! reading data for two phases for path coordinate tracing
         if (ind == 0) then ! define the parameters:
            Scell(SCN)%Na = Na ! Number of atoms is defined this way
            Scell(SCN)%Ne = SUM(matter%Atoms(:)%NVB*matter%Atoms(:)%percentage)/SUM(matter%Atoms(:)%percentage)*Scell(SCN)%Na
            Scell(SCN)%Ne_low = Scell(SCN)%Ne ! at the start, all electrons are low-energy
            Scell(SCN)%Ne_high = 0.0d0 ! no high-energy electrons at the start
            Scell(SCN)%Ne_emit = 0.0d0 ! no emitted electrons at the start
         else ! do not redefine the parameters, but check for consistency
            if (Na /= Scell(SCN)%Na) then
               write(Error_descript,'(a)') 'Inconsistent numbers of atoms in the files PHASE_1_atoms.dat and PHASE_2_atoms.dat, terminating XTANT'
               call Save_error_details(Err, 3, Error_descript)
               print*, trim(adjustl(Error_descript))
               goto 3417
            endif
         endif
      else
         Scell(SCN)%Na = Na ! Number of atoms is defined this way
         !Scell(SCN)%Ne = matter%Atoms(1)%Ne_shell(matter%Atoms(1)%sh)*Scell(SCN)%Na	! number of valence electrons
         Scell(SCN)%Ne = SUM(matter%Atoms(:)%NVB*matter%Atoms(:)%percentage)/SUM(matter%Atoms(:)%percentage)*Scell(SCN)%Na
         Scell(SCN)%Ne_low = Scell(SCN)%Ne ! at the start, all electrons are low-energy
         Scell(SCN)%Ne_high = 0.0d0 ! no high-energy electrons at the start
         Scell(SCN)%Ne_emit = 0.0d0 ! no emitted electrons at the start
      endif
      
      if (.not.allocated(MDAtoms)) allocate(MDAtoms(Scell(SCN)%Na))
      if (.not.allocated(Scell(SCN)%MDAtoms)) allocate(Scell(SCN)%MDAtoms(Scell(SCN)%Na))
      
      do i = 1, Na ! read atomic data:
         !read(FN,*,IOSTAT=Reason) Scell(SCN)%MDAtoms(i)%KOA, Scell(SCN)%MDAtoms(i)%S(:), Scell(SCN)%MDAtoms(i)%S0(:), Scell(SCN)%MDAtoms(i)%SV(:), Scell(SCN)%MDAtoms(i)%SV0(:)
         read(FN,*,IOSTAT=Reason) MDAtoms(i)%KOA, MDAtoms(i)%S(:), MDAtoms(i)%S0(:), MDAtoms(i)%SV(:), MDAtoms(i)%SV0(:)
         call read_file(Reason, count_lines, read_well)
!          print*, 'Read:', i, Scell(SCN)%MDAtoms(i)%S(:), Scell(SCN)%MDAtoms(i)%SV(:)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript)), Reason
            goto 3417
         endif
      enddo
      
      if (present(ind)) then
         if (ind == 0) then ! initial phase
            do i = 1, Na ! read atomic data:
               Scell(SCN)%MDAtoms(i)%KOA = MDAtoms(i)%KOA
               Scell(SCN)%MDAtoms(i)%S0(:) = MDAtoms(i)%S0(:)
               Scell(SCN)%MDAtoms(i)%SV0(:) = MDAtoms(i)%SV0(:)
            enddo
         else   ! final phase
            do i = 1, Na ! read atomic data:
               Scell(SCN)%MDAtoms(i)%S(:) = MDAtoms(i)%S(:)
               Scell(SCN)%MDAtoms(i)%SV(:) = MDAtoms(i)%SV(:)
            enddo
            call Coordinates_rel_to_abs(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"
            call velocities_abs_to_rel(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"
         endif
      else
         Scell(SCN)%MDAtoms = MDAtoms
         call Coordinates_rel_to_abs(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"
         !call velocities_abs_to_rel(Scell, SCN, if_old=.true.)	! from the module "Atomic_tools"
         call velocities_rel_to_abs(Scell, SCN) ! get absolute velocities, module "Atomic_tools"
      endif
      
      if (.not.allocated(MDAtoms)) deallocate(MDAtoms)

   case (3) ! coordinates in the XYZ file
      ! Replicate unit cell, if requested:
      unit_cell = Scell(SCN)%Supce  ! use it to resize the supercell
      call set_supercell_size_from_unitcells(Scell, SCN, matter, unit_cell, .true.)   ! below

      ! Set coordinates in the sueprcell:
      call set_initial_coords(matter, Scell, SCN, FN, File_name, Nat=Scell(SCN)%Na, INFO=INFO, Error_descript=Error_descript, XYZ=1)
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         goto 3417
      endif

   case default ! coordinates in the unit cell
      call set_initial_coords(matter, Scell, SCN, FN, File_name, INFO=INFO, Error_descript=Error_descript)
      if (INFO .NE. 0) then
         call Save_error_details(Err, INFO, Error_descript)
         goto 3417
      endif

   end select


!    print*, 'NVB_1 = ', Scell(SCN)%Ne
   ! Check for consistency of valence electrons in chemical formula and actual atoms in the supercell:
   do j = 1, size(perc)
      perc(j) = COUNT( Scell(SCN)%MDatoms(:)%KOA == j )
!       print*, j, matter%Atoms(j)%percentage, perc(j)
      if (perc(j) /= matter%Atoms(j)%percentage) then ! assume supercell gives the right number
!          print*, 'Overwriting the number of VB electrons with the consistent number from supercell data'
         matter%Atoms(j)%percentage = perc(j)
      endif
   enddo
   Scell(SCN)%Ne = SUM(matter%Atoms(:)%NVB*matter%Atoms(:)%percentage)/SUM(matter%Atoms(:)%percentage)*Scell(SCN)%Na
   Scell(SCN)%Ne_low = Scell(SCN)%Ne ! at the start, all electrons are low-energy

   if (numpar%verbose) then
      write(*, '(a)', advance='no') 'Number of valence electrons: '
      !write(*, '(a)') 'Number of valence electrons: '
      write(*,*) dble(matter%Atoms(:)%NVB)
      write(*,*) '(total: ', Scell(1)%Ne, 'per atom:', Scell(1)%Ne/Scell(1)%Na, ')'
   endif

!    print*, 'NVB_2 = ', Scell(SCN)%Ne
!    pause
   
   ! Save atomic coordinates at their equilibrium positions:
   if (present(ind)) then
      if (ind == 0) then ! initial phase
          do j = 1, Scell(SCN)%Na
            Scell(SCN)%MDatoms(j)%R_eq(:) = Scell(SCN)%MDatoms(j)%R0(:)	! save coords at equilibrium positions to ger mean square displacements later
            Scell(SCN)%MDatoms(j)%S_eq(:) = Scell(SCN)%MDatoms(j)%S0(:)	! save coords at equilibrium positions to ger mean square displacements later
         enddo ! j
      else   ! final phase
            ! has nothing to do here
      endif
   else
      do j = 1, Scell(SCN)%Na
         Scell(SCN)%MDatoms(j)%R_eq(:) = Scell(SCN)%MDatoms(j)%R(:)	! save coords at equilibrium positions to ger mean square displacements later
         Scell(SCN)%MDatoms(j)%S_eq(:) = Scell(SCN)%MDatoms(j)%S(:)	! save coords at equilibrium positions to ger mean square displacements later
      enddo ! j
   endif
   
   ! For Martyna algorithm (only start from zeros for now...):
   do i = 1, Scell(SCN)%Na
      Scell(SCN)%MDAtoms(i)%A = 0.0d0
      Scell(SCN)%MDAtoms(i)%A_tild(:) = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_F = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_J = 0.0d0
      Scell(SCN)%MDAtoms(i)%A0 = 0.0d0
      Scell(SCN)%MDAtoms(i)%A_tild0 = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_F0 = 0.0d0
      Scell(SCN)%MDAtoms(i)%v_J0 = 0.0d0
   enddo
   
   
3417 continue

!    do i = 1, Scell(SCN)%Na
!       if (present(ind)) then
!          write(6,'(i4,f,f,f,f,f,f)') ind, Scell(SCN)%MDAtoms(i)%S0(:), Scell(SCN)%MDAtoms(i)%S(:)
!       else
!          write(6,'(a)') 'Subroutine get_initial_atomic_coord is called without ind'
!       endif
!    enddo ! j

end subroutine get_initial_atomic_coord



subroutine set_initial_coords(matter,Scell,SCN,FN,File_name,Nat,INFO,Error_descript,XYZ)
   type(solid), intent(inout) :: matter	! materil parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: SCN ! number of supercell
   integer, intent(in) :: FN	! file number for reading initial coordinates inside the unit-cell
   character(*), intent(in) :: File_name ! file name where to read from
   integer, intent(inout), optional :: Nat ! number of atoms in the unit cell
   integer, intent(out) :: INFO ! did we read well from the file
   character(200), intent(inout) :: Error_descript
   integer, intent(in), optional :: XYZ   ! if XYZ file was given and read from
   !-----------------------------
   real(8) a, b, l2, x, y, z, RN, epsylon, coord_shift
   integer i, j, k, ik, nx, ny, nz, ncellx, ncelly, ncellz, Na
   !real(8), dimension(3,8) :: Relcoat
   real(8), dimension(:,:), allocatable :: Relcoat
   integer, dimension(:), allocatable :: KOA
   integer Reason, count_lines
   logical read_well, to_read

   epsylon = 1.0d-10    ! for tiny shift of coords

   if (present(XYZ)) then
      select case (XYZ)
      case (1)
         to_read = .false. ! no need to read, data were already provided in XYZ file
      case default
         to_read = .true.  ! coordinates to be read from unit-cell file
      end select
   else  ! default: read from file
      to_read = .true.  ! coordinates to be read from unit-cell file
   endif
   
   INFO = 0 ! at the start there is no errors
   if (present(Nat)) then
      Na = Nat ! given number of atoms in the unit cell
   else
      call Count_lines_in_file(FN, Na) ! that's how many atoms we have
   endif
   allocate(Relcoat(3,Na))
   allocate(KOA(Na))

   Scell(SCN)%Na = Na*matter%cell_x*matter%cell_y*matter%cell_z ! Number of atoms is defined this way
   !Scell(SCN)%Ne = matter%Atoms(1)%Ne_shell(matter%Atoms(1)%sh)*Scell(SCN)%Na	! number of valence electrons
   Scell(SCN)%Ne = SUM(matter%Atoms(:)%NVB*matter%Atoms(:)%percentage)/SUM(matter%Atoms(:)%percentage)*Scell(SCN)%Na
   Scell(SCN)%Ne_low = Scell(SCN)%Ne ! at the start, all electrons are low-energy
   
   if (.not.allocated(Scell(SCN)%MDatoms)) allocate(Scell(SCN)%MDatoms(Scell(SCN)%Na))

   if (to_read) then
      count_lines = 0
      do i = 1, Na
         !read(FN,*,IOSTAT=Reason) Scell(SCN)%MDatoms(i)%KOA, Relcoat(:,i)	! relative coordinates of atoms in the unit-cell
         read(FN,*,IOSTAT=Reason) KOA(i), Relcoat(:,i)	! relative coordinates of atoms in the unit-cell
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            INFO = 3
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            print*, trim(adjustl(Error_descript))
            goto 3417
         endif
         !matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma
      enddo
   else  ! No need to read from unit-cell file
      ! Unit-cell coordinates were already read from XYZ file:
      do j = 1, Na
         Relcoat(:,j) = Scell(SCN)%MDatoms(j)%S(:)
         KOA(j) = Scell(SCN)%MDatoms(j)%KOA
         !print*, 'KOA', j, KOA(j), Scell(SCN)%MDatoms(j)%KOA
      enddo
      ! Update the size of the supercell from the unit cell:
      if (size(Scell(SCN)%MDatoms) /= Scell(SCN)%Na) then
         deallocate(Scell(SCN)%MDatoms)
         allocate(Scell(SCN)%MDatoms(Scell(SCN)%Na))
      endif
   endif

   ! All atoms distribution (in the super-cell):
   j = 0 ! for the beginning
   ncellx = matter%cell_x
   ncelly = matter%cell_y
   ncellz = matter%cell_z
   do nx = 0, ncellx-1
      do ny = 0, ncelly-1
         do nz = 0, ncellz-1
            do k = 1,Na	! number of atoms in the unit cell is FIXED
               j = j + 1
               do i=1,3
                  if (i .EQ. 1) then
                     a = REAL(nx)*1.0d0
                     l2 = REAL(ncellx)
                  endif
                  if (i .EQ. 2) then
                     a = REAL(ny)*1.0d0
                     l2 =  REAL(ncelly)
                  endif
                  if (i .EQ. 3) then
                     a = REAL(nz)*1.0d0
                     l2 =  REAL(ncellz)
                  endif
                  
                  ! Add a tiny shift to coordinates of replicated atoms:
                  call random_number(RN)
                  if (RN > 0.5d0) then
                     coord_shift = epsylon*RN
                  else
                     coord_shift = -epsylon*RN
                  endif
                  
                  Scell(SCN)%MDatoms(j)%S(i) = (Relcoat(i,k) + a + coord_shift)/l2  ! relative coordinates of an atom
                  !Scell(SCN)%MDatoms(j)%KOA = Scell(SCN)%MDatoms(k)%KOA ! kind of atom
                  Scell(SCN)%MDatoms(j)%KOA = KOA(k) ! kind of atom
               enddo ! i
            enddo ! k
         enddo ! nz
      enddo ! ny
   enddo ! nx
   deallocate(Relcoat, KOA)

   call check_periodic_boundaries(matter, Scell, SCN)   ! module "Atomic_tools"
!    call Coordinates_rel_to_abs(Scell, SCN)	! from the module "Atomic_tools"

   do j = 1, Scell(SCN)%Na
      Scell(SCN)%MDatoms(j)%R0(:) = Scell(SCN)%MDatoms(j)%R(:)	! coords at "previous" time step
      Scell(SCN)%MDatoms(j)%S0(:) = Scell(SCN)%MDatoms(j)%S(:)	! relative coords at "previous" time step
!       write(*,'(i3,i2,f,f,f,f,f,f)') j, Scell(SCN)%MDatoms(j)%KOA, Scell(SCN)%MDatoms(j)%R(:), Scell(SCN)%MDatoms(j)%S(:)
   enddo ! j
!    pause 'INITIAL CONDITIONS'

   ! If we want to shift all atoms:
   ! call Shift_all_atoms(matter, atoms, 0.25d0, 0.25d0, 0.25d0) ! from the module "Atomic_tools"
3417 continue

end subroutine set_initial_coords


subroutine set_initial_velocities(matter, Scell, NSC, atoms, numpar, allow_rotation)
   type(solid), intent(inout) :: matter	! materil parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of the super-cell
   type(Atom), dimension(:), intent(inout) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   logical, intent(in) :: allow_rotation ! remove angular momentum or not?
   !-------------------------------------------------------
   real(8) :: xr, SCVol, Xcm, Ycm, Zcm, vx, vy, vz, BigL(3), BigI(3,3), BigIinv(3,3)
   real(8) :: x0(3),r1, v0(3), rxv(3), omeg(3), Na, V_temp, Mass
   real(8), dimension(:), allocatable :: indices ! working array of indices
   integer i, j

   ! Set random velocities for all atoms:
   do i = 1,Scell(NSC)%Na ! velociteis of all atoms
      !print*, 'Mass', i, Scell(NSC)%MDatoms(i)%KOA
      Mass = matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma
      ! Get initial velocity
      call Get_random_velocity(Scell(NSC)%TaeV, Mass, atoms(i)%V(1), atoms(i)%V(2), atoms(i)%V(3), 2) ! module "Atomic_tools"
   enddo

   ! Set initial velocities for the super-cell vectors:
   Scell(NSC)%Vsupce = 0.0d0 ! this part of energy is in atoms at first, so let them relax!
   Scell(NSC)%Vsupce0 = Scell(NSC)%Vsupce ! Derivatives of Super-cell vectors (velocities) on last time-step

   !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   ! Find unconnected fragments of which the material is constructed (can be one piece too):
   call get_fragments_indices(Scell, NSC, numpar, atoms, matter, indices) ! module "Atomic_tools"
   
!     print*, 'indices', indices
!     pause
   
   !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   ! Eliminate initial angular momentum of the super-cell and each fragment inside:
   
!    do i = 1,Scell(NSC)%Na ! velociteis of all atoms
!       print*, 'set_initial_velocities 0', atoms(i)%S(:), atoms(i)%R(:)
!    enddo
   
   if (.not.allow_rotation) call remove_angular_momentum(NSC, Scell, matter, atoms, indices) ! module "Atomic_tools"
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   ! Eliminate total momentum of the center-of-mass and each fragment inside:
   call remove_momentum(Scell, NSC, matter, atoms, indices) ! module "Atomic_tools"
   
!    do i = 1,Scell(NSC)%Na ! velociteis of all atoms
!       print*, 'set_initial_velocities 1', atoms(i)%S(:), atoms(i)%R(:)
!    enddo
   
   ! Set relative velocities according to the new absolute ones:
   call velocities_abs_to_rel(Scell, NSC)
!    print*, 'set_initial_velocities'
end subroutine set_initial_velocities



subroutine get_supercell_vectors(FN, File_name, Scell, SCN, which_one, matter, Err, ind)
   integer, intent(in) :: FN, which_one, SCN ! file number; type of file to read from (2=unit-cell, 1=super-cell); number of supercell
   character(*), intent(in) :: File_name ! file with the super-cell parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(in) :: matter	! all material parameters
   type(Error_handling), intent(inout) :: Err	! error save
   integer, intent(in), optional :: ind     ! index into which variable to save the data
   !=========================================
   real(8), dimension(3,3) :: unit_cell
   real(8), dimension(3) :: temp_vec
   integer Reason, count_lines, i
   character(200) Error_descript
   logical read_well
   count_lines = 0
   select case (which_one)
   case (1) ! super-cell
      do i = 1,15 ! read supercell data:
      select case (i)
      case (1:3)   ! supercell
         if (present(ind)) then
            if (ind == 1) then   ! to read this index or not
               read(FN,*,IOSTAT=Reason) temp_vec
               !Scell(SCN)%supce(i,:) =  temp_vec
               Scell(SCN)%supce(:,i) =  temp_vec
            else
               read(FN,*,IOSTAT=Reason) ! do not read this index
            endif
         else
            read(FN,*,IOSTAT=Reason) temp_vec
            !Scell(SCN)%supce(i,:) =  temp_vec
            Scell(SCN)%supce(:,i) =  temp_vec
         endif
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3417
         endif
      case (5:7)   ! supercell0
         if (present(ind)) then
            if (ind == 0) then
               read(FN,*,IOSTAT=Reason) temp_vec
               !Scell(SCN)%supce0(i-4,:) = temp_vec
               Scell(SCN)%supce0(:,i-4) = temp_vec
            else
               read(FN,*,IOSTAT=Reason) ! do not read this index
            endif
         else
            read(FN,*,IOSTAT=Reason) temp_vec
            !Scell(SCN)%supce0(i-4,:) = temp_vec
            Scell(SCN)%supce0(:,i-4) = temp_vec
         endif
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3417
         endif
      case (9:11)  ! supce vel
         if (present(ind)) then
            if (ind == 1) then
               read(FN,*,IOSTAT=Reason) temp_vec
               !Scell(SCN)%Vsupce(i-8,:) = temp_vec
               Scell(SCN)%Vsupce(:,i-8) = temp_vec
            else
               read(FN,*,IOSTAT=Reason) ! do not read this index
            endif
         else
            read(FN,*,IOSTAT=Reason) temp_vec
            !Scell(SCN)%Vsupce(i-8,:) = temp_vec
            Scell(SCN)%Vsupce(:,i-8) = temp_vec
         endif
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3417
         endif
      case (13:15)  ! supce vel0
       if (present(ind)) then
            if (ind == 0) then
               read(FN,*,IOSTAT=Reason) temp_vec
               !Scell(SCN)%Vsupce0(i-12,:) = temp_vec
               Scell(SCN)%Vsupce0(:,i-12) = temp_vec
            else
               read(FN,*,IOSTAT=Reason) ! do not read this index
            endif
         else
            read(FN,*,IOSTAT=Reason) temp_vec
            !Scell(SCN)%Vsupce0(i-12,:) = temp_vec
            Scell(SCN)%Vsupce0(:,i-12) = temp_vec
         endif
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3417
         endif
      case (4,8,12)
         read(FN,*,IOSTAT=Reason) ! skip the lines separating the data sets
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
            call Save_error_details(Err, 3, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 3417
         endif
      end select
      enddo
   case default ! unit-cell
      read(FN,*,IOSTAT=Reason) unit_cell
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i3,a,$)') 'Could not read line ', count_lines, ' in file '//trim(adjustl(File_name))
         call Save_error_details(Err, 3, Error_descript)
         print*, trim(adjustl(Error_descript))
         goto 3417
      endif
      ! Adjust size of the supercell given the number of unit cells:
      call set_supercell_size_from_unitcells(Scell, SCN, matter, unit_cell, .false.)  ! below

   end select
   
    if (present(ind)) then  ! if we define two phases
       if (ind == 1) then
          Scell(SCN)%SCforce%rep = 0.0d0
          Scell(SCN)%SCforce%att = 0.0d0
          Scell(SCN)%SCforce%total = 0.0d0
          call Det_3x3(Scell(SCN)%supce, Scell(SCN)%V) ! finding initial volume of the super-cell, module "Algebra_tools"
          call Reciproc(Scell(SCN)%supce, Scell(SCN)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
          Scell(SCN)%supce_eq = Scell(SCN)%supce	! [A] equilibrium lengths of super-cell
       else
          ! do nothing for this index
       endif
    else    ! if we define one supercell
       Scell(SCN)%SCforce%rep = 0.0d0
       Scell(SCN)%SCforce%att = 0.0d0
       Scell(SCN)%SCforce%total = 0.0d0
       call Det_3x3(Scell(SCN)%supce, Scell(SCN)%V) ! finding initial volume of the super-cell, module "Algebra_tools"
       call Reciproc(Scell(SCN)%supce, Scell(SCN)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
       Scell(SCN)%supce_eq = Scell(SCN)%supce	! [A] equilibrium lengths of super-cell
    endif

3417 continue
end subroutine get_supercell_vectors


subroutine set_supercell_size_from_unitcells(Scell, SCN, matter, unit_cell, def_par)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: SCN    ! index of supercell (=1)
   type(Solid), intent(in) :: matter	! all material parameters
   real(8), dimension(3,3), intent(in) :: unit_cell   ! unit cell vectors to construct the supercell
   logical, intent(in) :: def_par   ! flag to define other parameters of the supercell
   !-----------------------
   Scell(SCN)%supce(:,1) = matter%cell_x*unit_cell(:,1)   ! [A] length of super-cell X
   Scell(SCN)%supce(:,2) = matter%cell_y*unit_cell(:,2)   ! [A] length of super-cell Y
   Scell(SCN)%supce(:,3) = matter%cell_z*unit_cell(:,3)   ! [A] length of super-cell Z
   !Scell(SCN)%supce(1,:) = matter%cell_x*unit_cell(1,:)  ! [A] length of super-cell X
   !Scell(SCN)%supce(2,:) = matter%cell_y*unit_cell(2,:)  ! [A] length of super-cell Y
   !Scell(SCN)%supce(3,:) = matter%cell_z*unit_cell(3,:)  ! [A] length of super-cell Z
   Scell(SCN)%supce0 = Scell(SCN)%supce   ! [A] length of super-cell on the previous time-step
   Scell(SCN)%Vsupce = 0.0d0  ! initial velocity is 0
   Scell(SCN)%Vsupce0 = 0.0d0 ! initial velocity is 0
   if (def_par) then ! define volume, reciprocal, etc.
      Scell(SCN)%SCforce%rep = 0.0d0
      Scell(SCN)%SCforce%att = 0.0d0
      Scell(SCN)%SCforce%total = 0.0d0
      call Det_3x3(Scell(SCN)%supce, Scell(SCN)%V) ! finding initial volume of the super-cell, module "Algebra_tools"
      call Reciproc(Scell(SCN)%supce, Scell(SCN)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
      Scell(SCN)%supce_eq = Scell(SCN)%supce	! [A] equilibrium lengths of super-cell
   endif
end subroutine set_supercell_size_from_unitcells



END MODULE Initial_configuration
