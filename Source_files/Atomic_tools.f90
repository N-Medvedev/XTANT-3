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
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau
!   "PENELOPE-2014 A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2014)
! This module includes some tools for performing vector algebra operations:
MODULE Atomic_tools
use Universal_constants
use Objects
use Algebra_tools, only : Cross_Prod, Invers_3x3, Matrix_Vec_Prod, Transpose_M, d_detH_d_h_a_b, Two_Matr_mult, Det_3x3
use Little_subroutines, only : Find_in_array_monoton
implicit none
PRIVATE

! this interface finds by itself which of the two subroutine to use depending on the array passed:
interface shortest_distance
   module procedure shortest_distance_NEW ! for Scell as single object
   module procedure shortest_distance_OLD ! for Scell as an array
end interface shortest_distance


public :: define_subcells, Maxwell_int_shifted, Coordinates_rel_to_abs, velocities_abs_to_rel, make_time_step_supercell, &
get_energy_from_temperature, distance_to_given_cell, make_time_step_atoms, Rescale_atomic_velocities, save_last_timestep, &
get_interplane_indices, get_near_neighbours, get_number_of_image_cells, pair_correlation_function, get_fraction_of_given_sort, &
Reciproc_rel_to_abs, total_forces, Potential_super_cell_forces, super_cell_forces, Convert_reciproc_rel_to_abs, get_kinetic_energy_abs, &
get_mean_square_displacement, Cooling_atoms, Coordinates_abs_to_rel, get_Ekin, make_time_step_supercell_Y4, make_time_step_atoms_M, &
remove_angular_momentum, get_fragments_indices, remove_momentum, make_time_step_atoms_Y4, check_periodic_boundaries, &
Make_free_surfaces, Coordinates_abs_to_rel_single, velocities_rel_to_abs, check_periodic_boundaries_single, &
Coordinates_rel_to_abs_single, deflect_velosity, Get_random_velocity, shortest_distance, cell_vectors_defined_by_angles


!=======================================
! Yoshida parameters for 4th order MD integrator:
! https://en.wikipedia.org/wiki/Leapfrog_integration#Yoshida_algorithms
real(8) :: m_c1, m_c2, m_c3, m_c4, m_d1, m_d2, m_d3, m_w0, m_w1, m_cr2

parameter(m_cr2 = 2.0d0**(1.0d0/3.0d0))
parameter(m_w0 = -m_cr2/(2.0d0 - m_cr2))
parameter(m_w1 = 1.0d0/(2.0d0 - m_cr2))
parameter(m_c1 = m_w1/2.0d0)
parameter(m_c2 = (m_w0 + m_w1)/2.0d0)
parameter(m_c3 = m_c2)
parameter(m_c4 = m_c1)
parameter(m_d1 = m_w1)
parameter(m_d2 = m_w0)
parameter(m_d3 = m_d1)
!=======================================


 contains



subroutine cell_vectors_defined_by_angles(a, b, c, alpha, beta, gamm, a_vec, b_vec, c_vec, INFO)
   ! Definitioin from: http://gisaxs.com/index.php/Unit_cell
   real(8), intent(in) :: a, b, c   ! absolute values of the supercell vectors
   real(8), intent(in) :: alpha, beta, gamm   ! angles bebtween the supercell vectors
   real(8), dimension(3), intent(out) :: a_vec, b_vec, c_vec   ! cell vectors constructed
   integer, intent(out) :: INFO  ! flag if something is wrong
   !----------------------
   real(8) :: cos_alpha, cos_beta, cos_gamma, sin_beta, sin_gamma, eps, arg2, arg

   INFO = 0 ! to start with no errors
   eps = 1.0d-12  ! precision

   cos_alpha = cos(alpha)
   cos_beta = cos(beta)
   cos_gamma = cos(gamm)
   sin_beta = sin(beta)
   sin_gamma = sin(gamm)

   if (sin_gamma < eps) then
      INFO = 1 ! cannot construct two-dimentional cell
      return   ! exit the subroutine
   endif

   arg2 = (cos_alpha - cos_beta*cos_gamma)/sin_gamma
   arg = 1.0d0 - cos_beta**2 - arg2**2
   if (arg < eps) then
      INFO = 2 ! cannot construct imaginary cell
      return   ! exit the subroutine
   endif

   ! By defenition, let a_vec is aligned along X:
   a_vec = (/a, 0.0d0, 0.0d0/)
   ! Vector b is then:
   b_vec(1) = b*cos_gamma
   b_vec(2) = b*sin_gamma
   b_vec(3) = 0.0d0
   ! Vector c is then:
   c_vec(1) = c*cos_beta
   c_vec(2) = c*arg2
   c_vec(3) = c*sqrt(arg)
end subroutine cell_vectors_defined_by_angles



subroutine define_subcells(Scell, numpar)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters
   !-----------------------
   real(8) :: dm
   integer :: i
   
   select case (numpar%lin_scal)
   case default ! no linear scaling, do nothing
   case (1) ! use linear scaling -> define define_subcell
      ! Get the cut off radius:
      call get_near_neighbours(Scell, numpar, include_vdW=.true., dm=dm)  ! below
      ! Define the number of subcells by the interaction radius:
      numpar%N_subcels(1) = FLOOR( Scell(1)%supce(1,1) / dm )
      numpar%N_subcels(2) = FLOOR( Scell(1)%supce(2,2) / dm )
      numpar%N_subcels(3) = FLOOR( Scell(1)%supce(3,3) / dm )
      ! Define the subcell sizes:
      if (.not.allocated(numpar%Subcell_coord_sx)) allocate(numpar%Subcell_coord_sx(numpar%N_subcels(1)))
      if (.not.allocated(numpar%Subcell_coord_sy)) allocate(numpar%Subcell_coord_sy(numpar%N_subcels(2)))
      if (.not.allocated(numpar%Subcell_coord_sz)) allocate(numpar%Subcell_coord_sz(numpar%N_subcels(3)))
      ! Save the coordinates for each subcell:
      do i = 1, numpar%N_subcels(1)
         numpar%Subcell_coord_sx(i) = dble(i)/dble(numpar%N_subcels(1))
      enddo
      do i = 1, numpar%N_subcels(2)
         numpar%Subcell_coord_sy(i) = dble(i)/dble(numpar%N_subcels(2))
      enddo
      do i = 1, numpar%N_subcels(3)
         numpar%Subcell_coord_sz(i) = dble(i)/dble(numpar%N_subcels(3))
      enddo
      ! Allocate parameters of the subcells:
      if (.not.allocated(Scell(1)%Subcell)) then
         allocate(Scell(1)%Subcell(numpar%N_subcels(1),numpar%N_subcels(2),numpar%N_subcels(3)))
      endif
      
   end select
end subroutine define_subcells
 
 
 
pure subroutine Find_outermost_atom(Scell, axis_ind, which_surf, N_ind)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: axis_ind	! along which direction: 1=X, 2=Y, 3=Z?
   character(*), intent(in) :: which_surf	! at which surface: UP of DOWN? ("Up" also means "Right" or a positive in case of X or Y surface)
   integer, intent(out) :: N_ind	! number of atom that is the outermost one
   !---------------
   integer :: NScell
   NScell = 1	! so far, only 1 supercell
   select case (which_surf)	! Which surface, bottom or top?
   case ('UP', 'U', 'up', 'u', 'Up')	! maximal coordinate: "Top"
      N_ind = transfer(MAXLOC(Scell(NScell)%MDAtoms(:)%S(axis_ind)), 1)	! generic trasfer is used to convert from one-dimensional array with size 1 into an integer number
   case ('DOWN', 'D', 'down', 'd', 'Down')	! minimal coordinate: "Bottom"
      N_ind = transfer(MINLOC(Scell(NScell)%MDAtoms(:)%S(axis_ind)), 1)	! generic trasfer is used to convert from one-dimensional array with size 1 into an integer number
   end select
end subroutine Find_outermost_atom

 
 
pure subroutine Make_free_surfaces(Scell, numpar, matter)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   type(solid), intent(in) :: matter	! material parameters
   !---------------
   integer :: N_at, i_at
   real(8) :: factr
   factr = 50.0d0	! how much empty space to add around the finite size sample
   
   N_at = size(Scell(1)%MDAtoms) ! number of atoms
   
   if (.not.numpar%r_periodic(1)) then		! free surface along X
      ! Expand the simulation box along this direction:
      !Scell(1)%supce(:,1) = Scell(1)%supce(:,1)*factr   ! incorrect
      !Scell(1)%supce0(:,1) = Scell(1)%supce(:,1)
      Scell(1)%supce(1,:) = Scell(1)%supce(1,:)*factr ! correct
      Scell(1)%supce0(1,:) = Scell(1)%supce(1,:)
      ! Rescale relative coordinates of atoms and place atoms into the middle of the simulation box:
      do i_at = 1, N_at
         call Coordinate_rescaling(Scell(1), i_at, 1, factr)
         call Velocity_rescaling(Scell(1), i_at, 1, factr)
      enddo
      call Shift_all_atoms(matter, Scell, 1, shx=-0.5d0, N_start=1, N_end=N_at)
   endif
   
   if (.not.numpar%r_periodic(2)) then		! free surface along Y
      ! Expand the simulation box along this direction:
      !Scell(1)%supce(:,2) = Scell(1)%supce(:,2)*factr   ! incorrect
      !Scell(1)%supce0(:,2) = Scell(1)%supce(:,2)
      Scell(1)%supce(2,:) = Scell(1)%supce(2,:)*factr ! correct
      Scell(1)%supce0(2,:) = Scell(1)%supce(2,:)
      ! Rescale relative coordinates of atoms and place atoms into the middle of the simulation box:
      do i_at = 1, N_at
         call Coordinate_rescaling(Scell(1), i_at, 2, factr)
         call Velocity_rescaling(Scell(1), i_at, 2, factr)
      enddo
      call Shift_all_atoms(matter, Scell, 1, shy=-0.5d0, N_start=1, N_end=N_at)
   endif
   
   if (.not.numpar%r_periodic(3)) then		! free surface along Z
      ! Expand the simulation box along this direction:
      !Scell(1)%supce(:,3) = Scell(1)%supce(:,3)*factr ! incorrect
      !Scell(1)%supce0(:,3) = Scell(1)%supce(:,3)
      Scell(1)%supce(3,:) = Scell(1)%supce(3,:)*factr ! correct
      Scell(1)%supce0(3,:) = Scell(1)%supce(3,:)
      ! Rescale relative coordinates of atoms and place atoms into the middle of the simulation box:
      do i_at = 1, N_at
         call Coordinate_rescaling(Scell(1), i_at, 3, factr)
         call Velocity_rescaling(Scell(1), i_at, 3, factr)
      enddo
      call Shift_all_atoms(matter, Scell, 1, shz=-0.5d0, N_start=1, N_end=N_at)
   endif
   
end subroutine Make_free_surfaces


pure subroutine Coordinate_rescaling(Scell, i_at, ind, factr)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: i_at, ind	! number of atom, coordinate index (x,y,z)
   real(8), intent(in) :: factr
   Scell%MDAtoms(i_at)%S(ind) = Scell%MDAtoms(i_at)%S(ind)/factr
   Scell%MDAtoms(i_at)%S0(ind) = Scell%MDAtoms(i_at)%S(ind)
end subroutine Coordinate_rescaling
 
pure subroutine Velocity_rescaling(Scell, i_at, ind, factr)
   type(Super_cell), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: i_at, ind	! number of atom, coordinate index (x,y,z)
   real(8), intent(in) :: factr
   Scell%MDAtoms(i_at)%SV(ind) = Scell%MDAtoms(i_at)%SV(ind)/factr
   Scell%MDAtoms(i_at)%SV0(ind) = Scell%MDAtoms(i_at)%SV(ind)
end subroutine Velocity_rescaling


subroutine Get_random_velocity(T, Mass, Vx, Vy, Vz, ind)
   real(8), intent(in) :: T ! [eV] Temperature to set the velocities accordingly
   real(8), intent(in) :: Mass ! [kg] mass of the atom
   real(8), intent(out) :: Vx, Vy, Vz ! velocities [A/fs]
   integer, intent(in) :: ind ! which distribution to use
   select case (ind)
   case (1) ! linear distributino of random values
      call Random_RN(T, Mass, Vx, Vy, Vz)
   case (2) ! Maxwellian distribution
      call Maxwell_RN(T, Mass, Vx, vy, Vz)
   end select
end subroutine Get_random_velocity

 
! Sample according to Maxwell distribution:
subroutine Maxwell_RN(T, Mass, Vx, vy, Vz)
   real(8), intent(in) :: T ! [eV] Temperature to set the velocities accordingly
   real(8), intent(in) :: Mass ! [kg] mass of the atom
   real(8), intent(out) :: Vx, Vy, Vz ! velocities [A/fs]
   real(8) RN(6) ! random numbers
   real(8) E, V, Pi2, theta, phi, cos_phi, cos2
   integer i
   Pi2 = g_Pi/2.0d0
   E = 1d10 ! just to start
   do while (E > T*10.0d0) ! exclude too high energies
      do i = 1,size(RN)
         call random_number(RN(i))
      enddo
      cos2 = cos(RN(3)*Pi2)
      cos2 = cos2*cos2
      E = T*(-log(RN(1)) - log(RN(2))*cos2) ! [eV]
      E = 2.0d0*E ! halp of the energy is kinetic, half potential, so double it to get the right temperature
      !print*, 'E=', E, T, log(RN(1)), log(RN(2)), cos(RN(3))
   enddo
   V = sqrt(E*2.0d0*g_e/Mass)*1d-5 ! [A/fs] absolute value of velocity
   theta = 2.0d0*g_Pi*RN(4) ! angle
   phi = -Pi2 + g_Pi*RN(5)  ! second angle
   cos_phi = cos(phi)
   if (RN(6) <= 0.33) then
      Vx = V*cos_phi*cos(theta)
      Vy = V*cos_phi*sin(theta)
      Vz = V*sin(phi)
   elseif(RN(6) <=0.67) then
      Vz = V*cos_phi*cos(theta)
      Vx = V*cos_phi*sin(theta)
      Vy = V*sin(phi)
   else
      Vy = V*cos_phi*cos(theta)
      Vz = V*cos_phi*sin(theta)
      Vx = V*sin(phi)
   endif
end subroutine Maxwell_RN


subroutine Random_RN(T, Mass, Vx, Vy, Vz)
   real(8), intent(in) :: T ! [eV] Temperature to set the velocities accordingly
   real(8), intent(in) :: Mass ! [kg] mass of the atom
   real(8), intent(out) :: Vx, Vy, Vz ! velocities [A/fs]
   real(8) V_temp, xr
   V_temp = sqrt(2.0d0*3.0d0*T*g_e/Mass)*1d10/1d15 ! [A/fs] average velocity corresponding to the given temperature
   !V_temp = 1.2d0*sqrt(2.0d0*3.0d0*Scell(NSC)%TaeV*g_e/matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma)*1d10/1d15 ! [A/fs] factor of 1.2 is to account for energy of supercell
   call random_number(xr)
   Vx = V_temp*(-1.0d0 + 2.0d0*xr) ! [A/fs] atomic velocity X
   call random_number(xr)
   Vy = V_temp*(-1.0d0 + 2.0d0*xr) ! [A/fs] atomic velocity Y
   call random_number(xr)
   Vz = V_temp*(-1.0d0 + 2.0d0*xr) ! [A/fs] atomic velocity Z
end subroutine Random_RN

!NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
! create a list of nearest neighbors:
subroutine get_near_neighbours(Scell, numpar, include_vdW, dm)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   logical, intent (in), optional :: include_vdW ! include vdW into consideration
   real(8), intent(out), optional :: dm ! [A] get cut-off radius
   !------------------------
   integer NSC
   real(8) rm, r(1), d(1), X
   do NSC = 1, size(Scell) ! for all supercells
      if (allocated(Scell(NSC)%TB_Hamil)) then
       ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:)) ! attractive part
         select type(ARRAY)
         type is (TB_H_Pettifor)
            r = maxval(ARRAY(:,:)%rm)
            rm = r(1)
         type is (TB_H_Fu)
            r = maxval(ARRAY(:,:)%rm)
            rm = r(1)
         type is (TB_H_Molteni)
            r = maxval(ARRAY(:,:)%rcut)
            d = maxval(ARRAY(:,:)%d)
            rm = r(1) + 10.0d0*d(1)
          type is (TB_H_NRL)
            r = maxval(ARRAY(:,:)%Rc)
            rm = r(1)*g_au2A	! convert from a.u. to Angstroms
          type is (TB_H_DFTB)
            r = maxval(ARRAY(:,:)%rcut)
            d = maxval(ARRAY(:,:)%d)
            rm = r(1) + 10.0d0*d(1)
          type is (TB_H_3TB)
            r = maxval(ARRAY(:,:)%rcut)
            d = maxval(ARRAY(:,:)%d)
            rm = r(1) + 10.0d0*d(1)
          type is (TB_H_BOP)
            rm = maxval(ARRAY(:,:)%rcut + ARRAY(:,:)%dcut)
          type is (TB_H_xTB)
            r = maxval(ARRAY(:,:)%rcut)
            d = maxval(ARRAY(:,:)%d)
            rm = r(1) + 10.0d0*d(1)
         end select
       END ASSOCIATE
      else
         rm = 10.0d0
         print*, 'TB Hamiltonian parameters are undefined, using default value:', rm
      endif

      if (allocated(Scell(NSC)%TB_Repuls)) then
       ASSOCIATE (ARRAY => Scell(NSC)%TB_Repuls(:,:)) ! repulsive part
         select type(ARRAY)
         type is (TB_Rep_Pettifor)
            r = maxval(ARRAY(:,:)%dm)
            if (rm < r(1)) rm = r(1) ! in case repulsive potential is longer than the attractive
         type is (TB_Rep_Fu)
            r = maxval(ARRAY(:,:)%dm)
            if (rm < r(1)) rm = r(1) ! in case repulsive potential is longer than the attractive
         type is (TB_Rep_Molteni)
            r = maxval(ARRAY(:,:)%rcut)
            d = maxval(ARRAY(:,:)%d)
            !if (rm < r(1) + 4.0d0*d(1)) rm = r(1) + 4.0d0*d(1)
            X = r(1) + 10.0d0*d(1)
            if (rm < X) rm = X
          type is (TB_Rep_DFTB)
             if (ARRAY(1,1)%ToP == 0) then  ! type of parameterization: 0=polinomial, 1=spline
                r = maxval(ARRAY(:,:)%rcut)
             else
                r = maxval(ARRAY(:,:)%rcut_spline)
             endif
             if (rm < r(1)) rm = r(1) ! in case repulsive potential is longer than the attractive
         end select
       END ASSOCIATE
      else
         print*, 'TB Repulsive parameters are undefined, using default value:', rm
      endif
      
      if (present(include_vdW) .and. allocated(Scell(NSC)%TB_Waals)) then ! if we have vdW potential defined
         r = maxval(Scell(NSC)%TB_Waals(:,:)%d0_cut)
         d = maxval(Scell(NSC)%TB_Waals(:,:)%dd_cut)
         X = r(1) + 10.0d0*d(1)
         if (rm < X) rm = X
      endif
      
      if (present(dm)) then
         dm = rm
      else
         call Find_nearest_neighbours(Scell, rm, numpar) ! get list of nearest neighbours from "Atomic_tools"
!          print*, 'rm', rm
!          pause
      endif
   enddo
end subroutine get_near_neighbours


subroutine get_number_of_image_cells(Scell, NSC, atoms, R_cut, Nx, Ny, Nz)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   real(8), intent(in) :: R_cut ! [A] cut-off distance for chosen potential
   integer, intent(out) :: Nx, Ny, Nz ! number of super-cells images that contribute to total energy
   !--------------------------
   real(8), dimension(3) :: zb ! super cell indices
   real(8) :: R
   
   ! Get the number of image cells in X-direction:
   !R = 0.0d0 ! just to start
   Nx = 0
   zb = (/1.0d0, 0.0d0, 0.0d0/)
   call distance_to_given_cell(Scell, NSC, atoms, zb, 1, 1, R) ! module "Atomic_tools"
   do while (R < 2.0d0*R_cut) ! do X-direction:
      Nx = Nx + 1 ! check next image of the super-cell in X-derection
      zb(1) = dble(Nx) ! create a vector of super-cell image numbers
      ! Check whether the distance is without cut-off or not (chech for atom #1 distance to its own image),
      ! if it holds for arbitrary number of atom, it holds for all of them (arbitrary chosen atom #1):
      call distance_to_given_cell(Scell, NSC, atoms, zb, 1, 1, R) ! module "Atomic_tools"
   enddo
   
   ! Get the number of image cells in Y-direction:
   !R = 0.0d0 ! just to start
   Ny = 0
   zb = (/0.0d0, 1.0d0, 0.0d0/)
   call distance_to_given_cell(Scell, NSC, atoms, zb, 1, 1, R) ! module "Atomic_tools"
   do while (R < 2.0d0*R_cut) ! do Y-direction:
      Ny = Ny + 1 ! check next image of the super-cell in Y-derection
      zb(2) = dble(Ny) ! create a vector of super-cell image numbers
      ! Check whether the distance is without cut-off or not (chech for atom #1 distance to its own image),
      ! if it holds for arbitrary number of atom, it holds for all of them (arbitrary chosen atom #1):
      call distance_to_given_cell(Scell, NSC, atoms, zb, 1, 1, R) ! module "Atomic_tools"
   enddo
   
   ! Get the number of image cells in Z-direction:
   !R = 0.0d0 ! just to start
   Nz = 0
   zb = (/0.0d0, 0.0d0, 1.0d0/)
   call distance_to_given_cell(Scell, NSC, atoms, zb, 1, 1, R) ! module "Atomic_tools"
   do while (R < 2.0d0*R_cut) ! do Z-direction:
      Nz = Nz + 1 ! check next image of the super-cell in Z-derection
      zb(3) = dble(Nz) ! create a vector of super-cell image numbers
      ! Check whether the distance is without cut-off or not (chech for atom #1 distance to its own image),
      ! if it holds for arbitrary number of atom, it holds for all of them (arbitrary chosen atom #1):
      call distance_to_given_cell(Scell, NSC, atoms, zb, 1, 1, R) ! module "Atomic_tools"
   enddo
end subroutine get_number_of_image_cells


subroutine remove_momentum(Scell, NSC, matter, atoms, indices, print_out)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of the super-cell
   type(solid), intent(inout), target :: matter	! materil parameters
   type(Atom), dimension(:), intent(inout) :: atoms	! array of atoms in the supercell
   real(8), dimension(:), intent(in), optional :: indices ! working array of indices
   logical, optional, intent(in) :: print_out ! print our c-o-m momemntum
   real(8), pointer :: Mass
   real(8) :: vx, vy, vz, Masstot
   integer i, Na
   Na = size(atoms) ! number of atoms
   Masstot = SUM(matter%Atoms(atoms(:)%KOA)%Ma) ! net-mass of all atoms
   ! Center of mass velocities:
   vx = SUM(matter%Atoms(atoms(:)%KOA)%Ma * atoms(:)%V(1))/Masstot
   vy = SUM(matter%Atoms(atoms(:)%KOA)%Ma * atoms(:)%V(2))/Masstot
   vz = SUM(matter%Atoms(atoms(:)%KOA)%Ma * atoms(:)%V(3))/Masstot

   if (present(print_out)) then
      write(*,'(a,es25.16,es25.16,es25.16)') 'CM1:', vx, vy, vz
   endif
   ! Subtract the velocity of the center of mass:
   atoms(:)%V(1) = atoms(:)%V(1) - vx
   atoms(:)%V(2) = atoms(:)%V(2) - vy
   atoms(:)%V(3) = atoms(:)%V(3) - vz
   
   ! Remove angular momenta of individual unconnected fragments if any:
   if (present(indices).and. (maxval(indices(:)) > 1)) then ! special treatement for graphite
      call remove_plane_momenta(Scell, NSC, matter, atoms, indices)
   endif

   nullify(Mass)
end subroutine remove_momentum


subroutine remove_plane_momenta(Scell, NSC, matter, atoms, indices, print_out)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of the super-cell
   type(solid), intent(inout) :: matter	! materil parameters
   type(Atom), dimension(:), intent(inout) :: atoms	! array of atoms in the supercell
   real(8), dimension(:), intent(in) :: indices ! working array of indices
   logical, intent(in), optional :: print_out ! just in case user want to print out c-of-m-momenta
   real(8) :: vx, vy, vz, Mass_sum
   integer i, j, N(1), Nplane
   integer :: Na
   Na = size(atoms) ! number of atoms
   ! Locate atoms within planes:
!    indices = 0
!    call get_interplane_indices(Scell, NSC, atoms, matter, indices) ! function see below

   N = maxval(indices(:)) ! total number of planes
   Nplane = N(1)
   do i = 1, Nplane
      Mass_sum = sum(matter%Atoms(atoms(:)%KOA)%Ma, MASK = (indices(:) == i)) ! total mass of atoms within this plane
      vx = sum(atoms(:)%V(1)*matter%Atoms(atoms(:)%KOA)%Ma, MASK = (indices(:) == i)) ! Px of center of mass of this plane
      vy = sum(atoms(:)%V(2)*matter%Atoms(atoms(:)%KOA)%Ma, MASK = (indices(:) == i)) ! Py of center of mass of this plane
      vz = sum(atoms(:)%V(3)*matter%Atoms(atoms(:)%KOA)%Ma, MASK = (indices(:) == i)) ! Pz of center of mass of this plane
      vx = vx/Mass_sum ! Px -> Vx
      vy = vy/Mass_sum ! Py -> Vy
      vz = vz/Mass_sum ! Pz -> Vz
      if (present(print_out)) then
         write(*,'(a,i2,es25.16,es25.16,es25.16)') 'CM2:', i, vx, vy, vz
      endif
      do j = 1, Na ! for all atoms
        if (indices(j) == i) then! do only those within this plane
         ! subtract center-of-mass velocities:
         atoms(j)%V(1) = atoms(j)%V(1) - vx
         atoms(j)%V(2) = atoms(j)%V(2) - vy
         atoms(j)%V(3) = atoms(j)%V(3) - vz
        endif
      enddo
   enddo
end subroutine remove_plane_momenta


subroutine remove_angular_momentum(NSC, Scell, matter, atoms, indices, print_out)
   integer, intent(in) :: NSC ! number of the super-cell
   type(Super_cell), dimension(:), intent(inout), target :: Scell ! suoer-cell with all the atoms inside
   type(solid), intent(inout), target :: matter	! materil parameters
   type(Atom), dimension(:), intent(inout) :: atoms	! array of atoms in the supercell
   real(8), dimension(:), intent(in), optional :: indices ! working array of indices
   logical, optional, intent(in) :: print_out ! print our c-o-m momemntum
   !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   real(8), pointer :: Mass
   real(8) :: Masstot, Xcm, Ycm, Zcm, BigL(3), BigI(3,3), BigIinv(3,3), detB
   real(8) :: r1, v0(3), rxv(3), omeg(3)
   real(8), dimension(:,:), allocatable ::  x0
   integer i, j
   integer, pointer :: Na
   Na => Scell(NSC)%Na ! number of atoms
   allocate(x0(Na,3))
   x0 = 0.0d0
   
   ! Eliminate initial angular momentum of the super-cell:

   ! Finding the center of inertia NOT assuming equal masses of particles:
   Masstot = SUM(matter%Atoms(atoms(:)%KOA)%Ma) ! total mass
   Xcm = SUM(matter%Atoms(atoms(:)%KOA)%Ma * atoms(:)%R(1))/Masstot
   Ycm = SUM(matter%Atoms(atoms(:)%KOA)%Ma * atoms(:)%R(2))/Masstot
   Zcm = SUM(matter%Atoms(atoms(:)%KOA)%Ma * atoms(:)%R(3))/Masstot

   ! Calculating the total moment of inertia tensor:
   BigL = 0.0d0 ! 3d vector
   BigI = 0.0d0 ! (3x3) tensor
   do i = 1,Na
      ! Atomic positions in the center-of-mass frame of reference:
      x0(i,1) = atoms(i)%R(1) - Xcm !/dble(Na)
      x0(i,2) = atoms(i)%R(2) - Ycm !/dble(Na)
      x0(i,3) = atoms(i)%R(3) - Zcm !/dble(Na)
      r1 = DSQRT(x0(i,1)*x0(i,1) + x0(i,2)*x0(i,2) + x0(i,3)*x0(i,3))
      Mass => matter%Atoms(atoms(i)%KOA)%Ma
      rxv = 0.0d0
      call Cross_Prod(x0(i,:), atoms(i)%V(:), rxv) ! from MODULE "Algebra_tools"
      BigL(1) = BigL(1) + rxv(1)*Mass
      BigL(2) = BigL(2) + rxv(2)*Mass
      BigL(3) = BigL(3) + rxv(3)*Mass
      ! calculating total moment of inertia tensor:
      BigI(1,1) = BigI(1,1) + (r1*r1 - x0(i,1)*x0(i,1))*Mass ! xx
      BigI(2,2) = BigI(2,2) + (r1*r1 - x0(i,2)*x0(i,2))*Mass ! yy
      BigI(3,3) = BigI(3,3) + (r1*r1 - x0(i,3)*x0(i,3))*Mass ! zz
      BigI(1,2) = BigI(1,2) - x0(i,1)*x0(i,2)*Mass           ! xy
      BigI(1,3) = BigI(1,3) - x0(i,1)*x0(i,3)*Mass           ! xz
      BigI(2,3) = BigI(2,3) - x0(i,2)*x0(i,3)*Mass           ! yz
   enddo
   BigI(2,1) = BigI(1,2)                         ! yx
   BigI(3,1) = BigI(1,3)                         ! zx
   BigI(3,2) = BigI(2,3)                         ! zy

   call Det_3x3(BigI,detB) ! find determinant of A, module "Algebra_tools"
   if (detB > 0.0d0) then	! only if there is angular momentum
      call Invers_3x3(BigI, BigIinv, 'remove_angular_momentum') ! calculate inverse tensor of inertia ! from MODULE "Algebra_tools"

      call Matrix_Vec_Prod(BigIinv,BigL,omeg) ! find omega - the angular velocity ! from MODULE "Algebra_tools"
      rxv = 0.0d0
      do i = 1, Na ! to subtract the total angular velocity from each atom:
         call Cross_Prod(omeg, x0(i,:), rxv) ! from MODULE "Algebra_tools"
         atoms(i)%V(:) = atoms(i)%V(:) - rxv(:)
         if (present(print_out)) then
            write(*,'(a,i4,es25.16,es25.16,es25.16)') 'AM2:', i, rxv(:)
         endif
      enddo
   endif
   
   ! Remove angular momenta of individual unconnected fragments if any:
   if (present(indices) .and. (maxval(indices(:)) > 1)) then ! special treatement for graphite
      call remove_plane_angular_momenta(Scell, NSC, matter, atoms, indices)
   endif

   deallocate(x0)
   nullify(Mass, Na)
end subroutine remove_angular_momentum


subroutine remove_plane_angular_momenta(Scell, NSC, matter, atoms, indices, print_out)
   integer, intent(in) :: NSC ! number of the super-cell
   type(Super_cell), dimension(:), intent(inout), target :: Scell ! suoer-cell with all the atoms inside
   type(solid), intent(inout), target :: matter	! materil parameters
   type(Atom), dimension(:), intent(inout) :: atoms	! array of atoms in the supercell
   real(8), dimension(:), intent(in) :: indices ! working array of indices
   logical, optional, intent(in) :: print_out ! print our c-o-m momemntum
   !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   real(8), pointer :: Mass
   integer, pointer :: Na
   real(8) :: Masstot, Xcm, Ycm, Zcm, BigL(3), BigI(3,3), BigIinv(3,3)
   real(8) :: r1, v0(3), rxv(3), omeg(3), detB
   real(8), dimension(:,:), allocatable ::  x0
   integer i, j, N(1), i_plane, Nplane, Na_plane
   Na => Scell(NSC)%Na ! number of atoms
   allocate(x0(Na,3))
   x0 = -1.0d25
   
   N = maxval(indices(:)) ! total number of planes
   Nplane = N(1)
   PLANES:do i_plane = 1, Nplane ! do for each plane
      ! Eliminate initial angular momentum of the super-cell:

      Na_plane = count(MASK = (indices(:) == i_plane)) ! total mass of atoms within this plane
      ! Finding the center of inertia NOT assuming equal masses of particles:
      Xcm = 0.0d0 !SUM(atoms(:)%R(1)*matter%Atoms(atoms(:)%KOA)%Ma)/Masstot
      Ycm = 0.0d0 !SUM(atoms(:)%R(2)*matter%Atoms(atoms(:)%KOA)%Ma)/Masstot
      Zcm = 0.0d0 !SUM(atoms(:)%R(3)*matter%Atoms(atoms(:)%KOA)%Ma)/Masstot
      do i = 1,Na
         if (indices(i) == i_plane) then ! only atoms from this plane
            Mass => matter%Atoms(atoms(i)%KOA)%Ma
            Xcm = Xcm + atoms(i)%R(1)*Mass
            Ycm = Ycm + atoms(i)%R(2)*Mass
            Zcm = Zcm + atoms(i)%R(3)*Mass
         endif
      enddo
      Masstot = sum(matter%Atoms(atoms(:)%KOA)%Ma, MASK = (indices(:) == i_plane)) ! total mass of atoms within this plane
!       print*, i_plane, Nplane
!       print*, Masstot, matter%Atoms(atoms(:)%KOA)%Ma
!       pause 'Masstot'
      Xcm = Xcm/Masstot
      Ycm = Ycm/Masstot
      Zcm = Zcm/Masstot

      ! Calculating the total moment of inertia tensor:
      BigL = 0.0d0 ! 3d vector
      BigI = 0.0d0 ! (3x3) tensor
      do i = 1,Na
         if (indices(i) == i_plane) then ! only atoms from this plane
            x0(i,1) = atoms(i)%R(1) - Xcm
            x0(i,2) = atoms(i)%R(2) - Ycm
            x0(i,3) = atoms(i)%R(3) - Zcm
            if (x0(i,1) < -1.d10) print*, 'PANIC ATTACK!', i, x0(i,:)
            r1 = DSQRT(x0(i,1)*x0(i,1) + x0(i,2)*x0(i,2) + x0(i,3)*x0(i,3))
            rxv = 0.0d0
            call Cross_Prod(x0(i,:), atoms(i)%V(:), rxv) ! from MODULE "Algebra_tools"
            Mass => matter%Atoms(atoms(i)%KOA)%Ma
            BigL(:) = BigL(:) + rxv(:)*Mass
            ! calculating total moment of inertia tensor:
            BigI(1,1) = BigI(1,1) + (r1*r1 - x0(i,1)*x0(i,1))*Mass ! xx
            BigI(2,2) = BigI(2,2) + (r1*r1 - x0(i,2)*x0(i,2))*Mass ! yy
            BigI(3,3) = BigI(3,3) + (r1*r1 - x0(i,3)*x0(i,3))*Mass ! zz
            BigI(1,2) = BigI(1,2) - x0(i,1)*x0(i,2)*Mass           ! xy
            BigI(1,3) = BigI(1,3) - x0(i,1)*x0(i,3)*Mass           ! xz
            BigI(2,3) = BigI(2,3) - x0(i,2)*x0(i,3)*Mass           ! yz
         endif
      enddo
      BigI(2,1) = BigI(1,2)                         ! yx
      BigI(3,1) = BigI(1,3)                         ! zx
      BigI(3,2) = BigI(2,3)                         ! zy

      call Det_3x3(BigI,detB) ! find determinant of A, module "Algebra_tools"
      if (detB > 0.0d0) then	! only if there is angular momentum
         call Invers_3x3(BigI,BigIinv, 'remove_plane_angular_momenta') ! calculate inverse tensor of inertia ! from MODULE "Algebra_tools"

         call Matrix_Vec_Prod(BigIinv,BigL,omeg) ! find omega - the angular velocity ! from MODULE "Algebra_tools"
         rxv = 0.0d0
         do i = 1, Na ! to subtract the total angular velocity from each atom:
            if (indices(i) == i_plane) then ! only atoms from this plane
               call Cross_Prod(omeg, x0(i,:), rxv) ! from MODULE "Algebra_tools"
               atoms(i)%V(:) = atoms(i)%V(:) - rxv(:)
               if (present(print_out)) then
                  write(*,'(a,i4,i2,es25.16,es25.16,es25.16)') 'AM2:', i, i_plane, rxv(:)
               endif
            endif
         enddo
      endif
   enddo PLANES

   deallocate(x0)
   nullify(Mass, Na)
end subroutine remove_plane_angular_momenta



subroutine get_fragments_indices(Scell, NSC, numpar, atoms, matter, indices)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of the super-cell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   type(solid), intent(in) :: matter	! materil parameters
   type(Atom), dimension(:), intent(in) :: atoms   ! array of atoms in the supercell
   real(8), dimension(:), intent(inout), allocatable :: indices ! working array of indices
   real(8) a_r, dm
   integer i, j, coun, ind_i, ind_same, Na
   Na = size(atoms) ! corresponding to number of atoms
   if (.not.allocated(indices)) then 
      allocate(indices(Na))
      indices = 0.0d0
   endif
   coun = 0
   call get_near_neighbours(Scell, numpar, dm = dm) ! get cut-off radius
   
!    print*, 'dm', dm
   
   do i = 1, Na ! Check each atom - to which fragment it belongs
      if (indices(i)==0) then ! this atom is not allocated within any fragment yet
         coun = coun + 1   ! next fragment
         indices(i) = coun ! atom is within this fragment
      endif
      ind_i = indices(i)
      do j = 1, Na ! check if there are other atoms within the same fragment
         if (i /= j) then
            call shortest_distance(Scell, NSC, atoms, i, j, a_r)
            if (a_r < dm) then ! this atom is within the same fragment as the atom 'i'
               if (indices(j) == 0) then ! this atom has not been allocated to a fragment yet:
                  indices(j) = indices(i) ! mark this atoms with the same index
               else if (indices(j) == indices(i)) then ! this atom is already within this fragment
                  ! do nothing, it's already been done
               else ! it turns out, it's a part of the same fragment, who could've guessed?!
                  ! change indices of the second fragment to the first fragment, since it's actually the same:
                  ind_same = indices(j)
                  ind_i = MIN(ind_i,ind_same)
                  where(indices == ind_same) indices = ind_i ! renumber these atoms to this fragment
               endif
            endif
         endif
      enddo
   enddo
   
   !Check that there are no empty fragments left:
   do i = 1, maxval(indices(:))
      j = count(indices(:) == i)	! that's how many atoms are in this fragment
      do while (j < 1)
         where(indices >= i) indices = indices - 1	! shift them by one
         if (i >= maxval(indices(:))) exit	! all indices are shifted down to here, nothing to do more
         j = count(indices(:) == i)	! that's how many atoms are in this fragment
!       print*, 'i,j', i, j, maxval(indices(:))
      enddo
   enddo
   
!    pause
   
end subroutine get_fragments_indices


subroutine get_interplane_indices(Scell, NSC, numpar, atoms, matter, indices)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of the super-cell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   type(solid), intent(in) :: matter	! materil parameters
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   real(8), dimension(:), intent(inout) :: indices ! working array of indices
   real(8) a_r, x1, y1, z1, dm
   integer i, j, coun, Na
   Na = size(indices)
   coun = 0
   call get_near_neighbours(Scell, numpar, dm = dm) ! get cut-off distance
   do i = 1, Na-1 ! Check each atom - to which plane it belongs
      if (indices(i)==0) then ! this atom is not allocated within any plane yet
         coun = coun + 1   ! next plane
         indices(i) = coun ! atom is within this plane
         do j = i+1, Na ! check if there are other atoms within the same plane
            call shortest_distance(Scell, NSC, atoms, i, j, a_r, x1=x1, y1=y1, z1=z1)
            if (ABS(z1) < dm) then ! this atom is within the same plane as the atom 'i'
               indices(j) = indices(i) ! mark this atoms with the same index
            endif
         enddo
      endif
   enddo
end subroutine get_interplane_indices



subroutine get_fraction_of_given_sort(Scell, KOA, N_KOA)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: KOA	! index of kind of atom
   real(8), intent(out) :: N_KOA	! fraction of atoms of this kind are in the supercell normalized to the total number of atoms
   integer :: i1, Nat
   Nat = size(Scell(1)%MDatoms)	! total number of atoms
   i1 = COUNT(Scell(1)%MDatoms(:)%KOA == KOA)	! that's how many atoms of the given kind are in the supercell
   N_KOA = dble(i1)/dble(Nat)
end subroutine get_fraction_of_given_sort




subroutine save_last_timestep(Scell)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer i, NSC
   do NSC = 1, size(Scell)
      do i = 1,size(Scell(NSC)%MDatoms)
         Scell(NSC)%MDatoms(i)%R0 = Scell(NSC)%MDatoms(i)%R
         Scell(NSC)%MDatoms(i)%S0 = Scell(NSC)%MDatoms(i)%S
         Scell(NSC)%MDatoms(i)%V0 = Scell(NSC)%MDatoms(i)%V
         Scell(NSC)%MDatoms(i)%SV0 = Scell(NSC)%MDatoms(i)%SV
         Scell(NSC)%MDatoms(i)%A0 = Scell(NSC)%MDatoms(i)%A
         !atoms(i)%Ekin = atoms(i)%Ekin
      enddo
      Scell(NSC)%supce0 = Scell(NSC)%supce
      Scell(NSC)%Vsupce0 = Scell(NSC)%Vsupce
   enddo
end subroutine save_last_timestep


subroutine Cooling_atoms(numpar, matter, Scell, t_time, dt, t_start, incl)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(Solid), intent(inout) :: matter ! material parameters
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters
   real(8), intent(in) :: t_time ! [fs] current timestep
   real(8), intent(in) :: dt ! [fs] how often to cool it
   real(8), intent(in) :: t_start ! [fs] when to start cooling
   logical, intent(in) :: incl	! yes or no
   integer NSC
   do NSC = 1, size(Scell)
      call Cooling_atoms_SC(numpar, Scell(NSC)%nrg, matter, Scell, NSC, t_time, dt, t_start, incl) ! below
   enddo
end subroutine Cooling_atoms


subroutine Cooling_atoms_SC(numpar, nrg, matter, Scell, NSC, t_time, dt, t_start, incl)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(inout) :: matter ! material parameters
   type(Energies), intent(inout) :: nrg	! all energies
   type(Numerics_param), intent(inout) :: numpar	! numerical parameters
   real(8), intent(in) :: t_time ! [fs] current timestep
   real(8), intent(in) :: dt ! [fs] how often to cool it
   real(8), intent(in) :: t_start ! [fs] when to start cooling
   logical, intent(in) :: incl	! yes or no
   integer k
   if (incl) then ! do this
      if ((t_time .GE. t_start) .AND. (numpar%dt_cooling .GT. dt)) then ! it's time to cool atoms down:
         numpar%dt_cooling = 0.0d0
         ! Cooling atoms instantly to Ta=0:
         do k = 1, size(Scell(NSC)%MDatoms) ! for all atoms: instant cool-down
            Scell(NSC)%MDatoms(k)%SV(:) = 0.0d0
            Scell(NSC)%MDatoms(k)%SV0(:) = 0.0d0
         enddo
         ! Cooling supercell:
         Scell(NSC)%Vsupce(:,:) = 0.0d0
         Scell(NSC)%Vsupce0(:,:) = 0.0d0
         ! Adjust absolute velocities accordingly:
         call velocities_rel_to_abs(Scell, NSC)
         call Atomic_kinetic_energies(Scell, NSC, matter)
         call get_kinetic_energy_abs(Scell, NSC, matter, nrg)
      else
         numpar%dt_cooling = numpar%dt_cooling + numpar%dt ! [fs]
      endif
      if (numpar%verbose) print*, 'Quenching succesful : Cooling_atoms_SC'
   endif
end subroutine Cooling_atoms_SC


pure subroutine Shift_all_atoms(matter, Scell, NSC, shx, shy, shz, N_start, N_end)
   type(solid), intent(in) :: matter     ! materil parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   real(8), intent(in), optional :: shx, shy, shz  ! shift all atoms by this much
   integer, intent(in), optional :: N_start, N_end ! from this to that atom only
   integer k, Ms, Me
   
   if (present(N_start)) then
      Ms = N_start ! starting from this atom
   else
      Ms = 1 ! start from the first one by default
   endif
   if (present(N_end)) then
      Me = N_end ! end at this atom
   else
      Me = Scell(NSC)%Na ! all atoms by default
   endif
   
   do k = Ms, Me
      if (present(shx)) then
         ! Shift along x direction:
         Scell(NSC)%MDatoms(k)%S(1) = Scell(NSC)%MDatoms(k)%S(1) - shx
         Scell(NSC)%MDatoms(k)%S0(1) = Scell(NSC)%MDatoms(k)%S0(1) - shx
      endif
      if (present(shy)) then
         ! Shift along y direction:
         Scell(NSC)%MDatoms(k)%S(2) = Scell(NSC)%MDatoms(k)%S(2) - shy
         Scell(NSC)%MDatoms(k)%S0(2) = Scell(NSC)%MDatoms(k)%S0(2) - shy
      endif
      if (present(shz)) then
         ! Shift along z direction:
         Scell(NSC)%MDatoms(k)%S(3) = Scell(NSC)%MDatoms(k)%S(3) - shz
         Scell(NSC)%MDatoms(k)%S0(3) = Scell(NSC)%MDatoms(k)%S0(3) - shz
      endif
   enddo
   call check_periodic_boundaries(matter, Scell, NSC) ! put atoms back into the supercell
   call Coordinates_rel_to_abs(Scell, NSC)
end subroutine Shift_all_atoms



subroutine rotate_sample(matter, Scell, NSC, Nat_start, Nat_end, theta, phi)
   type(solid), intent(in) :: matter     ! materil parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   integer, intent(in), optional :: Nat_start, Nat_end ! atoms from start to end ones to rotate
   real(8), intent(in), optional :: theta ! angle counted within (X,Y)-plane from Y-axis clockwise (0:2*Pi)
   real(8), intent(in), optional :: phi   ! angle counted from Z axis (-Pi:Pi)
   ! Spherical coordinates are defined as follows:
   ! X = R*sin(phi)*sin(theta)
   ! Y = R*sin(phi)*cos(theta)
   ! Z = R*cos(phi)
   !------------------------------------------
   integer k, M0, Me
   real(8) teta, fi, R, t, p
   real(8), pointer :: X, Y, Z
   
   if (present(Nat_start)) then
      M0 = Nat_start ! starting from this atom
   else
      M0 = 1 ! start from the first one by default
   endif
   if (present(Nat_end)) then
      Me = Nat_end ! end at this atom
   else
      Me = Scell(NSC)%Na ! all atoms by default
   endif
   if (present(theta)) then ! rotate atoms by this theta angle
      teta = theta*g_Pi/180.0d0 ! rotate by this angle
   else
      teta = 0.0d0 ! no rotation in this angular direction
   endif
   if (present(phi)) then ! rotate atoms by this phi angle
      fi = phi*g_Pi/180.0d0 ! rotate by this angle
   else
      fi = 0.0d0 ! no rotation in this angular direction
   endif

   do k = M0, Me ! for those atoms
      X => Scell(NSC)%MDatoms(k)%S(1)
      Y => Scell(NSC)%MDatoms(k)%S(2)
      Z => Scell(NSC)%MDatoms(k)%S(3)
      R = SQRT(X*X + Y*Y + Z*Z)
      t = DACOS(Z/R) + teta ! rotate by this angle
      p = DATAN(X/Y) + fi   ! rotate by this angle
      X = R*dsin(p)*dsin(t)
      Y = R*dsin(p)*dcos(t)
      Z = R*dcos(p)
   enddo
   call Coordinates_rel_to_abs(Scell, NSC)
   nullify(Z,Y,X)
end subroutine rotate_sample


subroutine C60_crystal_construction(Scell, matter, N_start, N_end)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   integer, intent(in), optional :: N_start, N_end
   !------------------------------------------
   integer :: i
   
   call Shift_all_atoms(matter, Scell, 1, shx=0.5d0, shy=0.5d0, N_start=1, N_end=60) ! above
!    call Shift_all_atoms(matter, Scell, 1, shz=0.0252d0, N_start=61, N_end=120) ! above
   !call rotate_sample(matter, Scell, 1, Nat_start=61, Nat_end=120, theta=22.0d0, phi=22.0d0) ! above
   
    open(6543, file='C60_Pa3.dat')
    do i = 1, size(Scell(1)%MDatoms)
       write(6543, '(i8,e24.16,e24.16,e24.16)') Scell(1)%MDatoms(i)%KOA, Scell(1)%MDatoms(i)%S(:)
    enddo
    close(6543)
    pause 'C60_crystal_construction completed'
end subroutine C60_crystal_construction



subroutine Rescale_atomic_velocities(dE_nonadiabat, matter, Scell, NSC, nrg)
   real(8), intent(in) :: dE_nonadiabat     ! energy gain by atoms [eV]
   type(solid), intent(in) :: matter    ! material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Energies), intent(inout) :: nrg ! all energies
   !------------------------------
   real(8) :: Ekin_tot, alpha, b, eps, Ekin, Epot
   integer :: i, Nat
   
   eps = 1.0d-13   ! acceptable change of energy
   
   Nat = Scell(NSC)%Na ! number of atoms
   call Atomic_kinetic_energies(Scell, NSC, matter)
   Ekin_tot = SUM(Scell(NSC)%MDatoms(:)%Ekin)   ! total atomic kinetic energy
   Ekin_tot = Ekin_tot + Supce_kin_energy(Scell, NSC, matter%W_PR)           ! supercell kinetic energy, function below

   if ((ABS(Ekin_tot) > eps) .and. (ABS(dE_nonadiabat) > Ekin_tot * eps)) then ! makes sense to change it
      ! Test if it worked correctly:
!       print*, 'Rescale_atomic_velocities before', Ekin_tot, dE_nonadiabat, Ekin_tot + dE_nonadiabat

      ! Get the scaling coefficient:
      b = dE_nonadiabat/Ekin_tot   ! ration entering the scaling
      if (b < -1.0d0) then
         alpha = 0.0d0
      else
         alpha = sqrt(1.0d0 + b)
      endif

      ! Rescale atomic velocities:
      do i = 1, Nat
         Scell(NSC)%MDatoms(i)%V(:) = Scell(NSC)%MDatoms(i)%V(:) * alpha
      enddo

      ! Rescale supercell velocities:
      Scell(NSC)%Vsupce = Scell(NSC)%Vsupce * alpha

      ! Rescaling the relative velocities:
      call velocities_abs_to_rel(Scell, NSC) ! !New relative velocities
      ! And energies:
      call Atomic_kinetic_energies(Scell, NSC, matter)

      ! And supercell energy:
      Ekin = Supce_kin_energy(Scell, NSC, matter%W_PR)    ! function below
      Epot = (matter%p_ext*Scell(NSC)%V)*1d-30/g_e	! potential part of the energy of the supercell [eV]
      nrg%E_supce = (Ekin + Epot)/dble(Nat) 	! total energy of the supercell [eV/atom]
   
      ! Test if it worked correctly:
!       print*, 'Rescale_atomic_velocities after', SUM(Scell(NSC)%MDatoms(:)%Ekin), Ekin_tot*(alpha*alpha), alpha
   endif ! (abs(dE_nonadiabat) > Ekin_tot * eps)
end subroutine Rescale_atomic_velocities



subroutine Atomic_kinetic_energies(Scell, NSC, matter)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in), target :: matter
   integer i, N
   real(8), pointer :: Mass
   N = size(Scell(NSC)%MDatoms)
   do i = 1,N
      !atoms(i)%Ekin = 0.5d0*matter%Ma*SUM(atoms(i)%V(:)*atoms(i)%V(:))*1d10/g_e	! [eV]
      Mass => matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma
      !Scell(NSC)%MDatoms(i)%Ekin = 0.5d0*matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma*SUM(Scell(NSC)%MDatoms(i)%V(:)*Scell(NSC)%MDatoms(i)%V(:))*1d10/g_e	! [eV]
      Scell(NSC)%MDatoms(i)%Ekin = 0.5d0*Mass*SUM(Scell(NSC)%MDatoms(i)%V(:)*Scell(NSC)%MDatoms(i)%V(:))*1d10/g_e	! [eV]
   enddo
   nullify(Mass)
end subroutine Atomic_kinetic_energies



subroutine Rescale_atomic_velocities_OLD(dE_nonadiabat, matter, Scell, NSC, nrg)
   REAL(8), INTENT(in) :: dE_nonadiabat 	! energy gain by atoms
   type(solid), intent(in), target :: matter	! material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of supercell
   type(Energies), intent(inout) :: nrg	! all energies
   real(8), pointer :: Mass
   integer i, Nat, ik
   real(8) Etot_at, dE, V, sign_f, e, vi(3), E_at_temp, E_temp, V0(3), koef(3), E_proj, eps
   parameter (e = g_e)         ! Electron charge       [Coulomb]
   parameter (eps = 1.0d-12)    ! presision in energy [eV]
   
   Nat = Scell(NSC)%Na ! number of atoms
   call Atomic_kinetic_energies(Scell, NSC, matter)
   !Etot_at = sqrt(nrg%At_kin*real(matter%Na))
   Etot_at = SUM(sqrt(Scell(NSC)%MDatoms(:)%Ekin))

   ! Do only if there is anything to do:
   if ((ABS(Etot_at) > eps) .and. (ABS(dE_nonadiabat) > eps)) then
      ! Rescale velocities according to changed energy:
      E_at_temp = 0.0d0
      E_temp = 0.0d0
      do i = 1, Nat
         !dE = dE_nonadiabat*(atoms(i)%Ekin)/(Etot_at) ! fraction of energy given to this atom
         dE = dE_nonadiabat*sqrt(Scell(NSC)%MDatoms(i)%Ekin)/(Etot_at) ! fraction of energy given to this atom
         !V = ABS(atoms(i)%V(1)) + ABS(atoms(i)%V(2)) + ABS(atoms(i)%V(3)) ! to get partions - how to distribute energy
         !koef(:) = ABS(atoms(i)%V(:))/V
         V = Scell(NSC)%MDatoms(i)%V(1)*Scell(NSC)%MDatoms(i)%V(1) + Scell(NSC)%MDatoms(i)%V(2)*Scell(NSC)%MDatoms(i)%V(2) + Scell(NSC)%MDatoms(i)%V(3)*Scell(NSC)%MDatoms(i)%V(3)
         koef(:) = Scell(NSC)%MDatoms(i)%V(:)*Scell(NSC)%MDatoms(i)%V(:)/V
         V0(:) = 0.0d0
         Mass => matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma
         do ik = 1,3
            if (Scell(NSC)%MDatoms(i)%V(ik) .GT. 0.0d0) then
               sign_f = 1.0d0
            else
               sign_f = -1.0d0
            endif
            !E_proj = matter%Ma*(atoms(i)%V(ik)*atoms(i)%V(ik))*0.5d0/e*1d10
            !E_proj = matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma*(Scell(NSC)%MDatoms(i)%V(ik)*Scell(NSC)%MDatoms(i)%V(ik))*0.5d0/e*1d10
            E_proj = Mass*(Scell(NSC)%MDatoms(i)%V(ik)*Scell(NSC)%MDatoms(i)%V(ik))*0.5d0/e*1d10
            !V0(ik) = sign_f*SQRT(2.0d0/matter%Ma*(ABS(E_proj + koef(ik)*dE)*e))*1d-5
            !V0(ik) = sign_f*SQRT(2.0d0/matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma*(ABS(E_proj + koef(ik)*dE)*e))*1d-5
            V0(ik) = sign_f*SQRT(2.0d0/Mass*(ABS(E_proj + koef(ik)*dE)*e))*1d-5
         enddo ! ik

         do ik = 1,3	! check
            if (ABS((V0(ik)-Scell(NSC)%MDatoms(i)%V(ik))/Scell(NSC)%MDatoms(i)%V(ik)) .GT. 10.0d0) then ! too strong change, SMTH might be wrong...
               write(*,'(a,i12,i12,e25.16,e25.16,e25.16,e25.16,e25.16)') 'SPEED CHANGE:', i, ik, V0(ik), Scell(NSC)%MDatoms(i)%V(ik), E_proj, koef(ik), dE
               !print*, SUM(V0(:)*V0(:))*matter%Ma*0.5d0/e*1d10, atoms(i)%Ekin+dE, dE, dE_nonadiabat/matter%Na
               print*, SUM(V0(:)*V0(:))*Mass*0.5d0/e*1d10, Scell(NSC)%MDatoms(i)%Ekin+dE, dE, dE_nonadiabat/Scell(NSC)%Na
               print*, koef(:)
               print*, V0(:)
               print*, Scell(NSC)%MDatoms(i)%V(:)
            endif
         enddo
         Scell(NSC)%MDatoms(i)%V(:) = V0(:)
         !E_at_temp = E_at_temp + matter%Ma*((atoms(i)%V(1)*atoms(i)%V(1) + atoms(i)%V(2)*atoms(i)%V(2) + atoms(i)%V(3)*atoms(i)%V(3)))*0.5d0/e*1d10
         E_at_temp = E_at_temp + Mass*((Scell(NSC)%MDatoms(i)%V(1)*Scell(NSC)%MDatoms(i)%V(1) + Scell(NSC)%MDatoms(i)%V(2)*Scell(NSC)%MDatoms(i)%V(2) + Scell(NSC)%MDatoms(i)%V(3)*Scell(NSC)%MDatoms(i)%V(3)))*0.5d0/e*1d10
         E_temp = E_temp + dE*SUM(koef(:))
      enddo ! i
      call velocities_abs_to_rel(Scell, NSC) ! !New relative velocities
      call Atomic_kinetic_energies(Scell, NSC, matter)
   endif ! ((ABS(Etot_at) > eps) .and. (ABS(dE_nonadiabat) > eps))
   nullify(Mass)
end subroutine Rescale_atomic_velocities_OLD



function Maxwell_int(Ta, E)
   real(8), intent(in) :: Ta	! temperature [eV]
   real(8), intent(inout) :: E	! lower limit of integration, upper one is infinity [eV]
   real(8) :: Maxwell_int   ! normalized to 1
!    if ((E .LT. 0.0d0) .AND. (E .GT. -1.0d-14)) E = 0.0d0
   if (E < 0.0d0) E = 0.0d0
   if (E/Ta < huge(E)) then
      Maxwell_int = 2.0d0*dsqrt(E/(g_Pi*Ta))*dexp(-E/Ta) - (derf(dsqrt(E/Ta)) - 1.0d0)
   else
      Maxwell_int = 0.0d0
   endif
   if (isnan(Maxwell_int)) then
      print*, E, Ta, Maxwell_int
      pause 'Maxwell is NaN'
   endif
end function Maxwell_int



function Maxwell_int_shifted(Ta, hw) result(G)
   real(8), intent(in) :: Ta	! temperature [eV]
   real(8), intent(inout) :: hw	! [eV] shift og the Maxwell function
   real(8) :: G ! normalized to 1
   real(8) :: eps   ! tolerance around zero
   real(8) :: arg
   eps = 1.0d-10
   if (abs(hw) < eps) then
      G = 1.0d0 ! unshifted Maxwell integrated from 0 to infinity gives 1
   else ! shifted Maxwell gives exponent
      arg = hw/Ta
      if (arg > huge(eps)) then
         G = 0.0d0
      else
          G = exp(-hw/Ta)
      endif
   endif
end function Maxwell_int_shifted





! This subroutine calculates an input of the super-cell changing into the atomic motion
! (it's not directly related to Verlet, but to the Parrinello-Rahman method):
subroutine PR_sc_at_OLD(supce, Vsupce, SVco, k, ggs)
   REAL(8), DIMENSION(:,:), INTENT(in) ::  Supce    ! P-R super cell (3x3 matrix)
   REAL(8), DIMENSION(:,:), INTENT(in) ::  Vsupce   ! P-R speed of super cell changes (3x3 matrix)
   REAL(8), DIMENSION(:), INTENT(in) ::  SVco      ! Relative atomic coordinates
   integer, INTENT(in) :: k  ! number of atom
   REAL(8), DIMENSION(3), INTENT(out) :: ggs      ! part of acceleration, related to P-R super-cell changes
   integer i,j
   REAL(8), DIMENSION(3,3) :: TSupce, TVsupce, Ginv, Gdot, Gdot1, Gdot2, GG, GinvGdot
   REAL(8), DIMENSION(3) :: v0
   ! Get g^(-1) term:
   call Transpose_M(Supce,TSupce)      ! transpose super-cell matrix: h^T
   call Two_Matr_mult(TSupce,Supce,GG) ! construct G-matrix: g=h^T*h
   call Invers_3x3(GG, Ginv, 'PR_sc_at_OLD')           ! inverse G-matrix: g^(-1)
   ! Get dg/dt term:
   call Transpose_M(VSupce,TVSupce)           ! transpose super-cell velocity matrix: hdot^(T)
   call Two_Matr_mult(TVsupce, Supce, Gdot1)  ! first half of dg/dt: hdot^(T)*h
   call Two_Matr_mult(TSupce, Vsupce, Gdot2)  ! second half of dg/dt: h^(T)*hdot
   Gdot = Gdot2 + Gdot1 	! dg/dt = hdot^(T)*h + h^(T)*hdot

   ! Get g^(-1)*dg/dt*ds/dt term:
   call Two_Matr_mult(Ginv, Gdot, GinvGdot)  ! a factor before sdot
   call Matrix_Vec_Prod(GinvGdot, SVco, ggs) ! full g^(-1)*dg/dt*ds/dt
end subroutine PR_sc_at_OLD


subroutine get_GGdot(supce, Vsupce, GinvGdot)
   REAL(8), DIMENSION(:,:), INTENT(in) ::  Supce    ! P-R super cell (3x3 matrix)
   REAL(8), DIMENSION(:,:), INTENT(in) ::  Vsupce   ! P-R speed of super cell changes (3x3 matrix)
   REAL(8), DIMENSION(3,3), INTENT(out) :: GinvGdot ! part of acceleration, related to P-R super-cell changes
   REAL(8), DIMENSION(3,3) :: TSupce, TVsupce, Ginv, Gdot, Gdot1, Gdot2, GG

   ! Get g^(-1) term:
   call Transpose_M(Supce,TSupce)      ! transpose super-cell matrix: h^T
   call Two_Matr_mult(TSupce,Supce,GG) ! construct G-matrix: g=h^T*h
   call Invers_3x3(GG, Ginv, 'get_GGdot')           ! inverse G-matrix: g^(-1)
   ! Get dg/dt term:
   call Transpose_M(VSupce,TVSupce)           ! transpose super-cell velocity matrix: hdot^(T)
   call Two_Matr_mult(TVsupce, Supce, Gdot1)  ! first half of dg/dt: hdot^(T)*h
   call Two_Matr_mult(TSupce, Vsupce, Gdot2)  ! second half of dg/dt: h^(T)*hdot
   Gdot = Gdot2 + Gdot1 	! dg/dt = hdot^(T)*h + h^(T)*hdot

   ! Get g^(-1)*dg/dt term:
   call Two_Matr_mult(Ginv, Gdot, GinvGdot)  ! a factor before sdot
end subroutine get_GGdot


subroutine PR_sc_at(GinvGdot, SVco, ggs)
   REAL(8), DIMENSION(3,3), INTENT(in) :: GinvGdot      ! part of acceleration, related to
   REAL(8), DIMENSION(3), INTENT(in) ::  SVco      ! Relative atomic coordinates
   REAL(8), DIMENSION(3), INTENT(out) :: ggs      ! part of acceleration, related to P-R super-cell changes
   call Matrix_Vec_Prod(GinvGdot, SVco, ggs) ! full g^(-1)*dg/dt*ds/dt
end subroutine PR_sc_at



subroutine make_time_step_atoms(Scell, matter, numpar, ind)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer, intent(in) :: ind	! =1, or =2, first or second half of the velocity Verlet algorithm
   !=========================
   integer :: NSC ! number of super-cell

   do NSC = 1, size(Scell)
      ! Make MD step:
      call make_time_step_atoms_SC(Scell, NSC, matter, numpar, ind) ! see below
      ! Update absolute coordinates:
      call Coordinates_rel_to_abs(Scell, NSC)
      ! Update absolute velocities:
      call velocities_abs_to_rel(Scell, NSC) ! !New relative velocities
   enddo
end subroutine make_time_step_atoms


subroutine make_time_step_atoms_SC(Scell, NSC, matter, numpar, ind)     ! update coordinates and velocities of atoms via intermediate Verlet step
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in), target :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer, intent(in) :: ind	! =1, or =2, first or second half of the velocity Verlet algorithm
   !----------------------------------
   real(8), pointer :: Mass
   REAL(8), DIMENSION(3,3) :: GinvGdot ! matrix of GGdot contributing to velocities
   integer :: nat, k
   real(8) Fors_s(3,Scell(NSC)%Na), fctr
   real(8) dsupce(3,3), tempM(3,3), gg2(3,3), tt(3,3), tempV(3), x0(3)

   nat = Scell(NSC)%Na	! number of atoms in the supercell
   
   call get_GGdot(Scell(NSC)%supce, Scell(NSC)%Vsupce, GinvGdot) ! see above
   
   Fors_s = 0.0d0
   do k = 1,nat ! All atoms - calculating new coordinates:
      call Transpose_M(Scell(NSC)%supce,dsupce) ! transpose super-cell matrix, from "Algebra_tools" module
      call Two_Matr_mult(dsupce,Scell(NSC)%supce,tempM)  ! g from Eq.(2.15), H.Jeschke PHD thesis, p.36, module "Algebra_tools"
!       call Two_Matr_mult(Scell(NSC)%supce,dsupce,tempM)  ! g from Eq.(2.15), H.Jeschke PHD thesis, p.36, module "Algebra_tools" 
      call Invers_3x3(tempM, gg2, 'make_time_step_atoms_SC') ! to get g^(-1), module "Algebra_tools"
      tempV(:) = Scell(NSC)%MDatoms(k)%forces%total(:)
      call Matrix_Vec_Prod(gg2,tempV,x0) ! module "Algebra_tools"

      Mass => matter%Atoms(Scell(NSC)%MDatoms(k)%KOA)%Ma ! atomic mass
      fctr = -1.0d0/Mass*g_e*1d-10 ! to transfer forces into the proper units: -> [eV/A] -> acceleration
      Fors_s(:,k) = fctr * x0(:)	 ! make proper normalization

      !call PR_sc_at_OLD(Scell(NSC)%supce, Scell(NSC)%Vsupce, Scell(NSC)%MDatoms(k)%SV0, k, tempV) ! second term in acceleration coming from the super-cell
      call PR_sc_at(GinvGdot, Scell(NSC)%MDatoms(k)%SV0, tempV) ! see above

      Fors_s(:,k) = Fors_s(:,k) - tempV(:)
      ! Verlet step of coordinates:
      if (ind .EQ. 2) Scell(NSC)%MDatoms(k)%S(:) = Scell(NSC)%MDatoms(k)%S0(:) + numpar%dt*Scell(NSC)%MDatoms(k)%SV0(:) + Fors_s(:,k)*numpar%dtsqare ! new X-coordinates
      ! Verlet velocities second part:
      Scell(NSC)%MDatoms(k)%SV(:) = Scell(NSC)%MDatoms(k)%SV0(:) + Fors_s(:,k)*numpar%halfdt !dt/2.0e0
   enddo
   if (ind .EQ. 2) call check_periodic_boundaries(matter, Scell, NSC) ! and set the absolute coordinates out of the new relative ones
   call velocities_rel_to_abs(Scell, NSC) ! set the absolute velocities out of the new relative ones
   
   nullify(Mass)
end subroutine make_time_step_atoms_SC



subroutine make_time_step_atoms_M(Scell, matter, numpar, ind)   ! Martyna algorithm
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer, intent(in) :: ind	! step of Matryna algorithm
   !=========================
   integer :: NSC   ! number of super-cell
   integer :: k     ! atoms index

   do NSC = 1, size(Scell)
      ! Make MD step:
      select case (ind)
      case (1)  ! coordinate
         !$omp PARALLEL do private(k)
         do k = 1,Scell(NSC)%Na ! All atoms - calculating new coordinates:
            ! Martyna step of coordinates:
            Scell(NSC)%MDatoms(k)%S(:) = Scell(NSC)%MDatoms(k)%S0(:) + numpar%dt*Scell(NSC)%MDatoms(k)%SV0(:) + &
               numpar%dtsqare * 0.25d0*(5.0d0*Scell(NSC)%MDatoms(k)%A(:) - Scell(NSC)%MDatoms(k)%A_tild(:)) + &
               numpar%dt3 * Scell(NSC)%MDatoms(k)%v_F(:) + numpar%dt4 * Scell(NSC)%MDatoms(k)%v_J(:)
            ! Update old coordinates:
            Scell(NSC)%MDatoms(k)%S0(:) = Scell(NSC)%MDatoms(k)%S(:)
         enddo
         !$omp end parallel do

         call check_periodic_boundaries(matter, Scell, NSC) ! and set the absolute coordinates out of the new relative ones
         
         ! Update absolute coordinates:
         call Coordinates_rel_to_abs(Scell, NSC)
      
      case (2)  ! velocity
         ! Update accelerations from the new potential:
         call get_accelerations_M(Scell, NSC, matter, numpar)   ! below

         ! Update relative velosities:
         !$omp PARALLEL do private(k)
         do k = 1,Scell(NSC)%Na ! fro all atoms
            ! Make a step for relative velocities:
            Scell(NSC)%MDatoms(k)%SV(:) = Scell(NSC)%MDatoms(k)%SV0(:) + &
               (numpar%dt/8.0d0)*(4.0d0*Scell(NSC)%MDatoms(k)%A0(:) + Scell(NSC)%MDatoms(k)%A_tild(:) + 3.0d0*Scell(NSC)%MDatoms(k)%A(:)) +&
               (numpar%dtsqare/4.0d0)*Scell(NSC)%MDatoms(k)%v_F(:)
            ! Update old relative velocities:
            Scell(NSC)%MDatoms(k)%SV0(:) = Scell(NSC)%MDatoms(k)%SV(:)
         enddo
         !$omp end parallel do
         ! Update absolute velosities:
         call velocities_rel_to_abs(Scell, NSC) ! set the absolute velocities out of the new relative ones
         ! Update old absolute velocities:
         !$omp PARALLEL do private(k)
         do k = 1,Scell(NSC)%Na ! fro all atoms
            Scell(NSC)%MDatoms(k)%V0(:) = Scell(NSC)%MDatoms(k)%V(:)
         enddo
         !$omp end parallel do
      case (3)  ! effective force
         !$omp PARALLEL do private(k)
         do k = 1,Scell(NSC)%Na ! fro all atoms
            ! MAke a step for effective forces:
            Scell(NSC)%MDatoms(k)%A_tild(:) = 0.5d0*(2.0d0*Scell(NSC)%MDatoms(k)%A0(:) +&
                Scell(NSC)%MDatoms(k)%A(:) - Scell(NSC)%MDatoms(k)%A_tild0(:)) + numpar%dt*0.5d0*Scell(NSC)%MDatoms(k)%v_F(:)
            ! Update old accelerations:
            Scell(NSC)%MDatoms(k)%A0(:) = Scell(NSC)%MDatoms(k)%A(:)
         enddo
         !$omp end parallel do
      case (4)  ! effective force velocity 
         !$omp PARALLEL do private(k)
         do k = 1,Scell(NSC)%Na ! fro all atoms
            ! Make a step for effective force velocities:
            Scell(NSC)%MDatoms(k)%v_F(:) = -0.50d0 * Scell(NSC)%MDatoms(k)%v_F0(:) + &
               1.50d0/numpar%dt * (Scell(NSC)%MDatoms(k)%A(:) - Scell(NSC)%MDatoms(k)%A_tild0(:))
         enddo
         !$omp end parallel do
      case (5)  ! effective force acceleration
         !$omp PARALLEL do private(k)
         do k = 1,Scell(NSC)%Na ! fro all atoms
            ! Make a step for effective force accelerations:
            Scell(NSC)%MDatoms(k)%v_J(:) = -Scell(NSC)%MDatoms(k)%v_J0(:) - 3.0d0/numpar%dt * Scell(NSC)%MDatoms(k)%v_F0(:) +&
                1.5d0/(numpar%dt*numpar%dt) * (Scell(NSC)%MDatoms(k)%A(:) - Scell(NSC)%MDatoms(k)%A_tild0(:))
            ! Update old effective force:
            Scell(NSC)%MDatoms(k)%A_tild0(:) = Scell(NSC)%MDatoms(k)%A_tild(:)
            ! Update old effective force velocity and acceleration:
            Scell(NSC)%MDatoms(k)%v_F0(:) = Scell(NSC)%MDatoms(k)%v_F(:)
            Scell(NSC)%MDatoms(k)%v_J0(:) = Scell(NSC)%MDatoms(k)%v_J(:)
         enddo
         !$omp end parallel do
      endselect
   enddo
end subroutine make_time_step_atoms_M



subroutine get_accelerations_M(Scell, NSC, matter, numpar)     ! update coordinates and velocities of atoms via intermediate Verlet step
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in), target :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !----------------------------------
   real(8), pointer :: Mass
   REAL(8), DIMENSION(3,3) :: GinvGdot ! matrix of GGdot contributing to velocities
   integer :: k
   real(8) :: Fors_s(3), fctr
   real(8) :: dsupce(3,3), tempM(3,3), gg2(3,3), tt(3,3), tempV(3), x0(3)

   call get_GGdot(Scell(NSC)%supce, Scell(NSC)%Vsupce, GinvGdot) ! see above
   
   !Fors_s(:) = 0.0d0
!    !$omp parallel private(k, Fors_s, dsupce, tempM, gg2, tempV, x0, Mass, fctr)
!    !$omp do
   do k = 1,Scell(NSC)%Na ! All atoms - calculating new coordinates:
      call Transpose_M(Scell(NSC)%supce,dsupce) ! transpose super-cell matrix, from "Algebra_tools" module
      call Two_Matr_mult(dsupce,Scell(NSC)%supce,tempM)  ! g from Eq.(2.15), H.Jeschke PHD thesis, p.36, module "Algebra_tools"
      call Invers_3x3(tempM, gg2, 'get_accelerations_M') ! to get g^(-1), module "Algebra_tools"
      tempV(:) = Scell(NSC)%MDatoms(k)%forces%total(:)
      call Matrix_Vec_Prod(gg2,tempV,x0) ! module "Algebra_tools"

      Mass => matter%Atoms(Scell(NSC)%MDatoms(k)%KOA)%Ma ! atomic mass
      fctr = -1.0d0/Mass*g_e*1d-10 ! to transfer forces into the proper units: -> [eV/A] -> acceleration
      Fors_s(:) = fctr * x0(:)	 ! make proper normalization
      ! second term in acceleration coming from the super-cell
      call PR_sc_at(GinvGdot, Scell(NSC)%MDatoms(k)%SV0, tempV) ! see above
      ! Total acceleration:
      Fors_s(:) = Fors_s(:) - tempV(:)
      
      ! Get new accelerations:
      Scell(NSC)%MDatoms(k)%A(:) = Fors_s(:)
   enddo
!    !$omp end do
!    !$omp end parallel
   
   nullify(Mass)
end subroutine get_accelerations_M




subroutine make_time_step_atoms_Y4(Scell, matter, numpar, ind_step, ind_cv) ! Yoshida MD algorithm, 4th order
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   integer, intent(in) :: ind_step  ! steps of Yoshida algorithm 1 to 4
   integer, intent(in) :: ind_cv    ! is this steps of Yoshida algorithm for coordinates or velocities
   !=========================
   integer :: NSC ! number of super-cell

   do NSC = 1, size(Scell)
      ! Make MD step:
      call make_time_step_atoms_SC_Y4(Scell, NSC, matter, numpar, ind_step, ind_cv) ! see below
      ! Update absolute coordinates:
      call Coordinates_rel_to_abs(Scell, NSC)    ! below
      ! Update absolute velocities:
      call velocities_abs_to_rel(Scell, NSC) ! below
   enddo
end subroutine make_time_step_atoms_Y4


subroutine make_time_step_atoms_SC_Y4(Scell, NSC, matter, numpar, ind_step, ind_cv)     ! update coordinates and velocities of atoms
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in), target :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   integer, intent(in) :: ind_step  ! steps of Yoshida algorithm 1 to 4
   integer, intent(in) :: ind_cv    ! is this steps of Yoshida algorithm for coordinates or velocities
   !----------------------------------
   real(8), pointer :: Mass
   REAL(8), DIMENSION(3,3) :: GinvGdot ! matrix of GGdot contributing to velocities
   integer :: nat, k
   real(8) Fors_s(3,Scell(NSC)%Na), fctr
   real(8) dsupce(3,3), tempM(3,3), gg2(3,3), tt(3,3), tempV(3), x0(3)
   real(8) :: Coef_C, Coef_V

   nat = Scell(NSC)%Na	! number of atoms in the supercell
   call get_GGdot(Scell(NSC)%supce, Scell(NSC)%Vsupce, GinvGdot) ! see above
   
   ! Get Yoshida coeefs depending on the step number:
   select case(ind_step)
   case (1)
      Coef_C = m_c1 
      Coef_V = m_d1
   case (2)
      Coef_C = m_c2
      Coef_V = m_d2
   case (3)
      Coef_C = m_c3
      Coef_V = m_d3
   case (4)
      Coef_C = m_c4
      Coef_V = 0.0d0
   endselect
   
   Fors_s = 0.0d0
   do k = 1,nat ! All atoms - calculating new coordinates:
      if (ind_cv == 1) then ! Yoshida step for coordinates:
         Scell(NSC)%MDatoms(k)%S(:) = Scell(NSC)%MDatoms(k)%S0(:) + Scell(NSC)%MDatoms(k)%SV0(:) * numpar%dt * Coef_C
         Scell(NSC)%MDatoms(k)%S0(:) = Scell(NSC)%MDatoms(k)%S(:)   ! update old coords for the next step
      else ! Yoshida step for velocities:
         call Transpose_M(Scell(NSC)%supce,dsupce) ! transpose super-cell matrix, from "Algebra_tools" module
         call Two_Matr_mult(dsupce,Scell(NSC)%supce,tempM)  ! g from Eq.(2.15), H.Jeschke PHD thesis, p.36, module "Algebra_tools"
         call Invers_3x3(tempM, gg2, 'make_time_step_atoms_SC') ! to get g^(-1), module "Algebra_tools"
         tempV(:) = Scell(NSC)%MDatoms(k)%forces%total(:)
         call Matrix_Vec_Prod(gg2,tempV,x0) ! module "Algebra_tools"
         Mass => matter%Atoms(Scell(NSC)%MDatoms(k)%KOA)%Ma ! atomic mass
         fctr = -1.0d0/Mass*g_e*1d-10 ! to transfer forces into the proper units: -> [eV/A] -> acceleration
         Fors_s(:,k) = fctr * x0(:)	 ! make proper normalization
         call PR_sc_at(GinvGdot, Scell(NSC)%MDatoms(k)%SV0, tempV) ! see above
         Fors_s(:,k) = Fors_s(:,k) - tempV(:)
         Scell(NSC)%MDatoms(k)%SV(:) = Scell(NSC)%MDatoms(k)%SV0(:) + Fors_s(:,k) * numpar%dt * Coef_V
         Scell(NSC)%MDatoms(k)%SV0(:) = Scell(NSC)%MDatoms(k)%SV(:)   ! update old coords for the next step
      endif
   enddo
   call check_periodic_boundaries(matter, Scell, NSC) ! and set the absolute coordinates out of the new relative ones
   call velocities_rel_to_abs(Scell, NSC) ! set the absolute velocities out of the new relative ones
   
   nullify(Mass)
end subroutine make_time_step_atoms_SC_Y4


subroutine make_time_step_supercell_Y4(Scell, matter, numpar, ind_step, ind_cv)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer, intent(in) :: ind_step  ! steps of Yoshida algorithm 1 to 4
   integer, intent(in) :: ind_cv    ! is this steps of Yoshida algorithm for coordinates or velocities
   !=========================
   integer :: NSC ! number of super-cell

   do NSC = 1, size(Scell)
      ! Make a Verlet step for the super cell:
      call make_time_step_supercell_SC_Y4(Scell, NSC, matter, numpar, Scell(NSC)%SCforce, ind_step, ind_cv) ! see below
      ! Now update kinetic energies of atoms and supercell:
      call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg)
   enddo
end subroutine make_time_step_supercell_Y4


subroutine make_time_step_supercell_SC_Y4(Scell, NSC, matter, numpar, supce_forces, ind_step, ind_cv)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Supce_force), intent(inout) :: supce_forces
   integer, intent(in) :: ind_step  ! steps of Yoshida algorithm 1 to 4
   integer, intent(in) :: ind_cv    ! is this steps of Yoshida algorithm for coordinates or velocities
   integer i, j
   real(8) :: Coef_C, Coef_V

   if (numpar%p_const) then ! P = const, else V=const   
      ! Get Yoshida coeefs depending on the step number:
      select case(ind_step)
      case (1)
         Coef_C = m_c1 
         Coef_V = m_d1
      case (2)
         Coef_C = m_c2
         Coef_V = m_d2
      case (3)
         Coef_C = m_c3
         Coef_V = m_d3
      case (4)
         Coef_C = m_c4
         Coef_V = 0.0d0
      endselect
      do i = 1,3
         do j = 1,3
            if (ind_cv == 1) then ! Yoshida step for coordinates:
               Scell(NSC)%supce(i,j) = Scell(NSC)%supce0(i,j) + Scell(NSC)%Vsupce0(i,j) * numpar%dt * Coef_C
               Scell(NSC)%supce0(i,j) = Scell(NSC)%supce(i,j)   ! save for next step
            else ! Yoshida step for velocities:
               Scell(NSC)%Vsupce(i,j) = Scell(NSC)%Vsupce0(i,j) + supce_forces%total(i,j) * numpar%dt * Coef_V
               Scell(NSC)%Vsupce0(i,j) = Scell(NSC)%Vsupce(i,j)   ! save for next step
            endif
         enddo ! j
      enddo ! i
      if (ind_cv .EQ. 1) then ! Update inverse and GG:
         call Transpose_M(Scell(NSC)%supce,Scell(NSC)%supce_t) ! transpose super-cell matrix: h^T
         call Two_Matr_mult(Scell(NSC)%supce_t,Scell(NSC)%supce,Scell(NSC)%GG) ! construct G-matrix: g=h^T*h
      endif
   endif
end subroutine make_time_step_supercell_SC_Y4




subroutine make_time_step_supercell_M(Scell, matter, numpar, ind)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer, intent(in) :: ind	! =1, or =2, first or second half of the velocity Verlet algorithm
   !=========================
   integer :: NSC ! number of super-cell

   do NSC = 1, size(Scell)
      ! Make a Verlet step for the super cell:
      call make_time_step_supercell_SC(Scell, NSC, matter, numpar, Scell(NSC)%SCforce, ind) ! see below
      ! Now update kinetic energies of atoms and supercell:
      call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg)
   enddo
end subroutine make_time_step_supercell_M



subroutine make_time_step_supercell(Scell, matter, numpar, ind)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   integer, intent(in) :: ind	! =1, or =2, first or second half of the velocity Verlet algorithm
   !=========================
   integer :: NSC ! number of super-cell

   do NSC = 1, size(Scell)
      ! Make a Verlet step for the super cell:
      call make_time_step_supercell_SC(Scell, NSC, matter, numpar, Scell(NSC)%SCforce, ind) ! see below
      ! Now update kinetic energies of atoms and supercell:
      call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg)
   enddo
end subroutine make_time_step_supercell 



subroutine make_time_step_supercell_SC(Scell, NSC, matter, numpar, supce_forces, ind)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in) :: matter	! material parameters
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Supce_force), intent(inout) :: supce_forces
   integer, intent(in) :: ind	! =1, or =2, first or second half of the velocity Verlet algorithm
   integer i, j

   if (numpar%p_const) then ! P = const, else V=const
      do i = 1,3
         do j = 1,3
               if (ind .EQ. 2) then 
                  Scell(NSC)%supce(i,j) = Scell(NSC)%supce0(i,j) + Scell(NSC)%Vsupce0(i,j)*numpar%dt + supce_forces%total(i,j)*numpar%dtsqare !*dt*dt/2.0e0
               endif
               Scell(NSC)%Vsupce(i,j) = Scell(NSC)%Vsupce0(i,j) + supce_forces%total(i,j)*numpar%halfdt !*dt/2.0e0
!                ! FOR TEST ONLY:
!                if (ind .EQ. 2) Scell(NSC)%supce(i,j) = Scell(NSC)%supce0(i,j) + Scell(NSC)%Vsupce0(i,j)*numpar%dt + supce_forces%total(j,i)*numpar%dtsqare !*dt*dt/2.0e0
!                Scell(NSC)%Vsupce(i,j) = Scell(NSC)%Vsupce0(i,j) + supce_forces%total(j,i)*numpar%halfdt !*dt/2.0e0
         enddo ! j
      enddo ! i
      if (ind .EQ. 2) then ! Update inverse and GG:
         call Transpose_M(Scell(NSC)%supce,Scell(NSC)%supce_t) ! transpose super-cell matrix: h^T
         call Two_Matr_mult(Scell(NSC)%supce_t,Scell(NSC)%supce,Scell(NSC)%GG) ! construct G-matrix: g=h^T*h
      endif
   endif
end subroutine make_time_step_supercell_SC


subroutine total_forces(atoms)
   !type(Forces), intent(inout) :: forces	! all interatomic forces
   type(Atom), dimension(:), intent(inout) :: atoms	! array of atoms in the supercell
   integer i, j, N, M
   do i = 1, size(atoms)
      atoms(i)%forces%total(:) = atoms(i)%forces%rep(:) + atoms(i)%forces%att(:) ! total force
      !print*, 'total_forces', i, atoms(i)%forces%rep(:)
   enddo
end subroutine total_forces


subroutine Potential_super_cell_forces(numpar, Scell, NSC, matter)
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors   
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Solid), intent(in) :: matter ! material parameters
   real(8) fact
   integer i, j
   if (numpar%p_const) then	! calculate this for P=const Parrinello-Rahman MD
      Scell(NSC)%SCforce%total = 0.0d0
      fact = g_e/1d10/matter%W_PR
      Scell(NSC)%SCforce%att = Scell(NSC)%SCforce%att*fact !/2.0d0
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep*fact !/2.0d0
      do i = 1,3
         do j = 1,3
!              Scell(NSC)%SCforce%total(i,j) = Scell(NSC)%SCforce%att(i,j) + Scell(NSC)%SCforce%rep(j,i) ! correct
             Scell(NSC)%SCforce%total(j,i) = Scell(NSC)%SCforce%att(i,j) + Scell(NSC)%SCforce%rep(j,i) ! test correct
!              Scell(NSC)%SCforce%total(i,j) = Scell(NSC)%SCforce%att(i,j) + Scell(NSC)%SCforce%rep(i,j) ! not correct
         enddo ! j
      enddo ! i
   endif
end subroutine Potential_super_cell_forces


subroutine super_cell_forces(numpar, Scell, NSC, matter, supce_forces, Sigma_tens_OUT)
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   !type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(solid), intent(in), target :: matter	! material parameters
   type(Supce_force), intent(inout) :: supce_forces
   real(8), dimension(3,3), intent(inout), optional :: Sigma_tens_OUT	! if it is only for pressure calculations, exclude some terms and print out the sigma tensor
   !-----------------------------------
   real(8), pointer :: Mass
   real(8) KPRESint(3,3), a, dsupce(3,3)
   real(8) KPRES_INV(3,3), KPRES_VV(3,3), V, Sigma_tens(3,3)
   integer nat, k, ik, ij, i, j
   if (numpar%p_const) then	! calculate this for P=const Parrinello-Rahman MD
      nat = Scell(NSC)%Na	! number of atoms
      KPRESint = 0.0d0
      KPRES_VV = 0.0d0
      do k = 1,nat ! periodic boundary conditions:
         Mass => matter%Atoms(Scell(NSC)%MDatoms(k)%KOA)%Ma
         do ik = 1,3 
            do ij = 1,3
               ! Parrinello-Rahman version of kinetic term:
               KPRES_VV(ij,ik) = KPRES_VV(ij,ik) + Scell(NSC)%MDatoms(k)%V(ij)*Scell(NSC)%MDatoms(k)%V(ik)*Mass !*
            enddo ! ij
         enddo ! ik
      enddo ! k
      KPRES_VV = KPRES_VV/Scell(NSC)%V ! kinetic part, to be multiplied by sigma below

      ! Pressure calculation finishing:
      call Det_3x3(Scell(NSC)%supce,Scell(NSC)%V) ! determinant of the super-cell is the volume, module "Algebra_tools"
      
      ! Construct full Sigma tensor:
      Sigma_tens = 0.0d0
      do i=1,3
         do j = 1,3
            !call d_detH_d_h_a_b_OLD(Scell(NSC)%supce,i,j,Sigma_tens(i,j)) ! calculates the derivatives of the determinant of h-vectors, "sigma" from PR, module "Algebra_tools"
            !print*, 'ONE:', Sigma_tens
            call d_detH_d_h_a_b(Scell(NSC)%supce,i,j,Sigma_tens(i,j)) ! calculates the derivatives of the determinant of h-vectors, "sigma" from PR, module "Algebra_tools"
            !print*, 'TWO:', Sigma_tens
         enddo
      enddo
      dsupce = Sigma_tens/matter%W_PR ! [A^2/kg] (sigma/wPR) from H.Jeschke Eq.(2.16), p.36
      ! If we are calculating the stress tensor and pressure:
      if (present(Sigma_tens_OUT)) Sigma_tens_OUT = Sigma_tens
      
      ! Multipliers:
      do i=1,3
         do j = 1,3
            KPRESint(i,j) = SUM(KPRES_VV(:,i)*Sigma_tens(j,:)) !*
            !KPRESint(i,j) = SUM(KPRES_VV(i,:)*Sigma_tens(:,j)) ! does not work well
         enddo ! j
      enddo ! i
      KPRESint = KPRESint/matter%W_PR
      
      ! Final total force as sum of three terms:
      if (present(Sigma_tens_OUT)) then	! construction of stress tensor and pressure require to exclude external part:
         do i=1,3
            do j = 1,3
               supce_forces%total(i,j) = KPRESint(i,j) - supce_forces%total(i,j) !- (matter%p_ext*dsupce(i,j))*1d-40 ! 
            enddo ! j
         enddo ! i
      else ! Forces calculations require full expression:
         do i=1,3
            do j = 1,3
               supce_forces%total(i,j) = KPRESint(i,j) - supce_forces%total(i,j) - (matter%p_ext*dsupce(i,j))*1d-40 ! 
            enddo ! j
         enddo ! i
      endif
      
      supce_forces%total0 = supce_forces%total ! save it for futher calculations of velocities
   endif
   nullify(Mass)
end subroutine super_cell_forces



pure subroutine check_periodic_boundaries(matter, Scell, NSC)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in) :: matter	! material parameters
   integer k, nat
   nat = size(Scell(NSC)%MDatoms)
   do k = 1,nat ! periodic boundary conditions:
      if ( (Scell(NSC)%MDatoms(k)%S(1) .GT. 1.0d0) .or. (Scell(NSC)%MDatoms(k)%S(1) .LT. -0.0d0) ) then
!          print*, ' X', Scell(NSC)%MDatoms(k)%S(1), FLOOR(Scell(NSC)%MDatoms(k)%S(1))
         Scell(NSC)%MDatoms(k)%S(1) = Scell(NSC)%MDatoms(k)%S(1) - FLOOR(Scell(NSC)%MDatoms(k)%S(1))
         Scell(NSC)%MDatoms(k)%S0(1) = Scell(NSC)%MDatoms(k)%S0(1) - FLOOR(Scell(NSC)%MDatoms(k)%S(1))
      endif
      if ( (Scell(NSC)%MDatoms(k)%S(2) .GT. 1.0d0) .or. (Scell(NSC)%MDatoms(k)%S(2) .LT. -0.0d0) ) then
!          print*, ' Y', Scell(NSC)%MDatoms(k)%S(2), FLOOR(Scell(NSC)%MDatoms(k)%S(2))
         Scell(NSC)%MDatoms(k)%S(2) = Scell(NSC)%MDatoms(k)%S(2) - FLOOR(Scell(NSC)%MDatoms(k)%S(2))
         Scell(NSC)%MDatoms(k)%S0(2) = Scell(NSC)%MDatoms(k)%S0(2) - FLOOR(Scell(NSC)%MDatoms(k)%S(2))
      endif
      if ( (Scell(NSC)%MDatoms(k)%S(3) .GT. 1.0d0) .or. (Scell(NSC)%MDatoms(k)%S(3) .LT. -0.0d0) ) then
!          print*, ' Z', Scell(NSC)%MDatoms(k)%S(3), FLOOR(Scell(NSC)%MDatoms(k)%S(3))
         Scell(NSC)%MDatoms(k)%S(3) = Scell(NSC)%MDatoms(k)%S(3) - FLOOR(Scell(NSC)%MDatoms(k)%S(3))
         Scell(NSC)%MDatoms(k)%S0(3) = Scell(NSC)%MDatoms(k)%S0(3) - FLOOR(Scell(NSC)%MDatoms(k)%S(3))
      endif
   enddo ! k
   call Coordinates_rel_to_abs(Scell, NSC)
end subroutine check_periodic_boundaries


pure subroutine check_periodic_boundaries_single(matter, Scell, NSC, k)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in) :: matter	! material parameters
   integer, intent(in) :: k
   if ( (Scell(NSC)%MDatoms(k)%S(1) .GT. 1.0d0) .or. (Scell(NSC)%MDatoms(k)%S(1) .LT. -0.0d0) ) then
      Scell(NSC)%MDatoms(k)%S(1) = Scell(NSC)%MDatoms(k)%S(1) - FLOOR(Scell(NSC)%MDatoms(k)%S(1))
      Scell(NSC)%MDatoms(k)%S0(1) = Scell(NSC)%MDatoms(k)%S0(1) - FLOOR(Scell(NSC)%MDatoms(k)%S(1))
   endif
   if ( (Scell(NSC)%MDatoms(k)%S(2) .GT. 1.0d0) .or. (Scell(NSC)%MDatoms(k)%S(2) .LT. -0.0d0) ) then
      Scell(NSC)%MDatoms(k)%S(2) = Scell(NSC)%MDatoms(k)%S(2) - FLOOR(Scell(NSC)%MDatoms(k)%S(2))
      Scell(NSC)%MDatoms(k)%S0(2) = Scell(NSC)%MDatoms(k)%S0(2) - FLOOR(Scell(NSC)%MDatoms(k)%S(2))
   endif
   if ( (Scell(NSC)%MDatoms(k)%S(3) .GT. 1.0d0) .or. (Scell(NSC)%MDatoms(k)%S(3) .LT. -0.0d0) ) then
      Scell(NSC)%MDatoms(k)%S(3) = Scell(NSC)%MDatoms(k)%S(3) - FLOOR(Scell(NSC)%MDatoms(k)%S(3))
      Scell(NSC)%MDatoms(k)%S0(3) = Scell(NSC)%MDatoms(k)%S0(3) - FLOOR(Scell(NSC)%MDatoms(k)%S(3))
   endif
   call Coordinates_rel_to_abs_single(Scell, NSC, k, .true.)
end subroutine check_periodic_boundaries_single



subroutine get_Ekin(Scell, matter)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   type(solid), intent(in) :: matter	! materil parameters
   integer NSC
   do NSC = 1, size(Scell)
      call get_kinetic_energy_abs(Scell, NSC, matter, Scell(NSC)%nrg)
      !print*, Scell(NSC)%nrg%At_kin, Scell(NSC)%nrg%E_supce
   enddo
end subroutine get_Ekin

subroutine get_kinetic_energy_abs(Scell, NSC, matter, nrg)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(in), target :: matter	! materil parameters
   type(Energies), intent(inout) :: nrg	! all energies
   !----------------------------
   real(8), pointer :: Mass
   real(8) :: V2, Vv, Ekin, Epot
   integer i, N, Nat
   N = size(Scell(NSC)%MDatoms)	! number of atoms
   Scell(NSC)%MDatoms(:)%Ekin = 0.0d0
   do i = 1,N	! all atoms:
      V2 = SUM(Scell(NSC)%MDatoms(i)%V(:)*Scell(NSC)%MDatoms(i)%V(:))*1d10 ! abs value of velocity [A/fs]^2 -> [m/s]^2
      Mass => matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma ! atomic mass
      !atoms(i)%Ekin = matter%Ma*V2/2.0d0/g_e ! [eV]
      !Scell(NSC)%MDatoms(i)%Ekin = matter%Atoms(Scell(NSC)%MDatoms(i)%KOA)%Ma*V2/2.0d0/g_e ! [eV]
      Scell(NSC)%MDatoms(i)%Ekin = Mass*V2/2.0d0/g_e ! [eV]
   enddo
   nrg%At_kin = SUM(Scell(NSC)%MDatoms(:)%Ekin)/dble(N)	! total atomic kinetic energy [eV/atom] CORRECT
   if (N > 2) then
      Scell(NSC)%TaeV = 2.0d0/(3.0d0*dble(N) - 6.0d0)*SUM(Scell(NSC)%MDatoms(:)%Ekin) ! Temperature [eV], Eq.(2.62) from H.Jeschke PhD thesis, p.49
   else	! use the eq. for nonperiodic boundaries:
      Scell(NSC)%TaeV = 2.0d0/(3.0d0*dble(N))*SUM(Scell(NSC)%MDatoms(:)%Ekin) ! Temperature [eV], Eq.(2.62) from H.Jeschke PhD thesis, p.49
   endif
   Scell(NSC)%Ta = Scell(NSC)%TaeV*g_kb	! [K]
   ! Temperature of different sublattices:
   if (.not.allocated(Scell(NSC)%Ta_sub)) allocate(Scell(NSC)%Ta_sub(matter%N_KAO))
   ! For all elements:
   if (matter%N_KAO > 1) then
      do i = 1, matter%N_KAO
         Nat = COUNT(MASK = (Scell(NSC)%MDatoms(:)%KOA == i)) ! how many atoms of this kind
         if (Nat > 2) then
            Scell(NSC)%Ta_sub(i) = 2.0d0/(3.0d0*dble(Nat) - 6.0d0)*SUM(Scell(NSC)%MDatoms(:)%Ekin, MASK = (Scell(NSC)%MDatoms(:)%KOA == i))
         elseif (Nat <= 0) then  ! no atoms, no temperature
            Scell(NSC)%Ta_sub(i) = 0.0d0
         else	! use the eq. for nonperiodic boundaries:
            Scell(NSC)%Ta_sub(i) = 2.0d0/(3.0d0*dble(Nat))*SUM(Scell(NSC)%MDatoms(:)%Ekin, MASK = (Scell(NSC)%MDatoms(:)%KOA == i))
         endif
         Scell(NSC)%Ta_sub(i) = Scell(NSC)%Ta_sub(i)*g_kb	! [K]
      enddo
   else ! there is only one atomic element:
      Scell(NSC)%Ta_sub(1) = Scell(NSC)%Ta  ! [K]
   endif
   
   ! CORRECT KINETIC ENERGY OF THE SUPERCELL:
!    Ekin = matter%W_PR*(SUM(Scell(NSC)%Vsupce(1,:)*Scell(NSC)%Vsupce(1,:)) + &
!                        SUM(Scell(NSC)%Vsupce(2,:)*Scell(NSC)%Vsupce(2,:)) + &
!                        SUM(Scell(NSC)%Vsupce(3,:)*Scell(NSC)%Vsupce(3,:)))*1d10/2.0d0/g_e ! kinetic part of the energy of the supercell [eV]
   Ekin = Supce_kin_energy(Scell, NSC, matter%W_PR)    ! function below

!    Ekin = matter%W_PR*(SUM(Scell(NSC)%Vsupce(:,1)*Scell(NSC)%Vsupce(:,1)) + &
!                        SUM(Scell(NSC)%Vsupce(:,2)*Scell(NSC)%Vsupce(:,2)) + &
!                        SUM(Scell(NSC)%Vsupce(:,3)*Scell(NSC)%Vsupce(:,3)))*1d10/2.0d0/g_e ! kinetic part of the energy of the supercell [eV]

   Epot = (matter%p_ext*Scell(NSC)%V)*1d-30/g_e	! potential part of the energy of the supercell [eV]
   nrg%E_supce = (Ekin + Epot)/dble(N) 	! total energy of the supercell [eV/atom]
   nullify(Mass)
end subroutine get_kinetic_energy_abs


pure function Supce_kin_energy(Scell, NSC, W_PR) result(Ekin)
   real(8) Ekin ! [eV] Super cell kinetic energy in Parrinello-Rahman method
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   real(8), intent(in) :: W_PR  ! PR mass of the supercell

   ! CORRECT KINETIC ENERGY OF THE SUPERCELL:
   Ekin = W_PR*(SUM(Scell(NSC)%Vsupce(1,:)*Scell(NSC)%Vsupce(1,:)) + &
                SUM(Scell(NSC)%Vsupce(2,:)*Scell(NSC)%Vsupce(2,:)) + &
                SUM(Scell(NSC)%Vsupce(3,:)*Scell(NSC)%Vsupce(3,:)))*1d10/2.0d0/g_e ! kinetic part of the energy of the supercell [eV]
end function Supce_kin_energy




subroutine get_energy_from_temperature(Na, Ta, Ekin)
   real(8), intent(in) :: Na ! number of atoms
   real(8), intent(in) :: Ta ! temperature of atoms [eV]
   real(8), intent(out) :: Ekin ! total kinetic energy of atoms [eV]
   Ekin = (3.0d0*Na - 6.0d0)/2.0d0*Ta ! kinetic energy in a box with periodic boundary
end subroutine get_energy_from_temperature



subroutine get_temperature_from_energy(Na, Ekin, Ta)
   real(8), intent(in) :: Na ! number of atoms
   real(8), intent(in) :: Ekin ! total kinetic energy of atoms [eV]
   real(8), intent(out) :: Ta ! temperature of atoms [eV]
   Ta = 2.0d0/(3.0d0*Na - 6.0d0)*Ekin ! temperature in a box with periodic boundary
end subroutine get_temperature_from_energy


subroutine get_mean_square_displacement(Scell, matter, MSD, MSDP, MSD_power)	! currently, it calculates mean displacement, without sqaring it
   type(Super_cell), dimension(:), intent(inout), target :: Scell	! super-cell with all the atoms inside
   type(Solid), intent(in) :: matter     ! material parameters
   real(8), intent(out) :: MSD	! [A^MSD_power] mean displacements average over all atoms
   real(8), dimension(:), allocatable, intent(out) :: MSDP ! [A] mean displacement of atoms for each sort of atoms in a compound
   integer, intent(in) :: MSD_power ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>)
   !-------------------------
   integer :: N, iat, ik, i,  j, k, Nat
   integer, pointer :: KOA
   real(8) :: zb(3), x, y, z, a_r, r1, x0, y0, z0
   real(8), dimension(:), pointer :: S, S0
   
   if (.not.allocated(MSDP)) allocate(MSDP(matter%N_KAO))
   
   N = size(Scell(1)%MDAtoms)	! number of atoms
   
   ! Get equilibrium relative coordinates from given absolute coordinates inside of the current supercell:
   call get_coords_in_new_supce(Scell, 1)	! below (S_eq arae updated, R_eq do not change)
   
   MSD = 0.0d0	! to start with
   MSDP = 0.0d0 ! to start with
   do iat = 1, N	! for all atoms
      KOA => Scell(1)%MDatoms(iat)%KOA
      S => Scell(1)%MDAtoms(iat)%S(:)
      S0 => Scell(1)%MDAtoms(iat)%S_eq(:)
      a_r = 1.0d31	! just to start from
      ! For the case of periodical boundaries:
      do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
         ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell: 
         ! periodic boundary conditions
         zb(1) = dble(i)
         do j = -1,1
            zb(2) = dble(j)
            do k = -1,1
               zb(3) = dble(k)
               x0 = 0.0d0
               y0 = 0.0d0
               z0 = 0.0d0
               do ik = 1,3
!                   x0 = x0 + (S(ik) - S0(ik) + zb(ik))*Scell(1)%supce(1,ik) ! incorrect
!                   y0 = y0 + (S(ik) - S0(ik) + zb(ik))*Scell(1)%supce(2,ik)
!                   z0 = z0 + (S(ik) - S0(ik) + zb(ik))*Scell(1)%supce(3,ik)
                  x0 = x0 + (S(ik) - S0(ik) + zb(ik))*Scell(1)%supce(ik,1) ! correct
                  y0 = y0 + (S(ik) - S0(ik) + zb(ik))*Scell(1)%supce(ik,2)
                  z0 = z0 + (S(ik) - S0(ik) + zb(ik))*Scell(1)%supce(ik,3)
               enddo ! ik
               r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
               if (r1 .LT. a_r) then
                  x = x0
                  y = y0
                  z = z0
                  a_r = r1
               endif !  (r1 .LT. a_r)
            enddo ! k
         enddo ! j
      enddo ! i
      MSD = MSD + a_r**MSD_power ! mean displacement^N
      MSDP(KOA) = MSDP(KOA) + a_r**MSD_power    ! mean displacement^N
   enddo ! iat
   MSD = MSD/dble(N)	! averaged over all atoms
   ! For all elements:
   do i = 1, matter%N_KAO
      ! how many atoms of this kind are in the supercell:
      Nat = COUNT(MASK = (Scell(1)%MDatoms(:)%KOA == i))
      if (Nat > 0) then
         MSDP(i) = MSDP(i) / dble(Nat)
      else
         MSDP(i) = 0.0d0
      endif
   enddo
   
   nullify(S,S0,KOA)
end subroutine get_mean_square_displacement


subroutine get_coords_in_new_supce(Scell, NSC) !  (S_eq are updated, R_eq do not change)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   real(8) :: sx, sy, sz
   real(8), dimension(3,3) :: supce_inv
   integer j, N, ik
   N = size(Scell(NSC)%MDatoms)
   call Invers_3x3(Scell(NSC)%supce, supce_inv, 'get_coords_in_new_supce')	! module "Algebra_tools"
   do j = 1,N	! all atoms:
      sx = 0.0d0
      sy = 0.0d0
      sz = 0.0d0
      do ik = 1,3
!          sx = sx + Scell(NSC)%MDatoms(j)%R_eq(ik)*supce_inv(1,ik)
!          sy = sy + Scell(NSC)%MDatoms(j)%R_eq(ik)*supce_inv(2,ik)
!          sz = sz + Scell(NSC)%MDatoms(j)%R_eq(ik)*supce_inv(3,ik)
         sx = sx + Scell(NSC)%MDatoms(j)%R_eq(ik)*supce_inv(ik,1) ! correct
         sy = sy + Scell(NSC)%MDatoms(j)%R_eq(ik)*supce_inv(ik,2)
         sz = sz + Scell(NSC)%MDatoms(j)%R_eq(ik)*supce_inv(ik,3)
      enddo ! ik
      Scell(NSC)%MDatoms(j)%S_eq(1) = sx
      Scell(NSC)%MDatoms(j)%S_eq(2) = sy
      Scell(NSC)%MDatoms(j)%S_eq(3) = sz
   enddo ! j
end subroutine get_coords_in_new_supce




subroutine get_kinetic_energy_rel(Scell, NSC, matter, nrg)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(solid), intent(inout) :: matter	! materil parameters
   type(Energies), intent(inout) :: nrg	! all energies
   call velocities_rel_to_abs(Scell, NSC)	! convert velocities into the absolute ones
   call get_kinetic_energy_abs(Scell, NSC, matter, nrg)	! use subroutine for the absolute velocities
end subroutine get_kinetic_energy_rel


subroutine velocities_rel_to_abs(Scell, NSC)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   real(8) :: v(3)
   integer k, ik, N
   N = size(Scell(NSC)%MDatoms)
   do k = 1,N ! periodic boundary conditions:
      v = 0.0d0
      do ik = 1,3
!          v(:) = v(:) + Scell(NSC)%MDatoms(k)%SV(ik)*Scell(NSC)%supce(:,ik) ! + Sco(ik,k)*Vsupce(1,ik)
         v(:) = v(:) + Scell(NSC)%MDatoms(k)%SV(ik)*Scell(NSC)%supce(ik,:) ! + Sco(ik,k)*Vsupce(1,ik)
      enddo ! ik
      Scell(NSC)%MDatoms(k)%V(:) = v(:) ! [A/fs]
   enddo ! k
end subroutine velocities_rel_to_abs


subroutine velocities_abs_to_rel(Scell, NSC, if_old)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   logical, optional :: if_old ! then do it for the previous time-step too
   real(8) v(3), dsupce(3,3)
   integer i, ik, N
   N = size(Scell(NSC)%MDatoms)
   !Relative velocities:
   call Invers_3x3(Scell(NSC)%supce, dsupce, 'velocities_abs_to_rel') ! from module "Algebra_tools"
   do i = 1, N
      v = 0.0d0
      do ik = 1,3
!          v(:) = v(:) + Scell(NSC)%MDatoms(i)%V(ik)*dsupce(:,ik)
         v(:) = v(:) + Scell(NSC)%MDatoms(i)%V(ik)*dsupce(ik,:)
      enddo ! ik
      Scell(NSC)%MDatoms(i)%SV(:) = v(:)
   enddo
   if (present(if_old)) then
      call Invers_3x3(Scell(NSC)%supce0, dsupce, 'velocities_abs_to_rel (2)') ! from module "Algebra_tools"
      do i = 1, N
         v = 0.0d0
         do ik = 1,3
!             v(:) = v(:) + Scell(NSC)%MDatoms(i)%V0(ik)*dsupce(:,ik)
            v(:) = v(:) + Scell(NSC)%MDatoms(i)%V0(ik)*dsupce(ik,:)
         enddo ! ik
         Scell(NSC)%MDatoms(i)%SV0(:) = v(:)
      enddo
   endif
end subroutine velocities_abs_to_rel


pure subroutine Coordinates_rel_to_abs(Scell, NSC, if_old)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   logical, intent(in), optional :: if_old ! then do it for the previous time-step too
   real(8) x, y, z
   integer j, N, ik
   N = size(Scell(NSC)%MDatoms)
   do j = 1,N	! all atoms:
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
      do ik = 1,3
!          x = x + Scell(NSC)%MDatoms(j)%S(ik)*Scell(NSC)%supce(1,ik)
!          y = y + Scell(NSC)%MDatoms(j)%S(ik)*Scell(NSC)%supce(2,ik)
!          z = z + Scell(NSC)%MDatoms(j)%S(ik)*Scell(NSC)%supce(3,ik)
         x = x + Scell(NSC)%MDatoms(j)%S(ik)*Scell(NSC)%supce(ik,1)
         y = y + Scell(NSC)%MDatoms(j)%S(ik)*Scell(NSC)%supce(ik,2)
         z = z + Scell(NSC)%MDatoms(j)%S(ik)*Scell(NSC)%supce(ik,3)
      enddo ! ik
      Scell(NSC)%MDatoms(j)%R(1) = x
      Scell(NSC)%MDatoms(j)%R(2) = y
      Scell(NSC)%MDatoms(j)%R(3) = z
   enddo ! j
   if (present(if_old)) then
      do j = 1,N	! all atoms:
         x = 0.0d0
         y = 0.0d0
         z = 0.0d0
         do ik = 1,3
!             x = x + Scell(NSC)%MDatoms(j)%S0(ik)*Scell(NSC)%supce0(1,ik)
!             y = y + Scell(NSC)%MDatoms(j)%S0(ik)*Scell(NSC)%supce0(2,ik)
!             z = z + Scell(NSC)%MDatoms(j)%S0(ik)*Scell(NSC)%supce0(3,ik)
            x = x + Scell(NSC)%MDatoms(j)%S0(ik)*Scell(NSC)%supce0(ik,1)
            y = y + Scell(NSC)%MDatoms(j)%S0(ik)*Scell(NSC)%supce0(ik,2)
            z = z + Scell(NSC)%MDatoms(j)%S0(ik)*Scell(NSC)%supce0(ik,3)
         enddo ! ik
         Scell(NSC)%MDatoms(j)%R0(1) = x
         Scell(NSC)%MDatoms(j)%R0(2) = y
         Scell(NSC)%MDatoms(j)%R0(3) = z
      enddo ! j
   endif
end subroutine Coordinates_rel_to_abs



subroutine Coordinates_abs_to_rel(Scell, NSC, if_old)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   logical, optional :: if_old ! then do it for the previous time-step too
   real(8) S(3), dsupce(3,3)
   integer i, ik, N
   N = size(Scell(NSC)%MDatoms)
   !Relative velocities:
   call Invers_3x3(Scell(NSC)%supce, dsupce, 'Coordinates_abs_to_rel') ! from module "Algebra_tools"
   do i = 1, N
      S = 0.0d0
      do ik = 1,3
         S(:) = S(:) + Scell(NSC)%MDatoms(i)%R(ik)*dsupce(ik,:)
      enddo ! ik
      Scell(NSC)%MDatoms(i)%S(:) = S(:)
   enddo
   if (present(if_old)) then
      call Invers_3x3(Scell(NSC)%supce0, dsupce, 'Coordinates_abs_to_rel (2)') ! from module "Algebra_tools"
      do i = 1, N
         S = 0.0d0
         do ik = 1,3
            S(:) = S(:) + Scell(NSC)%MDatoms(i)%R0(ik)*dsupce(ik,:)
         enddo ! ik
         Scell(NSC)%MDatoms(i)%S0(:) = S(:)
      enddo
   endif
end subroutine Coordinates_abs_to_rel




pure subroutine Coordinates_rel_to_abs_single(Scell, NSC, i_in, if_old)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   integer, intent(in) :: i_in   ! single atom to do the coordinates
   logical, intent(in), optional :: if_old ! then do it for the previous time-step too
   real(8) :: x, y, z
   integer :: ik
   x = 0.0d0
   y = 0.0d0
   z = 0.0d0
   do ik = 1,3
      x = x + Scell(NSC)%MDatoms(i_in)%S(ik)*Scell(NSC)%supce(ik,1)
      y = y + Scell(NSC)%MDatoms(i_in)%S(ik)*Scell(NSC)%supce(ik,2)
      z = z + Scell(NSC)%MDatoms(i_in)%S(ik)*Scell(NSC)%supce(ik,3)
   enddo ! ik
   Scell(NSC)%MDatoms(i_in)%R(1) = x
   Scell(NSC)%MDatoms(i_in)%R(2) = y
   Scell(NSC)%MDatoms(i_in)%R(3) = z

   if (present(if_old)) then
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
      do ik = 1,3
         x = x + Scell(NSC)%MDatoms(i_in)%S0(ik)*Scell(NSC)%supce0(ik,1)
         y = y + Scell(NSC)%MDatoms(i_in)%S0(ik)*Scell(NSC)%supce0(ik,2)
         z = z + Scell(NSC)%MDatoms(i_in)%S0(ik)*Scell(NSC)%supce0(ik,3)
      enddo ! ik
      Scell(NSC)%MDatoms(i_in)%R0(1) = x
      Scell(NSC)%MDatoms(i_in)%R0(2) = y
      Scell(NSC)%MDatoms(i_in)%R0(3) = z
   endif
end subroutine Coordinates_rel_to_abs_single



subroutine Coordinates_abs_to_rel_single(Scell, NSC, i_in, if_old)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   integer, intent(in) :: i_in   ! single atom to do the coordinates
   logical, optional :: if_old ! then do it for the previous time-step too
   real(8) S(3), dsupce(3,3)
   integer ik
   !Relative velocities:
   call Invers_3x3(Scell(NSC)%supce, dsupce, 'Coordinates_abs_to_rel_single') ! from module "Algebra_tools"
   S = 0.0d0
   do ik = 1,3
      S(:) = S(:) + Scell(NSC)%MDatoms(i_in)%R(ik)*dsupce(ik,:)
   enddo ! ik
   Scell(NSC)%MDatoms(i_in)%S(:) = S(:)

   if (present(if_old)) then
      call Invers_3x3(Scell(NSC)%supce0, dsupce, 'Coordinates_abs_to_rel_single (2)') ! from module "Algebra_tools"
      S = 0.0d0
      do ik = 1,3
         S(:) = S(:) + Scell(NSC)%MDatoms(i_in)%R0(ik)*dsupce(ik,:)
      enddo ! ik
      Scell(NSC)%MDatoms(i_in)%S0(:) = S(:)
   endif
end subroutine Coordinates_abs_to_rel_single




subroutine Reciproc_rel_to_abs(ksx, ksy, ksz, Scell, NSC, kx, ky, kz)
   real(8), intent(in) :: ksx, ksy, ksz ! relative reciprocal vector
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   real(8), intent(out) :: kx, ky, kz   ! absolute reciprocal vector [1/A]
!    kx = ksx*Scell(NSC)%k_supce(1,1) + ksy*Scell(NSC)%k_supce(1,2) + ksz*Scell(NSC)%k_supce(1,3)
!    ky = ksx*Scell(NSC)%k_supce(2,1) + ksy*Scell(NSC)%k_supce(2,2) + ksz*Scell(NSC)%k_supce(2,3)
!    kz = ksx*Scell(NSC)%k_supce(3,1) + ksy*Scell(NSC)%k_supce(3,2) + ksz*Scell(NSC)%k_supce(3,3)
   kx = ksx*Scell(NSC)%k_supce(1,1) + ksy*Scell(NSC)%k_supce(2,1) + ksz*Scell(NSC)%k_supce(3,1)
   ky = ksx*Scell(NSC)%k_supce(1,2) + ksy*Scell(NSC)%k_supce(2,2) + ksz*Scell(NSC)%k_supce(3,2)
   kz = ksx*Scell(NSC)%k_supce(1,3) + ksy*Scell(NSC)%k_supce(2,3) + ksz*Scell(NSC)%k_supce(3,3)
end subroutine Reciproc_rel_to_abs


subroutine Convert_reciproc_rel_to_abs(ksx, ksy, ksz, k_supce, kx, ky, kz)
   real(8), intent(in) :: ksx, ksy, ksz ! relative reciprocal vector
   real(8), dimension(3,3), intent(in) :: k_supce
   real(8), intent(out) :: kx, ky, kz   ! absolute reciprocal vector [1/A]
   kx = ksx*k_supce(1,1) + ksy*k_supce(2,1) + ksz*k_supce(3,1)
   ky = ksx*k_supce(1,2) + ksy*k_supce(2,2) + ksz*k_supce(3,2)
   kz = ksx*k_supce(1,3) + ksy*k_supce(2,3) + ksz*k_supce(3,3)
end subroutine Convert_reciproc_rel_to_abs



subroutine deflect_velosity(u0, v0, w0, theta, phi, u, v, w)    ! Eq.(1.131), p.37 [1]
   real(8), intent(in) :: u0, v0, w0     ! cosine directions of the old velosity
   real(8), intent(in) :: theta, phi     ! polar (0,Pi) and azimuthal (0,2Pi) angles
   real(8), intent(out) :: u, v, w       ! new cosine directions
   real(8) :: sin_theta, cos_theta, sin_phi, cos_phi, one_w, eps, sin_t_w, temp
   eps = 1.0d-8 ! margin of acceptance of w being along Z
   sin_theta = sin(theta)
   cos_theta = cos(theta)
   sin_phi = sin(phi)
   cos_phi = cos(phi)
   if ( abs(abs(w0)-1.0d0) < eps ) then   ! motion parallel to Z
      u = w0*sin_theta*cos_phi
      v = w0*sin_theta*sin_phi
      w = w0*cos_theta
   else ! any other direction of motion
      one_w = sqrt(1.0d0 - w0*w0)
      sin_t_w = sin_theta/one_w
      u = u0*cos_theta + sin_t_w*(u0*w0*cos_phi - v0*sin_phi)
      v = v0*cos_theta + sin_t_w*(v0*w0*cos_phi + u0*sin_phi)
      w = w0*cos_theta - one_w*sin_theta*cos_phi
   endif

   temp = sqrt(u*u + v*v + w*w)
   if (abs(temp-1.0d0) > eps) then  ! renormalize it:
      u = u/temp
      v = v/temp
      w = w/temp
   endif
end subroutine deflect_velosity



subroutine shortest_distance_NEW(Scell, i1, j1, a_r, x1, y1, z1, sx1, sy1, sz1, cell_x, cell_y, cell_z)
   type(Super_cell), intent(in), target :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: i1, j1 ! atomic numbers
   real(8), intent(out) ::  a_r	! [A] shortest distance between the two atoms within supercell with periodic boundaries
   real(8), intent(out), optional :: x1, y1, z1		! [A] projections of the shortest distance
   real(8), intent(out), optional :: sx1, sy1, sz1 	! relative projections of the shortest distance
   integer, intent(out), optional :: cell_x, cell_y, cell_z ! cell numbers
   real(8) x, y, z, zb(3), r, x0, y0, z0, r1
   integer i, j, k, ik
   type(Atom), dimension(:), pointer :: atoms	! array of atoms in the supercell

   atoms => Scell%MDAtoms
   x = 0.0d0
   y = 0.0d0
   z = 0.0d0
   if (i1 == j1) then ! it's the same atom:
      a_r = 0.0d0
      ! save the shortest distance projections
      if (present(x1)) x1 = x
      if (present(y1)) y1 = y
      if (present(z1)) z1 = z
      ! save the relative shortest distance projections
      if (present(sx1)) sx1 = 0.0d0
      if (present(sy1)) sy1 = 0.0d0
      if (present(sz1)) sz1 = 0.0d0
      ! save the cell numbers
      if (present(cell_x)) cell_x = 0
      if (present(cell_y)) cell_y = 0
      if (present(cell_z)) cell_z = 0
   else
      ! For the case of periodic boundaries:
      do ik = 1,3
         x = x + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell%supce(ik,1)
         y = y + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell%supce(ik,2)
         z = z + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell%supce(ik,3)
      enddo ! ik
      a_r = DSQRT(x*x + y*y + z*z)
      if (present(x1)) x1 = x
      if (present(y1)) y1 = y
      if (present(z1)) z1 = z
      if (present(sx1)) sx1 = atoms(i1)%S(1) - atoms(j1)%S(1)
      if (present(sy1)) sy1 = atoms(i1)%S(2) - atoms(j1)%S(2)
      if (present(sz1)) sz1 = atoms(i1)%S(3) - atoms(j1)%S(3)
      if (present(cell_x)) cell_x = 0
      if (present(cell_y)) cell_y = 0
      if (present(cell_z)) cell_z = 0

      do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
         ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell:
         ! periodic boundary conditions
         zb(1) = dble(i)
         do j =-1,1
            zb(2) = dble(j)
            do k = -1,1
               zb(3) = dble(k)
               x0 = 0.0d0
               y0 = 0.0d0
               z0 = 0.0d0
               do ik = 1,3
                  x0 = x0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell%supce(ik,1)
                  y0 = y0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell%supce(ik,2)
                  z0 = z0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell%supce(ik,3)
               enddo ! ik
               r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
               if (r1 <= a_r) then
                  x = x0
                  y = y0
                  z = z0
                  a_r = r1
                  if (present(x1)) x1 = x
                  if (present(y1)) y1 = y
                  if (present(z1)) z1 = z
                  if (present(sx1)) sx1 = atoms(i1)%S(1) - atoms(j1)%S(1) + zb(1)
                  if (present(sy1)) sy1 = atoms(i1)%S(2) - atoms(j1)%S(2) + zb(2)
                  if (present(sz1)) sz1 = atoms(i1)%S(3) - atoms(j1)%S(3) + zb(3)
                  if (present(cell_x)) cell_x = i
                  if (present(cell_y)) cell_y = j
                  if (present(cell_z)) cell_z = k
               endif
            enddo ! k
         enddo ! j
      enddo ! i
   endif ! i1 = j1
end subroutine shortest_distance_NEW



subroutine shortest_distance_OLD(Scell, NSC, atoms, i1, j1, a_r, x1, y1, z1, sx1, sy1, sz1, cell_x, cell_y, cell_z)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   integer, intent(in) :: i1, j1 ! atomic numbers
   real(8), intent(out) ::  a_r	! [A] shortest distance between the two atoms within supercell with periodic boundaries
   real(8), intent(out), optional :: x1, y1, z1		! [A] projections of the shortest distance
   real(8), intent(out), optional :: sx1, sy1, sz1 	! relative projections of the shortest distance
   integer, intent(out), optional :: cell_x, cell_y, cell_z ! cell numbers
   real(8) x, y, z, zb(3), r, x0, y0, z0, r1
   integer i, j, k, ik

  x = 0.0d0
  y = 0.0d0
  z = 0.0d0
  if (i1 .EQ. j1) then ! it's the same atom:
   a_r = 0.0d0
   if (present(x1) .AND. present(y1) .AND. present(z1)) then	! save the shortest distance projections
      x1 = x
      y1 = y
      z1 = z
   endif
   if (present(sx1) .AND. present(sy1) .AND. present(sz1)) then	! save the shortest distance projections
      sx1 = 0.0d0
      sy1 = 0.0d0
      sz1 = 0.0d0
   endif
   if (present(cell_x) .AND. present(cell_y) .AND. present(cell_z)) then	! save the cell numbers:
      cell_x = 0
      cell_y = 0
      cell_z = 0
   endif
  else
   ! For the case of periodic boundaries:
   do ik = 1,3
!       x = x + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell(NSC)%supce(1,ik)
!       y = y + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell(NSC)%supce(2,ik)
!       z = z + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell(NSC)%supce(3,ik)
      x = x + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell(NSC)%supce(ik,1)   ! correct
      y = y + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell(NSC)%supce(ik,2)
      z = z + (atoms(i1)%S(ik) - atoms(j1)%S(ik))*Scell(NSC)%supce(ik,3)
   enddo ! ik
   a_r = DSQRT(x*x + y*y + z*z)
   if (present(x1) .AND. present(y1) .AND. present(z1)) then	! save the shortest distance projections
      x1 = x
      y1 = y
      z1 = z
   endif
   if (present(sx1) .AND. present(sy1) .AND. present(sz1)) then	! save the shortest distance projections
      sx1 = atoms(i1)%S(1) - atoms(j1)%S(1)
      sy1 = atoms(i1)%S(2) - atoms(j1)%S(2)
      sz1 = atoms(i1)%S(3) - atoms(j1)%S(3)
   endif
   if (present(cell_x) .AND. present(cell_y) .AND. present(cell_z)) then	! save the cell numbers:
      cell_x = 0
      cell_y = 0
      cell_z = 0
   endif
   do i = -1,1 ! if the distance between the atoms is more than a half of supercell, we account for
      ! interaction with the atom not from this, but from the neigbour ("mirrored") supercell:
      ! periodic boundary conditions.
      zb(1) = dble(i)
      do j =-1,1
         zb(2) = dble(j)
         do k = -1,1
            zb(3) = dble(k)
            x0 = 0.0d0
            y0 = 0.0d0
            z0 = 0.0d0
            do ik = 1,3
!                x0 = x0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(1,ik) ! incorrect
!                y0 = y0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(2,ik)
!                z0 = z0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(3,ik)
               x0 = x0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(ik,1) ! correct
               y0 = y0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(ik,2)
               z0 = z0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(ik,3)
            enddo ! ik
            r1 = DSQRT(x0*x0 + y0*y0 + z0*z0)
            if (r1 <= a_r) then
               x = x0
               y = y0
               z = z0
               a_r = r1
               if (present(x1) .AND. present(y1) .AND. present(z1)) then	! save the shortest distance projections
                  x1 = x
                  y1 = y
                  z1 = z
               endif
               if (present(sx1) .AND. present(sy1) .AND. present(sz1)) then	! save the shortest distance projections
                  sx1 = atoms(i1)%S(1) - atoms(j1)%S(1) + zb(1)
                  sy1 = atoms(i1)%S(2) - atoms(j1)%S(2) + zb(2)
                  sz1 = atoms(i1)%S(3) - atoms(j1)%S(3) + zb(3)
               endif
               if (present(cell_x) .AND. present(cell_y) .AND. present(cell_z)) then	! save the cell numbers:
                  cell_x = i
                  cell_y = j
                  cell_z = k
               endif
            endif
         enddo ! k
      enddo ! j
   enddo ! i
  endif ! i1 = j1
end subroutine shortest_distance_OLD



subroutine distance_to_given_cell(Scell, NSC, atoms, zb, i1, j1, R, x, y, z, sx, sy, sz)
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   real(8), dimension(:), intent(in) :: zb ! super cell indices
   integer, intent(in) :: i1, j1 ! indices of atoms distance between which we are looking for
   real(8), intent(out) :: R ! [A] distance between the given atoms
   real(8), intent(out), optional :: x, y, z ! [A] projection of that distance
   real(8), intent(out), optional :: sx, sy, sz ! relative projection of that distance
   real(8) :: x0, y0, z0
   integer :: ik
   x0 = 0.0d0
   y0 = 0.0d0
   z0 = 0.0d0
   do ik = 1,3
!       x0 = x0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(1,ik) ! distance in X
!       y0 = y0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(2,ik) ! distance in Y
!       z0 = z0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(3,ik) ! distance in Z
      x0 = x0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(ik,1) ! distance in X
      y0 = y0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(ik,2) ! distance in Y
      z0 = z0 + (atoms(i1)%S(ik) - atoms(j1)%S(ik) + zb(ik))*Scell(NSC)%supce(ik,3) ! distance in Z
   enddo ! ik
   R = DSQRT(x0*x0 + y0*y0 + z0*z0) ! [A] get out distance
   ! get out projections too:
   if (present(x)) x=x0
   if (present(y)) y=y0
   if (present(z)) z=z0
   ! save the shortest distance projections:
   if (present(sx)) sx = atoms(i1)%S(1) - atoms(j1)%S(1) + zb(1)
   if (present(sy)) sy = atoms(i1)%S(2) - atoms(j1)%S(2) + zb(2)
   if (present(sz)) sz = atoms(i1)%S(3) - atoms(j1)%S(3) + zb(3)
end subroutine distance_to_given_cell



subroutine Find_nearest_neighbours(Scell, rm, numpar)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! super-cell with all the atoms inside
   real(8), intent(in) :: rm ! potentials cut-off distance [A]
   type(Numerics_param), intent(in) :: numpar	! numerical parameters
   !==========================================================
   integer NSC ! number of super-cell
   integer n, i1, j1, ik, i, j, k, coun, coun_user
   real(8) x,y,z, a_r, x0, y0, z0, zb(3), r1, sx, sy, sz
   !rm = TB_H%rm ! [A], that's the range of potential, atoms with r>rm are not interacting
   do NSC = 1, size(Scell) ! for all supercells
      n = size(Scell(NSC)%MDatoms) ! number of atoms in this supercell
      Scell(NSC)%Near_neighbor_size = 0
      Scell(NSC)%Near_neighbor_list = 0
      Scell(NSC)%Near_neighbor_dist = 1d8
      do i1 = 1, n
         coun = 0
         coun_user = 0
         do j1=1, n
             if (j1 .NE. i1) then
                call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)

!                 Testing:
!                 call shortest_distance_slow(Scell(NSC), i1, j1, a_r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
!                 print*, 'OLD:', a_r, x, y, z, sx, sy, sz
!                 call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
!                 print*, 'NEW:', a_r, x, y, z, sx, sy, sz

!                 if (a_r .LE. rm) then   ! this atoms do interact:
                if (a_r < rm) then   ! this atoms do interact:
                   coun = coun + 1
                   Scell(NSC)%Near_neighbor_list(i1, coun) = j1   ! it interacts with this atom
                   Scell(NSC)%Near_neighbor_dist(i1, coun,1) = x  ! at this distance, X
                   Scell(NSC)%Near_neighbor_dist(i1, coun,2) = y  ! at this distance, Y
                   Scell(NSC)%Near_neighbor_dist(i1, coun,3) = z  ! at this distance, Z
                   Scell(NSC)%Near_neighbor_dist(i1, coun,4) = a_r  ! at this distance, R
                   Scell(NSC)%Near_neighbor_dist_s(i1, coun,1) = sx  ! at this distance, SX
                   Scell(NSC)%Near_neighbor_dist_s(i1, coun,2) = sy  ! at this distance, SY
                   Scell(NSC)%Near_neighbor_dist_s(i1, coun,3) = sz  ! at this distance, SZ
!                    Scell(NSC)%Near_neighbor_size(i1) = coun	! that's how many nearest neighbours there are for this atom
!                    if (ABS(Scell(NSC)%Near_neighbor_dist(i1, coun,4)) .LT. 1.0) then
!                       print*, 'NND1', Scell(NSC)%Near_neighbor_dist(i1, coun,1)
!                       print*, 'NND2', Scell(NSC)%Near_neighbor_dist(i1, coun,2)
!                       print*, 'NND3', Scell(NSC)%Near_neighbor_dist(i1, coun,3)
!                       print*, 'NND4', Scell(NSC)%Near_neighbor_dist(i1, coun,4)
!                    endif ! (ABS(Near_neighbor_dist(i1, coun,4)) .LT. 1.0)
                endif ! (a_r .LE. rm)
                
                if (numpar%save_NN) then   ! if user wants to study number of nearest neighbors within defined radius
                   if (a_r < numpar%NN_radius) then   ! this atoms are counted as nearest neighbors
                      coun_user = coun_user + 1
                   endif
                endif

            endif ! (j1 .NE. i1)
         enddo ! j1=1, n
         
         Scell(NSC)%Near_neighbor_size(i1) = coun   ! that's how many nearest neighbours there are for this atom
         if (numpar%save_NN) Scell(NSC)%Near_neighbors_user(i1) = coun_user ! number of neighbors within user-defined redius
         
         if (coun .LT. 1) then ! test: no nearest neighbours
!             print*, 'Atom', i1, ' has no near neighbours!'
         endif
      enddo! i1 = 1, n
   enddo !NSC
end subroutine Find_nearest_neighbours



subroutine pair_correlation_function(atoms, matter, Scell, NSC)
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(solid), intent(inout) :: matter	! materil parameters
   type(Super_cell), dimension(:), intent(in) :: Scell ! super-cell with all the atoms inside
   integer, intent(in) :: NSC ! number of super-cell
   real(8) r, dr, a_r
   integer n, i, m, k, j
   n = size(atoms)
   if (.not. allocated(matter%PCF)) then
      m = INT(Scell(NSC)%supce(1,1)*20)
      allocate(matter%PCF(2,m))
      r = 0.0d0
      dr = Scell(NSC)%supce(1,1)/dble(m)
      do i = 1,m
         r = r + dr 
         matter%PCF(1,i) = r
      enddo
   else
      m = size(matter%PCF,2)
   endif
   matter%PCF(2,:) = 0.0d0
   !$omp PARALLEL private(i,j,a_r,k)
   !$omp do schedule(dynamic)
   do i = 1, n	! trace all atoms
      do j = 1, n	! trace all neighbours for PCF
         if (i .NE. j) then
            call shortest_distance(Scell, NSC, atoms, i, j, a_r)
            if (a_r .GE. matter%PCF(1,m)) then
               k = m
            else
               call Find_in_array_monoton(matter%PCF, a_r, 1, k)	! module "Little_subroutines"
            endif
            matter%PCF(2,k) = matter%PCF(2,k) + 1.0d0
         endif
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp do
   do k = 1,m	! all points of the PCF
      matter%PCF(2,k) = matter%PCF(2,k)/(matter%PCF(1,k)*matter%PCF(1,k)) 
   enddo
   !$omp end do
   !$omp end parallel
   dr = matter%PCF(1,2) - matter%PCF(1,1)
   matter%PCF(2,:) = matter%PCF(2,:)/(4.0d0*g_Pi*dr)*Scell(NSC)%V/dble(Scell(NSC)%Na*Scell(NSC)%Na) ! normalizing per volume
end subroutine pair_correlation_function


! This subroutine for pressure calculations is NOT VALID for periodic boundary conditions,
! see e.g. [M.J. Louwerse, E.J. Baerends, Chemical Physics Letters 421 (2006) 138-141]
! https://doi.org/10.1016/j.cplett.2006.01.087
function Get_pressure_nonperiodic(Tat, V, Atoms) result(P)
   real(8) :: P	! Pressure, output
   real(8), intent(in) :: Tat	! Atomic temperature [K]
   real(8), intent(in) :: V	! Volume of the supercell [A^3]
   type(Atom), dimension(:), intent(in) :: Atoms	! type that contains coordinates and forces needed for pressure calculations
   !-------------------------------------
   real(8) :: Joul, Term1, Term2
   integer :: Nat, i
   Nat = size(Atoms)	! number of atoms in the simulation box
   Joul = Tat*g_kb_J	! Temperature [K] -> Energy [Joules]
   Term1 = dble(Nat)*Joul
   Term2 = 0.0d0
!$omp PARALLEL private(i)
!$omp do reduction( + : Term2)
   do i = 1, Nat
      Term2 = Term2 - SUM(Atoms(i)%R(:) * Atoms(i)%forces%total(:))
   enddo
!$omp end do
!$omp end parallel

   ! Convert units into SI units: R from [A] to [m]; Force from [eV/A] into [J/m].
   Term2 = 1.0d0/3.0d0*Term2*g_e
   P = (Term1 + Term2)/(V*1d-30)	! [Pa] total pressure
end function Get_pressure_nonperiodic
 

END MODULE Atomic_tools
