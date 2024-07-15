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
! References:
! [1] M. Graf and P. Vogl, Phys. Rev. B 51, 4940 (1995)
! [2] F. Trani, G. Cantele, D. Ninno, and G. Iadonisi, Phys. Rev. B 72, 075423 (2005)
! [3] P. Yeh "Optical Waves in Layered Media", (2005) ISBN: 978-0-471-73192-4
! [4] J. Graves, "Electronic and structural response of semiconductors to ultra-intense laser pulses" PhD thesis, Texas A&M University (1997)
! [5] Notes on Kubo-Greenwood CDF: https://www.openmx-square.org/tech_notes/Dielectric_Function_YTL.pdf
! [6] G. S. Demyanov et al., Phys. Rev. E 105 (2022) 035307 https://doi.org/10.1103/PhysRevE.105.035307


MODULE Optical_parameters
use Universal_constants
use Objects
use Algebra_tools, only : sym_diagonalize, numerical_delta
use Electron_tools, only : get_number_of_CB_electrons, set_Fermi, set_Erf_distribution, find_mu_from_N_T, get_Ce_and_mu, &
                           Diff_Fermi_E, Diff_Fermi_Te, update_cross_section
use TB_Fu, only : Complex_Hamil_tot_F
use TB_Pettifor, only : Complex_Hamil_tot
use TB_Molteni, only : Complex_Hamil_tot_Molteni
use TB_NRL, only : Complex_Hamil_NRL
use TB_DFTB, only : Complex_Hamil_DFTB, identify_DFTB_orbitals_per_atom
use TB, only : k_point_choice, construct_complex_Hamiltonian, get_coupling_matrix_elements
use Little_subroutines, only : deallocate_array, Find_in_array_monoton, d_Fermi_function
use MC_cross_sections, only : Mean_free_path, velosity_from_kinetic_energy
use Nonadiabatic, only : get_Mij2, get_nonadiabatic_Pij

#ifdef MPI_USED
   use MPI_subroutines, only : do_MPI_Allreduce, MPI_barrier_wrapper
#endif

#ifdef _OPENMP
   USE OMP_LIB, only : OMP_GET_THREAD_NUM
#endif

implicit none
PRIVATE

real(8), parameter :: m_gamm = 1.5d14  ! [1/s] gamma parameter

real(8), parameter :: m_inv_sqrt_pi = 1.0d0/sqrt(g_Pi)
real(8), parameter :: m_e_h = g_h/g_e
real(8), parameter :: m_prefac = g_e**2 / (g_h*g_me**2)



public :: get_optical_parameters, allocate_Eps_hw, get_Onsager_coeffs, get_Kubo_Greenwood_CDF, get_kappa_e_e, get_Onsager_dynamic


 contains


subroutine get_optical_parameters(numpar, matter, Scell, Err) ! optical coefficients, module "Optical_parameters"
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Solid), intent(inout) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err	! error save
   !=====================================
   integer NSC, i
   real(8), dimension(:), allocatable :: Re_CDF
   logical :: do_together

   do NSC = 1, size(Scell) ! for all super-cells
      ! Get the number of CB electrons:
      call get_number_of_CB_electrons(Scell, NSC) ! module "Electron_tools"

      ! Check if kappa and optical CDF can be calculated together:
      if (numpar%do_kappa .and. ((numpar%optic_model == 4) .or. (numpar%optic_model == 5)) ) then  ! yes, it's possible
         do_together = .true.
      else  ! no, separate subroutines required (more time consuming!)
         do_together = .false.
      endif

      if (do_together) then   ! One subroutine to get both, kappa and CDF:
         call get_Kubo_Greenwood_all_complex(numpar, matter, Scell, NSC, Scell(NSC)%eps%all_w, Err)    ! below

      else ! separate ones for kappa and CDF:
         ! CDF (if requested):
         select case (numpar%optic_model)
            case (1) ! within the Drude model
               call get_drude(numpar, Scell, NSC)  ! below
            case (2) ! Trani et al. PRB 72, 075423 (2005) -- This subroutine is TB-parameterization specific:
               call get_trani_all_complex(numpar, Scell, NSC, Scell(NSC)%eps%all_w, Err)    ! below
            case (3) ! Trani at the Gamma-point only
               call get_trani_all(numpar, Scell, NSC, Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%eps%all_w)   ! below
            case (4:5) ! Kubo-Greenwood Refs.[2] and [5]
               call get_Kubo_Greenwood_all_complex(numpar, matter, Scell, NSC, Scell(NSC)%eps%all_w, Err)    ! below
            case (-4) ! Graf-Vogl Ref.[2] (Not realy working...)
               call get_Graf_Vogl_all_complex(numpar, Scell, NSC, Scell(NSC)%eps%all_w, Err)    ! below
            case default ! no optical coefficients needed
            ! nothing to do
         end select ! (numpar%optic_model)

         ! Kappa (if requested): -- Obsolete. Kappa will be calculated separately in module "TB_complex"
         !if (numpar%do_kappa) call get_Kubo_Greenwood_all_complex(numpar, matter, Scell, NSC, Scell(NSC)%eps%all_w, Err)    ! below
      endif ! do_together
      
      !-------------------------------------
      ! Convergence with respect to the number of k-points is better for Im part than Re part of CDF,
      ! so you may use Kramers-Kronig relations to restore Re part:
      if (Scell(NSC)%eps%all_w .and. Scell(NSC)%eps%KK) then
!          do i = 1, size(Scell(NSC)%eps%Eps_hw(1,:))
!             write(*,'(f,es,es)') Scell(NSC)%eps%Eps_hw(1,i), Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i)
!          enddo
         call Kramers_Kronig_Re_from_Im(Scell(NSC)%eps%Eps_hw(1,:), Scell(NSC)%eps%Eps_hw(3,:), Re_CDF)	! function below
         Scell(NSC)%eps%Eps_hw(2,:) = Re_CDF(:)
         
         ! And the same for each component:
         call Kramers_Kronig_Re_from_Im(Scell(NSC)%eps%Eps_hw(1,:), Scell(NSC)%eps%Eps_hw(12,:), Re_CDF)	! function below
         Scell(NSC)%eps%Eps_hw(11,:) = Re_CDF(:)
         call Kramers_Kronig_Re_from_Im(Scell(NSC)%eps%Eps_hw(1,:), Scell(NSC)%eps%Eps_hw(14,:), Re_CDF)	! function below
         Scell(NSC)%eps%Eps_hw(13,:) = Re_CDF(:)
         call Kramers_Kronig_Re_from_Im(Scell(NSC)%eps%Eps_hw(1,:), Scell(NSC)%eps%Eps_hw(16,:), Re_CDF)	! function below
         Scell(NSC)%eps%Eps_hw(15,:) = Re_CDF(:)
         
         ! And then recalculate the optical properties:
         do i = 1, size(Scell(NSC)%eps%Eps_hw(1,:))
            call get_RTA_from_CDF(Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i), Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, Scell(NSC)%eps%Eps_hw(5,i), Scell(NSC)%eps%Eps_hw(6,i), Scell(NSC)%eps%Eps_hw(7,i))
            call get_n_k(Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i), Scell(NSC)%eps%Eps_hw(8,i), Scell(NSC)%eps%Eps_hw(9,i))	! transfer Re(e) and Im(e) into n and k
            call save_Eps_hw(Scell(NSC)%eps%Eps_hw, i, Scell(NSC)%eps%Eps_hw(1,i), Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i), Scell(NSC)%eps%Eps_hw(4,i), Scell(NSC)%eps%Eps_hw(5,i), Scell(NSC)%eps%Eps_hw(6,i), Scell(NSC)%eps%Eps_hw(7,i), Scell(NSC)%eps%Eps_hw(8,i), Scell(NSC)%eps%Eps_hw(9,i), cmplx(Scell(NSC)%eps%Eps_hw(11,i),Scell(NSC)%eps%Eps_hw(12,i)), cmplx(Scell(NSC)%eps%Eps_hw(13,i),Scell(NSC)%eps%Eps_hw(14,i)), cmplx(Scell(NSC)%eps%Eps_hw(15,i),Scell(NSC)%eps%Eps_hw(16,i)),  Scell(NSC)%eps%Eps_hw(10,i))
!             write(*,'(f,es,es)') Scell(NSC)%eps%Eps_hw(1,i), Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i)
         enddo
      endif
      !-------------------------------------
   enddo ! NSC = 1, size(Scell) ! for all super-cells
end subroutine get_optical_parameters



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kubo-GReenwood implementation in orthogonalized Hamiltonian (combined Refs.[2] and [5]):
subroutine get_Kubo_Greenwood_all_complex(numpar, matter, Scell, NSC, all_w, Err)  ! From Ref. [2]
   type (Numerics_param), intent(inout) :: numpar  ! numerical parameters, including drude-function
   type(Solid), intent(inout) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell   ! supercell with all the atoms as one object
   integer, intent(in) :: NSC    ! number of supercell
   logical, intent(in) :: all_w  ! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err   ! error save
   !--------------------
   !complex(8), dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   !complex(8), dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   complex, dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz    ! effective momentum operators
   complex, dimension(:,:), allocatable :: CHij    ! eigenvectors of the hamiltonian
   real(8) :: w, kx, ky, kz
   real(8), dimension(:), allocatable :: w_grid
   integer :: i, j, N, FN, ix, iy, iz, ixm, iym, izm, schem, Ngp, Nsiz
   real(8), dimension(:), allocatable :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable :: Eps_hw ! array of all eps vs hw
   real(8), dimension(:,:), allocatable :: Eps_hw_temp ! array of all eps vs hw
   real(8), dimension(:), allocatable :: kappa, kappa_ee, kappa_ee_temp, kappa_temp    ! electron heat conductivity vs Te
   real(8), dimension(:), allocatable :: kappa_mu_grid, kappa_mu_grid_temp, kappa_Ce_grid, kappa_Ce_grid_temp
   integer :: Nsiz_Te
   real(8) :: Te_min, Te_max, dTe

   ! Define electron temperature grid:
   Te_min = numpar%kappa_Te_min  ! [K]
   Te_max = numpar%kappa_Te_max  ! [K]
   dTe = numpar%kappa_dTe  ! [K]
   Nsiz_Te = (Te_max+dTe - Te_min)/dTe

   if (.not.allocated(Scell(NSC)%kappa_Te_grid)) allocate(Scell(NSC)%kappa_Te_grid(Nsiz_Te), source = 0.0d0)   ! Te [K]
   if (.not.allocated(Scell(NSC)%kappa_e_vs_Te)) allocate(Scell(NSC)%kappa_e_vs_Te(Nsiz_Te), source = 0.0d0)   ! kappa
   if (.not.allocated(Scell(NSC)%kappa_mu_grid)) allocate(Scell(NSC)%kappa_mu_grid(Nsiz_Te), source = 0.0d0)   ! mu [eV]
   if (.not.allocated(Scell(NSC)%kappa_Ce_grid)) allocate(Scell(NSC)%kappa_Ce_grid(Nsiz_Te), source = 0.0d0)   ! Ce [J/(m^3 K)]

   if (.not.allocated(kappa)) allocate(kappa(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_ee)) allocate(kappa_ee(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_mu_grid)) allocate(kappa_mu_grid(Nsiz_Te), source = 0.0d0)   ! mu [eV]
   if (.not.allocated(kappa_Ce_grid)) allocate(kappa_Ce_grid(Nsiz_Te), source = 0.0d0)   ! Ce [J/(m^3 K)]

   Scell(NSC)%kappa_Te_grid(1) = Te_min
   do i = 2, Nsiz_Te ! Set the electronic temperature grid
      Scell(NSC)%kappa_Te_grid(i) = Scell(NSC)%kappa_Te_grid(i-1) + dTe
   enddo

   ! Allocate the array of optical coefficients spectrum (if not allocated before):
   ! Define the grid size:
   if (Scell(NSC)%eps%all_w) then ! full spectrum:
      call allocate_Eps_hw(Scell(NSC)%eps%E_min, Scell(NSC)%eps%E_max, Scell(NSC)%eps%dE, Scell(NSC)%eps%Eps_hw) ! see below
      N = size(Scell(NSC)%eps%Eps_hw,2) ! use this size to define temporary arrays within this subroutine
      ! Grid of frequencies:
      allocate(w_grid(N),source = 0.0d0)
      w_grid(1) =  Scell(NSC)%eps%E_min*g_e/g_h ! [1/s] frequency starting point for the optical spectrum (from [eV])
      do i = 2, N
         w_grid(i) = w_grid(i-1) + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
      enddo
   else  ! single frequency
      N = 1
      allocate(w_grid(N),source = 0.0d0)
      if (.not.allocated(Scell(NSC)%eps%Eps_hw)) then
         allocate(Scell(NSC)%eps%Eps_hw(16,N)) ! different parameters saved for all energy grid points
      endif
      w_grid(1) = Scell(NSC)%eps%w
   endif
   if (.not.allocated(Eps_hw)) allocate(Eps_hw(16,N), source = 0.0d0) ! all are there

   ! For the user-defind number of k-points:
   ixm = numpar%ixm
   iym = numpar%iym
   izm = numpar%izm
   if (allocated(numpar%k_grid)) then
      schem = 1   ! user-defined grid is present
      Nsiz = size(numpar%k_grid,1)  ! size of the user provided grid
   else
      schem = 0   ! no user-defined grid, use default
      Nsiz = ixm*iym*izm
   endif

   kappa = 0.0d0  ! to start with
   kappa_ee = 0.0d0  ! to start with

#ifdef MPI_USED   ! use the MPI version
   if (.not.allocated(Eps_hw_temp)) allocate(Eps_hw_temp(16,N), source = 0.0d0) ! all are there
   if (.not.allocated(kappa_temp)) allocate(kappa_temp(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_ee_temp)) allocate(kappa_ee_temp(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_mu_grid_temp)) allocate(kappa_mu_grid_temp(Nsiz_Te), source = 0.0d0)   ! mu [eV]
   if (.not.allocated(kappa_Ce_grid_temp)) allocate(kappa_Ce_grid_temp(Nsiz_Te), source = 0.0d0)   ! Ce
   ! Do sequential do-cycle here, because the parallel regions are inside of hamiltonian calculations
   do Ngp = 1, Nsiz
      ! Split total index into 3 coordinates indices:
      ix = ceiling( dble(Ngp)/dble(iym*izm) )
      iy = ceiling( dble(Ngp - (ix-1)*iym*izm)/dble(izm) )
      iz = Ngp - (ix-1)*iym*izm - (iy-1)*izm

      ! k-points:
      call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid) ! module "TB"

      if (numpar%MPI_param%process_rank == 0) then ! only master process does that
         if (numpar%verbose) write(*,'(a,i4,a,i6,i3,i3,i3,f9.4,f9.4,f9.4,a)') '[MPI Process #', numpar%MPI_param%process_rank, &
                                    '] point #', Ngp, ix, iy, iz, kx, ky, kz, ' Kubo-Greenwood'
      endif

      ! Get the effective momentum and kinetic-energy-related operators:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
         select type(ARRAY)
         type is (TB_H_Pettifor) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_Molteni)  ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_Fu) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_NRL)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         type is (TB_H_DFTB)  ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         type is (TB_H_3TB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         type is (TB_H_xTB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         end select
      END ASSOCIATE

      !-------------------------------
      ! Get the parameters of the CDF:
      if ((numpar%optic_model == 4) .or. (numpar%optic_model == 5)) then ! if requested
         call get_Kubo_Greenwood_CDF(numpar, Scell, NSC, w_grid, cPRRx, cPRRy, cPRRz, Ei, Eps_hw_temp)   ! below
      else  ! skip it, if not reqired
         Eps_hw_temp = 0.0d0
      endif
      ! Save data:
      Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points

      !-------------------------------
      ! If required, do Onsager coefficients (for electronic heat conductivity):
      if (numpar%do_kappa) then  ! if requested
         call get_Onsager_coeffs(numpar, matter, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, kappa_temp, &
                              Scell(NSC)%kappa_Te_grid, kappa_mu_grid_temp, kappa_Ce_grid_temp)   ! below
      else  ! if not required
         kappa_temp = 0.0d0
         kappa_mu_grid_temp = 0.0d0
         kappa_Ce_grid_temp = 0.0d0
      endif

      !-----------------------------------------
      ! Add contribution of the electronic term:
      if (numpar%do_kappa) then  ! if requested
         call get_kappa_e_e(numpar, matter, Scell, NSC, Ei, kappa_mu_grid_temp, Scell(NSC)%kappa_Te_grid, kappa_ee_temp) ! below
      else
         kappa_ee_temp = 0.0d0
      endif

      !=========================================
      ! Combine terms:
      kappa = kappa + kappa_temp ! sum up at different k-points
      kappa_mu_grid = kappa_mu_grid + kappa_mu_grid_temp    ! average mu
      kappa_Ce_grid = kappa_Ce_grid + kappa_Ce_grid_temp    ! average Ce
      kappa_ee = kappa_ee + kappa_ee_temp ! sum up at different k-points
   enddo ! Ngp
   if (allocated(Eps_hw_temp)) deallocate(Eps_hw_temp)
   if (allocated(kappa_temp)) deallocate(kappa_temp)
   if (allocated(kappa_mu_grid_temp)) deallocate(kappa_mu_grid_temp)
   if (allocated(kappa_Ce_grid_temp)) deallocate(kappa_Ce_grid_temp)
   if (allocated(kappa_ee_temp)) deallocate(kappa_ee_temp)

#else ! use OpenMP instead
   !$omp PARALLEL private(ix, iy, iz, Ngp, kx, ky, kz, cPRRx, cPRRy, cPRRz, CHij, Ei, Eps_hw_temp, &
   !$omp                  kappa_temp, kappa_ee_temp, kappa_mu_grid_temp, kappa_Ce_grid_temp)
   if (.not.allocated(Eps_hw_temp)) allocate(Eps_hw_temp(16,N), source = 0.0d0) ! all are there
   if (.not.allocated(kappa_temp)) allocate(kappa_temp(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_ee_temp)) allocate(kappa_ee_temp(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_mu_grid_temp)) allocate(kappa_mu_grid_temp(Nsiz_Te), source = 0.0d0)   ! mu [eV]
   if (.not.allocated(kappa_Ce_grid_temp)) allocate(kappa_Ce_grid_temp(Nsiz_Te), source = 0.0d0)   ! Ce
   !$omp do schedule(dynamic) reduction( + : Eps_hw, kappa, kappa_ee, kappa_mu_grid, kappa_Ce_grid)
   do Ngp = 1, Nsiz
      ! Split total index into 3 coordinates indices:
      ix = ceiling( dble(Ngp)/dble(iym*izm) )
      iy = ceiling( dble(Ngp - (ix-1)*iym*izm)/dble(izm) )
      iz = Ngp - (ix-1)*iym*izm - (iy-1)*izm

      ! k-points:
      call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid) ! module "TB"

#ifdef _OPENMP
      if (numpar%verbose) write(*,'(a,i4,a,i6,i3,i3,i3,f9.4,f9.4,f9.4,a)') 'Thread #', OMP_GET_THREAD_NUM(), &
                                     ' point #', Ngp, ix, iy, iz, kx, ky, kz, ' Kubo-Greenwood'
#else
      if (numpar%verbose) write(*,'(a,i7,i3,i3,i3,f9.4,f9.4,f9.4,a)') ' point #', Ngp, ix, iy, iz, kx, ky, kz, ' Kubo-Greenwood'
#endif

      ! Get the effective momentum and kinetic-energy-related operators:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
         select type(ARRAY)
         type is (TB_H_Pettifor) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_Molteni)  ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_Fu) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_NRL)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         type is (TB_H_DFTB)  ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         type is (TB_H_3TB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         type is (TB_H_xTB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij) ! module "TB"
         end select
      END ASSOCIATE

      !-------------------------------
      ! Get the parameters of the CDF:
      if ((numpar%optic_model == 4) .or. (numpar%optic_model == 5)) then ! if requested
         call get_Kubo_Greenwood_CDF(numpar, Scell, NSC, w_grid, cPRRx, cPRRy, cPRRz, Ei, Eps_hw_temp)   ! below
      else  ! skip it, if not reqired
         Eps_hw_temp = 0.0d0
      endif
      ! Save data:
      Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points

      !-------------------------------
      ! If required, do Onsager coefficients (for electronic heat conductivity):
      if (numpar%do_kappa) then  ! if requested
         call get_Onsager_coeffs(numpar, matter, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, kappa_temp, &
                              Scell(NSC)%kappa_Te_grid, kappa_mu_grid_temp, kappa_Ce_grid_temp)   ! below
      else  ! if not required
         kappa_temp = 0.0d0
         kappa_mu_grid_temp = 0.0d0
         kappa_Ce_grid_temp = 0.0d0
      endif

      !-----------------------------------------
      ! Add contribution of the electronic term:
      if (numpar%do_kappa) then  ! if requested
         call get_kappa_e_e(numpar, matter, Scell, NSC, Ei, kappa_mu_grid_temp, Scell(NSC)%kappa_Te_grid, kappa_ee_temp) ! below
      else
         kappa_ee_temp = 0.0d0
      endif

      !=========================================
      ! Combine terms:
      kappa = kappa + kappa_temp ! sum up at different k-points
      kappa_mu_grid = kappa_mu_grid + kappa_mu_grid_temp    ! average mu
      kappa_Ce_grid = kappa_Ce_grid + kappa_Ce_grid_temp    ! average Ce
      kappa_ee = kappa_ee + kappa_ee_temp ! sum up at different k-points
   enddo ! Ngp
   !$omp end do
   if (allocated(Eps_hw_temp)) deallocate(Eps_hw_temp)
   if (allocated(kappa_temp)) deallocate(kappa_temp)
   if (allocated(kappa_mu_grid_temp)) deallocate(kappa_mu_grid_temp)
   if (allocated(kappa_Ce_grid_temp)) deallocate(kappa_Ce_grid_temp)
   if (allocated(kappa_ee_temp)) deallocate(kappa_ee_temp)
   !$omp end parallel
#endif

   ! Save the k-point averages:
   Eps_hw = Eps_hw/dble(Nsiz) ! normalize k-point summation
   if ((numpar%optic_model == 4) .or. (numpar%optic_model == 5)) then ! if requested
      Scell(NSC)%eps%Eps_hw = Eps_hw   ! all data for spectrum in array
      ! Get the values for the single value of the probe pulse:
      call Find_in_array_monoton(w_grid, Scell(NSC)%eps%w, i)  ! module "Little_subroutines"
      ! Use closest value on the grid:
      Scell(NSC)%eps%ReEps = Eps_hw(2,i)  ! real part of CDF
      Scell(NSC)%eps%ImEps = Eps_hw(3,i)  ! imaginary part of CDF
      Scell(NSC)%eps%R = Eps_hw(5,i)   ! reflectivity
      Scell(NSC)%eps%T = Eps_hw(6,i)   ! transmission
      Scell(NSC)%eps%A = Eps_hw(7,i)   ! absorption
      Scell(NSC)%eps%n = Eps_hw(8,i)   ! optical n
      Scell(NSC)%eps%k = Eps_hw(9,i)   ! optical k
      Scell(NSC)%eps%dc_cond = Eps_hw(10,i)  ! dc-conductivity
      Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,i), Eps_hw(12,i))  ! Re_E_xx and Im_E_xx
      Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,i), Eps_hw(14,i))  ! Re_E_yy and Im_E_yy
      Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,i), Eps_hw(16,i))  ! Re_E_zz and Im_E_zz
   endif

   ! Save electron heat conductivity, averaged over k-points:
   if (numpar%do_kappa) then  ! if requested
      Scell(NSC)%kappa_e = kappa(1)/dble(Nsiz)  ! room temperature
      Scell(NSC)%kappa_e_vs_Te = kappa/dble(Nsiz)  ! array vs Te
      Scell(NSC)%kappa_ee_vs_Te = kappa_ee/dble(Nsiz)  ! array vs Te
      Scell(NSC)%kappa_mu_grid = kappa_mu_grid/dble(Nsiz)  ! array of mu vs Te
      Scell(NSC)%kappa_Ce_grid = kappa_Ce_grid/dble(Nsiz)  ! array of Ce vs Te
   endif

   ! Clean up:
   if (allocated(cPRRx)) deallocate(cPRRx)
   if (allocated(cPRRy)) deallocate(cPRRy)
   if (allocated(cPRRz)) deallocate(cPRRz)
   if (allocated(w_grid)) deallocate(w_grid)
   if (allocated(Ei)) deallocate(Ei)
   if (allocated(CHij)) deallocate(CHij)
   if (allocated(Eps_hw)) deallocate(Eps_hw)
   if (allocated(kappa)) deallocate(kappa)
   if (allocated(kappa_mu_grid)) deallocate(kappa_mu_grid)
   if (allocated(kappa_Ce_grid)) deallocate(kappa_Ce_grid)
   if (allocated(kappa_ee)) deallocate(kappa_ee)
end subroutine get_Kubo_Greenwood_all_complex


subroutine get_kappa_e_e(numpar, matter, Scell, NSC, Ev, mu, Te_grid, kappa_ee)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function
   type(Solid), intent(inout) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: Ev   ! [eV] energy levels (molecular orbitals)
   real(8), dimension(:), intent(in) :: mu   ! chemical potential [eV]
   real(8), dimension(:), intent(in) :: Te_grid    ! electronic temperature grid [K]
   real(8), dimension(:), intent(out) :: kappa_ee  ! electron heat conductivity tensor [W/(K*m)] vs Te
   !----------------------------
   real(8), dimension(:,:), allocatable :: dfe_dT_grid
   real(8), dimension(:), allocatable :: dmu, v, A, B, C
   real(8) :: Vol, Ele, L, prefact, temp
   integer :: i, n, Nsiz, N_Te_grid

   if (numpar%do_kappa .or. numpar%do_kappa_dyn) then ! only if requested

      Nsiz = size(Ev)   ! number of energy levels
      N_Te_grid = size(Te_grid)


      if (.not. allocated(dmu)) then
         allocate(dmu(N_Te_grid), source = 0.0d0)
      else
         dmu = 0.0d0   ! to start with
      endif

      allocate(dfe_dT_grid(N_Te_grid, Nsiz), source = 0.0d0)
      allocate(v(Nsiz), source = 0.0d0)
      allocate(A(N_Te_grid), source = 0.0d0)
      allocate(B(N_Te_grid), source = 0.0d0)
      allocate(C(N_Te_grid), source = 0.0d0)


      ! Supercell volume:
      Vol = Scell(NSC)%V*1.0d-30 ! [m^3]

      ! Get dmu:
      dmu(1) = (mu(2) - mu(1)) / ((Te_grid(2)-Te_grid(1))*g_kb_EV)
      do i = 2, N_Te_grid
         dmu(i) = (mu(i) - mu(i-1)) / ((Te_grid(i)-Te_grid(i-1))*g_kb_EV)  !
      enddo

      ! Derivative of fe by Te:
      do i = 1, N_Te_grid
         do n = 1, Nsiz ! all energy points
            ! Derivative df/dT (wrong, should be by dE!)
            !dfe_dT_grid(i,n) = Diff_Fermi_Te(Te_grid(i)*g_kb_EV, mu(i), dmu(i), Ev(n))   ! module "Electron_tools"
            ! Derivative df/dE * 1/Te
            dfe_dT_grid(i,n) = 1.0d0/(Te_grid(i)*g_kb_EV) * Diff_Fermi_E(Te_grid(i)*g_kb_EV, mu(i), Ev(n))   ! module "Electron_tools"
         enddo
      enddo

      do n = 1, Nsiz ! all energy points
         ! Count energy from the bottom of VB (CB for metals); assume free-electron mass:
         Ele = Ev(n) - Ev(1)
         ! Alternatively, count from the fermi-level (produces wrong results!)
         !if (Ev(n) < Scell(NSC)%E_VB_top) then ! formally, it's VB
         !   Ele = abs(Ev(n) - Scell(NSC)%E_VB_top)  ! electron energy from fermi energy [eV]
         !else  ! formally, it's CB:
         !   Ele = abs(Ev(n) - Scell(NSC)%E_bottom)  ! electron energy from fermi energy [eV]
         !endif
         v(n) = velosity_from_kinetic_energy(Ele, g_me, afs=.false.)    ! [m/s] module "MC_cross_sections"
      enddo


      do i = 1, N_Te_grid  ! all temperatures

         ! update cross-sections for different temrpeatures:
         call update_cross_section(Scell(NSC), matter, Te_in=Te_grid(i))  ! module "Electron_tools"

         do n = 1, Nsiz ! all energy points
            ! Define electron energy to find mean free path:
            if (Ev(n) < Scell(NSC)%E_VB_top) then ! formally, it's VB
               Ele = abs(Ev(n) - Scell(NSC)%E_VB_top)  ! electron energy from fermi energy [eV]
            else  ! formally, it's CB:
               Ele = abs(Ev(n) - Scell(NSC)%E_bottom)  ! electron energy from fermi energy [eV]
            endif
            ! its mean free path:
            call Mean_free_path(Ele, matter%El_MFP_tot, L) ! [A] module "MC_cross_sections"

            if (L < 1.0d6) then ! exclude infinities
               !kappa_ee(i) = kappa_ee(i) + dfe_dT_grid(i,n) * (Ev(n) - mu(i))**2 * v(n) * L
               temp = dfe_dT_grid(i,n) * v(n) * L
               A(i) = A(i) + temp * (Ev(n) - mu(i))**2
               C(i) = C(i) + temp * (Ev(n) - mu(i))
               B(i) = B(i) + temp
            endif
            !print*, i, n, kappa_ee(i), dfe_dT_grid(i,n), (Ev(n) - mu(i)), v(n), L
         enddo

         kappa_ee(i) = (A(i) - C(i)**2/B(i))
         !pause 'get_kappa_e_e'
      enddo ! i

      ! Include porefactors:
      prefact = 1.0d0/(3.0d0 * Vol) * g_e/g_kb * 1.0d-10
      kappa_ee(:) = prefact * abs(kappa_ee(:))  ! -> [W/(K*m)]

      ! Restore the cross-section:
      call update_cross_section(Scell(NSC), matter)  ! module "Electron_tools"

      ! Clean up:
      if (allocated(dfe_dT_grid)) deallocate(dfe_dT_grid)
      if (allocated(dmu)) deallocate(dmu)
      if (allocated(v)) deallocate(v)
   else
      kappa_ee = 0.0d0
   endif
end subroutine get_kappa_e_e




subroutine get_Onsager_coeffs(numpar, matter, Scell, NSC, cPRRx, cPRRy, cPRRz, Ev, kappa, Te_grid, mu_grid, Ce_grid)
   ! Following Ref. [6] for evaluation of Onsager coefficients, Eqs.(6-7) (assuming w->0)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function
   type(Solid), intent(in) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !complex(8), dimension(:,:), intent(in) :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   complex, dimension(:,:), intent(in) :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   real(8), dimension(:), intent(in) :: Ev   ! [eV] energy levels (molecular orbitals)
   real(8), dimension(:), intent(out) :: kappa ! electron heat conductivity tensor [W/(K*m)] vs Te
   real(8), dimension(:), intent(in) :: Te_grid ! electronic temperature grid [K]
   real(8), dimension(:), intent(out) :: mu_grid ! chemical potential [eV]
   real(8), dimension(:), intent(out) :: Ce_grid   ! heat capacity [J/(m^3 K)]
   !----------------------------
   real(8), dimension(:,:), allocatable :: fe_on_Te_grid, dfe_dE_grid
   real(8) :: prec, Vol, pref, prefact, delta, Eij, eta, P2, A_cur, B_cur, C_cur, E_mu, f_nm, f_delt
   real(8) :: kappa_e, n_s, v_F, pref_ke, ne
   real(8), dimension(:), allocatable :: A, B, C, mu, Ce
   integer :: Nsiz, n, m, N_Te_grid, i

   Nsiz = size(Ev)   ! number of energy levels
   N_Te_grid = size(Te_grid)

   ! Define the distributions on the temperature grid:
   if (.not. allocated(mu)) then
      allocate(mu(N_Te_grid), source = 0.0d0)
   else
      mu = 0.0d0
   endif

   if (.not. allocated(Ce)) then
      allocate(Ce(N_Te_grid), source = 0.0d0)
   else
      Ce = 0.0d0
   endif

   if (.not. allocated(fe_on_Te_grid)) then
      allocate(fe_on_Te_grid(N_Te_grid, Nsiz), source = 0.0d0)
   else
      fe_on_Te_grid = 0.0d0   ! to start with
   endif

   if (.not. allocated(dfe_dE_grid)) then
      allocate(dfe_dE_grid(N_Te_grid, Nsiz), source = 0.0d0)
   else
      dfe_dE_grid = 0.0d0   ! to start with
   endif

   do i = 1, N_Te_grid
      !call find_mu_from_N_T( Ev, dble(Scell(NSC)%Ne), mu(i), Te_grid(i)*g_kb_EV) ! module "Electron_tools"
      ! Instead of setting population in k-points, assume constant electornic 'sea' level across entire k-landscape (mu=mu(k=0)):
      call get_Ce_and_mu(Scell, NSC, Te_grid(i), Ev, Ce(i), mu(i), .true.)   ! module "Electron_tools"
      call set_Fermi( Ev, Te_grid(i)*g_kb_EV, mu(i), fe_on_Te_grid(i,:) )   ! module "Electron_tools"

      ! Save chem.potential:
      mu_grid(i) = mu(i)
      ! Save heat capacity:
      Ce_grid(i) = Ce(i)
      ! Derivative of fe by Te:
      do n = 1, Nsiz ! all energy points
         dfe_dE_grid(i,n) = Diff_Fermi_E(Te_grid(i)*g_kb_EV, mu(i), Ev(n))   ! module "Electron_tools"
      enddo
   enddo

   prec = 1.0d-10 ! [eV] acceptance for degenerate levels
   eta = m_gamm * m_e_h ! finite width of delta function [eV]
   ! Supercell volume:
   Vol = Scell(NSC)%V*1.0d-30 ! [m^3]

   ! Calculate only if requested:
   if (numpar%do_kappa) then ! only if requested

      if (.not. allocated(A)) allocate(A(N_Te_grid)) ! to start with
      A = 0.0d0

      if (.not. allocated(B)) allocate(B(N_Te_grid)) ! to start with
      B = 0.0d0

      if (.not. allocated(C)) allocate(C(N_Te_grid)) ! to start with
      C = 0.0d0

      ! Get the Onsager coefficients:
      call get_Onsager_ABC(numpar, Ev, cPRRx, cPRRy, cPRRz, eta, mu, fe_on_Te_grid, dfe_dE_grid, A, B, C, model=numpar%kappa_model)   ! below

      ! Collect terms to calculate thermal conductivity:
      !prefact = g_Pi * g_h / (g_me**2 * Vol * Scell(NSC)%Te) ! prefactor in L22
      pref = g_Pi * g_h / (g_me**2 * Vol) ! prefactor in L22
      do i = 1, N_Te_grid
         prefact = pref / Te_grid(i)
         if (abs(B(i)) > 0.0d0) then
            kappa(i) = prefact * (A(i) - C(i)**2/B(i))   ! [W/(K*m)]
         else
            kappa(i) = 0.0d0
         endif
      enddo ! i

   endif ! numpar%do_kappa

   if (allocated(fe_on_Te_grid)) deallocate(fe_on_Te_grid)
   if (allocated(A)) deallocate(A)
   if (allocated(B)) deallocate(B)
   if (allocated(C)) deallocate(C)
   if (allocated(mu)) deallocate(mu)
end subroutine get_Onsager_coeffs


subroutine get_Onsager_ABC(numpar, Ev, cPRRx, cPRRy, cPRRz, eta, mu, fe_on_Te_grid, dfe_dE_grid, A, B, C, model)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   complex, dimension(:,:), intent(in) :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   real(8), dimension(:), intent(in) :: Ev   ! [eV] energy levels (molecular orbitals)
   real(8), intent(in) :: eta
   real(8), dimension(:), intent(in) :: mu
   real(8), dimension(:,:), intent(in) :: fe_on_Te_grid, dfe_dE_grid
   real(8), dimension(:), intent(inout) :: A, B, C
   integer, intent(in) :: model
   !---------------------
   real(8) :: prec, delta, Eij, P2, E_mu, A_cur, B_cur, C_cur, f_nm, f_delt
   integer :: n, m, i, N_Te_grid, Nsiz
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   prec = 1.0d-10 ! [eV] acceptance for degenerate levels
   N_Te_grid = size(mu)
   Nsiz = size(Ev)   ! number of energy levels


   select case (model)
   !-----------------------
   case default   ! get numerical df/dE
#ifdef MPI_USED   ! use MPI
      N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
      Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
      Nend = Nsiz
      ! Do the cycle (parallel) calculations:
      do n = Nstart, Nend, N_incr  ! each process does its own part
      !do n = 1, Nsiz ! all energy points
         do m = 1, Nsiz
            Eij = (Ev(m) - Ev(n))   ! [eV]

            ! Get approximate delta-function:
            delta = numerical_delta(Eij, eta)   ! [1/eV] module "Algebra_tools"

            if ( (n /= m) .and. (abs(Eij) > prec) .and. (abs(delta) > prec) ) then   ! nondegenerate levels
               ! Average momentum operator:
               P2 = ( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
                     dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
                     dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0
!                P2 = ( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) + aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
!                      dble(cPRRy(n,m)) * dble(cPRRy(m,n)) + aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
!                      dble(cPRRz(n,m)) * dble(cPRRz(m,n)) + aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0

               ! Collect terms (without prefactors)
               ! for a set of electronic temperatures:
               do i = 1, N_Te_grid
                  !f_nm = Scell(NSC)%fe(n) - Scell(NSC)%fe(m)
                  !E_mu = ( (Ev(n) + Ev(m))*0.5d0 - Scell(NSC)%mu )   ! [eV]

                  E_mu = ( (Ev(n) + Ev(m))*0.5d0 - mu(i) )   ! [eV]

                  f_nm = fe_on_Te_grid(i,n) - fe_on_Te_grid(i,m)
                  f_delt = f_nm * delta

                  B_cur = abs(P2) * f_delt / Eij
                  C_cur = B_cur * E_mu
                  A_cur = C_cur * E_mu

                  ! Sum up the terms:
                  A(i) = A(i) + A_cur
                  B(i) = B(i) + B_cur
                  C(i) = C(i) + C_cur
               enddo ! i
            endif ! (n /= m)

         enddo ! m
      enddo ! n
      error_part = 'Error in get_Onsager_ABC:case 0:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {A}', A) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {B}', B) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {C}', C) ! module "MPI_subroutines"

#else ! use OpenMP instead
      !$omp PARALLEL private(n, m, Eij, delta, P2, E_mu, i, A_cur, B_cur, C_cur, f_nm, f_delt)
      !$omp do schedule(dynamic)  reduction( + : A, B, C)
      do n = 1, Nsiz ! all energy points
         do m = 1, Nsiz
            Eij = (Ev(m) - Ev(n))   ! [eV]

            ! Get approximate delta-function:
            delta = numerical_delta(Eij, eta)   ! [1/eV] module "Algebra_tools"

            if ( (n /= m) .and. (abs(Eij) > prec) .and. (abs(delta) > prec) ) then   ! nondegenerate levels
               ! Average momentum operator:
               P2 = ( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
                     dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
                     dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0
!                P2 = ( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) + aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
!                      dble(cPRRy(n,m)) * dble(cPRRy(m,n)) + aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
!                      dble(cPRRz(n,m)) * dble(cPRRz(m,n)) + aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0

               ! Collect terms (without prefactors)
               ! for a set of electronic temperatures:
               do i = 1, N_Te_grid
                  !f_nm = Scell(NSC)%fe(n) - Scell(NSC)%fe(m)
                  !E_mu = ( (Ev(n) + Ev(m))*0.5d0 - Scell(NSC)%mu )   ! [eV]

                  E_mu = ( (Ev(n) + Ev(m))*0.5d0 - mu(i) )   ! [eV]

                  f_nm = fe_on_Te_grid(i,n) - fe_on_Te_grid(i,m)
                  f_delt = f_nm * delta

                  B_cur = abs(P2) * f_delt / Eij
                  C_cur = B_cur * E_mu
                  A_cur = C_cur * E_mu

                  ! Sum up the terms:
                  A(i) = A(i) + A_cur
                  B(i) = B(i) + B_cur
                  C(i) = C(i) + C_cur
               enddo ! i
            endif ! (n /= m)

         enddo ! m
      enddo ! n
      !$omp end do
      !$omp end parallel
#endif

   !-----------------------
   case (1) ! get analytical df/dE
#ifdef MPI_USED   ! use MPI
      N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
      Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
      Nend = Nsiz
      ! Do the cycle (parallel) calculations:
      do n = Nstart, Nend, N_incr  ! each process does its own part
      !do n = 1, Nsiz ! all energy points
         ! Get summ of P2:
         P2 = 0.0d0
         do m = 1, Nsiz
            if ( (n /= m) ) then   ! nondegenerate levels
               ! Average momentum operator:
               P2 = P2 + ( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
                     dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
                     dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0
!                P2 = P2 +( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) + aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
!                      dble(cPRRy(n,m)) * dble(cPRRy(m,n)) + aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
!                      dble(cPRRz(n,m)) * dble(cPRRz(m,n)) + aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0
            endif ! (n /= m)
         enddo ! m

         ! Collect terms (without prefactors)
         ! for a set of electronic temperatures:
         do i = 1, N_Te_grid
            E_mu = Ev(n) - mu(i) ! [eV]
            f_nm = -dfe_dE_grid(i,n)   ! [1/eV]

            B_cur = abs(P2) * f_nm
            C_cur = B_cur * E_mu
            A_cur = C_cur * E_mu

            ! Sum up the terms:
            A(i) = A(i) + A_cur
            B(i) = B(i) + B_cur
            C(i) = C(i) + C_cur

!             print*, n, i, P2, f_nm
         enddo ! i
      enddo ! n
      error_part = 'Error in get_Onsager_ABC:case 1:'
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {A}', A) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {B}', B) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {C}', C) ! module "MPI_subroutines"

#else ! use OpenMP instead
      !$omp PARALLEL private(n, m, P2, E_mu, i, A_cur, B_cur, C_cur, f_nm)
      !$omp do schedule(dynamic)  reduction( + : A, B, C)
      do n = 1, Nsiz ! all energy points
         ! Get summ of P2:
         P2 = 0.0d0
         do m = 1, Nsiz
            if ( (n /= m) ) then   ! nondegenerate levels
               ! Average momentum operator:
               P2 = P2 + ( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
                     dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
                     dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0
!                P2 = P2 +( dble(cPRRx(n,m)) * dble(cPRRx(m,n)) + aimag(cPRRx(n,m)) * aimag(cPRRx(m,n)) + &
!                      dble(cPRRy(n,m)) * dble(cPRRy(m,n)) + aimag(cPRRy(n,m)) * aimag(cPRRy(m,n)) + &
!                      dble(cPRRz(n,m)) * dble(cPRRz(m,n)) + aimag(cPRRz(n,m)) * aimag(cPRRz(m,n)) ) / 3.0d0
            endif ! (n /= m)
         enddo ! m

         ! Collect terms (without prefactors)
         ! for a set of electronic temperatures:
         do i = 1, N_Te_grid
            E_mu = Ev(n) - mu(i) ! [eV]
            f_nm = -dfe_dE_grid(i,n)   ! [1/eV]

            B_cur = abs(P2) * f_nm
            C_cur = B_cur * E_mu
            A_cur = C_cur * E_mu

            ! Sum up the terms:
            A(i) = A(i) + A_cur
            B(i) = B(i) + B_cur
            C(i) = C(i) + C_cur

!             print*, n, i, P2, f_nm
         enddo ! i
      enddo ! n
      !$omp end do
      !$omp end parallel
#endif
   end select
end subroutine get_Onsager_ABC



subroutine get_Onsager_dynamic(numpar, matter, Scell, NSC, kappa)
   ! Following Ref. [6] for evaluation of Onsager coefficients, Eqs.(6-7) (assuming w->0)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function
   type(Solid), intent(in) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(out) :: kappa ! electron heat conductivity tensor [W/(K*m)] vs Te
   !----------------------------
   real(8) :: Vol, prefact, eta
   real(8) :: A, B, C
   integer :: i

   kappa(:) = 0.0d0  ! to start with

   eta = m_gamm * m_e_h ! finite width of delta function [eV]
   ! Supercell volume:
   Vol = Scell(NSC)%V*1.0d-30 ! [m^3]

   ! Calculate only if requested:
   if (numpar%do_kappa_dyn) then ! only if requested
      A = 0.0d0
      B = 0.0d0
      C = 0.0d0

      ! Get the Onsager coefficients:
      call get_Onsager_ABC_dynamic(numpar, Scell, NSC, Scell(NSC)%Mij, eta, A, B, C) ! below

      ! Collect terms to calculate thermal conductivity:
      prefact = g_Pi * g_h / (g_me**2 * Vol * Scell(NSC)%Te) ! prefactor in L22
      if (abs(B) > 0.0d0) then
         kappa(1) = prefact * (A - C**2/B)   ! [W/(K*m)]
      endif
      !print*, 'get_Onsager_dynamic', kappa(1), A, B, C
   endif
end subroutine get_Onsager_dynamic


subroutine get_Onsager_ABC_dynamic(numpar, Scell, NSC, Mij, eta, A, B, C)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), allocatable, intent(inout) :: Mij ! matrix element for electron-ion coupling
   real(8), intent(in) :: eta ! [eV]
   real(8), intent(inout) :: A, B, C
   !---------------------
   real(8) :: prec, delta, Eij, P2, E_mu, A_cur, B_cur, C_cur, f_nm, f_delt, coef, coef_inv, dt_small, vn, vm
   integer :: n, m, i, j, Nsiz, N_orb, N_at
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   prec = 1.0d-10 ! [eV] acceptance for degenerate levels
   Nsiz = size(Scell(NSC)%Ei)   ! number of energy levels
!    coef = 2.0d0*g_e/g_h
!    coef_inv = numpar%M2_scaling*g_h/(g_e)   ! [s] with scaling factor added (e.g. 2 from d|a|^2/dt)
   dt_small = numpar%dt*(1d-15)  ! [fs]
   N_at = size(Scell(NSC)%MDatoms)  ! number of atoms
   N_orb = Nsiz / N_at  ! number of orbitals per atom

   ! Calculate Tully's nonadiabatic-coupling matrix element:
   call get_coupling_matrix_elements(1, Scell, NSC, 1, 1, 1, 1, 1, 1, Mij)  ! module "TB"

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Nsiz
   ! Do the cycle (parallel) calculations:
   do n = Nstart, Nend, N_incr  ! each process does its own part
   !do n = 1, Nsiz ! all energy points
      i = 1 + (n-1)/N_orb  ! atomic index
      !vn = sqrt( SUM(Scell(NSC)%MDatoms(i)%V(:)*Scell(NSC)%MDatoms(i)%V(:)) ) * 1.0d5  ! [(m/s)]

      do m = 1, Nsiz
         j = 1 + (m-1)/N_orb  ! atomic index
         !vm = sqrt( SUM(Scell(NSC)%MDatoms(j)%V(:)*Scell(NSC)%MDatoms(j)%V(:)) ) * 1.0d5  ! [(m/s)]

         !Eij = (Scell(NSC)%Ei(m) - Scell(NSC)%Ei(n))   ! [eV]
         Eij = (Scell(NSC)%Ei(n) - Scell(NSC)%Ei(m))   ! [eV] test

         ! Get approximate delta-function:
         delta = numerical_delta(Eij, eta)   ! [1/eV] module "Algebra_tools"

         if ( (n /= m) .and. (abs(Eij) > prec) .and. (abs(delta) > prec) ) then   ! nondegenerate levels
            ! Average momentum operator:
            !call get_nonadiabatic_Pij(n, m, Mij, dt_small, v, P2) ! module "Nonadiabatic"
            call get_nonadiabatic_Pij(n, m, Mij, dt_small, Scell(NSC)%MDatoms(i)%V*1.0d5, Scell(NSC)%MDatoms(j)%V*1.0d5, P2) ! module "Nonadiabatic"
            P2 = abs(P2)   ! [ (kg*m/s) ^2 ]
            !if (P2 > 0.0d0) print*, 'p', n, m, sqrt(P2)

            ! Collect terms (without prefactors)
            !f_nm = Scell(NSC)%fe(n) - Scell(NSC)%fe(m)
            f_nm = Scell(NSC)%fe(m) - Scell(NSC)%fe(n)   ! test
            E_mu = ( (Scell(NSC)%Ei(n) + Scell(NSC)%Ei(m))*0.5d0 - Scell(NSC)%mu )   ! [eV]
            f_delt = f_nm * delta

            B_cur = abs(P2) * f_delt / Eij
            C_cur = B_cur * E_mu
            A_cur = C_cur * E_mu

            ! Sum up the terms:
            A = A + A_cur
            B = B + B_cur
            C = C + C_cur
         endif ! (n /= m)
         !print*, 'dyn', n, m, A, B, C, P2, f_delt
      enddo ! m
   enddo ! n
   error_part = 'Error in get_Onsager_ABC_dynamic:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {A}', A) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {B}', B) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {C}', C) ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp PARALLEL private(n, m, Eij, delta, P2, E_mu, A_cur, B_cur, C_cur, f_nm, f_delt, i, j)
   !$omp do schedule(dynamic)  reduction( + : A, B, C)
   do n = 1, Nsiz ! all energy points
      i = 1 + (n-1)/N_orb  ! atomic index
      !vn = sqrt( SUM(Scell(NSC)%MDatoms(i)%V(:)*Scell(NSC)%MDatoms(i)%V(:)) ) * 1.0d5  ! [(m/s)]

      do m = 1, Nsiz
         j = 1 + (m-1)/N_orb  ! atomic index
         !vm = sqrt( SUM(Scell(NSC)%MDatoms(j)%V(:)*Scell(NSC)%MDatoms(j)%V(:)) ) * 1.0d5  ! [(m/s)]

         !Eij = (Scell(NSC)%Ei(m) - Scell(NSC)%Ei(n))   ! [eV]
         Eij = (Scell(NSC)%Ei(n) - Scell(NSC)%Ei(m))   ! [eV] test

         ! Get approximate delta-function:
         delta = numerical_delta(Eij, eta)   ! [1/eV] module "Algebra_tools"

         if ( (n /= m) .and. (abs(Eij) > prec) .and. (abs(delta) > prec) ) then   ! nondegenerate levels
            ! Average momentum operator:
            !call get_nonadiabatic_Pij(n, m, Mij, dt_small, v, P2) ! module "Nonadiabatic"
            call get_nonadiabatic_Pij(n, m, Mij, dt_small, Scell(NSC)%MDatoms(i)%V*1.0d5, Scell(NSC)%MDatoms(j)%V*1.0d5, P2) ! module "Nonadiabatic"
            P2 = abs(P2)   ! [ (kg*m/s) ^2 ]
            !if (P2 > 0.0d0) print*, 'p', n, m, sqrt(P2)

            ! Collect terms (without prefactors)
            !f_nm = Scell(NSC)%fe(n) - Scell(NSC)%fe(m)
            f_nm = Scell(NSC)%fe(m) - Scell(NSC)%fe(n)   ! test
            E_mu = ( (Scell(NSC)%Ei(n) + Scell(NSC)%Ei(m))*0.5d0 - Scell(NSC)%mu )   ! [eV]
            f_delt = f_nm * delta

            B_cur = abs(P2) * f_delt / Eij
            C_cur = B_cur * E_mu
            A_cur = C_cur * E_mu

            ! Sum up the terms:
            A = A + A_cur
            B = B + B_cur
            C = C + C_cur
         endif ! (n /= m)
         !print*, 'dyn', n, m, A, B, C, P2, f_delt
      enddo ! m
   enddo ! n
   !$omp end do
   !$omp end parallel
#endif
end subroutine get_Onsager_ABC_dynamic


subroutine get_Kubo_Greenwood_CDF(numpar, Scell, NSC, w_grid, cPRRx_in, cPRRy_in, cPRRz_in, Ev, Eps_hw_temp)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: w_grid ! frequency grid [1/s]
   !complex(8), dimension(:,:), intent(in) :: cPRRx_in, cPRRy_in, cPRRz_in  ! effective momentum operators
   complex, dimension(:,:), intent(in) :: cPRRx_in, cPRRy_in, cPRRz_in  ! effective momentum operators
   real(8), dimension(:), intent(in) :: Ev   ! [eV] energy levels (molecular orbitals)
   real(8), dimension(:,:), intent(inout) :: Eps_hw_temp ! all CDF data
   !----------------------------
   real(8), dimension(3,3) :: Re_eps_ij, Im_eps_ij, Re_small, Re_mid, Re_large, Im_small, Im_mid, Im_large
   real(8), dimension(:), allocatable :: fe_temp
   real(8), dimension(:,:), allocatable :: f_nm_w_nm
   real(8), dimension(:,:,:,:), allocatable :: A_sigma, B_sigma
   !complex(8), dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   !complex(8) :: Eps_xx, Eps_yy, Eps_zz   ! diagonal components of the complex dielectric tensor
   complex, dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   complex :: Eps_xx, Eps_yy, Eps_zz   ! diagonal components of the complex dielectric tensor
   real(8) :: Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, temp_1, temp_2, temp_3
   real(8) :: w, Vol, prefact, w_mn, g_sigma, w_sigma, denom, prec
   real(8) :: pxpx, pxpy, pxpz, pypx, pypy, pypz, pzpx, pzpy, pzpz
   integer :: i, j, Nsiz, m, n, N_w
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   prec = 1.0d-10 ! [eV] acceptance for degenerate levels
   Nsiz = size(Ev)   ! number of energy levels
   N_w = size(w_grid)   ! number of frequency grid points

   allocate(cPRRx(Nsiz,Nsiz))
   allocate(cPRRy(Nsiz,Nsiz))
   allocate(cPRRz(Nsiz,Nsiz))

   ! Supercell volume:
   Vol = Scell(NSC)%V*1.0d-30 ! [m^3]
   prefact = m_prefac*g_ke/Vol  ! prefactor of sigma

   cPRRx = sqrt(prefact)*cPRRx_in
   cPRRy = sqrt(prefact)*cPRRy_in
   cPRRz = sqrt(prefact)*cPRRz_in

   allocate(f_nm_w_nm(Nsiz,Nsiz), source = 0.0d0)
   allocate(A_sigma(Nsiz,Nsiz,3,3), source = 0.0d0)
   allocate(B_sigma(Nsiz,Nsiz,3,3), source = 0.0d0)
   Re_eps_ij = 0.0d0 ! to start with
   Im_eps_ij = 0.0d0 ! to start with


   ! Get the Fermi function, assuming the same Te and mu for all k-points
   ! #Note: it is extremely important to set proper distribution for k-points,
   ! accounting for the fact that the electronic 'sea' should have the same surface
   ! across the entire k-langscape. Setting wrong distribution affects the CDF drammatically!
   allocate(fe_temp(Nsiz), source = 0.0d0)
   call set_Fermi(Ev, Scell(NSC)%TeeV, Scell(NSC)%mu, fe_temp) ! module "Electron_tools"

   !-------------------
   ! 1) Get frequency-independent terms:
#ifdef MPI_USED   ! use the MPI version
   ! Make sure all the processes are here:
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Nsiz
   ! Do the cycle (parallel) calculations:
   do n = Nstart, Nend, N_incr  ! each process does its own part
   !do n = 1, Nsiz ! all energy points
      do m = 1, Nsiz
         w_mn = (Ev(n) - Ev(m))
         if ( (n /= m) .and. (abs(w_mn) > prec) ) then   ! nondegenerate levels
            w_mn = w_mn/m_e_h   ! [1/s] frequency point
            ! Keeping same distirbution for different k-points  does not work:
            !f_nm_w_nm(n,m) = (Scell(NSC)%fe(n) - Scell(NSC)%fe(m)) / w_mn**2
            ! Instead, use the distirbution with the same electornic surface across k-landscape:
            f_nm_w_nm(n,m) = (fe_temp(n) - fe_temp(m)) / w_mn**2

            select case(numpar%optic_model)
            case (4:5) ! full calculations
               ! (P) * (P) [2]:
               A_sigma(n,m,1,1) = dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n))
               A_sigma(n,m,2,2) = dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n))
               A_sigma(n,m,3,3) = dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n))

               B_sigma(n,m,1,1) = dble(cPRRx(n,m)) * aimag(cPRRx(m,n)) + aimag(cPRRx(n,m)) * dble(cPRRx(m,n))
               B_sigma(n,m,2,2) = dble(cPRRy(n,m)) * aimag(cPRRy(m,n)) + aimag(cPRRy(n,m)) * dble(cPRRy(m,n))
               B_sigma(n,m,3,3) = dble(cPRRz(n,m)) * aimag(cPRRz(m,n)) + aimag(cPRRz(n,m)) * dble(cPRRz(m,n))
            case (-5) ! only real part, to test
               A_sigma(n,m,1,1) = dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n))
               A_sigma(n,m,2,2) = dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n))
               A_sigma(n,m,3,3) = dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n))

               B_sigma(n,m,1,1) = 0.0d0
               B_sigma(n,m,2,2) = 0.0d0
               B_sigma(n,m,3,3) = 0.0d0
            end select
         endif
      enddo
   enddo
   error_part = 'Error in get_Kubo_Greenwood_CDF:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {A_sigma}', A_sigma) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {B_sigma}', B_sigma) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

#else ! nonparallel / OpenMP
!    !$omp PARALLEL private(n, m, w_mn)
!    !$omp do schedule(dynamic)
   do n = 1, Nsiz ! all energy points
      do m = 1, Nsiz
         w_mn = (Ev(n) - Ev(m))
         if ( (n /= m) .and. (abs(w_mn) > prec) ) then   ! nondegenerate levels
            w_mn = w_mn/m_e_h   ! [1/s] frequency point
            ! Keeping same distirbution for different k-points  does not work:
            !f_nm_w_nm(n,m) = (Scell(NSC)%fe(n) - Scell(NSC)%fe(m)) / w_mn**2
            ! Instead, use the distirbution with the same electornic surface across k-landscape:
            f_nm_w_nm(n,m) = (fe_temp(n) - fe_temp(m)) / w_mn**2

            select case(numpar%optic_model)
            case (4:5) ! full calculations
               ! (P) * (P) [2]:
               A_sigma(n,m,1,1) = dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n))
               A_sigma(n,m,2,2) = dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n))
               A_sigma(n,m,3,3) = dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n))

               B_sigma(n,m,1,1) = dble(cPRRx(n,m)) * aimag(cPRRx(m,n)) + aimag(cPRRx(n,m)) * dble(cPRRx(m,n))
               B_sigma(n,m,2,2) = dble(cPRRy(n,m)) * aimag(cPRRy(m,n)) + aimag(cPRRy(n,m)) * dble(cPRRy(m,n))
               B_sigma(n,m,3,3) = dble(cPRRz(n,m)) * aimag(cPRRz(m,n)) + aimag(cPRRz(n,m)) * dble(cPRRz(m,n))

!                A_sigma(n,m,1,1) = dcmplx(dconjg(cPRRx(n,m))) * dcmplx(cPRRx(n,m))
!                A_sigma(n,m,2,2) = dcmplx(dconjg(cPRRy(n,m))) * dcmplx(cPRRy(n,m))
!                A_sigma(n,m,3,3) = dcmplx(dconjg(cPRRz(n,m))) * dcmplx(cPRRz(n,m))
!
!                B_sigma(n,m,1,1) = 0.0d0   !aimag(dconjg(cPRRx(n,m)) * cPRRx(n,m))
!                B_sigma(n,m,2,2) = 0.0d0   !aimag(dconjg(cPRRx(n,m)) * cPRRx(n,m))
!                B_sigma(n,m,3,3) = 0.0d0   !aimag(dconjg(cPRRx(n,m)) * cPRRx(n,m))

            case (-5) ! only real part, to test
               A_sigma(n,m,1,1) = dble(cPRRx(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRx(m,n))
               A_sigma(n,m,2,2) = dble(cPRRy(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRy(m,n))
               A_sigma(n,m,3,3) = dble(cPRRz(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRz(m,n))

               B_sigma(n,m,1,1) = 0.0d0
               B_sigma(n,m,2,2) = 0.0d0
               B_sigma(n,m,3,3) = 0.0d0
            end select

            ! Off-diagonal (currently unused):
            !A_sigma(n,m,1,2) = dble(cPRRx(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRy(m,n))
            !A_sigma(n,m,1,3) = dble(cPRRx(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRx(n,m)) * aimag(cPRRz(m,n))
            !A_sigma(n,m,2,1) = dble(cPRRy(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRx(m,n))
            !A_sigma(n,m,2,3) = dble(cPRRy(n,m)) * dble(cPRRz(m,n)) - aimag(cPRRy(n,m)) * aimag(cPRRz(m,n))
            !A_sigma(n,m,3,1) = dble(cPRRz(n,m)) * dble(cPRRx(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRx(m,n))
            !A_sigma(n,m,3,2) = dble(cPRRz(n,m)) * dble(cPRRy(m,n)) - aimag(cPRRz(n,m)) * aimag(cPRRy(m,n))
            ! Off-diagonal (currently unused):
            !B_sigma(n,m,1,2) = dble(cPRRx(n,m)) * aimag(cPRRy(m,n)) + aimag(cPRRx(n,m)) * dble(cPRRy(m,n))
            !B_sigma(n,m,1,3) = dble(cPRRx(n,m)) * aimag(cPRRz(m,n)) + aimag(cPRRx(n,m)) * dble(cPRRz(m,n))
            !B_sigma(n,m,2,1) = dble(cPRRy(n,m)) * aimag(cPRRx(m,n)) + aimag(cPRRy(n,m)) * dble(cPRRx(m,n))
            !B_sigma(n,m,2,3) = dble(cPRRy(n,m)) * aimag(cPRRz(m,n)) + aimag(cPRRy(n,m)) * dble(cPRRz(m,n))
            !B_sigma(n,m,3,1) = dble(cPRRz(n,m)) * aimag(cPRRx(m,n)) + aimag(cPRRz(n,m)) * dble(cPRRx(m,n))
            !B_sigma(n,m,3,2) = dble(cPRRz(n,m)) * aimag(cPRRy(m,n)) + aimag(cPRRz(n,m)) * dble(cPRRy(m,n))
         endif
      enddo
   enddo
!    !$omp end do
!    !$omp end parallel
#endif

   !-------------------
   ! 2) Frequency-dependent values:
   do i = 1, N_w  ! all frequencies
      w = w_grid(i)  ! frequency grid point [1/s]

      ! Get the real and imaginary parts of CDF [1]:
      Re_eps_ij = 0.0d0 ! to start wirh
      Im_eps_ij = 0.0d0 ! to start wirh
      ! To sum up different parts, avoiding small-large numbers problem:
      Re_small = 0.0d0
      Re_mid = 0.0d0
      Re_large = 0.0d0
      Im_small = 0.0d0
      Im_mid = 0.0d0
      Im_large = 0.0d0

#ifdef MPI_USED   ! use the MPI version
      N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
      Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
      Nend = Nsiz
      ! Do the cycle (parallel) calculations:
      do n = Nstart, Nend, N_incr  ! each process does its own part
      !do n = 1, Nsiz ! all energy points
         do m = 1, Nsiz ! all energy points
            w_mn = (Ev(n) - Ev(m))   ! [eV] frequency point
            if ( (n /= m) .and. (abs(w_mn) > prec) ) then   ! nondegenerate levels
               w_mn = w_mn / m_e_h   ! [1/s] frequency point

               denom = (w_mn + w)**2 + m_gamm**2
               !denom = (-w_mn + w)**2 + m_gamm**2   ! test [M]
               g_sigma = m_gamm / denom
               w_sigma = (w_mn + w) / denom
               !w_sigma = (-w_mn + w) / denom  ! test [M]

               ! Optical conductivity / w:
               ! [ D ]
               Re_eps_ij(1,1) = Re_eps_ij(1,1) + f_nm_w_nm(n,m) * (A_sigma(n,m,1,1) * g_sigma - B_sigma(n,m,1,1) * w_sigma)
               Re_eps_ij(2,2) = Re_eps_ij(2,2) + f_nm_w_nm(n,m) * (A_sigma(n,m,2,2) * g_sigma - B_sigma(n,m,2,2) * w_sigma)
               Re_eps_ij(3,3) = Re_eps_ij(3,3) + f_nm_w_nm(n,m) * (A_sigma(n,m,3,3) * g_sigma - B_sigma(n,m,3,3) * w_sigma)

               Im_eps_ij(1,1) = Im_eps_ij(1,1) + f_nm_w_nm(n,m) * (A_sigma(n,m,1,1) * w_sigma + B_sigma(n,m,1,1) * g_sigma)
               Im_eps_ij(2,2) = Im_eps_ij(2,2) + f_nm_w_nm(n,m) * (A_sigma(n,m,2,2) * w_sigma + B_sigma(n,m,2,2) * g_sigma)
               Im_eps_ij(3,3) = Im_eps_ij(3,3) + f_nm_w_nm(n,m) * (A_sigma(n,m,3,3) * w_sigma + B_sigma(n,m,3,3) * g_sigma)
            endif
         enddo ! m
      enddo ! n
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_eps_ij}', Re_eps_ij) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_eps_ij}', Im_eps_ij) ! module "MPI_subroutines"
      call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

#else ! nonparallel / OpenMP
!       !$omp PARALLEL private(n, m, w_mn, denom, g_sigma, w_sigma)
!       !$omp do schedule(dynamic) reduction( + : Re_eps_ij, Im_eps_ij)
      do n = 1, Nsiz ! all energy points
         do m = 1, Nsiz ! all energy points
            w_mn = (Ev(n) - Ev(m))   ! [eV] frequency point
            if ( (n /= m) .and. (abs(w_mn) > prec) ) then   ! nondegenerate levels
               w_mn = w_mn / m_e_h   ! [1/s] frequency point

               denom = (w_mn + w)**2 + m_gamm**2
               !denom = (-w_mn + w)**2 + m_gamm**2   ! test [M]
               g_sigma = m_gamm / denom
               w_sigma = (w_mn + w) / denom
               !w_sigma = (-w_mn + w) / denom  ! test [M]

               ! Optical conductivity / w:
               ! [ D ]
               Re_eps_ij(1,1) = Re_eps_ij(1,1) + f_nm_w_nm(n,m) * (A_sigma(n,m,1,1) * g_sigma - B_sigma(n,m,1,1) * w_sigma)
               Re_eps_ij(2,2) = Re_eps_ij(2,2) + f_nm_w_nm(n,m) * (A_sigma(n,m,2,2) * g_sigma - B_sigma(n,m,2,2) * w_sigma)
               Re_eps_ij(3,3) = Re_eps_ij(3,3) + f_nm_w_nm(n,m) * (A_sigma(n,m,3,3) * g_sigma - B_sigma(n,m,3,3) * w_sigma)

               Im_eps_ij(1,1) = Im_eps_ij(1,1) + f_nm_w_nm(n,m) * (A_sigma(n,m,1,1) * w_sigma + B_sigma(n,m,1,1) * g_sigma)
               Im_eps_ij(2,2) = Im_eps_ij(2,2) + f_nm_w_nm(n,m) * (A_sigma(n,m,2,2) * w_sigma + B_sigma(n,m,2,2) * g_sigma)
               Im_eps_ij(3,3) = Im_eps_ij(3,3) + f_nm_w_nm(n,m) * (A_sigma(n,m,3,3) * w_sigma + B_sigma(n,m,3,3) * g_sigma)

!                temp_1 = f_nm_w_nm(n,m) * (A_sigma(n,m,1,1) * g_sigma)
!                temp_2 = f_nm_w_nm(n,m) * (A_sigma(n,m,2,2) * g_sigma)
!                temp_3 = f_nm_w_nm(n,m) * (A_sigma(n,m,3,3) * g_sigma)
!                call split_into_orders(temp_1, Re_small(1,1), Re_mid(1,1), Re_large(1,1))   ! below
!                call split_into_orders(temp_2, Re_small(2,2), Re_mid(2,2), Re_large(2,2))   ! below
!                call split_into_orders(temp_3, Re_small(3,3), Re_mid(3,3), Re_large(3,3))   ! below
!
!                temp_1 = f_nm_w_nm(n,m) * (-A_sigma(n,m,1,1) * w_sigma)   ! test [M]
!                temp_2 = f_nm_w_nm(n,m) * (-A_sigma(n,m,2,2) * w_sigma)   ! test [M]
!                temp_3 = f_nm_w_nm(n,m) * (-A_sigma(n,m,3,3) * w_sigma)   ! test [M]
!                call split_into_orders(temp_1, Im_small(1,1), Im_mid(1,1), Im_large(1,1))   ! below
!                call split_into_orders(temp_2, Im_small(2,2), Im_mid(2,2), Im_large(2,2))   ! below
!                call split_into_orders(temp_3, Im_small(3,3), Im_mid(3,3), Im_large(3,3))   ! below
            endif
         enddo ! m
      enddo ! n
!       !$omp end do
!       !$omp end parallel
#endif

      ! Get the prefactors: (Prefactors were already included in cPRR above!)
      !Re_eps_ij = prefact * Re_eps_ij
      !Im_eps_ij = prefact * Im_eps_ij
!       Re_eps_ij = Re_small + Re_mid + Re_large  ! sum app different-orders terms
!       Im_eps_ij = Im_small + Im_mid + Im_large  ! sum app different-orders terms

      ! Convert conductivity into CDF:
      !Eps_xx = dcmplx( 1.0d0 - 4.0d0*g_Pi*Im_eps_ij(1,1) / w, 4.0d0*g_Pi*Re_eps_ij(1,1) / w )
      !Eps_yy = dcmplx( 1.0d0 - 4.0d0*g_Pi*Im_eps_ij(2,2) / w, 4.0d0*g_Pi*Re_eps_ij(2,2) / w )
      !Eps_zz = dcmplx( 1.0d0 - 4.0d0*g_Pi*Im_eps_ij(3,3) / w, 4.0d0*g_Pi*Re_eps_ij(3,3) / w )
      ! Convert (sigma/w) into CDF:
      Eps_xx = dcmplx( 1.0d0 - 4.0d0*g_Pi*Im_eps_ij(1,1), 4.0d0*g_Pi*Re_eps_ij(1,1))
      Eps_yy = dcmplx( 1.0d0 - 4.0d0*g_Pi*Im_eps_ij(2,2), 4.0d0*g_Pi*Re_eps_ij(2,2))
      Eps_zz = dcmplx( 1.0d0 - 4.0d0*g_Pi*Im_eps_ij(3,3), 4.0d0*g_Pi*Re_eps_ij(3,3))

      ! Average CDF:
      Im_eps = aimag(Eps_xx + Eps_yy + Eps_zz) / 3.0d0
      Re_eps = dble (Eps_xx + Eps_yy + Eps_zz) / 3.0d0

      ! DC-conductivity:
      dc_cond = Im_eps*w*g_e0     ! averaged over x, y, and z

      ! Get optical coefficients:
      call get_RTA_from_CDF(Re_eps, Im_eps, Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, R, T, A) ! below
      call get_n_k(Re_eps, Im_eps, opt_n, opt_k)  ! transfer Re(e) and Im(e) into n and k

      ! Write them all into array:
      call save_Eps_hw(Eps_hw_temp, i, w*m_e_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, &
               opt_n, opt_k, Eps_xx, Eps_yy, Eps_zz, dc_cond=dc_cond) ! below

   enddo ! i

   ! Clean up:
   deallocate(f_nm_w_nm, fe_temp)
   deallocate(cPRRx, cPRRy, cPRRz)
end subroutine get_Kubo_Greenwood_CDF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines from Ref.[2]:
subroutine get_Graf_Vogl_all_complex(numpar, Scell, NSC, all_w, Err)  ! From Ref. [2]
   type (Numerics_param), intent(inout) :: numpar  ! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(inout) :: Scell   ! supercell with all the atoms as one object
   integer, intent(in) :: NSC    ! number of supercell
   logical, intent(in) :: all_w  ! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err   ! error save
   !--------------------
   real(8), dimension(:,:,:), allocatable :: m_eff ! [1/me] inverse effective mass (in units of electron mass)
   !complex(8), dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   !complex(8), dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   !complex(8) :: Eps_xx, Eps_yy, Eps_zz   ! diagonal components of the complex dielectric tensor
   complex, dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   complex :: Eps_xx, Eps_yy, Eps_zz   ! diagonal components of the complex dielectric tensor
   real(8), dimension(:,:,:,:), allocatable :: cTnn ! kinetic energy-related [dimensionless]
   real(8), dimension(:,:), allocatable :: Eps
   real(8) :: Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   real(8) :: w, kx, ky, kz
   real(8), dimension(:), allocatable :: w_grid
   integer :: i, j, N, FN, ix, iy, iz, ixm, iym, izm, schem, Ngp, Nsiz, Nx, Ny, Nz
   real(8), dimension(:), allocatable :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable :: Eps_hw ! array of all eps vs hw
   real(8), dimension(:,:), allocatable :: Eps_hw_temp ! array of all eps vs hw
   !--------------------

   ! Allocate the array of optical coefficients spectrum (if not allocated before):
   ! Define the grid size:
   if (Scell(NSC)%eps%all_w) then ! full spectrum:
      call allocate_Eps_hw(Scell(NSC)%eps%E_min, Scell(NSC)%eps%E_max, Scell(NSC)%eps%dE, Scell(NSC)%eps%Eps_hw) ! see below
      N = size(Scell(NSC)%eps%Eps_hw,2) ! use this size to define temporary arrays within this subroutine
      ! Grid of frequencies:
      allocate(w_grid(N),source = 0.0d0)
      w_grid(1) =  Scell(NSC)%eps%E_min*g_e/g_h ! [1/s] frequency starting point for the optical spectrum [eV]
      do i = 2, N
         w_grid(i) = w_grid(i-1) + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
      enddo
   else
      N = 1
      if (.not.allocated(Scell(NSC)%eps%Eps_hw)) then
         allocate(Scell(NSC)%eps%Eps_hw(16,N)) ! different parameters saved for all energy grid points
      endif
   endif

   if (.not.allocated(Eps_hw)) allocate(Eps_hw(16,N), source = 0.0d0) ! all are there

   ! For the user-defind number of k-points:
   ixm = numpar%ixm
   iym = numpar%iym
   izm = numpar%izm
   if (allocated(numpar%k_grid)) then
      schem = 1   ! user-defined grid is present
      Nsiz = size(numpar%k_grid,1)  ! size of the user provided grid
   else
      schem = 0   ! no user-defined grid, use default
      Nsiz = ixm*iym*izm
   endif


#ifdef MPI_USED   ! use the MPI version [tested]
   if (.not.allocated(Eps_hw_temp)) allocate(Eps_hw_temp(16,N), source = 0.0d0) ! all are there
   ! Do sequential do-cycle here, because parallel calculations are inside of Hamiltonian etc.
   do Ngp = 1, Nsiz
      ! Split total index into 3 coordinates indices:
      ix = ceiling( dble(Ngp)/dble(iym*izm) )
      iy = ceiling( dble(Ngp - (ix-1)*iym*izm)/dble(izm) )
      iz = Ngp - (ix-1)*iym*izm - (iy-1)*izm

      ! k-points:
      call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid) ! module "TB"

      if (numpar%MPI_param%process_rank == 0) then ! only master process does that
         if (numpar%verbose) write(*,'(a,i4,a,i7,i3,i3,i3,f9.4,f9.4,f9.4,a)') '[MPI process #', numpar%MPI_param%process_rank, &
                                     '] point #', Ngp, ix, iy, iz, kx, ky, kz, ' Graf-Vogl'
      endif

      ! Get the effective momentum and kinetic-energy-related operators:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
         select type(ARRAY)
         type is (TB_H_Pettifor) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, cTnn=cTnn)  ! module "TB"
         type is (TB_H_Molteni)  ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, cTnn=cTnn)  ! module "TB"
         type is (TB_H_Fu) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, cTnn=cTnn)  ! module "TB"
         type is (TB_H_NRL)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         type is (TB_H_DFTB)  ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         type is (TB_H_3TB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         type is (TB_H_xTB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         end select
      END ASSOCIATE

      ! Get the inverse effective mass (in units of m_e):
      call inv_effective_mass(numpar, cPRRx, cPRRy, cPRRz, cTnn, Ei, m_eff)  ! below

      ! Get the parameters of the CDF:
      if (Scell(NSC)%eps%all_w) then ! full spectrum:
         do i = 1, N
            w = w_grid(i)  ! frequency
            ! Get CDF:
            call get_Graf_Vogl_CDF(numpar, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, m_eff, w, &
               Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, Eps_xx, Eps_yy, Eps_zz)   ! below

            ! Write them all into array:
            call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, &
                  opt_n, opt_k, Eps_xx, Eps_yy, Eps_zz, dc_cond=dc_cond) ! below
         enddo
         Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points

      else ! only for given probe:
         w = Scell(NSC)%eps%w
         ! Get CDF:
         call get_Graf_Vogl_CDF(numpar, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, m_eff, w, &
               Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, Eps_xx, Eps_yy, Eps_zz)   ! below

         ! Write them all into array:
         call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, &
               opt_n, opt_k, Eps_xx, Eps_yy, Eps_zz, dc_cond=dc_cond) ! below
         Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points
      endif ! (Scell(NSC)%eps%all_w)
   enddo ! Ngp
   if (allocated(Eps_hw_temp)) deallocate(Eps_hw_temp)

#else ! use OpenMP instead
   !$omp PARALLEL private(ix, iy, iz, Ngp, kx, ky, kz, cPRRx, cPRRy, cPRRz, cTnn, CHij, Ei, m_eff, w, i, &
   !$omp                  Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, Eps_hw_temp, Eps_xx, Eps_yy, Eps_zz)
   if (.not.allocated(Eps_hw_temp)) allocate(Eps_hw_temp(16,N), source = 0.0d0) ! all are there
   !$omp do schedule(dynamic) reduction( + : Eps_hw)
   do Ngp = 1, Nsiz
      ! Split total index into 3 coordinates indices:
      ix = ceiling( dble(Ngp)/dble(iym*izm) )
      iy = ceiling( dble(Ngp - (ix-1)*iym*izm)/dble(izm) )
      iz = Ngp - (ix-1)*iym*izm - (iy-1)*izm

      ! k-points:
      call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid) ! module "TB"

#ifdef _OPENMP
      if (numpar%verbose) write(*,'(a,i4,a,i7,i3,i3,i3,f9.4,f9.4,f9.4,a)') 'Thread #', OMP_GET_THREAD_NUM(), &
                                     ' point #', Ngp, ix, iy, iz, kx, ky, kz, ' Graf-Vogl'
#else
      if (numpar%verbose) write(*,'(a,i7,i3,i3,i3,f9.4,f9.4,f9.4,a)') ' point #', Ngp, ix, iy, iz, kx, ky, kz, ' Graf-Vogl'
#endif


      ! Get the effective momentum and kinetic-energy-related operators:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
         select type(ARRAY)
         type is (TB_H_Pettifor) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, cTnn=cTnn)  ! module "TB"
         type is (TB_H_Molteni)  ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, cTnn=cTnn)  ! module "TB"
         type is (TB_H_Fu) ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, cTnn=cTnn)  ! module "TB"
         type is (TB_H_NRL)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         type is (TB_H_DFTB)  ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         type is (TB_H_3TB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         type is (TB_H_xTB)   ! nonorthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz, Sij=Scell(NSC)%Sij, cTnn=cTnn) ! module "TB"
         end select
      END ASSOCIATE

      ! Get the inverse effective mass (in units of m_e):
      call inv_effective_mass(numpar, cPRRx, cPRRy, cPRRz, cTnn, Ei, m_eff)  ! below

      ! Get the parameters of the CDF:
      if (Scell(NSC)%eps%all_w) then ! full spectrum:
         do i = 1, N
            w = w_grid(i)  ! frequency
            ! Get CDF:
            call get_Graf_Vogl_CDF(numpar, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, m_eff, w, &
               Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, Eps_xx, Eps_yy, Eps_zz)   ! below

            ! Write them all into array:
            call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, &
                  opt_n, opt_k, Eps_xx, Eps_yy, Eps_zz, dc_cond=dc_cond) ! below
         enddo
         Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points

      else ! only for given probe:
         w = Scell(NSC)%eps%w
         ! Get CDF:
         call get_Graf_Vogl_CDF(numpar, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, m_eff, w, &
               Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, Eps_xx, Eps_yy, Eps_zz)   ! below

         ! Write them all into array:
         call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, &
               opt_n, opt_k, Eps_xx, Eps_yy, Eps_zz, dc_cond=dc_cond) ! below
         Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points
      endif ! (Scell(NSC)%eps%all_w)
   enddo ! Ngp
   !$omp end do
   if (allocated(Eps_hw_temp)) deallocate(Eps_hw_temp)
   !$omp end parallel
#endif

   ! Save the k-point avreages:
   Eps_hw = Eps_hw/dble(Nsiz) ! normalize k-point summation
   Scell(NSC)%eps%Eps_hw = Eps_hw   ! all data for spectrum in array
   ! Get the values for the single value of the probe pulse:
   call Find_in_array_monoton(w_grid, Scell(NSC)%eps%w, i)  ! module "Little_subroutines"
   ! Use closest value on the grid:
   Scell(NSC)%eps%ReEps = Eps_hw(2,i)  ! real part of CDF
   Scell(NSC)%eps%ImEps = Eps_hw(3,i)  ! imaginary part of CDF
   Scell(NSC)%eps%R = Eps_hw(5,i)   ! reflectivity
   Scell(NSC)%eps%T = Eps_hw(6,i)   ! transmission
   Scell(NSC)%eps%A = Eps_hw(7,i)   ! absorption
   Scell(NSC)%eps%n = Eps_hw(8,i)   ! optical n
   Scell(NSC)%eps%k = Eps_hw(9,i)   ! optical k
   Scell(NSC)%eps%dc_cond = Eps_hw(10,i)  ! dc-conductivity
   Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,i), Eps_hw(12,i))  ! Re_E_xx and Im_E_xx
   Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,i), Eps_hw(14,i))  ! Re_E_yy and Im_E_yy
   Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,i), Eps_hw(16,i))  ! Re_E_zz and Im_E_zz

   ! Clean up:
   if (allocated(m_eff)) deallocate(m_eff)
   if (allocated(cPRRx)) deallocate(cPRRx)
   if (allocated(cPRRy)) deallocate(cPRRy)
   if (allocated(cPRRz)) deallocate(cPRRz)
   if (allocated(cTnn)) deallocate(cTnn)
   if (allocated(Eps)) deallocate(Eps)
   if (allocated(w_grid)) deallocate(w_grid)
   if (allocated(Ei)) deallocate(Ei)
   if (allocated(CHij)) deallocate(CHij)
   if (allocated(Eps_hw)) deallocate(Eps_hw)
end subroutine get_Graf_Vogl_all_complex



subroutine get_Graf_Vogl_CDF(numpar, Scell, NSC, cPRRx, cPRRy, cPRRz, Ev, m_eff, w, &
                        Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond, Eps_xx, Eps_yy, Eps_zz)
   type (Numerics_param), intent(in) :: numpar ! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !complex(8), dimension(:,:), intent(in) :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   complex, dimension(:,:), intent(in) :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   real(8), dimension(:), intent(in) :: Ev   ! [eV] energy levels (molecular orbitals)
   real(8), dimension(:,:,:), intent(in) :: m_eff ! [1/me] inverse effective mass (in units of electron mass)
   real(8), intent(in) :: w   ! frequency [1/s]
   real(8), intent(inout) :: Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   !complex(8), intent(inout) :: Eps_xx, Eps_yy, Eps_zz   ! diagonal components of the complex dielectric tensor
   complex, intent(inout) :: Eps_xx, Eps_yy, Eps_zz   ! diagonal components of the complex dielectric tensor
   !----------------------------
   real(8), dimension(3,3) :: Re_eps_ij, Im_eps_ij, Re_eps_ij_term1
   real(8) :: delt, eta, Vol, w_mn, prefact, Re_prefact, temp, term_1_SI, term_2_SI, Im_term_SI, f_P_term, ww, prec
   real(8) :: pxpx, pxpy, pxpz, pypx, pypy, pypz, pzpx, pzpy, pzpz
   integer :: i, j, Nsiz, m, n

   Nsiz = size(Ev)   ! number of energy levels

   ! Precision for the principal value:
   prec = 1.0d-2

   ! Approximate the delta-function with a gaussian, Chapter VII in [4]
   !eta = 0.3d0 ! [eV]
   eta = 3.0d0 * m_gamm*m_e_h  ! [eV] gamma parameter ~0.3 eV

   ! Supercell volume:
   Vol = Scell(NSC)%V*1.0d-30 ! [m^3]
   temp = g_e**2/Vol ! prefactor in Re(Eps)
   Im_term_SI = g_ke             ! to convert Im(Eps) into SI units (Coulomb-law)
   term_1_SI = Im_term_SI/g_me   ! to convert first term in Re(Eps) into SI units
   term_2_SI = term_1_SI/g_me    ! to convert second term in Re(Eps) into SI units

   ! Get the imaginary part of CDF [1]:
   Re_eps_ij = 0.0d0 ! to start wirh
   Im_eps_ij = 0.0d0 ! to start wirh
   Re_eps_ij_term1 = 0.0d0 ! to start wirh

   do n = 1, Nsiz
      ! First term in Re_eps:
      do i = 1, 3
         do j = 1,3
            !Re_eps_ij_term1(i,j) = Re_eps_ij_term1(i,j) - temp/((w)**2) * Scell(NSC)%fe(n) * dble(m_eff(n,i,j)) * term_1_SI
            Re_eps_ij_term1(i,j) = Re_eps_ij_term1(i,j) - Scell(NSC)%fe(n) * m_eff(n,i,j)
         enddo ! j
      enddo ! i
      ! Im_CDF and second term in Re_CDF:
      do m = 1, Nsiz
         if (m /= n) then
            ! Imaginary part:
            w_mn = (Ev(n) - Ev(m))*g_e/g_h ! [1/s] frequency point
            delt = m_inv_sqrt_pi * exp( -((w - w_mn)/(eta/m_e_h))**2 ) / (eta/m_e_h)   ! delta-function approx. [s]
            prefact = -delt * (Scell(NSC)%fe(n) - Scell(NSC)%fe(m))

!             pxpx = conjg(cPRRx(n,m)) * cPRRx(m,n)
!             pxpy = conjg(cPRRx(n,m)) * cPRRy(m,n)
!             pxpz = conjg(cPRRx(n,m)) * cPRRz(m,n)
!             pypx = conjg(cPRRy(n,m)) * cPRRx(m,n)
!             pypy = conjg(cPRRy(n,m)) * cPRRy(m,n)
!             pypz = conjg(cPRRy(n,m)) * cPRRz(m,n)
!             pzpx = conjg(cPRRz(n,m)) * cPRRx(m,n)
!             pzpy = conjg(cPRRz(n,m)) * cPRRy(m,n)
!             pzpz = conjg(cPRRz(n,m)) * cPRRz(m,n)

            pxpx = dble( cPRRx(n,m) * cPRRx(m,n) )
            pxpy = dble( cPRRx(n,m) * cPRRy(m,n) )
            pxpz = dble( cPRRx(n,m) * cPRRz(m,n) )
            pypx = dble( cPRRy(n,m) * cPRRx(m,n) )
            pypy = dble( cPRRy(n,m) * cPRRy(m,n) )
            pypz = dble( cPRRy(n,m) * cPRRz(m,n) )
            pzpx = dble( cPRRz(n,m) * cPRRx(m,n) )
            pzpy = dble( cPRRz(n,m) * cPRRy(m,n) )
            pzpz = dble( cPRRz(n,m) * cPRRz(m,n) )

            Im_eps_ij(1,1) = Im_eps_ij(1,1) + prefact * pxpx
            Im_eps_ij(1,2) = Im_eps_ij(1,2) + prefact * pxpy
            Im_eps_ij(1,3) = Im_eps_ij(1,3) + prefact * pxpz
            Im_eps_ij(2,1) = Im_eps_ij(2,1) + prefact * pypx
            Im_eps_ij(2,2) = Im_eps_ij(2,2) + prefact * pypy
            Im_eps_ij(2,3) = Im_eps_ij(2,3) + prefact * pypz
            Im_eps_ij(3,1) = Im_eps_ij(3,1) + prefact * pzpx
            Im_eps_ij(3,2) = Im_eps_ij(3,2) + prefact * pzpy
            Im_eps_ij(3,3) = Im_eps_ij(3,3) + prefact * pzpz

            ! Real part:
            ww = (w_mn-w)
            if (abs(ww) > prec) then
               Re_prefact = (Scell(NSC)%fe(n) - Scell(NSC)%fe(m))/((w_mn**2)*g_h*ww)
               Re_eps_ij(1,1) = Re_eps_ij(1,1) + Re_prefact * pxpx
               Re_eps_ij(1,2) = Re_eps_ij(1,2) + Re_prefact * pxpy
               Re_eps_ij(1,3) = Re_eps_ij(1,3) + Re_prefact * pxpz
               Re_eps_ij(2,1) = Re_eps_ij(2,1) + Re_prefact * pypx
               Re_eps_ij(2,2) = Re_eps_ij(2,2) + Re_prefact * pypy
               Re_eps_ij(2,3) = Re_eps_ij(2,3) + Re_prefact * pypz
               Re_eps_ij(3,1) = Re_eps_ij(3,1) + Re_prefact * pzpx
               Re_eps_ij(3,2) = Re_eps_ij(3,2) + Re_prefact * pzpy
               Re_eps_ij(3,3) = Re_eps_ij(3,3) + Re_prefact * pzpz
            endif ! (abs(ww) > prec)

         endif ! (m /= n)
      enddo ! j
   enddo ! i

   ! Include prefactors:
   Im_eps_ij = m_prefac/(w*w * Vol) * Im_eps_ij * Im_term_SI   ! -> SI units
   Re_eps_ij_term1 = Re_eps_ij_term1 * temp/((w-eta/m_e_h)**2) * term_1_SI
   Re_eps_ij = Re_eps_ij * temp*term_2_SI + Re_eps_ij_term1   ! combine terms

   ! Save into output variables:
   Eps_xx = dcmplx( 1.0d0 + 4.0d0*g_Pi*Re_eps_ij(1,1),  4.0d0*g_Pi*Im_eps_ij(1,1) )
   Eps_yy = dcmplx( 1.0d0 + 4.0d0*g_Pi*Re_eps_ij(2,2),  4.0d0*g_Pi*Im_eps_ij(2,2) )
   Eps_zz = dcmplx( 1.0d0 + 4.0d0*g_Pi*Re_eps_ij(3,3),  4.0d0*g_Pi*Im_eps_ij(3,3) )

   ! Convert from conductivity to CDF:
   Im_eps = aimag(Eps_xx + Eps_yy + Eps_zz) / 3.0d0
   Re_eps = dble (Eps_xx + Eps_yy + Eps_zz) / 3.0d0

   ! DC-conductivity:
   dc_cond = Im_eps*w*g_e0     ! averaged over x, y, and z

   ! Get optical coefficients:
   call get_RTA_from_CDF(Re_eps, Im_eps, Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, R, T, A) ! below
   call get_n_k(Re_eps, Im_eps, opt_n, opt_k)  ! transfer Re(e) and Im(e) into n and k
end subroutine get_Graf_Vogl_CDF



subroutine inv_effective_mass(numpar, cPRRx, cPRRy, cPRRz, cTnn, Ev, m_eff)  ! Ref.[2], Eq(8)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   complex, dimension(:,:), intent(in) :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   real(8), dimension(:,:,:,:), intent(in) :: cTnn ! kinetic energy-related [dimensionless]
   real(8), dimension(:), intent(in) :: Ev   ! [eV] energy levels (molecular orbitals)
   real(8), dimension(:,:,:), intent(inout), allocatable :: m_eff ! [1/me] inverse effective mass (in units of electron mass)
   !-----------------------------
   integer :: Nsiz, n, m, i, j
   real(8) :: fact, me_inv
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   Nsiz = size(cTnn,1)
   me_inv = 1.0d0/g_me

   if (.not.allocated(m_eff)) allocate(m_eff(Nsiz,3,3))
   m_eff = 0.0d0 ! to start with

#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Nsiz
   ! Do the cycle (parallel) calculations:
   do n = Nstart, Nend, N_incr  ! each process does its own part
   !do n = 1, Nsiz
      do i = 1, 3
         do j = 1, 3
            ! Term 1, kinetic energy (diagonal) contribution:
            m_eff(n,i,j) = cTnn(n,n,i,j)  ! 1/me is excluded, due to units choice
!             print*, n, i, j, cTnn(n,n,i,j)
         enddo ! j
      enddo ! i

      ! Term 2, momentum contribution:
      do m = 1, Nsiz
         if (m /= n) then ! off-diagonal terms:
            fact = me_inv/(Ev(n) - Ev(m)) ! 1/(me*(En-Em))
            m_eff(n,1,1) = m_eff(n,1,1) + fact * ( conjg(cPRRx(n,m))*cPRRx(m,n) + conjg(cPRRx(m,n))*cPRRx(n,m) ) ! i = x, j = x
            m_eff(n,1,2) = m_eff(n,1,2) + fact * ( conjg(cPRRx(n,m))*cPRRy(m,n) + conjg(cPRRx(m,n))*cPRRy(n,m) ) ! i = x, j = y
            m_eff(n,1,3) = m_eff(n,1,3) + fact * ( conjg(cPRRx(n,m))*cPRRz(m,n) + conjg(cPRRx(m,n))*cPRRz(n,m) ) ! i = x, j = z
            m_eff(n,2,1) = m_eff(n,2,1) + fact * ( conjg(cPRRy(n,m))*cPRRx(m,n) + conjg(cPRRy(m,n))*cPRRx(n,m) ) ! i = y, j = x
            m_eff(n,2,2) = m_eff(n,2,2) + fact * ( conjg(cPRRy(n,m))*cPRRy(m,n) + conjg(cPRRy(m,n))*cPRRy(n,m) ) ! i = y, j = y
            m_eff(n,2,3) = m_eff(n,2,3) + fact * ( conjg(cPRRy(n,m))*cPRRz(m,n) + conjg(cPRRy(m,n))*cPRRz(n,m) ) ! i = y, j = z
            m_eff(n,3,1) = m_eff(n,3,1) + fact * ( conjg(cPRRz(n,m))*cPRRx(m,n) + conjg(cPRRz(m,n))*cPRRx(n,m) ) ! i = z, j = x
            m_eff(n,3,2) = m_eff(n,3,2) + fact * ( conjg(cPRRz(n,m))*cPRRy(m,n) + conjg(cPRRz(m,n))*cPRRy(n,m) ) ! i = z, j = y
            m_eff(n,3,3) = m_eff(n,3,3) + fact * ( conjg(cPRRz(n,m))*cPRRz(m,n) + conjg(cPRRz(m,n))*cPRRz(n,m) ) ! i = z, j = z
         endif ! (m /= n)
      enddo ! m
   enddo ! n
   error_part = 'Error in Optical_parameters: inv_effective_mass:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {m_eff}', m_eff) ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp parallel
   !$omp do private(n, m, i, j, fact)
   do n = 1, Nsiz
      do i = 1, 3
         do j = 1, 3
            ! Term 1, kinetic energy (diagonal) contribution:
            m_eff(n,i,j) = cTnn(n,n,i,j)  ! 1/me is excluded, due to units choice
!             print*, n, i, j, cTnn(n,n,i,j)
         enddo ! j
      enddo ! i

      ! Term 2, momentum contribution:
      do m = 1, Nsiz
         if (m /= n) then ! off-diagonal terms:
            fact = me_inv/(Ev(n) - Ev(m)) ! 1/(me*(En-Em))
            m_eff(n,1,1) = m_eff(n,1,1) + fact * ( conjg(cPRRx(n,m))*cPRRx(m,n) + conjg(cPRRx(m,n))*cPRRx(n,m) ) ! i = x, j = x
            m_eff(n,1,2) = m_eff(n,1,2) + fact * ( conjg(cPRRx(n,m))*cPRRy(m,n) + conjg(cPRRx(m,n))*cPRRy(n,m) ) ! i = x, j = y
            m_eff(n,1,3) = m_eff(n,1,3) + fact * ( conjg(cPRRx(n,m))*cPRRz(m,n) + conjg(cPRRx(m,n))*cPRRz(n,m) ) ! i = x, j = z
            m_eff(n,2,1) = m_eff(n,2,1) + fact * ( conjg(cPRRy(n,m))*cPRRx(m,n) + conjg(cPRRy(m,n))*cPRRx(n,m) ) ! i = y, j = x
            m_eff(n,2,2) = m_eff(n,2,2) + fact * ( conjg(cPRRy(n,m))*cPRRy(m,n) + conjg(cPRRy(m,n))*cPRRy(n,m) ) ! i = y, j = y
            m_eff(n,2,3) = m_eff(n,2,3) + fact * ( conjg(cPRRy(n,m))*cPRRz(m,n) + conjg(cPRRy(m,n))*cPRRz(n,m) ) ! i = y, j = z
            m_eff(n,3,1) = m_eff(n,3,1) + fact * ( conjg(cPRRz(n,m))*cPRRx(m,n) + conjg(cPRRz(m,n))*cPRRx(n,m) ) ! i = z, j = x
            m_eff(n,3,2) = m_eff(n,3,2) + fact * ( conjg(cPRRz(n,m))*cPRRy(m,n) + conjg(cPRRz(m,n))*cPRRy(n,m) ) ! i = z, j = y
            m_eff(n,3,3) = m_eff(n,3,3) + fact * ( conjg(cPRRz(n,m))*cPRRz(m,n) + conjg(cPRRz(m,n))*cPRRz(n,m) ) ! i = z, j = z
         endif ! (m /= n)
      enddo ! m

!       print*, n, m_eff(n,:,:)
   enddo ! n
   !$omp end do
   !$omp end parallel
#endif
!    pause 'inv_effective_mass'
end subroutine inv_effective_mass



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complex-value Trani's subroutines from Ref.[1]:

!subroutine get_trani_all_complex(numpar, Scell, NSC, fe, all_w, Err)
subroutine get_trani_all_complex(numpar, Scell, NSC, all_w, Err)  ! From Ref. [2]
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
!    real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err	! error save
   !-----------------------------------------
   real(8), dimension(:,:), allocatable :: Fnnx, Fnny, Fnnz
   real(8), dimension(:,:), allocatable :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), dimension(:,:), allocatable :: Eps
   real(8), dimension(:), allocatable :: w_grid
   real(8) Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   real(8) w, kx, ky, kz
   integer i, j, N, FN, ix, iy, iz, ixm, iym, izm, schem, Ngp, Nsiz
   !integer*4 today(3), now(3)
   real(8), dimension(:), allocatable :: Ei	! energy levels [eV]
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   real(8), dimension(:,:), allocatable :: Eps_hw ! array of all eps vs hw
   real(8), dimension(:,:), allocatable :: Eps_hw_temp ! array of all eps vs hw

!    N = 1000
!    if (.not.allocated(Scell(NSC)%eps%Eps_hw)) then
!       allocate(Scell(NSC)%eps%Eps_hw(9,N)) ! all are there
!    endif
   ! Allocate the array of optical coefficients spectrum (if not allocated before):
   if (all_w) then ! full spectrum:
      call allocate_Eps_hw(Scell(NSC)%eps%E_min, Scell(NSC)%eps%E_max, Scell(NSC)%eps%dE, Scell(NSC)%eps%Eps_hw) ! see below
      N = size(Scell(NSC)%eps%Eps_hw,2) ! use this size to define temporary arrays within this subroutine
   else
      N = 1
      if (.not.allocated(Scell(NSC)%eps%Eps_hw)) then
         allocate(Scell(NSC)%eps%Eps_hw(16,N)) ! different parameters saved for all energy grid points
      endif     
   endif

   if (.not.allocated(Eps_hw)) allocate(Eps_hw(16,N)) ! all are there
   Eps_hw = 0.0d0
   if (.not.allocated(Eps_hw_temp)) then
      allocate(Eps_hw_temp(16,N)) ! all are there
      Eps_hw_temp = 0.0d0
   endif
   allocate(w_grid(N))

   ixm = numpar%ixm
   iym = numpar%iym
   izm = numpar%izm
   if (allocated(numpar%k_grid)) then
      schem = 1	! user-defined grid is present
      Nsiz = size(numpar%k_grid,1)	! size of the user provided grid
   else
      schem = 0	! no user-defined grid
      Nsiz = ixm*iym*izm
   endif

!    !$omp PARALLEL private(ix,iy,iz,i,Eps_hw_temp,CHij,Ei,Fnnx,Fnny,Fnnz,w,Re_eps,Im_eps,R,T,A,opt_n,opt_k) shared(Eps_hw)
!    !$omp do reduction( + : Eps_hw)
   do ix = 1, ixm
      do iy = 1, iym
         do iz = 1, izm
         
            Ngp = (ix-1)*iym*ixm + (iy-1)*ixm + iz	! number of grid point for user defined grid
            if (Ngp > Nsiz) goto 3457
            
!             ! k-points choice from [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]
!             kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
!             ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
!             kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
            call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid)	! module "TB"
            if (numpar%verbose) write(*,'(i3,i3,i3,f9.4,f9.4,f9.4,a)') ix, iy, iz, kx, ky, kz, ' Trani'
            
            call get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
            ! With off-diagonal elements:
!             call get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnxx, Fnnyy, Fnnzz, kx=kx, ky=ky, kz=kz, Err=Err, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy) ! see below
!             call get_Fnn_complex(numpar, matter, atoms, TB_Hamil, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx=real(ix)/real(ixm), ky=real(iy)/real(iym), kz=real(iz)/real(izm), Err=Err) ! see below

            if (all_w) then ! full spectrum:
               w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] frequency starting point for the optical spectrum [eV]
               w_grid(1) = w  ! save for the grid
               do i = 1, N
                  call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
                  ! With off-diagonal elements:
!                   call get_trani(numpar, Scell, NSC, Fnnxx, Fnnyy, Fnnzz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy)
                  !call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
                  call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), &
                  R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy
                  w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
                  if (i<N) w_grid(i+1) = w   ! save for the grid
               enddo
               Eps_hw = Eps_hw + Eps_hw_temp
            else ! only for given probe:
               w = Scell(NSC)%eps%w
               w_grid(1) = w  ! save for the grid
               call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
               ! With off-diagonal elements:
!                call get_trani(numpar, Scell, NSC, Fnnxx, Fnnyy, Fnnzz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy)
               !call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
               call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), &
               R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy
               Eps_hw = Eps_hw + Eps_hw_temp
            endif
!             print*, 'Fnn', Fnnx, Fnny, Fnnz
!             print*, w, Re_eps, Im_eps
         enddo ! iz
      enddo ! iy
   enddo ! ix

!    !$omp end do
!    !$omp end parallel
!    if (all_w) then ! full spectrum:
3457 continue
!       !Eps_hw = Eps_hw/real((ixm+1)*(iym+1)*(izm+1))
!       Eps_hw = Eps_hw/dble(ixm*iym*izm)
     Eps_hw = Eps_hw/dble(Nsiz)
     Scell(NSC)%eps%Eps_hw = Eps_hw     

     ! Get the values for the single value of the probe pulse:
     call Find_in_array_monoton(w_grid, Scell(NSC)%eps%w, i)  ! module "Little_subroutines"

     Scell(NSC)%eps%ReEps = Eps_hw(2,i)	! real part of CDF
     Scell(NSC)%eps%ImEps = Eps_hw(3,i)	! imaginary part of CDF
     Scell(NSC)%eps%R = Eps_hw(5,i)	! reflectivity
     Scell(NSC)%eps%T = Eps_hw(6,i)	! transmission
     Scell(NSC)%eps%A = Eps_hw(7,i)	! absorption
     Scell(NSC)%eps%n = Eps_hw(8,i)	! optical n
     Scell(NSC)%eps%k = Eps_hw(9,i)	! optical k
     Scell(NSC)%eps%dc_cond = Eps_hw(10,i)	! dc-conductivity
     Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,i), Eps_hw(12,i))  ! Re_E_xx and Im_E_xx
     Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,i), Eps_hw(14,i))  ! Re_E_yy and Im_E_yy
     Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,i), Eps_hw(16,i))  ! Re_E_zz and Im_E_zz

!       ! Save parameters for the selected probe-pulse:
!       w = Scell(NSC)%eps%w
!       call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, fe, w, Re_eps, Im_eps, dc_cond=dc_cond)
!    endif
end subroutine get_trani_all_complex



subroutine allocate_Eps_hw(E_min, E_max, dE, Eps_hw)
   real(8), intent(inout) :: E_min, E_max, dE ! min, max, and step for the array of Eps_hw (given in input file)
   real(8), dimension(:,:), allocatable, intent(inout) :: Eps_hw ! array of all eps vs hw
   integer N ! number of grid points in the array Eps_hw to be defined
   if (min(E_min, E_max, dE) <= 0.0d0) then ! parameters are ill-defined, use default instead:
      N = 1000
      E_min = 0.05d0
      E_max = 50.0d0
      dE = 0.05d0
   else ! use user-defined grid
      N = CEILING(ABS(E_max - E_min)/dE) + 1
   endif

   if (.not.allocated(Eps_hw)) then
      allocate(Eps_hw(16,N)) ! different parameters saved for all energy grid points
   endif
!    These 16 parameters are:
!    Eps_hw(1,i) = hw     ! energy
!    Eps_hw(2,i) = Re_eps ! real part of CDF
!    Eps_hw(3,i) = Im_eps ! imaginary part of CDF
!    Eps_hw(4,i) = LF  ! loss function
!    Eps_hw(5,i) = R   ! reflectivity
!    Eps_hw(6,i) = T   ! transmission
!    Eps_hw(7,i) = A   ! absorption
!    Eps_hw(8,i) = n   ! optical n
!    Eps_hw(9,i) = k   ! optical k
!    Eps_hw(10,i) = sigma   ! dc-conductivity
!    Eps_hw(11,i) = Re_E_xx
!    Eps_hw(12,i) = Im_E_xx
!    Eps_hw(13,i) = Re_E_yy
!    Eps_hw(14,i) = Im_E_yy
!    Eps_hw(15,i) = Re_E_zz
!    Eps_hw(16,i) = Im_E_zz
end subroutine allocate_Eps_hw


subroutine get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx ,ky, kz, Err, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy) ! Ref. [2]
!subroutine get_Fnn_complex(numpar, matter, atoms, TB, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx ,ky, kz, Err)
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), allocatable, intent(out) :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable, intent(out) :: Fnnx, Fnny, Fnnz
   real(8), dimension(:,:), allocatable, intent(out), optional :: Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy
   !real(8), dimension(:,:), allocatable, intent(out) :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), intent(in) :: kx ,ky, kz ! k-point to get Hamiltonian at
   type(Error_handling), intent(inout) :: Err	! error save
   !------------------------------------------------------------------------
   ! complex, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   ! complex, dimension(:,:), allocatable :: PBx, PBy, PBz
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the complex hamiltonian
   complex, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   complex, dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz
   complex, dimension(:,:), allocatable :: PBx, PBy, PBz, CSij
   integer i, j, n, nn, FN, n_orb, nat, Nsiz
   character(200) :: Error_descript
   real(8), dimension(:,:), allocatable :: Hij
   real(8), dimension(:), allocatable :: Ei_r, Norm1
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

!    real(8), dimension(:), allocatable :: Ei_test
   
   Error_descript = ''

!    nat = size(Scell(NSC)%MDatoms) ! number of atoms
   Nsiz = size(Scell(NSC)%Ha,1) ! total number of orbitals
!    n_orb =  Nsiz/nat ! orbitals per atom

   if (.not.allocated(CHij)) allocate(CHij(Nsiz,Nsiz))
   if (.not.allocated(Ei)) allocate(Ei(Nsiz))
!    if (.not.allocated(Ei_test)) allocate(Ei_test(Nsiz))
   if (.not.allocated(Fnnx)) allocate(Fnnx(Nsiz,Nsiz))
   if (.not.allocated(Fnny)) allocate(Fnny(Nsiz,Nsiz))
   if (.not.allocated(Fnnz)) allocate(Fnnz(Nsiz,Nsiz))
   if (present(Fnnxy)) then ! off-diagonal
      if (.not.allocated(Fnnxy)) allocate(Fnnxy(Nsiz,Nsiz))
      if (.not.allocated(Fnnxz)) allocate(Fnnxz(Nsiz,Nsiz))
      if (.not.allocated(Fnnyx)) allocate(Fnnyx(Nsiz,Nsiz))
      if (.not.allocated(Fnnyz)) allocate(Fnnyz(Nsiz,Nsiz))
      if (.not.allocated(Fnnzx)) allocate(Fnnzx(Nsiz,Nsiz))
      if (.not.allocated(Fnnzy)) allocate(Fnnzy(Nsiz,Nsiz))
   endif
   
!    if (.not.allocated(Hij)) allocate(Hij(Nsiz,Nsiz))
!    if (.not.allocated(Ei_r)) allocate(Ei_r(Nsiz))
    if (.not.allocated(Norm1)) allocate( Norm1(Nsiz) )
    Norm1 = 0.0d0 ! to start with
   
!    if (.not.allocated(Fnnxx)) allocate(Fnnxx(Nsiz,Nsiz))
!    if (.not.allocated(Fnnxy)) allocate(Fnnxy(Nsiz,Nsiz))
!    if (.not.allocated(Fnnxz)) allocate(Fnnxz(Nsiz,Nsiz))
!    if (.not.allocated(Fnnyx)) allocate(Fnnyx(Nsiz,Nsiz))
!    if (.not.allocated(Fnnyy)) allocate(Fnnyy(Nsiz,Nsiz))
!    if (.not.allocated(Fnnyz)) allocate(Fnnyz(Nsiz,Nsiz))
!    if (.not.allocated(Fnnzx)) allocate(Fnnzx(Nsiz,Nsiz))
!    if (.not.allocated(Fnnzy)) allocate(Fnnzy(Nsiz,Nsiz))
!    if (.not.allocated(Fnnzz)) allocate(Fnnzz(Nsiz,Nsiz))
   allocate(Fnn_temp_x(Nsiz,Nsiz))
   allocate(Fnn_temp_y(Nsiz,Nsiz))
   allocate(Fnn_temp_z(Nsiz,Nsiz))
   allocate(PBx(Nsiz,Nsiz))
   allocate(PBy(Nsiz,Nsiz))
   allocate(PBz(Nsiz,Nsiz))

   ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
      select type(ARRAY)
      type is (TB_H_Pettifor)	! orthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz) ! module "TB"
      type is (TB_H_Molteni)	! orthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz)  ! module "TB"
      type is (TB_H_Fu)		! orthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz) ! module "TB"
      type is (TB_H_NRL)	! nonorthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz, Scell(NSC)%Sij, CSij) ! module "TB"
      type is (TB_H_DFTB) 	! nonorthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz, Scell(NSC)%Sij, CSij) ! module "TB"
      type is (TB_H_3TB) 	! nonorthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz, Scell(NSC)%Sij, CSij) ! module "TB"
      type is (TB_H_xTB) 	! nonorthogonal
         call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, Err, cPRRx, cPRRy, cPRRz, Scell(NSC)%Sij, CSij) ! module "TB"
      end select
   END ASSOCIATE
   
   ! test:
!    print*, ' kx, ky, kz ',  kx, ky, kz    
!     do j = 250,size(CHij,1)
!        do i = 250,size(CHij,2)
!          if ( ABS(dble(CHij(j,i)))  > 1.0d-14 ) print*, j, i, CHij(j,i) !, CHij_test(j,i)
!        enddo
!     enddo
!     PAUSE ' CHij - 2'

   do j = 1,size(CHij,1)
      do i = 1,size(CHij,2)
         if (isnan(REAL(CHij(j,i))) .or. isnan(AIMAG(CHij(j,i)))) then
             Error_descript = 'Module Optical_parameters: subroutine get_Fnn_complex got NaNs'
             call Save_error_details(Err, 8, Error_descript)
             print*, trim(adjustl(Error_descript)), j, i, CHij(j,i)
         endif
      enddo
   enddo


#ifdef MPI_USED
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = Nsiz
   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, Nsiz ! ensure WF normalization to 1
      Norm1(i) = SQRT( SUM( conjg(CHij(:,i)) * CHij(:,i) ) )
   enddo
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_Fnn_complex:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Norm1', Norm1) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, Nsiz
      do nn = 1, Nsiz
         PBx(i,nn) = SUM(cPRRx(i,:)*CHij(:,nn)) / Norm1(nn)
         PBy(i,nn) = SUM(cPRRy(i,:)*CHij(:,nn)) / Norm1(nn)
         PBz(i,nn) = SUM(cPRRz(i,:)*CHij(:,nn)) / Norm1(nn)
      enddo ! j
   enddo ! i
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'PBx', PBx) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'PBy', PBy) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'PBz', PBz) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

   do n = Nstart, Nend, N_incr  ! each process does its own part
   !do n = 1, Nsiz
      do nn = 1, Nsiz
         Fnn_temp_x(n,nn) = SUM(conjg(CHij(:,n))*PBx(:,nn)) / Norm1(n)
         Fnn_temp_y(n,nn) = SUM(conjg(CHij(:,n))*PBy(:,nn)) / Norm1(n)
         Fnn_temp_z(n,nn) = SUM(conjg(CHij(:,n))*PBz(:,nn)) / Norm1(n)
      enddo ! nn
   enddo ! n
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnn_temp_x', Fnn_temp_x) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnn_temp_y', Fnn_temp_y) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnn_temp_z', Fnn_temp_z) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, Nsiz
      do j = 1, Nsiz
!          Fnnx(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_x(i,j)
!          Fnny(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_y(i,j)
!          Fnnz(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_z(i,j)

!        ! Must convert it to double complex, otherwise it gives wrong results!
         Fnnx(i,j) = DCMPLX(conjg(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j))
         Fnny(i,j) = DCMPLX(conjg(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         Fnnz(i,j) = DCMPLX(conjg(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_z(i,j))
         if (present(Fnnxy)) then ! off-diagonal
            Fnnxy(i,j) = DCMPLX(conjg(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_y(i,j))
            Fnnxz(i,j) = DCMPLX(conjg(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_z(i,j))
            Fnnyx(i,j) = DCMPLX(conjg(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_x(i,j))
            Fnnyz(i,j) = DCMPLX(conjg(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_z(i,j))
            Fnnzx(i,j) = DCMPLX(conjg(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_x(i,j))
            Fnnzy(i,j) = DCMPLX(conjg(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         endif
      enddo
   enddo
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnx', Fnnx) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnny', Fnny) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnz', Fnnz) ! module "MPI_subroutines"
   if (present(Fnnxy)) then ! off-diagonal
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnxy', Fnnxy) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnxz', Fnnxz) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnyx', Fnnyx) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnyz', Fnnyz) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnzx', Fnnzx) ! module "MPI_subroutines"
      call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//'Fnnzy', Fnnzy) ! module "MPI_subroutines"
   endif
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp PARALLEL private(i,nn,n,j)
   !$omp do
   do i = 1, Nsiz ! ensure WF normalization to 1
      Norm1(i) = SQRT( SUM( conjg(CHij(:,i)) * CHij(:,i) ) )
   enddo
   !$omp end do
   !$OMP BARRIER
   !$omp do
   do i = 1, Nsiz
      do nn = 1, Nsiz
         PBx(i,nn) = SUM(cPRRx(i,:)*CHij(:,nn)) / Norm1(nn)
         PBy(i,nn) = SUM(cPRRy(i,:)*CHij(:,nn)) / Norm1(nn)
         PBz(i,nn) = SUM(cPRRz(i,:)*CHij(:,nn)) / Norm1(nn)
      enddo ! j
   enddo ! i
   !$omp end do
   !$OMP BARRIER
   !$omp do
   do n = 1, Nsiz
      do nn = 1, Nsiz
         Fnn_temp_x(n,nn) = SUM(conjg(CHij(:,n))*PBx(:,nn)) / Norm1(n)
         Fnn_temp_y(n,nn) = SUM(conjg(CHij(:,n))*PBy(:,nn)) / Norm1(n)
         Fnn_temp_z(n,nn) = SUM(conjg(CHij(:,n))*PBz(:,nn)) / Norm1(n)
      enddo ! nn
   enddo ! n
   !$omp end do
   !$OMP BARRIER
   !$omp do
   do i = 1, Nsiz
      do j = 1, Nsiz
!          Fnnx(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_x(i,j)
!          Fnny(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_y(i,j)
!          Fnnz(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_z(i,j)

!        ! Must convert it to double complex, otherwise it gives wrong results!
         Fnnx(i,j) = DCMPLX(conjg(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j))
         Fnny(i,j) = DCMPLX(conjg(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         Fnnz(i,j) = DCMPLX(conjg(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_z(i,j))
         if (present(Fnnxy)) then ! off-diagonal
            Fnnxy(i,j) = DCMPLX(conjg(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_y(i,j))
            Fnnxz(i,j) = DCMPLX(conjg(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_z(i,j))
            Fnnyx(i,j) = DCMPLX(conjg(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_x(i,j))
            Fnnyz(i,j) = DCMPLX(conjg(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_z(i,j))
            Fnnzx(i,j) = DCMPLX(conjg(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_x(i,j))
            Fnnzy(i,j) = DCMPLX(conjg(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         endif
      enddo
   enddo
   !$omp end do
   !$omp end parallel
#endif

!   open(NEWUNIT=FN, FILE = 'OUTPUT_PBx.dat')
!      do i = 1, size(PBx,1)
!         do j = 1, size(PBx,2)
!            write(FN, '(i,i,e)') i, j, PBx(i,j)
!         enddo
!      enddo
!   close(FN)
!   pause 'get_Fnn'

   ! Cleaning up:
   call deallocate_array(Ei_r)     ! module "Little_subroutines"
   call deallocate_array(Norm1) ! module "Little_subroutines"
   if (allocated(CHij)) deallocate(CHij)
   if (allocated(Fnn_temp_x)) deallocate(Fnn_temp_x)
   if (allocated(Fnn_temp_y)) deallocate(Fnn_temp_y)
   if (allocated(Fnn_temp_z)) deallocate(Fnn_temp_z)
   if (allocated(PBx)) deallocate(PBx)
   if (allocated(PBy)) deallocate(PBy)
   if (allocated(PBz)) deallocate(PBz)
   if (allocated(cPRRx)) deallocate(cPRRx)
   if (allocated(cPRRy)) deallocate(cPRRy)
   if (allocated(cPRRz)) deallocate(cPRRz)
!    if (allocated(Fnnxy)) deallocate(Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy)
end subroutine get_Fnn_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real-value Trani's subroutines (Gamma-point):
subroutine get_trani_all(numpar, Scell, NSC, Ei, Ha, fe, all_w) ! Ref. [2]
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: Ei	! energy levels [eV]
   real(8), dimension(:,:), intent(in) :: Ha	! eigenvectors of the hamiltonian
   real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   real(8), dimension(size(Ei),size(Ei)) :: Fnnx, Fnny, Fnnz
!    real(8), allocatable, dimension(:,:) :: Fnnx, Fnny, Fnnz
   real(8), dimension(size(Ei),size(Ei)) :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), dimension(:,:), allocatable :: Eps
   real(8), dimension(:), allocatable :: Re_CDF
   real(8) Re_eps, Im_eps, dc_cond, loss_f
   real(8) w
   integer i, N, FN, SEi
   integer*4 today(3), now(3)

   SEi = size(Ei)
!    allocate(Fnnx(SEi,SEi))
!    allocate(Fnny(SEi,SEi))
!    allocate(Fnnz(SEi,SEi))

!    call print_time_step('Start opt time:', g_time, msec=.true.)   ! module "Little_subroutines"call print_
   
   call get_Fnn(numpar, Scell, NSC, Ha, Ei, Fnnx, Fnny, Fnnz) ! see below
    ! With off-diagonal elements:
!     call get_Fnn(Scell, NSC, Ha, Ei, Fnnxx, Fnnyy, Fnnzz, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy)

!    call print_time_step('get_Fnn:', g_time, msec=.true.)   ! module "Little_subroutines"call print_

   if (all_w) then ! full spectrum:
!       wmax = 50.0d0*g_e/g_h	! [1/s] maximum frequency to trace
!       N = 1000
!       allocate(Eps(2,N))
!       if (.not.allocated(Scell(NSC)%eps%Eps_hw)) then
!          allocate(Scell(NSC)%eps%Eps_hw(9,N)) ! all are there
!       endif
      
      ! Allocate the array of optical coefficients spectrum (if not allocated before):
      call allocate_Eps_hw(Scell(NSC)%eps%E_min, Scell(NSC)%eps%E_max, Scell(NSC)%eps%dE, Scell(NSC)%eps%Eps_hw) ! see below
      N = size(Scell(NSC)%eps%Eps_hw,2) ! use this size to define temporary arrays within this subroutine
      allocate(Eps(2,N))

      !w = 0.0d0
      w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] maximum frequency to trace, starting point for the optical spectrum [eV]
      do i = 1, N
!          w = real(i)/real(N)*50.0d0*g_e/g_h
         call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, dc_cond=dc_cond)
         ! With off-diagonal elements:
         !call get_trani(numpar, Scell, NSC, Fnnxx, Fnnyy, Fnnzz, Ei, w, Re_eps, Im_eps, dc_cond=dc_cond, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy)
!          call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps)
         Eps(1,i) = Re_eps
         Eps(2,i) = Im_eps
         !write(*,'(f,es,es,es)') w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps)
         !write(FN,'(f,es,es,es,es,es,es,es,es)') w/g_e*g_h, Eps(1,i), Eps(2,i), Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), Scell(NSC)%eps%R, Scell(NSC)%eps%T, Scell(NSC)%eps%A, Scell(NSC)%eps%n, Scell(NSC)%eps%k
         loss_f = Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps)
         call save_Eps_hw(Scell(NSC)%eps%Eps_hw, i, w/g_e*g_h, Eps(1,i), Eps(2,i), loss_f, Scell(NSC)%eps%R, Scell(NSC)%eps%T, Scell(NSC)%eps%A, Scell(NSC)%eps%n, Scell(NSC)%eps%k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond)
         w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] next frequency grid point [eV]
      enddo
      !close(FN)
 
      deallocate(Eps)
      !pause 'All spectrum is done'
   else ! only for given probe:
      w = Scell(NSC)%eps%w
      ! Diagonal elements only:
      call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, dc_cond=dc_cond)
      ! With off-diagonal elements:
!       call get_trani(numpar, Scell, NSC, Fnnxx, Fnnyy, Fnnzz, Ei, w, Re_eps, Im_eps, dc_cond=dc_cond, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy)
   endif
   
!    deallocate(Fnnx, Fnny, Fnnz)
!    deallocate(Fnnxx, Fnnyy, Fnnzz, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy)
   
!    call print_time_step('get_trani_all:', g_time, msec=.true.)   ! module "Little_subroutines"call print_
end subroutine get_trani_all



subroutine get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, n, k, dc_cond, &
                     Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy) ! Ref. [1]
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: Ei	! energy levels [eV]
!    real(8), dimension(:), intent(in), target :: fe	! electron distribution function
   real(8), intent(in) :: w	! [1/s] frequency of the probe pulse
   real(8), intent(out) :: Re_eps, Im_eps ! real and imaginary parts of the dielectric function
   real(8), dimension(:,:), intent(in) :: Fnnx, Fnny, Fnnz
   real(8), dimension(:,:), intent(in), optional :: Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy
!    real(8), dimension(:,:), intent(in) :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), intent(out), optional :: R, T, A, n, k, dc_cond
   type(Error_handling) :: Err	! errors
   !------------------------------------------------------
   real(8)  A_const, V, fij, Eij, Eijw, whe, e2, he, g, hg, hg2, temp, hw, temp2, temp3
   real(8), pointer :: fi, fj
   real(8), dimension(size(Ei)), target :: fe_temp	! electron distribution function
   real(8) Re(3), Im(3) 	! xx, yy, zz
   real(8) Re_off(6), Im_off(6) 	! xy, xz, yx, yz, zx, zy
   real(8) Re_mid(3), Re_small(3), Re_large(3) 	! x, y, z
   real(8) Im_mid(3), Im_small(3), Im_large(3)	! x, y, z
   real(8) Re_mid_off(6), Re_small_off(6), Re_large_off(6)	! xy, xz, yx, yz, zx, zy
   real(8) Im_mid_off(6), Im_small_off(6), Im_large_off(6)	    ! xy, xz, yx, yz, zx, zy
!    real(8), dimension(:,:), allocatable :: Re_array, Im_array
   real(8) dc(3) ! x, y, z
!    real(8) Ret(3,3), Imt(3,3) ! x, y, z
!    real(8) Re_eps_t(3,3), Im_eps_t(3,3)
   integer i, j, NEL, FN
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part


   NEL = size(Ei)
   
!    allocate(Re_array(NEL,3))
!    allocate(Im_array(NEL,3))
   
!    g = 1.0d13	! [1/s] gamma parameter
!    g = 1.5d14	! [1/s] gamma parameter
   g = m_gamm  ! [1/s] gamma parameter

!    V = Scell(NSC)%V*1d-30
   V = Scell(NSC)%V*1d-10	! power of 20 cancels with the same power in Fnn'
   
!    he = g_h/g_e
   e2 = g_e*g_e
   hw = w*g_h
!    A_const = e2*g_h*g_h/(g_me*g_me*V*g_e0)
! mass and Plank constant cancel out in the final expression:
   A_const = e2/(V*g_e0)

   hg = g_h*g
   hg2 = hg*hg
   Re = 0.0d0
   Im = 0.0d0
   dc = 0.0d0
   
!    Re_array = 0.0d0
!    Im_array = 0.0d0
   Re_mid = 0.0d0
   Re_small = 0.0d0
   Re_large = 0.0d0
   Im_mid = 0.0d0
   Im_small = 0.0d0
   Im_large = 0.0d0
!    Ret = 0.0d0
!    Imt = 0.0d0

   call set_Fermi(Ei,Scell(NSC)%TeeV,Scell(NSC)%mu,fe_temp)	! module "Electron_tools"
   !-----------------------------------------------
   ! Project-specific routine:
   ! Testing Erf-function distribution:
!    call set_Erf_distribution(Ei,Scell(NSC)%TeeV,Scell(NSC)%mu,fe_temp)	! module "Electron_tools"
   !-----------------------------------------------

#ifdef MPI_USED   ! use the MPI version [tested]
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = NEL
   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1,NEL
      fi => fe_temp(i)
      do j = 1,NEL
         fj => fe_temp(j)
         fij = (fi - fj)
         Eij = (Ei(i) - Ei(j))
         if ((abs(Eij) > 1d-13) .and. (abs(fij) > 1.0d-14)) then
            Eijw = hw - Eij*g_e	! [J]
            temp = fij/(Eij*Eij)/(Eijw*Eijw + hg2)
            temp2 = temp*Eijw

            ! To avoid the problem of large and small numbers, save them independently:
            temp3 = temp2*Fnnx(i,j) ! xx
            call split_into_orders(temp3, Re_small(1), Re_mid(1), Re_large(1))   ! below
            temp3 = temp*Fnnx(i,j)  ! xx
            call split_into_orders(temp3, Im_small(1), Im_mid(1), Im_large(1))   ! below

            temp3 = temp2*Fnny(i,j) ! yy
            call split_into_orders(temp3, Re_small(2), Re_mid(2), Re_large(2))   ! below
            temp3 = temp*Fnny(i,j)  ! yy
            call split_into_orders(temp3, Im_small(2), Im_mid(2), Im_large(2))   ! below

            temp3 = temp2*Fnnz(i,j) ! zz
            call split_into_orders(temp3, Re_small(3), Re_mid(3), Re_large(3))   ! below
            temp3 = temp*Fnnz(i,j)   ! zz
            call split_into_orders(temp3, Im_small(3), Im_mid(3), Im_large(3))   ! below

            ! Off-diagonal elements:
            if (present(Fnnxy)) then
               temp3 = temp2*Fnnxy(i,j) ! xy
               call split_into_orders(temp3, Re_small_off(1), Re_mid_off(1), Re_large_off(1))   ! below
               temp3 = temp*Fnnxy(i,j)  ! xy
               call split_into_orders(temp3, Im_small_off(1), Im_mid_off(1), Im_large_off(1))   ! below

               temp3 = temp2*Fnnxz(i,j) ! xz
               call split_into_orders(temp3, Re_small_off(2), Re_mid_off(2), Re_large_off(2))   ! below
               temp3 = temp*Fnnxz(i,j)  ! xz
               call split_into_orders(temp3, Im_small_off(2), Im_mid_off(2), Im_large_off(2))   ! below

               Re_small_off(3) = Re_small_off(1)    ! yx
               Re_large_off(3) = Re_large_off(1)
               Re_mid_off(3) = Re_mid_off(1)
               Im_small_off(3) = Im_small_off(1)    ! yx
               Im_large_off(3) = Im_large_off(1)
               Im_mid_off(3) = Im_mid_off(1)

               temp3 = temp2*Fnnyz(i,j) ! yz
               call split_into_orders(temp3, Re_small_off(4), Re_mid_off(4), Re_large_off(4))   ! below
               temp3 = temp*Fnnyz(i,j)  ! yz
               call split_into_orders(temp3, Im_small_off(4), Im_mid_off(4), Im_large_off(4))   ! below

               Re_small_off(5) = Re_small_off(2)    ! zx
               Re_large_off(5) = Re_large_off(2)
               Re_mid_off(5) = Re_mid_off(2)
               Im_small_off(5) = Im_small_off(2)    ! zx
               Im_large_off(5) = Im_large_off(2)
               Im_mid_off(5) = Im_mid_off(2)

               Re_small_off(6) = Re_small_off(4)    ! zy
               Re_large_off(6) = Re_large_off(4)
               Re_mid_off(6) = Re_mid_off(4)
               Im_small_off(6) = Im_small_off(4)    ! zy
               Im_large_off(6) = Im_large_off(4)
               Im_mid_off(6) = Im_mid_off(4)
            endif
         endif
      enddo
   enddo
   error_part = 'Error in get_trani:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_mid}', Re_mid) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_small}', Re_small) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_large}', Re_large) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_mid}', Im_mid) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_small}', Im_small) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_large}', Im_large) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_small_off}', Re_small_off) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_mid_off}', Re_mid_off) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Re_large_off}', Re_large_off) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_small_off}', Im_small_off) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_mid_off}', Im_mid_off) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Im_large_off}', Im_large_off) ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp PARALLEL private(i,j,fi,fj,fij,Eij,Eijw,temp,temp2, temp3)
!    !$omp do reduction( + : Re,Im)
   !$omp do reduction( + : Re_mid, Re_small, Re_large, Im_mid, Im_small, Im_large, Re_small_off, Re_mid_off, Re_large_off, Im_small_off, Im_mid_off, Im_large_off)
   do i = 1,NEL
      !fi = fe(i)
      fi => fe_temp(i)
      do j = 1,NEL
         fj => fe_temp(j)
         fij = (fi - fj)
         Eij = (Ei(i) - Ei(j))
         !if ((i .NE. j) .and. (abs(Eij) > 1d-14) .and. (abs(fij) > 1.0d-14)) then
         if ((abs(Eij) > 1d-13) .and. (abs(fij) > 1.0d-14)) then
            !fj = fe(j)
!             Eij = Eij*g_e
            Eijw = hw - Eij*g_e	! [J]
!             print*, 'TEMP:', fij, Ei(i), Ei(j)
!             print*, Eij
!             print*, Eijw
!             print*, (Eijw*Eijw + hg2)
!             print*, 'test 0'
            temp = fij/(Eij*Eij)/(Eijw*Eijw + hg2)
            temp2 = temp*Eijw
            
            ! To avoid the problem of large and small numbers, save them independently:
            temp3 = temp2*Fnnx(i,j) ! xx
            call split_into_orders(temp3, Re_small(1), Re_mid(1), Re_large(1))   ! below 
            temp3 = temp*Fnnx(i,j)  ! xx
            call split_into_orders(temp3, Im_small(1), Im_mid(1), Im_large(1))   ! below 
            
            temp3 = temp2*Fnny(i,j) ! yy
            call split_into_orders(temp3, Re_small(2), Re_mid(2), Re_large(2))   ! below 
            temp3 = temp*Fnny(i,j)  ! yy
            call split_into_orders(temp3, Im_small(2), Im_mid(2), Im_large(2))   ! below 
            
            temp3 = temp2*Fnnz(i,j) ! zz
            call split_into_orders(temp3, Re_small(3), Re_mid(3), Re_large(3))   ! below 
            temp3 = temp*Fnnz(i,j)   ! zz
            call split_into_orders(temp3, Im_small(3), Im_mid(3), Im_large(3))   ! below 
            
            ! Off-diagonal elements:
            if (present(Fnnxy)) then
               temp3 = temp2*Fnnxy(i,j) ! xy
               call split_into_orders(temp3, Re_small_off(1), Re_mid_off(1), Re_large_off(1))   ! below
               temp3 = temp*Fnnxy(i,j)  ! xy
               call split_into_orders(temp3, Im_small_off(1), Im_mid_off(1), Im_large_off(1))   ! below
               
               temp3 = temp2*Fnnxz(i,j) ! xz
               call split_into_orders(temp3, Re_small_off(2), Re_mid_off(2), Re_large_off(2))   ! below 
               temp3 = temp*Fnnxz(i,j)  ! xz
               call split_into_orders(temp3, Im_small_off(2), Im_mid_off(2), Im_large_off(2))   ! below 
            
               Re_small_off(3) = Re_small_off(1)    ! yx
               Re_large_off(3) = Re_large_off(1)
               Re_mid_off(3) = Re_mid_off(1)
               Im_small_off(3) = Im_small_off(1)    ! yx
               Im_large_off(3) = Im_large_off(1)
               Im_mid_off(3) = Im_mid_off(1)
               
               temp3 = temp2*Fnnyz(i,j) ! yz
               call split_into_orders(temp3, Re_small_off(4), Re_mid_off(4), Re_large_off(4))   ! below
               temp3 = temp*Fnnyz(i,j)  ! yz
               call split_into_orders(temp3, Im_small_off(4), Im_mid_off(4), Im_large_off(4))   ! below
               
               Re_small_off(5) = Re_small_off(2)    ! zx
               Re_large_off(5) = Re_large_off(2)
               Re_mid_off(5) = Re_mid_off(2)
               Im_small_off(5) = Im_small_off(2)    ! zx
               Im_large_off(5) = Im_large_off(2)
               Im_mid_off(5) = Im_mid_off(2)
               
               Re_small_off(6) = Re_small_off(4)    ! zy
               Re_large_off(6) = Re_large_off(4)
               Re_mid_off(6) = Re_mid_off(4)
               Im_small_off(6) = Im_small_off(4)    ! zy
               Im_large_off(6) = Im_large_off(4)
               Im_mid_off(6) = Im_mid_off(4)
            endif
            
            ! dc-conductivity
!             dc(1) = dc(1) + fij*Fnnx(i,j)
!             dc(2) = dc(2) + fij*Fnny(i,j)
!             dc(3) = dc(3) + fij*Fnnz(i,j)
!             Ret(1,1) = Ret(1,1) + temp2*Fnnxx(i,j) ! xx
!             Imt(1,1) = Imt(1,1) + temp*Fnnxx(i,j)
!             Ret(1,2) = Ret(1,2) + temp2*Fnnxy(i,j) ! xy
!             Imt(1,2) = Imt(1,2) + temp*Fnnxy(i,j)
!             Ret(1,3) = Ret(1,3) + temp2*Fnnxz(i,j) ! xz
!             Imt(1,3) = Imt(1,3) + temp*Fnnxz(i,j)
!             Ret(2,1) = Ret(2,1) + temp2*Fnnyx(i,j) ! yx
!             Imt(2,1) = Imt(2,1) + temp*Fnnyx(i,j)
!             Ret(2,2) = Ret(2,2) + temp2*Fnnyy(i,j) ! yy
!             Imt(2,2) = Imt(2,2) + temp*Fnnyy(i,j)
!             Ret(2,3) = Ret(2,3) + temp2*Fnnyz(i,j) ! yz
!             Imt(2,3) = Imt(2,3) + temp*Fnnyz(i,j)
!             Ret(3,1) = Ret(3,1) + temp2*Fnnzx(i,j) ! zx
!             Imt(3,1) = Imt(3,1) + temp*Fnnzx(i,j)
!             Ret(3,2) = Ret(3,2) + temp2*Fnnzy(i,j) ! zy
!             Imt(3,2) = Imt(3,2) + temp*Fnnzy(i,j)
!             Ret(3,3) = Ret(3,3) + temp2*Fnnzz(i,j) ! zz
!             Imt(3,3) = Imt(3,3) + temp*Fnnzz(i,j)
         endif
      enddo
   enddo
   !$omp end do
   !$omp end parallel
#endif

   ! test:
!    write(*,'(a,es,es,es,es,es,es)') 'Re(1)', Re(1), Re_mid(1), Re_small(1), Im(1), Im_mid(1), Im_small(1)
!    write(*,'(a,es,es,es,es,es,es)') 'Re(2)', Re(2), Re_mid(2), Re_small(2), Im(2), Im_mid(2), Im_small(2)
!    write(*,'(a,es,es,es,es,es,es)') 'Re(3)', Re(3), Re_mid(3), Re_small(3), Im(3), Im_mid(3), Im_small(3)

   ! Sum up large and small numbers contributions with the normal ones:
   Re(:) = Re_mid(:) + Re_small(:) + Re_large(:)
   Im(:) = Im_mid(:) + Im_small(:) + Im_large(:)
   
   ! Off-diagonals:
   if (present(Fnnxy)) then
      Re_off(:) = Re_mid_off(:) + Re_small_off(:) + Re_large_off(:)
      Im_off(:) = Im_mid_off(:) + Im_small_off(:) + Im_large_off(:)
   endif
   
!    print*, 'S', Re_small(:) 
!    print*, 'L', Re_large(:)
!    print*, 'M', Re_mid(:)
!    print*, 'IS', Im_small(:) 
!    print*, 'IL', Im_large(:)
!    print*, 'IM', Im_mid(:)

   ! Aferage over all directions:
   !Re_eps = 1.0d0 + A_const*SUM(Re)/3.0d0 ! average over x, y, and z
   !Im_eps = -hg*A_const*SUM(Im)/3.0d0     ! average over x, y, and z
   
   ! With off-diagonals:
   Re_eps = 1.0d0 + A_const*SUM(Re)/3.0d0  ! average over x, y, and z
   Im_eps = -hg*A_const*SUM(Im)/3.0d0         ! average over x, y, and z
   ! Off-diagonal elements:
   if (present(Fnnxy)) then
      Re_eps = Re_eps + A_const*SUM(Re_off)/6.0d0
      Im_eps = Im_eps - hg*A_const*SUM(Im_off)/6.0d0
   endif
   
   ! Save components of the calculated CDF:
   Scell(NSC)%eps%Eps_xx = dcmplx(1.0d0 + A_const*Re(1), -hg*A_const*Im(1))
   Scell(NSC)%eps%Eps_yy = dcmplx(1.0d0 + A_const*Re(2), -hg*A_const*Im(2))
   Scell(NSC)%eps%Eps_zz = dcmplx(1.0d0 + A_const*Re(3), -hg*A_const*Im(3))
   if (present(Fnnxy)) then ! off-diagonal elements:
      Scell(NSC)%eps%Eps_xy = dcmplx(A_const*Re_off(1), -hg*A_const*Im_off(1))
      Scell(NSC)%eps%Eps_xz = dcmplx(A_const*Re_off(2), -hg*A_const*Im_off(2))
      Scell(NSC)%eps%Eps_yx = dcmplx(A_const*Re_off(3), -hg*A_const*Im_off(3))
      Scell(NSC)%eps%Eps_yz = dcmplx(A_const*Re_off(4), -hg*A_const*Im_off(4))
      Scell(NSC)%eps%Eps_zx = dcmplx(A_const*Re_off(5), -hg*A_const*Im_off(5))
      Scell(NSC)%eps%Eps_zy = dcmplx(A_const*Re_off(6), -hg*A_const*Im_off(6))
   endif

!    print*, 'w', w, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, Scell(NSC)%eps%Eps_xy, Scell(NSC)%eps%Eps_xz, Scell(NSC)%eps%Eps_yx, Scell(NSC)%eps%Eps_yz, Scell(NSC)%eps%Eps_zx, Scell(NSC)%eps%Eps_zy

   if (present(dc_cond)) then	! save the dc-conductivity                                                            
      dc_cond = Im_eps*w*g_e0     ! averaged over x, y, and z
      Scell(NSC)%eps%dc_cond = dc_cond
!       print*, 'dc_cond', dc_cond
!       print*, 'cd_con2', Im_eps*w*g_e0
   endif

!    ! Average a tensor:
!    Re_eps = 0.0d0
!    Im_eps = 0.0d0
!    Re_eps_t = 0.0d0
!    Im_eps_t = 0.0d0
!    do i = 1,3
!       do j = 1,3
!          if (i .EQ. j) then 
!             Re_eps_t(i,j) = 1.0d0 + A_const*Ret(i,j) ! average over x, y, and z
!          else
!             Re_eps_t(i,j) = A_const*Ret(i,j) ! average over x, y, and z
!          endif
!          Im_eps_t(i,j) = -hg*A_const*Imt(i,j)     ! average over x, y, and z
!       enddo
!    enddo
!    call sym_diagonalize(Re_eps_t, Re, Err) ! "TB_Hamiltonian" modeule
!    call sym_diagonalize(Im_eps_t, Im, Err) ! "TB_Hamiltonian" modeule
!    Re_eps = SUM(Re)/3.0d0
!    Im_eps = SUM(Im)/3.0d0

   ! Don't average, save only one:
!     Re_eps = 1.0d0 + A_const*Re(3) ! x, y, or z
!     Im_eps = -hg*A_const*Im(3)     ! x, y, or z
   ! Average for graphite:
!     Re_eps = 1.0d0 + A_const*(Re(1)+Re(2))/2.0d0   ! average over x and y (graphite planes)
!     Im_eps = -hg*A_const*(Im(1)+Im(2))/2.0d0     ! average over x and y (graphite planes)

   Scell(NSC)%eps%ReEps = Re_eps
   Scell(NSC)%eps%ImEps = Im_eps

   !call get_RTA_from_CDF(Scell(NSC)%eps%ReEps, Scell(NSC)%eps%ImEps, Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, Scell(NSC)%eps%R, Scell(NSC)%eps%T, Scell(NSC)%eps%A)
   !call get_n_k(Scell(NSC)%eps%ReEps, Scell(NSC)%eps%ImEps, Scell(NSC)%eps%n, Scell(NSC)%eps%k)	! transfer Re(e) and Im(e) into n and k
   if (present(R)) then
      call get_RTA_from_CDF(Re_eps, Im_eps, Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, R, T, A)
      Scell(NSC)%eps%R = R
      Scell(NSC)%eps%T = T
      Scell(NSC)%eps%A = A
      call get_n_k(Re_eps, Im_eps, n, k)	! transfer Re(e) and Im(e) into n and k
      Scell(NSC)%eps%n = n
      Scell(NSC)%eps%k = k
   else
      call get_RTA_from_CDF(Re_eps, Im_eps, Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, Scell(NSC)%eps%R, Scell(NSC)%eps%T, Scell(NSC)%eps%A)
      call get_n_k(Re_eps, Im_eps, Scell(NSC)%eps%n, Scell(NSC)%eps%k)	! transfer Re(e) and Im(e) into n and k
   endif

! print*, 'n, K:', Scell(NSC)%eps%n, Scell(NSC)%eps%k
! print*, 'Full:', Re_eps, Im_eps
! print*, 'Re', 1.0d0 + A_const*(Re)
! print*, 'Im', -hg*A_const*(Im)
   nullify (fi, fj)
end subroutine get_trani


subroutine split_into_orders(temp, V_small, V_mid, V_large)
   real(8), intent(in) :: temp  ! value to be sorted
   real(8), intent(inout) :: V_small, V_mid, V_large    ! values to be sorted into
   if (ABS(V_mid) > 1.0d-10) then
      if (ABS(temp) < 1d-6*ABS(V_mid)) then	! too small values
         V_small = V_small + temp
      else if (ABS(temp) > 1d6*ABS(V_mid)) then	! too large values
         V_large = V_large + temp
      else ! normal values
         V_mid = V_mid + temp
      endif
   else
      V_mid = V_mid + temp
   endif
end subroutine split_into_orders


subroutine get_Fnn(numpar, Scell, NSC, Ha, Ei, Fnnx, Fnny, Fnnz, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy) ! Ref. [2]
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), intent(in) :: Ei	! energy levels [eV]
   real(8), dimension(:,:), intent(in) :: Ha	! eigenvectors of the hamiltonian
   real(8), dimension(:,:), intent(out) :: Fnnx, Fnny, Fnnz ! diagonal elements
   real(8), dimension(:,:), intent(out), optional :: Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy   ! off-diagonal elements
   real(8), allocatable, dimension(:,:) :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   real(8), allocatable, dimension(:,:) :: PBx, PBy, PBz
   real(8), allocatable, dimension(:) :: Norm1
!    real(8), dimension(size(Fnnxx,1),size(Fnnxx,2)) :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
!    real(8), dimension(size(Fnnxx,1),size(Fnnxx,2)) :: PBx, PBy, PBz
   integer i, j, m, n, nn, FN
   integer :: N_incr, Nstart, Nend
   character(100) :: error_part

   m = size(Ha,1)
   allocate(Fnn_temp_x(size(Fnnx,1),size(Fnnx,2)), source = 0.0d0)
   allocate(Fnn_temp_y(size(Fnnx,1),size(Fnnx,2)), source = 0.0d0)
   allocate(Fnn_temp_z(size(Fnnx,1),size(Fnnx,2)), source = 0.0d0)
   allocate(PBx(size(Fnnx,1),size(Fnnx,2)), source = 0.0d0)
   allocate(PBy(size(Fnnx,1),size(Fnnx,2)), source = 0.0d0)
   allocate(PBz(size(Fnnx,1),size(Fnnx,2)), source = 0.0d0)
   allocate(Norm1(size(Fnnx,1)), source = 0.0d0)

   ! Do it with mkl:
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Scell(NSC)%PRRx, m, Ha, m, 0.0d0, PBx, m) ! Original
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Scell(NSC)%PRRy, m, Ha, m, 0.0d0, PBy, m)
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Scell(NSC)%PRRz, m, Ha, m, 0.0d0, PBz, m)
   
   CALL dgemm ('N','N', m, m, m, 1.0d0, Scell(NSC)%PRRx, m, Ha, m, 0.0d0, PBx, m)
   CALL dgemm ('N','N', m, m, m, 1.0d0, Scell(NSC)%PRRy, m, Ha, m, 0.0d0, PBy, m)
   CALL dgemm ('N','N', m, m, m, 1.0d0, Scell(NSC)%PRRz, m, Ha, m, 0.0d0, PBz, m)
   
#ifdef MPI_USED   ! use the MPI version
   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank   ! starting point for each process
   Nend = m
   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, m ! ensure WF normalization to 1
      Norm1(i) = DSQRT(SUM( Ha(i,:) * Ha(i,:) ))
   enddo
   ! Collect information from all processes into the master process, and distribute the final arrays to all processes:
   error_part = 'Error in get_Fnn:'
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Norm1}', Norm1) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, m
      PBx(:,i) = PBx(:,i) / Norm1(i)
      PBy(:,i) = PBy(:,i) / Norm1(i)
      PBz(:,i) = PBz(:,i) / Norm1(i)
   enddo ! i
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {PBx}', PBx) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {PBy}', PBy) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {PBz}', PBz) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

   ! Do it with mkl:
   CALL dgemm ('T','N', m, m, m, 1.0d0, PBx, m, Ha, m, 0.0d0, Fnn_temp_x, m)
   CALL dgemm ('T','N', m, m, m, 1.0d0, PBy, m, Ha, m, 0.0d0, Fnn_temp_y, m)
   CALL dgemm ('T','N', m, m, m, 1.0d0, PBz, m, Ha, m, 0.0d0, Fnn_temp_z, m)

   ! Do the cycle (parallel) calculations:
   do i = Nstart, Nend, N_incr  ! each process does its own part
   !do i = 1, m
      Fnn_temp_x(:,i) = Fnn_temp_x(:,i) / Norm1(i)
      Fnn_temp_y(:,i) = Fnn_temp_y(:,i) / Norm1(i)
      Fnn_temp_z(:,i) = Fnn_temp_z(:,i) / Norm1(i)
   enddo ! i
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Fnn_temp_x}', Fnn_temp_x) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Fnn_temp_y}', Fnn_temp_y) ! module "MPI_subroutines"
   call do_MPI_Allreduce(numpar%MPI_param, trim(adjustl(error_part))//' {Fnn_temp_z}', Fnn_temp_z) ! module "MPI_subroutines"
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"

#else ! use OpenMP instead
   !$omp parallel private(i)
   !$omp do
   do i = 1, m ! ensure WF normalization to 1
      Norm1(i) = DSQRT(SUM( Ha(i,:) * Ha(i,:) ))
!       Norm1(i) = DSQRT(SUM( Ha(:,i) * Ha(:,i) ))
   enddo
   !$omp end do
   !$omp do
   do i = 1, m
      PBx(:,i) = PBx(:,i) / Norm1(i)
      PBy(:,i) = PBy(:,i) / Norm1(i)
      PBz(:,i) = PBz(:,i) / Norm1(i)
   enddo ! i
   !$omp end do
   !$omp end parallel

   ! Do it with mkl:
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Ha, m, PBx, m, 0.0d0, Fnn_temp_x, m)
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Ha, m, PBy, m, 0.0d0, Fnn_temp_y, m)
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Ha, m, PBz, m, 0.0d0, Fnn_temp_z, m)
   CALL dgemm ('T','N', m, m, m, 1.0d0, PBx, m, Ha, m, 0.0d0, Fnn_temp_x, m)
   CALL dgemm ('T','N', m, m, m, 1.0d0, PBy, m, Ha, m, 0.0d0, Fnn_temp_y, m)
   CALL dgemm ('T','N', m, m, m, 1.0d0, PBz, m, Ha, m, 0.0d0, Fnn_temp_z, m)

   !$omp parallel private(i)
   !$omp do
   do i = 1, m
      Fnn_temp_x(:,i) = Fnn_temp_x(:,i) / Norm1(i)
      Fnn_temp_y(:,i) = Fnn_temp_y(:,i) / Norm1(i)
      Fnn_temp_z(:,i) = Fnn_temp_z(:,i) / Norm1(i)
   enddo ! i
   !$omp end do
   !$omp end parallel
#endif

   ! Do whole arrays at once:
   Fnnx = Fnn_temp_x*Fnn_temp_x
   Fnny = Fnn_temp_y*Fnn_temp_y
   Fnnz = Fnn_temp_z*Fnn_temp_z
   
   if (present(Fnnxy)) then ! do off-diagonal elements
      Fnnxy = Fnn_temp_x*Fnn_temp_y
      Fnnxz = Fnn_temp_x*Fnn_temp_z
      Fnnyx = Fnnxy
      Fnnyz = Fnn_temp_y*Fnn_temp_z
      Fnnzx = Fnnxz
      Fnnzy = Fnnyz
   endif
   
!    do i = 1, m
!       do j = 1, m
!          if (Fnnx(i,j) > 1d10) print*, 'get_Fnn', i,j, Fnnx(i,j), Fnn_temp_x(i,j)
!       enddo
!    enddo
!    print*,  'get_Fnn', MINVAL(Fnnx), MAXVAL(Fnnx)
!    print*,  'get_Fnn', MINVAL(Fnny), MAXVAL(Fnny)
!    print*,  'get_Fnn', MINVAL(Fnnz), MAXVAL(Fnnz)
   
   deallocate(Fnn_temp_x, Fnn_temp_y, Fnn_temp_z, PBx, PBy, PBz)
end subroutine get_Fnn


subroutine save_Eps_hw(Eps_hw, i, hw, Re_eps, Im_eps, LF, R, T, A, n, k, Eps_xx, Eps_yy, Eps_zz, dc_cond)
   real(8), dimension(:,:), intent(inout) :: Eps_hw
   integer, intent(in) :: i
   real(8), intent(in) :: hw, Re_eps, Im_eps, LF, R, T, A, n, k
   !complex(8), intent(in) :: Eps_xx, Eps_yy, Eps_zz
   complex, intent(in) :: Eps_xx, Eps_yy, Eps_zz
   real(8), intent(in), optional :: dc_cond
   Eps_hw(1,i) = hw  ! energy
   Eps_hw(2,i) = Re_eps ! real part of CDF
   Eps_hw(3,i) = Im_eps ! imaginary part of CDF
   Eps_hw(4,i) = LF  ! loss function
   Eps_hw(5,i) = R   ! reflectivity
   Eps_hw(6,i) = T   ! transmission
   Eps_hw(7,i) = A   ! absorption
   Eps_hw(8,i) = n   ! optical n
   Eps_hw(9,i) = k   ! optical k
   if (present(dc_cond)) Eps_hw(10,i) = dc_cond    ! cd-conductivity
   Eps_hw(11,i) = dble(Eps_xx)
   Eps_hw(12,i) = aimag(Eps_xx)
   Eps_hw(13,i) = dble(Eps_yy)
   Eps_hw(14,i) = aimag(Eps_yy)
   Eps_hw(15,i) = dble(Eps_zz)
   Eps_hw(16,i) = aimag(Eps_zz)
end subroutine save_Eps_hw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Drude subroutines:
subroutine get_drude(numpar, Scell, NSC)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   !===============================================
   real(8) :: wp, wph, me_eff, mh_eff, wl, vv, tth

   if (numpar%do_drude) then
      if (Scell(NSC)%eps%ReEps0 .LE. 0.0d0) call get_Re_and_Im(Scell(NSC)%eps%ReEps0, Scell(NSC)%eps%ImEps0, Scell(NSC)%eps%n, Scell(NSC)%eps%k)

      me_eff = Scell(NSC)%eps%me_eff
      mh_eff = Scell(NSC)%eps%mh_eff

      wl = Scell(NSC)%eps%w
      vv = 1.0d15/Scell(NSC)%eps%tau_e	! [1/fs] frequency of electron scattering
      if (Scell(NSC)%eps%tau_h .NE. Scell(NSC)%eps%tau_e) then
         tth = 1.0d15/Scell(NSC)%eps%tau_e ! [1/sec] collisional frequency for holes
      else
         tth = vv/mh_eff*me_eff    ! [1/sec] collisional frequency for holes
      endif
      wp = w_plasma(Scell(NSC)%Ne_CB, Scell(NSC)%V, me_eff) ! plasma frequency, function below
      wph = w_plasma(Scell(NSC)%Ne_CB, Scell(NSC)%V, mh_eff) ! plasma frequency for holes, function below

      ! real part of the dielectric constant:
      Scell(NSC)%eps%ReEps = Scell(NSC)%eps%ReEps0 - (wp/wl)*(wp/wl)/(1.0d0+vv*vv/(wl*wl)) - (wph/wl)*(wph/wl)/(1.0d0+tth*tth/(wl*wl)) 
      ! imaginary part of the dielectric constant
      Scell(NSC)%eps%ImEps = Scell(NSC)%eps%ImEps0 + (wp/wl)*(wp/wl)*vv/wl/(1.0d0+vv*vv/(wl*wl)) + (wph/wl)*(wph/wl)*tth/wl/(1.0d0+tth*tth/(wl*wl))

      call get_RTA_from_CDF(Scell(NSC)%eps%ReEps, Scell(NSC)%eps%ImEps, Scell(NSC)%eps%l, Scell(NSC)%eps%dd, Scell(NSC)%eps%teta, numpar%drude_ray, Scell(NSC)%eps%R, Scell(NSC)%eps%T, Scell(NSC)%eps%A)
      call get_n_k(Scell(NSC)%eps%ReEps, Scell(NSC)%eps%ImEps, Scell(NSC)%eps%n, Scell(NSC)%eps%k)	! transfer Re(e) and Im(e) into n and k
   endif
end subroutine get_drude




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kramers-Kronig relations for CDF:
subroutine Kramers_Kronig_Re_from_Im(hw, Im_CDF, Re_CDF)
   real(8), dimension(:), intent(in), target :: hw			! energy grid
   real(8), dimension(:), intent(in) :: Im_CDF		! imaginary part of the complex dielectric function
   real(8), dimension(:), intent(out), allocatable :: Re_CDF	! real part of the complex dielectric function
   !-------------------------------------------------------------
   integer :: Nsiz, i, j
   real(8) :: dw
   real(8), pointer :: w
   Nsiz = size(hw)
   if (.not.allocated(Re_CDF)) allocate(Re_CDF(Nsiz))
   
   ! Integrate:
   Re_CDF = 0.0d0	! to start with
   do i = 1, Nsiz
      w => hw(i)
      do j = 1, Nsiz
         if (j == 1) then	! start with
            if ( (hw(j+1) >= w) .and. (hw(j) <= w) ) then		! Cauchy principal value
               ! just exclude the special point
            else
               dw = hw(j)
               Re_CDF(i) = Re_CDF(i) + ( Im_CDF(j)*hw(j)/(hw(j)*hw(j) - w*w)*2.0d0 - Im_CDF(j+1)*hw(j+1)/(hw(j+1)*hw(j+1) - w*w) ) * 0.5d0 * dw
            endif
         else
            if ( (hw(j) >= w) .and. (hw(j-1) <= w) ) then		! Cauchy principal value
               ! just exclude the special point
            else
               dw = hw(j) - hw(j-1)
               Re_CDF(i) = Re_CDF(i) + ( Im_CDF(j)*hw(j)/(hw(j)*hw(j) - w*w) + Im_CDF(j)*hw(j-1)/(hw(j-1)*hw(j-1) - w*w) ) * 0.5d0 * dw
            endif
         endif
      enddo ! j
   enddo ! i
   Re_CDF = 1.0d0 + 2.0d0/g_Pi * Re_CDF
   
!    do i = 1, Nsiz
!       print*, 'KK:',  hw(i), Re_CDF(i)
!    enddo
!    pause 'Kramers_Kronig_Re_from_Im'
   
   nullify(w)
end subroutine Kramers_Kronig_Re_from_Im




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! General optical coefficients subroutines:
subroutine get_RTA_from_CDF(ReEps, ImEps, ll, dd, teta0, drude_ray, R, T, A)  ! Ref. [3]
   real(8), intent(in) :: ReEps, ImEps	! Re(e) and Im(e) of dielectric function
   real(8), intent(in) :: ll  ! [nm] probe-pulse wavelength
   real(8), intent(in) :: dd  ! [nm] experimental layer thickness
   real(8), intent(in) :: teta0 ! angle of the probe pulse [radians]
   integer, intent(in) :: drude_ray	! over how many rays to sum up: (1)=1 or (>1)=all
   real(8), intent(out) :: R, T, A	! reflection, transmittion and absorption coefficients
   complex :: cxn1, cxn2, ckax, ckaz
   complex cs(5), cr12, ct12, cr23, ct23, cxki, cd, cref, ctran
   real(8) :: A2, B2, n, k, alpha, temp
   real(8) :: sai, xkaz, xkax, nm1, np1, exp2ak, dak, dan

   call get_n_k(ReEps, ImEps, n, k)	! transfer Re(e) and Im(e) into n and k
   alpha = 2.0d0*g_Pi*dd/ll
   nm1 = n - 1.0d0
   np1 = n + 1.0d0
   dak = 2.0d0*alpha*k
   if (dak .GT. 300) then
      exp2ak = 1.0d200
   else
      exp2ak = exp(dak)
   endif

   dan = 2.0d0*alpha*n
!    A2 = ((nm1*(1.0d0-exp2ak*cos(dan)-k*exp2ak*sin(dan)))**2+(k*(1.0d0-exp2ak)*cos(dan)+nm1*sin(dan)*exp2ak)**2)/(nm1**2+k**2)
!    B2 = (((np1*np1+k*k)**2-exp2ak*(((n*n+k*k-1.0d0)**2-4.0d0*k*k)*cos(dan)+4.0d0*k*((n*n-1.0d0)+k*k)*sin(dan)))**2 + (exp2ak*(cos(dan)*(-4.0d0*k*(n*n-1.0d0+k*k))+((n*n+k*k-1.0d0)**2-4.0d0*k*k)*sin(dan)))**2 )/(np1*np1+k*k)**4

   cxn1 = cmplx(n,-k)
   sai   = 1.0d0	! air
   cxn2 = cmplx(sai,0.0d0) ! air
   ! Eq (4.4.3) from "multilayer book", p.97:
   temp = 2.0d0*g_Pi*sai/ll
   xkaz = temp*sin(teta0)
   xkax = temp*cos(teta0)
   ckax = cmplx(xkax,0.0d0)
   ckaz = cmplx(xkaz,0.0d0)
   cxki  = 2.0d0*g_Pi*cxn1/ll
   cs(3) = sqrt(cxki*cxki-ckaz*ckaz)

   cd    = cmplx(dd,0.0d0)
   cs(2) = cphase(cs(3),cd) ! function below
   ckaz  = cmplx(xkaz,0.0d0)
   cxki  = 2.0d0*g_Pi*cxn2/ll
   cs(5) = sqrt(cxki*cxki-ckaz*ckaz)

   cr12 = crefl(ckax,cs(3))
   ct12 = ctrans(ckax,cs(3))

   cr23 = crefl(cs(3),cs(5))
   ct23 = ctrans(cs(3),cs(5))
   if (drude_ray .GT. 1) then  ! total summed up:
      cref = (cr12+cr23*cs(2)*cs(2))/(1.0d0+cr12*cr23*cs(2)*cs(2))
      ctran = (ct12*ct23*cs(2))/(1.0d0+cr12*cr23*cs(2)*cs(2))
      R = cabs(cref)*cabs(cref)
      T = cabs(ctran)*cabs(ctran)*cs(5)/ckax
      A = 1.0d0 - R - T
   else ! first ray propagation:
      R = cabs(cr12)*cabs(cr12)
      T = cabs(ct12*ct23*cs(2))
      T = T*T*cs(5)/ckax
      A = 1.0d0 - R - T
   endif
end subroutine get_RTA_from_CDF


function w_plasma(Ne, V, me_eff)   ! plasma frequency [1/s]
   real(8), intent(in) :: Ne, V ! number of electrons in the supercell; volume of the supercell [anstroms]
   real(8), intent(in), optional :: me_eff	! effective electron mass [kg]
   real(8) :: w_plasma
   if (present(me_eff)) then
      w_plasma = SQRT((Ne/(V*1d-30))/(g_e0)*g_e*g_e/me_eff) ! plasma frequency for effective electron mass
   else
      w_plasma = SQRT((Ne/(V*1d-30))/(g_e0)*g_e*g_e/g_me)   ! plasma frequency of a free-electron
   endif
end function w_plasma


subroutine get_n_k(ReEps, ImEps, n, k)
   real(8), intent(in) :: ReEps, ImEps	! real and imagenary parts of dialectric function 
   real(8), intent(out) :: n, k	! n and k optical coefficients
   !------------------
   complex(8) :: CDF, sqrtCDF
   CDF = cmplx(ReEps,ImEps)
   sqrtCDF = sqrt(CDF)
   n =  dble( sqrtCDF )
   k = aimag( sqrtCDF )
end subroutine get_n_k


subroutine get_n_k_old(ReEps, ImEps, n, k)
   real(8), intent(in) :: ReEps, ImEps	! real and imagenary parts of dialectric function
   real(8), intent(out) :: n, k	! n and k optical coefficients
   n = sqrt((ReEps + sqrt(ReEps*ReEps+ImEps*ImEps))*0.5d0)
   k =sqrt((-ReEps + sqrt(ReEps*ReEps+ImEps*ImEps))*0.5d0)
end subroutine get_n_k_old


subroutine get_Re_and_Im(ReEps, ImEps, n, k)
   real(8), intent(out) :: ReEps, ImEps	! real and imagenary parts of dialectric function 
   real(8), intent(in) :: n, k	! n and k optical coefficients
   ReEps = n*n - k*k
   ImEps = 2.0d0*n*k
end subroutine get_Re_and_Im


! reflection coefficient between two layers as function of two wave vectors:
function crefl(cxki,cxkj)  ! Ref. [3]
   complex :: crefl, cxki, cxkj
   crefl = (cxki-cxkj)/(cxki+cxkj)
end function crefl

! transmission coefficient between two layers as function of two wave vectors:
function ctrans(cxki,cxkj) ! Ref. [3]
   complex :: ctrans, cxki, cxkj
   ctrans = 2.0d0*cxki/(cxki+cxkj)
end function ctrans

! phase factor within the layer i:
function cphase(cxkx,cd)   ! Ref. [3]
   complex :: cphase, cxkx, cd
   complex cc, carg
   cc = (0.0d0,1.0d0)
   carg = -cc*cxkx*cd
   cphase=exp(carg)
end function cphase


END MODULE Optical_parameters
