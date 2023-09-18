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
! This module assembles subroutines that require complex Hamiltonian in TB

MODULE TB_complex
use Universal_constants
use Objects
use TB, only : k_point_choice, construct_complex_Hamiltonian
use Optical_parameters, only : allocate_Eps_hw, get_Onsager_coeffs, get_Kubo_Greenwood_CDF, get_kappa_e_e, get_Onsager_dynamic
use Electron_tools, only : get_DOS_sort
use Little_subroutines, only : Find_in_array_monoton, linear_interpolation

USE OMP_LIB, only : OMP_GET_THREAD_NUM

implicit none
PRIVATE

real(8), parameter :: m_gamm = 1.5d14  ! [1/s] gamma parameter

public :: use_complex_Hamiltonian


 contains


subroutine use_complex_Hamiltonian(numpar, matter, Scell, NSC, Err)  ! From Ref. [2]
   type (Numerics_param), intent(inout) :: numpar  ! numerical parameters, including drude-function
   type(Solid), intent(inout) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell   ! supercell with all the atoms as one object
   integer, intent(in) :: NSC    ! number of supercell
   type(Error_handling), intent(inout) :: Err   ! error save
   !--------------------
   !complex(8), dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   !complex(8), dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   complex, dimension(:,:), allocatable :: cPRRx, cPRRy, cPRRz  ! effective momentum operators
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
   real(8) :: w, kx, ky, kz
   real(8), dimension(:), allocatable :: w_grid
   integer :: i, j, N, FN, ix, iy, iz, ixm, iym, izm, schem, Ngp, Nsiz, N_wgrid
   real(8), dimension(:), allocatable :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable :: Eps_hw ! array of all eps vs hw
   real(8), dimension(:,:), allocatable :: Eps_hw_temp ! array of all eps vs hw
   real(8), dimension(:), allocatable :: kappa, kappa_ee, kappa_temp, kappa_ee_temp    ! electron heat conductivity vs Te
   real(8), dimension(:), allocatable :: kappa_mu_grid, kappa_mu_grid_temp, kappa_Ce_grid, kappa_Ce_grid_temp
   integer :: Nsiz_Te, Nsiz_DOS_1, Nsiz_DOS_2, Nsiz_DOS_3
   real(8) :: Te_min, Te_max, dTe, kap_temp
   real(8), dimension(:,:), allocatable :: DOS, DOS_temp  ! [eV] grid; [a.u.] DOS
   real(8), dimension(:,:,:), allocatable :: DOS_partial, DOS_partial_temp ! partial DOS made of each orbital type
   logical :: anything_to_do

   !-----------------------------------------------
   ! Check if there is anything to do with the complex Hamiltonian:
   anything_to_do = (numpar%save_DOS .or. numpar%do_kappa .or. numpar%do_kappa_dyn .or. &
                     (numpar%optic_model == 4) .or. (numpar%optic_model == 5))
   if (.not.anything_to_do) return  ! nothing to do, exit


   !-----------------------------------------------
   ! Define the total number of k-points:
   call create_DOS_arrays(numpar, Scell(NSC), matter, DOS, DOS_partial)  ! below
   Nsiz_DOS_1 = size(DOS_partial,1)
   Nsiz_DOS_2 = size(DOS_partial,2)
   Nsiz_DOS_3 = size(DOS_partial,3)

   !-----------------------------------------------
   ! For electron-temperature-dependent quantities:
   ! Define electron temperature grid:
   call set_temperature_grid(Scell(NSC), numpar, kappa, kappa_ee, kappa_mu_grid, kappa_Ce_grid) ! below
   Nsiz_Te = size(kappa)

   !-----------------------------------------------
   ! Allocate the array of optical coefficients spectrum (if not allocated before):
   ! Define the grid size:
   call set_frequency_grid(Scell(NSC), w_grid, Eps_hw)  ! below
   N_wgrid = size(w_grid)

   !-----------------------------------------------
   ! Define the total number of k-points:
   call get_total_num_of_k_points(numpar, ixm, iym, izm, schem, Nsiz)   ! below

   !-----------------------------------------------
   ! Calculate what's required for all k-points:
   !$omp PARALLEL private(ix, iy, iz, Ngp, kx, ky, kz, cPRRx, cPRRy, cPRRz, CHij, Ei, Eps_hw_temp, &
   !$omp                  kappa_temp, kappa_ee_temp, kappa_mu_grid_temp, kappa_Ce_grid_temp, DOS_temp, DOS_partial_temp)
   if (.not.allocated(Eps_hw_temp)) allocate(Eps_hw_temp(16,N_wgrid), source = 0.0d0) ! all are there
   if (.not.allocated(kappa_temp)) allocate(kappa_temp(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_ee_temp)) allocate(kappa_ee_temp(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_mu_grid_temp)) allocate(kappa_mu_grid_temp(Nsiz_Te), source = 0.0d0)   ! mu [eV]
   if (.not.allocated(kappa_Ce_grid_temp)) allocate(kappa_Ce_grid_temp(Nsiz_Te), source = 0.0d0)   ! Ce
   if (.not.allocated(DOS_temp)) allocate(DOS_temp(2,Nsiz_DOS_3), source = DOS)   ! DOS
   if (.not.allocated(DOS_partial_temp)) allocate(DOS_partial_temp(Nsiz_DOS_1,Nsiz_DOS_2,Nsiz_DOS_3), source = 0.0d0)   ! DOS_partial
   !$omp do schedule(dynamic) reduction( + : Eps_hw, kappa, kappa_ee, kappa_mu_grid, kappa_Ce_grid, DOS, DOS_partial)
   do Ngp = 1, Nsiz
      ! Split total index into 3 coordinates indices:
      ix = ceiling( dble(Ngp)/dble(iym*izm) )
      iy = ceiling( dble(Ngp - (ix-1)*iym*izm)/dble(izm) )
      iz = Ngp - (ix-1)*iym*izm - (iy-1)*izm

      !-------------------------------
      ! k-points:
      call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid) ! module "TB"
      if (numpar%verbose) write(*,'(a,i4,a,i6,i3,i3,i3,f9.4,f9.4,f9.4,a)') 'Thread #', OMP_GET_THREAD_NUM(), &
                                     ' point #', Ngp, ix, iy, iz, kx, ky, kz, ' k-points'

      !-------------------------------
      ! Get the parameters of the complex Hamiltonian:
      ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
         select type(ARRAY)
         type is (TB_H_Pettifor)    ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_Molteni)     ! orthogonal
            call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij, Ei, kx, ky, kz, &
            cPRRx=cPRRx, cPRRy=cPRRy, cPRRz=cPRRz)  ! module "TB"
         type is (TB_H_Fu)          ! orthogonal
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
      ! Get DOS:
      if (numpar%save_DOS) then ! if required
         call get_DOS_on_k_points(numpar, Ei, CHij, DOS_temp, DOS_partial_temp) ! below
      else  ! skip DOS calculations
         DOS_temp = 0.0d0
         DOS_partial_temp = 0.0d0
      endif
      ! Save DOS data:
      DOS(2,:) = DOS(2,:) + DOS_temp(2,:)
      DOS_partial = DOS_partial + DOS_partial_temp

      !-------------------------------
      ! Get the parameters of the CDF:
      if ((numpar%optic_model == 4) .or. (numpar%optic_model == 5)) then ! if requested
         call get_Kubo_Greenwood_CDF(numpar, Scell, NSC, w_grid, cPRRx, cPRRy, cPRRz, Ei, Eps_hw_temp)   ! module "Optical_parameters"
      else  ! skip it, if not reqired
         Eps_hw_temp = 0.0d0
      endif
      ! Save CDF data:
      Eps_hw = Eps_hw + Eps_hw_temp ! sum data at different k-points

      !-------------------------------
      ! If required, do Onsager coefficients (for electronic heat conductivity):
      if (numpar%do_kappa .or. numpar%do_kappa_dyn) then  ! if requested
         ! Electron-phonon contribution:
         if (numpar%do_kappa) then ! (static, KG)
            call get_Onsager_coeffs(numpar, matter, Scell, NSC, cPRRx, cPRRy, cPRRz, Ei, kappa_temp, &
                              Scell(NSC)%kappa_Te_grid, kappa_mu_grid_temp, kappa_Ce_grid_temp)   ! module "Optical_parameters"

            ! Contribution of the electronic term:
            call get_kappa_e_e(numpar, matter, Scell, NSC, Ei, kappa_mu_grid_temp, &
                              Scell(NSC)%kappa_Te_grid, kappa_ee_temp) ! module "Optical_parameters"
         endif

         if (numpar%do_kappa_dyn) then ! (dynamic)
            ! Only gamma-point calculations, so only 1 time:
            if (Ngp == 1) then
               call get_Onsager_dynamic(numpar, matter, Scell, NSC, kappa_temp) ! module "Optical_parameters"

               ! Contribution of the electronic term:
               call get_kappa_e_e(numpar, matter, Scell, NSC, Scell(NSC)%Ei, kappa_mu_grid_temp, &
                              Scell(NSC)%kappa_Te_grid, kappa_ee_temp) ! module "Optical_parameters"

               !print*, 'kappa_temp', kappa_temp
               !print*, 'kappa_ee_temp', kappa_ee_temp
            else
               kappa_temp = 0.0d0
               kappa_ee_temp = 0.0d0
            endif
         endif

      else  ! if not required
         kappa_temp = 0.0d0
         kappa_mu_grid_temp = 0.0d0
         kappa_Ce_grid_temp = 0.0d0
         kappa_ee_temp = 0.0d0
      endif

      !-------------------------------
      ! Save kappa - electronic heat conductivity data:
      kappa = kappa + kappa_temp ! sum up at different k-points
      kappa_ee = kappa_ee + kappa_ee_temp ! sum up at different k-points
      kappa_mu_grid = kappa_mu_grid + kappa_mu_grid_temp    ! average mu
      kappa_Ce_grid = kappa_Ce_grid + kappa_Ce_grid_temp    ! average Ce
   enddo ! Ngp
   !$omp end do
   if (allocated(Eps_hw_temp)) deallocate(Eps_hw_temp)
   if (allocated(kappa_temp)) deallocate(kappa_temp)
   if (allocated(kappa_ee_temp)) deallocate(kappa_ee_temp)
   if (allocated(kappa_mu_grid_temp)) deallocate(kappa_mu_grid_temp)
   if (allocated(kappa_Ce_grid_temp)) deallocate(kappa_Ce_grid_temp)
   if (allocated(DOS_temp)) deallocate(DOS_temp)
   if (allocated(DOS_partial_temp)) deallocate(DOS_partial_temp)
   !$omp end parallel

   !-----------------------------------------------
   ! Save the k-point averages:
   !-------------------------------
   ! CDF:
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

   !-------------------------------
   ! Save electron heat conductivity, averaged over k-points:
   if (numpar%do_kappa) then  ! (static, if requested)
      Scell(NSC)%kappa_e = kappa(1)/dble(Nsiz)  ! transient temperature
      Scell(NSC)%kappa_e_vs_Te = kappa/dble(Nsiz)       ! e-ph array vs Te
      Scell(NSC)%kappa_ee_vs_Te = kappa_ee/dble(Nsiz)   ! e-e array vs Te
      Scell(NSC)%kappa_mu_grid = kappa_mu_grid/dble(Nsiz)  ! array of mu vs Te
      Scell(NSC)%kappa_Ce_grid = kappa_Ce_grid/dble(Nsiz)  ! array of Ce vs Te
   endif
   if (numpar%do_kappa_dyn) then ! (dynamic)
      Scell(NSC)%kappa_e_vs_Te = kappa       ! e-ph
      Scell(NSC)%kappa_ee_vs_Te = kappa_ee   ! e-e
      Scell(NSC)%kappa_e = 0.0d0
      if (kappa(1) > 0.0d0) Scell(NSC)%kappa_e = Scell(NSC)%kappa_e + 1.0d0/kappa(1)
      ! e-e term for the given temperature:
      call Find_in_array_monoton(Scell(NSC)%kappa_Te_grid, Scell(NSC)%Te, i)  ! module "Little_subroutines"
      call linear_interpolation(Scell(NSC)%kappa_Te_grid, kappa_ee, Scell(NSC)%Te, kap_temp, i) ! module "Little_subroutines"
      ! Save it to print out:
      Scell(NSC)%kappa_ee_vs_Te(1) = kap_temp
      if (kap_temp > 0.0d0) Scell(NSC)%kappa_e = Scell(NSC)%kappa_e + 1.0d0/kap_temp
      ! Total kappa:
      if (Scell(NSC)%kappa_e > 0.0d0) Scell(NSC)%kappa_e = 1.0d0 / Scell(NSC)%kappa_e ! total
      !print*, 'use_complex_Hamiltonian:', kappa(1), kap_temp, Scell(NSC)%kappa_e
   endif

   !-------------------------------
   ! Save DOS:
   if (numpar%save_DOS) then
      Scell(NSC)%DOS(2,:) = DOS(2,:)/dble(Nsiz)
      select case (numpar%DOS_splitting)
      case (1)
         Scell(NSC)%partial_DOS = DOS_partial/dble(Nsiz)
      end select
   endif

   !-------------------------------
   ! Clean up:
   if (allocated(cPRRx)) deallocate(cPRRx)
   if (allocated(cPRRy)) deallocate(cPRRy)
   if (allocated(cPRRz)) deallocate(cPRRz)
   if (allocated(w_grid)) deallocate(w_grid)
   if (allocated(Ei)) deallocate(Ei)
   if (allocated(CHij)) deallocate(CHij)
   if (allocated(Eps_hw)) deallocate(Eps_hw)
   if (allocated(kappa)) deallocate(kappa)
   if (allocated(kappa_ee)) deallocate(kappa_ee)
   if (allocated(kappa_mu_grid)) deallocate(kappa_mu_grid)
   if (allocated(kappa_Ce_grid)) deallocate(kappa_Ce_grid)
   if (allocated(DOS)) deallocate(DOS)
   if (allocated(DOS_partial)) deallocate(DOS_partial)
end subroutine use_complex_Hamiltonian



subroutine get_DOS_on_k_points(numpar, Ei, CHij, DOS_temp, DOS_partial_temp)
   type (Numerics_param), intent(in) :: numpar  ! numerical parameters
   real(8), dimension(:), intent(in) :: Ei  ! [eV] energy levels
   complex, dimension(:,:), intent(in) :: CHij ! eigenvectors of the hamiltonian
   real(8), dimension(:,:), intent(inout) :: DOS_temp  ! [eV] grid; [a.u.] DOS
   real(8), dimension(:,:,:), intent(inout) :: DOS_partial_temp ! partial DOS made of each orbital type
   !-------------------------

   call get_DOS_sort(Ei, DOS_temp, numpar%Smear_DOS, DOS_partial_temp, numpar%mask_DOS, CHij = CHij)  ! module "Electron_tools"
end subroutine get_DOS_on_k_points



subroutine create_DOS_arrays(numpar, Scell, matter, DOS, DOS_partial)
   type (Numerics_param), intent(in) :: numpar  ! numerical parameters
   type (Super_cell), intent(inout) :: Scell   ! supercell with all the atoms as one object
   type(Solid), intent(in) :: matter     ! material parameters
   real(8), dimension(:,:), allocatable, intent(inout) :: DOS  ! [eV] grid; [a.u.] DOS
   real(8), dimension(:,:,:), allocatable, intent(inout) :: DOS_partial ! partial DOS made of each orbital type
   !----------------------
   real(8) :: dE, Emax, Estart
   integer :: Ei_siz, Nsiz, n_types, i

   ! Set grid for DOS:
   if (.not.allocated(DOS)) then
      dE = 0.1d0  ! [eV] uniform energy grid step for DOS
      Ei_siz = size(Scell%Ei) ! the uppermost energy level at the start
      Estart = dble(FLOOR(Scell%Ei(1) - 50.0d0*numpar%Smear_DOS))  ! [eV]
      Emax = min(Scell%Ei(Ei_siz) + 50.0d0*numpar%Smear_DOS,100.0)   ! no need to trace levels higher than 100 eV
      Nsiz = CEILING( (Emax - Estart)/dE )
      allocate(DOS(2,Nsiz), source = 0.0d0)
      do i = 1, Nsiz ! set energy grid [eV]
         DOS(1,i) = Estart + dE*dble(i-1)
      enddo

      if (.not.allocated(Scell%DOS)) then ! it's the first time, set it:
         allocate(Scell%DOS(2,Nsiz), source = 0.0d0)
         Scell%DOS(1,:) = DOS(1,:)  ! save energy grid
      endif

      n_types = size(numpar%mask_DOS,2)
      if (.not.allocated(DOS_partial)) then
         allocate(DOS_partial(matter%N_KAO, n_types, Nsiz), source = 0.0d0)
      endif
      ! Partial DOS if needed:
      select case (numpar%DOS_splitting)
      case (1)
         if (.not. allocated(Scell%partial_DOS)) then
            allocate(Scell%partial_DOS(matter%N_KAO, n_types, Nsiz), source = 0.0d0)
         endif
      case default
         ! No need to sort DOS per orbitals
      endselect

   endif
end subroutine create_DOS_arrays



subroutine get_total_num_of_k_points(numpar, ixm, iym, izm, schem, Nsiz)
   type (Numerics_param), intent(in) :: numpar  ! numerical parameters
   integer, intent(out) :: ixm, iym, izm, schem, Nsiz
   !----------------------
   ! For the user-defind number of k-points along axis::
   ixm = numpar%ixm  ! x
   iym = numpar%iym  ! y
   izm = numpar%izm  ! z
   if (allocated(numpar%k_grid)) then
      schem = 1   ! user-defined grid is present
      Nsiz = size(numpar%k_grid,1)  ! size of the user provided grid
   else
      schem = 0   ! no user-defined grid, use default
      Nsiz = ixm*iym*izm
   endif
end subroutine get_total_num_of_k_points



subroutine set_frequency_grid(Scell, w_grid, Eps_hw)
   type (Super_cell), intent(inout) :: Scell   ! supercell with all the atoms as one object
   real(8), dimension(:), allocatable :: w_grid
   real(8), dimension(:,:), allocatable :: Eps_hw
   !------------------
   integer :: N, i

   if (Scell%eps%all_w) then ! full spectrum:
      call allocate_Eps_hw(Scell%eps%E_min, Scell%eps%E_max, Scell%eps%dE, Scell%eps%Eps_hw) ! module "Optical_parameters"
      N = size(Scell%eps%Eps_hw,2) ! use this size to define temporary arrays within this subroutine
      ! Grid of frequencies:
      allocate(w_grid(N), source = 0.0d0)
      w_grid(1) =  Scell%eps%E_min*g_e/g_h ! [1/s] frequency starting point for the optical spectrum (from [eV])
      do i = 2, N
         w_grid(i) = w_grid(i-1) + Scell%eps%dE*g_e/g_h  ! [1/s] frequency next grid point [eV]
      enddo
   else  ! single frequency
      N = 1
      allocate(w_grid(N), source = 0.0d0)
      if (.not.allocated(Scell%eps%Eps_hw)) then
         allocate(Scell%eps%Eps_hw(16,N)) ! different parameters saved for all energy grid points
      endif
      w_grid(1) = Scell%eps%w
   endif
   if (.not.allocated(Eps_hw)) allocate(Eps_hw(16,N), source = 0.0d0) ! all are there
end subroutine set_frequency_grid



subroutine set_temperature_grid(Scell, numpar, kappa, kappa_ee, kappa_mu_grid, kappa_Ce_grid)
   type (Super_cell), intent(inout) :: Scell   ! supercell with all the atoms as one object
   type (Numerics_param), intent(in) :: numpar  ! numerical parameters
   real(8), dimension(:), allocatable :: kappa, kappa_ee, kappa_mu_grid, kappa_Ce_grid
   !-----------------------------
   real(8) :: Te_min, Te_max, dTe
   integer :: Nsiz_Te, i

   ! Define temperature grid:
   Te_min = numpar%kappa_Te_min  ! [K]
   Te_max = numpar%kappa_Te_max  ! [K]
   dTe = numpar%kappa_dTe  ! [K]
   Nsiz_Te = (Te_max+dTe - Te_min)/dTe

   ! Allocate the arrays:
   if (.not.allocated(kappa)) allocate(kappa(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_ee)) allocate(kappa_ee(Nsiz_Te), source = 0.0d0)
   if (.not.allocated(kappa_mu_grid)) allocate(kappa_mu_grid(Nsiz_Te), source = 0.0d0)   ! mu [eV]
   if (.not.allocated(kappa_Ce_grid)) allocate(kappa_Ce_grid(Nsiz_Te), source = 0.0d0)   ! Ce [J/(m^3 K)]

   if (numpar%do_kappa .or. numpar%do_kappa_dyn) then
      if (.not.allocated(Scell%kappa_Te_grid)) allocate(Scell%kappa_Te_grid(Nsiz_Te), source = 0.0d0)   ! Te [K]
      if (.not.allocated(Scell%kappa_e_vs_Te)) allocate(Scell%kappa_e_vs_Te(Nsiz_Te), source = 0.0d0)   ! kappa
      if (.not.allocated(Scell%kappa_mu_grid)) allocate(Scell%kappa_mu_grid(Nsiz_Te), source = 0.0d0)   ! mu [eV]
      if (.not.allocated(Scell%kappa_Ce_grid)) allocate(Scell%kappa_Ce_grid(Nsiz_Te), source = 0.0d0)   ! Ce [J/(m^3 K)]
      ! Set the grid:
      Scell%kappa_Te_grid(1) = Te_min
      do i = 2, Nsiz_Te ! Set the electronic temperature grid
         Scell%kappa_Te_grid(i) = Scell%kappa_Te_grid(i-1) + dTe
      enddo
   endif
end subroutine set_temperature_grid


END MODULE TB_complex
