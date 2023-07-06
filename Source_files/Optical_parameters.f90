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
MODULE Optical_parameters
use Universal_constants
use Objects
use Algebra_tools, only : sym_diagonalize
use Electron_tools, only : get_number_of_CB_electrons, set_Fermi, set_Erf_distribution
use TB_Fu, only : Complex_Hamil_tot_F
use TB_Pettifor, only : Complex_Hamil_tot
use TB_Molteni, only : Complex_Hamil_tot_Molteni
use TB_NRL, only : Complex_Hamil_NRL
use TB_DFTB, only : Complex_Hamil_DFTB, identify_DFTB_orbitals_per_atom
use TB, only : k_point_choice, construct_complex_Hamiltonian
use Little_subroutines, only : deallocate_array

implicit none
PRIVATE

public :: get_optical_parameters


 contains


subroutine get_optical_parameters(numpar, matter, Scell, Err) ! optical coefficients, module "Optical_parameters"
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Solid), intent(inout) :: matter  ! Material parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   type(Error_handling), intent(inout) :: Err	! error save
   !=====================================
   integer NSC, i
   real(8), dimension(:), allocatable :: Re_CDF

   do NSC = 1, size(Scell) ! for all super-cells
      ! Get the number of CB electrons:
      call get_number_of_CB_electrons(Scell, NSC) ! module "Electron_tools"

      select case (numpar%optic_model)
         case (1)	! within the Drude model
            call get_drude(numpar, Scell, NSC)  ! below
         case (2)	! Trani et al. PRB 72, 075423 (2005) -- This subroutine is TB-parameterization specific:
            !call get_trani_all_complex(numpar, Scell, NSC, Scell(NSC)%fe, Scell(NSC)%eps%all_w, Err)    ! below
            call get_trani_all_complex(numpar, Scell, NSC, Scell(NSC)%eps%all_w, Err)    ! below
            
!             ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
!             select type(ARRAY)
!                type is (TB_H_Pettifor)
!                   call get_trani_all_complex(numpar, Scell, NSC, Scell(NSC)%fe, Scell(NSC)%eps%all_w, Err)
!                type is (TB_H_Molteni)
!                   call get_trani_all_complex_M(numpar, Scell, NSC, Scell(NSC)%MDatoms, ARRAY, Scell(NSC)%fe, Scell(NSC)%eps%all_w, Err)
!                type is (TB_H_Fu)
!                    call get_trani_all_complex_F(numpar, Scell, NSC, Scell(NSC)%MDatoms, ARRAY, Scell(NSC)%fe, Scell(NSC)%eps%all_w, Err)
!                type is (TB_H_NRL)
!                    call get_trani_all_complex_NRL(numpar, Scell, NSC, ARRAY, Scell(NSC)%fe, Scell(NSC)%eps%all_w, Err)
!                type is (TB_H_DFTB)
!                    call get_trani_all_complex_DFTB(numpar, Scell, NSC, Scell(NSC)%fe, Scell(NSC)%eps%all_w, Err)
!                end select	!  type of Hamiltonain parameterization
!             END ASSOCIATE
         case (3)	! Trani at the Gamma-point only
            call get_trani_all(numpar, Scell, NSC, Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%eps%all_w)   ! below
         case default ! no optical coefficients needed
            !call get_trani_all(numpar, Scell, NSC, Scell(NSC)%Ei, Scell(NSC)%Ha, Scell(NSC)%fe, Scell(NSC)%eps%all_w)
      end select ! (numpar%optic_model)
      
      !-------------------------------------
      ! Convergence with respect to the number of k-points is better for Im part than Re part of CDF, so use Kramers-Kronig relations to restore Re part:
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
            call save_Eps_hw(Scell(NSC)%eps%Eps_hw, i, Scell(NSC)%eps%Eps_hw(1,i), Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i), Scell(NSC)%eps%Eps_hw(4,i), Scell(NSC)%eps%Eps_hw(5,i), Scell(NSC)%eps%Eps_hw(6,i), Scell(NSC)%eps%Eps_hw(7,i), Scell(NSC)%eps%Eps_hw(8,i), Scell(NSC)%eps%Eps_hw(9,i), dcmplx(Scell(NSC)%eps%Eps_hw(11,i),Scell(NSC)%eps%Eps_hw(12,i)), dcmplx(Scell(NSC)%eps%Eps_hw(13,i),Scell(NSC)%eps%Eps_hw(14,i)), dcmplx(Scell(NSC)%eps%Eps_hw(15,i),Scell(NSC)%eps%Eps_hw(16,i)),  Scell(NSC)%eps%Eps_hw(10,i))
!             write(*,'(f,es,es)') Scell(NSC)%eps%Eps_hw(1,i), Scell(NSC)%eps%Eps_hw(2,i), Scell(NSC)%eps%Eps_hw(3,i)
         enddo
      endif
      !-------------------------------------
   enddo ! NSC = 1, size(Scell) ! for all super-cells
end subroutine get_optical_parameters



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complex-value Trani's subroutines:

!subroutine get_trani_all_complex(numpar, Scell, NSC, fe, all_w, Err)
subroutine get_trani_all_complex(numpar, Scell, NSC, all_w, Err)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
!    real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err	! error save

   real(8), dimension(:,:), allocatable :: Fnnx, Fnny, Fnnz
   real(8), dimension(:,:), allocatable :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), dimension(:,:), allocatable :: Eps
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
            write(*,'(i3,i3,i3,f9.4,f9.4,f9.4,a)') ix, iy, iz, kx, ky, kz, ' OPT'
            
            call get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
            ! With off-diagonal elements:
!             call get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnxx, Fnnyy, Fnnzz, kx=kx, ky=ky, kz=kz, Err=Err, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy) ! see below
!             call get_Fnn_complex(numpar, matter, atoms, TB_Hamil, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx=real(ix)/real(ixm), ky=real(iy)/real(iym), kz=real(iz)/real(izm), Err=Err) ! see below

            if (all_w) then ! full spectrum:
               w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] frequency starting point for the optical spectrum [eV]
               do i = 1, N
                  call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
                  ! With off-diagonal elements:
!                   call get_trani(numpar, Scell, NSC, Fnnxx, Fnnyy, Fnnzz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy)
                  !call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
                  call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy
                  w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
               enddo
               Eps_hw = Eps_hw + Eps_hw_temp
            else ! only for given probe:
               w = Scell(NSC)%eps%w
               call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
               ! With off-diagonal elements:
!                call get_trani(numpar, Scell, NSC, Fnnxx, Fnnyy, Fnnzz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond, Fnnxy=Fnnxy, Fnnxz=Fnnxz, Fnnyx=Fnnyx, Fnnyz=Fnnyz, Fnnzx=Fnnzx, Fnnzy=Fnnzy)
               !call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
               call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy   
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
     
     Scell(NSC)%eps%ReEps = Eps_hw(2,1)	! real part of CDF
     Scell(NSC)%eps%ImEps = Eps_hw(3,1)	! imaginary part of CDF
     Scell(NSC)%eps%R = Eps_hw(5,1)	! reflectivity
     Scell(NSC)%eps%T = Eps_hw(6,1)	! transmission
     Scell(NSC)%eps%A = Eps_hw(7,1)	! absorption
     Scell(NSC)%eps%n = Eps_hw(8,1)	! optical n
     Scell(NSC)%eps%k = Eps_hw(9,1)	! optical k
     Scell(NSC)%eps%dc_cond = Eps_hw(10,1)	! dc-conductivity
     Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,1), Eps_hw(12,1))  ! Re_E_xx and Im_E_xx
     Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,1), Eps_hw(14,1))  ! Re_E_yy and Im_E_yy
     Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,1), Eps_hw(16,1))  ! Re_E_zz and Im_E_zz


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


subroutine get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx ,ky, kz, Err, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy)
!subroutine get_Fnn_complex(numpar, matter, atoms, TB, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx ,ky, kz, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
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
subroutine get_trani_all(numpar, Scell, NSC, Ei, Ha, fe, all_w)
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
   
   call get_Fnn(Scell, NSC, Ha, Ei, Fnnx, Fnny, Fnnz) ! see below
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



subroutine get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, n, k, dc_cond, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy)
! subroutine get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, n, k)
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
   NEL = size(Ei)
   
!    allocate(Re_array(NEL,3))
!    allocate(Im_array(NEL,3))
   
!    g = 1.0d13	! [1/s] gamma parameter
   g = 1.5d14	! [1/s] gamma parameter

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


subroutine get_Fnn(Scell, NSC, Ha, Ei, Fnnx, Fnny, Fnnz, Fnnxy, Fnnxz, Fnnyx, Fnnyz, Fnnzx, Fnnzy)
! subroutine get_Fnn(matter, Ha, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz)
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
   m = size(Ha,1)
   allocate(Fnn_temp_x(size(Fnnx,1),size(Fnnx,2)))
   allocate(Fnn_temp_y(size(Fnnx,1),size(Fnnx,2)))
   allocate(Fnn_temp_z(size(Fnnx,1),size(Fnnx,2)))
   allocate(PBx(size(Fnnx,1),size(Fnnx,2)))
   allocate(PBy(size(Fnnx,1),size(Fnnx,2)))
   allocate(PBz(size(Fnnx,1),size(Fnnx,2)))
   allocate(Norm1(size(Fnnx,1)))

   ! Do it with mkl:
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Scell(NSC)%PRRx, m, Ha, m, 0.0d0, PBx, m) ! Original
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Scell(NSC)%PRRy, m, Ha, m, 0.0d0, PBy, m)
!    CALL dgemm ('T','N', m, m, m, 1.0d0, Scell(NSC)%PRRz, m, Ha, m, 0.0d0, PBz, m)
   
   CALL dgemm ('N','N', m, m, m, 1.0d0, Scell(NSC)%PRRx, m, Ha, m, 0.0d0, PBx, m)
   CALL dgemm ('N','N', m, m, m, 1.0d0, Scell(NSC)%PRRy, m, Ha, m, 0.0d0, PBy, m)
   CALL dgemm ('N','N', m, m, m, 1.0d0, Scell(NSC)%PRRz, m, Ha, m, 0.0d0, PBz, m)
   
   
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
   
   ! Do whole arrays at once:
    !$OMP WORKSHARE
   Fnnx = Fnn_temp_x*Fnn_temp_x
   Fnny = Fnn_temp_y*Fnn_temp_y
   Fnnz = Fnn_temp_z*Fnn_temp_z
   !$OMP END WORKSHARE
   
   if (present(Fnnxy)) then ! do off-diagonal elements
      !$OMP WORKSHARE
      Fnnxy = Fnn_temp_x*Fnn_temp_y
      Fnnxz = Fnn_temp_x*Fnn_temp_z
      Fnnyx = Fnnxy
      Fnnyz = Fnn_temp_y*Fnn_temp_z
      Fnnzx = Fnnxz
      Fnnzy = Fnnyz
      !$OMP END WORKSHARE
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
   complex(8), intent(in) :: Eps_xx, Eps_yy, Eps_zz
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
subroutine get_RTA_from_CDF(ReEps, ImEps, ll, dd, teta0, drude_ray, R, T, A)
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
   n = sqrt((ReEps + sqrt(ReEps*ReEps+ImEps*ImEps))*0.5d0)
   k =sqrt((-ReEps + sqrt(ReEps*ReEps+ImEps*ImEps))*0.5d0)
end subroutine get_n_k

subroutine get_Re_and_Im(ReEps, ImEps, n, k)
   real(8), intent(out) :: ReEps, ImEps	! real and imagenary parts of dialectric function 
   real(8), intent(in) :: n, k	! n and k optical coefficients
   ReEps = n*n - k*k
   ImEps = 2.0d0*n*k
end subroutine get_Re_and_Im


!c creflection coefficient between two layers as function of two wave vectors:
function crefl(cxki,cxkj)
   complex :: crefl, cxki, cxkj
   crefl = (cxki-cxkj)/(cxki+cxkj)
end function crefl

!c transmission coefficient between two layers as function of two wave vectors:
function ctrans(cxki,cxkj)
   complex :: ctrans, cxki, cxkj
   ctrans = 2.0d0*cxki/(cxki+cxkj)
end function ctrans

!c phase factor within the layer i:
function cphase(cxkx,cd)
   complex :: cphase, cxkx, cd
   complex cc, carg
   cc = (0.0d0,1.0d0)
   carg = -cc*cxkx*cd
   cphase=exp(carg)
end function cphase





!======================================================
! Obsolete subs:


subroutine get_trani_all_complex_M(numpar, Scell, NSC, atoms, TB_Hamil, fe, all_w, Err)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err	! error save

   real(8), dimension(:,:), allocatable :: Fnnx, Fnny, Fnnz
   !real(8), dimension(:,:), allocatable :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), dimension(:,:), allocatable :: Eps
   real(8) Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   real(8) w, wmax, kx, ky, kz
   integer i, j, N, FN, ix, iy, iz, ixm, iym, izm, iw
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
   call allocate_Eps_hw(Scell(NSC)%eps%E_min, Scell(NSC)%eps%E_max, Scell(NSC)%eps%dE, Scell(NSC)%eps%Eps_hw) ! see below
   N = size(Scell(NSC)%eps%Eps_hw,2) ! use this size to define temporary arrays within this subroutine
   
   if (.not.allocated(Eps_hw)) allocate(Eps_hw(16,N)) ! all are there
   Eps_hw = 0.0d0
   if (.not.allocated(Eps_hw_temp)) then
      allocate(Eps_hw_temp(16,N)) ! all are there
      Eps_hw_temp = 0.0d0
   endif
!    wmax =50.0d0*g_e/g_h ! [1/s] maximum frequency to trace
   ixm = numpar%ixm
   iym = numpar%iym
   izm = numpar%izm
   
!    !$omp PARALLEL private(ix,iy,iz,i,Eps_hw_temp,CHij,Ei,Fnnx,Fnny,Fnnz,w,Re_eps,Im_eps,R,T,A,opt_n,opt_k) shared(Eps_hw)
!    !$omp do reduction( + : Eps_hw)
   do ix = 1, ixm
      do iy = 1, iym
         do iz = 1, izm
            
            print*, ix, iy, iz, ixm
            
            ! k-points choice from [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]
            kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
            ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
            kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))

            call get_Fnn_complex_M(numpar, Scell, NSC, atoms, TB_Hamil, CHij, Ei, Fnnx, Fnny, Fnnz, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
!             call get_Fnn_complex(numpar, matter, atoms, TB_Hamil, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx=real(ix)/real(ixm), ky=real(iy)/real(iym), kz=real(iz)/real(izm), Err=Err) ! see below

            if (all_w) then ! full spectrum:
! !$omp PARALLEL private(i,w,Eps_hw_temp,Re_eps,Im_eps,R,T,A,opt_n,opt_k)
! !$omp do reduction( + : Eps_hw)
               w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] frequency starting point for the optical spectrum [eV]
               do i = 1, N
!                   w = real(i)/real(N)*wmax
                  call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
!                 call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps)
                  call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz) ! Eps_hw_temp(1,i) = hw  ! energy
                  w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
               enddo
! !$omp end do
               Eps_hw = Eps_hw + Eps_hw_temp
! !$omp end parallel
            else ! only for given probe:
               w = Scell(NSC)%eps%w
               call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, dc_cond=dc_cond)
!                call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps)
            endif
         enddo ! iz
      enddo ! iy
   enddo ! ix

!    !$omp end do
!    !$omp end parallel
!    if (all_w) then ! full spectrum:
!       !Eps_hw = Eps_hw/real((ixm+1)*(iym+1)*(izm+1))
      Eps_hw = Eps_hw/dble(ixm*iym*izm)
      Scell(NSC)%eps%Eps_hw = Eps_hw

     Scell(NSC)%eps%ReEps = Eps_hw(2,1)	! real part of CDF
     Scell(NSC)%eps%ImEps = Eps_hw(3,1)	! imaginary part of CDF
     Scell(NSC)%eps%R = Eps_hw(5,1)	! reflectivity
     Scell(NSC)%eps%T = Eps_hw(6,1)	! transmission
     Scell(NSC)%eps%A = Eps_hw(7,1)	! absorption
     Scell(NSC)%eps%n = Eps_hw(8,1)	! optical n
     Scell(NSC)%eps%k = Eps_hw(9,1)	! optical k
     Scell(NSC)%eps%dc_cond = Eps_hw(10,1)	! dc-conductivity
     Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,1), Eps_hw(12,1))  ! Re_E_xx and Im_E_xx
     Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,1), Eps_hw(14,1))  ! Re_E_yy and Im_E_yy
     Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,1), Eps_hw(16,1))  ! Re_E_zz and Im_E_zz
!       ! Save parameters for the selected probe-pulse:
!       w = Scell(NSC)%eps%w
!       call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, fe, w, Re_eps, Im_eps, dc_cond=dc_cond)
!    endif
   
end subroutine get_trani_all_complex_M



subroutine get_trani_all_complex_F(numpar, Scell, NSC, atoms, TB_Hamil, fe, all_w, Err)
   type (Numerics_param), intent(inout) :: numpar ! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Fu), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err	! error save

   real(8), dimension(:,:), allocatable :: Fnnx, Fnny, Fnnz
   !real(8), dimension(:,:), allocatable :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), dimension(:,:), allocatable :: Eps
   real(8) Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   real(8) w, kx, ky, kz
   integer i, j, N, FN, ix, iy, iz, ixm, iym, izm
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
         allocate(Scell(NSC)%eps%Eps_hw(16,N)) ! 10 different parameters saved for all energy grid points
      endif     
   endif

   if (.not.allocated(Eps_hw)) allocate(Eps_hw(16,N)) ! all are there
   Eps_hw = 0.0d0
   if (.not.allocated(Eps_hw_temp)) then
      allocate(Eps_hw_temp(16,N)) ! all are there
      Eps_hw_temp = 0.0d0
   endif
!    wmax =50.0d0*g_e/g_h ! [1/s] maximum frequency to trace
   ixm = numpar%ixm
   iym = numpar%iym
   izm = numpar%izm

!    !$omp PARALLEL private(ix,iy,iz,i,Eps_hw_temp,CHij,Ei,Fnnx,Fnny,Fnnz,w,Re_eps,Im_eps,R,T,A,opt_n,opt_k) shared(Eps_hw)
!    !$omp do reduction( + : Eps_hw)
   do ix = 1, ixm
      do iy = 1, iym
         do iz = 1, izm
         
            print*, ix, iy, iz, ixm, 'P'
         
            ! k-points choice from [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]
            kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
            ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
            kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
            call get_Fnn_complex_F(numpar, Scell, NSC, atoms, TB_Hamil, CHij, Ei, Fnnx, Fnny, Fnnz, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
!             call get_Fnn_complex(numpar, matter, atoms, TB_Hamil, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx=real(ix)/real(ixm), ky=real(iy)/real(iym), kz=real(iz)/real(izm), Err=Err) ! see below

            if (all_w) then ! full spectrum:
               w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] frequency starting point for the optical spectrum [eV]
               do i = 1, N
                  !w = real(i)/real(N)*wmax
                  call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
                  !call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
                  call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy
                  w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
               enddo
               Eps_hw = Eps_hw + Eps_hw_temp
            else ! only for given probe:
               w = Scell(NSC)%eps%w
               !call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, fe, w, Re_eps, Im_eps, dc_cond=dc_cond)
!                call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps)
               call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
               !call get_trani(numpar, matter, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, Ei, fe, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k)
               call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy   
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
!       !Eps_hw = Eps_hw/real((ixm+1)*(iym+1)*(izm+1))
      Eps_hw = Eps_hw/dble(ixm*iym*izm)
      Scell(NSC)%eps%Eps_hw = Eps_hw     
     
     Scell(NSC)%eps%ReEps = Eps_hw(2,1)	! real part of CDF
     Scell(NSC)%eps%ImEps = Eps_hw(3,1)	! imaginary part of CDF
     Scell(NSC)%eps%R = Eps_hw(5,1)	! reflectivity
     Scell(NSC)%eps%T = Eps_hw(6,1)	! transmission
     Scell(NSC)%eps%A = Eps_hw(7,1)	! absorption
     Scell(NSC)%eps%n = Eps_hw(8,1)	! optical n
     Scell(NSC)%eps%k = Eps_hw(9,1)	! optical k
     Scell(NSC)%eps%dc_cond = Eps_hw(10,1)	! dc-conductivity
     Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,1), Eps_hw(12,1))  ! Re_E_xx and Im_E_xx
     Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,1), Eps_hw(14,1))  ! Re_E_yy and Im_E_yy
     Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,1), Eps_hw(16,1))  ! Re_E_zz and Im_E_zz
     
!       ! Save parameters for the selected probe-pulse:
!       w = Scell(NSC)%eps%w
!       call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, fe, w, Re_eps, Im_eps, dc_cond=dc_cond)
!    endif
end subroutine get_trani_all_complex_F


subroutine get_trani_all_complex_NRL(numpar, Scell, NSC, TB_Hamil, fe, all_w, Err)
 type (Numerics_param), intent(inout) :: numpar	! numerical parameters, including drude-function 
   type(Super_cell), dimension(:), intent(inout) :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   type(TB_H_NRL), dimension(:,:), intent(in) :: TB_Hamil	! parameters of the Hamiltonian of TB
   real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err	! error save
   !----------------------------------------------------------------------------
   real(8), dimension(:,:), allocatable :: Fnnx, Fnny, Fnnz
   real(8), dimension(:,:), allocatable :: Eps
   real(8) Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   real(8) w, kx, ky, kz
   integer i, j, N, FN, ix, iy, iz, ixm, iym, izm, schem, Ngp, Nsiz
   real(8), dimension(:), allocatable :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable :: Eps_hw	! array of all eps vs hw
   real(8), dimension(:,:), allocatable :: Eps_hw_temp	! array of all eps vs hw

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

   do ix = 1, ixm
      do iy = 1, iym
         do iz = 1, izm
         
            print*, ix, iy, iz, ixm, 'NRL'
            Ngp = (ix-1)*iym*ixm + (iy-1)*ixm + iz	! number of grid point for user defined grid
            if (Ngp > Nsiz) goto 3456
            
!             ! k-points choice from [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]
!             kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
!             ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
!             kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
            call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid)	! module "TB"
            
            call get_Fnn_complex_NRL(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
            
            if (all_w) then ! full spectrum:
               w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] frequency starting point for the optical spectrum [eV]
               do i = 1, N
                  call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
                  call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy
                  w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
               enddo
               Eps_hw = Eps_hw + Eps_hw_temp
            else ! only for given probe:
               w = Scell(NSC)%eps%w
               call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
               call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy   
               Eps_hw = Eps_hw + Eps_hw_temp
            endif
         enddo ! iz
      enddo ! iy
   enddo ! ix

3456  continue
   !Eps_hw = Eps_hw/dble(ixm*iym*izm)
   Eps_hw = Eps_hw/dble(Nsiz)
   Scell(NSC)%eps%Eps_hw = Eps_hw     
     
   Scell(NSC)%eps%ReEps = Eps_hw(2,1)	! real part of CDF
   Scell(NSC)%eps%ImEps = Eps_hw(3,1)	! imaginary part of CDF
   Scell(NSC)%eps%R = Eps_hw(5,1)	! reflectivity
   Scell(NSC)%eps%T = Eps_hw(6,1)	! transmission
   Scell(NSC)%eps%A = Eps_hw(7,1)	! absorption
   Scell(NSC)%eps%n = Eps_hw(8,1)	! optical n
   Scell(NSC)%eps%k = Eps_hw(9,1)	! optical k
   Scell(NSC)%eps%dc_cond = Eps_hw(10,1)	! dc-conductivity
   Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,1), Eps_hw(12,1))  ! Re_E_xx and Im_E_xx
   Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,1), Eps_hw(14,1))  ! Re_E_yy and Im_E_yy
   Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,1), Eps_hw(16,1))  ! Re_E_zz and Im_E_zz
end subroutine get_trani_all_complex_NRL


subroutine get_trani_all_complex_DFTB(numpar, Scell, NSC, fe, all_w, Err)
   type (Numerics_param), intent(inout) :: numpar	! numerical parameters, including drude-function
   type(Super_cell), dimension(:), intent(inout) :: Scell	! supercell with all the atoms as one object
   integer, intent(in) :: NSC	! number of supercell
   real(8), dimension(:), intent(in) :: fe	! electron distribution function
   logical, intent(in) :: all_w	! get all spectrum of hv, or only for given probe wavelength
   type(Error_handling), intent(inout) :: Err	! error save
   !----------------------------------------------------------------------------
   real(8), dimension(:,:), allocatable :: Fnnx, Fnny, Fnnz
   real(8), dimension(:,:), allocatable :: Fnnx_test, Fnny_test, Fnnz_test
   real(8), dimension(:,:), allocatable :: Eps
   real(8) Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond
   real(8) w, kx, ky, kz
   integer i, j, N, FN, ix, iy, iz, ixm, iym, izm, schem, Ngp, Nsiz
   real(8), dimension(:), allocatable :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable :: Eps_hw	! array of all eps vs hw
   real(8), dimension(:,:), allocatable :: Eps_hw_temp	! array of all eps vs hw

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

   do ix = 1, ixm
      do iy = 1, iym
         do iz = 1, izm
         
            print*, ix, iy, iz, ixm, 'DFTB'
            Ngp = (ix-1)*iym*ixm + (iy-1)*ixm + iz	! number of grid point for user defined grid
            if (Ngp > Nsiz) goto 3456
            
!             ! k-points choice from [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]
!             kx = (2.0d0*real(ix) - real(ixm) - 1.0d0)/(2.0d0*real(ixm))
!             ky = (2.0d0*real(iy) - real(iym) - 1.0d0)/(2.0d0*real(iym))
!             kz = (2.0d0*real(iz) - real(izm) - 1.0d0)/(2.0d0*real(izm))
            call k_point_choice(schem, ix, iy, iz, ixm, iym, izm, kx, ky, kz, numpar%k_grid)	! module "TB"
            
            call get_Fnn_complex_DFTB(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
            
!             call get_Fnn_complex(numpar, Scell, NSC, Ei, Fnnx_test, Fnny_test, Fnnz_test, kx=kx, ky=ky, kz=kz, Err=Err) ! see below
!             do i = 1, size(Fnnx,1)
!                do j = 1, size(Fnnx,2)
!                   if (ABS(Fnnx(i,j)) > 1.0d-13) print*, i, j, Fnnx(i,j), Fnnx_test(i,j) ! the same, CHECKED
!                enddo
!             enddo
!             pause 'Test get_Fnn_complex'
            
            if (all_w) then ! full spectrum:
               w = Scell(NSC)%eps%E_min*g_e/g_h	! [1/s] frequency starting point for the optical spectrum [eV]
               do i = 1, N
                  call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
                  call save_Eps_hw(Eps_hw_temp, i, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy
                  w = w + Scell(NSC)%eps%dE*g_e/g_h	! [1/s] frequency next grid point [eV]
               enddo
               Eps_hw = Eps_hw + Eps_hw_temp
            else ! only for given probe:
               w = Scell(NSC)%eps%w
               call get_trani(numpar, Scell, NSC, Fnnx, Fnny, Fnnz, Ei, w, Re_eps, Im_eps, R, T, A, opt_n, opt_k, dc_cond=dc_cond)
               call save_Eps_hw(Eps_hw_temp, 1, w/g_e*g_h, Re_eps, Im_eps, Im_eps/(Re_eps*Re_eps + Im_eps*Im_eps), R, T, A, opt_n, opt_k, Scell(NSC)%eps%Eps_xx, Scell(NSC)%eps%Eps_yy, Scell(NSC)%eps%Eps_zz, dc_cond=dc_cond) ! Eps_hw_temp(1,i) = hw  ! energy   
               Eps_hw = Eps_hw + Eps_hw_temp
            endif
         enddo ! iz
      enddo ! iy
   enddo ! ix

3456  continue
   !Eps_hw = Eps_hw/dble(ixm*iym*izm)
   Eps_hw = Eps_hw/dble(Nsiz)
   Scell(NSC)%eps%Eps_hw = Eps_hw     
     
   Scell(NSC)%eps%ReEps = Eps_hw(2,1)	! real part of CDF
   Scell(NSC)%eps%ImEps = Eps_hw(3,1)	! imaginary part of CDF
   Scell(NSC)%eps%R = Eps_hw(5,1)	! reflectivity
   Scell(NSC)%eps%T = Eps_hw(6,1)	! transmission
   Scell(NSC)%eps%A = Eps_hw(7,1)	! absorption
   Scell(NSC)%eps%n = Eps_hw(8,1)	! optical n
   Scell(NSC)%eps%k = Eps_hw(9,1)	! optical k
   Scell(NSC)%eps%dc_cond = Eps_hw(10,1)	! dc-conductivity
   Scell(NSC)%eps%Eps_xx = dcmplx(Eps_hw(11,1), Eps_hw(12,1))  ! Re_E_xx and Im_E_xx
   Scell(NSC)%eps%Eps_yy = dcmplx(Eps_hw(13,1), Eps_hw(14,1))  ! Re_E_yy and Im_E_yy
   Scell(NSC)%eps%Eps_zz = dcmplx(Eps_hw(15,1), Eps_hw(16,1))  ! Re_E_zz and Im_E_zz
end subroutine get_trani_all_complex_DFTB





subroutine get_Fnn_complex_NRL(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx, ky, kz, Err) ! see below
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), allocatable, intent(out) :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable, intent(out) :: Fnnx, Fnny, Fnnz
   real(8), intent(in) :: kx ,ky, kz ! k-point to get Hamiltonian at
   type(Error_handling), intent(inout) :: Err	! error save
   !----------------------------------------------------------------------------
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
!    complex, dimension(:,:), allocatable :: CHij_temp	! eigenvectors of the hamiltonian
   complex, dimension(:,:), allocatable :: CSij	! overlap matrix of the nonorthogonal hamiltonian
   complex*8, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   complex*8, dimension(:,:), allocatable :: PBx, PBy, PBz
   integer i, j, n, nn, Ne, m
   character(200) :: Error_descript
   
   real(8), dimension(:,:), allocatable :: Hij
   real(8), dimension(:), allocatable :: Ei_r
   
   Error_descript = ''

   ! Allocate Hamiltonan matrices:
   Ne = 9.0d0*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
   m = Ne
   
   if (.not.allocated(Ei)) allocate(Ei(Ne))
   if (.not.allocated(Fnnx)) allocate(Fnnx(Ne,Ne))
   if (.not.allocated(Fnny)) allocate(Fnny(Ne,Ne))
   if (.not.allocated(Fnnz)) allocate(Fnnz(Ne,Ne))
   allocate(CHij(Ne,Ne))
   allocate(CSij(Ne,Ne))
   allocate(Fnn_temp_x(Ne,Ne))
   allocate(Fnn_temp_y(Ne,Ne))
   allocate(Fnn_temp_z(Ne,Ne))
   allocate(PBx(Ne,Ne))
   allocate(PBy(Ne,Ne))
   allocate(PBz(Ne,Ne))

   ! Construct complex Hamiltonian from the real one for the given k-point:
   call Complex_Hamil_NRL(numpar, Scell, NSC, CHij, CSij, Ei, kx, ky, kz, Err) ! "TB_NRL"
   
!    ! TEST:
!    do j = 1,size(CHij,1)
!       do i = 1,size(CHij,2)
!          if ( ABS(dble(CSij(j,i)) - Scell(NSC)%Sij(j,i)) > 1.0d-14 ) print*, j, i, CSij(j,i), Scell(NSC)%Sij(j,i)
!       enddo
!    enddo
!    PAUSE 'dble(CSij(j,i)) - '

!    do j = 1,size(CHij,1)
!       do i = 1,size(CHij,2)
!          if (isnan(REAL(CHij(j,i))) .or. isnan(AIMAG(CHij(j,i)))) then
!              Error_descript = 'Module Optical_parameters: subroutine get_Fnn_complex_NRL got NaNs'
!              call Save_error_details(Err, 8, Error_descript)
!              print*, trim(adjustl(Error_descript)), j, i, CHij(j,i)
!          endif
!       enddo
!    enddo
   
   !$omp PARALLEL private(i,nn,n,j)
   !$omp do
   do i = 1, m
      do nn = 1, m
         PBx(i,nn) = SUM(Scell(NSC)%cPRRx(i,:)*CHij(:,nn))
         PBy(i,nn) = SUM(Scell(NSC)%cPRRy(i,:)*CHij(:,nn))
         PBz(i,nn) = SUM(Scell(NSC)%cPRRz(i,:)*CHij(:,nn))
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp do
   do n = 1, m
      do nn = 1, m
         Fnn_temp_x(n,nn) = SUM(CONJG(CHij(:,n))*PBx(:,nn))
         Fnn_temp_y(n,nn) = SUM(CONJG(CHij(:,n))*PBy(:,nn))
         Fnn_temp_z(n,nn) = SUM(CONJG(CHij(:,n))*PBz(:,nn))
      enddo ! nn
   enddo ! n
   !$omp end do
   !$omp do
   do i = 1, m
      do j = 1, m
!        ! Must convert it to double complex, otherwise it gives wrong results!
         Fnnx(i,j) = DCMPLX(CONJG(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j))
         Fnny(i,j) = DCMPLX(CONJG(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         Fnnz(i,j) = DCMPLX(CONJG(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_z(i,j))
      enddo
   enddo
   !$omp end do
   !$omp end parallel
   
   deallocate (Fnn_temp_x,Fnn_temp_y,Fnn_temp_z,PBx,PBy,PBz,CHij,CSij)
end subroutine get_Fnn_complex_NRL



subroutine get_Fnn_complex_DFTB(numpar, Scell, NSC, Ei, Fnnx, Fnny, Fnnz, kx, ky, kz, Err) ! see below
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:), allocatable, intent(out) :: Ei	! energy levels [eV]
   real(8), dimension(:,:), allocatable, intent(out) :: Fnnx, Fnny, Fnnz
   real(8), intent(in) :: kx ,ky, kz ! k-point to get Hamiltonian at
   type(Error_handling), intent(inout) :: Err	! error save
   !----------------------------------------------------------------------------
   complex, dimension(:,:), allocatable :: CHij	! eigenvectors of the hamiltonian
!    complex, dimension(:,:), allocatable :: CHij_temp	! eigenvectors of the hamiltonian
   complex, dimension(:,:), allocatable :: CSij	! overlap matrix of the nonorthogonal hamiltonian
   complex*8, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   complex*8, dimension(:,:), allocatable :: PBx, PBy, PBz
   integer i, j, n, nn, Ne, m, n_orb
   character(200) :: Error_descript
   
   real(8), dimension(:,:), allocatable :: Hij
   real(8), dimension(:), allocatable :: Ei_r
   
   complex, dimension(:,:), allocatable :: CHij_test
   real(8), dimension(:), allocatable :: Ei_test	   ! energy levels [eV]
   
   Error_descript = ''

   ! Size of the basis set per atom:
   n_orb =  identify_DFTB_orbitals_per_atom(numpar%N_basis_size)    ! module "TB_DFTB"
   Ne = n_orb*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
   m = Ne
   ! Allocate Hamiltonan matrices:   
   if (.not.allocated(Ei)) allocate(Ei(Ne))
   if (.not.allocated(Ei_test)) allocate(Ei_test(Ne))
   if (.not.allocated(Fnnx)) allocate(Fnnx(Ne,Ne))
   if (.not.allocated(Fnny)) allocate(Fnny(Ne,Ne))
   if (.not.allocated(Fnnz)) allocate(Fnnz(Ne,Ne))
   allocate(CHij(Ne,Ne))
   allocate(CHij_test(Ne,Ne))
   allocate(CSij(Ne,Ne))
   allocate(Fnn_temp_x(Ne,Ne))
   allocate(Fnn_temp_y(Ne,Ne))
   allocate(Fnn_temp_z(Ne,Ne))
   allocate(PBx(Ne,Ne))
   allocate(PBy(Ne,Ne))
   allocate(PBz(Ne,Ne))

   ! Construct complex Hamiltonian from the real one for the given k-point:
   call Complex_Hamil_DFTB(numpar, Scell, NSC, CHij, CSij, Ei, kx, ky, kz, Err) ! "TB_DFTB"
   
!    call construct_complex_Hamiltonian(numpar, Scell, NSC, Scell(NSC)%H_non, CHij_test, Ei_test, kx, ky, kz, Err, Scell(NSC)%cPRRx, Scell(NSC)%cPRRy, Scell(NSC)%cPRRz, Scell(NSC)%Sij, CSij)
   
!    print*, ' kx, ky, kz ',  kx, ky, kz
   
!    do i = 1,size(CHij,1)
!       print*, i, Ei(i), Ei_test(i)
!    enddo
!    pause 'Test construct_complex_Hamiltonian'
   
!     ! test:
!     do j = 250,size(CHij,1)
!        do i = 250,size(CHij,2)
! !           if ( ABS(dble(CSij(j,i)) - Scell(NSC)%Sij(j,i)) > 1.0d-14 ) print*, j, i, CSij(j,i), Scell(NSC)%Sij(j,i)
!          if ( ABS(dble(CHij(j,i)))  > 1.0d-14 ) print*, j, i, CHij(j,i) , CHij_test(j,i)  ! The same, CHECKED
!        enddo
!     enddo
!     PAUSE ' CHij '

!    do j = 1,size(CHij,1)
!       do i = 1,size(CHij,2)
!          if (isnan(REAL(CHij(j,i))) .or. isnan(AIMAG(CHij(j,i)))) then
!              Error_descript = 'Module Optical_parameters: subroutine get_Fnn_complex_NRL got NaNs'
!              call Save_error_details(Err, 8, Error_descript)
!              print*, trim(adjustl(Error_descript)), j, i, CHij(j,i)
!          endif
!       enddo
!    enddo
   
   !$omp PARALLEL private(i,nn,n,j)
   !$omp do
   do i = 1, m
      do nn = 1, m
         PBx(i,nn) = SUM(Scell(NSC)%cPRRx(i,:)*CHij(:,nn))
         PBy(i,nn) = SUM(Scell(NSC)%cPRRy(i,:)*CHij(:,nn))
         PBz(i,nn) = SUM(Scell(NSC)%cPRRz(i,:)*CHij(:,nn))
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp do
   do n = 1, m
      do nn = 1, m
         Fnn_temp_x(n,nn) = SUM(CONJG(CHij(:,n))*PBx(:,nn))
         Fnn_temp_y(n,nn) = SUM(CONJG(CHij(:,n))*PBy(:,nn))
         Fnn_temp_z(n,nn) = SUM(CONJG(CHij(:,n))*PBz(:,nn))
      enddo ! nn
   enddo ! n
   !$omp end do
   !$omp do
   do i = 1, m
      do j = 1, m
!        ! Must convert it to double complex, otherwise it gives wrong results!
         Fnnx(i,j) = DCMPLX(CONJG(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j))
         Fnny(i,j) = DCMPLX(CONJG(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         Fnnz(i,j) = DCMPLX(CONJG(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_z(i,j))
      enddo
   enddo
   !$omp end do
   !$omp end parallel
   
   deallocate (Fnn_temp_x,Fnn_temp_y,Fnn_temp_z,PBx,PBy,PBz,CHij,CSij)
end subroutine get_Fnn_complex_DFTB



subroutine get_Fnn_complex_F(numpar, Scell, NSC, atoms, TB, CHij, Ei, Fnnx, Fnny, Fnnz, kx ,ky, kz, Err)
!subroutine get_Fnn_complex(numpar, matter, atoms, TB, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx ,ky, kz, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Fu), dimension(:,:), intent(in) :: TB   ! parameters of the Hamiltonian of TB
   real(8), dimension(:), allocatable, intent(out) :: Ei	! energy levels [eV]
   complex, dimension(:,:), allocatable, intent(out) :: CHij	! eigenvectors of the hamiltonian
   real(8), dimension(:,:), allocatable, intent(out) :: Fnnx, Fnny, Fnnz
   !real(8), dimension(:,:), allocatable, intent(out) :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), intent(in) :: kx ,ky, kz ! k-point to get Hamiltonian at
   type(Error_handling), intent(inout) :: Err	! error save
   ! complex, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   ! complex, dimension(:,:), allocatable :: PBx, PBy, PBz
   complex*8, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   complex*8, dimension(:,:), allocatable :: PBx, PBy, PBz
   integer i, j, m, n, nn, FN, Ne
   character(200) :: Error_descript
   
   
   real(8), dimension(:,:), allocatable :: Hij
   real(8), dimension(:), allocatable :: Ei_r
   
   
   Error_descript = ''
   
   Ne = Scell(NSC)%Ne  ! initial number of VB electrons
   m = Ne

   ! Allocate Hamiltonan matrices:
   ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
      select type(ARRAY)
      type is (TB_H_Pettifor)
         m = size(ARRAY(NSC,NSC)%V0)*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
      type is (TB_H_Molteni)
         m = size(ARRAY(NSC,NSC)%V0)*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
      type is (TB_H_Fu)
         m = size(ARRAY(NSC,NSC)%V0)*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
      end select
   END ASSOCIATE
   Ne = m


   if (.not.allocated(CHij)) allocate(CHij(Ne,Ne))
   if (.not.allocated(Ei)) allocate(Ei(Ne))
   if (.not.allocated(Fnnx)) allocate(Fnnx(Ne,Ne))
   if (.not.allocated(Fnny)) allocate(Fnny(Ne,Ne))
   if (.not.allocated(Fnnz)) allocate(Fnnz(Ne,Ne))
   
   if (.not.allocated(Hij)) allocate(Hij(Ne,Ne))
   if (.not.allocated(Ei_r)) allocate(Ei_r(Ne))
   
!    if (.not.allocated(Fnnxx)) allocate(Fnnxx(Ne,Ne))
!    if (.not.allocated(Fnnxy)) allocate(Fnnxy(Ne,Ne))
!    if (.not.allocated(Fnnxz)) allocate(Fnnxz(Ne,Ne))
!    if (.not.allocated(Fnnyx)) allocate(Fnnyx(Ne,Ne))
!    if (.not.allocated(Fnnyy)) allocate(Fnnyy(Ne,Ne))
!    if (.not.allocated(Fnnyz)) allocate(Fnnyz(Ne,Ne))
!    if (.not.allocated(Fnnzx)) allocate(Fnnzx(Ne,Ne))
!    if (.not.allocated(Fnnzy)) allocate(Fnnzy(Ne,Ne))
!    if (.not.allocated(Fnnzz)) allocate(Fnnzz(Ne,Ne))
   allocate(Fnn_temp_x(Ne,Ne))
   allocate(Fnn_temp_y(Ne,Ne))
   allocate(Fnn_temp_z(Ne,Ne))
   allocate(PBx(Ne,Ne))
   allocate(PBy(Ne,Ne))
   allocate(PBz(Ne,Ne))

   call Complex_Hamil_tot_F(numpar, Scell, NSC, atoms, TB, CHij=CHij, ksx=kx, ksy=ky, ksz=kz) ! "TB_Fu"
   
   call sym_diagonalize(CHij, Ei, Err%Err_descript) ! module "Algebra_tools"
!    pause 'Complex_Hamil_tot TEST'

   do j = 1,size(CHij,1)
      do i = 1,size(CHij,2)
         if (isnan(REAL(CHij(j,i))) .or. isnan(AIMAG(CHij(j,i)))) then
             Error_descript = 'Module Optical_parameters: subroutine get_Fnn_complex got NaNs'
             call Save_error_details(Err, 8, Error_descript)
             print*, trim(adjustl(Error_descript)), j, i, CHij(j,i)
         endif
      enddo
   enddo
   
!    call Hamil_tot(numpar, Scell, NSC, TB, Hij)
!    call sym_diagonalize(Hij, Ei_r, Err%Err_descript) ! modeule "Algebra_tools"
!    do i = 1, size(Ei) ! compare Real and Complex Hamiltonians eigenvalues:
!       write(*,'(i,e,e,f,f)') i, Ei(i), Ei_r(i), CHij(i,1)
!    enddo
!    PAUSE 'Test CHij'

   !$omp PARALLEL private(i,nn,n,j)
   !$omp do
   do i = 1, m
      do nn = 1, m
         PBx(i,nn) = SUM(Scell(NSC)%cPRRx(i,:)*CHij(:,nn))
         PBy(i,nn) = SUM(Scell(NSC)%cPRRy(i,:)*CHij(:,nn))
         PBz(i,nn) = SUM(Scell(NSC)%cPRRz(i,:)*CHij(:,nn))
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp do
   do n = 1, m
      do nn = 1, m
         Fnn_temp_x(n,nn) = SUM(CONJG(CHij(:,n))*PBx(:,nn))
         Fnn_temp_y(n,nn) = SUM(CONJG(CHij(:,n))*PBy(:,nn))
         Fnn_temp_z(n,nn) = SUM(CONJG(CHij(:,n))*PBz(:,nn))
      enddo ! nn
   enddo ! n
   !$omp end do
   !$omp do
   do i = 1, m
      do j = 1, m
!          Fnnx(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_x(i,j)
!          Fnny(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_y(i,j)
!          Fnnz(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_z(i,j)

!        ! Must convert it to double complex, otherwise it gives wrong results!
         Fnnx(i,j) = DCMPLX(CONJG(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j))
         Fnny(i,j) = DCMPLX(CONJG(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         Fnnz(i,j) = DCMPLX(CONJG(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_z(i,j))
         
!          Fnnxx(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_x(i,j)
!          Fnnxy(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_y(i,j)
!          Fnnxz(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_z(i,j)
!          Fnnyx(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_x(i,j)
!          Fnnyy(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_y(i,j)
!          Fnnyz(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_z(i,j)
!          Fnnzx(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_x(i,j)
!          Fnnzy(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_y(i,j)
!          Fnnzz(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_z(i,j)
      enddo
   enddo
   !$omp end do
   !$omp end parallel


!   open(NEWUNIT=FN, FILE = 'OUTPUT_PBx.dat')
!      do i = 1, size(PBx,1)
!         do j = 1, size(PBx,2)
!            write(FN, '(i,i,e)') i, j, PBx(i,j)
!         enddo
!      enddo
!   close(FN)
!   pause 'get_Fnn'

end subroutine get_Fnn_complex_F


subroutine get_Fnn_complex_M(numpar, Scell, NSC, atoms, TB, CHij, Ei, Fnnx, Fnny, Fnnz, kx ,ky, kz, Err)
!subroutine get_Fnn_complex(numpar, matter, atoms, TB, CHij, Ei, Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz, kx ,ky, kz, Err)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB   ! parameters of the Hamiltonian of TB
   real(8), dimension(:), allocatable, intent(out) :: Ei	! energy levels [eV]
   complex, dimension(:,:), allocatable, intent(out) :: CHij	! eigenvectors of the hamiltonian
   real(8), dimension(:,:), allocatable, intent(out) :: Fnnx, Fnny, Fnnz
   !real(8), dimension(:,:), allocatable, intent(out) :: Fnnxx, Fnnxy, Fnnxz, Fnnyx, Fnnyy, Fnnyz, Fnnzx, Fnnzy, Fnnzz
   real(8), intent(in) :: kx ,ky, kz ! k-point to get Hamiltonian at
   type(Error_handling), intent(inout) :: Err	! error save
   complex*8, dimension(:,:), allocatable :: Fnn_temp_x, Fnn_temp_y, Fnn_temp_z
   complex*8, dimension(:,:), allocatable :: PBx, PBy, PBz
   integer i, j, m, n, nn, FN, Ne
   character(200) :: Error_descript
   Error_descript = ''

   ! Allocate Hamiltonan matrices:
   ASSOCIATE (ARRAY => Scell(NSC)%TB_Hamil(:,:))
      select type(ARRAY)
      type is (TB_H_Pettifor)
         m = size(ARRAY(NSC,NSC)%V0)*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
      type is (TB_H_Molteni)
         m = size(ARRAY(NSC,NSC)%V0)*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
       type is (TB_H_Fu)
         m = size(ARRAY(NSC,NSC)%V0)*Scell(NSC)%Na ! number of energy levels is defined by the number of TB parameters included
      end select
   END ASSOCIATE
   Ne = m

   if (.not.allocated(CHij)) allocate(CHij(Ne,Ne))
   if (.not.allocated(Ei)) allocate(Ei(Ne))
   if (.not.allocated(Fnnx)) allocate(Fnnx(Ne,Ne))
   if (.not.allocated(Fnny)) allocate(Fnny(Ne,Ne))
   if (.not.allocated(Fnnz)) allocate(Fnnz(Ne,Ne))
!    if (.not.allocated(Fnnxx)) allocate(Fnnxx(Ne,Ne))
!    if (.not.allocated(Fnnxy)) allocate(Fnnxy(Ne,Ne))
!    if (.not.allocated(Fnnxz)) allocate(Fnnxz(Ne,Ne))
!    if (.not.allocated(Fnnyx)) allocate(Fnnyx(Ne,Ne))
!    if (.not.allocated(Fnnyy)) allocate(Fnnyy(Ne,Ne))
!    if (.not.allocated(Fnnyz)) allocate(Fnnyz(Ne,Ne))
!    if (.not.allocated(Fnnzx)) allocate(Fnnzx(Ne,Ne))
!    if (.not.allocated(Fnnzy)) allocate(Fnnzy(Ne,Ne))
!    if (.not.allocated(Fnnzz)) allocate(Fnnzz(Ne,Ne))
   allocate(Fnn_temp_x(Ne,Ne))
   allocate(Fnn_temp_y(Ne,Ne))
   allocate(Fnn_temp_z(Ne,Ne))
   allocate(PBx(Ne,Ne))
   allocate(PBy(Ne,Ne))
   allocate(PBz(Ne,Ne))

   call Complex_Hamil_tot_Molteni(numpar, Scell, NSC, atoms, TB, CHij=CHij, ksx=kx, ksy=ky, ksz=kz) ! "TB_Hamiltonian"

   do j = 1,size(CHij,1)
      do i = 1,size(CHij,2)
         if (isnan(REAL(CHij(j,i))) .or. isnan(AIMAG(CHij(j,i)))) then
             !print*, j,i,CHij(j,i)
             Error_descript = 'Module Optical_parameters: subroutine get_Fnn_complex_M got NaNs'
             call Save_error_details(Err, 8, Error_descript)
             print*, trim(adjustl(Error_descript)), j, i, CHij(j,i)
         endif
      enddo
   enddo
   
   call sym_diagonalize(CHij, Ei, Err%Err_descript) ! module "Algebra_tools"

!    do i = 1, size(Ei) ! compare Real and Complex Hamiltonians eigenvalues:
!       write(*,'(i,f,f,f)') i, Ei(i), CHij(i,1)
!    enddo

   !$omp PARALLEL private(i,nn,n,j)
   !$omp do
   do i = 1, m
      do nn = 1, m
         PBx(i,nn) = SUM(Scell(NSC)%cPRRx(i,:)*CHij(:,nn))
         PBy(i,nn) = SUM(Scell(NSC)%cPRRy(i,:)*CHij(:,nn))
         PBz(i,nn) = SUM(Scell(NSC)%cPRRz(i,:)*CHij(:,nn))
!          write(*,'(es,es,es,es,es,es)') PBx(i,nn), PBy(i,nn), PBz(i,nn)
      enddo ! j
   enddo ! i
   !$omp end do
   !$omp do
   do n = 1, m
      do nn = 1, m
         Fnn_temp_x(n,nn) = SUM(CONJG(CHij(:,n))*PBx(:,nn))
         Fnn_temp_y(n,nn) = SUM(CONJG(CHij(:,n))*PBy(:,nn))
         Fnn_temp_z(n,nn) = SUM(CONJG(CHij(:,n))*PBz(:,nn))
!          write(*,'(es,es,es,es,es,es)') Fnn_temp_x(n,nn), Fnn_temp_y(n,nn), Fnn_temp_z(n,nn)
      enddo ! nn
   enddo ! n
   !$omp end do
   
   !$omp do
   do i = 1, m
      do j = 1, m
!          Fnnx(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_x(i,j)
!          Fnny(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_y(i,j)
!          Fnnz(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_z(i,j)

         ! Must convert it to double complex, otherwise it gives wrong results!
         Fnnx(i,j) = DCMPLX(CONJG(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j))
         Fnny(i,j) = DCMPLX(CONJG(Fnn_temp_y(i,j)))*DCMPLX(Fnn_temp_y(i,j))
         Fnnz(i,j) = DCMPLX(CONJG(Fnn_temp_z(i,j)))*DCMPLX(Fnn_temp_z(i,j))

!          Fnnxx(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_x(i,j)
!          Fnnxy(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_y(i,j)
!          Fnnxz(i,j) = CONJG(Fnn_temp_x(i,j))*Fnn_temp_z(i,j)
!          Fnnyx(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_x(i,j)
!          Fnnyy(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_y(i,j)
!          Fnnyz(i,j) = CONJG(Fnn_temp_y(i,j))*Fnn_temp_z(i,j)
!          Fnnzx(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_x(i,j)
!          Fnnzy(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_y(i,j)
!          Fnnzz(i,j) = CONJG(Fnn_temp_z(i,j))*Fnn_temp_z(i,j)
         !write(*,'(es,es,es,es,es,es)') Fnnx(i,j), Fnny(i,j), Fnnz(i,j)
!          write(*,'(es,es,es,es,es,es)') Fnn_temp_x(i,j), DCMPLX(CONJG(Fnn_temp_x(i,j)))*DCMPLX(Fnn_temp_x(i,j)), Fnnx(i,j)
      enddo
   enddo
   !$omp end do
   !$omp end parallel

!   open(NEWUNIT=FN, FILE = 'OUTPUT_PBx.dat')
!      do i = 1, size(PBx,1)
!         do j = 1, size(PBx,2)
!            write(FN, '(i,i,e)') i, j, PBx(i,j)
!         enddo
!      enddo
!   close(FN)
!   pause 'get_Fnn'

end subroutine get_Fnn_complex_M



END MODULE Optical_parameters
