! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2016-2021 Nikita Medvedev
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
! This module contains subroutines to calculate the electron-ion coupling:

MODULE Nonadiabatic
use Universal_constants
use Objects
!use Variables
use Algebra_tools
use Electron_tools
use Atomic_tools


implicit none

 contains


subroutine Electron_ion_collision_int(Scell, numpar, nrg, Mij, wr, wr0, distre, dE_nonadiabat, kind_M, DOS_weights, G_ei_partial)
   type(Super_cell), intent(inout) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(Energies), intent(inout) :: nrg	! all energies
   real(8), dimension(:,:), intent(in) :: Mij  ! Matrix element coupling electron energy levels via ion motion
   real(8), dimension(:), intent(in) :: wr  ! electron energy levels [eV]
   real(8), dimension(:), intent(in) :: wr0  ! electron energy levels on last timestep [eV]
   real(8), dimension(:), intent(inout) :: distre   ! electron distribution function
   real(8), intent(out) :: dE_nonadiabat ! energy exchanged [eV]
   integer, intent(in) :: kind_M	! index of the model used for the matri element calculation
   real(8), dimension(:,:,:), intent(in) :: DOS_weights   ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:,:), intent(inout) :: G_ei_partial ! partial contributions into coupling from shells of atoms
   !---------------------------------------
   integer :: i, N_orb, N_steps, i_t, i_pc
   real(8) :: coef, coef_inv, dt_small, N_tot_cut, E_tot_cur, alpha, A_pc, B_pc, AB, A2AB, arg, eps, eps_N
   real(8), dimension(:), allocatable :: dfdt_e ! change of the distribution function
   real(8), dimension(:), allocatable :: dfdt_e2 ! predictotr change of the distribution function
   real(8), dimension(:), allocatable :: distre_temp  ! electron distribution
   real(8), dimension(:), allocatable :: distre_temp2  ! predictor electron distribution
   real(8), dimension(:), allocatable :: dfdt_A, dfdt_B  ! factors for explicit scheme
   real(8), dimension(:), allocatable :: dfdt_A2, dfdt_B2  ! predictor factors for explicit scheme
!    real(8), dimension(:), allocatable :: inv_t_e_ph ! inverse electron-phonon scattering time [1/fs]
   logical :: pc_iterations

   alpha = 0.35d0  ! mixing parameter in predictor-corrector
   eps = 1.0d-10   ! precision in A-B
   eps_N = 1.0d-6  ! precision in number of electrons
   
   
   N_orb = size(wr) ! how many energy levels
   allocate(dfdt_e(N_orb), source = 0.0d0)
   allocate(dfdt_e2(N_orb), source = 0.0d0)
   allocate(dfdt_A(N_orb), source = 0.0d0)
   allocate(dfdt_B(N_orb), source = 0.0d0)
   allocate(dfdt_A2(N_orb), source = 0.0d0)
   allocate(dfdt_B2(N_orb), source = 0.0d0)
   allocate(distre_temp(N_orb))
   allocate(distre_temp2(N_orb))
!    allocate(inv_t_e_ph(N_orb))
   N_steps = 1

   ! Deal with atomic distribution:
   N_tot_cut = 0.0d0
   coef = 2.0d0*g_e/g_h
   !coef_inv =g_h/( 2.0d0*g_e)   ! [s]
   coef_inv = numpar%M2_scaling*g_h/(g_e)   ! [s] with scaling factor added (e.g. 2 from d|a|^2/dt)
110 continue
   dt_small = numpar%dt*(1d-15)/dble(N_steps) ! [s] smaller time-steps to iterate with better precision
   dE_nonadiabat = 0.0d0
   distre_temp = distre ! to start
   distre_temp2 = distre_temp ! to start
   !E_tot_cur = Scell%E_tot
   E_tot_cur = nrg%El_low

   TIME_STEPS:do i_t = 1, N_steps
      ! initializing:
      dfdt_e = 0.0d0
      dfdt_e2 = 0.0d0
      dfdt_A = 0.0d0
      dfdt_B = 0.0d0
      dfdt_A2 = 0.0d0
      dfdt_B2 = 0.0d0
      i_pc = 1 ! counter of predictor-corrector iterations
      pc_iterations = .true.  ! predictor-corrector iterations needed
!     inv_t_e_ph = 0.0d0

      P_C:do while (pc_iterations)
         G_ei_partial = 0.0d0 ! to start with

         ! Get the collision integral parameters:
         !call get_el_ion_kernel(Scell, numpar, Mij, wr, wr0, dt_small, kind_M, distre_temp2, dfdt_e2, dfdt_A2, dfdt_B2, &
         !                        DOS_weights, G_ei_partial) ! below
         call get_el_ion_kernel(Scell, numpar, Mij, wr, wr0, numpar%dt*(1d-15), kind_M, distre_temp2, dfdt_e2, dfdt_A2, dfdt_B2, &
                                 DOS_weights, G_ei_partial) ! below

         ! first step, no iterative data yet, use what we have:
         if (i_pc == 1) then
            dfdt_e = dfdt_e2
            dfdt_A = dfdt_A2
            dfdt_B = dfdt_B2
         endif
         pc_iterations = .false. ! assume this iteration of predictor-corrector is the last one (to check below)

         ! Having the collision integral, update the distribution function:
         do i = 1, N_orb
            ! Explicit scheme:
            distre_temp2(i) = distre_temp(i) + dfdt_e(i)*dt_small
            ! Implicit scheme:
            !distre_temp(i) = (distre_temp(i) + 2.0d0*dfdt_A(i)*dt_small)/(1.0d0 + (dfdt_B(i) + dfdt_A(i))*dt_small)
            ! Explicit predictor-corrector:
            !distre_temp2(i) = distre_temp(i) + ( dfdt_e(i)*(1.0d0-alpha) + dfdt_e2(i)*alpha ) * dt_small
            ! Implicit predictor-corrector:
            !A_pc = dfdt_A(i)*(1.0d0-alpha) + dfdt_A2(i)*alpha
            !B_pc = dfdt_B(i)*(1.0d0-alpha) + dfdt_B2(i)*alpha
            !AB = A_pc + B_pc
            !if (abs(AB) < eps) then
            !   A2AB = 0.0d0
            !else
            !   A2AB = 2.0d0*A_pc / AB
            !endif
            !arg = dt_small*AB
            !distre_temp2(i) = A2AB - (A2AB - distre_temp(i))* exp(-dt_small*AB)
         
            ! Check if another iteration of predictor-corrector is required:
            !if ((distre_temp2(i) > 2.0d0) .OR. (distre_temp2(i) < 0.0d0)) then   ! obviously
               !if (i < 5) print*, 'Predictor-corrector problem:', i_pc, i, distre_temp2(i)
            !   pc_iterations = .true.
            !endif

            ! Reducing time-steps if needed:
            if ((distre_temp2(i) > 2.0d0) .OR. (distre_temp2(i) < 0.0d0)) then
               N_steps = INT(N_steps*2.0d0)  ! decrease time step in Boltzmann equation and recalculate
               !print*, 'N_steps=', N_steps, 'SOMETHING MIGHT BE WRONG (1)'
               !write(*,'(i5, e25.16, e25.16, e25.16, e25.16, e25.16, e25.16)') i, dt_small, coef, coef_inv, distre_temp2(i), dfdt_e(i), dfdt_e(i)*dt_small
               ! Recalculate with smaller time steps:
               goto 110
            endif
         enddo
         !i_pc = i_pc + 1   ! one more iteration
         ! for the next iteration:
         dfdt_e = dfdt_e2
         dfdt_A = dfdt_A2
         dfdt_B = dfdt_B2

         ! Predictor-corrector check:
         !N_tot_cut = SUM(distre_temp2(:))
         !if (ABS(N_tot_cut - Scell%Ne_low)/Scell%Ne_low > eps_N) then
         !   pc_iterations = .true.  ! another iteration required
         !   print*, 'PC_steps=', i_pc, 'SOMETHING MIGHT BE WRONG (1):', N_tot_cut, Scell%Ne_low
         !endif
         !if (i_pc > 10000) then ! maximal number of p-c iterations
         !   print*, 'Too many iterations in P-C'
            !print*, N_tot_cut, Scell%Ne_low
         !   pc_iterations = .false.
         !endif
      enddo P_C
      ! When we have an acceptable solution:
      distre_temp(:) = distre_temp2(:)

      G_ei_partial = G_ei_partial*dt_small

      N_tot_cut = SUM(distre_temp(:))
      !E_tot_cur = SUM(distre_temp(:)*wr(:))

      if (ABS(N_tot_cut - Scell%Ne_low)/Scell%Ne_low > eps_N) then
         N_steps = INT(N_steps*2.0d0)
         print*, 'Trouble in Electron_ion_collision_int:'
         print*, 'N_steps=', N_steps, 'SOMETHING MIGHT BE WRONG (2):'
         print*, N_tot_cut, Scell%Ne_low
!          print*, dfdt_e(:)
         goto 110
      endif
   enddo TIME_STEPS

   ! Update electron energy:
   nrg%El_low = SUM(distre_temp(:)*wr(:)) ! new total electron energy [eV]
   dE_nonadiabat = SUM( (distre(:) - distre_temp(:)) * wr(:) )    ! [eV] transfered energy
   !print*, 'dE:', dE_nonadiabat

   ! Save the collision integral, if needed:
   if (numpar%do_kappa) then
      Scell%I_ij(:) = (distre(:) - distre_temp(:))/(numpar%dt * 1d-15)  ! [1/s]
   endif

   distre = distre_temp ! update electron distribution function after energy exchange with the atomic system
   
   deallocate(dfdt_e, dfdt_A, dfdt_B, distre_temp)
   deallocate(dfdt_e2, dfdt_A2, dfdt_B2, distre_temp2)
!    deallocate(inv_t_e_ph)
end subroutine Electron_ion_collision_int




subroutine get_el_ion_kernel(Scell, numpar, Mij, wr, wr0, dt_small, kind_M, distre_temp, &
                              dfdt_e, dfdt_A, dfdt_B, DOS_weights, G_ei_partial)
   type(Super_cell), intent(in) :: Scell ! supercell with all the atoms as one object
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:), intent(in) :: Mij  ! Matrix element coupling electron energy levels via ion motion
   real(8), dimension(:), intent(in) :: wr, wr0  ! electron energy levels [eV], on current and last timesteps
   real(8), intent(in) :: dt_small  ! timestep [fs]
   integer, intent(in) :: kind_M	! index of the model used for the matri element calculation
   real(8), dimension(:), intent(in) :: distre_temp	! electron distribution
   real(8), dimension(:), intent(inout) :: dfdt_e	! change of the distribution function
   real(8), dimension(:), intent(inout) :: dfdt_A, dfdt_B ! factors for explicit scheme
   real(8), dimension(:,:,:), intent(in) :: DOS_weights     ! weigths of the particular type of orbital on each energy level
   real(8), dimension(:,:), intent(inout) :: G_ei_partial     ! partial contributions into coupling from different shells of
   !------------------
   real(8) :: dfdt_i, wij, ViVj, Mij2, distr_at_in, distr_at_fin, dtij, dfdt, dfdt_aij, dfdt_bij, coef, coef_inv
   integer :: i, j

   coef = 2.0d0*g_e/g_h
   coef_inv = numpar%M2_scaling*g_h/(g_e)   ! [s] with scaling factor added (e.g. 2 from d|a|^2/dt) OLD
   !coef_inv = numpar%M2_scaling   ! scaling factor added (e.g. 2 from d|a|^2/dt)
   ! Initializing for summation:
   dfdt_e = 0.0d0
   dfdt_A = 0.0d0
   dfdt_B = 0.0d0

!$omp PARALLEL private(i,j,dfdt_i,ViVj,Mij2,distr_at_in,distr_at_fin, dtij, wij, dfdt, dfdt_aij, dfdt_bij) ! Implicit scheme
!$omp do reduction( + : dfdt_e, G_ei_partial, dfdt_A, dfdt_B) ! explicit scheme
    do i = 1, size(wr) ! all energy levels
      dfdt_i = 0.0d0  ! implicit
      do j = 1, size(wr) ! all pairs of energy levels
          ! All transitions except degenerate levels
          wij = abs( wr(i) - wr(j) )    ! [eV] energy difference between the levels = transferred energy
          ! Exclude too large jumps, and quasidegenerate levels:
          large_jumps:if ( (wij < numpar%acc_window) .and. (wij > numpar%degeneracy_eV) ) then
               ! Select which scheme for the matrix element to use:
               select case (kind_M)
               case (-1)    ! Landau nonperturbative expression (with ad hoc correction):
                  dtij = coef_inv !/ wij
                  Mij2 = Mij(i,j)*Mij(i,j) / dt_small * dtij/dt_small
               case(0)
                  Mij2 = 0.0d0 ! no coupling
               case(2)  ! Fermi's golden rule:
                  ViVj = (wr(i)+wr0(i) - (wr(j)+wr0(j)))/2.0d0 ! for matrix element:
                  Mij2 = Mij(i,j)*Mij(i,j)*ViVj*ViVj*g_Pi*coef
               case(3)  ! incomplete Fermi golden rule:
                  ViVj = (wr(i)+wr0(i) - (wr(j)+wr0(j)))/2.0d0 ! for matrix element:
                  Mij2 = Mij(i,j)*Mij(i,j)*ViVj*sin(ViVj*g_e/g_h*dt_small)*coef
               case default ! Tully's expression:
                  Mij2 = Mij(i,j)*Mij(i,j)/dt_small !/dble(N_steps)
               end select

               ! analytical integration of maxwell function:
               if (i > j) then
                  distr_at_fin = Maxwell_int_shifted(Scell%TaeV, wij)  ! from module "Atomic_tools"
                  distr_at_in = 1.0d0	! all atoms can absorb this energy
               else
                  distr_at_fin = 1.0d0	! all atoms can absorb this energy
                  distr_at_in = Maxwell_int_shifted(Scell%TaeV, wij)  ! from module "Atomic_tools"
               endif

               ! Explicit scheme (reqiures small timestep):
               dfdt = Mij2*( distre_temp(j)*distr_at_fin*(2.0d0-distre_temp(i)) - distre_temp(i)*distr_at_in*(2.0d0-distre_temp(j)) )

               ! Implicit scheme (unconditionally stable, but does not conserve number of electrons!):
               dfdt_aij = Mij2*distr_at_fin
               dfdt_bij = Mij2*distr_at_in
               ! Sum the contributions up:
               dfdt_A(i) = dfdt_A(i) + dfdt_aij*distre_temp(j)
               dfdt_B(i) = dfdt_B(i) + dfdt_bij*(2.0d0-distre_temp(j))

               !if (abs(dfdt * dt_small) < 2.0d0) then
                   dfdt_i = dfdt_i + dfdt
                   ! Get contributions from different shells:
                   call Get_Gei_shells_contrib(i, j, DOS_weights, dfdt, wr, G_ei_partial)    ! below
               !endif ! (abs(dfdt * dt_small) < 2.0d0)


          endif large_jumps
      enddo ! j
      ! Implicit scheme:
      dfdt_e(i) = dfdt_i
    enddo ! i
!$omp end do
!$omp end parallel

end subroutine get_el_ion_kernel



! Get contributions from different shells analogously to Mulliken charge calculations:
pure subroutine Get_Gei_shells_contrib(i, j, DOS_weights, dfdt, wr, G_ei_partial)
   integer, intent(in) :: i, j
   real(8), dimension(:,:,:), intent(in) :: DOS_weights
   real(8), intent(in) :: dfdt
   real(8), dimension(:), intent(in) :: wr
   real(8), dimension(:,:), intent(out) :: G_ei_partial
   !-----------------------------
   integer :: i_at, N_at, i_types, N_types, i_at2, i_types2, i_G1, i_G2
   real(8) :: G_temp
   
   N_at = size(DOS_weights,1)    ! number of kinds of atoms
   N_types = size(DOS_weights,2) ! number of atomic shells (basis set size)
!    Nsiz = size(DOS_weights,3)    ! number of energy levels

   do i_at = 1, N_at
      do i_types = 1, N_types
         i_G1 = (i_at-1) * N_types + i_types
         G_temp = dfdt * DOS_weights(i_at, i_types,i) * wr(i)
         do i_at2 = 1, N_at
            do i_types2 = 1, N_types
               i_G2 = (i_at2-1) * N_types + i_types2
               G_ei_partial(i_G1,i_G2) = G_ei_partial(i_G1,i_G2) - G_temp * DOS_weights(i_at2, i_types2,j)  ! sign defined as energy loss
            enddo   ! i_types2
         enddo ! i_at2
      enddo   ! i_types
   enddo ! i_at
   
end subroutine Get_Gei_shells_contrib


subroutine Electron_ion_coupling_Mij(wr, Ha, Ha0, Mij, kind_M, Sij) ! calculates the electron-ion coupling matrix element
   real(8), dimension(:), intent(in) ::  wr ! electron energy levels [eV], after diagonalization of Ha
   real(8), dimension(:,:), intent(in) :: Ha  ! diagonilized Hamiltonian Ha, eigenvectors
   real(8), dimension(:,:), intent(in) :: Ha0 ! eigenvectors on the previous time-step
   real(8), dimension(:,:), allocatable, intent(inout) :: Mij  ! Matrix element coupling electron WF via ion motion
!    real(8), dimension(:), intent(out) ::  Norm_WF   ! WF normalization coefficient
   integer, intent(in), optional :: kind_M ! what kind of matrix element to use
   real(8), dimension(:,:), intent(in), optional :: Sij ! overlap matrix, in case of non-orthogonal basis set
   !----------------------------
   integer i, j, N, k
   real(8), dimension(:), allocatable :: psi0 ! WF_0
   real(8), dimension(:), allocatable :: psi_back ! WF
   real(8), dimension(size(wr))  :: Norm0, Norm1   ! normalization factors
   real(8), dimension(size(wr),size(wr)) :: tij
   real(8) :: eps, Norm_val0, Norm_val1
   
!    print*, 'Electron_ion_coupling_Mij start'
   
   N = size(wr)
   if (.not. allocated(Mij)) allocate(Mij(N,N))
   Mij = 0.0d0
   eps = 1.0d-13    ! acceptable error
   
!    print*, 'Electron_ion_coupling_Mij 0'
   
   if (present (kind_M)) then
      select case (kind_M)
      case(-1) ! the full probability by Landau
         
!          print*, 'Electron_ion_coupling_Mij 1'
         
         !$omp PARALLEL private(i, j, k, psi0, psi_back, Norm_val1) shared(Norm0, Norm1, tij)
         allocate(psi0(N))
         allocate(psi_back(N))
         
!          print*, 'Electron_ion_coupling_Mij 2'
         
         ! Ensure WF normalization to 1: 
         !$omp do
         do i = 1, N 
            Norm0(i) = DSQRT(SUM( Ha0(:,i) * Ha0(:,i) ))
            Norm1(i) = DSQRT(SUM( Ha(:,i) * Ha(:,i) ))
         enddo
         !$omp end do
         !$OMP BARRIER
         
!          print*, 'Electron_ion_coupling_Mij 3'
         
         if (present(Sij)) then ! non-orthogonal
!             print*, 'Test 0 '
            !$omp do
            do i = 1, N ! all energy levels
               do k = 1, N
                  if (Norm1(i) < eps) then
                     tij(i,k) = 0.0d0
                  else
                     tij(i,k) = SUM(Sij(k,:)*Ha(:,i) / Norm1(i)) ! correct
                  endif
               enddo ! k = 1, N
            enddo !  i = 1, N
            !$omp end do
            !$OMP BARRIER
!             print*, 'Electron_ion_coupling_Mij 4'
         endif

         !$omp do
         do i = 1, N ! all energy levels
!             psi0(:) = Ha0(:,i)
            if (Norm0(i) < eps) then
               psi0(:) = 0.0d0
            else
               psi0(:) = Ha0(:,i) / Norm0(i)
            endif
            
!             print*, 'Electron_ion_coupling_Mij 5'
            
!             psi_back(:) = Ha(:,i)
            do j = i, N ! only upper triangle, using symmetry
!             do j = 1, N
               if (i .NE. j) then
                  if ( Norm1(j) < eps ) then
                     Norm_val1 = 0.0d0
                  else 
                     Norm_val1 = 1.0d0/Norm1(j)
                  endif
                  
!                   print*, 'Electron_ion_coupling_Mij 6'
                  
                  if (present(Sij)) then ! non-orthogonal
!                      print*, 'Test 1 ', i, j
                     Mij(i,j) = SUM( psi0(:)* tij(j,:) )
                  else ! orthogonal
                     Mij(i,j) = SUM( psi0(:) * Ha(:,j)*Norm_val1 )
                  endif ! (present(Sij))
                  
!                   print*, 'Electron_ion_coupling_Mij 7'

                  ! Use symmetry:
                  Mij(j,i) = Mij(i,j)
               endif
            enddo
         enddo
         !$omp end do !!end PARALLEL do
         deallocate(psi0,psi_back)
         !$omp end parallel
         
!          print*, 'Electron_ion_coupling_Mij 8'
         
      case default ! Tully's expression
         !$omp PARALLEL private(i, j, psi0, psi_back, Norm_val1, Norm_val0) shared(Norm0, Norm1, tij)
         allocate(psi0(N))
         allocate(psi_back(N))
         ! Ensure WF normalization to 1: 
         !$omp do
         do i = 1, N 
            Norm0(i) = DSQRT(SUM( Ha0(:,i) * Ha0(:,i) ))
            Norm1(i) = DSQRT(SUM( Ha(:,i) * Ha(:,i) ))
         enddo
         !$omp end do
         !$OMP BARRIER

         if (present(Sij)) then ! non-orthogonal
            !$omp do
            do i = 1, N ! all energy levels
               do k = 1, N
                  if (Norm1(i) < eps) then
                     tij(i,k) = 0.0d0
                  else
                     tij(i,k) = SUM(Sij(k,:)*Ha(:,i) / Norm1(i)) ! correct
                  endif
               enddo ! k = 1, N
            enddo !  i = 1, N
            !$omp end do
            !$OMP BARRIER
         endif
         
         !$omp do
         do i = 1, N ! all energy levels
            if (Norm0(i) < eps) then
               psi0(:) = 0.0d0
            else
               psi0(:) = Ha0(:,i) / Norm0(i)
            endif
            if (Norm1(i) < eps) then
               psi_back(:) = 0.0d0
            else
               psi_back(:) = Ha(:,i) / Norm1(i)
            endif
            
            do j = i, N
               if (i .NE. j) then !First Born approximation:
                  if ( Norm1(j) < eps ) then
                     Norm_val1 = 0.0d0
                  else
                     Norm_val1 = 1.0d0/Norm1(j)
                  endif
                  if ( Norm0(j) < eps ) then
                     Norm_val0 = 0.0d0
                  else
                     Norm_val0 = 1.0d0/Norm0(j)
                  endif
                  if (present(Sij)) then ! non-orthogonal
                     Mij(i,j) = SUM( psi0(:)* tij(j,:)  - (Ha0(:,j)*Norm_val0)* tij(i,:) )
                     Mij(j,i) = Mij(i,j)
                  else ! orthogonal
                     Mij(i,j) = SUM(psi0(:)*(Ha(:,j)*Norm_val1) - (Ha0(:,j)*Norm_val0)*psi_back(:) )
                     Mij(j,i) = Mij(i,j)
                  endif ! (present(Sij))
!                    if (j==129) print*, 'M', i, j, Mij(i,j)
               endif ! (i .NE. j)
            enddo
         enddo
         !$omp end do !!end PARALLEL do
         deallocate(psi0, psi_back)
         !$omp end parallel
         Mij = Mij*0.5d0
      end select
      
!       Norm_WF = Norm1   ! save normalization coefficient
   endif ! present(kind_M)

!    print*, 'Electron_ion_coupling_Mij end'

end subroutine Electron_ion_coupling_Mij





subroutine Electron_ion_coupling_Mij_OLD(wr, Ha, Ha0, Mij, kind_M, Sij) ! calculates the electron-ion coupling matrix element
   real(8), dimension(:), intent(in) ::  wr ! electron energy levels [eV], after diagonalization of Ha
   real(8), dimension(:,:), intent(in) :: Ha  ! diagonilized Hamiltonian Ha, eigenvectors
   real(8), dimension(:,:), intent(in) :: Ha0 ! eigenvectors on the previous time-step
   real(8), dimension(:,:), allocatable, intent(out) :: Mij  ! Matrix element coupling electron WF via ion motion
   integer, intent(in), optional :: kind_M ! what kind of matrix element to use
   real(8), dimension(:,:), intent(in), optional :: Sij ! overlap matrix, in case of non-orthogonal basis set
   !----------------------------
   integer i, j, N
   real(8), dimension(:), allocatable :: psi0 ! WF_0
   real(8), dimension(:), allocatable :: psi_back ! WF
   real(8), dimension(size(wr))  :: Norm0, Norm1   ! normalization factors
   real(8), dimension(:), allocatable :: tij, til
   real(8) :: eps, Norm_val0, Norm_val1
   N = size(wr)
   if (.not. allocated(Mij)) allocate(Mij(N,N))
   Mij = 0.0d0
   eps = 1.0d-13    ! acceptable error
   
   if (present (kind_M)) then
      select case (kind_M)
      case(-1) ! the full probability by Landau
         !$omp PARALLEL private(i,j,psi0, psi_back)
         allocate(psi0(N))
         allocate(psi_back(N))
         !$omp do
         do i = 1, N ! all energy levels
            psi0(:) = Ha0(i,:)
            psi_back(:) = Ha(i,:)
            do j = i, N
               if (i .NE. j) then !First Born approximation:
                  Mij(i,j) = SUM(psi0(:)*Ha(j,:))
                  Mij(i,j) = ABS(Mij(i,j))*(Mij(i,j) -  SUM(Ha0(j,:)*psi_back(:)))
                  Mij(j,i) = Mij(i,j)
               endif
            enddo
         enddo
         !$omp end do !!end PARALLEL do
         deallocate(psi0)
         deallocate(psi_back)
         !$omp end parallel
      case default ! Tully's expression
         !$omp PARALLEL private(i, j, psi0, psi_back, Norm_val1, Norm_val0) shared(Norm0, Norm1)
         
         allocate(psi0(N))
         allocate(psi_back(N))
         ! Ensure WF normalization to 1: 
         !$omp do
         do i = 1, N 
            Norm0(i) = DSQRT(SUM( Ha0(i,:) * Ha0(i,:) ))
            Norm1(i) = DSQRT(SUM( Ha(i,:) * Ha(i,:) ))
         enddo
         !$omp end do
         !$OMP BARRIER

         !$omp do
         do i = 1, N ! all energy levels
!             psi0(:) = Ha0(i,:)
!             psi_back(:) = Ha(i,:)
            if (Norm0(i) < eps) then
               psi0(:) = 0.0d0
            else
               psi0(:) = Ha0(i,:) / Norm0(i)
            endif
            if (Norm1(i) < eps) then
               psi_back(:) = 0.0d0
            else
               psi_back(:) = Ha(i,:) / Norm1(i)
            endif
            
            do j = i, N
               if (i .NE. j) then !First Born approximation:
                  if ( Norm1(j) < eps ) then
                     Norm_val1 = 0.0d0
                  else 
                     Norm_val1 = 1.0d0/Norm1(j)
                  endif
                  if ( Norm0(j) < eps ) then
                     Norm_val0 = 0.0d0
                  else 
                     Norm_val0 = 1.0d0/Norm0(j)
                  endif
!                   Mij(i,j) = SUM(psi0(:)*(Ha(j,:)/Norm1(j)) - (Ha0(j,:)/Norm0(j))*psi_back(:))
                  Mij(i,j) = SUM(psi0(:)*(Ha(j,:)*Norm_val1) - (Ha0(j,:)*Norm_val0)*psi_back(:) )
                  Mij(j,i) = Mij(i,j)
               endif
            enddo
         enddo
         !$omp end do !!end PARALLEL do
         deallocate(psi0, psi_back)
         !$omp end parallel
         Mij = Mij*0.5d0
      end select
   endif ! present(kind_M)

end subroutine Electron_ion_coupling_Mij_OLD


subroutine Electron_ion_coupling_Mij_complex(Ha_non, CHij, CHij0, Mij,  ix, iy, iz, ixm, iym, izm, kind_M) ! calculates the electron-ion coupling matrix element
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha_non	! non-diagonilized Hamiltonian Ha
   complex, dimension(:,:,:,:,:), allocatable, INTENT(inout) :: CHij	! Complex hamiltonian matrix for each (kx, ky, kz) points
   complex, dimension(:,:,:,:,:), allocatable, INTENT(inout) :: CHij0	! Complex hamiltonian matrix on the last step
   REAL(8), DIMENSION(:,:), allocatable, INTENT(out) :: Mij	! Matrix element coupling electron WF via ion motion
   integer, intent(in) ::  ix, iy, iz, ixm, iym, izm	! current indices of k-points at which the Hamiltonian and ensuing things are calculated
   integer, intent(in), optional :: kind_M ! what kind of matrix element to use
   !----------------------------
   real(8) :: kx, ky, kz
   integer :: N
   ! k-points choice from [J. Moreno, J. M. Soler, PRB 45, 13891 (1992)]
   kx = (2.0d0*dble(ix) - dble(ixm) - 1.0d0)/(2.0d0*dble(ixm))
   ky = (2.0d0*dble(iy) - dble(iym) - 1.0d0)/(2.0d0*dble(iym))
   kz = (2.0d0*dble(iz) - dble(izm) - 1.0d0)/(2.0d0*dble(izm))
   print*, 'K:', ix, iy, iz, kx, ky, kz
   
   ! size of the arrays:
   N = size(Ha_non,1)
   if (.not. allocated(Mij)) allocate(Mij(N,N))
   Mij = 0.0d0
   
   if (.not.allocated(CHij0)) then ! that means, it the first run, no coupling may be calculated here, set it to zero:
      allocate(CHij0(ixm, iym, izm, N, N))
      ! Get the values of the complex Hamiltonian to use at the next time-step:
      
      return		! no need to continue, no coupling can be calculated at the first step
   endif
   
   if (.not.allocated(CHij)) then
      allocate(CHij(ixm, iym, izm, N, N))   
   endif
   
   CHij0 = CHij	! save it for the next step
end subroutine Electron_ion_coupling_Mij_complex



subroutine get_G_ei(Scell, NSC, numpar, dE)
   type(Super_cell), dimension(:), intent(inout) :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), intent(in) :: dE	! [eV] energy transferred from electrons to ions
   real(8) :: temp
   temp = g_e/((Scell(NSC)%Te - Scell(NSC)%Ta)*Scell(NSC)%V*1d-30*(numpar%dt*1d-15))	! [W/(m^3 K)]
   !Scell(NSC)%G_ei = dE*g_e/((Scell(NSC)%Te - Scell(NSC)%Ta)*Scell(NSC)%V*1d-30*(numpar%dt*1d-15))	! [W/(m^3 K)]
   Scell(NSC)%G_ei = dE*temp	! [W/(m^3 K)]
   if (numpar%DOS_splitting == 1) then ! get partial G
      Scell(NSC)%G_ei_partial = Scell(NSC)%G_ei_partial*temp
   endif
end subroutine get_G_ei 




!pppppppppppppppppppppppp
! Project specific subroutines:
subroutine Landau_vs_Sakurei_test(Mij, Ha, Ha0, Ha_non, Ha_non0, wr, wr0)
   REAL(8), DIMENSION(:,:), INTENT(in) :: Mij  ! Matrix element coupling electron WF via ion motion
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha   ! diagonilized Hamiltonian Ha, eigenvectors
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha0  ! eigenvectors on the previous time-step
   REAL(8), DIMENSION(:,:), INTENT(in) :: Ha_non	! non-diagonilized Hamiltonian
   REAL(8), DIMENSION(:,:), INTENT(in), optional :: Ha_non0	! non-diagonilized Hamiltonian at last time-step
   REAL(8), DIMENSION(:), INTENT(in) ::  wr ! electron energy levels [eV], after diagonalization of Ha
   REAL(8), DIMENSION(:), INTENT(in) ::  wr0 ! electron energy levels [eV], after diagonalization of Ha, last time step
   real(8) :: diff, Ev
   real(8), dimension(size(Ha,1), size(Ha,2)) :: Ha_1, Ha_new, Ha_non1, Vmn
   real(8), dimension(size(Ha,1)) :: Evec, E_2
   character(200) :: Error_descript, File_name
   integer :: N_lev, n, m, i, k, FN
   logical :: file_opened
   Ha_1 = 0.0d0
   N_lev = size(Ha,1)
   
   if (present(Ha_non0)) then
!       Ha_non1 = Ha_non0
      ! Diagonalize the Hamiltonian to get electron energy levels:
!       call sym_diagonalize(Ha_non1, wr0, Error_descript, check_M=.true.)
     
     File_name = 'OUTPUT_Landau_energy_shifts.dat'
     FN = 7755
     inquire(file=trim(adjustl(File_name)),opened=file_opened)
     if (.not.file_opened) open(UNIT=FN, FILE = trim(adjustl(File_name)), status = 'unknown', position="append", action="write")
   
     
      do n = 1, N_lev
         do i = 1, N_lev
            if (n < 2) then
               Ha_1(n,i) = SUM((Ha_non(n+1:,i) - Ha_non0(n+1:,i))*Ha0(n+1:,i)/(wr0(n) - wr0(n+1:)))
            else if (n > N_lev - 1) then
               Ha_1(n,i) = SUM((Ha_non(1:n-1,i) - Ha_non0(1:n-1,i))*Ha0(1:n-1,i)/(wr0(n) - wr0(1:n-1)))
            else
               Ha_1(n,i) = SUM((Ha_non(1:n-1,i) - Ha_non0(1:n-1,i))*Ha0(1:n-1,i)/(wr0(n) - wr0(1:n-1))) + SUM((Ha_non(n+1:,i) - Ha_non0(n+1:,i))*Ha0(n+1:,i)/(wr0(n) - wr0(n+1:)))
            endif
         enddo
      enddo
      
      Ha_new = Ha0 + Ha_1
      
      do n = 1, N_lev
         diff = (SUM(Ha_new(n,:)*Ha_new(n,:)))
         Ha_new(n,:) = Ha_new(n,:)/SQRT(diff)	! normalizing to 1
      enddo
      
!       Ev = 0.0d0
!       do i = 1,N_lev
!          do k = 1,N_lev
!             !Evec(k) = SUM(Ha_non0(k,:)*Ha_new(:,i)) ! correct
!             Evec(k) = SUM((Ha_non(k,:) - Ha_non0(k,:))*Ha_new(:,i)) ! correct
!          enddo
!          Ev = SUM(Evec(:)*Ha_new(:,i))
!          write(*,'(i4,f,f,f,f)') i, (wr(i)-wr0(i))/abs(wr0(i)), Ev/abs(wr0(i)), Ev/(wr(i)-wr0(i))
!       enddo
      
      E_2 = 0.0d0 ! second order correction to the energy
      do i = 1,N_lev
         do k = 1,N_lev
            Vmn(k,i) = SUM( (Ha_non(k,:) - Ha_non0(k,:))*Ha(i,:)*Ha(k,:) )
         enddo
      enddo
      do i = 1,N_lev
         E_2(i) = SUM( Vmn(i,:)*Vmn(i,:)/(wr0(i) - wr0(:)), MASK = (wr0(:) /= wr0(i)) )
      enddo
      
      
      Ev = 0.0d0 ! first order correction to the energy
      do i = 1,N_lev
         do k = 1,N_lev
            !Evec(k) = SUM(Ha_non0(k,:)*Ha_new(:,i)) ! correct
            Evec(k) = SUM((Ha_non(k,:) - Ha_non0(k,:))*Ha(:,i)) ! correct
         enddo
         Ev = SUM(Evec(:)*Ha(:,i))
         write(*,'(i4,f25.16,f25.16,f25.16,f25.16)') i, (wr(i)-wr0(i))/abs(wr0(i)), Ev/abs(wr0(i)), Ev/(wr(i)-wr0(i)), (Ev+E_2(i))/(wr(i)-wr0(i))  
         write(FN,'(f25.16)') Ev/(wr(i)-wr0(i))
      enddo
      write(FN,'(a)') '' ! break
      
      PAUSE 'Landau_vs_Sakurei_test'
         
   else ! do the test with precalculated matrix elements:
      diff = 0.0d0
      do n = 1, N_lev
         diff = 0.0d0
         do m = 1, N_lev
            if (n < 2) then
               Ha_1(n,m) = SUM(Mij(n+1:,m)*Ha0(n+1:,m) )
            else if (n > N_lev - 1) then
               Ha_1(n,m) = SUM(Mij(1:n-1,m)*Ha0(1:n-1,m) )
            else
               Ha_1(n,m) = SUM(Mij(1:n-1,m)*Ha0(1:n-1,m) ) + SUM(Mij(n+1:,m)*Ha0(n+1:,m) )
            endif
            diff = diff + ( Ha(n,m) - (Ha0(n,m)+Ha_1(n,m)) )
         enddo ! j
!          print*, n, diff
      enddo ! i
      Ha_new = Ha0 + Ha_1
      do n = 1, N_lev
         diff = (SUM(Ha_new(n,:)*Ha_new(n,:)))
         Ha_new(n,:) = Ha_new(n,:)/SQRT(diff)		! normalizing to 1
      enddo
      
      Ev = 0.0d0
      do i = 1,N_lev
         do k = 1,N_lev
            Evec(k) = SUM(Ha_non(k,:)*Ha_new(:,i)) ! correct
         enddo
         Ev = SUM(Evec(:)*Ha_new(:,i))
         write(*,'(i4,f25.16,f25.16,f25.16)') i, (wr(i)-wr0(i))/abs(wr0(i)), (Ev-wr0(i))/abs(wr0(i)), (wr(i)-wr0(i))/(Ev-wr0(i))
      enddo
   
      PAUSE 'Landau_vs_Sakurei_test'
   endif
end subroutine Landau_vs_Sakurei_test


END MODULE Nonadiabatic
