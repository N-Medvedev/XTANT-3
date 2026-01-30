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
MODULE Electron_electron_scattering
use Universal_constants
use Objects
use Algebra_tools, only : Two_Vect_Matr
use Little_subroutines, only : Find_in_array_monoton, Fermi_interpolation, linear_interpolation, Gaussian
use Atomic_tools, only : shortest_distance


implicit none
PRIVATE

public :: get_Boltzmann_alpha_beta, Boltzmann_solution, test_change_of_fe, Electron_electron_scattering_Kij


 contains


!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! Nonequilibrium electron kinetics


subroutine get_Boltzmann_alpha_beta(Scell, i, Ev, K_ij, fe, dt, alpha_ij, beta_ij)
   type(Super_cell), intent(in), target :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: i   ! current energy level
   real(8), dimension(:), intent(in), target :: Ev    ! [eV] electron energy levels
   real(8), dimension(:,:), intent(in) :: K_ij        ! electron-electron scattering matrix elements (probabilities)
   real(8), dimension(:), intent(in) :: fe            ! electron distribution function on the last time-step
   real(8), intent(in) :: dt        ! time step [fs]
   real(8), dimension(:,:), intent(inout) :: alpha_ij, beta_ij
   !--------------------------
   integer :: j, i_f, j_g, N, N_basis_set, j_cur, i_f_cur
   integer :: i_at, j_at, i_f_at, j_g_at, atom_2, atom_3
   real(8) :: Ej_g, fe_g, alpha_ij_temp, beta_ij_temp, w_ijfg
   real(8) :: rij, rif, rjg, r_ave
   logical :: within_range ! check if final energy level is within possible range
   real(8), pointer :: Ei, Ej, Ei_f
   integer, pointer :: m, mf

   N = size(fe)   ! Number of electrons (normalization of fe)
   Ei => Ev(i)    ! Electron energy levels [eV]
   ! Size of the basis set, to identify atoms:
   N_basis_set = N/Scell%Na
   ! Index of the atom, to which the level "i" belongs:
   if (i < N) then
      i_at = 1+i/N_basis_set
   else
      i_at = i/N_basis_set
   endif
   ! Use symmetry: only upper triangle of (i,j)-matrix
   !E_FROM:do j = i, N ! all levels from where second electron can scatter off
   m => Scell%Near_neighbor_size(i_at) ! how many atoms we can interact with
   E_FROM:do atom_2 = 0,m ! do only for atoms close to that one, which can interact
      ! Find the distance between atoms i and j, between the two scattering electrons:
      if (atom_2 == 0) then ! self atom
         j_at = i_at
         rij = 0.0d0
      else ! different atoms
         j_at = Scell%Near_neighbor_list(i_at,atom_2) ! this is the list of such close atoms
         rij = Scell%Near_neighbor_dist(i_at,atom_2,4) ! at this distance, R [A]
      endif
      do j_cur = 1, N_basis_set  ! all orbitals of this atom
         j = (j_at-1)*N_basis_set + j_cur

         Ej => Ev(j) ! [eV]

         ! Get the inner sum over j (and g):
         alpha_ij_temp = 0.0d0   ! to start with
         beta_ij_temp = 0.0d0    ! to start with

         !E_TO:do i_f = 1, N ! all levels where first electron can jump to
         mf => Scell%Near_neighbor_size(i_at) ! how many atoms we can interact with
         E_TO:do atom_3 = 0,mf ! do only for atoms close to that one, which can interact
            ! Find the distance between atoms j and g, between which the 1st electron scatteres:
            if (atom_3 == 0) then ! self atom
               i_f_at = i_at
               rif = 0.0d0
            else ! different atoms
               i_f_at = Scell%Near_neighbor_list(i_at,atom_3) ! this is the list of such close atoms
               rif = Scell%Near_neighbor_dist(i_at,atom_3,4) ! at this distance, R [A]
            endif

            do i_f_cur = 1, N_basis_set  ! all orbitals of this atom
               ! Corresponding energy level:
               i_f = (i_f_at-1)*N_basis_set + i_f_cur

               Ei_f => Ev(i_f) ! [eV] final level of 1st electron
               Ej_g = Ej - (Ei_f - Ei) ! [eV] final energy (effective level) of 2d electron

               ! Check if such a transition is possible:
               within_range = .true. ! to start with
               ! Exclude transitions outside of energy levels given within TB:
               if (Ei_f < Ev(1) .or. Ei_f > Ev(N)) then
                  within_range = .false.
               endif
               ! Same for the second electron:
               if (Ej_g < Ev(1) .or. Ej_g > Ev(N)) then
                  within_range = .false.
               endif

               if (within_range) then ! if transition is possible
                  ! Find the value of distribution in-between the energy levels (final state of 2d electron):
                  call mean_distribution(Ev, fe, Ej_g, fe_g, i_f1=j_g, method=5) ! see below

                  if (j_g > 0) then ! scattering is allowed:
                     if (j_g < N) then
                        j_g_at = 1 + (j_g/N_basis_set)
                     else
                        j_g_at = (j_g/N_basis_set)
                     endif

                     ! Find the distance between atoms j and g, between which the 2d electron scatteres:
                     call shortest_distance(Scell, i_f_at, j_g_at, rjg)  ! module "Atomic_tools"

                     ! Scattering probability (i,j) -> (f,g):
                     !w_ijfg = 1.0d-3   ! [eV] just to test
                     r_ave = (rij+rif+rjg)
                     if (r_ave < 1.0d-5) r_ave = g_a0   ! same atom, eliminate zero distance
                     w_ijfg = ( (3.0d0*g_e*g_ke/(r_ave)*1.0d10) * K_ij(i,i_f)*K_ij(j,j_g) )**2    ! [eV] just to test

                     ! testing:
                     !if (w_ijfg > 1.0d-20) print*, i, i_f, j, j_g, (3.0d0*g_e*g_ke/(r_ave)*1.0d10), w_ijfg, K_ij(i,i_f), K_ij(j,j_g)

                     ! coefficients:
                     alpha_ij_temp = alpha_ij_temp + fe(i_f) * fe_g * w_ijfg
                     beta_ij_temp = beta_ij_temp + (2.0d0 - fe(i_f)) * (2.0d0 - fe_g) * w_ijfg
                  endif

               endif ! within_range

            enddo ! i_f_cur
         enddo E_TO
      enddo ! j_cur
      alpha_ij(i,j) = alpha_ij_temp ! save the coefficient for given (i,j)
      beta_ij(i,j) = beta_ij_temp   ! save the coefficient for given (i,j)
   enddo E_FROM

   nullify(Ei, Ej, Ei_f, m, mf)

   !pause 'get_Boltzmann_alpha_beta'
end subroutine get_Boltzmann_alpha_beta


! Electron-electron scatterign integral:
subroutine Boltzmann_solution(i, fe_temp, fe, dt, alpha_ij, beta_ij, scheme)
   integer, intent(in) :: i   ! current energy level
   real(8), dimension(:), intent(inout) :: fe_temp    ! electron distribution function on the curent time-step
   real(8), dimension(:), intent(in) :: fe            ! electron distribution function on the last time-step
   real(8), intent(in) :: dt        ! time step [fs]
   real(8), dimension(:,:), intent(in) :: alpha_ij, beta_ij ! Boltzmann alpha, beta coeffs
   integer, intent(in) :: scheme    ! which integration scheme in the Boltzmann equation to use: 0=full, 1=1st order
   !--------------------------
   integer :: j, N
   real(8) :: A, B, prefac, eps, A_cur, B_cur, a_com, b_com
   real(8) :: A_cur_0, B_cur_0, A_cur_1, B_cur_1, A_cur_2, B_cur_2, A_temp, B_temp
   logical :: within_range ! check if final energy level is within possible range
   
   eps = 1.0d-20  ! scattering allowed window [1/fs]

   N = size(fe)   ! Number of electrons (normalization of fe)

   prefac = 2.0d0 * g_e * 1.0d-15 / (g_h)  ! [1/fs] prefactor

   A = 0.0d0 ! to start with
   B = 0.0d0 ! to start with
   A_cur = 0.0d0 ! to start with
   B_cur = 0.0d0 ! to start with
   ! To avoid the problem of small and large numbers:
   A_cur_0 = 0.0d0 ! to start with
   B_cur_0 = 0.0d0 ! to start with
   A_cur_1 = 0.0d0 ! to start with
   B_cur_1 = 0.0d0 ! to start with
   A_cur_2 = 0.0d0 ! to start with
   B_cur_2 = 0.0d0 ! to start with

   ! Use precalculated alpha nad  beta:
   do j = 1, N ! all levels from where second electron can scatter off
      ! Naive addition:
      !A_cur = A_cur + (2.0d0 - fe(j)) * alpha_ij(i,j)
      !B_cur = B_cur + fe(j) * beta_ij(i,j)

      ! A little more involved addition:
      ! Temporary value for the given j:
      A_temp = (2.0d0 - fe(j)) * alpha_ij(i,j)
      B_temp = fe(j) * beta_ij(i,j)

      ! Add it together, accounting for possibilities of small and large numbers:
      if (abs(A_temp) > 1.0d-6) then  ! largest contributions
         A_cur_0 = A_cur_0 + A_temp
      elseif (abs(A_temp) > 1.0d-14) then   ! medium contributions
         A_cur_1 = A_cur_1 + A_temp
      else ! smallest contribution, not to lose
         A_cur_2 = A_cur_2 + A_temp
      endif
      ! Same for beta:
      if (abs(B_temp) > 1.0d-6) then  ! largest contributions
         B_cur_0 = B_cur_0 + B_temp
      elseif (abs(B_temp) > 1.0d-14) then   ! medium contributions
         B_cur_1 = B_cur_1 + B_temp
      else ! smallest contribution, not to lose
         B_cur_2 = B_cur_2 + B_temp
      endif
      ! Combine them together:

      ! Test symmetry:
!       a_com = max(alpha_ij(i,j),alpha_ij(j,i))
!       b_com = max(beta_ij(i,j),beta_ij(j,i))
!       if ( (a_com > 1.0d-10) .and. (abs(alpha_ij(i,j)-alpha_ij(j,i)) > 1.0d-6*a_com) ) then
!          print*, 'Problem in Boltzmann_solution (1a):', i, j, alpha_ij(i,j), alpha_ij(j,i)
!       endif
!       if ( (b_com > 1.0d-10) .and. (abs(beta_ij(i,j) - beta_ij(j,i)) > 1.0d-6*b_com) ) then
!          print*, 'Problem in Boltzmann_solution (1b):', i, j, beta_ij(i,j), beta_ij(j,i)
!       endif
   enddo ! j
   ! Add the final contributions together:
   !A = A_cur
   !B = B_cur
   A = A_cur_0 + A_cur_1 + A_cur_2
   B = B_cur_0 + B_cur_1 + B_cur_2

   !write(*,'(i0,f,es,es,es,f,f)') i, exp(-prefac*(A+B)*dt), A, B, A/(A+B), fe_temp(i), fe(i)

   if ( abs(A+B) > eps ) then ! scattering possible:
      select case (scheme) ! scheme of integration
      case default   ! analytic solution
         fe_temp(i) = (fe(i) - 2.0d0*A/(A+B))*exp(-prefac*(A+B)*dt) + 2.0d0*A/(A+B)
      case (1) ! 1st order, a.k.a. explicit solution for differential equation
         fe_temp(i) = fe(i) + 2.0d0*prefac*((2.0d0 - fe(i))*A - fe(i)*B)*dt
      case (2)       ! 1st order implicit
         fe_temp(i) = (fe(i) + 2.0d0*A*prefac*dt) / (1.0d0 + prefac*(A+B)*dt)
      end select
      !print*, i, prefac, A, B, A/(A+B), fe_temp(i), fe(i)
   else     ! without scattering, nothing changes:
      fe_temp(i) = fe(i)
   endif

   !pause 'Boltzmann_solution'
end subroutine Boltzmann_solution


subroutine mean_distribution(Ev, fe, Ef, fe_f, i_f1, method)
   real(8), dimension(:), intent(in) :: Ev ! [eV] electron energy levels
   real(8), dimension(:), intent(in) :: fe ! electron distribution function
   real(8), intent(in) :: Ef    ! given energy (in between energy levels)
   real(8), intent(out) :: fe_f  ! distribution function at the given point between the energy levels
   integer, intent(out), optional :: i_f1  ! number of energy level the closest one
   integer, intent(in), optional :: method ! what kind of average do we use
   real(8) Fm, sigma, mu, Gaus, three_sigma
   integer i_f, i_cur, N, i_start, i_end
   
   N = size(eV) ! number of energy levels
   sigma = 0.1d0 ! [eV] width of each energy level assigned
   three_sigma = 3.0d0*sigma
   
   ! Find the closest energy level (from above) of the second scattering electron
   call Find_in_array_monoton(Ev, Ef, i_f) ! module "Little_subroutines"
   if (present(i_f1)) i_f1 = i_f ! for the output

   !if (i_f > 1) print*, i_f, Ev(i_f-1), Ev(i_f), Ef
   !pause 'mean_distribution'

   select case (method)
   case (1) ! try Fermi interpolation :: POOR CONSERVATION PROPERTIES
      call Fermi_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
   case (2) ! try "geometric" average :: DOESN'T WORK
      Fm = fe(i_f-1)*(Ef - Ev(i_f-1)) + fe(i_f)*(Ev(i_f) - Ef)
      fe_f = fe(i_f-1)*fe(i_f)*(Ev(i_f) - Ev(i_f-1))/Fm
   case (3) ! do gaussian width of energy levels :: DOESN'T WORK
      fe_f = 0.0d0
      call Find_in_array_monoton(Ev, Ef-three_sigma, i_start) ! module "Little_subroutines"
      call Find_in_array_monoton(Ev, Ef+three_sigma, i_end) ! module "Little_subroutines"
      do i_cur = i_start, i_end ! within three sigma range
         call Gaussian(Ev(i_cur), sigma, Ef, Gaus)
         fe_f = fe_f + Gaus !*(2.0d0*sqrt(2.0d0*g_Pi)*sigma) ! we use this kind of nuormalization to make 2 electrons exctly at a level
         !print*, 'Gaus', i_cur, fe_f, Gaus, (2.0d0*sqrt(2.0d0*g_Pi)*sigma)
      enddo
   case (4) ! try linear interpolation + gaussian width :: DOESN'T WORK
      if (i_f > 1) then
         call linear_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
         if ((Ef - Ev(i_f-1)) > (Ev(i_f) - Ef)) then
            call Gaussian(Ev(i_f), sigma, Ef, Gaus, 2.0d0)
         else 
            call Gaussian(Ev(i_f-1), sigma, Ef, Gaus, 2.0d0)
         endif
         fe_f = fe_f*Gaus
      else
         fe_f = fe(1)
      endif
   case (5) ! linear interpolation with exclusion of too-far-lying levels (Not a very good conservation)
      if (i_f > 1) then
         if ((Ev(i_f) - Ev(i_f-1)) < three_sigma) then
            call linear_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
         else
            fe_f = 0.0d0
            if (present(i_f1)) i_f1 = -1 ! set scattering as impossible
         endif
      else
         fe_f = fe(1)
      endif

   case (6) ! Fermi interpolation with exclusion of too-far-lying levels: POOR CONSERVATION

      if (i_f > 1) then
         if ((Ev(i_f) - Ev(i_f-1)) < three_sigma) then

            call Fermi_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"

         endif ! ((Ev(i_f) - Ev(i_f-1)) < three_sigma)
      endif ! i_f

   case (7) ! Try just the average of the two :: POOR CONSERVATION
      if (i_f > 1) then
         if ((Ev(i_f) - Ev(i_f-1)) < three_sigma) then
            fe_f = 0.5d0*(fe(i_f) + fe(i_f-1))
         else
            fe_f = 0.0d0
            if (present(i_f1)) i_f1 = -1 ! set scattering as impossible
         endif
      else
         fe_f = fe(1)
      endif

   case default ! linear interpolation
      if (i_f > 1) then
         call linear_interpolation(Ev, fe, Ef, fe_f, i_f) ! module "Little_subroutines"
      else
         fe_f = fe(1)
      endif
   endselect

end subroutine mean_distribution


subroutine find_average_level(fe, Ev, i_f1, E_mean, fe_mean)
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   real(8), dimension(:), intent(in) :: fe  ! electron distribution function (ocupation numbers of the energy levels Ei)
   integer, intent(in) :: i_f1 ! number of the level
   real(8), intent(in) :: E_mean ! [eV] exact energy we want to find corresponding occupation to
   real(8), intent(out) :: fe_mean ! average occupation number
   integer :: i_plus
   if (i_f1 == size(Ev)) then ! borderline:
      fe_mean = fe(i_f1)
   else ! within the borders:
      i_plus = i_f1+1
      fe_mean = (fe(i_f1)*(Ev(i_plus) - E_mean) + fe(i_plus)*(E_mean - Ev(i_f1))) / (Ev(i_plus) - Ev(i_f1))
   endif
end subroutine find_average_level



subroutine test_change_of_fe(fe, fe_test)
   real(8), dimension(:), intent(in) :: fe
   real(8), dimension(:), intent(out) :: fe_test
   real(8) :: RN
   integer :: i, N
   N = size(fe)
   do i = 1, N
      call random_number(RN)
      fe_test(i) = fe(i) + 0.001d0*(RN-0.5d0)*fe(i)
   enddo

   where (fe_test(:) > 2.0d0) ! cannot be
      fe_test(:) = 2.0d0
   end where
   where (fe_test(:) < 0.0d0) ! cannot be
      fe_test(:) = 0.0d0
   end where
end subroutine test_change_of_fe


subroutine share_energy(Ev, i, j, Eprime, a, b) 
! Shares energy between levels i and j so that total energy Eprime is conserved
   real(8), dimension(:), intent(in) :: Ev     ! [eV] electron energy levels
   integer, intent(in) :: i, j ! levels to which we give parts of the energy Eprime
   real(8), intent(in) :: Eprime ! [eV] energy that is somewhere in between levels Ev(i) and Ev(j)
   real(8), intent(out) :: a, b  ! coefficients, how much energy must be in the level i and j to conserve total energy Eprime
   real(8) :: Ediff
   Ediff = (Ev(i)-Ev(j))
   a = (Eprime - Ev(j))/Ediff
   b = (Ev(i) - Eprime)/Ediff
!    print*, 'a,b', a, b, Eprime, Ev(i), Ev(j)
!    pause 
end subroutine share_energy



subroutine Electron_electron_scattering_Kij(Ha, Mij, Sij) ! calculates factor in electron-electron scattering matrix element
   real(8), dimension(:,:), intent(in) :: Ha  ! diagonilized Hamiltonian Ha, eigenvectors
   real(8), dimension(:,:), intent(inout) :: Mij  ! Matrix element coupling electron WF via ion motion:
   real(8), dimension(:,:), intent(in), optional :: Sij ! overlap matrix, in case of non-orthogonal basis set
   !----------------------------
   ! For the first and simplest approximation, assume that Mij=Sij:
   if (present(Sij)) then ! non-orthogonal
      Mij = Sij
   else  ! having no idea about the overlap matrix, just use some test-value:
      Mij = 0.01d0
   endif
end subroutine Electron_electron_scattering_Kij




! This subroutine is not ready, it is here only for testing:
subroutine Electron_electron_scattering_Kij_TEST(Ha, Mij, Sij) ! calculates factor in electron-electron scattering matrix element
   real(8), dimension(:,:), intent(in) :: Ha  ! diagonilized Hamiltonian Ha, eigenvectors
   real(8), dimension(:,:), intent(inout) :: Mij  ! Matrix element coupling electron WF via ion motion:
   real(8), dimension(:,:), intent(in), optional :: Sij ! overlap matrix, in case of non-orthogonal basis set
   !----------------------------
   integer i, j, N, k
   real(8), dimension(:), allocatable :: psi0 ! WF_0
   real(8), dimension(size(Ha,1))  :: Norm1   ! normalization factors
   real(8), dimension(size(Ha,1),size(Ha,2)) :: tij
   real(8) :: eps, Norm_val0, Norm_val1

   N = size(Ha,1)

   Mij = 0.0d0 ! to start with
   eps = 1.0d-13    ! acceptable error

!$omp PARALLEL private(i, j, k, psi0, Norm_val1) shared(Norm1, tij)
   allocate(psi0(N))

   ! Ensure WF normalization to 1:
!$omp do
   do i = 1, N
      Norm1(i) = DSQRT(SUM( Ha(:,i) * Ha(:,i) ))
      !print*, i, Norm1(i)
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

            !print*, i, k, tij(i,k)
         enddo ! k = 1, N
      enddo !  i = 1, N
!$omp end do
!$OMP BARRIER
   endif

   !pause 'test 007'

!$omp do
   do i = 1, N ! all energy levels
      if (Norm1(i) < eps) then
         psi0(:) = 0.0d0
      else
         psi0(:) = Ha(:,i) / Norm1(i)
      endif

      !print*, i, psi0(:)
      !pause

      do j = i, N ! only upper triangle, using symmetry
         if (i .NE. j) then

            if (present(Sij)) then ! non-orthogonal
               Mij(i,j) = SUM( psi0(:) * tij(j,:) )
            else ! orthogonal
               if ( Norm1(j) < eps ) then
                  Norm_val1 = 0.0d0
               else
                  Norm_val1 = 1.0d0/Norm1(j)
               endif
               Mij(i,j) = SUM( psi0(:) * Ha(:,j)*Norm_val1 )
            endif ! (present(Sij))

            !print*, i, j, Mij(i,j), Sij(i,j)

            ! Use symmetry:
            Mij(j,i) = Mij(i,j)
         endif
      enddo
   enddo
!$omp end do !!end PARALLEL do

   deallocate(psi0)
!$omp end parallel

   !pause 'Electron_electron_scattering_Kij'
end subroutine Electron_electron_scattering_Kij_TEST




END MODULE Electron_electron_scattering
