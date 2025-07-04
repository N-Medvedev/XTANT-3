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
! This module contains subroutines to deal with ZBL potential
! https://en.wikipedia.org/wiki/Stopping_power_(particle_radiation)#Repulsive_interatomic_potentials
!
! And its improved version, the NLH potential:
! https://doi.org/10.1103/PhysRevA.111.032818



MODULE ZBL_potential
use Universal_constants
use Objects
use Atomic_tools, only : shortest_distance
#ifdef MPI_USED
use MPI_subroutines, only : do_MPI_Reduce
#endif

implicit none
PRIVATE

! Modular parameters:
real(8), parameter :: m_a_u = 0.8854d0
real(8), parameter :: m_phi1 = 0.1818d0
real(8), parameter :: m_phi2 = 0.5099d0
real(8), parameter :: m_phi3 = 0.2802d0
real(8), parameter :: m_phi4 = 0.02817d0
real(8), parameter :: m_exp1 = -3.2d0
real(8), parameter :: m_exp2 = -0.9423d0
real(8), parameter :: m_exp3 = -0.4028d0
real(8), parameter :: m_exp4 = -0.2016d0
real(8), parameter :: m_k = 1.0d0/(4.0d0 * g_Pi * g_e0)


character(30), parameter :: m_INPUT_NLH = 'nlh_coeffs.dat'

!==============================================
! For reading the NLH potential coefficients:
type NLH_coeffs    ! file with the coefficients, "nlh_coeffs.dat"
   integer :: Z1, Z2     ! atomic number Z1 and Z2 of the pair of elements
   real(8), dimension(3) :: a ! miltipliers
   real(8), dimension(3) :: b ! exponential coefficients
endtype NLH_coeffs
!==============================================


public :: ZBL_pot, d_ZBL_pot, get_total_ZBL, d2_ZBL_pot, get_NLH_coefficients, get_total_NLH, NLH_pot, d_NLH_pot, d2_NLH_pot


 contains



!======================================================
! NLH potential:


subroutine get_NLH_coefficients(Z1, Z2, a, b, pathname, error_message, filename_in)
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), dimension(3), intent(out) :: a, b    ! exponential coefficients: [-], [1/A]
   character(*), intent(in) :: pathname ! where to find the file with the coefficients
   character(*), intent(inout) :: error_message ! error message
   character(*), intent(in), optional :: filename_in     ! file with the coefficients, if differs from default
   !---------------------------------
   integer :: FN, Reason, i, j
   character(500) :: Filename
   logical :: file_exists, file_opened
   type(NLH_coeffs), dimension(92,92) :: NLH_data

   error_message = ''   ! no error to start with
   a(:) = 0.0d0 ! to start with
   b(:) = 0.0d0 ! to start with


   ! Find the file with NLH coefficients:
   if (present(filename_in)) then
      Filename = trim(adjustl(pathname))//trim(adjustl(filename_in))
   else ! default
      Filename = trim(adjustl(pathname))//trim(adjustl(m_INPUT_NLH))
   endif

   inquire(file=trim(adjustl(Filename)),exist=file_exists) ! check if input file is there
   exists:if (file_exists) then
      ! Read all the coefficients from this file:
      FN = 101
      open (unit=FN, file=trim(adjustl(Filename)))
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) then ! read data from it
         read(FN,*,IOSTAT=Reason) ! skip first line with comments
         if (Reason /= 0) then
            write(error_message,*) 'Could not read input file: ', trim(adjustl(Filename))
            return ! nothing else to do
         endif

         do i = 1, 92 ! number of lines in the file defined by number of elements included in it
            do j = 1, 92
               if (i <= j) then
                  read(FN,*,IOSTAT=Reason) NLH_data(i,j)%Z1, NLH_data(i,j)%Z2, NLH_data(i,j)%a(1), NLH_data(i,j)%b(1), &
                                                                               NLH_data(i,j)%a(2), NLH_data(i,j)%b(2), &
                                                                               NLH_data(i,j)%a(3), NLH_data(i,j)%b(3)
                  if (Reason /= 0) then
                     write(error_message,*) 'Could not read input file: ', trim(adjustl(Filename))
                     return ! nothing else to do
                  endif
               else ! Use symmetry to fill in the lower triangle:
                  NLH_data(i,j)%Z1 = NLH_data(j,i)%Z1
                  NLH_data(i,j)%Z2 = NLH_data(j,i)%Z2
                  NLH_data(i,j)%a(:) = NLH_data(j,i)%a(:)
                  NLH_data(i,j)%b(:) = NLH_data(j,i)%a(:)
               endif
            enddo ! j
         enddo ! i
         close(FN)
      else
         write(error_message,*) 'Could not open the file: ', trim(adjustl(Filename))
         write(*,'(a)') trim(adjustl(error_message))
      endif
   else exists ! not
      write(error_message,*) 'Could not find the file: ', trim(adjustl(Filename))
      write(*,'(a)') trim(adjustl(error_message))
   endif exists

   ! Get out the data for the required elements:
   a(:) = NLH_data(Z1,Z2)%a(:)
   b(:) = NLH_data(Z1,Z2)%b(:)
end subroutine get_NLH_coefficients


pure function NLH_pot(Z1, Z2, r, a, b) result(V_NLH)
   real(8) V_NLH    ! NLH potential [eV]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8), dimension(3), intent(in) :: a, b    ! exponential coefficients: [-], [1/A]
   !----------
   real(8) :: phi

   phi = SUM(a(:) * exp(-b(:) * r ))

   V_NLH = m_k * Z1 * Z2 * g_e/(1.0d-10*r) * phi    ! [eV] NLH potential
end function NLH_pot


pure function d_NLH_pot(Z1, Z2, r, a, b) result(d_V_NLH)
   real(8) d_V_NLH  ! derivative of NLH potential [eV/A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8), dimension(3), intent(in) :: a, b    ! exponential coefficients: [-], [1/A]
   !--------------------
   real(8) phi, d_phi, V_NLH_part

   V_NLH_part = m_k * Z1 * Z2 * g_e/(1.0d-10*r) ! part without phi [eV]

   phi      = SUM(a(:) * exp(-b(:) * r ))
   d_phi    = SUM(-b(:) * a(:) * exp(-b(:) * r ))

   d_V_NLH  = V_NLH_part * (d_phi - phi/r)  ! [eV/A] derivative of NLH potential
end function d_NLH_pot


pure function d2_NLH_pot(Z1, Z2, r, a, b) result(d2_V_NLH)
   real(8) d2_V_NLH  ! derivative of NLH potential [eV/A^2]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8), dimension(3), intent(in) :: a, b    ! exponential coefficients: [-], [1/A]
   !---------------------
   real(8) phi, d_phi, V_NLH_part, d2_phi, d_V_NLH

   V_NLH_part = m_k * Z1 * Z2 * g_e/(1.0d-10*r) ! part without phi [eV]

   phi      = SUM(a(:) * exp(-b(:) * r ))
   d_phi    = SUM(-b(:) * a(:) * exp(-b(:) * r ))
   d2_phi    = SUM(b(:)**2 * a(:) * exp(-b(:) * r ))

   d_V_NLH  = V_NLH_part * (d_phi - phi/r)  ! [eV/A] derivative of NLH potential

   d2_V_NLH  = d_V_NLH*(d_phi - phi/r) + V_NLH_part*(d2_phi - d_phi/r + phi/r**2)     ! [eV/A^2] second derivative of NLH potential
end function d2_NLH_pot




subroutine get_total_NLH(Scell, NSC, matter, numpar, a)   ! NLH energy
! This subroutine is only used for comparison of the interlayer NLH energy with other works
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(solid), intent(in) :: matter   ! material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a  ! total NLH repulsive energy [eV]
   !=====================================================
   real(8) :: sum_a, a_r, Z1, Z2, NLH_a(3), NLH_b(3)
   INTEGER(4) i1, j1, m, atom_2
   integer :: Nstart, Nend, N_incr
   integer, pointer :: KOA1, KOA2

   sum_a = 0.0d0

#ifdef MPI_USED   ! only does anything if the code is compiled with MPI

   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank
   Nend = Scell(NSC)%Na

   ! Do the cycle (parallel) calculations:
   do i1 = Nstart, Nend, N_incr  ! each process does its own part
      m = Scell(NSC)%Near_neighbor_size(i1)
      KOA1 => Scell(NSC)%MDatoms(i1)%KOA   ! kind of atom #1
      Z1 = matter%Atoms(KOA1)%Z  ! Z of element #1
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         KOA2 => Scell(NSC)%MDatoms(j1)%KOA   ! kind of atom #2
         Z2 = matter%Atoms(KOA2)%Z  ! Z of element #2
         if ( j1 /= i1 ) then ! count only interplane energy:
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"

            ASSOCIATE (ARRAY => Scell(1)%TB_Expwall(:,:)) ! this is the sintax we have to use to check the class of defined types
               select type (ARRAY)
               type is (TB_Short_Rep)
                  NLH_a(:) = ARRAY(KOA1,KOA2)%f_NLH%a(:)
                  NLH_b(:) = ARRAY(KOA1,KOA2)%f_NLH%b(:)
               endselect
            END ASSOCIATE
            sum_a = sum_a + NLH_pot(Z1, Z2, a_r, NLH_a, NLH_b)    ! function above
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1

   !-----------
   ! Collect information from all processes into the master process:
   ! https://rookiehpc.org/mpi/docs/mpi_reduce/index.html
   call do_MPI_Reduce(numpar%MPI_param, 'Error in get_total_NLH:', sum_a) ! module "MPI_subroutines"

   a = sum_a * 0.5d0

#else ! OpenMP is used instead
   !$omp PARALLEL private(i1,j1,m,KOA1,Z1,atom_2,KOA2,Z2,a_r, NLH_a, NLH_b)
   !$omp do reduction( + : sum_a)
   do i1 = 1, Scell(NSC)%Na ! all atoms
      m = Scell(NSC)%Near_neighbor_size(i1)
      KOA1 => Scell(NSC)%MDatoms(i1)%KOA   ! kind of atom #1
      Z1 = matter%Atoms(KOA1)%Z  ! Z of element #1
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         KOA2 => Scell(NSC)%MDatoms(j1)%KOA   ! kind of atom #2
         Z2 = matter%Atoms(KOA2)%Z  ! Z of element #2
         if ( j1 /= i1 ) then ! count only interplane energy:
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"

            ASSOCIATE (ARRAY => Scell(1)%TB_Expwall(:,:)) ! this is the sintax we have to use to check the class of defined types
               select type (ARRAY)
               type is (TB_Short_Rep)
                  NLH_a(:) = ARRAY(KOA1,KOA2)%f_NLH%a(:)
                  NLH_b(:) = ARRAY(KOA1,KOA2)%f_NLH%b(:)
               endselect
            END ASSOCIATE

            sum_a = sum_a + NLH_pot(Z1, Z2, a_r, NLH_a, NLH_b)    ! function above
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1
   !$omp end do
   nullify(KOA1, KOA2)
   !$omp end parallel
   a = sum_a * 0.5d0
#endif
   nullify(KOA1, KOA2)
end subroutine get_total_NLH



!======================================================
! ZBL potential:


subroutine get_total_ZBL(Scell, NSC, matter, numpar, a)   ! ZBL energy
! This subroutine is only used for comparison of the interlayer ZBL energy with other works
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(solid), intent(in) :: matter   ! material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a  ! total ZBL repulsive energy [eV]
   !=====================================================
   real(8) :: sum_a, a_r, Z1, Z2
   INTEGER(4) i1, j1, m, atom_2
   integer :: Nstart, Nend, N_incr
   integer, pointer :: KOA1, KOA2

   sum_a = 0.0d0

#ifdef MPI_USED   ! only does anything if the code is compiled with MPI

   N_incr = numpar%MPI_param%size_of_cluster    ! increment in the loop
   Nstart = 1 + numpar%MPI_param%process_rank
   Nend = Scell(NSC)%Na

   ! Do the cycle (parallel) calculations:
   !do i1 = 1, Scell(NSC)%Na, N_incr ! all atoms
   do i1 = Nstart, Nend, N_incr  ! each process does its own part
      m = Scell(NSC)%Near_neighbor_size(i1)
      KOA1 => Scell(NSC)%MDatoms(i1)%KOA   ! kind of atom #1
      Z1 = matter%Atoms(KOA1)%Z  ! Z of element #1
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         KOA2 => Scell(NSC)%MDatoms(j1)%KOA   ! kind of atom #2
         Z2 = matter%Atoms(KOA2)%Z  ! Z of element #2
         if ( j1 /= i1 ) then ! count only interplane energy:
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
            sum_a = sum_a + ZBL_pot(Z1, Z2, a_r)    ! function below
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1

   !-----------
   ! Collect information from all processes into the master process:
   ! https://rookiehpc.org/mpi/docs/mpi_reduce/index.html
   call do_MPI_Reduce(numpar%MPI_param, 'Error in get_total_ZBL:', sum_a) ! module "MPI_subroutines"

   a = sum_a * 0.5d0

#else ! OpenMP is used instead
   !$omp PARALLEL private(i1,j1,m,KOA1,Z1,atom_2,KOA2,Z2,a_r)
   !$omp do reduction( + : sum_a)
   do i1 = 1, Scell(NSC)%Na ! all atoms
      m = Scell(NSC)%Near_neighbor_size(i1)
      KOA1 => Scell(NSC)%MDatoms(i1)%KOA   ! kind of atom #1
      Z1 = matter%Atoms(KOA1)%Z  ! Z of element #1
      do atom_2 = 1, m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         KOA2 => Scell(NSC)%MDatoms(j1)%KOA   ! kind of atom #2
         Z2 = matter%Atoms(KOA2)%Z  ! Z of element #2
         if ( j1 /= i1 ) then ! count only interplane energy:
            call shortest_distance(Scell, NSC, Scell(NSC)%MDatoms, i1, j1, a_r) ! module "Atomic_tools"
            sum_a = sum_a + ZBL_pot(Z1, Z2, a_r)    ! function below
         endif ! (j1 .NE. i1)
      enddo ! j1
   enddo ! i1
   !$omp end do
   nullify(KOA1, KOA2)
   !$omp end parallel
   a = sum_a * 0.5d0
#endif
   nullify(KOA1, KOA2)
end subroutine get_total_ZBL


pure function ZBL_pot(Z1, Z2, r) result(V_ZBL)
   real(8) V_ZBL    ! ZBL potential [eV]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) phi
   phi = ZBL_phi(Z1, Z2, r) ! see below
   V_ZBL = m_k * Z1 * Z2 * g_e/(1.0d-10*r) * phi    ! [eV] ZBL potential
end function ZBL_pot


pure function d_ZBL_pot(Z1, Z2, r) result(d_V_ZBL)
   real(8) d_V_ZBL  ! derivative of ZBL potential [eV/A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) phi, d_phi, V_ZBL_part
   
   V_ZBL_part = m_k * Z1 * Z2 * g_e/(1.0d-10*r) ! part without phi [eV]
   phi      = ZBL_phi(Z1, Z2, r) ! see below
   d_phi    = d_ZBL_phi(Z1, Z2, r) ! see below
   d_V_ZBL  = V_ZBL_part * (d_phi - phi/r)  ! [eV/A] derivative of ZBL potential
end function d_ZBL_pot


pure function d2_ZBL_pot(Z1, Z2, r) result(d2_V_ZBL)
   real(8) d2_V_ZBL  ! derivative of ZBL potential [eV/A^2]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) phi, d_phi, V_ZBL_part, d2_phi, d_V_ZBL

   V_ZBL_part = m_k * Z1 * Z2 * g_e/(1.0d-10*r) ! part without phi [eV]
   phi      = ZBL_phi(Z1, Z2, r) ! see below
   d_phi    = d_ZBL_phi(Z1, Z2, r) ! see below
   d_V_ZBL  = V_ZBL_part * (d_phi - phi/r)  ! [eV/A] derivative of ZBL potential
   d2_phi   = d2_ZBL_phi(Z1, Z2, r) ! see below

   d2_V_ZBL  = d_V_ZBL*(d_phi - phi/r) + V_ZBL_part*(d2_phi - d_phi/r + phi/r**2)     ! [eV/A^2] second derivative of ZBL potential
end function d2_ZBL_pot



pure function d_ZBL_phi(Z1, Z2, r) result(d_phi)
   real(8) d_phi ! [1/A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) :: a, x
   a = ZBL_a(Z1, Z2)    ! see below
   x = r / a    ! normalized distance
   d_phi = (m_phi1 * m_exp1 * exp(m_exp1 * x) + &
            m_phi2 * m_exp2 * exp(m_exp2 * x) + &
            m_phi3 * m_exp3 * exp(m_exp3 * x) + &
            m_phi4 * m_exp4 * exp(m_exp4 * x)) / a
end function d_ZBL_phi


pure function d2_ZBL_phi(Z1, Z2, r) result(d_phi)
   real(8) d_phi ! [1/A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) :: a, x
   a = ZBL_a(Z1, Z2)    ! see below
   x = r / a    ! normalized distance
   d_phi = (m_phi1 * m_exp1**2 * exp(m_exp1 * x) + &
            m_phi2 * m_exp2**2 * exp(m_exp2 * x) + &
            m_phi3 * m_exp3**2 * exp(m_exp3 * x) + &
            m_phi4 * m_exp4**2 * exp(m_exp4 * x)) / a**2
end function d2_ZBL_phi


pure function ZBL_phi(Z1, Z2, r) result(phi)
   real(8) phi
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   real(8), intent(in) :: r ! interatomic distance [A]
   real(8) :: a, x
   a = ZBL_a(Z1, Z2)    ! see below
   x = r / a    ! normalized distance
   phi = m_phi1 * exp(m_exp1 * x) + &
         m_phi2 * exp(m_exp2 * x) + &
         m_phi3 * exp(m_exp3 * x) + &
         m_phi4 * exp(m_exp4 * x)
end function ZBL_phi


pure function ZBL_a(Z1, Z2) result (a)
   real(8) a    ! [A]
   real(8), intent(in) :: Z1, Z2    ! atomic nuclear charges (atomic numbers)
   a = m_a_u * g_a0/(Z1**0.23d0 + Z2**0.23d0)
end function ZBL_a




END MODULE ZBL_potential
