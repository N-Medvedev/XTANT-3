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
! This module contains subroutines to deal with Coulomb potential
!
! For Ewalds summation, it follows the example of Allen and Tildesley:
! [1] https://github.com/Allen-Tildesley/examples/blob/master/ewald_module.f90
! Which accompanies the book "Computer Simulation of Liquids", second edition, 2017
! [2] http://micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf
!
! For Wolf's et al. treatement of truncated Coulomb, see:
! [3] D. Wolf, et al., J. Chem. Phys. 110, 8254 (1999); https://doi.org/10.1063/1.478738
! [4] C. J. Fennell and J. D. Gezelter, J. Chem. Phys. 124, 234104 (2006)



MODULE Coulomb
use Universal_constants
use Objects
use Little_subroutines, only : count_3d
use Atomic_tools, only : get_number_of_image_cells, distance_to_given_cell

implicit none
PRIVATE

real(8) :: m_k, m_2Pi2, m_sqrtPi

parameter(m_k = g_ke * g_e * 1.0d10)  ! Constant in the Coulomb law, converting potential into [eV]
parameter(m_2Pi2 = g_2Pi*g_2Pi)     ! 2*Pi^2
parameter(m_sqrtPi = sqrt(g_Pi))    ! sqrt(Pi)

public :: get_Coulomb_Wolf_s, f_cut_L_C, d_f_cut_L_C, m_sqrtPi, Coulomb_Wolf_pot, Coulomb_Wolf_self_term, cut_off_distance, &
m_k, Construct_B_C, get_Coulomb_s, d_Coulomb_Wolf_pot



 contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Wolf's method of Coulomb trunkation:

subroutine get_Coulomb_Wolf_s(Scell, NSC, matter, E_coulomb, gam_ij)   ! Coulomb energy
! This subroutine calcualtes the full Coulomb energy following Wolf's truncation method
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(solid), intent(in), target :: matter   ! material parameters
   real(8), intent(out) :: E_coulomb  ! Total Coulomb energy of all atoms [eV]
   real(8), dimension(:,:), intent(in), optional :: gam_ij  ! effective energy values [eV]
   !=====================================================
   real(8) :: alpha, sum_a, Coul_pot, r_cut, q1, q2
   integer :: j, nat, atom_2, i
   real(8), pointer :: r
   integer, pointer :: m, KOA1, KOA2

   nat = Scell(NSC)%Na ! number of atoms
   ! Get the cut-off distance
   r_cut = cut_off_distance(Scell(NSC)) ! below
   sum_a = 0.0d0  ! to start with

   !$omp PARALLEL private( j, KOA1, q1, m, atom_2, i, KOA2, q2, r, alpha, Coul_pot )
   !$omp do reduction( + : sum_a)
   do j = 1, nat  ! atom #1
      KOA1 => Scell(NSC)%MDatoms(j)%KOA   ! kind of atom #1
      q1 = matter%Atoms(KOA1)%NVB - matter%Atoms(KOA1)%mulliken_Ne ! Mulliken charge of atom #1
      m => Scell(NSC)%Near_neighbor_size(j)  ! number of nearest neighbors of atom #1
      do atom_2 = 1, m  ! only for atoms close to #1
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! atom #2
         KOA2 => Scell(NSC)%MDatoms(i)%KOA   ! kind of atom #2
         q2 = matter%Atoms(KOA2)%NVB - matter%Atoms(KOA2)%mulliken_Ne ! Mulliken charge of atom #2

         r => Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! distance between atoms #1 and #2, R [A]
         alpha = 3.0d0/(4.0d0*r_cut) ! it's chosen according to optimal value from [4]

         Coul_pot = Coulomb_Wolf_pot(q1, q2, r, r_cut, alpha)    ! function below
         if (present(gam_ij)) then  ! Effective value according to Hubbard model:
            ! Renormalize Coulomb potential to a.u., and then to the given Hubbard parameters:
            sum_a = sum_a + Coul_pot * gam_ij(j,i) * g_ev2au
         else  ! Bare Coulomb:
            sum_a = sum_a + Coul_pot
         endif
      enddo ! atom_2
      ! Add the Wolf's self-term:
      sum_a = sum_a + Coulomb_Wolf_self_term(q1, r_cut, alpha) ! below
   enddo ! j
   !$omp end do
   !$omp end parallel
   ! Total Coulomb energy, excluding double-counting:
   E_coulomb = sum_a * 0.5d0   ! [eV]

   nullify(KOA1, KOA2, r)
end subroutine get_Coulomb_Wolf_s



pure function Coulomb_Wolf_pot(q1, q2, r, Rc, alpha) result(WPot) ! truncated Coulomb potential [4]
   real(8) :: WPot  ! [eV]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r, Rc ! [A] interatomic distance, truncation distance
   real(8), intent(in) :: alpha ! truncation parameter
   real(8) :: Rc2, a2, term
   if (r < Rc) then
      Rc2 = Rc*Rc
      a2 = alpha*alpha
      term = erfc(alpha*Rc)/Rc
      WPot = m_k*q1*q2 * (erfc(alpha*r)/r - term + &
         (term/Rc + 2.0d0*alpha/m_sqrtPi*exp(-a2*Rc2)/Rc)*(r-Rc) )  ! [eV]
   else
      WPot = 0.0d0
   endif
end function Coulomb_Wolf_pot


pure function Coulomb_Wolf_self_term(q1, Rc, alpha) result(SelfPot)   ! Self-term, Eq.(5.13) [3]
   real(8) :: SelfPot  ! [eV]
   real(8), intent(in) :: q1    ! [e] charges
   real(8), intent(in) :: Rc    ! [A] interatomic distance, truncation distance
   real(8), intent(in) :: alpha ! truncation parameter
   SelfPot = m_k*q1*q1*(erfc(alpha*Rc)/Rc + 2.0d0*alpha/m_sqrtPi)    ! [eV]
end function Coulomb_Wolf_self_term


pure function d_Coulomb_Wolf_pot(q1, q2, r, Rc, alpha) result(dWPot) ! derivative truncated Coulomb potential, Eq.(5.22) [3]
   real(8) :: dWPot ! [eV/A]
   real(8), intent(in) :: q1, q2    ! [e] charges
   real(8), intent(in) :: r, Rc ! [A] interatomic distance, truncation distance
   real(8), intent(in) :: alpha ! truncation parameter
   real(8) :: r2, Rc2, a2
   if (r < Rc) then
      r2 = r*r
      Rc2 = Rc*Rc
      a2 = alpha*alpha
      dWPot = -m_k*q1*q2*( erfc(alpha*r)/r2 - erfc(alpha*Rc)/Rc2 + &
               2.0d0*alpha/m_sqrtPi*(exp(-a2*r2)/r - exp(-a2*Rc2)/Rc) )  ! [eV/A]
   else
      dWPot = 0.0d0
   endif
end function d_Coulomb_Wolf_pot



function cut_off_distance(Scell) result(d_cut)
   real(8) d_cut  ! [A] cut off distance defined for the given TB Hamiltonian
   type(Super_cell), intent(in) :: Scell  ! supercell with all the atoms as one object
   real(8) :: r(1), Sup_Cell

   ! Cut off of the potential/Hamiltonian:
   ASSOCIATE (ARRAY => Scell%TB_Hamil(:,:)) ! attractive part
      select type(ARRAY)
      type is (TB_H_Pettifor)
         r = maxval(ARRAY(:,:)%rm)
         d_cut = r(1)
      type is (TB_H_Fu)
         r = maxval(ARRAY(:,:)%rm)
         d_cut = r(1)
      type is (TB_H_Molteni)
         r = maxval(ARRAY(:,:)%rcut)
         d_cut = r(1)
      type is (TB_H_NRL)
         r = maxval(ARRAY(:,:)%Rc)
         d_cut = r(1)*g_au2A	! [a.u.] -> [A]
      type is (TB_H_DFTB)
         r = maxval(ARRAY(:,:)%rcut)
         d_cut = r(1)
      type is (TB_H_3TB)
         r = maxval(ARRAY(:,:)%rcut)
         d_cut = r(1)
      type is (TB_H_BOP)
         d_cut = maxval(ARRAY(:,:)%rcut + ARRAY(:,:)%dcut)
      type is (TB_H_xTB)
         r = maxval(ARRAY(:,:)%rcut)
         d_cut = r(1)
      end select
   END ASSOCIATE

   ! Make sure it is not larger than half of the supercell:
   !Sup_Cell = min( sqrt( sum(Scell%supce(1,:)**2) ) , sqrt( sum(Scell%supce(2,:)**2) ), sqrt( sum(Scell%supce(3,:)**2) ) )
   !d_cut = min( d_cut, Sup_Cell*0.5d0 )
end function cut_off_distance




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Softly trunkated Coulomb potential:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_mirror_cell_num_C(Scell, NSC, TB_Coul, atoms, Nx, Ny, Nz)
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Coulomb_cut), dimension(:,:), intent(in) :: TB_Coul ! Coulomb parameters
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   integer, intent(out) :: Nx, Ny, Nz ! number of super-cells images that contribute to total energy
   real(8), dimension(3) :: zb ! super cell indices
   real(8) :: R_cut

   R_cut = TB_Coul(1,1)%dm ! [A] cut-off distance
   call get_number_of_image_cells(Scell, NSC, atoms, R_cut, Nx, Ny, Nz) ! module "Atomic_tools"

   ! It has to be at least one image cell:
   if (Nx < 1) Nx = 1
   if (Ny < 1) Ny = 1
   if (Nz < 1) Nz = 1
end subroutine get_mirror_cell_num_C


subroutine get_Coulomb_s(TB_Coul, Scell, NSC, numpar, a)   ! Coulomb energy
! This subroutine calcualtes the full Coulomb energy
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Coulomb_cut), dimension(:,:), intent(in) :: TB_Coul ! Coulomb parameters
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a	! [eV]
   !=====================================================
   real(8) :: sum_a, a_r, tot_pot, Coul_pot
   integer :: Nx, Ny, Nz, zb(3), N_ind
   INTEGER(4) i1, j1, m, atom_2, x_cell, y_cell, z_cell
   logical :: origin_cell
   
   ! Find how many image cells along each direction we need to include:
   call get_mirror_cell_num_C(Scell, NSC, TB_Coul, Scell(NSC)%MDatoms, Nx, Ny, Nz) ! subroutine above
   
!    call Find_outermost_atom(Scell, 3, 'UP', N_ind)	! module "Atomic_tools"
   tot_pot = 0.0d0
   
   sum_a = 0.0d0
   !$omp PARALLEL private(i1,j1,a_r,x_cell,y_cell,z_cell,zb,origin_cell, Coul_pot)
!   !$omp do reduction( + : sum_a, tot_pot)
   !$omp do reduction( + : sum_a)
   XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
      YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
         ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
            !zb = (/dble(Nx),dble(Ny),dble(Nx)/) ! vector of image of the super-cell
            zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
            origin_cell = ALL(zb==0) ! if it is the origin cell
            do i1 = 1, Scell(NSC)%Na ! all atoms
               do j1 = 1, Scell(NSC)%Na ! all pairs of atoms
                  if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                     call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, a_r) ! module "Atomic_tools"
                     Coul_pot = Coulomb_pot(Scell, NSC, TB_Coul, i1, j1, a_r)    ! function below
                     sum_a = sum_a + Coul_pot
!                      if (i1 == N_ind) tot_pot = tot_pot + Coul_pot
                  endif ! (j1 .NE. i1)
               enddo ! j1
            enddo ! i1
         enddo ZC
      enddo YC
   enddo XC
   !$omp end do
   !$omp end parallel
   a = sum_a	! [eV]
!    write(*,'(a,i4,f,f,f,f)') 'tot_pot', N_ind, Scell(NSC)%MDAtoms(N_ind)%S(3), tot_pot, tot_pot + Coulomb_pot(Scell, NSC, TB_Coul, N_ind, N_ind, 0.25d0), sum_a/dble(Scell(NSC)%Na)
end subroutine get_Coulomb_s


function Coulomb_pot(Scell, NSC, TB_Coul, i1, j1, a_r) ! function according to Girifalco, based on lenard-Jones potential + cut-offs
   !real(8) :: vdW_Girifalco ! van der Waals potential energy [eV]
   real(8) :: Coulomb_pot ! coulomb potential energy [eV]
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Coulomb_cut), dimension(:,:), intent(in), target :: TB_Coul ! Coulomb parameters with TB
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   integer, intent(in) :: i1, j1 ! atoms
   !====================================
   real(8), pointer :: k	! 
   real(8), pointer :: dm 	! [A] radius where to switch to polinomial
   real(8), pointer :: dd 	! [A] cut-off radius at large r
   real(8), pointer :: Q	! mean charge
   integer, pointer :: KOAi, KOAj
   !====================================
   real(8) :: f_cut_large
   
   KOAi => Scell(NSC)%MDatoms(i1)%KOA ! kind of first atom
   KOAj => Scell(NSC)%MDatoms(j1)%KOA ! kind of second atom
   Q => Scell(NSC)%Q
   k => TB_Coul(KOAi, KOAj)%k
   dd => TB_Coul(KOAi, KOAj)%dd
   dm => TB_Coul(KOAi, KOAj)%dm

   if (a_r > dm+dd*10.0d0) then ! anything beyond cut-offs is zero:
      Coulomb_pot = 0.0d0
   else ! at shorter distances we use proper potential:
      f_cut_large = f_cut_L_C(a_r, dm, dd) ! function below
      Coulomb_pot = Coulomb_potential(k, Q, Q, a_r)*f_cut_large
   endif
!    print*, 'Coulomb_pot:', k, Q, a_r, Coulomb_pot
!    pause 'test'
   nullify(KOAi, KOAj, k, dd, dm, Q)
end function Coulomb_pot



pure function Coulomb_potential(k, Q1, Q2, r) ! Coulomb potential
   real(8) :: Coulomb_potential
   real(8), intent(in) :: k, Q1, Q2, r
   Coulomb_potential = k*Q1*Q2/r ! [eV]
end function Coulomb_potential


pure function d_Coulomb_potential(k, Q1, Q2, r) ! derivative of Coulomb potential
   real(8) :: d_Coulomb_potential
   real(8), intent(in) :: k, Q1, Q2, r
   d_Coulomb_potential = -k*Q1*Q2/(r*r)
end function d_Coulomb_potential


pure function f_cut_L_C(a_r, r_L, d_L) ! cut-off function at large distances
   real(8) :: f_cut_L_C
   real(8), intent(in) :: a_r, r_L, d_L
   f_cut_L_C = 1.0d0/(1.0d0 + dexp((a_r - r_L)/d_L))
end function f_cut_L_C


pure function d_f_cut_L_C(a_r, r_L, d_L) ! derivative of cut-off function at large distances
   real(8) :: d_f_cut_L_C
   real(8), intent(in) :: a_r, r_L, d_L
   real(8) :: exp_r, exp_r2, arg
   arg = (a_r - r_L)/d_L
   if (arg >= log(HUGE(a_r))) then
      d_f_cut_L_C = 0.0d0
   else
      exp_r = dexp(arg)
      exp_r2 = 1.0d0 + exp_r
      d_f_cut_L_C = -exp_r/(d_L*exp_r2*exp_r2)
   endif
end function d_f_cut_L_C



function dCoulomb(Scell, TB_Coul, a_r)
   real(8) :: dCoulomb ! derivative of the trunkated Coulomb
   type(Super_cell), intent(in), target :: Scell  ! supercell as one object
   type(TB_Coulomb_cut), intent(in), target :: TB_Coul   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_r ! [A] distance between the atoms
   !====================================
   real(8), pointer :: k	! 
   real(8), pointer :: dm 	! [A] radius where to switch to polinomial
   real(8), pointer :: dd 	! [A] cut-off radius at large r
   real(8), pointer :: Q	! mean charge
   !====================================
   real(8) :: E_C, f_cut_large, d_f_large, d_E_C

   Q => Scell%Q
   k => TB_Coul%k
   dd => TB_Coul%dd
   dm => TB_Coul%dm

   if (a_r > dm+dd*10.0d0) then ! anything beyond cut-offs is zero:
      dCoulomb = 0.0d0
   else ! at shorter distances we use proper potential:
      f_cut_large = f_cut_L_C(a_r, dm, dd) ! function above
      E_C = Coulomb_potential(k, Q, Q, a_r) ! function above
      d_E_C = d_Coulomb_potential(k, Q, Q, a_r) ! derivative of  Coulomb part of the potential
      d_f_large = d_f_cut_L_C(a_r, dm, dd) ! derivative of cut-off function, function above
      dCoulomb = d_E_C*f_cut_large + E_C*d_f_large
!       print*, 'Coulomb_pot:', k, Q, a_r, dCoulomb, d_E_C, f_cut_large, E_C, d_f_large
   endif
!    pause 'test'
   nullify(k, dd, dm, Q)
end function dCoulomb



! Construct matrices of multipliers often used later to calculate derivatives of vdW:
subroutine Construct_B_C(TB_Coul, Scell, NSC, atoms, Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz)
   type(TB_Coulomb_cut), dimension(:,:), intent(in)   :: TB_Coul ! van der Waals parameters within TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   real(8), dimension(:,:,:), allocatable, intent(out) :: Bij, A_rij, XijSupce, YijSupce, ZijSupce, Xij, Yij, Zij, SXij, SYij, SZij ! intermediate calculations matrices
   integer, intent(out) :: Nx, Ny, Nz ! number of super-cell images to consider
   real(8) :: x, y, z, sx, sy, sz
   integer :: x_cell, y_cell, z_cell, coun_cell
   integer :: i1, j1, n
   integer, DIMENSION(3) :: zb
   logical :: origin_cell
   
   n = size(atoms) ! total number of atoms
  
  ! Find how many image cells along each direction we need to include:
   call get_mirror_cell_num_C(Scell, NSC, TB_Coul, atoms, Nx, Ny, Nz) ! subroutine above

   coun_cell = (2*Nx+1)*(2*Ny+1)*(2*Nz+1) ! total number of image cells

   if (.not.allocated(Bij))   allocate(Bij(coun_cell,n,n))
   if (.not.allocated(A_rij)) allocate(A_rij(coun_cell,n,n))
   if (.not.allocated(Xij))   allocate(Xij(coun_cell,n,n))
   if (.not.allocated(Yij))   allocate(Yij(coun_cell,n,n))
   if (.not.allocated(Zij))   allocate(Zij(coun_cell,n,n))
   if (.not.allocated(SXij))  allocate(SXij(coun_cell,n,n))
   if (.not.allocated(SYij))  allocate(SYij(coun_cell,n,n))
   if (.not.allocated(SZij))  allocate(SZij(coun_cell,n,n))
   if (.not.allocated(XijSupce))   allocate(XijSupce(coun_cell,n,n))
   if (.not.allocated(YijSupce))   allocate(YijSupce(coun_cell,n,n))
   if (.not.allocated(ZijSupce))   allocate(ZijSupce(coun_cell,n,n))
   Bij  = 0.0d0
   A_rij = 1.0d30
   Xij  = 0.0d0
   Yij  = 0.0d0
   Zij  = 0.0d0
   SXij = 0.0d0
   SYij = 0.0d0
   SZij = 0.0d0
   XijSupce = 0.0d0
   YijSupce = 0.0d0
   ZijSupce = 0.0d0
   
   CHARGE:if (Scell(NSC)%Q <= 1.0d-12) then ! no charge, no force
      ! do nothing then, it's all zero
   else CHARGE ! there is charge, get Coulomb
      !$omp PARALLEL private(i1,j1,x,y,z,sx,sy,sz,x_cell,y_cell,z_cell,zb,origin_cell,coun_cell)
      !$omp do
      do i1 = 1,n
         do j1 = 1,n
            XC:do x_cell = -Nx, Nx ! all images of the super-cell along X
               YC:do y_cell = -Ny, Ny ! all images of the super-cell along Y
                  ZC:do z_cell = -Nz, Nz ! all images of the super-cell along Z
                     coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"
                     zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
                     origin_cell = ALL(zb==0) ! if it is the origin cell
                     if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                        call distance_to_given_cell(Scell, NSC, Scell(NSC)%MDatoms, dble(zb), i1, j1, A_rij(coun_cell,i1,j1), x=x, y=y, z=z, sx=sx,sy=sy,sz=sz) ! module "Atomic_tools"
                        !call shortest_distance(Scell, NSC, atoms, i1, j1, A_rij(coun_cell,i1,j1), x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz) ! module "Atomic_tools"

                        ! Save necessary elements:
                        Xij(coun_cell,i1,j1) = x
                        Yij(coun_cell,i1,j1) = y
                        Zij(coun_cell,i1,j1) = z
                        SXij(coun_cell,i1,j1) = sx
                        SYij(coun_cell,i1,j1) = sy
                        SZij(coun_cell,i1,j1) = sz

                        XijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
                        YijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
                        ZijSupce(coun_cell,i1,j1) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)
                        !Bij(coun_cell,i1,j1) = dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1)) ! function above
                        Bij(coun_cell,i1,j1) = dCoulomb(Scell(NSC), TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA), A_rij(coun_cell,i1,j1)) ! function above

!                         print*, 'TEST', coun_cell,i1,j1, Xij(coun_cell,i1,j1), Yij(coun_cell,i1,j1), Zij(coun_cell,i1,j1), A_rij(coun_cell,i1,j1), Bij(coun_cell,i1,j1)
                     else ! No self-interaction
                        A_rij(coun_cell,i1,j1) = 1.0d30
                        Xij(coun_cell,i1,j1) = 0.0d0
                        Yij(coun_cell,i1,j1) = 0.0d0
                        Zij(coun_cell,i1,j1) = 0.0d0
                        Bij(coun_cell,i1,j1) = 0.0d0
                     endif
                  enddo ZC
               enddo YC
            enddo XC
         enddo ! j1
      enddo ! i1
      !$omp end do
      !$omp end parallel
   endif CHARGE
!    pause 'Construct_B'
end subroutine Construct_B_C



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unused subroutines:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Derivatives of the coulomb energy by s:
subroutine dCoulombdr_s(atoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, Nx, Ny, Nz) 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: Bij, A_rij, Xij, Yij, Zij ! intermediate calculations matrices
   integer, intent(in) :: Nx, Ny, Nz ! number of super-cells to go through
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   REAL(8), DIMENSION(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, x, y, z, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, j1, ian, dik, djk, n, m, atom_2
   integer :: x_cell, y_cell, z_cell, coun_cell
   integer, DIMENSION(3) :: zb
   real(8), dimension(:,:), allocatable :: Erx_s
   logical :: origin_cell
   n = size(atoms) ! total number of atoms
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

   !$omp PARALLEL private(i1,j1,ian,dik,djk,dpsi,x1,b,a_r,ddlta,b_delta, x_cell, y_cell, z_cell, coun_cell, zb, origin_cell)
   !$omp DO
   do ian = 1, n  ! Forces for all atoms
      !Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with   
      do i1 = 1, n ! contribution from all atoms
         if (ian == i1) then
            dik = 1
         else
            dik = 0
         endif
         dpsi = 0.0d0
         do j1 = 1,n ! for each pair of atoms
            if (ian == j1) then
               djk = 1
            else
               djk = 0
            endif
            cos_if:if ((dik-djk) /= 0) then ! without it, it gives ERROR
               XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
                  YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
                     ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
                        zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
                        origin_cell = ALL(zb==0) ! if it is the origin cell
                        cell_if:if ((j1 /= i1) .or. (.not.origin_cell)) then ! exclude self-interaction only within original super cell
                           ! contribution from this image of the cell:
                           coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"
                     
                           x1(1) = Xij(coun_cell,i1,j1) ! x*supce(1,1) + y*supce(1,2) + z*supce(1,3)
                           x1(2) = Yij(coun_cell,i1,j1) ! x*supce(2,1) + y*supce(2,2) + z*supce(2,3)
                           x1(3) = Zij(coun_cell,i1,j1) ! x*supce(3,1) + y*supce(3,2) + z*supce(3,3)
                           b = Bij(coun_cell,i1,j1) ! dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))
                           a_r = A_rij(coun_cell,i1,j1) ! call shortest_distance(Scell, NSC, atoms, i1, j1, A_rij(coun_cell,i1,j1), x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz) 

                           ddlta = real(dik - djk)/a_r
                           b_delta = b*ddlta
                           dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                           dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                           dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
                        endif cell_if
                     enddo ZC2
                  enddo YC2
               enddo XC2
            endif cos_if
         enddo ! j1

         Erx_s(1,ian) = Erx_s(1,ian) + dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      ! Add van der Waals force to already calculated other forces:
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Scell(NSC)%MDatoms(ian)%forces%rep(:) + Erx_s(:,ian)
!       print*, 'NEW',ian, Erx_s(:,ian)
   enddo ! ian
   !$omp end do
   !$omp end parallel

   deallocate(Erx_s)
!    pause 'dCoulombdr_s'
END subroutine dCoulombdr_s



! Derivatives of the coulomb energy by h:
subroutine dCoulombdr_Pressure_s(atoms, Scell, NSC, numpar, Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij, Nx, Ny, Nz) 
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   real(8), dimension(:,:,:), intent(in) :: Bij, A_rij, Xij, Yij, Zij, SXij, SYij, SZij ! intermediate calculations matrices
   integer, intent(in) :: Nx, Ny, Nz ! number of super-cells to go through
   !===============================================
   REAL(8), DIMENSION(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i,j,k,l,n
   integer x_cell, y_cell, z_cell, zb(3), coun_cell
   real(8) r, rcur(3), scur(3), PForce(3,3)
   real(8) df_psy, psi, dpsy 
   logical :: origin_cell

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = size(atoms)

      PForce = 0.0d0 ! to start with
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         do j = 1, n ! do for all pairs of atoms
            XC2:do x_cell = -Nx, Nx ! all images of the super-cell along X
               YC2:do y_cell = -Ny, Ny ! all images of the super-cell along Y
                  ZC2:do z_cell = -Nz, Nz ! all images of the super-cell along Z
                     zb = (/x_cell,y_cell,z_cell/) ! vector of image of the super-cell
                     origin_cell = ALL(zb==0) ! if it is the origin cell
                     cell_if:if ((j /= i) .or. (.not.origin_cell)) then ! exclude self-interaction within original super cell
                        ! contribution from this image of the cell:
                        coun_cell = count_3d(Nx,Ny,Nz,x_cell,y_cell,z_cell) ! module "Little_sobroutine"

                        rcur(1) = Xij(coun_cell,i,j) ! x
                        rcur(2) = Yij(coun_cell,i,j) ! y
                        rcur(3) = Zij(coun_cell,i,j) ! z
                        scur(1) = SXij(coun_cell,i,j) ! SX
                        scur(2) = SYij(coun_cell,i,j) ! SY
                        scur(3) = SZij(coun_cell,i,j) ! SZ
                        r = A_rij(coun_cell,i,j)     ! R
                        dpsy = Bij(coun_cell,i,j)    ! dvdW(TB_Coul(Scell(NSC)%MDatoms(j1)%KOA,Scell(NSC)%MDatoms(i1)%KOA),A_rij(coun_cell,i1,j1))

                        do k = 1,3 ! supce indices: a,b,c
                           do l = 1,3  ! supce indices: x,y,z
                              Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                           enddo ! l
                        enddo ! k
                     endif cell_if
                  enddo ZC2
               enddo YC2
            enddo XC2
         enddo ! j

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               !Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + Rep_Pr(l,k) !*0.5d0
               PForce(l,k) = PForce(l,k) + Rep_Pr(l,k)
            enddo ! l
         enddo ! k
      enddo ! i
      Scell(NSC)%SCforce%rep = Scell(NSC)%SCforce%rep + PForce ! add vdW part to existing TB part
!       print*, 'NEW P', PForce
   endif
!    pause 'dvWdr_Pressure_s'
end subroutine dCoulombdr_Pressure_s


END MODULE Coulomb
