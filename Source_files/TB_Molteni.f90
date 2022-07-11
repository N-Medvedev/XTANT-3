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
! This module contains subroutines to deal with TB hamiltonian in the Pettifor parametrization

MODULE TB_Molteni
use Universal_constants
use Objects
! use Variables
use Algebra_tools
use Little_subroutines
use Atomic_tools
use Electron_tools
use Nonadiabatic

implicit none

 contains

subroutine construct_TB_H_Molteni(numpar, matter, TB_Hamil, Scell, NSC, Ha, Err) ! module "TB_Molteni"
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(solid), intent(inout) :: matter	! materil parameters
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:), intent(out) :: Ha ! Constructed Hamiltonian
   type(Error_handling), intent(inout) :: Err	! error save
   character(200) :: Error_descript
   integer i, j
!    complex, dimension(:,:), allocatable :: Ha_c
!    complex, dimension(:), allocatable :: Ei_c

   Error_descript = ''

   ! Get lists of nearest neighbors:
!    call get_near_neighbours(Scell) 	! see "Atomic_tools"
   Scell(NSC)%Ha0 = Scell(NSC)%Ha 	! save Hamiltonial from previous time-step
   Scell(NSC)%H_non0 = Scell(NSC)%H_non	! save non-diagonalized Hamiltonian from last time-step
   Scell(NSC)%Ei0 = Scell(NSC)%Ei ! save energy levels for the next timestep
   
   ! Construct TB Hamiltonian (with Molteni parameters):
   call Hamil_tot_M(numpar, Scell, NSC, TB_Hamil, Ha) ! see below
   Scell(NSC)%H_non = Ha ! save non-diagonalized Hamiltonian
   !call check_symmetry(Ha) ! just for testing, must be hermitian Hamiltonian

   ! Diagonalize symmetric Hamiltonian to get electron energy levels:
   call sym_diagonalize(Ha, Scell(NSC)%Ei, Error_descript, check_M=.true.)

   if (LEN(trim(adjustl(Error_descript))) .GT. 0) then
      Error_descript = 'Subroutine construct_TB_H_Molteni: '//trim(adjustl(Error_descript))
      call Save_error_details(Err, 6, Error_descript)
      print*, trim(adjustl(Error_descript))
   endif

   ! Get band gap:
   call find_band_gap(Scell(NSC)%Ei, Scell(NSC), matter, numpar) ! module "Electron_tools"
end subroutine construct_TB_H_Molteni


subroutine Complex_Hamil_tot_Molteni(numpar, Scell, NSC, atoms, TB, Hij, CHij, ksx, ksy, ksz)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB   ! parameters of the Hamiltonian of TB
   REAL(8), DIMENSION(:,:), allocatable, INTENT(inout), optional :: Hij	 ! Real Hamiltonian at Gamma-point
   real(8), intent(in) :: ksx, ksy, ksz ! k-point to get Hamiltonian [relative]
   complex, DIMENSION(:,:), allocatable, INTENT(inout), optional :: CHij ! Complex Hamiltonian at k-point
   integer nat, Ne, i, j, j1, i1, k, l, m, atom_2, NumTB, Nsiz, n1
   integer :: cell_x, cell_y, cell_z ! cell numbers
   real(8) :: x, y, z, r, sx, sy, sz, temp
   real(8) :: kx, ky, kz ! k-point to get Hamiltonian at, [1/A]
   REAL(8), DIMENSION(5,5) :: Hij1
   complex :: C_I, expfac

   nat = size(atoms) ! number of atoms
   Ne = Scell(NSC)%Ne  ! initial number of VB electrons
   if (present(Hij)) then
      Nsiz = size(Hij,1)
   else
      Nsiz = size(TB(1,1)%V0)*Scell(1)%Na
   endif

   if (numpar%optic_model .EQ. 2) then
!       temp = g_me*g_e/g_h*1d-10
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0
      if (.not.allocated(Scell(NSC)%cPRRx)) allocate(Scell(NSC)%cPRRx(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%cPRRy)) allocate(Scell(NSC)%cPRRy(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%cPRRz)) allocate(Scell(NSC)%cPRRz(Nsiz,Nsiz))
      Scell(NSC)%cPRRx = (0.0d0,0.0d0)
      Scell(NSC)%cPRRy = (0.0d0,0.0d0)
      Scell(NSC)%cPRRz = (0.0d0,0.0d0)
   endif

   call Reciproc(Scell(NSC)%supce, Scell(NSC)%k_supce) ! create reciprocal super-cell, module "Algebra_tools"
   call Reciproc_rel_to_abs(ksx, ksy, ksz, Scell, NSC, kx, ky, kz) ! get absolute k-values

   if (.not.present(Hij) .AND. .not.present(CHij)) then
      print*, 'Neither real, nor complex Hamiltonian is present,'
      print*, 'so, what do you expect this subroutine to do?'
      print*, 'The subroutine "Complex_Hamil_tot" from the module TB_Hamiltonian'
      print*, 'requires at least one of them to be present.'
   endif

   if (present(Hij) .AND. .not.present(CHij)) then ! just get real Hamiltonian:
      !call Hamil_tot(numpar, Scell, NSC, TB, Hij) ! see below
      call Hamil_tot_M(numpar, Scell, NSC, TB, Hij) ! see below
   endif

   if (.not.present(Hij) .AND. present(CHij)) then ! get complex Hamiltonian and the real one as well:
      
      C_I = dcmplx(0.0d0,1.0d0)	! complex unity
      if (.not.allocated(CHij)) then
         allocate(CHij(Nsiz,Nsiz))
      endif
      CHij = dcmplx(0.0d0,0.0d0)
      !if (.not.allocated(CHij1)) allocate(CHij1(Ne,Ne))
      do j = 1,nat	! all atoms
         m = Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one  
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            endif ! (atom_2 .EQ. 0)
            if (i .GT. 0) then
               !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)
               !call Hamilton_one(Scell, NSC, i, j, atoms, TB(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA), Hij1, x1=x, y1=y, z1=z, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z)
               call Hamilton_one_M(Scell, NSC, i, j, atoms, TB, Hij1,  x1=x, y1=y, z1=z, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z) ! below
               !expfac = exp(2.0d0*g_Pi*C_I*(kx*real(cell_x) + ky*real(cell_y) + kz*real(cell_z)))
               if ((abs(kx) .LT. 1.0d-14) .AND. (abs(ky) .LT. 1.0d-14) .AND. (abs(kz) .LT. 1.0d-14)) then
                  expfac = cmplx(1.0d0,0.0d0)
               else
!                   expfac = cexp(C_I*(kx*x + ky*y + kz*z))
                  expfac = exp(C_I*(kx*x + ky*y + kz*z))
               endif
               
               n1 = size(TB(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)%V0) ! that's how many orbitals per atom

               do j1 = 1,n1 ! all orbitals
                  l = (j-1)*n1+j1
                  do i1 = 1,n1 ! all orbitals
                     k = (i-1)*n1+i1

                     CHij(k,l) = DCMPLX(Hij1(j1,i1),0.0d0)*expfac ! make complex number of two real ones
                     
                     if (isnan(REAL(CHij(k,l))) .or. isnan(AIMAG(CHij(k,l)))) then
                        print*, k,l,CHij(k,l)
                        pause 'CHij pause'
                     endif

                     if (numpar%optic_model .EQ. 2) then ! create matrix element:
                        if (i .EQ. j) then ! for the dielectric function, according to Trani:
                           Scell(NSC)%cPRRx(k,l) = DCMPLX(temp*0.270d0,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = DCMPLX(temp*0.270d0,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = DCMPLX(temp*0.270d0,0.0d0)*CHij(k,l)
                        else
                           Scell(NSC)%cPRRx(k,l) = DCMPLX(temp*x,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = DCMPLX(temp*y,0.0d0)*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = DCMPLX(temp*z,0.0d0)*CHij(k,l)
                        endif
                        if (real(Scell(NSC)%cPRRx(k,l)) .GT. 1d10) print*, 'cPRRx', i, j, Scell(NSC)%cPRRx(k,l)
                        if (real(Scell(NSC)%cPRRy(k,l)) .GT. 1d10) print*, 'cPRRy',i, j, Scell(NSC)%cPRRy(k,l)
                        if (real(Scell(NSC)%cPRRz(k,l)) .GT. 1d10) print*, 'cPRRz',i, j, Scell(NSC)%cPRRz(k,l)
                     endif

                  enddo ! i1
               enddo !j1
            endif ! (i .GT. 0)
         enddo ! atom_2
      enddo ! j
   endif

   if (present(Hij) .AND. present(CHij)) then ! get complex Hamiltonian from the real one:
      
      C_I = dcmplx(0.0d0,1.0d0)	! complex unity
      if (.not.allocated(CHij)) then
         allocate(CHij(Nsiz,Nsiz))
      endif
      CHij = dcmplx(0.0d0,0.0d0)
      do j = 1,nat	! all atoms
         m = Scell(NSC)%Near_neighbor_size(j)
         do atom_2 = 0,m ! do only for atoms close to that one  
            if (atom_2 .EQ. 0) then
               i = j
            else
               i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
            endif ! (atom_2 .EQ. 0)
            if (i .GT. 0) then
               call shortest_distance(Scell, NSC, atoms, i, j, r, x1=x, y1=y, z1=z, cell_x=cell_x, cell_y=cell_y, cell_z=cell_z)
               !expfac = exp(2.0d0*g_Pi*C_I*(kx*real(cell_x) + ky*real(cell_y) + kz*real(cell_z)))
               if ((abs(kx) .LT. 1.0d-14) .AND. (abs(ky) .LT. 1.0d-14) .AND. (abs(kz) .LT. 1.0d-14)) then
                  expfac = dcmplx(1.0d0,0.0d0)
               else
!                   expfac = cexp(C_I*(kx*x + ky*y + kz*z))
                  expfac = exp(C_I*(kx*x + ky*y + kz*z))
               endif
               
               n1 = size(TB(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)%V0) ! that's how many orbitals per atom

               do j1 = 1,n1 ! all orbitals
                  l = (j-1)*n1+j1
                  do i1 = 1,n1 ! all orbitals
                     k = (i-1)*n1+i1

                     CHij(k,l) = DCMPLX(Hij(k,l),0.0d0)*expfac

                     if ((isnan(real(CHij(k,l)))) .OR. isnan(aimag(CHij(k,l)))) then
                        print*, i, j, k, l, CHij(k,l)
                        pause 'ISNAN'
                     endif

                     if (numpar%optic_model .EQ. 2) then ! create matrix element:
                        if (i .EQ. j) then ! for the dielectric function, according to Trani:
                           Scell(NSC)%cPRRx(k,l) = temp*0.270d0*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = temp*0.270d0*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = temp*0.270d0*CHij(k,l)
                        else
                           Scell(NSC)%cPRRx(k,l) = temp*x*CHij(k,l)
                           Scell(NSC)%cPRRy(k,l) = temp*y*CHij(k,l)
                           Scell(NSC)%cPRRz(k,l) = temp*z*CHij(k,l)
                        endif
                        if (real(Scell(NSC)%cPRRx(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRx(k,l)
                        if (real(Scell(NSC)%cPRRy(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRy(k,l)
                        if (real(Scell(NSC)%cPRRz(k,l)) .GT. 1d10) print*, i, j, Scell(NSC)%cPRRz(k,l)
                     endif
                  enddo ! i1
               enddo !j1
            endif ! (i .GT. 0)
         enddo ! atom_2
      enddo ! j
   endif

end subroutine Complex_Hamil_tot_Molteni


! hamiltonian for atoms:
subroutine Hamil_tot_M(numpar, Scell, NSC, TB_Hamil, Hij)
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   REAL(8), DIMENSION(:,:), INTENT(inout) :: Hij	 ! Hamiltonian
   REAL(8), DIMENSION(5,5) :: Hij1
   REAL(8), DIMENSION(5) :: Vij1
   integer :: nat, Ne, NumTB, Nsiz
   integer :: i, j, j1, i1, k, l, atom_2, m, FN, n1
   real(8) :: r, x, y, z, temp
   nat = Scell(NSC)%Na
   Ne = Scell(NSC)%Ne
   Nsiz = size(Hij,1)

   if (numpar%optic_model .EQ. 3) then
!       temp = g_me*g_e/g_h*1d-10
      ! mass and Plank constant cancel out in the final expression (subroutine get_Trani, module "Optical_parameters")
      temp = 1.0d0
      if (.not.allocated(Scell(NSC)%PRRx)) allocate(Scell(NSC)%PRRx(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%PRRy)) allocate(Scell(NSC)%PRRy(Nsiz,Nsiz))
      if (.not.allocated(Scell(NSC)%PRRz)) allocate(Scell(NSC)%PRRz(Nsiz,Nsiz))
      Scell(NSC)%PRRx = 0.0d0
      Scell(NSC)%PRRy = 0.0d0
      Scell(NSC)%PRRz = 0.0d0
   endif

   Hij = 0.0d0

   ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis:
   do j = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(j)
      do atom_2 = 0,m ! do only for atoms close to that one  
         if (atom_2 .EQ. 0) then
            i = j
         else
            i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms
         endif
         if (i .GT. 0) then ! if there really is a nearest neighbor

            ! Which kinds of atoms are these (which TB parameters to use):
            if (numpar%optic_model .EQ. 3) then ! create matrix element:
               !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z) ! module "Atomic_tools"
               call Hamilton_one_M(Scell, NSC, i, j, Scell(NSC)%MDatoms, TB_Hamil, Hij1, x1=x, y1=y, z1=z) ! block-hamiltonian and the shortest distances
            else
               call Hamilton_one_M(Scell, NSC, i, j, Scell(NSC)%MDatoms, TB_Hamil, Hij1) ! this calles the block-hamiltonian
            endif

            n1 = size(TB_Hamil(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)%V0) ! that's how many orbitals per atom

            do j1 = 1,n1 ! all orbitals
               l = (j-1)*n1+j1
               do i1 = 1,n1 ! all orbitals
                  k = (i-1)*n1+i1
                  !Hij((i-1)*4+i1, (j-1)*4+j1) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  Hij(k,l) = Hij1(i1,j1) ! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)

                  if (numpar%optic_model .EQ. 3) then ! create matrix element:
                     if (i .EQ. j) then ! for the dielectric function, according to Trani:
                        Scell(NSC)%PRRx(k,l) = temp*0.270d0*Hij(k,l)
                        Scell(NSC)%PRRy(k,l) = temp*0.270d0*Hij(k,l)
                        Scell(NSC)%PRRz(k,l) = temp*0.270d0*Hij(k,l)
                     else
                        Scell(NSC)%PRRx(k,l) = temp*x*Hij(k,l)
                        Scell(NSC)%PRRy(k,l) = temp*y*Hij(k,l)
                        Scell(NSC)%PRRz(k,l) = temp*z*Hij(k,l)
                     endif
                     if (Scell(NSC)%PRRx(k,l) .GT. 1d10) print*, 'Hamil_tot_M: Scell(NSC)%PRRx(k,l)', i, j, Scell(NSC)%PRRx(k,l)
                     if (Scell(NSC)%PRRy(k,l) .GT. 1d10) print*, 'Hamil_tot_M: Scell(NSC)%PRRx(k,l)', i, j, Scell(NSC)%PRRy(k,l)
                     if (Scell(NSC)%PRRz(k,l) .GT. 1d10) print*, 'Hamil_tot_M: Scell(NSC)%PRRx(k,l)', i, j, Scell(NSC)%PRRz(k,l)
!                      write(*,'(i4,i4,e,e,e)') k,l, matter%PRRx(k,l), matter%PRRy(k,l), matter%PRRz(k,l)
                  endif
               enddo ! i1
            enddo ! j1
         endif ! (i > 0)
      enddo ! j
   enddo ! i

!   open(NEWUNIT=FN, FILE = 'OUTPUT_Hamiltonain_Si.dat')
!      do i = 1, size(Hij,1)
!         do j = 1, size(Hij,2)
!            write(FN, '(i,i,e)') i, j, Hij(i,j)
!         enddo
!      enddo
!   close(FN)
!    pause 'Hamil_tot'

end subroutine Hamil_tot_M



subroutine Hamilton_one_M(Scell, NSC, i, j, atoms, TB, Hij,  x1, y1, z1, sx1, sy1, sz1, cell_x, cell_y, cell_z)
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   integer(4), INTENT(IN) :: i, j
   type(Atom), dimension(:), intent(in), target :: atoms	! array of atoms in the supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB	! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: Hij  ! hamiltonian
   real(8), intent(out), optional :: x1, y1, z1    ! shortest distances
   real(8), intent(out), optional :: sx1, sy1, sz1 ! shortest distance in relative coordinates
   integer, intent(out), optional :: cell_x, cell_y, cell_z ! cell numbers
   integer(4) ki, kj, i1, j1, k1, ik, n1
   REAL(8), DIMENSION(5,5) :: ts ! hopping integrals
   real(8) x,y,z,r,r1, x0, y0, z0, sx, sy, sz !  interatomic distance projections, and the total distance
   integer :: cell_x1, cell_y1, cell_z1 ! cell numbers
   integer, pointer :: KOA1, KOA2
   real(8), dimension(3) :: zb

   KOA1 => atoms(i)%KOA ! kind of first atom
   KOA2 => atoms(j)%KOA ! kind of second atom

   n1 = size(TB(KOA1,KOA2)%V0) ! that's how many orbitals per atom
   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      Hij = 0.0d0   ! Nondiagonals are zeros
      Hij(1,1) = TB(KOA1,KOA2)%Es ! diagonals are equal to onsite energies
      do ki = 2,4
         Hij(ki,ki) = TB(KOA1,KOA2)%Ep ! diagonals are equal to onsite energies 
      enddo
      Hij(5,5) = TB(KOA1,KOA2)%Esa ! diagonals are equal to onsite energies
      if (present(x1)) x1 = 0.0d0
      if (present(y1)) y1 = 0.0d0
      if (present(z1)) z1 = 0.0d0
      if (present(sx1)) sx1 = 0.0d0
      if (present(sy1)) sy1 = 0.0d0
      if (present(sz1)) sz1 = 0.0d0
      if (present(cell_x)) cell_x = 0
      if (present(cell_y)) cell_y = 0
      if (present(cell_z)) cell_z = 0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call shortest_distance(Scell, NSC, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz, cell_x=cell_x1, cell_y=cell_y1, cell_z=cell_z1) ! module "Atomic_tools"
      if (present(x1)) x1 = x
      if (present(y1)) y1 = y
      if (present(z1)) z1 = z
      if (present(sx1)) sx1 = sx
      if (present(sy1)) sy1 = sy
      if (present(sz1)) sz1 = sz
      if (present(cell_x)) cell_x = cell_x1
      if (present(cell_y)) cell_y = cell_y1
      if (present(cell_z)) cell_z = cell_z1
      !call Hopping_old(TB, ts, -x, -y, -z)
      call Hopping_M(TB, KOA1, KOA2, ts, x, y, z)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1,n1
         do kj = 1,n1
            Hij(kj,ki) = ts(kj,ki)   ! Hopping Integrals
         enddo ! kj
      enddo  ! ki
   endif
   nullify(KOA1, KOA2)
end subroutine Hamilton_one_M



subroutine Hopping_M(TB, N1, N2, ts, x, y, z)
! subroutine making the hopping integrals
! ts --- the Hopping Integrals themselves
! x,y,z --- the distances between naighbour atoms we are analizing
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB	! all tight binding parameters
   integer, intent(in) :: N1, N2 ! kind of atom 1 and 2
   REAL(8), DIMENSION(5,5), INTENT(out) :: ts
   real(8), intent(in) :: x,y,z	! relative distances (projections) between the atoms
   !=============================================
   INTEGER(4) :: i,j, k
   INTEGER(4) :: i1,j1
   real(8) r, r1, xr, yr, zr ! Vsr(5), Vsr2, Vsr5
   real(8) :: Vsr2, Vsr5
   real(8), dimension(5) :: Vsr
   r = SQRT(x*x + y*y + z*z)   ! total distance between the two atoms within the supercell
   ! All of the hopping integrals ts for each pair of atoms i and j,
   ! according to the Eqs.(2.23)-(2.29) from PhD Thesis of H.Jeschke, Page 38
   xr = x/r
   yr = y/r
   zr = z/r
   Vsr(1) = Vs_M(TB,1,N1,N2,r)	! calling the function given below, Vss
   Vsr(2) = Vs_M(TB,2,N1,N2,r)	! calling the function given below, Vsp
   Vsr2   = Vs_M(TB,2,N2,N1,r)	! calling the function given below, Vps
   Vsr(3) = Vs_M(TB,3,N1,N2,r)	! calling the function given below, VppSigma
   Vsr(4) = Vs_M(TB,4,N1,N2,r)	! calling the function given below, VppPi
   Vsr(5) = Vs_M(TB,5,N1,N2,r)	! calling the function given below, Vs*p
   Vsr5   = Vs_M(TB,5,N2,N1,r)	! calling the function given below, Vps*

!    Vsr(1) = Vs_M(TB,1,N1,N2,r,.false.)	! calling the function given below, Vss
!    Vsr(2) = Vs_M(TB,2,N1,N2,r,.false.)	! calling the function given below, Vsp
!    Vsr2   = Vs_M(TB,2,N2,N1,r,.false.)	! calling the function given below, Vps
!    Vsr(3) = Vs_M(TB,3,N1,N2,r,.false.)	! calling the function given below, VppSigma
!    Vsr(4) = Vs_M(TB,4,N1,N2,r,.false.)	! calling the function given below, VppPi
!    Vsr(5) = Vs_M(TB,5,N1,N2,r,.false.)	! calling the function given below, Vs*p
!    Vsr5   = Vs_M(TB,5,N2,N1,r,.false.)	! calling the function given below, Vps*


   ts(1,1) = Vsr(1)
   ts(1,2) = -xr*Vsr(2)
   ts(1,3) = -yr*Vsr(2)
   ts(1,4) = -zr*Vsr(2)
   ts(1,5) = 0.0d0

   ts(2,1) = xr*Vsr2
   ts(2,2) = xr*xr*(Vsr(3) - Vsr(4)) + Vsr(4)
   ts(2,3) = xr*yr*(Vsr(3) - Vsr(4))
   ts(2,4) = xr*zr*(Vsr(3) - Vsr(4))
   ts(2,5) = xr*Vsr5

   ts(3,1) = yr*Vsr2
   ts(3,2) = xr*yr*(Vsr(3) - Vsr(4))
   ts(3,3) = yr*yr*(Vsr(3) - Vsr(4)) + Vsr(4)
   ts(3,4) = yr*zr*(Vsr(3) - Vsr(4))
   ts(3,5) = yr*Vsr5

   ts(4,1) = zr*Vsr2
   ts(4,2) = xr*zr*(Vsr(3) - Vsr(4))
   ts(4,3) = yr*zr*(Vsr(3) - Vsr(4))
   ts(4,4) = zr*zr*(Vsr(3) - Vsr(4)) + Vsr(4)
   ts(4,5) = zr*Vsr5

   ts(5,1) = 0.0d0
   ts(5,2) = -xr*Vsr(5)
   ts(5,3) = -yr*Vsr(5)
   ts(5,4) = -zr*Vsr(5)
   ts(5,5) = 0.0d0

end subroutine Hopping_M


function Vs_M(TB,i,N1,N2,r,cut_off) ! functions Vsssigma, /p sigma, /p pi
! for the Hopping Integrals according to Eqs.(2.23)-(2.29), (2.32), 
! and the Appendix E with all parameters,
! from H.Jeschke PhD thesis.
! i   ---  stands for index of V, i=1 is "s s sigma", i=2 is "s p sigma", i=3 is "p p sigma", i=4 is "p p pi"
! x, y, z --- are the coordinates of interatomic distance (positive or negative!)
! r  --- absolute interatomic distance
   real(8) Vs_M ! function itself (Eqs.(2.23)-(2.27))
   type(TB_H_Molteni), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   INTEGER(4), intent(in) :: i ! index of V
   INTEGER(4), intent(in) :: N1, N2  ! kinds of 1st and 2d atoms
   real(8), intent(in) :: r  ! absolute distance between the atoms [A]
   logical, intent(in), optional :: cut_off ! do or don't do cut off?
   !--------------------------------------------
   real(8), dimension(:), pointer :: V0
   real(8), pointer :: nc, rc, rcut, d, r0, n
   real(8) :: r01
   integer j

   V0 => TB(N1,N2)%V0
   n => TB(N1,N2)%n
   rcut => TB(N1,N2)%rcut
   d => TB(N1,N2)%d
   r0 => TB(N1,N2)%r0
   nc => TB(N1,N2)%nc
   rc => TB(N1,N2)%rc

   ! Unconvenional scaling:
   !Vs_M = V0(i)*((r0/r)**n)*exp(-n*((r/rc)**nc - (r0/rc)**nc))
   ! Harrison scaling [Molteni et al., J. Phys.: Condens. Matter 6 (1994) 5243]:
   Vs_M = V0(i)*((r0/r)**n)

   ! Introduce cut-off function:
   if (present(cut_off)) then
      if (cut_off) Vs_M = Vs_M/(1.0d0+exp((r-rcut)*2.0d0/d)) ! cutting off long-distance atoms
   else
      Vs_M = Vs_M/(1.0d0+exp((r-rcut)*2.0d0/d)) ! cutting off long-distance atoms
   endif

   if (isnan(Vs_M) .or. (abs(Vs_M) > 1e25)) then 
      print*, 'Vs_M is:', Vs_M, 1.0d0/(1.0d0+exp((r-rcut)*2.0d0/d)), i, V0(i), r0, r, n, (r0/r)**n
      print*, N1, N2, 'TB:', TB(N1,N2)
      pause 'Vs_M'
   endif
   nullify(V0, nc, rc, rcut, d, r0, n)
end function Vs_M


!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

! Derivative of the atomic TB Hamiltonian:
subroutine dHij_s_M(TB_Hamil, atoms, Scell, NSC, numpar, Aij, M_x1, M_xrr) ! attractive forces for atoms, module "TB_Hamiltonian"
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors   
   REAL(8), DIMENSION(:,:), INTENT(in) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis
   REAL(8), DIMENSION(:,:,:), INTENT(in) :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in) :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   !=========================================
   REAL(8), DIMENSION(3) :: Eelectr_s
   real(8), dimension(:,:), allocatable :: dHijx_s_all ! allocatable
   real(8), dimension(:,:), allocatable :: dHijy_s_all ! allocatable
   real(8), dimension(:,:), allocatable :: dHijz_s_all ! allocatable
   real(8), dimension(:,:,:), allocatable :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), allocatable :: M_dVs ! matrix of functions dVs
   integer k, nat, nat4, my_id, OMP_GET_THREAD_NUM
   nat = size(atoms)
   nat4 = size(Scell(NSC)%Ha,1) ! size of the Hamiltonian depends on the basis set
   
   ! Construct array of functions Vs and dVs for all pairs of atoms to use for forces:
   call Construct_M_Vs_M(Scell, NSC, TB_Hamil, M_Vs, M_dVs) ! subroitine below

!$omp PARALLEL private(k, Eelectr_s, dHijx_s_all, dHijy_s_all, dHijz_s_all) 
   if (.not.allocated(dHijx_s_all)) allocate(dHijx_s_all(nat4,nat4))
   if (.not.allocated(dHijy_s_all)) allocate(dHijy_s_all(nat4,nat4))
   if (.not.allocated(dHijz_s_all)) allocate(dHijz_s_all(nat4,nat4))   
!$omp do
   do k = 1,nat ! initial conditions for atoms:
!      my_id = OMP_GET_THREAD_NUM() ! identify which thread it is
      Scell(NSC)%MDatoms(k)%forces%att(:) = 0.0d0
      
      call dHamil_tot_s_M(dHijx_s_all, dHijy_s_all, dHijz_s_all, TB_Hamil, Scell, NSC, numpar, k, M_x1, M_xrr, M_Vs, M_dVs) ! see below

!       call check_symmetry(dHijx_s_all) ! just for testing, must be hermitian Hamiltonian
!       call check_symmetry(dHijy_s_all) ! just for testing, must be hermitian Hamiltonian
!       call check_symmetry(dHijz_s_all) ! just for testing, must be hermitian Hamiltonian

      call Attract_TB_H3_near_M(Aij, dHijx_s_all, dHijy_s_all, dHijz_s_all, Scell, NSC, Eelectr_s) ! see below
      Scell(NSC)%MDatoms(k)%forces%att(:) = Eelectr_s(:) ! save attractive forces

   enddo ! k
!$omp end do 
   if (allocated(dHijx_s_all)) deallocate(dHijx_s_all)
   if (allocated(dHijy_s_all)) deallocate(dHijy_s_all)
   if (allocated(dHijz_s_all)) deallocate(dHijz_s_all)
!$omp end parallel

   if (allocated(M_Vs)) deallocate(M_Vs)
   if (allocated(M_dVs)) deallocate(M_dVs)
end subroutine dHij_s_M




subroutine Construct_M_Vs_M(Scell, NSC, TB_Hamil, M_Vs, M_dVs)
   type(Super_cell), dimension(:), intent(in), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil   ! parameters of the Hamiltonian of TB
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), allocatable, intent(out) :: M_dVs ! matrix of functions dVs
   !---------------------------
   real(8), pointer :: r
   integer, pointer :: m,  j, nat
   integer i, atom_2
   nat => Scell(NSC)%Na ! number of atoms
   if (.not.allocated(M_Vs)) allocate(M_Vs(7,nat,nat))
   if (.not.allocated(M_dVs)) allocate(M_dVs(7,nat,nat))
   M_Vs = 0.0d0
   M_dVs = 0.0d0
   !$omp PARALLEL DO private(i, m, atom_2, j, r) 
   do i = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one  
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R
            M_dVs(1,i,j) = dVs_M(TB_Hamil,1,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r,1) ! function below
            M_dVs(2,i,j) = dVs_M(TB_Hamil,2,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r,1) ! function below
            M_dVs(3,i,j) = dVs_M(TB_Hamil,2,Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA,r,1) ! function below
            M_dVs(4,i,j) = dVs_M(TB_Hamil,3,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r,1) ! function below
            M_dVs(5,i,j) = dVs_M(TB_Hamil,4,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r,1) ! function below
            M_dVs(6,i,j) = dVs_M(TB_Hamil,5,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r,1) ! function below
            M_dVs(7,i,j) = dVs_M(TB_Hamil,5,Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA,r,1) ! function below
            
            M_Vs(1,i,j) = Vs_M(TB_Hamil,1,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r) ! function below
            M_Vs(2,i,j) = Vs_M(TB_Hamil,2,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r) ! function below
            M_Vs(3,i,j) = Vs_M(TB_Hamil,2,Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA,r) ! function below
            M_Vs(4,i,j) = Vs_M(TB_Hamil,3,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r) ! function below
            M_Vs(5,i,j) = Vs_M(TB_Hamil,4,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r) ! function below
            M_Vs(6,i,j) = Vs_M(TB_Hamil,5,Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA,r) ! function below
            M_Vs(7,i,j) = Vs_M(TB_Hamil,5,Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA,r) ! function below
         endif
      enddo
   enddo
   !$omp END PARALLEL DO
   nullify(m, j, nat, r)
end subroutine Construct_M_Vs_M



subroutine Attract_TB_H3_near_M(Aij, dHijx_s, dHijy_s, dHijz_s, Scell, NSC, Eelectr_s)
   REAL(8), DIMENSION(:,:), intent(in) :: dHijx_s, dHijy_s, dHijz_s
   REAL(8), DIMENSION(:,:), INTENT(in), target :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD thesis 
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   REAL(8), DIMENSION(:), intent(out)  :: Eelectr_s ! part of the forces

   
   integer, pointer :: m, j1
   real(8), pointer :: Aijij, Aijij1, Aijij2, Aijij3, Aijij4
   integer i, j, k, n, i2, ste
   real(8) E_sum(3)

   n = size(Aij,1)
   Eelectr_s = 0.0d0

   i2 = 0
   ste = 1
   do i = 1, n ! all orbitals
       if (i .GE. ste) then
          i2 = i2 + 1
          ste = ste + 5
       endif
       m => Scell(NSC)%Near_neighbor_size(i2)
       do k = 1, m
          j1 => Scell(NSC)%Near_neighbor_list(i2,k) ! this is the list of such close atoms
          j = (j1-1)*5 + 1
          if (j .GT. 0) then
             E_sum = 0.0d0
             Aijij => Aij(i,j)
             Aijij1 => Aij(i,j+1)
             Aijij2 => Aij(i,j+2)
             Aijij3 => Aij(i,j+3)
             Aijij4 => Aij(i,j+4)

             E_sum(3) = (dHijz_s(i,j)*Aijij) + (dHijz_s(i,j+1)*Aijij1) + (dHijz_s(i,j+2)*Aijij2) + (dHijz_s(i,j+3)*Aijij3) + (dHijz_s(i,j+4)*Aijij4)

             E_sum(2) = (dHijy_s(i,j)*Aijij) + (dHijy_s(i,j+1)*Aijij1) + (dHijy_s(i,j+2)*Aijij2) + (dHijy_s(i,j+3)*Aijij3) + (dHijy_s(i,j+4)*Aijij4)

             !PROBLEM IN dHijx_s(i,j+3)???

             E_sum(1) = (dHijx_s(i,j)*Aijij) + (dHijx_s(i,j+1)*Aijij1) + (dHijx_s(i,j+2)*Aijij2) + (dHijx_s(i,j+3)*Aijij3) + (dHijx_s(i,j+4)*Aijij4)

             Eelectr_s(:) = Eelectr_s(:) + E_sum(:) ! sum for all 3 coordinates

             if (minval(ABS(E_sum(:))) .GE. 1.0d6) then
                write(*,'(a)') 'Trouble in subroutine Attract_TB_H3_near_M, too large attractive force:'
                write(*,'(i12,i12,f25.16,f25.16,f25.16)') i, j, E_sum(:)
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j), dHijy_s(i,j), dHijz_s(i,j), Aijij
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+1), dHijy_s(i,j+1), dHijz_s(i,j+1), Aijij1
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+2), dHijy_s(i,j+2), dHijz_s(i,j+2), Aijij2
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+3), dHijy_s(i,j+3), dHijz_s(i,j+3), Aijij3
                write(*,'(i12,i12,f25.16,f25.16,f25.16,e25.16)') i, j, dHijx_s(i,j+4), dHijy_s(i,j+4), dHijz_s(i,j+4), Aijij4
                write(*,'(e25.16,e25.16,e25.16)') Eelectr_s(:)
             endif

!              if (minval(ABS(E_sum(:))) .GE. 1.0d-12) write(*,'(i4,i4,f,f,f,f,f,f)') i, j, E_sum(:), Eelectr_s(:)

          endif
       enddo
   enddo
   nullify(m, j1, Aijij, Aijij1, Aijij2, Aijij3, Aijij4)
end subroutine Attract_TB_H3_near_M



subroutine dHamil_tot_s_M(dHijx, dHijy, dHijz, TB_Hamil, Scell, NSC, numpar, k, M_x1, M_xrr, M_Vs, M_dVs)
! construct the whole Hamilton matrix:
! it uses input of Hij, the Hamiltonain block,
! input Es and Ep - the onsite atomic energies
! X - the array of coordinates of all atoms
! nat - the number of atoms
! supce - the size of the supercell, used for periodical boundary conditions
! k - is the atom for which we calculate the forces 
! (with respect to which Rk we take the derivatives, Appendix F of H.Jeschke PhD Thesis)
   integer, INTENT(IN) :: k ! number of atom
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   REAL(8), DIMENSION(:,:), INTENT(out) :: dHijx, dHijy, dHijz
   REAL(8), DIMENSION(:,:,:), INTENT(in) :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in) :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   !-----------------------------------------------
   REAL(8), DIMENSION(5,5) :: dHijx1, dHijy1, dHijz1
   integer i, j1, i1, ki, atom_2, l, k1, n1
   integer, pointer :: nat, m, j

   nat => Scell(NSC)%Na
   dHijx = 0.0d0
   dHijy = 0.0d0
   dHijz = 0.0d0

   do i = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one
         j => Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            dHijx1 = 0.0d0
            dHijy1 = 0.0d0
            dHijz1 = 0.0d0
            call dHamilton_one_s_M(TB_Hamil, i, j, k, dHijz1, dHijy1, dHijx1, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)

            n1 = size(TB_Hamil(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)%V0) ! that's how many orbitals per atom
            do j1 = 1,n1 ! all orbitals
               l = (j-1)*n1+j1
               do i1 = 1,n1 ! all orbitals
                  k1 = (i-1)*n1+i1

                  dHijx(k1,l) = dHijx1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  dHijy(k1,l) = dHijy1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
                  dHijz(k1,l) = dHijz1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)

                  ! This is WRONG (!!!):
!                 dHijx(l,k1) = dHijx1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
!                 dHijy(l,k1) = dHijy1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)
!                 dHijz(l,k1) = dHijz1(i1,j1)	! construct the total Hamiltonian from the blocks of one-atom Hamiltonian, Eq.(2.40)

!                    if (isnan(dHijx(k1,l))) write(*,'(a,i4,i4,i4,i4,es)') 'x ', k1, l, j1, i1, dHijx(k1,l)
!                    if (isnan(dHijy(k1,l))) write(*,'(a,i4,i4,i4,i4,es)') 'y ', k1, l, j1, i1, dHijy(k1,l)
!                    if (isnan(dHijz(k1,l))) write(*,'(a,i4,i4,i4,i4,es)') 'z ', k1, l, j1, i1, dHijz(k1,l)
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0)
      enddo ! j
   enddo ! i
   nullify(nat, m, j)
end subroutine dHamil_tot_s_M



subroutine dHamilton_one_s_M(TB, i, j, k, dHijz, dHijy, dHijx, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB ! all tight binding parameters
   integer, INTENT(IN) :: i, j, k, atom_2
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   REAL(8), DIMENSION(:,:), INTENT(out) :: dHijx, dHijy, dHijz
   REAL(8), DIMENSION(:,:,:), INTENT(in) :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in) :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), intent(in) :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in) :: M_dVs ! matrix of functions dVs
   !---------------------------------------
   integer(4) ki, kj, i1, j1, k1, n1
   REAL(8), DIMENSION(5,5) :: dtsx, dtsy, dtsz ! hopping integrals
   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      dHijz = 0.0d0   ! Nondiagonals are zeros
      dHijy = 0.0d0   ! Nondiagonals are zeros
      dHijx = 0.0d0   ! Nondiagonals are zeros
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_s_M(dtsz, dtsy, dtsx, i, j, k, TB, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      n1 = size(TB(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)%V0) ! that's how many orbitals per atom
      do ki = 1,n1
         do kj = 1,n1 ! Checked, correct:
              dHijx(kj,ki) = dtsx(ki,kj)   ! Hopping Integrals
              dHijy(kj,ki) = dtsy(ki,kj)   ! Hopping Integrals
              dHijz(kj,ki) = dtsz(ki,kj)   ! Hopping Integrals
!              if (ABS(dHijx(kj,ki)) .LE. 1d-14) then
!                 dHijx(kj,ki) = 0.0d0
!              endif
!              if (ABS(dHijy(kj,ki)) .LE. 1d-14) then
!                 dHijy(kj,ki) = 0.0d0
!              endif
!              if (ABS(dHijz(kj,ki)) .LE. 1d-14) then
!                 dHijz(kj,ki) = 0.0d0
!              endif

!              dHijz(ki,kj) = dtsz(ki,kj)   ! Hopping Integrals, checked WRONG (!!!)
!              dHijy(ki,kj) = dtsy(ki,kj)   ! Hopping Integrals
!              dHijx(ki,kj) = dtsx(ki,kj)   ! Hopping Integrals
!             if (ABS(dHijx(ki,kj)) .LE. 1d-14) then
!                dHijx(ki,kj) = 0.0d0
!             endif
!             if (ABS(dHijy(ki,kj)) .LE. 1d-14) then
!                dHijy(ki,kj) = 0.0d0
!             endif
!             if (ABS(dHijz(ki,kj)) .LE. 1d-14) then
!                dHijz(ki,kj) = 0.0d0
!             endif

         enddo ! kj
      enddo  ! ki
   endif
end subroutine dHamilton_one_s_M



subroutine dHopping_s_M(dtsz, dtsy, dtsx, i, j, k, TB, Scell, NSC, atom_2, M_x1, M_xrr, M_Vs, M_dVs)
! subroutine making the derivatives of the hopping integrals
! dtsz, dtsy, dtsx --- the derivatives Hopping Integrals in Z,Y,X axis
! Xco --- atomic coordinates of all the atoms, 
! supce --- the supercell size
! i, j - a pair of atoms, with respect to which we calculate a force acting on 
! the atom "k".
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB ! all tight binding parameters
   REAL(8), DIMENSION(:,:), INTENT(out) :: dtsx, dtsy, dtsz
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   INTEGER(4), intent(in) :: i,j,k, atom_2
   REAL(8), DIMENSION(:,:,:), INTENT(in), target :: M_x1  ! Matrix of x1 elements, used for forces
   real(8), dimension(:,:,:), intent(in), target :: M_xrr ! matrix of coefficients xrr, yrr, zrr
   real(8), dimension(:,:,:), intent(in), target :: M_Vs  ! matrix of functions Vs
   real(8), dimension(:,:,:), intent(in), target :: M_dVs ! matrix of functions dVs
   !-----------------------------------------------
   integer, pointer :: KOA_i, KOA_j
   real(8), pointer :: x,y,z,r, xrr, yrr, zrr
   REAL(8), DIMENSION(:), pointer :: x1
   real(8), pointer :: dVs1, dVs2, dVs21, dVs3, dVs4, dVs5, dVs51
   real(8), pointer :: Vs1, Vs2, Vs21, Vs3, Vs4, Vs5, Vs51
   
   real(8) x0, y0, z0, sx, sy, sz	! relative distances (projections) between the atoms 
   real(8) r1, temp, temp1, temp2, z_temp
   INTEGER(4) :: i1, j1, k1, dik, djk, dij, ik
   REAL(8), DIMENSION(3) :: dlds, dmds, dnds, zb
   real(8) drdsx, drdsy, drdsz, b, xr, yr, zr

   if (i .EQ. k) then
      dik = 1
   else
      dik = 0
   endif
   if (k .EQ. j) then
      djk = 1
   else
      djk = 0
   endif

   if ((dik - djk) /= 0) then ! only then it is non-zero

      KOA_i => Scell(NSC)%MDatoms(i)%KOA ! kind of atom for atom i
      KOA_j => Scell(NSC)%MDatoms(j)%KOA ! kind of atom for atom j

      x => Scell(NSC)%Near_neighbor_dist(i,atom_2,1) ! at this distance, X
      y => Scell(NSC)%Near_neighbor_dist(i,atom_2,2) ! at this distance, Y
      z => Scell(NSC)%Near_neighbor_dist(i,atom_2,3) ! at this distance, Z
      r => Scell(NSC)%Near_neighbor_dist(i,atom_2,4) ! at this distance, R

      b = dble(dik - djk)/r
      
      x1 => M_x1(:,i,j) !x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
      
!       x1(1) = M_x1(1,i,j) !x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
!       x1(2) = M_x1(2,i,j) !x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
!       x1(3) = M_x1(3,i,j) !x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)

      
      drdsx = b*x1(1)
      drdsy = b*x1(2)
      drdsz = b*x1(3)
      xrr => M_xrr(1,i,j) !x/(r*r)
      yrr => M_xrr(2,i,j) !y/(r*r)
      zrr => M_xrr(3,i,j) !z/(r*r)

      dlds(1) = b*Scell(NSC)%supce(1,1) - xrr*drdsx
      dlds(2) = b*Scell(NSC)%supce(2,1) - xrr*drdsy
      dlds(3) = b*Scell(NSC)%supce(3,1) - xrr*drdsz
      dmds(1) = b*Scell(NSC)%supce(1,2) - yrr*drdsx
      dmds(2) = b*Scell(NSC)%supce(2,2) - yrr*drdsy
      dmds(3) = b*Scell(NSC)%supce(3,2) - yrr*drdsz
      dnds(1) = b*Scell(NSC)%supce(1,3) - zrr*drdsx
      dnds(2) = b*Scell(NSC)%supce(2,3) - zrr*drdsy
      dnds(3) = b*Scell(NSC)%supce(3,3) - zrr*drdsz
      ! All of the hopping integrals ts for each pair of atoms i and j,
      ! according to the Eqs.(2.23)-(2.29) from PhD Thesis of H.Jeschke, Page 38
      xr = (x/r)
      yr = (y/r)
      zr = (z/r)

      dVs1  => M_dVs(1,i,j) !dVs_M(TB,1,KOA_i,KOA_j,r,1)
      dVs2  => M_dVs(2,i,j) !dVs_M(TB,2,KOA_i,KOA_j,r,1)
      dVs21 => M_dVs(3,i,j) !dVs_M(TB,2,KOA_j,KOA_i,r,1)
      dVs3  => M_dVs(4,i,j) !dVs_M(TB,3,KOA_i,KOA_j,r,1)
      dVs4  => M_dVs(5,i,j) !dVs_M(TB,4,KOA_i,KOA_j,r,1)
      dVs5  => M_dVs(6,i,j) !dVs_M(TB,5,KOA_i,KOA_j,r,1)
      dVs51 => M_dVs(7,i,j) !dVs_M(TB,5,KOA_j,KOA_i,r,1)
      
!       if ( abs(dVs1 - dVs_M(TB,1,KOA_i,KOA_j,r,1)) > 1.0d-14) print*, 'dVs1', dVs1, dVs_M(TB,1,KOA_i,KOA_j,r,1)
!       if ( abs(dVs2 - dVs_M(TB,2,KOA_i,KOA_j,r,1)) > 1.0d-14) print*, 'dVs2', dVs2, dVs_M(TB,2,KOA_i,KOA_j,r,1)
!       if ( abs(dVs21 - dVs_M(TB,2,KOA_j,KOA_i,r,1)) > 1.0d-14) print*, 'dVs21', dVs21, dVs_M(TB,2,KOA_j,KOA_i,r,1)
!       if ( abs(dVs3 - dVs_M(TB,3,KOA_i,KOA_j,r,1)) > 1.0d-14) print*, 'dVs3', dVs3, dVs_M(TB,3,KOA_i,KOA_j,r,1)
!       if ( abs(dVs4 - dVs_M(TB,4,KOA_i,KOA_j,r,1)) > 1.0d-14) print*, 'dVs4', dVs4, dVs_M(TB,4,KOA_i,KOA_j,r,1)
!       if ( abs(dVs5 - dVs_M(TB,5,KOA_i,KOA_j,r,1)) > 1.0d-14) print*, 'dVs5', dVs5, dVs_M(TB,5,KOA_i,KOA_j,r,1)
!       if ( abs(dVs51 - dVs_M(TB,5,KOA_j,KOA_i,r,1)) > 1.0d-14) print*, 'dVs51', dVs51, dVs_M(TB,5,KOA_j,KOA_i,r,1)
      
      Vs1  => M_Vs(1,i,j) !Vs_M(TB,1,KOA_i,KOA_j,r)
      Vs2  => M_Vs(2,i,j) !Vs_M(TB,2,KOA_i,KOA_j,r)
      Vs21 => M_Vs(3,i,j) !Vs_M(TB,2,KOA_j,KOA_i,r)
      Vs3  => M_Vs(4,i,j) !Vs_M(TB,3,KOA_i,KOA_j,r)
      Vs4  => M_Vs(5,i,j) !Vs_M(TB,4,KOA_i,KOA_j,r)
      Vs5  => M_Vs(6,i,j) !Vs_M(TB,5,KOA_i,KOA_j,r)
      Vs51 => M_Vs(7,i,j) !Vs_M(TB,5,KOA_j,KOA_i,r)

!       if ( abs(Vs1 - Vs_M(TB,1,KOA_i,KOA_j,r)) > 1.0d-14) print*, 'Vs1', Vs1, Vs_M(TB,1,KOA_i,KOA_j,r)
!       if ( abs(Vs2 - Vs_M(TB,2,KOA_i,KOA_j,r)) > 1.0d-14) print*, 'Vs2', Vs2, Vs_M(TB,2,KOA_i,KOA_j,r)
!       if ( abs(Vs21 - Vs_M(TB,2,KOA_j,KOA_i,r)) > 1.0d-14) print*, 'Vs21', Vs21, Vs_M(TB,2,KOA_j,KOA_i,r)
!       if ( abs(Vs3 - Vs_M(TB,3,KOA_i,KOA_j,r)) > 1.0d-14) print*, 'Vs3', Vs3, Vs_M(TB,3,KOA_i,KOA_j,r)
!       if ( abs(Vs4 - Vs_M(TB,4,KOA_i,KOA_j,r)) > 1.0d-14) print*, 'Vs4', Vs4, Vs_M(TB,4,KOA_i,KOA_j,r)
!       if ( abs(Vs5 - Vs_M(TB,5,KOA_i,KOA_j,r)) > 1.0d-14) print*, 'Vs5', Vs5, Vs_M(TB,5,KOA_i,KOA_j,r)
!       if ( abs(Vs51 - Vs_M(TB,5,KOA_j,KOA_i,r)) > 1.0d-14) print*, 'Vs51', Vs51, Vs_M(TB,5,KOA_j,KOA_i,r)
      

      dtsx(1,1) = dVs1*drdsx ! s s
      dtsy(1,1) = dVs1*drdsy ! s s
      dtsz(1,1) = dVs1*drdsz ! s s

      temp = xr*dVs2
      dtsx(2,1) = -(dlds(1)*Vs2 + temp*drdsx) ! -s px
      dtsy(2,1) = -(dlds(2)*Vs2 + temp*drdsy) ! -s px
      dtsz(2,1) = -(dlds(3)*Vs2 + temp*drdsz) ! -s px

      temp = xr*dVs21
      dtsx(1,2) = (dlds(1)*Vs21 + temp*drdsx) ! px s
      dtsy(1,2) = (dlds(2)*Vs21 + temp*drdsy) ! px s
      dtsz(1,2) = (dlds(3)*Vs21 + temp*drdsz) ! px s

      temp = yr*dVs2
      dtsx(3,1) = -(dmds(1)*Vs2 + temp*drdsx) ! -s py
      dtsy(3,1) = -(dmds(2)*Vs2 + temp*drdsy) ! -s py
      dtsz(3,1) = -(dmds(3)*Vs2 + temp*drdsz) ! -s py

      temp = yr*dVs21
      dtsx(1,3) = (dmds(1)*Vs21 + temp*drdsx) ! s py
      dtsy(1,3) = (dmds(2)*Vs21 + temp*drdsy) ! s py
      dtsz(1,3) = (dmds(3)*Vs21 + temp*drdsz) ! s py

      temp = zr*dVs2
      dtsx(4,1) = -(dnds(1)*Vs2 + temp*drdsx)  ! -s pz
      dtsy(4,1) = -(dnds(2)*Vs2 + temp*drdsy)  ! -s pz
      dtsz(4,1) = -(dnds(3)*Vs2 + temp*drdsz)  ! -s pz

      temp = zr*dVs21
      dtsx(1,4) = (dnds(1)*Vs21 + temp*drdsx)  ! s pz
      dtsy(1,4) = (dnds(2)*Vs21 + temp*drdsy)  ! s pz
      dtsz(1,4) = (dnds(3)*Vs21 + temp*drdsz)  ! s pz

      temp = 2.0d0*xr*(Vs3 - Vs4)
      temp1 = (xr*xr*(dVs3 - dVs4) + dVs4)
      dtsx(2,2) = temp*dlds(1) + temp1*drdsx ! px px
      dtsy(2,2) = temp*dlds(2) + temp1*drdsy ! px px
      dtsz(2,2) = temp*dlds(3) + temp1*drdsz ! px px

      temp = Vs3 - Vs4
      temp1 = xr*yr*(dVs3 - dVs4)
      dtsx(3,2) = (xr*dmds(1) + yr*dlds(1))*temp + temp1*drdsx ! px py
      dtsy(3,2) = (xr*dmds(2) + yr*dlds(2))*temp + temp1*drdsy ! px py
      dtsz(3,2) = (xr*dmds(3) + yr*dlds(3))*temp + temp1*drdsz ! px py

      dtsx(2,3) = dtsx(3,2)
      dtsy(2,3) = dtsy(3,2)
      dtsz(2,3) = dtsz(3,2)

      temp1 = xr*zr*(dVs3 - dVs4)
      dtsx(4,2) = (xr*dnds(1) + zr*dlds(1))*temp + temp1*drdsx ! px pz
      dtsy(4,2) = (xr*dnds(2) + zr*dlds(2))*temp + temp1*drdsy ! px pz
      dtsz(4,2) = (xr*dnds(3) + zr*dlds(3))*temp + temp1*drdsz ! px pz

      dtsx(2,4) = dtsx(4,2)
      dtsy(2,4) = dtsy(4,2)
      dtsz(2,4) = dtsz(4,2)

      temp2 = 2.0d0*yr*(Vs3 - Vs4)
      temp1 = yr*yr*(dVs3 - dVs4) + dVs4
      dtsx(3,3) = temp2*dmds(1) + temp1*drdsx ! py py
      dtsy(3,3) = temp2*dmds(2) + temp1*drdsy ! py py
      dtsz(3,3) = temp2*dmds(3) + temp1*drdsz ! py py

      temp1 = yr*zr*(dVs3 - dVs4)
      dtsx(4,3) = (yr*dnds(1) + zr*dmds(1))*temp + temp1*drdsx ! py pz
      dtsy(4,3) = (yr*dnds(2) + zr*dmds(2))*temp + temp1*drdsy ! py pz
      dtsz(4,3) = (yr*dnds(3) + zr*dmds(3))*temp + temp1*drdsz ! py pz

      dtsx(3,4) = dtsx(4,3)
      dtsy(3,4) = dtsy(4,3)
      dtsz(3,4) = dtsz(4,3)

      temp1 = zr*zr*(dVs3 - dVs4) + dVs4
      z_temp = 2.0d0*zr*temp
      dtsx(4,4) = z_temp*dnds(1) + temp1*drdsx ! pz pz
      dtsy(4,4) = z_temp*dnds(2) + temp1*drdsy ! pz pz
      dtsz(4,4) = z_temp*dnds(3) + temp1*drdsz ! pz pz
!    dtsx(4,4) = 2.0d0*zr*dnds(1)*temp + temp1*drdsx ! pz pz
!    dtsy(4,4) = 2.0d0*zr*dnds(2)*temp + temp1*drdsy ! pz pz
!    dtsz(4,4) = 2.0d0*zr*dnds(3)*temp + temp1*drdsz ! pz pz

      dtsx(1,5) = 0.0d0 ! s* s
      dtsy(1,5) = 0.0d0 ! s* s
      dtsz(1,5) = 0.0d0 ! s* s

      dtsx(5,1) = 0.0d0 ! s s*
      dtsy(5,1) = 0.0d0 ! s s*
      dtsz(5,1) = 0.0d0 ! s s*

      temp = xr*dVs5
      dtsx(2,5) = -(dlds(1)*Vs5 + temp*drdsx) ! -s* px
      dtsy(2,5) = -(dlds(2)*Vs5 + temp*drdsy) ! -s* px
      dtsz(2,5) = -(dlds(3)*Vs5 + temp*drdsz) ! -s* px

      temp = xr*dVs51
      dtsx(5,2) = (dlds(1)*Vs51 + temp*drdsx) ! px s*
      dtsy(5,2) = (dlds(2)*Vs51 + temp*drdsy) ! px s*
      dtsz(5,2) = (dlds(3)*Vs51 + temp*drdsz) ! px s*

      temp = yr*dVs5
      dtsx(3,5) = -(dmds(1)*Vs5 + temp*drdsx) ! -s* py
      dtsy(3,5) = -(dmds(2)*Vs5 + temp*drdsy) ! -s* py
      dtsz(3,5) = -(dmds(3)*Vs5 + temp*drdsz) ! -s* py

      temp = yr*dVs51
      dtsx(5,3) = (dmds(1)*Vs51 + temp*drdsx) ! s* py
      dtsy(5,3) = (dmds(2)*Vs51 + temp*drdsy) ! s* py
      dtsz(5,3) = (dmds(3)*Vs51 + temp*drdsz) ! s* py

      temp = zr*dVs5
      dtsx(4,5) = -(dnds(1)*Vs5 + temp*drdsx)  ! -s* pz
      dtsy(4,5) = -(dnds(2)*Vs5 + temp*drdsy)  ! -s* pz
      dtsz(4,5) = -(dnds(3)*Vs5 + temp*drdsz)  ! -s* pz

      temp = zr*dVs51
      dtsx(5,4) = (dnds(1)*Vs51 + temp*drdsx)  ! s* pz
      dtsy(5,4) = (dnds(2)*Vs51 + temp*drdsy)  ! s* pz
      dtsz(5,4) = (dnds(3)*Vs51 + temp*drdsz)  ! s* pz

      dtsx(5,5) = 0.0d0 ! s* s*
      dtsy(5,5) = 0.0d0 ! s* s*
      dtsz(5,5) = 0.0d0 ! s* s*

      do i1 = 1,5
         do j1 = 1,5
            if (isnan(dtsx(i1,j1))) print*, 'dHopping_s_M: dtsx', i1,j1, dtsx(i1,j1)
            if (isnan(dtsy(i1,j1))) print*, 'dHopping_s_M: dtsy', i1,j1, dtsy(i1,j1)
            if (isnan(dtsz(i1,j1))) print*, 'dHopping_s_M: dtsz', i1,j1, dtsz(i1,j1)
         enddo
      enddo
      
      nullify(KOA_i, KOA_j, x,y,z,r, xrr, yrr, zrr, x1)
      nullify(dVs1, dVs2, dVs21, dVs3, dVs4, dVs5, dVs51)
      nullify(Vs1, Vs2, Vs21, Vs3, Vs4, Vs5, Vs51)
   else ! it's all zeros:
      dtsx = 0.0d0
      dtsy = 0.0d0
      dtsz = 0.0d0
   endif
end subroutine dHopping_s_M



function dVs_M(TB,i,N1,N2,r,scaling) ! functions Vsssigma, /p sigma, /p pi
! for the Hopping Integrals according to Eqs.(2.23)-(2.29), (2.32), 
! and the Appendix E with all parameters,
! from H.Jeschke PhD thesis.
! i   ---  stands for index of V, i=1 is "s s sigma", i=2 is "s p sigma", i=3 is "p p sigma", i=4 is "p p pi"
! x, y, z --- are the coordinates of interatomic distance (positive or negative!)
! r  --- absolute interatomic distance
   type(TB_H_Molteni), dimension(:,:), intent(in), target :: TB	! all tight binding parameters
   INTEGER(4), intent(in) :: i ! index of V
   integer, intent(in) :: N1,N2 ! kinds of atoms
   real(8), intent(in) :: r  ! absolute distance between the atoms	
   real(8) dVs_M ! function itself (Eqs.(2.23)-(2.27))
   integer(4), intent(in) :: scaling ! which scaling rule to use
   !--------------------------------------------
   real(8), dimension(:), pointer :: V0
   real(8), pointer :: nc, rc
   real(8), pointer :: r0, n
   real(8), pointer :: rcut, d
   real(8) :: Vs, denom, p2, r01
   integer j

   V0 => TB(N1,N2)%V0
   n => TB(N1,N2)%n
   rcut => TB(N1,N2)%rcut
   d => TB(N1,N2)%d
   r0 => TB(N1,N2)%r0
   nc => TB(N1,N2)%nc
   rc => TB(N1,N2)%rc

   ! Unconvenional scaling:
   !dVs_M = -n*Vs(TB,i,N1,N2,r,cut_off=.false.)/r*(1.0d0 + nc*(r/rc)**nc)

   ! Harrison scaling [Molteni et al., J. Phys.: Condens. Matter 6 (1994) 5243]:
   Vs = Vs_M(TB,i,N1,N2,r,.false.) ! function above
   dVs_M = -n/r*Vs

   select case (scaling)
   case (0) ! Fermi-like smoothing (inconsistent!):
      ! Introduce cut-off function:
      dVs_M = dVs_M/(1.0d0+exp((r-rcut)*2.0d0/d)) ! cutting off long-distance atoms
   case default ! Consistent smoothing:
      denom = (1.0d0+exp((r-rcut)*2.0d0/d))
      p2 = 2.0d0/d*exp((r-rcut)*2.0d0/d)/(denom*denom)
      dVs_M = dVs_M/denom - Vs*p2 ! cutting off long-distance atoms
   end select

   if (isnan(dVs_M) .or. (abs(dVs_M) > 1d25)) print*, 'dVs_M is:', dVs_M, 1.0d0/(1.0d0+exp((r-rcut)*2.0d0/d)), r, Vs_M(TB,i,N1,N2,r,.false.), Vs, p2, i, N1, N2

   nullify( V0, nc, rc, r0, n, rcut, d)
end function dVs_M


!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! Super-cell forces:
subroutine Attract_TB_Forces_Press_M(TB_Hamil, atoms, Scell, NSC, numpar, Aij)
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   REAL(8), DIMENSION(:,:), INTENT(in) :: Aij  ! Matrix Aij, Eq. (H.3), p.153 Jeschke's PhD
   !REAL(8), DIMENSION(9,size(Aij,1)) :: dwr_press
   !REAL(8), DIMENSION(9,size(Aij,1),size(Aij,2)) :: dHij
   REAL(8), allocatable, DIMENSION(:,:) :: dwr_press
   REAL(8), allocatable, DIMENSION(:,:,:) :: dHij
   integer i, j, k, n
   if (numpar%p_const) then	! calculate this for P=const Parrinello-Rahman MD
      n = size(Aij,1)
      allocate(dwr_press(9,n))
      allocate(dHij(9,n,size(Aij,2)))
      dHij = 0.0d0
      dwr_press = 0.0d0
      call dHamil_tot_Press_M(Scell(NSC)%MDatoms, Scell, NSC, numpar, TB_Hamil, dHij)
      !$omp PARALLEL DO private(i,j)
      do i = 1, n
          do j = 1,9
            dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(i,:)) ! old, tested, good
            !dwr_press(j,i) = dwr_press(j,i) + SUM(dHij(j,i,:)*Aij(:,i)) ! wrong!!
         enddo ! j
      enddo ! i
      !$OMP END PARALLEL DO
      Scell(NSC)%SCforce%att = 0.0d0
      do i = 1,3
         do k = 1,3
            Scell(NSC)%SCforce%att(k, i) = Scell(NSC)%SCforce%att(k, i) + SUM(dwr_press((i-1)*3+k,:))
         enddo ! k
      enddo ! i
      deallocate(dwr_press)
      deallocate(dHij)
   endif
end subroutine Attract_TB_Forces_Press_M




subroutine dHamil_tot_Press_M(atoms, Scell, NSC, numpar, TB_Hamil, dHij)
! construct the whole Hamilton matrix:
! (with respect to which Rk we take the derivatives, Appendix F of H.Jeschke PhD Thesis)
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB_Hamil ! all tight binding parameters
   REAL(8), DIMENSION(:,:,:), INTENT(inout) :: dHij
   REAL(8), DIMENSION(9,5,5) :: dHij1
   integer :: i, j, j1, i1, ki, atom_2, m, nat, k, i2, j2, NumTB
   !nat = matter%Na
   nat = size(atoms)
   dHij = 0.0d0
   do i = 1,nat	! all atoms
      m = Scell(NSC)%Near_neighbor_size(i)
      do atom_2 = 1,m ! do only for atoms close to that one
         j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
         if (j .GT. 0) then
            !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(i)%KOA, Scell(NSC)%MDatoms(j)%KOA)
            call dHamilton_one_Press_M(i, j, atom_2, Scell(NSC)%MDatoms, Scell, NSC, TB_Hamil, dHij1)
            ! Eqs. (2.41), (2.42), Page 40 in H.Jeschke PhD thesis.
            do j1 = 1,5 ! all orbitals
               do i1 = 1,5 ! all orbitals
                  i2 = (i-1)*5+i1
                  j2 = (j-1)*5+j1
                  dHij(:,i2,j2) = dHij1(:,i1,j1)	! construct the total Hamiltonian from
               enddo ! i1
            enddo ! j1
         endif ! (j .GT. 0) then
      enddo ! j
   enddo ! i
end subroutine dHamil_tot_Press_M ! CHECKED


subroutine dHamilton_one_Press_M(i, j, atom_2, atoms, Scell, NSC, TB, dHij_press)
! Create a Hamiltonain-matrix, a block
! of which the total Hamiltonain is constructed
! See H.Jeschke PhD thesis, Eq.(2.40) and its description, Page 40
   integer(4), INTENT(IN) :: i, j, atom_2
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) ::  TB	! all tight binding parameters   
   REAL(8), DIMENSION(:,:,:), INTENT(out) :: dHij_press
   integer :: ki, kj
   REAL(8), DIMENSION(9,5,5) :: dts_press ! hopping integrals
   dHij_press = 0.0d0
!    j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
   if (i .EQ. j) then ! Onsite contributions, according to H.Jeschke, PhD thesis, Eq.(2.41), Page 40
      dHij_press = 0.0d0
   else
      ! For pairs of atoms, fill the hamiltonain with Hopping Integrals.
      ! Call the subroutine to calculate these hoppings: 
      call dHopping_Press_h_M(i, j, atom_2, atoms, Scell, NSC, TB, dts_press)
      ! Then fill the hamiltonain with correspongin hopping integrals, as in
      ! H.Jeschke PhD thesis, Eq.(2.42), Page 40:
      do ki = 1,5
         do kj = 1,5
            dHij_press(:,kj,ki) = dts_press(:,ki,kj)   ! Hopping Integrals, TESTED, GOOD
            !dHij_press(:,ki,kj) = dts_press(:,ki,kj)   ! Hopping Integrals, WRONG!!!
         enddo ! kj
      enddo  ! ki
   endif
end subroutine dHamilton_one_Press_M


subroutine dHopping_Press_h_M(i, j, atom_2, atoms, Scell, NSC, TB, dts_press)
! subroutine making the derivatives of the hopping integrals
   INTEGER, intent(in) :: i,atom_2 ! atoms
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_H_Molteni), dimension(:,:), intent(in) :: TB	! all tight binding parameters
   REAL(8), DIMENSION(:,:,:), INTENT(out) :: dts_press
   REAL(8), DIMENSION(9,5,5) :: dtsx
   integer :: k, ik
   real(8) :: x,y,z, x0, y0, z0, sx, sy, sz	! relative distances (projections) between the atoms 
   real(8) r, r1, diff_dircos_l, diff_dircos_m, diff_dircos_n, drijdhgd
   INTEGER(4) :: j,i1,j1,k1, dik, djk, dij
   integer KOA_i, KOA_j
   real(8), dimension(3) :: zb
   real(8) sxr, syr, szr, Vs2r, Vs2r2, Vs3r, Vs4r, Vs5r, Vs5r2, x_r, y_r, z_r
   real(8) dVs2r, dVs2r2, dVs3r, dVs4r, dVs5r, dVs5r2

   KOA_i = Scell(NSC)%MDatoms(i)%KOA ! kind of atom for atom i
   KOA_j = Scell(NSC)%MDatoms(j)%KOA ! kind of atom for atom j

   dts_press = 0.0d0

   !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
!    j = Scell(NSC)%Near_neighbor_list(i,atom_2)   ! it interacts with this atom
   x = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
   y = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
   z = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
   r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
   sx = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
   sy = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
   sz = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

   sxr = sx/r
   syr = sy/r
   szr = sz/r

   x_r = x/r
   y_r = y/r
   z_r = z/r

!    Vs2r = Vs_M(TB,2,KOA_i,KOA_j,r,.false.)
!    Vs2r2 = Vs_M(TB,2,KOA_j,KOA_i,r,.false.)
!    Vs3r = Vs_M(TB,3,KOA_i,KOA_j,r,.false.)
!    Vs4r = Vs_M(TB,4,KOA_i,KOA_j,r,.false.)
!    Vs5r = Vs_M(TB,5,KOA_i,KOA_j,r,.false.)
!    Vs5r2 = Vs_M(TB,5,KOA_j,KOA_i,r,.false.)
   Vs2r = Vs_M(TB,2,KOA_i,KOA_j,r)
   Vs2r2 = Vs_M(TB,2,KOA_j,KOA_i,r)
   Vs3r = Vs_M(TB,3,KOA_i,KOA_j,r)
   Vs4r = Vs_M(TB,4,KOA_i,KOA_j,r)
   Vs5r = Vs_M(TB,5,KOA_i,KOA_j,r)
   Vs5r2 = Vs_M(TB,5,KOA_j,KOA_i,r)

   dVs2r = dVs_M(TB,2,KOA_i,KOA_j,r,1)
   dVs2r2 = dVs_M(TB,2,KOA_j,KOA_i,r,1)
   dVs3r = dVs_M(TB,3,KOA_i,KOA_j,r,1)
   dVs4r = dVs_M(TB,4,KOA_i,KOA_j,r,1)
   dVs5r = dVs_M(TB,5,KOA_i,KOA_j,r,1)
   dVs5r2 = dVs_M(TB,5,KOA_j,KOA_i,r,1)

   ! All of the hopping integrals ts for each pair of atoms i and j,
   ! according to the Eqs.(2.23)-(2.29) from PhD Thesis of H.Jeschke, Page 38
  do k1 = 1,9 ! all the components of the h_alpha_beta(3,3):
   diff_dircos_l = 0.0d0
   diff_dircos_m = 0.0d0
   diff_dircos_n = 0.0d0
   if (k1 .EQ. 1) then ! (1,1) ax
      drijdhgd = x*sxr
      diff_dircos_l = sxr
   !endif
   else if (k1 .EQ. 2) then ! (2,1) ay
      drijdhgd = y*sxr
      diff_dircos_m = sxr
   !endif
   else if (k1 .EQ. 3) then ! (3,1) az
      drijdhgd = z*sxr
      diff_dircos_n = sxr
   !endif
   else if (k1 .EQ. 4) then ! (1,2) bx 
      drijdhgd = x*syr
      diff_dircos_l = syr
   !endif
   else if (k1 .EQ. 5) then ! (2,2) by
      drijdhgd = y*syr
      diff_dircos_m = syr
   !endif
   else if (k1 .EQ. 6) then ! (3,2) bz
      drijdhgd = z*syr
      diff_dircos_n = syr
   !endif
   else if (k1 .EQ. 7) then ! (1,3) cx
      drijdhgd = x*szr
      diff_dircos_l = szr
   !endif
   else if (k1 .EQ. 8) then ! (2,3) cy
      drijdhgd = y*szr
      diff_dircos_m = szr
   !endif
   else if (k1 .EQ. 9) then ! (3,3) cz
      drijdhgd = z*szr
      diff_dircos_n = szr
   endif

   diff_dircos_l = diff_dircos_l - x*drijdhgd/(r*r)
   diff_dircos_m = diff_dircos_m - y*drijdhgd/(r*r)
   diff_dircos_n = diff_dircos_n - z*drijdhgd/(r*r)

   dtsx(k1,1,1) = dVs_M(TB,1,KOA_i,KOA_j,r,1)*drijdhgd ! s s
   dtsx(k1,2,1) = -(diff_dircos_l*Vs2r + x_r*dVs2r*drijdhgd) ! -px s
   dtsx(k1,3,1) = -(diff_dircos_m*Vs2r + y_r*dVs2r*drijdhgd)  ! -py s
   dtsx(k1,4,1) = -(diff_dircos_n*Vs2r + z_r*dVs2r*drijdhgd)  ! -pz s
   dtsx(k1,5,1) = 0.0d0 ! s* s

   dtsx(k1,1,2) = (diff_dircos_l*Vs2r2 + x_r*dVs2r2*drijdhgd) ! s px
   dtsx(k1,2,2) = 2.0d0*x_r*(diff_dircos_l)*(Vs3r - Vs4r) + (x_r*x_r*(dVs3r - dVs4r) + dVs4r)*drijdhgd !px px
   dtsx(k1,3,2) =(x_r*(diff_dircos_m)+y_r*(diff_dircos_l))*(Vs3r-Vs4r)+x_r*y_r*(dVs3r-dVs4r)*drijdhgd ! py px
   dtsx(k1,4,2) =(x_r*(diff_dircos_n)+z_r*(diff_dircos_l))*(Vs3r-Vs4r)+x_r*z_r*(dVs3r-dVs4r)*drijdhgd ! pz px
   dtsx(k1,5,2) = (diff_dircos_l*Vs5r2 + x_r*dVs5r2*drijdhgd)  ! s* px

   dtsx(k1,1,3) = (diff_dircos_m*Vs2r2 + y_r*dVs2r2*drijdhgd)  ! s py
   dtsx(k1,2,3) = dtsx(k1,3,2)
   dtsx(k1,3,3) = 2.0d0*y_r*(diff_dircos_m)*(Vs3r-Vs4r)+(y_r*y_r*(dVs3r-dVs4r)+dVs4r)*drijdhgd ! pypy
   dtsx(k1,4,3) = (y_r*(diff_dircos_n)+z_r*(diff_dircos_m))*(Vs3r-Vs4r)+y_r*z_r*(dVs3r-dVs4r)*drijdhgd ! pypz
   dtsx(k1,5,3) = (diff_dircos_m*Vs5r2 + y_r*dVs5r2*drijdhgd)  ! s* py

   dtsx(k1,1,4) = (diff_dircos_n*Vs2r2 + z_r*dVs2r2*drijdhgd) ! s pz
   dtsx(k1,2,4) = dtsx(k1,4,2)
   dtsx(k1,3,4) = dtsx(k1,4,3)
   dtsx(k1,4,4) = 2.0d0*z_r*(diff_dircos_n)*(Vs3r-Vs4r)+(z_r*z_r*(dVs3r-dVs4r)+dVs4r)*drijdhgd ! pzpz
   dtsx(k1,5,4) = (diff_dircos_n*Vs5r2 + z_r*dVs5r2*drijdhgd)  ! s* pz

   dtsx(k1,1,5) = 0.0d0 ! s s*
   dtsx(k1,2,5) = -(diff_dircos_l*Vs5r + x_r*dVs5r*drijdhgd)  ! -px s*
   dtsx(k1,3,5) = -(diff_dircos_m*Vs5r + y_r*dVs5r*drijdhgd)  ! -py s*
   dtsx(k1,4,5) = -(diff_dircos_n*Vs5r + z_r*dVs5r*drijdhgd)  ! -pz s*
   dtsx(k1,5,5) = 0.0d0 ! s* s*

  enddo ! k1 = 1,9 ! all the components of the h_alpha_beta(3,3):
  dts_press = dtsx
end subroutine dHopping_Press_h_M



!RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
! Repulsive part of the potential energy:

subroutine get_Erep_s_M(TB_Repuls, Scell, NSC, numpar, a)   ! repulsive energy, module "TB_Molteni"
   type(Super_cell), dimension(:), intent(in) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(TB_Rep_Molteni), dimension(:,:), intent(in)   :: TB_Repuls
   type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
   real(8), intent(out) :: a
   !=====================================================
   real(8) a_r, b
   INTEGER(4) i1, j1, m, atom_2, NumTB
   real(8), DIMENSION(3) :: zb
   a = 0.0d0
   do i1 = 1, Scell(NSC)%Na
      b = 0.0d0
      m = Scell(NSC)%Near_neighbor_size(i1)
      do atom_2 = 1,m ! do only for atoms close to that one
         j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
         if (j1 .NE. i1) then
           if (j1 .GT. 0) then ! if there really is a nearest neighbor
            !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA)
            !call shortest_distance(matter, atoms, i1, j1, a_r) ! from module Atomic_tools
            a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4)  ! at this distance, R
            if (a_r <= 0.0d0) a_r = 1d25 ! exclude by nullifying the potential
            b = b + phi_M(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r)
!             b = b + phi_M(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r,.false.)
           endif ! nearest neighbor exists
         endif ! (j1 .NE. i1)
      enddo ! j1
      !a = a + fx(TB,b)
      a = a + b
!       if (isnan(b)) print*, 'b', b
!       if (isnan(a)) print*, 'a', a
   enddo ! i1
   a = a/2.0d0 ! it was doubled
end subroutine get_Erep_s_M



function phi_M(TB_Repuls,a_x,cut_off)	! repulsive potential
   real(8) :: phi_M ! the function itself
   type(TB_Rep_Molteni), intent(in) :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_x ! distance between the atoms
   logical, intent(in), optional :: cut_off ! do or don't do cut off?
   integer :: NP, m ! Np = number of potential: 1=exp (Table 3), 2=rational (Table 4)
   real(8) :: phi1, phi2, r0, alpha
   real(8) :: a, b, c, r2, r4, inv_r
   real(8) :: rcut, d

   NP = TB_Repuls%NP
   phi1 = TB_Repuls%phi1
   phi2 = TB_Repuls%phi2
   r0 = TB_Repuls%r0
   rcut = TB_Repuls%rcut
   d = TB_Repuls%d

   select case(NP)
   case (1) ! exp (Molteni):
      alpha = TB_Repuls%alpha
      phi_M = phi1*exp(-(a_x - r0)/alpha) + phi2*r0/a_x
   case (2) ! rational (Molteni):
      m = TB_Repuls%m
      phi_M = phi1*(r0/a_x)**m + phi2*r0/a_x
   case default ! Allen's potential:
      a = TB_Repuls%a ! alpha
      b = TB_Repuls%b ! beta
      c = TB_Repuls%c ! gamma
      inv_r = 1.0d0/a_x
      r2 = inv_r*inv_r ! (1/r)^2
      r4 = r2*r2       ! (1/r)^4
      phi_M = a*r4 + b*r4*r2 + c*r4*r4
   end select
!    if (isnan(phi_M)) print*, 'phi', phi_M, phi1, exp(-(a_x-r0)/alpha), a_x

   ! Introduce cut-off function:
   if (present(cut_off)) then
      if (cut_off) phi_M = phi_M/(1.0d0+exp((a_x-rcut)*2.0d0/d)) ! cutting off long-distance atoms
   else ! do cut off by default
      phi_M = phi_M/(1.0d0+exp((a_x-rcut)*2.0d0/d)) ! cutting off long-distance atoms
   endif
end function phi_M



function dphi_M(TB_Repuls,a_x,scaling)	! repulsive potential
   real(8) :: dphi_M ! the function itself
   type(TB_Rep_Molteni), intent(in) :: TB_Repuls   ! parameters of the repulsive part of TB-H
   real(8), intent(in) :: a_x ! distance between the atoms
   integer(4), intent(in) :: scaling ! which scaling law to use: with or without concistent cut off
   integer :: NP, m ! Np = number of potential: 1=exp (Table 3), 2=rational (Table 4)
   real(8) :: phi1, phi2, r0, alpha, dphi_2, denom, p2
   real(8) :: a, b, c, r2, r4, inv_r
   real(8) :: rcut, d

   NP = TB_Repuls%NP
   phi1 = TB_Repuls%phi1
   phi2 = TB_Repuls%phi2
   r0 = TB_Repuls%r0
   rcut = TB_Repuls%rcut
   d = TB_Repuls%d

   select case(NP)
   case (1) ! exp (Molteni):
      alpha = TB_Repuls%alpha
      dphi_M = -phi1*exp(-(a_x - r0)/alpha)/alpha - phi2*r0/(a_x*a_x)
   case (2) ! rational (Molteni):
      m = TB_Repuls%m
      dphi_M = -m/a_x*phi1*(r0/a_x)**m - phi2*r0/(a_x*a_x)
   case default ! Allen's potential:
      a = TB_Repuls%a ! alpha
      b = TB_Repuls%b ! beta
      c = TB_Repuls%c ! gamma
      inv_r = 1.0d0/a_x
      r2 = inv_r*inv_r ! (1/r)^2
      r4 = r2*r2       ! (1/r)^4
      dphi_M = -a*4.0d0*r4*inv_r - b*6.0d0*r4*r2*inv_r - c*8.0d0*r4*r4*inv_r
   end select
!    if (isnan(phi_M)) print*, 'phi', phi_M, phi1, exp(-(a_x-r0)/alpha), a_x

   ! Introduce cut-off function:
   select case (scaling)
   case (0) ! Fermi-like cut off (inconsistent!!! keep for comparison and historical reasons):
      dphi_M = dphi_M/(1.0d0+exp((a_x-rcut)*2.0d0/d)) ! cutting off long-distance atoms
   case default
      dphi_2 = phi_M(TB_Repuls,a_x,.false.) ! see function above
      denom = 1.0d0 + exp((a_x-rcut)*2.0d0/d)
      p2 = 2.0d0/d*exp((a_x-rcut)*2.0d0/d)/(denom*denom)
      dphi_M = dphi_M/denom - dphi_2*p2 ! cutting off long-distance atoms
   end select
end function dphi_M



subroutine dErdr_s_M(TB_Repuls, atoms, Scell, NSC, numpar) ! derivatives of the repulsive energy by s
   type(TB_Rep_Molteni), dimension(:,:), intent(in)   :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms	! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar	! numerical parameters, including lists of earest neighbors
   !type(Forces), dimension(:,:), intent(inout) :: forces1	! all interatomic forces
   REAL(8), DIMENSION(3) :: x1  ! for coordinates of all atoms (X,Y,Z)-for all atoms
   real(8) dpsi(3), psi, a_r, x, y, z, r1, x0, y0, z0, a, b, ddlta, b_delta
   integer i, j, k, ik, i1, j1, ian, dik, djk, n, m, atom_2, NumTB
   real(8), DIMENSION(3) :: zb
   real(8), dimension(:,:), allocatable :: Erx_s
   n = size(atoms)
   allocate(Erx_s(3,n)) ! x,y,z-forces for each atoms
   Erx_s = 0.0d0

   !$omp PARALLEL DO private(ian,i1,dpsi,m,j1,x,y,z,a_r,dik,djk,x1,b,ddlta,b_delta,atom_2)
   do ian = 1, n  ! Forces for all atoms
     Scell(NSC)%MDatoms(ian)%forces%rep(:) = 0.0d0 ! just to start with
     do i1 = 1, n
         dpsi = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i1)
         do atom_2 = 1,m ! do only for atoms close to that one
            j1 = Scell(NSC)%Near_neighbor_list(i1,atom_2) ! this is the list of such close atoms
            if (j1 .GT. 0) then
                if (ian .EQ. i1) then
                   dik = 1
                else
                   dik = 0
                endif
                if (ian .EQ. j1) then
                   djk = 1
                else
                   djk = 0
                endif

                if ((j1 .NE. i1) .OR. ((dik-djk) .NE. 0)) then ! without it, it gives ERROR
                  x = Scell(NSC)%Near_neighbor_dist(i1,atom_2,1) ! at this distance, X
                  y = Scell(NSC)%Near_neighbor_dist(i1,atom_2,2) ! at this distance, Y
                  z = Scell(NSC)%Near_neighbor_dist(i1,atom_2,3) ! at this distance, Z
                  a_r = Scell(NSC)%Near_neighbor_dist(i1,atom_2,4) ! at this distance, R

                  x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(1,2) + z*Scell(NSC)%supce(1,3)
                  x1(2) = x*Scell(NSC)%supce(2,1) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(2,3)
                  x1(3) = x*Scell(NSC)%supce(3,1) + y*Scell(NSC)%supce(3,2) + z*Scell(NSC)%supce(3,3)
!                   x1(1) = x*Scell(NSC)%supce(1,1) + y*Scell(NSC)%supce(2,1) + z*Scell(NSC)%supce(3,1)
!                   x1(2) = x*Scell(NSC)%supce(1,2) + y*Scell(NSC)%supce(2,2) + z*Scell(NSC)%supce(3,2)
!                   x1(3) = x*Scell(NSC)%supce(1,3) + y*Scell(NSC)%supce(2,3) + z*Scell(NSC)%supce(3,3)

                  b = dphi_M(TB_Repuls(Scell(NSC)%MDatoms(j1)%KOA, Scell(NSC)%MDatoms(i1)%KOA),a_r,1)

                  ddlta = real(dik - djk)/a_r
                  b_delta = b*ddlta
                  dpsi(1) = dpsi(1) + b_delta*x1(1) ! X, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(2) = dpsi(2) + b_delta*x1(2) ! Y, Eq.(F21), H.Jeschke PhD Thesis
                  dpsi(3) = dpsi(3) + b_delta*x1(3) ! Z, Eq.(F21), H.Jeschke PhD Thesis
               endif
            endif
         enddo ! j1

         Erx_s(1,ian) = Erx_s(1,ian) + dpsi(1) ! repulsive part in X-coordinate
         Erx_s(2,ian) = Erx_s(2,ian) + dpsi(2) ! repulsive part in Y-coordinate
         Erx_s(3,ian) = Erx_s(3,ian) + dpsi(3) ! repulsive part in Z-coordinate
      enddo ! i1
      Scell(NSC)%MDatoms(ian)%forces%rep(:) = Erx_s(:,ian)*0.5d0 ! all repulsive forces
   enddo ! ian
   !$OMP END PARALLEL DO
   deallocate(Erx_s)
END subroutine dErdr_s_M



subroutine dErdr_Pressure_s_M(TB_Repuls, atoms, Scell, NSC, numpar)! derivatives of the repulsive energy by h
   type(TB_Rep_Molteni), dimension(:,:), intent(in) :: TB_Repuls ! repulsive TB parameters
   type(Super_cell), dimension(:), intent(inout) :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   type(Atom), dimension(:), intent(in) :: atoms ! array of atoms in the supercell
   type(Numerics_param), intent(in) :: numpar ! numerical parameters, including lists of earest neighbors
   !===============================================
   REAL(8), DIMENSION(3,3) :: Rep_Pr  ! for dErep/dr for (X,Y,Z) for all atoms
   integer i,j,k,l, i1, i2, i3, ik, m, atom_2, n, NumTB
   real(8) x,y,z,sx,sy,sz,r, rcur(3), scur(3)
   real(8) df_psy, psi, dpsy

   if (numpar%p_const) then ! calculate this for P=const Parrinello-Rahman MD
      n = size(atoms)
      Scell(NSC)%SCforce%rep = 0.0d0
      do i = 1, n ! Forces from all atoms
         Rep_Pr = 0.0d0 ! to start
         dpsy = 0.0d0
         m = Scell(NSC)%Near_neighbor_size(i)
         do atom_2 = 1,m ! do only for atoms close to that one
            j = Scell(NSC)%Near_neighbor_list(i,atom_2) ! this is the list of such close atoms
            if ((i .NE. j) .AND. (j .GT. 0)) then
               !NumTB = numpar%El_num_ij(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA)
               !call shortest_distance(matter, atoms, i, j, r, x1=x, y1=y, z1=z, sx1=sx, sy1=sy, sz1=sz)
               rcur(1) = Scell(NSC)%Near_neighbor_dist(i,atom_2,1)  ! at this distance, X
               rcur(2) = Scell(NSC)%Near_neighbor_dist(i,atom_2,2)  ! at this distance, Y
               rcur(3) = Scell(NSC)%Near_neighbor_dist(i,atom_2,3)  ! at this distance, Z
               r = Scell(NSC)%Near_neighbor_dist(i,atom_2,4)  ! at this distance, R
               scur(1) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,1) ! at this distance, SX
               scur(2) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,2) ! at this distance, SY
               scur(3) = Scell(NSC)%Near_neighbor_dist_s(i,atom_2,3) ! at this distance, SZ

!              psi = psi + phi(TB_Repuls(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA),r) ! Eq.(F21), p.147 H.Jeschke PhD Thesis
!              dpsy = dphi(TB_Repuls(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA),r)
               dpsy = dphi_M(TB_Repuls(Scell(NSC)%MDatoms(j)%KOA, Scell(NSC)%MDatoms(i)%KOA),r,1)
               do k = 1,3 ! supce indices: a,b,c
                  do l = 1,3  ! supce indices: x,y,z
                     Rep_Pr(l,k) = Rep_Pr(l,k) + dpsy*rcur(k)*scur(l)/r
                  enddo ! l
               enddo ! k
            endif ! i=j
         enddo ! j

         do k = 1,3 ! supce indices
            do l = 1,3  ! supce indices
               !Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + df_psy*Rep_Pr(l,k)  ! checked for Pettifor...
               Scell(NSC)%SCforce%rep(l,k) = Scell(NSC)%SCforce%rep(l,k) + Rep_Pr(l,k)*0.5d0
            enddo ! l
         enddo ! k
      enddo ! i
   endif
end subroutine dErdr_Pressure_s_M



END MODULE TB_Molteni
