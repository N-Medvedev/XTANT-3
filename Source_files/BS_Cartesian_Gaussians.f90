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
! This module contains subroutines to deal with Gaussian functions in Carthesian coordinates,
! such as calculating overlap integrals following Obara–Saika (OS) recurrence relations
! see e.g. http://esqc.ups-tlse.fr/11/lectures/KlopperLecture4-1x2.pdf

module BS_Cartesian_Gaussians

use Universal_constants

implicit none
PRIVATE
 
public :: Overlap_of_high_order_gaussians

 contains
 
!=============================================
subroutine Overlap_of_high_order_gaussians(Ax, Ay, Az, alpha, AMa, Norm_A, Bx, By, Bz, beta, AMb, Norm_B, Sij, Sijx, Sijy, Sijz)
   real(8), intent(in) :: Ax, Ay, Az, Bx, By, Bz	! [au] coordinates of the centers of the two gaussians
   real(8), intent(in) :: alpha, beta	! [1/au] exponential factors of the two gaussians
   integer, dimension(3), intent(in) :: AMa, AMb	! angular momenta along each direction for the two GTO
   real(8), intent(in) :: Norm_A, Norm_B	! the normalization coefficients for each direction
   real(8), intent(out) :: Sij	! overlap of the two given GTOs
   real(8), intent(out) :: Sijx, Sijy, Sijz	! overlap of the two given GTOs in each direction
   !--------------------------
   ! Since gaussians are separable, we calculate each direction individually,
   ! starting from unnormalized gaussians (to do normalization later):
   ! Along X:
   Sijx = overlap_Sij(Ax, alpha, AMa(1), Bx, beta, AMb(1)) ! see below
   ! Along Y:
   Sijy = overlap_Sij(Ay, alpha, AMa(2), By, beta, AMb(2)) ! see below
   ! Along Z:
   Sijz = overlap_Sij(Az, alpha, AMa(3), Bz, beta, AMb(3)) ! see below
   ! Total:
   Sij = Norm_A*Norm_B*Sijx*Sijy*Sijz ! all directions of gaussians together and normalized properly
end subroutine Overlap_of_high_order_gaussians


!=============================================
! Obara–Saika (OS) recurrence relations to construct overlap integrals of higher angular momenta:
! Using the following concept from http://esqc.ups-tlse.fr/11/lectures/KlopperLecture4-1x2.pdf
!              S00
!             /   \
!           |_     _|
!          S10     S01
!         /  \     /  \
!       |_    _| |_    _|
!     S20      S11     S02
!   ...        ...       ...
function overlap_Sij(RAx, alpha, AMa, RBx, beta, AMb) result(Sijx)
   real(8), intent(in) :: RAx, RBx	! [au] coordinates of the centers of the two gaussians
   real(8), intent(in) :: alpha, beta	! [1/au] exponential factors of the two gaussians
   integer, intent(in) :: AMa, AMb	! angular momenta along each direction for the two GTO
   real(8) :: Sijx
   !--------------------
   real(8), dimension(:,:), allocatable :: Sij_intermed	! for reusing intemediate elements of the recursion
   real(8) :: eta, XPA, XPB
   integer :: N_i, N_j, i, j, i_stop, j_stop, N_i_j
   N_i = AMa + 1	! indices are counted from 1, whereas AMs are counted from 0, so add 1
   N_j = AMb + 1	! indices are counted from 1, whereas AMs are counted from 0, so add 1
   allocate(Sij_intermed(N_i,N_j))
   Sij_intermed = 0.0d0 ! to start with
   
   eta = alpha + beta
   
   ! Get element S00 to start from
   call Overlap_primitive_gauss_1d(RAx, alpha, RBx, beta, Sij_intermed(1,1)) ! see below

   ! The shortest way to (N_i,N_j) is first along direction (i,0) or (0,j) (depending on N_i>=N_j or not),
   ! and then from (N_i,0) -> (N_i,N_j) {or (0,N_j) -> (N_i,N_j)}. That's what we do below:
   ! find along which direction to go first, go along it, and then go along the perpendicular direction   
   
   ALONG_I:if (N_i >= N_j) then ! the shortest way is along i first
      RECUR_I:if (N_i > 1) then ! there is something to do next in the recurrence:
         XPA = -beta/eta*(RAx - RBx)
         XPB = alpha/eta*(RAx - RBx)
   
         ! First, go along Si0 line as far as possible (and only later go inside of the triangle of recursion),
         ! to avoid calculating unnesessary elements:
         N_i_j = N_i - N_j
         i_stop = N_i_j + 1
         PRECICLE:if (i_stop >= 2) then
            do i = 2, i_stop ! go along this direction
               if (i == 2) then ! only one term of recursion
                  Sij_intermed(i,1) = recurrence_i_1j(i, 1, XPA, eta, Sij_intermed(i-1,1), 0.0d0, 0.0d0) ! below
               else ! both terms enter here
                  Sij_intermed(i,1) = recurrence_i_1j(i, 1, XPA, eta, Sij_intermed(i-1,1), Sij_intermed(i-2,1), 0.0d0) ! below
               endif
            enddo
         endif PRECICLE

         ! Then, go inside of the triangle, only if nesessary:
         SUBRECUR_J:if (N_j > 1) then
            do i = i_stop+1, N_i ! for these elements, we need to go inside of the triangle
               ! Next step along inreasing i:
               if (i == 2) then ! only one term of recursion
                  Sij_intermed(i,1) = recurrence_i_1j(i, 1, XPA, eta, Sij_intermed(i-1,1), 0.0d0, 0.0d0) ! below
               else ! both terms enter here
                  Sij_intermed(i,1) = recurrence_i_1j(i, 1, XPA, eta, Sij_intermed(i-1,1), Sij_intermed(i-2,1), 0.0d0) ! below
               endif
               do j = 2, i-N_i_j ! that's only how far we need to go inside
                  ! Next step along inreasing j:
                  if (j == 2) then ! only one term of recursion
                     Sij_intermed(i,j) = recurrence_j_1i(i, j, XPB, eta, Sij_intermed(i,j-1), 0.0d0, Sij_intermed(i-1,j-1)) ! below
                  else ! both terms enter here
                     Sij_intermed(i,j) = recurrence_j_1i(i, j, XPB, eta, Sij_intermed(i,j-1), Sij_intermed(i,j-2), Sij_intermed(i-1,j-1)) ! below
                  endif
!                   Sij_intermed(i,j) = recurrence_j_1i(i, j, XPB, eta, Sij_intermed(i,j-1), Sij_intermed(i,j-2), Sij_intermed(i-1,j-1)) ! below
               enddo
            enddo
         endif SUBRECUR_J
      endif RECUR_I
   else ALONG_I ! the shortest way is along j first
      RECUR_J:if (N_j > 1) then ! there is something to do next in the recurrence:
         XPA = -beta/eta*(RAx - RBx)
         XPB = alpha/eta*(RAx - RBx)
   
         ! First, go along Si0 line as far as possible (and only later go inside of the triangle of recursion),
         ! to avoid calculating unnesessary elements:
         N_i_j = N_j - N_i
         j_stop = N_i_j + 1
         PRECICLE2:if (j_stop >= 2) then
            do j = 2, j_stop ! go along this direction
               if (j == 2) then ! only one term of recursion
                  Sij_intermed(1,j) = recurrence_j_1i(1, j, XPB, eta, Sij_intermed(1,j-1), 0.0d0, 0.0d0) ! below
               else ! both terms enter here
                  Sij_intermed(1,j) = recurrence_j_1i(1, j, XPB, eta, Sij_intermed(1,j-1), Sij_intermed(1,j-2), 0.0d0) ! below
               endif
            enddo
         endif PRECICLE2

         ! Then, go inside of the triangle, only if nesessary:
         SUBRECUR_I:if (N_i > 1) then
            do j = j_stop+1, N_j ! for these elements, we need to go inside of the triangle
               ! Next step along inreasing i:
               if (j == 2) then ! only one term of recursion
                  Sij_intermed(1,j) = recurrence_j_1i(1, j, XPB, eta, Sij_intermed(1,j-1), 0.0d0, 0.0d0) ! below
               else ! both terms enter here
                  Sij_intermed(1,j) = recurrence_j_1i(1, j, XPB, eta, Sij_intermed(1,j-1), Sij_intermed(1,j-2), 0.0d0) ! below
               endif
               do i = 2, j-N_i_j ! that's only how far we need to go inside
                  ! Next step along inreasing j:
                  if (i == 2) then ! only one term of recursion
                     Sij_intermed(i,j) = recurrence_i_1j(i, j, XPA, eta, Sij_intermed(i-1,j), 0.0d0, Sij_intermed(i-1,j-1))  ! below
                  else ! both terms enter here
                     Sij_intermed(i,j) = recurrence_i_1j(i, j, XPA, eta, Sij_intermed(i-1,j), Sij_intermed(i-2,j), Sij_intermed(i-1,j-1))  ! below
                  endif
               enddo
            enddo
         endif SUBRECUR_I
      endif RECUR_J
   endif ALONG_I
   
   Sijx = Sij_intermed(N_i,N_j)	! reqired overlap integral obtined by recurrence
   deallocate(Sij_intermed)
end function overlap_Sij


function recurrence_i_1j(i, j, XPA, eta, Si_1j, Si_2j, Si_1j_1) result(Sijx)
   real(8), intent(in) :: XPA, eta, Si_1j, Si_2j, Si_1j_1
   integer, intent(in) :: i, j
   real(8) :: Sijx
   !===============
   select case (j)
   case (1) ! AM=0 for j=1
      if (i == 2) then ! only one term of recursion
         Sijx = XPA*Si_1j
      else ! both terms enter here
         Sijx = XPA*Si_1j + dble(i-2)/(2.0d0*eta)*Si_2j
      endif
   case default ! AM>0 for j/=1
      Sijx = XPA*Si_1j + ( dble(i-2)*Si_2j + dble(j-1)*Si_1j_1 )/(2.0d0*eta)
   endselect
end function recurrence_i_1j


function recurrence_j_1i(i, j, XPB, eta, Sij_1, Sij_2, Si_1j_1) result(Sijx)
   real(8), intent(in) :: XPB, eta, Sij_1, Sij_2, Si_1j_1
   integer, intent(in) :: i, j
   real(8) :: Sijx
   !===============
   select case (i)
   case (1) ! AM=0 for i=1
      if (j == 2) then ! only one term of recursion
         Sijx = XPB*Sij_1
      else ! both terms enter here
         Sijx = XPB*Sij_1 + dble(j-2)/(2.0d0*eta)*Sij_2
      endif
   case default ! AM>0 for i/=1
      Sijx = XPB*Sij_1 + ( dble(i-1)*Si_1j_1 + dble(j-2)*Sij_2 )/(2.0d0*eta)
   endselect
end function recurrence_j_1i
!=============================================


!=============================================
! Unnormalized overlap of two GTO as given in 
! [Szabo, Ostlund, "Modern Quantum Chemistry" 1996, p.412, Eq.(A.9)]
function Overlap_primitive_gauss_3d(RAx, RAy, RAz, alpha, RBx, RBy, RBz, beta) result(S)
   real(8) :: S	! output, ovelap for 1s-gaussian integral
   real(8), intent(in) :: RAx, RAy, RAz, RBx, RBy, RBz	! [au] coordinates of the centers of the two gaussians
   real(8), intent(in) :: alpha, beta	! [1/au] exponential factors of the two gaussians
   !----------------------------------
   real(8) :: Sx, Sy, Sz, prefactor, expfactor
   prefactor = -1.0d0 ! just to start
   expfactor = -1.0d0 ! just to start
   ! Along X:
   call Overlap_primitive_gauss_1d(RAx, alpha, RBx, beta, Sx, prefactor, expfactor) ! see below
   ! Along Y:
   call Overlap_primitive_gauss_1d(RAy, alpha, RBy, beta, Sy, prefactor, expfactor) ! see below
   ! Along Z:
   call Overlap_primitive_gauss_1d(RAz, alpha, RBz, beta, Sz, prefactor, expfactor) ! see below
   ! Total:
   S = Sx*Sy*Sz
end function Overlap_primitive_gauss_3d


! 1s overlap integral along one component:
subroutine Overlap_primitive_gauss_1d(RAx, alpha, RBx, beta, S, prefactor, expfactor)
   real(8), intent(out) :: S	! output, ovelap for 1s-gaussian integral
   real(8), intent(in) :: RAx, RBx	! [au] coordinates of the centers of the two gaussians
   real(8), intent(in) :: alpha,  beta	! [1/au] exponential factors of the two gaussians
   real(8), intent(inout), optional :: prefactor, expfactor ! if prefactors have already been precalculated, reuse them
   !----------------------------------
   real(8) :: eta, ieta, pre, expfac, R2, x

   x = RAx - RBx
   R2 = x*x
   if (R2 < 1.0d-5) then ! it is the same atom
      eta = alpha + beta
      ieta = 1.0d0/eta
      pre = DSQRT(g_Pi*ieta)
      S = pre
   else ! it is two different atoms:
      if (present(prefactor) .and. present(expfactor)) then
         if (prefactor > 0.0d0) then 
            pre = prefactor	! it was already precalculated, reuse it now
            expfac = expfactor	! it was already precalculated, reuse it now
         else
            eta = alpha + beta
            ieta = 1.0d0/eta
            pre = DSQRT(g_Pi*ieta)
            prefactor = pre	! output it to reuse later
            expfac = alpha*beta*ieta
            expfactor = expfac	! output it to reuse later
         endif
      else
         eta = alpha + beta
         ieta = 1.0d0/eta
         pre = DSQRT(g_Pi*ieta)
         expfac = alpha*beta*ieta
      endif
      S = pre*dexp(-expfac*R2)
   endif
end subroutine Overlap_primitive_gauss_1d



 
end module BS_Cartesian_Gaussians
