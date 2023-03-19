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
! This module contains all global variables

MODULE Variables
use Objects

 implicit none

! Program parameters:
real(8) g_as1	! duration of execution
real(8) g_time	! time step [fs]
real(8) g_dt_save	! time until saving output data [fs]
integer g_t_int	! number of time-step
integer g_ctim(8), g_c1(8) ! time stamps

! Derived datatypes (see description in module "Objects"):
type(Solid) :: g_matter           ! All material parameters
type(Error_handling) :: g_Err     ! error save
type(Pulse), dimension(:), allocatable :: g_laser       ! Laser pulse parameters
type(Numerics_param) :: g_numpar  ! all numerical parameters
type(Super_cell), dimension(:), allocatable :: g_Scell  ! supercell with all the atoms as one object

!All MC variables:
type(MC_data), dimension(:), allocatable :: g_MC        ! all MC arrays for photons, electrons and holes

! Minor and testing variables:
logical :: file_exists, file_opened, file_named  ! checking if file exists or opened well
character(200) chtest   ! file names and so on
character(3) chtest2   ! numbers so on
real(8) testvar
integer i_test, j_test, FN_test        ! for any regular counter in cycles and testing

! Temporary variables user for testing:
real(8) :: Tconf
! real(8), dimension(:,:), allocatable :: F, dF

 contains

! Deallocates all global variables:
subroutine deallocate_all()
   if (allocated(g_Scell)) deallocate(g_Scell)
   if (allocated(g_laser)) deallocate(g_laser)
   if (allocated(g_MC)) deallocate(g_MC)
   if (allocated(g_matter%Atoms)) deallocate(g_matter%Atoms)
   if (allocated(g_matter%PCF)) deallocate(g_matter%PCF)
   if (allocated(g_numpar%mask_DOS)) deallocate(g_numpar%mask_DOS)
   if (allocated(g_numpar%DOS_weights)) deallocate(g_numpar%DOS_weights)
end subroutine deallocate_all


END MODULE Variables
