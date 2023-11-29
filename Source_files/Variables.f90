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
logical :: file_exists, file_opened, file_named, logitest  ! checking if file exists or opened well
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
   ! All allocatables in the g_numpar must be reset:
   if (allocated(g_numpar%mask_DOS)) deallocate(g_numpar%mask_DOS)
   if (allocated(g_numpar%DOS_weights)) deallocate(g_numpar%DOS_weights)
   if (allocated(g_numpar%Subcell_coord_sx)) deallocate(g_numpar%Subcell_coord_sx)
   if (allocated(g_numpar%Subcell_coord_sy)) deallocate(g_numpar%Subcell_coord_sy)
   if (allocated(g_numpar%Subcell_coord_sz)) deallocate(g_numpar%Subcell_coord_sz)
   if (allocated(g_numpar%fe_input)) deallocate(g_numpar%fe_input)
   if (allocated(g_numpar%high_DOS)) deallocate(g_numpar%high_DOS)
   if (allocated(g_numpar%dt_MD_reset_grid)) deallocate(g_numpar%dt_MD_reset_grid)
   if (allocated(g_numpar%dt_MD_grid)) deallocate(g_numpar%dt_MD_grid)
   if (allocated(g_numpar%At_bath_reset_grid)) deallocate(g_numpar%At_bath_reset_grid)
   if (allocated(g_numpar%At_bath_grid_Ta)) deallocate(g_numpar%At_bath_grid_Ta)
   if (allocated(g_numpar%At_bath_grid_tau)) deallocate(g_numpar%At_bath_grid_tau)
   if (allocated(g_numpar%El_bath_reset_grid)) deallocate(g_numpar%El_bath_reset_grid)
   if (allocated(g_numpar%El_bath_grid_Ta)) deallocate(g_numpar%El_bath_grid_Ta)
   if (allocated(g_numpar%El_bath_grid_tau)) deallocate(g_numpar%El_bath_grid_tau)
   if (allocated(g_numpar%k_grid)) deallocate(g_numpar%k_grid)
end subroutine deallocate_all


END MODULE Variables
