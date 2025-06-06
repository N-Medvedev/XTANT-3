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
! This module contains subroutines to deal with MPI calls,
! such as wrappers and MPI derived types for particular purposes
!**************************************************************


module MPI_subroutines

#ifdef MPI_USED
use mpi
#endif

use Objects

implicit none
PRIVATE  ! hides items not listed on public statement


! Interface to automatically choose from possible variables in MPI broadcasting:
interface broadcast_variable
   module procedure broadcast_variable_int
   module procedure broadcast_variable_real
   module procedure broadcast_variable_complex
   module procedure broadcast_variable_logic
   module procedure broadcast_variable_char
end interface broadcast_variable

interface broadcast_array
   module procedure broadcast_array_logic
   module procedure broadcast_array_int
   module procedure broadcast_array_real
   module procedure broadcast_array_char
   module procedure broadcast_2d_array_real
   module procedure broadcast_3d_array_real
   module procedure broadcast_4d_array_real
   module procedure broadcast_2d_array_complex
end interface broadcast_array

interface broadcast_allocatable_array
   module procedure broadcast_allocatable_char_1d_array
   module procedure broadcast_allocatable_int_1d_array      ! for integer 1d arrays
   module procedure broadcast_allocatable_int_2d_array      ! for integer 2d arrays
   module procedure broadcast_allocatable_real_1d_array     ! for read 1d arrays
   module procedure broadcast_allocatable_real_2d_array     ! for read 2d arrays
   module procedure broadcast_allocatable_real_3d_array     ! for read 3d arrays
   module procedure broadcast_allocatable_real_4d_array     ! for real 4d arrays
   module procedure broadcast_allocatable_logic_1d_array    ! for logical 1d arrays
   module procedure broadcast_allocatable_logic_3d_array    ! for logical 3d arrays
end interface broadcast_allocatable_array


interface do_MPI_Allreduce ! it uses MPI_SUM
   module procedure do_MPI_Allreduce_real_variable
   module procedure do_MPI_Allreduce_real_1d_array
   module procedure do_MPI_Allreduce_real_2d_array
   module procedure do_MPI_Allreduce_real_3d_array
   module procedure do_MPI_Allreduce_real_4d_array
   module procedure do_MPI_Allreduce_complex_2d_array
end interface do_MPI_Allreduce



interface do_MPI_Reduce ! it uses MPI_SUM
   module procedure do_MPI_Reduce_real_variable
   module procedure do_MPI_Reduce_real_1d_array
   module procedure do_MPI_Reduce_real_2d_array
   module procedure do_MPI_Reduce_real_3d_array
   module procedure do_MPI_Reduce_real_4d_array
end interface do_MPI_Reduce



public :: initialize_MPI, initialize_random_seed, Initialize_ScaLAPACK, get_MPI_lapsed_time, MPI_barrier_wrapper, MPI_fileopen_wrapper, &
            MPI_fileclose_wrapper, MPI_error_wrapper, MPI_share_Read_Input_Files, MPI_share_matter, MPI_share_numpar, &
            MPI_share_initial_configuration, MPI_share_electron_MFPs, MPI_share_photon_attenuation, MPI_share_add_data, &
            do_MPI_Reduce, do_MPI_Allreduce, broadcast_allocatable_array, MPI_share_Ritchi_CDF, broadcast_variable, broadcast_array




contains




subroutine Initialize_ScaLAPACK(MPI_param, Err)
   ! https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(Error_handling), intent(inout) :: Err   ! error save
   !----------------------
   integer :: info, nproc, nprow, npcol, myid, myrow, mycol, ctxt, ctxt_sys, ctxt_all
   integer :: i, n_cur

#ifdef MPI_USED
   ! Determine rank of calling process and the size of processor set for BLACS:
   ! https://www.ibm.com/docs/de/pessl/5.4?topic=blacs-pinfo-routine
   call BLACS_PINFO(myid, nproc) ! library ScaLAPACK (must be linked in compilation; MKL)
   ! Consistency check:
   if ( (myid /= MPI_param%process_rank) .or. (nproc /= MPI_param%size_of_cluster) ) then
      write(6,*) '[MPI process', MPI_param%process_rank, '] Problem with BLACS ranks:', myid, MPI_param%process_rank, nproc, MPI_param%size_of_cluster
   endif

   ! Get the internal default context:
   ! https://www.ibm.com/docs/en/pessl/5.4?topic=blacs-get-routine
   call BLACS_GET( -1, 0, ctxt_sys ) ! library ScaLAPACK
   MPI_param%BLACS_icontxt = ctxt_sys  ! save into the object variable

   !write(6,*) '[MPI process', MPI_param%process_rank, ']: grid 0:', ctxt, MPI_param%BLACS_icontxt, myid, myrow, mycol, MPI_param%BLACS_nprow, MPI_param%BLACS_npcol

   ! Set up a process grid of the chosen size:
   ! https://www.ibm.com/docs/en/pessl/5.5?topic=blacs-gridinit-routine
   ctxt = MPI_param%BLACS_icontxt
   ! Choice of the grid:
   MPI_param%BLACS_nprow = floor(sqrt(dble(MPI_param%size_of_cluster)))
   MPI_param%BLACS_npcol = MPI_param%size_of_cluster / MPI_param%BLACS_nprow
   n_cur = MPI_param%BLACS_nprow * MPI_param%BLACS_npcol
   i = 0 ! to start with
   do while (n_cur < MPI_param%size_of_cluster)
      i = i + 1
      n_cur = MPI_param%BLACS_nprow * (MPI_param%BLACS_npcol + i)
      if (n_cur > MPI_param%size_of_cluster) then
         i = i - 1   ! go back to the one that wasn't above the limit
         exit
      endif
   enddo
   MPI_param%BLACS_npcol = MPI_param%BLACS_npcol + i

   ! Alternative choise of the grid: 1xN:
   !MPI_param%BLACS_nprow = MPI_param%size_of_cluster
   !MPI_param%BLACS_npcol = 1

   ! Set the grid:
   call BLACS_GRIDINIT( ctxt, 'R', MPI_param%BLACS_nprow, MPI_param%BLACS_npcol )   ! library ScaLAPACK
   MPI_param%BLACS_icontxt = ctxt

   !write(6,*) '[MPI process', MPI_param%process_rank, ']: grid 1:', ctxt, MPI_param%BLACS_icontxt, myid, myrow, mycol, MPI_param%BLACS_nprow, MPI_param%BLACS_npcol

   MPI_param%BLACS_myrow = -1 ! to start with
   MPI_param%BLACS_mycol = -1 ! to start with

   ! Processes not belonging to ctxt have nothing to do
   if (ctxt < 0) return

   ! Get the process coordinates in the grid
   call BLACS_GRIDINFO( ctxt, MPI_param%BLACS_nprow, MPI_param%BLACS_npcol, myrow, mycol )   ! library ScaLAPACK
   MPI_param%BLACS_myrow = myrow
   MPI_param%BLACS_mycol = mycol

   !write(6,*) '[MPI process', MPI_param%process_rank, ']: grid 2:', ctxt, MPI_param%BLACS_icontxt, myid, myrow, mycol, MPI_param%BLACS_nprow, MPI_param%BLACS_npcol
   !  now ScaLAPACK or PBLAS procedures can be used
#endif
!pause 'Initialize_ScaLAPACK'
end subroutine Initialize_ScaLAPACK



subroutine MPI_share_photon_attenuation(matter, numpar, Err)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-----------------------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   character(100) :: error_part
   integer :: Nshl, N_grid, Nsiz, i, j

#ifdef MPI_USED
   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------

   ! First of, allocate the TB parameterization array:
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_photon_attenuation' ! part of the error message

   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Err%Err}', Err%Err) ! below

   ! Allocate arrays for non-master processes:
   if (MPI_param%process_rank == 0) then   ! MPI master processes
      N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
   else
      N_grid = 0   ! to start with
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {N_grid}', N_grid) ! below

   if (MPI_param%process_rank /= 0) then   ! MPI non-master processes
      ! Total mean free paths:
      if (.not.allocated(matter%Ph_MFP_tot%L)) allocate(matter%Ph_MFP_tot%L(N_grid))
      if (.not.allocated(matter%Ph_MFP_tot%E)) allocate(matter%Ph_MFP_tot%E(N_grid))

      do i = 1, size(matter%Atoms) ! for all atoms
         Nshl = size(matter%Atoms(i)%Ip)
         do j = 1, Nshl ! for all shells of this atom
            if (.not.allocated(matter%Atoms(i)%Ph_MFP(j)%E)) allocate(matter%Atoms(i)%Ph_MFP(j)%E(N_grid))
            if (.not.allocated(matter%Atoms(i)%Ph_MFP(j)%L)) allocate(matter%Atoms(i)%Ph_MFP(j)%L(N_grid))
         enddo
      enddo
   endif ! (MPI_param%process_rank /= 0)
   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------

   ! Once allocated, they are ready to accept the broadcasting variables:
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Ph_MFP_tot%E}', matter%Ph_MFP_tot%E) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Ph_MFP_tot%L}', matter%Ph_MFP_tot%L) ! below

   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl ! for all shells of this atom
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Ph_MFP(j)%E}', matter%Atoms(i)%Ph_MFP(j)%E) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Ph_MFP(j)%L}', matter%Atoms(i)%Ph_MFP(j)%L) ! below

         !print*, '[MPI process #', MPI_param%process_rank, '] Ph:', i, matter%Atoms(i)%Ph_MFP(j)%L(:)
      enddo
   enddo

   nullify(MPI_param)
#endif
end subroutine MPI_share_photon_attenuation




subroutine MPI_share_electron_MFPs(matter, numpar, Err)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-----------------------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   character(100) :: error_part
   integer :: N_Te_points, N_grid, Nsiz, i, j, k

#ifdef MPI_USED
   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------

   ! First of, allocate the TB parameterization array:
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_electron_MFPs' ! part of the error message
   !print*, '[MPI process #', MPI_param%process_rank, '] MPI_share_electron_MFPs test 1'


   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Err%Err}', Err%Err) ! below

   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%hw_plasma}', matter%hw_plasma) ! below

   ! Allocate all MFP arrays in non-maaster processes:
   if (MPI_param%process_rank == 0) then   ! MPI non-master processes
      N_Te_points = size(matter%Atoms(1)%El_MFP_vs_T)
   else
      N_Te_points = 0   ! to start with
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {N_Te_points}', N_Te_points) ! below

   if (MPI_param%process_rank == 0) then   ! MPI non-master processes
      N_grid = size(matter%El_MFP_tot%L)
   else
      N_grid = 0   ! to start with
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {N_grid}', N_grid) ! below
   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------

   if (MPI_param%process_rank /= 0) then   ! MPI non-master processes
      ! Total mean free paths:
      if (.not.allocated(matter%El_MFP_tot%L)) allocate(matter%El_MFP_tot%L(N_grid))
      if (.not.allocated(matter%El_MFP_tot%E)) allocate(matter%El_MFP_tot%E(N_grid))
      if (.not.allocated(matter%El_EMFP_tot%L)) allocate(matter%El_EMFP_tot%L(N_grid))
      if (.not.allocated(matter%El_EMFP_tot%E)) allocate(matter%El_EMFP_tot%E(N_grid))
   endif

   ! For atom and shell:
   do i = 1, size(matter%Atoms)
      if (MPI_param%process_rank == 0) then   ! MPI non-master processes
         Nsiz = size(matter%Atoms(i)%El_MFP)
      else
         Nsiz = 0   ! to start with
      endif
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Nsiz}', Nsiz) ! below

      if (MPI_param%process_rank /= 0) then   ! MPI non-master processes
         if (.not.allocated(matter%Atoms(i)%El_MFP)) allocate(matter%Atoms(i)%El_MFP(Nsiz))
         if (.not.allocated(matter%Atoms(i)%El_EMFP%E)) allocate(matter%Atoms(i)%El_EMFP%E(N_grid))
         if (.not.allocated(matter%Atoms(i)%El_EMFP%L)) allocate(matter%Atoms(i)%El_EMFP%L(N_grid))

         do k = 1, Nsiz
            if (.not.allocated(matter%Atoms(i)%El_MFP(k)%E)) allocate(matter%Atoms(i)%El_MFP(k)%E(N_grid))
            if (.not.allocated(matter%Atoms(i)%El_MFP(k)%L)) allocate(matter%Atoms(i)%El_MFP(k)%L(N_grid))
         enddo ! k = 1, Nsiz

         ! Temperature dependence (for the valence/conduction band only):
         if (.not. allocated(matter%Atoms(i)%El_MFP_vs_T)) then
            allocate(matter%Atoms(i)%El_MFP_vs_T(N_Te_points))
            do j = 1, N_Te_points
               allocate(matter%Atoms(i)%El_MFP_vs_T(j)%L(N_grid))
               allocate(matter%Atoms(i)%El_MFP_vs_T(j)%E(N_grid))
            enddo ! j = 1, N_Te_points
         endif
      endif

   enddo ! i = 1, size(matter%Atoms)
   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------

   ! Once allocated, they are ready to accept the broadcasting variables:
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%El_MFP_tot%E}', matter%El_MFP_tot%E) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%El_MFP_tot%L}', matter%El_MFP_tot%L) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%El_EMFP_tot%E}', matter%El_EMFP_tot%E) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%El_EMFP_tot%L}', matter%El_EMFP_tot%L) ! below
   do i = 1, size(matter%Atoms)
      call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%El_EMFP%E}', matter%Atoms(i)%El_EMFP%E) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%El_EMFP%L}', matter%Atoms(i)%El_EMFP%L) ! below
      Nsiz = size(matter%Atoms(i)%El_MFP)
      do k = 1, Nsiz
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%El_MFP(k)%E}', matter%Atoms(i)%El_MFP(k)%E) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%El_MFP(k)%L}', matter%Atoms(i)%El_MFP(k)%L) ! below
         !print*, i, k, matter%Atoms(1)%El_MFP
      enddo ! k = 1, Nsiz
      if (i == 1) then ! CS vs Te is only defined for the VB, represented as the first element (others are undefined):
         do j = 1, N_Te_points
            call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%El_MFP_vs_T(j)%E}', matter%Atoms(i)%El_MFP_vs_T(j)%E) ! below
            call broadcast_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%El_MFP_vs_T(j)%L}', matter%Atoms(i)%El_MFP_vs_T(j)%L) ! below
         enddo ! j = 1, N_Te_points
      endif
   enddo ! i = 1, size(matter%Atoms)

   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------
   !pause 'MPI_share_electron_MFPs'
   nullify(MPI_param)
#endif
end subroutine MPI_share_electron_MFPs




subroutine MPI_share_initial_configuration(Scell, matter, numpar, laser, MC, Err)   ! module "MPI_subroutines"
   type(Super_cell), dimension(:), intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(MC_data), dimension(:), allocatable, intent(inout) :: MC ! all MC parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !-----------------------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   character(100) :: error_part
   integer :: Nsiz, N_arr_siz(3), i, j, n1
   logical :: array_is_allocated, array_is_allocated2

#ifdef MPI_USED
   ! First of, allocate the TB parameterization array:
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_initial_configuration' ! part of the error message

   !-----------------------------------------
   ! If BOP creating potential is used:
   ! If file with BOP repulsive potential does not exist, create it:
   if (numpar%create_BOP_repulse) then
      do i = 1, size(Scell(1)%TB_Repuls,1)
         do j = 1, size(Scell(1)%TB_Repuls,2)
            ASSOCIATE (TB_Repuls => Scell(1)%TB_Repuls)
               select type(TB_Repuls)
               type is (TB_Rep_BOP)
                  ASSOCIATE (TB_Hamil => Scell(1)%TB_Hamil)
                     select type(TB_Hamil)
                     type is (TB_H_BOP)   ! it can be various basis sets:
                        call MPI_share_BOP_TB_Params(MPI_param, TB_Hamil, TB_Repuls)   ! below
                     endselect
                  END ASSOCIATE
               endselect
            END ASSOCIATE
         enddo ! j
      enddo ! i
   endif

   !-----------------------------------------
   ! Share updated numpar parameters:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{numpar%do_path_coordinate}', numpar%do_path_coordinate) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{numpar%save_files_used}', numpar%save_files_used) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{numpar%vel_from_file}', numpar%vel_from_file) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{numpar%fe_filename}', numpar%fe_filename) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{numpar%fe_input_exists}', numpar%fe_input_exists) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{numpar%mask_DOS}', numpar%mask_DOS) ! below

   ! Updated matter parameters:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{matter%W_PR}', matter%W_PR) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{matter%dens}', matter%dens) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{matter%At_dens}', matter%At_dens) ! below

   ! Update laser parameters:
   do j = 1,size(laser) ! for each pulse:
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{laser(j)%Fabs}', laser(j)%Fabs) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{laser(j)%Nph}', laser(j)%Nph) ! below
   enddo

   ! Supercell parameters:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%V}', Scell(1)%V) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%supce}', Scell(1)%supce) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%supce0}', Scell(1)%supce0) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Vsupce}', Scell(1)%Vsupce) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Vsupce0}', Scell(1)%Vsupce0) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%supce_eq}', Scell(1)%supce_eq) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%k_supce}', Scell(1)%k_supce) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%supce_t}', Scell(1)%supce_t) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%GG}', Scell(1)%GG) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%SCforce%total}', Scell(1)%SCforce%total) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%SCforce%att}', Scell(1)%SCforce%att) ! below
   call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%SCforce%rep}', Scell(1)%SCforce%rep) ! below

   ! Atomic and electronic supercell parameters:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Na}', Scell(1)%Na) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ne}', Scell(1)%Ne) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ne_low}', Scell(1)%Ne_low) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ne_high}', Scell(1)%Ne_high) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ne_emit}', Scell(1)%Ne_emit) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ta}', Scell(1)%Ta) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%TaeV}', Scell(1)%TaeV) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Te}', Scell(1)%Te) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%TeeV}', Scell(1)%TeeV) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Nph}', Scell(1)%Nph) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Nh}', Scell(1)%Nh) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%G_ei_partial}', Scell(1)%G_ei_partial) ! below

   if (MPI_param%process_rank /= 0) then   ! MPI non-master processes
      if (.not.allocated(Scell(1)%MDAtoms)) allocate(Scell(1)%MDAtoms(Scell(1)%Na))
   endif

   do i = 1, Scell(1)%Na   ! for all atoms
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%KOA}', Scell(1)%MDAtoms(i)%KOA) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%Ekin}', Scell(1)%MDAtoms(i)%Ekin) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%Epot}', Scell(1)%MDAtoms(i)%Epot) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%q}', Scell(1)%MDAtoms(i)%q) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%Nx_subcel}', Scell(1)%MDAtoms(i)%Nx_subcel) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%Ny_subcel}', Scell(1)%MDAtoms(i)%Ny_subcel) ! below
      call broadcast_variable(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%Nz_subcel}', Scell(1)%MDAtoms(i)%Nz_subcel) ! below

      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%R}', Scell(1)%MDAtoms(i)%R) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%S}', Scell(1)%MDAtoms(i)%S) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%V}', Scell(1)%MDAtoms(i)%V) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%SV}', Scell(1)%MDAtoms(i)%SV) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%A}', Scell(1)%MDAtoms(i)%A) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%R0}', Scell(1)%MDAtoms(i)%R0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%S0}', Scell(1)%MDAtoms(i)%S0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%V0}', Scell(1)%MDAtoms(i)%V0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%SV0}', Scell(1)%MDAtoms(i)%SV0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%A0}', Scell(1)%MDAtoms(i)%A0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%R_eq}', Scell(1)%MDAtoms(i)%R_eq) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%S_eq}', Scell(1)%MDAtoms(i)%S_eq) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%A_tild}', Scell(1)%MDAtoms(i)%A_tild) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%v_F}', Scell(1)%MDAtoms(i)%v_F) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%v_J}', Scell(1)%MDAtoms(i)%v_J) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%A_tild0}', Scell(1)%MDAtoms(i)%A_tild0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%v_F0}', Scell(1)%MDAtoms(i)%v_F0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%v_J0}', Scell(1)%MDAtoms(i)%v_J0) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%accel}', Scell(1)%MDAtoms(i)%accel) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%accel0}', Scell(1)%MDAtoms(i)%accel0) ! below

      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%forces%rep}', Scell(1)%MDAtoms(i)%forces%rep) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%forces%att}', Scell(1)%MDAtoms(i)%forces%att) ! below
      call broadcast_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%MDAtoms(i)%forces%total}', Scell(1)%MDAtoms(i)%forces%total) ! below

      !print*, '[MPI process #', MPI_param%process_rank, ']', i, Scell(1)%MDAtoms(i)%R, Scell(1)%MDAtoms(i)%A_tild
   enddo


   ! Nearest neighbors lists:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Near_neighbor_list}', Scell(1)%Near_neighbor_list) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Near_neighbor_dist}', Scell(1)%Near_neighbor_dist) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Near_neighbor_dist_s}', Scell(1)%Near_neighbor_dist_s) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Near_neighbor_size}', Scell(1)%Near_neighbor_size) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Near_neighbors_user}', Scell(1)%Near_neighbors_user) ! below

   ! Hamiltonian arrays to allocate:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ha}', Scell(1)%Ha) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Sij}', Scell(1)%Sij) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Hij}', Scell(1)%Hij) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Hij_sol}', Scell(1)%Hij_sol) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ha0}', Scell(1)%Ha0) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%H_non}', Scell(1)%H_non) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%H_non0}', Scell(1)%H_non0) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ei}', Scell(1)%Ei) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ei0}', Scell(1)%Ei0) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Mij}', Scell(1)%Mij) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Aij}', Scell(1)%Aij) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%Ei_scc_part}', Scell(1)%Ei_scc_part) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%eigen_S}', Scell(1)%eigen_S) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%fe}', Scell(1)%fe) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//'{Scell(1)%fe_eq}', Scell(1)%fe_eq) ! below


   ! Share MC data:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(Scell(1)%MChole)) then
         array_is_allocated = .true.
         Nsiz = size(Scell(1)%MChole)
      else
         array_is_allocated = .false.
         Nsiz = 0
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Nsiz}', Nsiz) ! below

   if (array_is_allocated) then  ! broadcast the TB parameterization
      if (MPI_param%process_rank /= 0) then   ! MPI non-master processes
         allocate(Scell(1)%MChole(Nsiz))

         do i = 1, Nsiz ! for each kind of atoms:
            allocate(Scell(1)%MChole(i)%Noh(matter%Atoms(i)%sh))
            Scell(1)%MChole(i)%Noh(:) = 0.0d0 ! no holes in any shell
         enddo
         if (numpar%NMC > 0) then
            allocate(MC(numpar%NMC))	! all MC arrays for photons, electrons and holes

            do i = 1, size(MC)
               MC(i)%noe = 0.0d0
               MC(i)%noe_emit = 0.0d0
               MC(i)%noh_tot = 0.0d0
               allocate(MC(i)%electrons(Scell(1)%Ne))
               MC(i)%electrons(:)%E = 0.0d0
               MC(i)%electrons(:)%ti = 1d25
               MC(i)%electrons(:)%colls = 0
               allocate(MC(i)%holes(Scell(1)%Ne))
               MC(i)%holes(:)%E = 0.0d0
               MC(i)%holes(:)%ti = 1d26
            enddo
         endif !(size(MC) > 0)
      endif
   endif

   !print*, '[MPI process #', MPI_param%process_rank, '] Ha:', numpar%NMC, allocated(MC), size(MC)

   !-----------------------------------------
   ! Synchronize all processes:
   call MPI_barrier_wrapper(numpar%MPI_param)   ! below
   !-----------------------------------------
   nullify(MPI_param)

   !pause 'MPI_share_initial_configuration'
#endif
end subroutine MPI_share_initial_configuration





subroutine MPI_share_Read_Input_Files(matter, numpar, laser, Scell, Err)
   type(Solid), intent(inout) :: matter	! all material parameters
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   type(Error_handling), intent(inout) :: Err	! error save
   !-------------------------------------

#ifdef MPI_USED

   !--------------------------------------------------------------------------
   ! Before going any further, check if master process detected any error:
   call broadcast_variable(numpar%MPI_param, 'MPI_share_Read_Input_Files{Err%Err}', Err%Err) ! below
   if (Err%Err) return   ! if there was an error in the input files, cannot continue, go to the end...

   call broadcast_variable(numpar%MPI_param, 'MPI_share_Read_Input_Files{Err%Stopsignal}', Err%Stopsignal) ! below
   if (Err%Stopsignal) return     ! if the USER does not want to run the calculations, stop

   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"
   !print*, '[MPI process #', numpar%MPI_param%process_rank, '] Read_Input:0'

   !--------------------------------------------------------------------------
   ! Master process shares entire matter:
   call MPI_share_matter(numpar, matter) ! below
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"
   !print*, '[MPI process #', numpar%MPI_param%process_rank, '] Read_Input:1'

   !--------------------------------------------------------------------------
   ! Master process shares entire numpar:
   call MPI_share_numpar(numpar) ! below
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"
   !print*, '[MPI process #', numpar%MPI_param%process_rank, '] Read_Input:2'

   !--------------------------------------------------------------------------
   ! Master process shares entire Scell:
   call MPI_share_Scell(numpar, Scell) ! below
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"
   !print*, '[MPI process #', numpar%MPI_param%process_rank, '] Read_Input:3'

   !--------------------------------------------------------------------------
   ! Master process shares entire laser:
   call MPI_share_laser(numpar, laser) ! below
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"
   !print*, '[MPI process #', numpar%MPI_param%process_rank, '] Read_Input:4'

   !--------------------------------------------------------------------------
   ! Master process shares TB parameters:
   ! For now, assume a single supercell (as it is throughout most of he code...):
   call MPI_share_TB_parameters(matter, numpar, Scell(1)%TB_Repuls, Scell(1)%TB_Hamil, &
                                Scell(1)%TB_Waals, Scell(1)%TB_Coul, Scell(1)%TB_Expwall) ! above
   call MPI_barrier_wrapper(numpar%MPI_param)  ! module "MPI_subroutines"
  !print*, '[MPI process #', numpar%MPI_param%process_rank, '] Read_Input:5'

#endif
end subroutine MPI_share_Read_Input_Files




subroutine MPI_share_TB_parameters(matter, numpar, TB_Repuls, TB_Hamil, TB_Waals, TB_Coul, TB_Expwall)
   type(Solid), intent(in) :: matter	! all material parameters
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   class(TB_repulsive), dimension(:,:), allocatable, intent(inout) :: TB_Repuls   ! parameters of the repulsive part of TB
   class(TB_Hamiltonian), dimension(:,:), allocatable, intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   class(TB_vdW),  dimension(:,:), allocatable, intent(inout) :: TB_Waals         ! parameters of the van der Waals for TB
   class(TB_Coulomb),  dimension(:,:), allocatable, intent(inout) :: TB_Coul	! parameters of the Coulomb together with TB
   class(TB_Exp_wall),  dimension(:,:), allocatable, intent(inout) :: TB_Expwall	! parameters of the exponential wall with TB
   !--------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   integer :: Nsiz, N_arr_siz(3), i, j
   character(100) :: error_part, TB_param_name, TB_repulse_name, TB_Waals_name, TB_Coul_name, TB_Expwall_name
   logical :: array_is_allocated, array_is_allocated2


#ifdef MPI_USED

   ! First of, allocate the TB parameterization array:
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_TB_parameters' ! part of the error message
   TB_param_name = ''   ! to start with
   TB_repulse_name = ''   ! to start with
   TB_Waals_name = ''   ! to start with
   TB_Coul_name = ''   ! to start with
   TB_Expwall_name = ''   ! to start with



   !===============================================================
   ! 1) TB Hamiltonian part:
   ! Make sure it is defined:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(TB_Hamil)) then
         array_is_allocated = .true.
         N_arr_siz(1) = size(TB_Hamil,1)
         N_arr_siz(2) = size(TB_Hamil,2)
         N_arr_siz(3) = 0
         TB_param_name = TB_Hamil(1,1)%Param
      else
         array_is_allocated = .false.
         N_arr_siz = 0
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_param_name}', TB_param_name) ! below

   if (array_is_allocated) then  ! broadcast the TB parameterization
      ! Make the TB parameters of a selected class, depending on what is read in the file:
      select case (trim(adjustl(TB_param_name)))
      case ('Pettifor')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_Pettifor::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Pettifor parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
      case ('Molteni')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_Molteni::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Molteni parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
      case ('Fu')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_Fu::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Fu parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
      case ('Mehl', 'NRL')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_NRL::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for Mehl parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
      case ('DFTB')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_DFTB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
         ! DFTB skf files contain parameters for both Hamiltonian and Repulsive potential, allocate both of them here:
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_DFTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_param_name))
      case ('DFTB_no_repulsion', 'DFTB_no_repulsive', 'DFTB_no_rep')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_DFTB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
         ! DFTB skf files does not contain Repulsive potential, allocate special case:
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_DFTB_no::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_param_name))
      case ('3TB')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_3TB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for 3TB parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
      case ('BOP')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_BOP::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for BOP parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
         ! BOP files contain parameters for both Hamiltonian and Repulsive potential, allocate both of them here:
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_BOP::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for BOP parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_param_name))
      case ('xTB', 'GFN', 'GFN0')
         if (.not.allocated(TB_Hamil)) then
            allocate(TB_H_xTB::TB_Hamil(matter%N_KAO,matter%N_KAO)) ! make it for xTB parametrization
            TB_Hamil%Param = ''
         endif
         TB_Hamil(:,:)%Param = trim(adjustl(TB_param_name))
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_xTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for xTB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_param_name))
      case default
         print*, '[MPI process #', MPI_param%process_rank, '] Wrong TB-Hamiltonian parametrization class ' &
                  //trim(adjustl(TB_param_name))//' encounterred in MPI_share_TB_parameters subroutine'
         return
      end select
   endif ! array_is_allocated


   ! Now, once it is allocated, get the parameters of the Hamiltonian:
   select type (TB_Hamil)
   type is (TB_H_Pettifor)
      call MPI_share_Pettifor_TB_Hamiltonian(MPI_param, TB_Hamil)

   type is (TB_H_Molteni)
      call MPI_share_Molteni_TB_Hamiltonian(MPI_param, TB_Hamil)

   type is (TB_H_Fu)
      call MPI_share_Fu_TB_Hamiltonian(MPI_param, TB_Hamil)

   type is (TB_H_NRL)
      call MPI_share_Mehl_TB_Hamiltonian(MPI_param, TB_Hamil)

   type is (TB_H_DFTB) !in this case, read both Hamiltonian and Repulsive parts together:
      select type (TB_Repuls)  ! to confirm that repulsive part is consistent with the Hamiltonian
      type is (TB_Rep_DFTB)   ! repulsive parameters provided in the skf-file
         call MPI_share_DFTB_TB_Params(MPI_param, TB_Repuls, TB_Hamil)
      type is (TB_Rep_DFTB_no)   ! no repulsive parameters in skf-file
         call MPI_share_DFTB_TB_Params_no_rep(MPI_param, TB_Repuls, TB_Hamil)
      endselect

   type is (TB_H_3TB)
      call MPI_share_3TB_TB_Params(MPI_param, TB_Hamil)

   type is (TB_H_BOP) !in this case, read both Hamiltonian and Repulsive parts together:
      select type (TB_Repuls)  ! to confirm that repulsive part is consistent with the Hamiltonian
      type is (TB_Rep_BOP)
         call MPI_share_BOP_TB_Params(MPI_param, TB_Hamil, TB_Repuls)
      endselect

   type is (TB_H_xTB) !in this case, read both Hamiltonian and Repulsive parts together:
      select type (TB_Repuls)  ! to confirm that repulsive part is consistent with the Hamiltonian
      type is (TB_Rep_xTB)
         call MPI_share_xTB_Params(MPI_param, TB_Hamil, TB_Repuls)
      endselect
   end select

   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%basis_size_ind}', numpar%basis_size_ind) ! below



   !===============================================================
   ! 2) TB Repulsive part:
   ! Make the TB parameters of a selected class, depending on what is read in the file:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(TB_Repuls)) then
         array_is_allocated = .true.
         TB_repulse_name = TB_Repuls(1,1)%Param
      else
         array_is_allocated = .false.
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated#2}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_repulse_name}', TB_repulse_name) ! below

   if (array_is_allocated) then
      select case (trim(adjustl(TB_repulse_name)))
      case ('Pettifor')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_Pettifor::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for Pettifor parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('Molteni')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_Molteni::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for Molteni parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('Fu')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_Fu::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for Fu parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('Mehl', 'NRL')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_NRL::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for NRL parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('DFTB')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_DFTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('DFTB_no_repulsion', 'DFTB_no_repulsive', 'DFTB_no_rep')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_DFTB_no::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for DFTB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('3TB')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_3TB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for 3TB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('BOP')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_BOP::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for BOP parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case ('xTB')
         if (.not.allocated(TB_Repuls)) then
            allocate(TB_Rep_xTB::TB_Repuls(matter%N_KAO,matter%N_KAO)) ! make it for xTB parametrization
            TB_Repuls%Param = ''
         endif
         TB_Repuls(:,:)%Param = trim(adjustl(TB_repulse_name))
      case default
         print*, '[MPI process #', MPI_param%process_rank, '] Wrong TB-repulsive parametrization class ' &
               //trim(adjustl(TB_param_name))//' encounterred in MPI_share_TB_parameters subroutine'
         return
      end select
   endif

   ! Prior to use TB parameters, we now always have to find out which class they belong to:
   select type (TB_Repuls)
   type is (TB_Rep_Pettifor)
      call MPI_share_Pettifor_TB_repulsive(MPI_param, TB_Repuls)

   type is (TB_Rep_Molteni)
      call MPI_share_Molteni_TB_repulsive(MPI_param, TB_Repuls)

   type is (TB_Rep_Fu)
      call read_Fu_TB_repulsive(MPI_param, TB_Repuls)

   type is (TB_Rep_NRL)
      ! There is no repulsive part in NRL

   type is (TB_Rep_DFTB)
      call MPI_share_DFTB_TB_repulsive(MPI_param, TB_Repuls)

   type is (TB_Rep_DFTB_no)
      ! Nothing to read since repulsive potential is not provided in skf-file

   type is (TB_Rep_3TB)
      ! There is no repulsive part in 3TB

   type is (TB_Rep_BOP)
      ! Nothing to do yet with repulsive part in BOP parameterization

   type is (TB_Rep_xTB)
      ! Repulsive part has already been read above together with Hamiltonian parameters
   end select



   !===============================================================
   ! 3) Classical potential part:
   ! Make the TB vdW parameters of a selected class, depending on what is read in the file:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(TB_Waals)) then
         array_is_allocated = .true.
         TB_Waals_name = TB_Waals(1,1)%Param
      else
         array_is_allocated = .false.
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated#3}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals_name}', TB_Waals_name) ! below


   if (array_is_allocated) then
      select case (trim(adjustl(TB_Waals_name)))
      case ('Girifalco')
         if (.not.allocated(TB_Waals)) then
            allocate(TB_vdW_Girifalco::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for Girifalco parametrization
         endif
         TB_Waals(:,:)%Param = trim(adjustl(TB_Waals_name))

      case('LJ', 'Lennar-Jones', 'Lennard_Jones', 'lj', 'Lj', 'LENNARD-JONES', 'LENNARD_JONES', 'lennard-jones', 'lennard_jones')
         if (.not.allocated(TB_Waals)) then
            allocate(TB_vdW_LJ_cut::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for LJ parametrization
         endif
         TB_Waals(:,:)%Param = trim(adjustl(TB_Waals_name))

      case('ILJ', 'ilj', 'Improved_LJ', 'Ilj', 'improved_lj', 'I_LJ', 'i_lj')
         if (.not.allocated(TB_Waals)) then
            allocate(TB_vdW_ILJ_cut::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for LJ parametrization
         endif
         TB_Waals(:,:)%Param = trim(adjustl(TB_Waals_name))

      case ('Dumitrica') ! UNFINISHED, DO NOT USE
         if (.not.allocated(TB_Waals)) then
            allocate(TB_vdW_Dumitrica::TB_Waals(matter%N_KAO,matter%N_KAO)) ! make it for Dumitrica parametrization
         endif
         TB_Waals(:,:)%Param = trim(adjustl(TB_Waals_name))
      case default
         ! ignore unknown parameterization
      end select
   endif

   ! Read the parameters of the dispersion correction (vdW-type additional potential):
   select type (TB_Waals)
   type is (TB_vdW_Girifalco)
      call MPI_share_vdW_Girifalco_TB(MPI_param, TB_Waals)

   type is (TB_vdW_LJ_cut) ! Lennard-Jones (smoothly cut at short and large sitances)
      call MPI_share_vdW_LJ_TB(MPI_param, TB_Waals)

   type is (TB_vdW_ILJ_cut) ! Improved Lennard-Jones (smoothly cut at short and large sitances)
      call MPI_share_vdW_ILJ_TB(MPI_param, TB_Waals)

   type is (TB_vdW_Dumitrica) ! UNFINISHED, DO NOT USE
      call MPI_share_vdW_Dumitrica_TB(MPI_param, TB_Waals)
   end select


   !===============================================================
   ! 4) Coulomb potential part:
   ! Make the Coulomb parameters of a selected class, depending on what is read in the file:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(TB_Coul)) then
         array_is_allocated = .true.
         TB_Coul_name = TB_Coul(1,1)%Param
      else
         array_is_allocated = .false.
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated#4}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Coul_name}', TB_Coul_name) ! below

   if (array_is_allocated) then
      select case (trim(adjustl(TB_Coul_name)))
      case ('Coulomb_cut')
         if (.not.allocated(TB_Coul)) then
            allocate(TB_Coulomb_cut::TB_Coul(matter%N_KAO,matter%N_KAO)) ! make it for Coulomb parametrization
            TB_Coul%Param = ''
         endif
         TB_Coul(:,:)%Param = trim(adjustl(TB_Coul_name))
      case ('Cutie') ! testing ONLY
         if (.not.allocated(TB_Coul)) then
            allocate(Cutie::TB_Coul(matter%N_KAO,matter%N_KAO)) ! make it for Coulomb parametrization
            TB_Coul%Param = ''
         endif
         TB_Coul(:,:)%Param = trim(adjustl(TB_Coul_name))
      case default
         ! ignore unknown parameterization
      end select
   endif

   select type (TB_Coul)
   type is (TB_Coulomb_cut)
      call MPI_share_Coulomb_cut_TB(MPI_param, TB_Coul)
   end select



   !===============================================================
   ! 5) Additional short-range potential:
   ! Make the exponential wall  parameters of a selected class, depending on what is read in the file:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(TB_Expwall)) then
         array_is_allocated = .true.
         TB_Expwall_name = TB_Expwall(1,1)%Param
      else
         array_is_allocated = .false.
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated#5}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall_name}', TB_Expwall_name) ! below

   if (array_is_allocated) then
      select case (trim(adjustl(TB_Expwall_name)))
      case ('Simple_wall', 'SIMPLE_WALL', 'simple_wall')
         if (.not.allocated(TB_Expwall)) then
            allocate(TB_Exp_wall_simple::TB_Expwall(matter%N_KAO,matter%N_KAO)) ! make it for exponential wall  parametrization
            TB_Expwall%Param = ''
         endif
         TB_Expwall(:,:)%Param = trim(adjustl(TB_Expwall_name))

      case ('General', 'general', 'GENERAL')
         if (.not.allocated(TB_Expwall)) then
            allocate(TB_Short_Rep::TB_Expwall(matter%N_KAO,matter%N_KAO)) ! make it for exponential wall  parametrization
            TB_Expwall%Param = ''
         endif
         TB_Expwall(:,:)%Param = trim(adjustl(TB_Expwall_name))

      case default
         ! nothing to do here
      end select
   endif

   ! Prior to use Exponential wall parameters, we now always have to find out which class the belong to:
   select type (TB_Expwall)
   type is (TB_Exp_wall_simple)
      call MPI_share_Exponential_wall_TB(MPI_param, TB_Expwall)   ! below
   type is (TB_Short_Rep)
      call MPI_share_Short_Rep_TB(MPI_param, TB_Expwall)   ! below
   end select

   !print*, '[MPI process #', MPI_param%process_rank, '] test 0:', allocated(TB_Hamil), trim(adjustl(TB_param_name)), trim(adjustl(TB_repulse_name))
   !pause 'MPI_share_TB_parameters'
   nullify(MPI_param)
#endif
end subroutine MPI_share_TB_parameters




subroutine MPI_share_Short_Rep_TB(MPI_param, TB_Expwall)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Short_Rep), dimension(:,:), intent(inout) ::  TB_Expwall ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j, N1, k
   character(100) :: error_part
   logical :: array_is_allocated

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Short_Rep_TB' ! part of the error message

   do i = 1, size(TB_Expwall,1)
      do j = 1, size(TB_Expwall,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_cut%d0}', TB_Expwall(i,j)%f_cut%d0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_cut%dd}', TB_Expwall(i,j)%f_cut%dd) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_cut_inv%use_it}', TB_Expwall(i,j)%f_cut_inv%use_it) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_cut_inv%d0}', TB_Expwall(i,j)%f_cut_inv%d0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_cut_inv%dd}', TB_Expwall(i,j)%f_cut_inv%dd) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_inv_exp%use_it}', TB_Expwall(i,j)%f_inv_exp%use_it) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_inv_exp%C}', TB_Expwall(i,j)%f_inv_exp%C) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_inv_exp%r0}', TB_Expwall(i,j)%f_inv_exp%r0) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_exp%use_it}', TB_Expwall(i,j)%f_exp%use_it) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_exp%Phi}', TB_Expwall(i,j)%f_exp%Phi) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_exp%r0}', TB_Expwall(i,j)%f_exp%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_exp%a}', TB_Expwall(i,j)%f_exp%a) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_ZBL%use_it}', TB_Expwall(i,j)%f_ZBL%use_it) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_NLH%use_it}', TB_Expwall(i,j)%f_NLH%use_it) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_NLH%a}', TB_Expwall(i,j)%f_NLH%a) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_NLH%b}', TB_Expwall(i,j)%f_NLH%b) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%use_it}', TB_Expwall(i,j)%f_tab%use_it) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%use_spline}', TB_Expwall(i,j)%f_tab%use_spline) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%R}', TB_Expwall(i,j)%f_tab%R) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%E}', TB_Expwall(i,j)%f_tab%E) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%E}', TB_Expwall(i,j)%f_tab%E) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%a}', TB_Expwall(i,j)%f_tab%a) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%b}', TB_Expwall(i,j)%f_tab%b) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%c}', TB_Expwall(i,j)%f_tab%c) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_tab%d}', TB_Expwall(i,j)%f_tab%d) ! below

         ! Check if there is incomming pulse:
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            if (allocated(TB_Expwall(i,j)%f_pow)) then
               array_is_allocated = .true.
               N1 = size(TB_Expwall(i,j)%f_pow)
            else
               array_is_allocated = .false.
               N1 = 0
            endif
         endif
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated}', array_is_allocated) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {N1}', N1) ! below
         if (array_is_allocated) then
            if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
               if (.not.allocated(TB_Expwall(i,j)%f_pow)) allocate(TB_Expwall(i,j)%f_pow(N1))
            endif
            do k = 1, N1
               call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_pow(k)%use_it}', TB_Expwall(i,j)%f_pow(k)%use_it) ! below
               call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_pow(k)%Phi}', TB_Expwall(i,j)%f_pow(k)%Phi) ! below
               call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_pow(k)%r0}', TB_Expwall(i,j)%f_pow(k)%r0) ! below
               call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%f_pow(k)%m}', TB_Expwall(i,j)%f_pow(k)%m) ! below
            enddo
         endif
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Short_Rep_TB



subroutine MPI_share_Exponential_wall_TB(MPI_param, TB_Expwall)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Exp_wall_simple), dimension(:,:), intent(inout) ::  TB_Expwall ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Exponential_wall_TB' ! part of the error message

   do i = 1, size(TB_Expwall,1)
      do j = 1, size(TB_Expwall,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%C}', TB_Expwall(i,j)%C) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%r0}', TB_Expwall(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%d0}', TB_Expwall(i,j)%d0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Expwall(i,j)%dd}', TB_Expwall(i,j)%dd) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Exponential_wall_TB



subroutine MPI_share_Coulomb_cut_TB(MPI_param, TB_Coul)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Coulomb_cut), dimension(:,:), intent(inout) ::  TB_Coul ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Coulomb_cut_TB' ! part of the error message

   do i = 1, size(TB_Coul,1)
      do j = 1, size(TB_Coul,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Coul(i,j)%ind}', TB_Coul(i,j)%ind) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Coul(i,j)%k}', TB_Coul(i,j)%k) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Coul(i,j)%dm}', TB_Coul(i,j)%dm) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Coul(i,j)%dd}', TB_Coul(i,j)%dd) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Coul(i,j)%alpha}', TB_Coul(i,j)%alpha) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Coulomb_cut_TB


subroutine MPI_share_vdW_Dumitrica_TB(MPI_param, TB_Waals)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_vdW_Dumitrica), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_vdW_Dumitrica_TB' ! part of the error message

   do i = 1, size(TB_Waals,1)
      do j = 1, size(TB_Waals,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%C6}', TB_Waals(i,j)%C6) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%alpha}', TB_Waals(i,j)%alpha) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_vdW_Dumitrica_TB



subroutine MPI_share_vdW_ILJ_TB(MPI_param, TB_Waals)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_vdW_ILJ_cut), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_vdW_ILJ_TB' ! part of the error message

   do i = 1, size(TB_Waals,1)
      do j = 1, size(TB_Waals,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%eps}', TB_Waals(i,j)%eps) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%r0}', TB_Waals(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%n}', TB_Waals(i,j)%n) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%m}', TB_Waals(i,j)%m) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%d0_short}', TB_Waals(i,j)%d0_short) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%dd_short}', TB_Waals(i,j)%dd_short) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_vdW_ILJ_TB



subroutine MPI_share_vdW_LJ_TB(MPI_param, TB_Waals)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_vdW_LJ_cut), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_vdW_LJ_TB' ! part of the error message

   do i = 1, size(TB_Waals,1)
      do j = 1, size(TB_Waals,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%eps}', TB_Waals(i,j)%eps) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%r0}', TB_Waals(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%n}', TB_Waals(i,j)%n) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%d0_short}', TB_Waals(i,j)%d0_short) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%dd_short}', TB_Waals(i,j)%dd_short) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_vdW_LJ_TB



subroutine MPI_share_vdW_Girifalco_TB(MPI_param, TB_Waals)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_vdW_Girifalco), dimension(:,:), intent(inout) ::  TB_Waals ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_vdW_Girifalco_TB' ! part of the error message

   do i = 1, size(TB_Waals,1)
      do j = 1, size(TB_Waals,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%C12}', TB_Waals(i,j)%C12) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%C6}', TB_Waals(i,j)%C6) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%r_L}', TB_Waals(i,j)%r_L) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%d_L}', TB_Waals(i,j)%d_L) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%r_S}', TB_Waals(i,j)%r_S) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%d_S}', TB_Waals(i,j)%d_S) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%r_LJ}', TB_Waals(i,j)%r_LJ) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%dm}', TB_Waals(i,j)%dm) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%d_cut}', TB_Waals(i,j)%d_cut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%a}', TB_Waals(i,j)%a) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%b}', TB_Waals(i,j)%b) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%c}', TB_Waals(i,j)%c) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%d}', TB_Waals(i,j)%d) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%dsm}', TB_Waals(i,j)%dsm) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%ds_cut}', TB_Waals(i,j)%ds_cut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%as}', TB_Waals(i,j)%as) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%bs}', TB_Waals(i,j)%bs) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%cs}', TB_Waals(i,j)%cs) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%ds}', TB_Waals(i,j)%ds) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%es}', TB_Waals(i,j)%es) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Waals(i,j)%fs}', TB_Waals(i,j)%fs) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_vdW_Girifalco_TB




subroutine MPI_share_DFTB_TB_repulsive(MPI_param, TB_Repuls)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Rep_DFTB), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_DFTB_TB_repulsive' ! part of the error message

   do i = 1, size(TB_Repuls,1)
      do j = 1, size(TB_Repuls,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%param_name}', TB_Repuls(i,j)%param_name) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%ToP}', TB_Repuls(i,j)%ToP) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%rcut}', TB_Repuls(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%rcut_spline}', TB_Repuls(i,j)%rcut_spline) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%c}', TB_Repuls(i,j)%c) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%a}', TB_Repuls(i,j)%a) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%R}', TB_Repuls(i,j)%R) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%V_rep}', TB_Repuls(i,j)%V_rep) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_DFTB_TB_repulsive


subroutine read_Fu_TB_repulsive(MPI_param, TB_Repuls)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Rep_Fu), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in read_Fu_TB_repulsive' ! part of the error message

   do i = 1, size(TB_Repuls,1)
      do j = 1, size(TB_Repuls,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%E0_TB}', TB_Repuls(i,j)%E0_TB) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%phi0}', TB_Repuls(i,j)%phi0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%d0}', TB_Repuls(i,j)%d0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%m}', TB_Repuls(i,j)%m) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%mc}', TB_Repuls(i,j)%mc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%dc}', TB_Repuls(i,j)%dc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%d1}', TB_Repuls(i,j)%d1) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%dm}', TB_Repuls(i,j)%dm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%a0}', TB_Repuls(i,j)%a0) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%c0}', TB_Repuls(i,j)%c0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%C_a}', TB_Repuls(i,j)%C_a) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine read_Fu_TB_repulsive




subroutine MPI_share_Molteni_TB_repulsive(MPI_param, TB_Repuls)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Rep_Molteni), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Molteni_TB_repulsive' ! part of the error message

   do i = 1, size(TB_Repuls,1)
      do j = 1, size(TB_Repuls,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%NP}', TB_Repuls(i,j)%NP) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%m}', TB_Repuls(i,j)%m) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%phi1}', TB_Repuls(i,j)%phi1) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%phi2}', TB_Repuls(i,j)%phi2) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%r0}', TB_Repuls(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%alpha}', TB_Repuls(i,j)%alpha) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%rcut}', TB_Repuls(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%d}', TB_Repuls(i,j)%d) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%a}', TB_Repuls(i,j)%a) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%b}', TB_Repuls(i,j)%b) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%c}', TB_Repuls(i,j)%c) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Molteni_TB_repulsive




subroutine MPI_share_Pettifor_TB_repulsive(MPI_param, TB_Repuls)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_Rep_Pettifor), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Pettifor_TB_repulsive' ! part of the error message

   do i = 1, size(TB_Repuls,1)
      do j = 1, size(TB_Repuls,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%E0_TB}', TB_Repuls(i,j)%E0_TB) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%phi0}', TB_Repuls(i,j)%phi0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%d0}', TB_Repuls(i,j)%d0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%m}', TB_Repuls(i,j)%m) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%mc}', TB_Repuls(i,j)%mc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%dc}', TB_Repuls(i,j)%dc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%d1}', TB_Repuls(i,j)%d1) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%dm}', TB_Repuls(i,j)%dm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%a0}', TB_Repuls(i,j)%a0) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%c0}', TB_Repuls(i,j)%c0) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Pettifor_TB_repulsive



subroutine MPI_share_xTB_Params(MPI_param, TB_Hamil, TB_Repuls)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_xTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_xTB), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_xTB_Params' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%param_name}', TB_Hamil(i,j)%param_name) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%AO_names}', TB_Hamil(i,j)%AO_names) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rcut}', TB_Hamil(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%d}', TB_Hamil(i,j)%d) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Nprim}', TB_Hamil(i,j)%Nprim) ! below
         ! This part is not ready yet:
         !type(Basis_set_STO), dimension(:), allocatable :: STO ! STO parameteres (via GTO)
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_xTB_Params



subroutine MPI_share_BOP_TB_Params(MPI_param, TB_Hamil, TB_Repuls)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_BOP), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_BOP), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_BOP_TB_Params' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%dcut}', TB_Hamil(i,j)%dcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rcut}', TB_Hamil(i,j)%rcut) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%H_ci}', TB_Hamil(i,j)%H_ci) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%H_li}', TB_Hamil(i,j)%H_li) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%H_ni}', TB_Hamil(i,j)%H_ni) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%E_ci}', TB_Hamil(i,j)%E_ci) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%E_li}', TB_Hamil(i,j)%E_li) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%E_ni}', TB_Hamil(i,j)%E_ni) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%S_ci}', TB_Hamil(i,j)%S_ci) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%S_li}', TB_Hamil(i,j)%S_li) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%S_ni}', TB_Hamil(i,j)%S_ni) ! below

         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%R}', TB_Repuls(i,j)%R) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%V_rep}', TB_Repuls(i,j)%V_rep) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_BOP_TB_Params



subroutine MPI_share_3TB_TB_Params(MPI_param, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_3TB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_3TB_TB_Params' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rcut}', TB_Hamil(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%d}', TB_Hamil(i,j)%d) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rc}', TB_Hamil(i,j)%rc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%include_3body}', TB_Hamil(i,j)%include_3body) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%nullify_diag_cf}', TB_Hamil(i,j)%nullify_diag_cf) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ed}', TB_Hamil(i,j)%Ed) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ep}', TB_Hamil(i,j)%Ep) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Es}', TB_Hamil(i,j)%Es) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Hhavg}', TB_Hamil(i,j)%Hhavg) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Hhcf}', TB_Hamil(i,j)%Hhcf) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Hh3bdy}', TB_Hamil(i,j)%Hh3bdy) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Vrfx}', TB_Hamil(i,j)%Vrfx) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Srfx}', TB_Hamil(i,j)%Srfx) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%V3bdy}', TB_Hamil(i,j)%V3bdy) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_3TB_TB_Params


subroutine MPI_share_DFTB_TB_Params_no_rep(MPI_param, TB_Repuls, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_DFTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_DFTB_no), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_DFTB_TB_Params_no_rep' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%param_name}', TB_Hamil(i,j)%param_name) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rcut}', TB_Hamil(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%d}', TB_Hamil(i,j)%d) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ed}', TB_Hamil(i,j)%Ed) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ep}', TB_Hamil(i,j)%Ep) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Es}', TB_Hamil(i,j)%Es) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ud}', TB_Hamil(i,j)%Ud) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Up}', TB_Hamil(i,j)%Up) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Us}', TB_Hamil(i,j)%Us) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Rr}', TB_Hamil(i,j)%Rr) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Vr}', TB_Hamil(i,j)%Vr) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Sr}', TB_Hamil(i,j)%Sr) ! below
      enddo ! j
   enddo ! i
#endif

end subroutine MPI_share_DFTB_TB_Params_no_rep



subroutine MPI_share_DFTB_TB_Params(MPI_param, TB_Repuls, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_DFTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   type(TB_Rep_DFTB), dimension(:,:), intent(inout) ::  TB_Repuls ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_DFTB_TB_Params' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%param_name}', TB_Hamil(i,j)%param_name) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rcut}', TB_Hamil(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%d}', TB_Hamil(i,j)%d) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ed}', TB_Hamil(i,j)%Ed) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ep}', TB_Hamil(i,j)%Ep) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Es}', TB_Hamil(i,j)%Es) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ud}', TB_Hamil(i,j)%Ud) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Up}', TB_Hamil(i,j)%Up) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Us}', TB_Hamil(i,j)%Us) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Rr}', TB_Hamil(i,j)%Rr) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Vr}', TB_Hamil(i,j)%Vr) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Sr}', TB_Hamil(i,j)%Sr) ! below

         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%param_name}', TB_Repuls(i,j)%param_name) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%ToP}', TB_Repuls(i,j)%ToP) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%rcut}', TB_Repuls(i,j)%rcut) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%c}', TB_Repuls(i,j)%c) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%rcut_spline}', TB_Repuls(i,j)%rcut_spline) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%a}', TB_Repuls(i,j)%a) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%R}', TB_Repuls(i,j)%R) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {TB_Repuls(i,j)%V_rep}', TB_Repuls(i,j)%V_rep) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_DFTB_TB_Params



subroutine MPI_share_Mehl_TB_Hamiltonian(MPI_param, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_NRL), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Mehl_TB_Hamiltonian' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%ind_split}', TB_Hamil(i,j)%ind_split) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%ind_overlap}', TB_Hamil(i,j)%ind_overlap) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%lambd}', TB_Hamil(i,j)%lambd) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Rc}', TB_Hamil(i,j)%Rc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%lden}', TB_Hamil(i,j)%lden) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%al}', TB_Hamil(i,j)%al) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%bl}', TB_Hamil(i,j)%bl) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%cl}', TB_Hamil(i,j)%cl) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%dl}', TB_Hamil(i,j)%dl) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%ellm}', TB_Hamil(i,j)%ellm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%fllm}', TB_Hamil(i,j)%fllm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%gllm}', TB_Hamil(i,j)%gllm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%hllm}', TB_Hamil(i,j)%hllm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%pllm}', TB_Hamil(i,j)%pllm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%qllm}', TB_Hamil(i,j)%qllm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rllm}', TB_Hamil(i,j)%rllm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%sllm}', TB_Hamil(i,j)%sllm) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Mehl_TB_Hamiltonian



subroutine MPI_share_Molteni_TB_Hamiltonian(MPI_param, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_Molteni), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Molteni_TB_Hamiltonian' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Es}', TB_Hamil(i,j)%Es) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ep}', TB_Hamil(i,j)%Ep) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Esa}', TB_Hamil(i,j)%Esa) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%nc}', TB_Hamil(i,j)%nc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rc}', TB_Hamil(i,j)%rc) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%r0}', TB_Hamil(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%n}', TB_Hamil(i,j)%n) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rcut}', TB_Hamil(i,j)%rcut) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%d}', TB_Hamil(i,j)%d) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%V0}', TB_Hamil(i,j)%V0) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Molteni_TB_Hamiltonian


subroutine MPI_share_Fu_TB_Hamiltonian(MPI_param, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_Fu), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Pettifor_TB_Hamiltonian' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Es}', TB_Hamil(i,j)%Es) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ep}', TB_Hamil(i,j)%Ep) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%r0}', TB_Hamil(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%n}', TB_Hamil(i,j)%n)   ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%r1}', TB_Hamil(i,j)%r1) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rm}', TB_Hamil(i,j)%rm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%nc}', TB_Hamil(i,j)%nc) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rc}', TB_Hamil(i,j)%rc) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c0}', TB_Hamil(i,j)%c0) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c1}', TB_Hamil(i,j)%c1) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c2}', TB_Hamil(i,j)%c2) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c3}', TB_Hamil(i,j)%c3) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%V0}', TB_Hamil(i,j)%V0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%C_a}', TB_Hamil(i,j)%C_a) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Fu_TB_Hamiltonian



subroutine MPI_share_Pettifor_TB_Hamiltonian(MPI_param, TB_Hamil)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(TB_H_Pettifor), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !--------------------
   integer :: i, j
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Pettifor_TB_Hamiltonian' ! part of the error message

   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,1)
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Es}', TB_Hamil(i,j)%Es) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%Ep}', TB_Hamil(i,j)%Ep) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%r0}', TB_Hamil(i,j)%r0) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%n}', TB_Hamil(i,j)%n)   ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%r1}', TB_Hamil(i,j)%r1) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rm}', TB_Hamil(i,j)%rm) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%nc}', TB_Hamil(i,j)%nc) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%rc}', TB_Hamil(i,j)%rc) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c0}', TB_Hamil(i,j)%c0) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c1}', TB_Hamil(i,j)%c1) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c2}', TB_Hamil(i,j)%c2) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%c3}', TB_Hamil(i,j)%c3) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {TB_Hamil(i,j)%V0}', TB_Hamil(i,j)%V0) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Pettifor_TB_Hamiltonian



subroutine MPI_share_laser(numpar, laser)
   type(Numerics_param), intent(inout), target :: numpar ! all numerical parameters
   type(Pulse), dimension(:), allocatable, intent(inout) :: laser	! Laser pulse parameters
   !--------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   integer :: Nsiz, N1, N2, N3, i
   character(100) :: error_part
   logical :: array_is_allocated, array_is_allocated2

#ifdef MPI_USED
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_laser' ! part of the error message

   ! Check if there is incomming pulse:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(laser)) then
         array_is_allocated = .true.
         N1 = size(laser)
      else
         array_is_allocated = .false.
         N1 = 0
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {N1#1}', N1) ! below

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(laser)) allocate(laser(N1))
      endif

      do i = 1, N1   ! for all pulses
         ! kind of pulse: 0 = flat-top, 1 = Gaussian, 2 = SASE:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%KOP}', laser(i)%KOP) ! below

         ! [eV] photon energy:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%hw}', laser(i)%hw) ! below

         ! [eV] distribution of photon energy spectrum (assumed gaussian):
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%FWHM_hw}', laser(i)%FWHM_hw) ! below

         ! [fs] pulse duration:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%t}', laser(i)%t) ! below

         ! [fs] pulse center position:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%t0}', laser(i)%t0) ! below

         ! [J/cm^2]  incoming fluence:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%F_in}', laser(i)%F_in) ! below

         ! [eV/atom] absorbed dose:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%F}', laser(i)%F) ! below

         ! [eV] total absorbed energy per simulation box:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%Fabs}', laser(i)%Fabs) ! below

         ! number of absorbed photons
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {laser(i)%Nph}', laser(i)%Nph) ! below

         !print*, '[MPI process #', MPI_param%process_rank, '] test:', i, laser(i)%KOP, laser(i)%hw, laser(i)%FWHM_hw, laser(i)%t, laser(i)%t0, laser(i)%F_in, laser(i)%F, laser(i)%Fabs, laser(i)%Nph
      enddo
   endif ! array_is_allocated



! type Pulse
!    integer :: KOP    ! kind of pulse: 0 = flat-top, 1 = Gaussian, 2 = SASE
!    real(8) :: hw     ! [eV] photon energy
!    real(8) :: FWHM_hw   ! [eV] distribution of photon energy spectrum (assumed gaussian)
!    real(8) :: t      ! [fs] pulse duration
!    real(8) :: t0     ! [fs] pulse maximum position
!    real(8) :: F_in   ! [J/cm^2]  incoming fluence
!    real(8) :: F      ! [eV/atom] absorbed dose
!    real(8) :: Fabs   ! [eV] total absorbed energy per simulation box
!    real(8) :: Nph    ! number of absorbed photons
! end type Pulse

   nullify(MPI_param)
#endif
end subroutine MPI_share_laser




subroutine MPI_share_Scell(numpar, Scell)
   type(Numerics_param), intent(inout), target :: numpar ! all numerical parameters
   type(Super_cell), dimension(:), allocatable, intent(inout) :: Scell ! suoer-cell with all the atoms inside
   !--------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   integer :: Nsiz, N1, N2, N3, i
   character(100) :: error_part
   logical :: array_is_allocated, array_is_allocated2

#ifdef MPI_USED
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_Scell' ! part of the error message

   !print*, '[MPI process #', MPI_param%process_rank, '] test 0:', allocated(Scell)

   ! Only share the data the master process got from the input files (and correspondingly defined):
   if (MPI_param%process_rank /= 0) then   ! only MPI master process does it
      if (.not.allocated(Scell)) then  ! initialize it
         allocate(Scell(1)) ! for the moment, only one super-cell
         ! Undefined yet variables:
         ! Scell(1)%Subcell, Scell(1)%MCholeScell(1)%Ta_sub, Scell(1)%MSDP
         !
         ! Unused yet variables, set default:
         Scell(1)%Na = 0
         Scell(1)%Ne = 0
         Scell(1)%Ne_low = 0.0d0
         Scell(1)%Ne_high = 0.0d0
         Scell(1)%Ne_emit = 0.0d0
         Scell(1)%Nh = 0.0d0
         Scell(1)%Nph = 0.0d0
         Scell(1)%Q = 0.0d0
         Scell(1)%Ne_CB = 0.0d0
         Scell(1)%Ta_var(:) = 0.0d0
         Scell(1)%Ta_r_var(:) = 0.0d0
         Scell(1)%Ta_conf_run_average(:) = 0.0d0
         Scell(1)%Fv = 0.0d0
         Scell(1)%Tconf = 0.0d0
         Scell(1)%Tconf2 = 0.0d0
         Scell(1)%Pressure = 0.0d0
         Scell(1)%Stress = 0.0d0
         Scell(1)%Pot_Pressure = 0.0d0
         Scell(1)%Pot_Stress = 0.0d0
         Scell(1)%MSD = 0.0d0
      endif
   endif

   ! Variables defined, to share:
   ! Electronic temperature:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Te}', Scell(1)%Te) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%TeeV}', Scell(1)%TeeV) ! below

   ! Atomic temperature:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ta}', Scell(1)%Ta) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%TaeV}', Scell(1)%TaeV) ! below


   ! [A] mean displacements for atoms masked:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(Scell(1)%Displ)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(Scell(1)%Displ)
      else
         N1 = 0
      endif
   endif
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {array_is_allocated}', array_is_allocated) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {N1#1}', N1) ! below

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(Scell(1)%Displ)) allocate(Scell(1)%Displ(N1))
      endif

      do i = 1, N1   ! for all allocated bits
         ! name of the mask:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%mask_name}', Scell(1)%Displ(i)%mask_name) ! below

         ! power of the mean displacement for this particular analysis:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%MSD_power}', Scell(1)%Displ(i)%MSD_power) ! below


         ! [A^MSD_power] mean displacement of all atoms:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%mean_disp}', Scell(1)%Displ(i)%mean_disp) ! below

         ! [A^MSD_power] mean displacement of atoms by sort:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%mean_disp_sort}', Scell(1)%Displ(i)%mean_disp_sort) ! below

         ! along which exis the user requirested the analysis:
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%print_r}', Scell(1)%Displ(i)%print_r) ! below

         ! [A^MSD_power] mean displacements along X, Y, Z axes:
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%mean_disp_r}', Scell(1)%Displ(i)%mean_disp_r) ! below

         ! [A^MSD_power] mean displacements along X, Y, Z axes by sort:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%mean_disp_r_sort}', Scell(1)%Displ(i)%mean_disp_r_sort) ! below

         ! atomic mask to be used:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%Atomic_mask}', Scell(1)%Displ(i)%Atomic_mask) ! below

         ! which one is used:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%logical_and}', Scell(1)%Displ(i)%logical_and) ! below
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%logical_or}', Scell(1)%Displ(i)%logical_or) ! below

         ! index of the axis for the section:
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%axis_ind}', Scell(1)%Displ(i)%axis_ind) ! below

         ! section starting and ending points:
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%r_start}', Scell(1)%Displ(i)%r_start) ! below
         call broadcast_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Displ(i)%r_end}', Scell(1)%Displ(i)%r_end) ! below
      enddo ! i
   endif ! array_is_allocated

   ! Chemical potential of electrons:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%mu}', Scell(1)%mu) ! below

   ! Electronic entropy:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Se}', Scell(1)%Se) ! below

   ! Equivalent electronic entropy:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Se_eq}', Scell(1)%Se_eq) ! below

   ! electron distribution:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fe}', Scell(1)%fe) ! below

   ! equivalent Fermi electron distribution:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fe_eq}', Scell(1)%fe_eq) ! below

   ! energy grid for atomic distribution:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ea_grid}', Scell(1)%Ea_grid) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ea_grid_out}', Scell(1)%Ea_grid_out) ! below

   ! atomic distribution and equivalent Maxwell distribution:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa}', Scell(1)%fa) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_eq}', Scell(1)%fa_eq) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_out}', Scell(1)%fa_out) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_eq_out}', Scell(1)%fa_eq_out) ! below

   ! energy grid for atomic distribution
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ea_pot_grid}', Scell(1)%Ea_pot_grid) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ea_pot_grid_out}', Scell(1)%Ea_pot_grid_out) ! below

   ! atomic distribution of potential energies
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_pot}', Scell(1)%fa_pot) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_eq_pot}', Scell(1)%fa_eq_pot) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_pot_out}', Scell(1)%fa_pot_out) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_eq_pot_out}', Scell(1)%fa_eq_pot_out) ! below

   ! total-energy grid for atomic distribution:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ea_tot_grid}', Scell(1)%Ea_tot_grid) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%Ea_tot_grid_out}', Scell(1)%Ea_tot_grid_out) ! below

   ! atomic distribution of total energies
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_tot}', Scell(1)%fa_tot) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {Scell(1)%fa_tot_out}', Scell(1)%fa_tot_out) ! below

   ! Share optical parameters:
   call MPI_share_eps(numpar, Scell(1)%eps)  ! below

   nullify(MPI_param)
#endif
end subroutine MPI_share_Scell



subroutine MPI_share_eps(numpar, eps)
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   type(Drude), intent(inout) :: eps	! epsylon, dielectric function and its parameters
   !-------------------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   character(100) :: error_part

#ifdef MPI_USED
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_eps' ! part of the error message

   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%ReEps}', eps%ReEps) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%ReEps0}', eps%ReEps0) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%ImEps}', eps%ImEps) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%ImEps0}', eps%ImEps0) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_xx}', eps%Eps_xx) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_yy}', eps%Eps_yy) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_zz}', eps%Eps_zz) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_xy}', eps%Eps_xy) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_xz}', eps%Eps_xz) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_yx}', eps%Eps_yx) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_yz}', eps%Eps_yz) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_zx}', eps%Eps_zx) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%Eps_zy}', eps%Eps_zy) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%n}', eps%n) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%k}', eps%k) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%R}', eps%R) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%T}', eps%T) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%A}', eps%A) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%dc_cond}', eps%dc_cond) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%w}', eps%w) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%l}', eps%l) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%tau}', eps%tau) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%me_eff}', eps%me_eff) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%mh_eff}', eps%mh_eff) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%tau_e}', eps%tau_e) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%tau_h}', eps%tau_h) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%teta}', eps%teta) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%dd}', eps%dd) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%all_w}', eps%all_w) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%KK}', eps%KK) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%E_min}', eps%E_min) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%E_max}', eps%E_max) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {eps%dE}', eps%dE) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {eps%Eps_hw}', eps%Eps_hw) ! below

   nullify(MPI_param)
#endif
end subroutine MPI_share_eps



subroutine MPI_share_numpar(numpar)
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   !-------------------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   integer :: Nsiz, N1, N2, N3, i
   character(100) :: error_part
   logical :: array_is_allocated

#ifdef MPI_USED
   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_numpar' ! part of the error message


   ! use linear scaling TB (1), or not (0)
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%lin_scal}', numpar%lin_scal) ! below

   ! number of subcells along each axis: X, Y, Z
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {numpar%N_subcels}', numpar%N_subcels) ! below

   ! pair correlation function (if required):
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%Subcell_coord_sx}', numpar%Subcell_coord_sx) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%Subcell_coord_sy}', numpar%Subcell_coord_sy) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%Subcell_coord_sz}', numpar%Subcell_coord_sz) ! below

   ! number of input file used (for using more then one sequentially):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%which_input}', numpar%which_input) ! below

   ! verbose options:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%verbose}', numpar%verbose) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%nonverbose}', numpar%nonverbose) ! below

   ! flag to recalculate mean free paths, if needed:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%redo_MFP}', numpar%redo_MFP) ! below

   ! flag to printout mean free paths:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%print_MFP}', numpar%print_MFP) ! below

   ! flag to use the distribution from a file:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%fe_input_exists}', numpar%fe_input_exists) ! below

   ! file name with user-provided initial electronic distribution:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%fe_filename}', numpar%fe_filename) ! below

   ! initial distribution function:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%fe_input}', numpar%fe_input) ! below

   ! DOS for high-energy electron distribution, if required
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%high_DOS}', numpar%high_DOS) ! below

   ! number of time-steps over which to average the distribution on the grid:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%fe_aver_num}', numpar%fe_aver_num) ! below

   ! [fs] characteristic relaxation time of ALL electrons (used for the relaxation time approximation):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%tau_fe}', numpar%tau_fe) ! below

   ! [fs] characteristic relaxation time of CB electrons (used for the relaxation time approximation):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%tau_fe_CB}', numpar%tau_fe_CB) ! below

   ! [fs] characteristic relaxation time of VB electrons (used for the relaxation time approximation):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%tau_fe_VB}', numpar%tau_fe_VB) ! below

   ! flag for band-resolved thermalization:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_partial_thermal}', numpar%do_partial_thermal) ! below

   ! [fs] time-step for MD:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%dt}', numpar%dt) ! below

   ! dt/2, often used:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%halfdt}', numpar%halfdt) ! below

   ! dt*dt/2, often used:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%dtsqare}', numpar%dtsqare) ! below

   ! dt^3/6, often used:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%dt3}', numpar%dt3) ! below

   ! dt^4/48, often used:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%dt4}', numpar%dt4) ! below

   ! 0=Verlet (2d order); 1=Yoshida (4th order); 2=Martyna(4th order):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%MD_algo}', numpar%MD_algo) ! below

   ! grid, when to change the MD timestep:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%dt_MD_reset_grid}', numpar%dt_MD_reset_grid) ! below

   ! grid, when to change the MD timestep:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%dt_MD_grid}', numpar%dt_MD_grid) ! below

   ! which timestep from the array "dt_MD_grid" to use now:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%i_dt}', numpar%i_dt) ! below

   ! filename with MD time step grid
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%MD_step_grid_file}', numpar%MD_step_grid_file) ! below

   ! flag for various atomic temperature definitions:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%print_Ta}', numpar%print_Ta) ! below

   ! index to set starting velocity distribution: 1=linear; 2=Maxwellian
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%ind_starting_V}', numpar%ind_starting_V) ! below

   ! index to mark whether velocities were read from a file or set:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%vel_from_file}', numpar%vel_from_file) ! below

   ! EPICS data: EADL, EPDL, EEDL:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%EADL_file}', numpar%EADL_file) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%EPDL_file}', numpar%EPDL_file) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%EEDL_file}', numpar%EEDL_file) ! below

   ! grid, when to change the Atomic bath parameters:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%At_bath_reset_grid}', numpar%At_bath_reset_grid) ! below

   ! Atomic bath temperatures array [K]:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%At_bath_grid_Ta}', numpar%At_bath_grid_Ta) ! below

   ! Atomic bath characteristic times array [fs]:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%At_bath_grid_tau}', numpar%At_bath_grid_tau) ! below

   ! which timestep from the array to use now:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%i_At_bath_dt}', numpar%i_At_bath_dt) ! below

   ! filename with Atomic bath parameters:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%At_bath_step_grid_file}', numpar%At_bath_step_grid_file) ! below

   ! grid, when to change the Electronic bath parameters:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%El_bath_reset_grid}', numpar%El_bath_reset_grid) ! below

   ! Electronic bath temperatures array [K]:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%El_bath_grid_Ta}', numpar%El_bath_grid_Ta) ! below

   ! Electronic bath characteristic times array [fs]:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%El_bath_grid_tau}', numpar%El_bath_grid_tau) ! below

   ! which timestep from the array to use now:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%i_El_bath_dt}', numpar%i_El_bath_dt) ! below

   ! filename with Electronic bath parameters:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%El_bath_step_grid_file}', numpar%El_bath_step_grid_file) ! below

   !-----------------
   ! [fs] starting time of simulation:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%t_start}', numpar%t_start) ! below

   ! [fs] total time of simulation:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%t_total}', numpar%t_total) ! below

   ! time when we switch from Te=const, to Ee=const [fs] / negative value, when not using it:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%t_Te_Ee}', numpar%t_Te_Ee) ! below

   ! time when we switch on nonadiabatic terms [fs]:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%t_NA}', numpar%t_NA) ! below

   ! time-step how often to cool down the atoms [fs]:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%dt_cooling}', numpar%dt_cooling) ! below

   ! P=const, otherwise V=const for Parinello-Rahman MD simulations:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%p_const}', numpar%p_const) ! below

   ! true=included / false=excluded nonadiabatic coupling:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Nonadiabat}', numpar%Nonadiabat) ! below

   ! non-adiabatic model kind:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%NA_kind}', numpar%NA_kind) ! below

   ! model for atomic distribution in the nonadiabatic coupling:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%ind_at_distr}', numpar%ind_at_distr) ! below

   ! include self-consistent charge corrections in TB?:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%scc}', numpar%scc) ! below

   ! index for the model for gamma in scc term:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%scc_gam_ind}', numpar%scc_gam_ind) ! below

   ! mixing factor for scc calculations:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%scc_mix}', numpar%scc_mix) ! below

   ! index for the size of the basis set used: s=0; sp3=1; sp3d5=2; sp3s*=3; sp3d5s*=4;:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%basis_size_ind}', numpar%basis_size_ind) ! below

   ! basis set size (orbitals per atom):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%N_basis_size}', numpar%N_basis_size) ! below

   ! true=included / false=excluded ; for atoms:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Transport}', numpar%Transport) ! below

   ! true=included / false=excluded ; for electrons:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Transport_e}', numpar%Transport_e) ! below

   ! [eV] cut-off energy, separating low-energy-electrons from high-energy-electrons:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%E_cut}', numpar%E_cut) ! below

   ! include evolution of E_cut due to changes of the band struture:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%E_cut_dynamic}', numpar%E_cut_dynamic) ! below

   ! [eV] cut-off energy (work function) above which electrons are emitted from the material:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%E_work}', numpar%E_work) ! below

   ! [eV] Coulomb energy attracting electrons back to the material:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%E_Coulomb}', numpar%E_Coulomb) ! below

   ! number of MC iterations:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%NMC}', numpar%NMC) ! below

   ! Number of threads for openmp:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%NOMP}', numpar%NOMP) ! below

   ! save data into files every 'dt_save' [fs]:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%dt_save}', numpar%dt_save) ! below

   ! multiplication factor often used in MD:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%fctr}', numpar%fctr) ! below

   ! [eV] acceptance window for nonadiabatic coupling:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%acc_window}', numpar%acc_window) ! below

   ! [eV] window to exclude quasidegenerate levels in nonadiabatic coupling:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%degeneracy_eV}', numpar%degeneracy_eV) ! below

   ! scaling factor for electron-ion coupling matrix element:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%M2_scaling}', numpar%M2_scaling) ! below

   ! [fs], cooling of atoms: when to start, and how often to do:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%at_cool_dt}', numpar%at_cool_dt) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%at_cool_start}', numpar%at_cool_start) ! below

   ! smearing to use for DOS calculations:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Smear_DOS}', numpar%Smear_DOS) ! below

   ! user-defined radius to count nearest neighbors (for printout only, not for MD calculations!):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%NN_radius}', numpar%NN_radius) ! below

   ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%MSD_power}', numpar%MSD_power) ! below

   ! model to split DOS:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%DOS_splitting}', numpar%DOS_splitting) ! below

   ! to identify and separate different bands:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%mask_DOS}', numpar%mask_DOS) ! below

   ! use weights on energy levels accourding to DOS_weights:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%weigthed_DOS}', numpar%weigthed_DOS) ! below

   ! to identify and separate different bands:
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%DOS_weights}', numpar%DOS_weights) ! below

   ! input folder address:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%input_path}', numpar%input_path) ! below

   ! output folder address:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%output_path}', numpar%output_path) ! below

   ! This variable was already defined in all processes:
!    character(1) :: path_sep   ! path separator

   ! where to take atomic data from (EADL, CDF, XATOM...):
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%At_base}', numpar%At_base) ! below

   ! [eV] to start with, user-defined gap:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%user_defined_E_gap}', numpar%user_defined_E_gap) ! below

   ! flag for embedding in water:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%embed_water}', numpar%embed_water) ! below

   ! how many water molecules to use:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%N_water_mol}', numpar%N_water_mol) ! below

   ! File numbers:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_temperatures}', numpar%FN_temperatures) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_energies}', numpar%FN_energies) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_atoms_R}', numpar%FN_atoms_R) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_atoms_S}', numpar%FN_atoms_S) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_supercell}', numpar%FN_supercell) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_electron_properties}', numpar%FN_electron_properties) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_numbers}', numpar%FN_numbers) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_all_w}', numpar%FN_all_w) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_deep_holes}', numpar%FN_deep_holes) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Ei}', numpar%FN_Ei) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_fe}', numpar%FN_fe) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_PCF}', numpar%FN_PCF) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_optics}', numpar%FN_optics) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_parameters}', numpar%FN_parameters) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_communication}', numpar%FN_communication) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_cif}', numpar%FN_cif) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_pressure}', numpar%FN_pressure) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_DOS}', numpar%FN_DOS) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_coupling}', numpar%FN_coupling) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_neighbors}', numpar%FN_neighbors) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Ce}', numpar%FN_Ce) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_kappa}', numpar%FN_kappa) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_kappa_dyn}', numpar%FN_kappa_dyn) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Se}', numpar%FN_Se) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_fe_on_grid}', numpar%FN_fe_on_grid) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Te}', numpar%FN_Te) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_mu}', numpar%FN_mu) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_orb_resolved}', numpar%FN_orb_resolved) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_fa}', numpar%FN_fa) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Sa}', numpar%FN_Sa) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Ta}', numpar%FN_Ta) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_fa_pot}', numpar%FN_fa_pot) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_Ta_part}', numpar%FN_Ta_part) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_fa_tot}', numpar%FN_fa_tot) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%FN_testmode}', numpar%FN_testmode) ! below
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%FN_displacements}', numpar%FN_displacements) ! below

   ! time when the communication.txt file was last modified:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%MOD_TIME}', numpar%MOD_TIME) ! below

   ! Optical model:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%drude_ray}', numpar%drude_ray) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%optic_model}', numpar%optic_model) ! below

   ! Coupling scheme:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%el_ion_scheme}', numpar%el_ion_scheme) ! below

   ! number of k-points in each direction for eps-calculations:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%ixm}', numpar%ixm) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%iym}', numpar%iym) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%izm}', numpar%izm) ! below

   ! for the case of user-provided grid for k-space (for CDF and DOS calculations):
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {numpar%k_grid}', numpar%k_grid) ! below

   ! periodic boundaries in each of the three spatial dimensions:
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {numpar%r_periodic}', numpar%r_periodic) ! below

   ! Flags to save output:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_Ei}', numpar%save_Ei) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_fe}', numpar%save_fe) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_PCF}', numpar%save_PCF) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_XYZ}', numpar%save_XYZ) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_drude}', numpar%do_drude) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_cool}', numpar%do_cool) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_atoms}', numpar%do_atoms) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%change_size}', numpar%change_size) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%allow_rotate}', numpar%allow_rotate) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_fe_grid}', numpar%save_fe_grid) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_fe_orb}', numpar%save_fe_orb) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_fa}', numpar%save_fa) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_testmode}', numpar%save_testmode) ! below

   ! Flags of save_XYZ_extra indices: (1) atomic mass; (2) atomic charge; (3) kinetic energy
   call broadcast_array(MPI_param, trim(adjustl(error_part))//' {numpar%save_XYZ_extra}', numpar%save_XYZ_extra) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_elastic_MC}', numpar%do_elastic_MC) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_path_coordinate}', numpar%do_path_coordinate) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_kappa}', numpar%do_kappa) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_DOS}', numpar%do_DOS) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%do_kappa_dyn}', numpar%do_kappa_dyn) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_CIF}', numpar%save_CIF) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_pressure}', numpar%save_pressure) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_DOS}', numpar%save_DOS) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_raw}', numpar%save_raw) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_NN}', numpar%save_NN) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_CDF}', numpar%save_CDF) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Mulliken_model}', numpar%Mulliken_model) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%ind_fig_extention}', numpar%ind_fig_extention) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%change_size_step}', numpar%change_size_step) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%change_size_max}', numpar%change_size_max) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%change_size_min}', numpar%change_size_min) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%kappa_Te_min}', numpar%kappa_Te_min) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%kappa_Te_max}', numpar%kappa_Te_max) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%kappa_dTe}', numpar%kappa_dTe) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%kappa_model}', numpar%kappa_model) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%input_CDF_file}', numpar%input_CDF_file) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%fig_extention}', numpar%fig_extention) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%create_BOP_repulse}', numpar%create_BOP_repulse) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%BOP_Folder_name}', numpar%BOP_Folder_name) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Filename_communication}', numpar%Filename_communication) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%BOP_bond_length}', numpar%BOP_bond_length) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%power_b}', numpar%power_b) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%Cell_filename}', numpar%Cell_filename) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%save_files_used}', numpar%save_files_used) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%numpar_in_input}', numpar%numpar_in_input) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%output_extra_name}', numpar%output_extra_name) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%output_name_add}', numpar%output_name_add) ! below

   nullify(MPI_param)
#endif
end subroutine MPI_share_numpar


subroutine MPI_share_matter(numpar, matter)
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   type(Solid), intent(inout) :: matter	! all material parameters
   !---------------------------
   type(Used_MPI_parameters), pointer :: MPI_param
   integer :: Nsiz, N1, N2, i
   character(100) :: error_part
   logical :: array_is_allocated

#ifdef MPI_USED

   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_matter' ! part of the error message

   ! Material name
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Name}', matter%Name) ! below

   ! Chemical formula:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Chem}', matter%Chem) ! below

   ! number of different kinds of atoms in this compound:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%N_KAO}', matter%N_KAO) ! below

   ! pair correlation function (if required):
   call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%PCF}', matter%PCF) ! below

   ! number of unit-cells in x, y and z directions:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%cell_x}', matter%cell_x) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%cell_y}', matter%cell_y) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%cell_z}', matter%cell_z) ! below

   ! [kg] Parinello_Rahman super-cell mass:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%W_PR}', matter%W_PR) ! below

   ! [Pa] external pressure applied:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%p_ext}', matter%p_ext) ! below

   ! [g/cm^3] density of the material:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%dens}', matter%dens) ! below

   ! [1/cm^3] atomic density of the material:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%At_dens}', matter%At_dens) ! below

   ! [eV] plasmon energy:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%hw_plasma}', matter%hw_plasma) ! below

   ! [K] Bath temeprature for Berendsen thermostat for atoms:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%T_bath}', matter%T_bath) ! below

   ! [K] Bath temeprature for Berendsen thermostat for electrons:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%T_bath_e}', matter%T_bath_e) ! below

   ! [fs] time constant of cooling via thermostat for atoms:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%tau_bath}', matter%tau_bath) ! below

   ! [fs] time constant of cooling via thermostat for electrons:
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%tau_bath_e}', matter%tau_bath_e) ! below

   ! all kinds of atoms of the compound:
   if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(matter%Atoms)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(matter%Atoms)
      else
         N1 = 0
      endif
   endif
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {array_is_allocated#2}') ! module "MPI_subroutines"
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N1#2}') ! module "MPI_subroutines"
   if (numpar%MPI_param%process_rank /= 0) then   ! MPI non-master process
      if (array_is_allocated) then ! in the MASTER process it is allocated, so allocate it in non-master processes too
         allocate(matter%Atoms(N1))
      endif
   endif

   if (array_is_allocated) then ! in the MASTER process it is allocated, so allocate it in non-master processes too
      do i = 1, N1   ! for all atoms
         ! Chemical element name:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Name}', matter%Atoms(i)%Name) ! below

         ! [kg] atomic mass:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Ma}', matter%Atoms(i)%Ma) ! below

         ! atomic number:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Z}', matter%Atoms(i)%Z) ! below

         ! number of shells of the element:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%sh}', matter%Atoms(i)%sh) ! below

         ! contribution to the compound:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%percentage}', matter%Atoms(i)%percentage) ! below

         ! number of shells of the element:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%NVB}', matter%Atoms(i)%NVB) ! below

         ! [eV] Hubbard U for SCC:
         call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Hubbard_U}', matter%Atoms(i)%Hubbard_U) ! below

         ! These variables are currently unused, so skip broadcasting:
!        type(Basis_set), dimension(:), allocatable :: Cart_Basis   ! Cartesian GTO basis set functions for this element (if xTB is used)
!        type(Basis_set), dimension(:), allocatable :: Spher_Basis  ! Spherical (pure) GTO basis set functions for this element (if xTB is used)

         ! These parameters are undefined at the beginning, to be caluclated later (skip them):
!        real(8) :: mulliken_Ne   ! electron population according to Mulliken analysis
!        real(8) :: mulliken_q   ! charge according to Mulliken analysis

         ! EADL shell designator:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Shl_dsgnr}', matter%Atoms(i)%Shl_dsgnr) ! below

         ! Shell name:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Shell_name}', matter%Atoms(i)%Shell_name) ! below

         ! [eV] ionization potentials for all shells:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Ip}', matter%Atoms(i)%Ip) ! below

         ! [eV] mean kinetic energy of all shells:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Ek}', matter%Atoms(i)%Ek) ! below

         ! number of electrons in each shell:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Ne_shell}', matter%Atoms(i)%Ne_shell) ! below

         ! [fs] Auger-decay times for all shells:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Auger}', matter%Atoms(i)%Auger) ! below

         ! type of electron scattering cross-section used for each shell:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%TOCS}', matter%Atoms(i)%TOCS) ! below

         ! type of photon scattering cross-section used for each shell:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%TOCSph}', matter%Atoms(i)%TOCSph) ! below

         ! current number of deep-shell holes in each shell:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Nh_shell}', matter%Atoms(i)%Nh_shell) ! below

         ! Number of CDF functions:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%N_CDF}', matter%Atoms(i)%N_CDF) ! below

         ! EADL shell designator for atomic shells without VB:
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {matter%Atoms(i)%Shl_dsgnr_atomic}', matter%Atoms(i)%Shl_dsgnr_atomic) ! below

         ! Derived types have to be processes individually:
         ! CDF:
         if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            if (allocated(matter%Atoms(i)%CDF)) then
               array_is_allocated = .true.
            else
               array_is_allocated = .false.
            endif
            if (array_is_allocated) then
               N1 = size(matter%Atoms(i)%CDF)
            else
               N1 = 0
            endif
         endif
         call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {array_is_allocated#3}') ! module "MPI_subroutines"
         call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N1#3}') ! module "MPI_subroutines"
         if (numpar%MPI_param%process_rank /= 0) then   ! MPI non-master process
            if (array_is_allocated) then ! in the MASTER process it is allocated, so allocate it in non-master processes too
               allocate(matter%Atoms(i)%CDF(N1))
            endif
         endif


         if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            if (allocated(matter%Atoms(i)%El_MFP)) then
               array_is_allocated = .true.
            else
               array_is_allocated = .false.
            endif
            if (array_is_allocated) then
               N1 = size(matter%Atoms(i)%El_MFP)
            else
               N1 = 0
            endif
         endif
         call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {array_is_allocated#4}') ! module "MPI_subroutines"
         call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N1#4}') ! module "MPI_subroutines"
         if (numpar%MPI_param%process_rank /= 0) then   ! MPI non-master process
            if (array_is_allocated) then ! in the MASTER process it is allocated, so allocate it in non-master processes too
               allocate(matter%Atoms(i)%El_MFP(N1))
            endif
         endif


         if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            if (allocated(matter%Atoms(i)%Ph_MFP)) then
               array_is_allocated = .true.
            else
               array_is_allocated = .false.
            endif
            if (array_is_allocated) then
               N1 = size(matter%Atoms(i)%Ph_MFP)
            else
               N1 = 0
            endif
         endif
         call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {array_is_allocated#5}') ! module "MPI_subroutines"
         call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N1#5}') ! module "MPI_subroutines"
         if (numpar%MPI_param%process_rank /= 0) then   ! MPI non-master process
            if (array_is_allocated) then ! in the MASTER process it is allocated, so allocate it in non-master processes too
               allocate(matter%Atoms(i)%Ph_MFP(N1))
            endif
         endif


         if (numpar%MPI_param%process_rank == 0) then   ! only MPI master process does it
            if (allocated(matter%Atoms(i)%El_MFP_vs_T)) then
               array_is_allocated = .true.
            else
               array_is_allocated = .false.
            endif
            if (array_is_allocated) then
               N1 = size(matter%Atoms(i)%El_MFP_vs_T)
            else
               N1 = 0
            endif
         endif
         call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {array_is_allocated#6}') ! module "MPI_subroutines"
         call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N1#6}') ! module "MPI_subroutines"
         if (numpar%MPI_param%process_rank /= 0) then   ! MPI non-master process
            if (array_is_allocated) then ! in the MASTER process it is allocated, so allocate it in non-master processes too
               if (.not.allocated(matter%Atoms(i)%El_MFP_vs_T)) allocate(matter%Atoms(i)%El_MFP_vs_T(N1))
            endif
         endif

         ! These functions are not defined yet, to be read from files later:
!        type(Ritchi), dimension(:), allocatable :: CDF  ! coefficients of CDF
!        type(MFP), dimension(:), allocatable :: El_MFP  ! electron inelastic mean free paths for each shell (inverted [1/A])
!        type(MFP), dimension(:), allocatable :: Ph_MFP  ! photon mean free paths for each shell (inverted [1/A])
!        type(MFP), dimension(:), allocatable :: El_MFP_vs_T ! electron MFP for (inverted [1/A]) for different Te
!        type(MFP) :: El_EMFP ! electron elastic mean free paths (inverted [1/A])

         !print*, '[MPI process #', MPI_param%process_rank, '] :', allocated(matter%Atoms(i)%El_MFP), allocated(matter%Atoms(i)%Ph_MFP), allocated(matter%Atoms(i)%El_MFP_vs_T), allocated(matter%Atoms(i)%El_EMFP%E)

      enddo
   endif

   ! These functions are not defined yet, to be read from files later:
!    type(MFP) :: El_MFP_tot  ! Total electron inelastic mean free paths (inverted [1/A])
!    type(MFP) :: El_EMFP_tot ! Total electron elastic mean free paths (inverted [1/A])
!    type(MFP) :: Ph_MFP_tot  ! Total photon mean free paths for each shell (inverted [1/A])

!print*, '[MPI process #', MPI_param%process_rank, '] :', allocated(matter%El_MFP_tot%E), allocated(matter%El_EMFP_tot%E), allocated(matter%Ph_MFP_tot%E)

   nullify(MPI_param)
#endif
end subroutine MPI_share_matter



subroutine MPI_share_Ritchi_CDF(MPI_param, matter)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(Solid), intent(inout) :: matter                  ! all material parameters
   !-----------------------
   character(100) :: error_part, error_report
   logical :: array_is_allocated
   integer :: N1, i, j, Nat, Nsh

#ifdef MPI_USED
   error_part = 'ERROR in MPI_share_Ritchi_CDF' ! part of the error message

   Nat = size(matter%Atoms)

   ! allocate CDF arrays for non-master processes:
   do i = 1, Nat  ! for all atoms
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         if (allocated(matter%Atoms(i)%CDF)) then
            array_is_allocated = .true.
            N1 = size(matter%Atoms(i)%CDF)
         else
            array_is_allocated = .false.
            N1 = 0
         endif
      endif

      error_report = trim(adjustl(error_part))//' {array_is_allocated}'
      call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

      error_report = trim(adjustl(error_part))//' {N1}'
      call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

      if (array_is_allocated) then  ! allocate it also in non-master MPI processes:
         if (MPI_param%process_rank /= 0) then ! only do that for the master process
            if (.not.allocated(matter%Atoms(i)%CDF)) allocate(matter%Atoms(i)%CDF(N1))
         endif
      endif
   enddo ! i

   ! allocate the arrays in the CDF-objects:
   do i = 1, Nat   ! for all CDF functions
      Nsh = size(matter%Atoms(i)%Ip)
      do j = 1, Nsh  ! for all shells
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {CDF(i)%A}', matter%Atoms(i)%CDF(j)%A) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {CDF(i)%E0}', matter%Atoms(i)%CDF(j)%E0) ! below
         call broadcast_allocatable_array(MPI_param, trim(adjustl(error_part))//' {CDF(i)%G}', matter%Atoms(i)%CDF(j)%G) ! below
      enddo ! j
   enddo ! i
#endif
end subroutine MPI_share_Ritchi_CDF



!===============================================
! Basic MPI subroutines wrappers:


subroutine initialize_MPI(MPI_param, Err_data)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   type(Error_handling), intent(inout) :: Err_data  ! save data about error if any
   !-----------------------
   character(100) :: Error_message

#ifdef MPI_USED
   ! Initialize MPI:
   call MPI_INIT(MPI_param%ierror)
   if (MPI_param%ierror /= 0) then
      write(Error_message, *) 'Error initializing MPI!'
      call MPI_Save_error_details(Err_data, -1, Error_message, MPI_param)   ! above
      write(6, '(a)') trim(adjustl(Error_message))
      return
   endif

   ! Determine the size of the group associated with a communicator (cluster size, number of processes):
   call MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_param%size_of_cluster, MPI_param%ierror)
   if (MPI_param%ierror /= 0) then
      write(Error_message, *) 'Error getting MPI cluster size (number of processes)!'
      call MPI_Save_error_details(Err_data, -1, Error_message, MPI_param)   ! above
      write(6, '(a)') trim(adjustl(Error_message))
      return
   endif

   ! Determine the rank of the calling process in the communicator:
   call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_param%process_rank, MPI_param%ierror)
   if (MPI_param%ierror /= 0) then
      write(Error_message, *) 'Error getting MPI process rank!'
      call MPI_Save_error_details(Err_data, -1, Error_message, MPI_param)   ! bove
      write(6, '(a)') trim(adjustl(Error_message))
      return
   endif

   if (MPI_param%process_rank == 0) then ! only do that for the master process
      ! initialize MPI time counter:
      call get_MPI_lapsed_time(MPI_param%Wt0) ! below
   endif
#else
   ! No MPI, so only one process is there:
   MPI_param%process_rank = 0    ! index of the master process
   MPI_param%size_of_cluster = 1 ! total number of processes: 1 if no MPI is used
   MPI_param%ierror =0           ! error handler (no errors)
#endif
   write(MPI_param%rank_ch,'(i0)') MPI_param%process_rank
end subroutine initialize_MPI



subroutine initialize_random_seed(MPI_param)
   type(Used_MPI_parameters), intent(in) :: MPI_param
   !-----------------
   integer :: RN_seed
#ifdef MPI_USED
   ! Initialize different random seed for each process:
   CALL SYSTEM_CLOCK(count=RN_seed)
   RN_seed = RN_seed/100000 + MPI_param%process_Rank*100000
   call random_seed(put = (/RN_seed/) ) ! standard FORTRAN seeding of random numbers
#else
   ! Without MPI, use the default random seed:
   call random_seed() ! standard FORTRAN seeding of random numbers
#endif
end subroutine initialize_random_seed



subroutine MPI_share_add_data(numpar, Err)
   type(Numerics_param), intent(inout), target :: numpar 	! all numerical parameters
   type(Error_handling), intent(inout) :: Err      ! save data about error if any
   !--------------------
   character(100) :: error_part
   type(Used_MPI_parameters), pointer :: MPI_param

#ifdef MPI_USED

   MPI_param => numpar%MPI_param ! shorthand notation
   error_part = 'ERROR in MPI_share_add_data' ! part of the error message

   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%path_sep}', numpar%path_sep) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%change_size}', numpar%change_size) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%allow_rotate}', numpar%allow_rotate) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%verbose}', numpar%verbose) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {numpar%nonverbose}', numpar%nonverbose) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Err%Err}', Err%Err) ! below
   call broadcast_variable(MPI_param, trim(adjustl(error_part))//' {Err%Stopsignal}', Err%Stopsignal) ! below

   nullify(MPI_param)
#endif
end subroutine MPI_share_add_data




subroutine MPI_barrier_wrapper(MPI_param)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------

#ifdef MPI_USED
   call MPI_BARRIER( MPI_COMM_WORLD, MPI_param%ierror) ! module "MPI"
#endif
end subroutine MPI_barrier_wrapper



subroutine MPI_error_wrapper(process_rank, ierror, error_message)
   integer, intent(inout) :: process_rank, ierror
   character(*), intent(in) :: error_message   ! message to print about error
   !--------------------------
   if (ierror /= 0) then
        write(*,'(a,i0,a,i0)') '[MPI process #', process_rank, '] '//trim(adjustl(error_message)), ierror
        ! Cannot continue if the calculations are wrong:
#ifdef MPI_USED
        call MPI_Abort(MPI_COMM_WORLD, -1, ierror)   ! module "MPI"
#endif
    endif
end subroutine MPI_error_wrapper




subroutine MPI_fileopen_wrapper(MPI_param, File_name, FN, Error_message, readonly, err_msg, MPI_master)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: File_name
   integer, intent(inout) :: FN
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(in), optional :: readonly, MPI_master
   character(*), intent(in), optional :: err_msg   ! optional message to print together with error
   !------------------------
   logical :: if_read_only, MPI_master_only
   integer :: access_mode, ierr
   character(200) :: Error_descript
   !------------------------

   ! To chekc if the file is in read_only mode
   if (present(readonly)) then
      if_read_only = readonly
   else  ! by default, it is not
      if_read_only = .false.
   endif

   ! Initialize possible error message
   if (present(err_msg)) then
      Error_descript = trim(adjustl(err_msg))
   else
      Error_descript = ''  ! nothing yet
   endif

   ! To check if all MPI threads need access to the file, or only the master thread:
   if (present(MPI_master)) then
      MPI_master_only = MPI_master
   else
      MPI_master_only = .false.
   endif


   ! Opening file for (possibly) parallel i/o in it:
#ifdef MPI_USED
   ! https://rookiehpc.org/mpi/docs/mpi_file_open/index.html
   if (MPI_master_only) then ! only master thread opens the file:
      if (MPI_param%process_rank == 0) then  ! master thread has rank = 0 by default!
         call non_MPI_fileopen(trim(adjustl(File_name)), FN, Error_message, Error_descript, if_read_only, MPI_param) ! below
      endif
   else  ! open file for all MPI processes
      if (if_read_only) then
         access_mode = MPI_MODE_RDONLY ! With read-only access
      else
         access_mode = MPI_MODE_CREATE ! Create the file if it does not exist
         !access_mode = access_mode + MPI_MODE_EXCL ! The file must not exist, to avoid mistakenly erasing a file
         access_mode = access_mode + MPI_MODE_RDWR ! With read-write access
      endif

      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(adjustl(File_name)), access_mode, MPI_INFO_NULL, FN, MPI_param%ierror) ! module "MPI"

      if (MPI_param%ierror /= MPI_SUCCESS) then
         write(Error_descript, '(A,I0,A)') '[MPI process #', MPI_param%process_rank, '] Failure in opening file: '//trim(adjustl(File_name))
         call MPI_Save_error_details(Error_message, 1, Error_descript, MPI_param) ! above
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
         call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
      endif
   endif
#else
   call non_MPI_fileopen(trim(adjustl(File_name)), FN, Error_message, Error_descript, if_read_only, MPI_param) ! below
#endif
end subroutine MPI_fileopen_wrapper



subroutine non_MPI_fileopen(File_name, FN, Error_message, err_msg, read_only, MPI_param)
   character(*), intent(in) :: File_name
   integer, intent(inout) :: FN
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   character(*), intent(in), optional :: err_msg   ! optional message to print together with error
   logical, intent(in) :: read_only
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------
   character(200) :: Error_descript
   integer :: ierr
   !------------------

   ! Initialize possible error message
   if (present(err_msg)) then
      Error_descript = trim(adjustl(err_msg))
   else
      Error_descript = ''  ! nothing yet
   endif

   if (read_only) then  ! readonly option for existing file
      open(newunit=FN, FILE = trim(adjustl(File_name)), status = 'old', readonly, IOSTAT = ierr)
   else
      open(newunit=FN, FILE = trim(adjustl(File_name)), IOSTAT = ierr)
   endif

   if (ierr /= 0) then ! error opening the file
      write(Error_descript, '(A)') 'Failure in opening file: '//trim(adjustl(File_name))
      call MPI_Save_error_details(Error_message, 1, Error_descript, MPI_param) ! above
      print*, trim(adjustl(Error_descript)) ! print it also on the sreen
   endif
end subroutine non_MPI_fileopen



subroutine MPI_fileclose_wrapper(MPI_param, FN, Error_message, delete_file, err_msg, MPI_master)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   integer, intent(inout) :: FN
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(in), optional :: delete_file, MPI_master
   character(*), intent(in), optional :: err_msg   ! optional message to print together with error
   !-----------------
   logical :: file_to_be_deleted, MPI_master_only
   character(200) :: Error_descript

   ! Check if the file is to be deleted or not:
   if (present(delete_file)) then
      file_to_be_deleted = delete_file
   else  ! by default, don't delete it
      file_to_be_deleted = .false.
   endif

   ! Initialize possible error message
   if (present(err_msg)) then
      Error_descript = trim(adjustl(err_msg))
   else
      Error_descript = ''  ! nothing yet
   endif

   ! To check if all MPI threads need access to the file, or only the master thread:
   if (present(MPI_master)) then
      MPI_master_only = MPI_master
   else
      MPI_master_only = .false.
   endif

   ! Closing the opened file:
#ifdef MPI_USED
   ! https://www.open-mpi.org/doc/v3.0/man3/MPI_File_close.3.php
   if (MPI_master_only) then ! only master thread opens the file:
      if (MPI_param%process_rank == 0) then  ! master thread has rank = 0 by default!
         call non_MPI_fileclose(FN, file_to_be_deleted, Error_message, MPI_param=MPI_param)   ! below
      endif
   else  ! open file for all MPI processes
      call MPI_FILE_CLOSE(FN, MPI_param%ierror)    ! module "MPI"
      if (MPI_param%ierror /= MPI_SUCCESS) then
         write(Error_descript, '(A,I0,A)') '[MPI process #', MPI_param%process_rank, '] Failure in closing file'
         call MPI_Save_error_details(Error_message, 1, Error_descript, MPI_param) ! module "Objects"
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
         call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
      endif
   endif
#else
   call non_MPI_fileclose(FN, file_to_be_deleted, Error_message, MPI_param=MPI_param)   ! below
#endif
end subroutine MPI_fileclose_wrapper


subroutine non_MPI_fileclose(FN, delete_file, Error_message, err_msg, MPI_param)
   integer, intent(in) :: FN
   logical, intent(in) :: delete_file
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   character(*), intent(in), optional :: err_msg   ! optional message to print together with error
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !---------------------
   logical :: file_opened
   integer :: ierr
   character(200) :: Error_descript
   !---------------------

      ! Initialize possible error message
   if (present(err_msg)) then
      Error_descript = trim(adjustl(err_msg))
   else
      Error_descript = ''  ! nothing yet
   endif

   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) then
      if (delete_file) then
         close(FN, status='delete', IOSTAT = ierr)
      else
         close(FN, IOSTAT = ierr)
      endif

      if (ierr /= 0) then ! error opening the file
         write(Error_descript, '(A)') trim(adjustl(Error_descript))//'Failure in closing file'
         call MPI_Save_error_details(Error_message, 1, Error_descript, MPI_param) ! module "Objects"
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
      endif
   endif
end subroutine non_MPI_fileclose




subroutine get_MPI_lapsed_time(last_time, current_time, lasped_time, printout_text)
    real(8), intent(inout) :: last_time       ! previous time point
    real(8), intent(out), optional :: current_time    ! currect time point
    real(8), intent(out), optional :: lasped_time       ! lapsed time to be calculated
    character(*), intent(in), optional :: printout_text ! text to printout together with the lapsed time
    !-------------------
    real(8) :: a_lasped_time, a_current_time
    character(100) :: time_duration_string

#ifdef MPI_USED
    ! Define the current time point:
    if ( present(current_time) .or. (present(lasped_time)) .or. (present(printout_text)) ) then
        ! Get the current time:
        a_current_time = MPI_Wtime()    ! module "MPI"
        if (present(current_time)) current_time = a_current_time    ! printout, if requested

        ! Get the lapsed time:
        a_lasped_time = a_current_time - last_time
        if (present(lasped_time)) lasped_time = a_lasped_time   ! printout, if requested

        ! If user requests, printout the message about the time lapsed:
        if (present(printout_text)) then
            call pars_MPI_lasped_time(a_lasped_time, time_duration_string)   ! below
            print*, trim(adjustl(printout_text))//' '//trim(adjustl(time_duration_string))
        endif
    else    ! initialization: only one variable provided, save time point into it:
        last_time = MPI_Wtime()    ! module "MPI"
    endif
#endif
end subroutine get_MPI_lapsed_time



subroutine pars_MPI_lasped_time(sec, time_string)
   real(8), intent(inout) :: sec ! time interval in [sec]
   character(*), intent(out) :: time_string ! split it into mins, hours, days...
   !-----------------------
   character(100) :: temp
   real(8) :: days, hours, mins, msec
   days = 0.0d0     ! to start with
   hours = 0.0d0    ! to start with
   mins = 0.0d0     ! to start with
   msec = 0.0d0     ! to start with

   if (sec < 1.0d0) then    ! msec
      msec = sec * 1.0d3
   else if (sec .GE. 60.0d0) then   ! there are minutes
      mins = FLOOR(sec/60.0d0)  ! minutes
      sec = sec - mins*60.0d0   ! update seconds
      if (mins .GT. 60.0d0) then    ! there are hours
         hours = FLOOR(mins/60.0d0) ! hours
         mins = mins - hours*60.0d0 ! update minutes
         if (hours .GT. 24.0d0) then    ! there are days
            days = FLOOR(hours/24.0d0)  ! days
            hours = hours - days*24.0d0 ! hourse
         endif
      endif
   endif
   time_string = '' ! to start with
   temp = ''        ! to start with

   ! Write #days in the string
   if (days .GT. 1.0d0) then
      write(temp, '(i9)') int(days)
      write(time_string, '(a,a)') trim(adjustl(temp)), ' days'
   else if (days .GT. 0.0d0) then
      write(temp, '(i9)') int(days)
      write(time_string, '(a,a)') trim(adjustl(temp)), ' day'
   endif

   ! Write #hours in the string
   if (hours .GT. 1.0d0) then
      write(temp, '(i9)') int(hours)
      write(time_string, '(a,a,a)') trim(adjustl(time_string)), ' '//trim(adjustl(temp)), ' hours'
   else if (hours .GT. 0.0d0) then
      write(temp, '(i9)') int(hours)
      write(time_string, '(a,a,a)') trim(adjustl(time_string)), ' '//trim(adjustl(temp)), ' hour'
   endif

   ! Write #minutes in the string
   if (mins .GT. 1.0d0) then
      write(temp, '(i9)') int(mins)
      write(time_string, '(a,a,a)') trim(adjustl(time_string)), ' '//trim(adjustl(temp)), ' mins'
   else if (mins .GT. 0.0d0) then
      write(temp, '(i9)') int(mins)
      write(time_string, '(a,a,a)') trim(adjustl(time_string)), ' '//trim(adjustl(temp)), ' min'
   endif

   ! Write msec in the string
   if (msec > 1.0d-10) then
      write(temp, '(f15.6)') msec
      write(time_string, '(a,a,a)') trim(adjustl(time_string)), ' '//trim(adjustl(temp)), ' msec'
   else
      write(temp, '(f15.6)') sec
      write(time_string, '(a,a,a)') trim(adjustl(time_string)), ' '//trim(adjustl(temp)), ' sec'
   endif
end subroutine pars_MPI_lasped_time



subroutine MPI_Save_error_details(Err_name, Err_num, Err_data, MPI_param)
   class(Error_handling) :: Err_name    ! object containing all details
   integer, intent(in) :: Err_num       ! number of error asigned
   character(*), intent(in) :: Err_data   ! description of the error
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !---------------------
   integer :: FN, Nsiz   ! number of file where to write error log

   if (MPI_param%process_Rank == 0) then ! it is a master thread, it has access to the file
      call Save_error_details(Err_name, Err_num, Err_data)  ! module "Objects"
      return   ! that is it, wrote into the file, nothing else to do
   else  ! it is not a master thread, send info to the master to write into a file
#ifdef MPI_USED
      FN = Err_name%File_Num   ! this number is provided in the Err_name object
      Err_name%Err = .true.    ! error occured, mark it as "true"
      Err_name%Err_Num = Err_num   ! number of error we asign to it
      Err_name%Err_descript = Err_data ! descriptino of an error

      if (MPI_param%process_Rank /= 0) then  ! non-master process
         call MPI_SEND(FN, 1, MPI_INTEGER, 0, MPI_param%process_Rank, MPI_COMM_WORLD, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error issuing send request {Save_error_details:FN} from process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
         call MPI_SEND(Err_name%Err, 1, MPI_LOGICAL, 0, MPI_param%process_Rank, MPI_COMM_WORLD, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error issuing send request {Save_error_details:Err} from process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
         call MPI_SEND(Err_name%Err_Num, 1, MPI_INTEGER, 0, MPI_param%process_Rank, MPI_COMM_WORLD, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error issuing send request {Save_error_details:Err_Num}  from process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
         Nsiz = LEN(Err_name%Err_descript)
         call MPI_SEND(Nsiz, 1, MPI_INTEGER, 0, MPI_param%process_Rank, MPI_COMM_WORLD, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error issuing send request {Save_error_details:Nsiz}  from process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
      else  ! master process
         call MPI_RECV(FN, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error in receiving request {Save_error_details:FN} by process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
         call MPI_RECV(Err_name%Err, 1, MPI_LOGICAL, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error in receiving request {Save_error_details:Err} by process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
         call MPI_RECV(Err_name%Err_Num, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error in receiving request {Save_error_details:Err_Num} by process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
         call MPI_RECV(Nsiz, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error in receiving request {Save_error_details:Nsiz} by process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
      endif

      !------------------------------------------------------
      ! Synchronize MPI processes: make sure master process got the length of the message, to recieve all the text
      call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
      !------------------------------------------------------

      if (MPI_param%process_Rank /= 0) then  ! non-master process
         call MPI_SEND(Err_name%Err_descript, Nsiz, MPI_CHARACTER, 0, MPI_param%process_Rank, MPI_COMM_WORLD, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error issuing send request {Save_error_details:Err_descript} from process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
      else
         call MPI_RECV(Err_name%Err_descript, Nsiz, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, MPI_param%ierror)
         if (MPI_param%ierror /= 0) then
            write(*, *) 'Error in receiving request {Save_error_details:Err_descript} by process: ', MPI_param%process_Rank, MPI_param%ierror
            call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
         endif
      endif

      write(FN, '(a,i2,1x,a)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript))   ! write it all into the file
#else ! use the nonMPI version of the error saving
      call Save_error_details(Err_name, Err_num, Err_data)  ! module "Objects"
#endif
   endif
end subroutine MPI_Save_error_details





!===============================================
! Broadcast wrappers:

subroutine broadcast_variable_int(MPI_param, error_message, var_int)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   integer, intent(inout) :: var_int
   !---------------------------
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {integer}'
   call mpi_bcast(var_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_variable_int



subroutine broadcast_variable_logic(MPI_param, error_message, var_logic)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   logical, intent(inout) :: var_logic
   !---------------------------
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {logical}'
   call mpi_bcast(var_logic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_variable_logic


subroutine broadcast_variable_real(MPI_param, error_message, var_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   real(8), intent(inout) :: var_real
   !---------------------------
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {real}'
   call mpi_bcast(var_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_variable_real


subroutine broadcast_variable_complex(MPI_param, error_message, var_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   complex, intent(inout) :: var_real
   !---------------------------
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {complex}'
   call mpi_bcast(var_real, 1, MPI_COMPLEX, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_variable_complex



subroutine broadcast_variable_char(MPI_param, error_message, var_char)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   character(*), intent(inout) :: var_char
   !---------------------------
   integer :: Nsiz
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {character}'
   Nsiz = LEN(var_char)
   call mpi_bcast(var_char, Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_variable_char



subroutine broadcast_array_logic(MPI_param, error_message, array_logic)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   logical, dimension(:), intent(inout) :: array_logic
   !---------------------------
   integer :: Nsiz
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {logical_array}'
   Nsiz = size(array_logic)
   call mpi_bcast(array_logic, Nsiz, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_array_logic



subroutine broadcast_array_int(MPI_param, error_message, array_int)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   integer, dimension(:), intent(inout), optional :: array_int
   !---------------------------
   integer :: Nsiz
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {integer_array}'
   Nsiz = size(array_int)
   call mpi_bcast(array_int, Nsiz, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_array_int



subroutine broadcast_array_real(MPI_param, error_message, array_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   real(8), dimension(:), intent(inout) :: array_real
   !---------------------------
   integer :: Nsiz
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {real_array}'
   Nsiz = size(array_real)
   call mpi_bcast(array_real, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_array_real



subroutine broadcast_2d_array_real(MPI_param, error_message, array_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   real(8), dimension(:,:), intent(inout) :: array_real
   !---------------------------
   integer :: N1, N2
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {real_2d_array}'
   N1 = size(array_real,1)
   N2 = size(array_real,2)
   call mpi_bcast(array_real, N1*N2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_2d_array_real



subroutine broadcast_2d_array_complex(MPI_param, error_message, array_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   complex, dimension(:,:), intent(inout) :: array_real
   !---------------------------
   integer :: N1, N2
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {complex_2d_array}'
   N1 = size(array_real,1)
   N2 = size(array_real,2)
   call mpi_bcast(array_real, N1*N2, MPI_COMPLEX, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_2d_array_complex



subroutine broadcast_3d_array_real(MPI_param, error_message, array_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   real(8), dimension(:,:,:), intent(inout) :: array_real
   !---------------------------
   integer :: N1, N2, N3
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {real_3d_array}'
   N1 = size(array_real,1)
   N2 = size(array_real,2)
   N3 = size(array_real,3)
   call mpi_bcast(array_real, N1*N2*N3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_3d_array_real



subroutine broadcast_4d_array_real(MPI_param, error_message, array_real)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   real(8), dimension(:,:,:,:), intent(inout) :: array_real
   !---------------------------
   integer :: N1, N2, N3, N4
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {real_4d_array}'
   N1 = size(array_real,1)
   N2 = size(array_real,2)
   N3 = size(array_real,3)
   N4 = size(array_real,4)
   call mpi_bcast(array_real, N1*N2*N3*N4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_4d_array_real




subroutine broadcast_array_char(MPI_param, error_message, array_char)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   ! Non-allocatable, or already allocated, arrays only:
   character(*), dimension(:), intent(inout) :: array_char
   !---------------------------
   integer :: Nsiz, N1
   character(300) :: error_report

#ifdef MPI_USED
   error_report = trim(adjustl(error_message))//' {character_array}'
   Nsiz = LEN(array_char(1))
   N1 = size(array_char)
   call mpi_bcast(array_char, N1*Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report)) ) ! module "MPI_subroutines"
#endif
end subroutine broadcast_array_char




subroutine broadcast_allocatable_char_1d_array(MPI_param, error_message, array)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   character(*), dimension(:), allocatable, intent(inout) :: array
   !--------------------------
   logical :: array_is_allocated
   integer :: N1, Nsiz
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array)
      else
         N1 = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N1}'
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1))
      endif

      if (size(array) > 0) then
         Nsiz = LEN(array(1))
      else
         Nsiz = 0
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1*Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif

#endif
end subroutine broadcast_allocatable_char_1d_array




subroutine broadcast_allocatable_logic_1d_array(MPI_param, error_message, array)
   logical, dimension(:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array)
      else
         N1 = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N1}'
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif

#endif
end subroutine broadcast_allocatable_logic_1d_array



subroutine broadcast_allocatable_logic_3d_array(MPI_param, error_message, array)
   logical, dimension(:,:,:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1, N2, N3, N3d(3)
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array,1)
         N2 = size(array,2)
         N3 = size(array,3)
      else
         N1 = 0
         N2 = 0
         N3 = 0
      endif
      N3d = (/N1,N2,N3/)
   endif


   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N3d}'
   call mpi_bcast(N3d, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   N1 = N3d(1)
   N2 = N3d(2)
   N3 = N3d(3)

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1, N2, N3))
      endif

      error_report = trim(adjustl(error_message))//' {array(logical)}'
      call mpi_bcast(array, N1*N2*N3, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif
#endif
end subroutine broadcast_allocatable_logic_3d_array




subroutine broadcast_allocatable_int_1d_array(MPI_param, error_message, array)
   integer, dimension(:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array)
      else
         N1 = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N1}'
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif

#endif
end subroutine broadcast_allocatable_int_1d_array




subroutine broadcast_allocatable_int_2d_array(MPI_param, error_message, array)
   integer, dimension(:,:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1, N2
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array,1)
         N2 = size(array,2)
      else
         N1 = 0
         N2 = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N1}'
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N2}'
   call mpi_bcast(N2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1,N2))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1*N2, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif

#endif
end subroutine broadcast_allocatable_int_2d_array



subroutine broadcast_allocatable_real_1d_array(MPI_param, error_message, array)
   real(8), dimension(:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
         N1 = size(array)
      else
         array_is_allocated = .false.
         N1 = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N1}'
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif
#endif
end subroutine broadcast_allocatable_real_1d_array




subroutine broadcast_allocatable_real_2d_array(MPI_param, error_message, array)
   real(8), dimension(:,:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1, N2
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array,1)
         N2 = size(array,2)
      else
         N1 = 0
         N2 = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N1}'
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N2}'
   call mpi_bcast(N2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1, N2))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1*N2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif
#endif
end subroutine broadcast_allocatable_real_2d_array




subroutine broadcast_allocatable_real_3d_array(MPI_param, error_message, array)
   real(8), dimension(:,:,:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N1, N2, N3, N3d(3)
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N1 = size(array,1)
         N2 = size(array,2)
         N3 = size(array,3)
      else
         N1 = 0
         N2 = 0
         N3 = 0
      endif
      N3d = (/N1,N2,N3/)
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N3d}'
   call mpi_bcast(N3d, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   N1 = N3d(1)
   N2 = N3d(2)
   N3 = N3d(3)

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N1, N2, N3))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N1*N2*N3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif
#endif
end subroutine broadcast_allocatable_real_3d_array




subroutine broadcast_allocatable_real_4d_array(MPI_param, error_message, array)
   real(8), dimension(:,:,:,:), allocatable, intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   logical :: array_is_allocated
   integer :: N4d(4)
   character(300) :: error_report

#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (allocated(array)) then
         array_is_allocated = .true.
      else
         array_is_allocated = .false.
      endif
      if (array_is_allocated) then
         N4d(1) = size(array,1)
         N4d(2) = size(array,2)
         N4d(3) = size(array,3)
         N4d(4) = size(array,4)
      else
         N4d = 0
      endif
   endif

   error_report = trim(adjustl(error_message))//' {array_is_allocated}'
   call mpi_bcast(array_is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   error_report = trim(adjustl(error_message))//' {N3d}'
   call mpi_bcast(N4d, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"

   if (array_is_allocated) then
      if (MPI_param%process_rank /= 0) then   ! also allocate it in all other processes
         if (.not.allocated(array)) allocate(array(N4d(1), N4d(2), N4d(3), N4d(4)))
      endif

      error_report = trim(adjustl(error_message))//' {array}'
      call mpi_bcast(array, N4d(1)*N4d(2)*N4d(3)*N4d(4), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
   endif
#endif
end subroutine broadcast_allocatable_real_4d_array



subroutine do_MPI_Allreduce_real_variable(MPI_param, error_message, var)
   real(8), intent(inout) :: var
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   character(300) :: error_report
#ifdef MPI_USED
   ! https://rookiehpc.org/mpi/docs/mpi_allreduce/index.html
   CALL MPI_Allreduce(MPI_IN_PLACE, var, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"

   error_report = trim(adjustl(error_message))//' {Allreduce_real_variable}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Allreduce_real_variable



subroutine do_MPI_Allreduce_real_1d_array(MPI_param, error_message, array)
   real(8), dimension(:), intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   integer :: Nsiz
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz = size(array)

   CALL MPI_Allreduce(MPI_IN_PLACE, array, Nsiz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"

   error_report = trim(adjustl(error_message))//' {Allreduce_real_1d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Allreduce_real_1d_array



subroutine do_MPI_Allreduce_real_2d_array(MPI_param, error_message, array)
   real(8), dimension(:,:), intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   integer :: Nsiz(2)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)

   CALL MPI_Allreduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"

   error_report = trim(adjustl(error_message))//' {Allreduce_real_2d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Allreduce_real_2d_array



subroutine do_MPI_Allreduce_complex_2d_array(MPI_param, error_message, array)
   complex, dimension(:,:), intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   integer :: Nsiz(2)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)

   CALL MPI_Allreduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2), MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"

   error_report = trim(adjustl(error_message))//' {Allreduce_complex_2d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Allreduce_complex_2d_array


subroutine do_MPI_Allreduce_real_3d_array(MPI_param, error_message, array)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   real(8), dimension(:,:,:), intent(inout) :: array
   !--------------------------
   integer :: Nsiz(3)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)
   Nsiz(3) = size(array,3)

   CALL MPI_Allreduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2)*Nsiz(3), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"

   error_report = trim(adjustl(error_message))//' {Allreduce_real_3d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Allreduce_real_3d_array


subroutine do_MPI_Allreduce_real_4d_array(MPI_param, error_message, array)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   real(8), dimension(:,:,:,:), intent(inout) :: array
   !--------------------------
   integer :: Nsiz(4)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)
   Nsiz(3) = size(array,3)
   Nsiz(4) = size(array,4)

   CALL MPI_Allreduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2)*Nsiz(3)*Nsiz(4), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"

   error_report = trim(adjustl(error_message))//' {Allreduce_real_4d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Allreduce_real_4d_array




subroutine do_MPI_Reduce_real_variable(MPI_param, error_message, var)
   real(8), intent(inout) :: var
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   character(300) :: error_report
#ifdef MPI_USED
   if (MPI_param%process_rank == 0) then   ! MPI master process gets the value from itself, hence MPI_IN_PLACE
      call MPI_Reduce(MPI_IN_PLACE, var, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   else
      call MPI_Reduce(var, var, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   endif

   error_report = trim(adjustl(error_message))//' {Reduce_real_variable}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Reduce_real_variable



subroutine do_MPI_Reduce_real_1d_array(MPI_param, error_message, array)
   real(8), dimension(:), intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   integer :: Nsiz
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz = size(array)

   if (MPI_param%process_rank == 0) then   ! MPI master process gets the value from itself, hence MPI_IN_PLACE
      call MPI_Reduce(MPI_IN_PLACE, array, Nsiz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   else
      call MPI_Reduce(array, array, Nsiz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   endif

   error_report = trim(adjustl(error_message))//' {Reduce_real_1d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Reduce_real_1d_array



subroutine do_MPI_Reduce_real_2d_array(MPI_param, error_message, array)
   real(8), dimension(:,:), intent(inout) :: array
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   !--------------------------
   integer :: Nsiz(2)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)

   if (MPI_param%process_rank == 0) then   ! MPI master process gets the value from itself, hence MPI_IN_PLACE
      call MPI_Reduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   else
      call MPI_Reduce(array, array, Nsiz(1)*Nsiz(2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   endif

   error_report = trim(adjustl(error_message))//' {Reduce_real_2d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Reduce_real_2d_array



subroutine do_MPI_Reduce_real_3d_array(MPI_param, error_message, array)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   real(8), dimension(:,:,:), intent(inout) :: array
   !--------------------------
   integer :: Nsiz(3)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)
   Nsiz(3) = size(array,3)

   if (MPI_param%process_rank == 0) then   ! MPI master process gets the value from itself, hence MPI_IN_PLACE
      call MPI_Reduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2)*Nsiz(3), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   else
      call MPI_Reduce(array, array, Nsiz(1)*Nsiz(2)*Nsiz(3), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   endif

   error_report = trim(adjustl(error_message))//' {Reduce_real_3d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Reduce_real_3d_array



subroutine do_MPI_Reduce_real_4d_array(MPI_param, error_message, array)
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   character(*), intent(in) :: error_message
   real(8), dimension(:,:,:,:), intent(inout) :: array
   !--------------------------
   integer :: Nsiz(4)
   character(300) :: error_report
#ifdef MPI_USED
   Nsiz(1) = size(array,1)
   Nsiz(2) = size(array,2)
   Nsiz(3) = size(array,3)
   Nsiz(4) = size(array,4)

   if (MPI_param%process_rank == 0) then   ! MPI master process gets the value from itself, hence MPI_IN_PLACE
      call MPI_Reduce(MPI_IN_PLACE, array, Nsiz(1)*Nsiz(2)*Nsiz(3)*Nsiz(4), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   else
      call MPI_Reduce(array, array, Nsiz(1)*Nsiz(2)*Nsiz(3)*Nsiz(4), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, MPI_param%ierror)   ! module "mpi"
   endif

   error_report = trim(adjustl(error_message))//' {Reduce_real_4d_array}'
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_report))) ! module "MPI_subroutines"
#endif
end subroutine do_MPI_Reduce_real_4d_array





end module MPI_subroutines
