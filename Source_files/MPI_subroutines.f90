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
! This module contains subroutines to deal with MPI calls,
! such as wrappers and MPI derived types for particular purposes
!**************************************************************


module MPI_subroutines

#ifdef MPI_USED
use mpi
#endif

use Objects, only : Error_handling, Used_MPI_parameters, Save_error_details

implicit none
PRIVATE  ! hides items not listed on public statement

public :: get_MPI_lapsed_time, initialize_MPI, initialize_random_seed, MPI_barrier_wrapper, MPI_fileopen_wrapper, &
            MPI_fileclose_wrapper, MPI_error_wrapper


contains


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






subroutine initialize_random_seed(MPI_param)
   type(Used_MPI_parameters), intent(in) :: MPI_param
   !-----------------
   integer :: RN_seed
   ! Initialize different random seed for each process:
#ifdef MPI_USED
   CALL SYSTEM_CLOCK(count=RN_seed)
   RN_seed = RN_seed/100000 + MPI_param%process_Rank*100000
   call random_seed(put = (/RN_seed/) ) ! standard FORTRAN seeding of random numbers
#else
   ! Without MPI, use the default random seed:
   call random_seed() ! standard FORTRAN seeding of random numbers
#endif
end subroutine initialize_random_seed




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


end module MPI_subroutines
