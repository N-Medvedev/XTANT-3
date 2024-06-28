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
! This module contains subroutines to deal with cross-sections.
! References used in the module:
! [1]  F. Salvat, J. M. Fernandez-Varea, E. Acosta, J. Sempau
!   "PENELOPE-2014 A Code System for Monte Carlo Simulation of Electron and Photon Transport", OECD (2014)
! [2] R. Rymzhanov et al. Phys. Status Solidi B 252, 159-164 (2015) / DOI 10.1002/pssb.201400130
! [3] M. Azzolini et al. J. Phys. Condens. Matter 31, 055901 (2019)
! [4] N. Medvedev et al., Advanced Theory and Simulations 5, 2200091 (2022)



MODULE MC_cross_sections
use Objects
use Universal_constants
use Little_subroutines, only : Find_in_array_monoton, linear_interpolation, print_progress, &
                              interpolate_data_on_grid, find_order_of_number, exclude_doubles
use Dealing_with_files, only : Count_lines_in_file, read_file
use Dealing_with_EADL, only : m_EPDL_file, READ_EPDL_TYPE_FILE_real, next_designator, Read_EPDL_data
use Algebra_tools, only : sort_array

implicit none
PRIVATE


real(8), parameter :: m_two_third = 2.0d0/3.0d0


public :: Mean_free_path, which_atom, which_shell, NRG_transfer_elastic_atomic, Electron_energy_transfer_inelastic
public :: get_MFPs, Get_photon_attenuation, TotIMFP, velosity_from_kinetic_energy


 contains

!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
!General functions:

subroutine Mean_free_path(E, mfps, MFP_cur, inversed) ! finds total mean free path in [A]
! from arrays, when they are already precalculated
   REAL(8), INTENT(in) ::  E      ! energy of the traced particle [eV]
   type(MFP), intent(in) :: mfps	! electron mean free paths for each shell
   REAL(8), INTENT(out) ::  MFP_cur   ! mean free path of the traced particle [A]
   logical, optional, intent(in) :: inversed ! do we need it streight or inversed?
   real(8) Ecur, MFP_temp
   integer i,j, N_temmp, Msh, N_last
   call Find_in_array_monoton(mfps%E, E, N_temmp) ! module "Little_subroutines"
   call linear_interpolation(mfps%E, mfps%L, E, MFP_temp, N_temmp)	! module "Little_subroutines"
   if (.not.present(inversed)) then
      if (MFP_temp < 1.0d-12) then ! avoid divide by zero
         MFP_cur = 1.3d21  ! infinity
      else
         MFP_cur = 1.0d0/MFP_temp ! [A] inverse it back
      endif
   else
      MFP_cur = MFP_temp ! give inverse MFP as output [1/A]
   endif
end subroutine Mean_free_path



subroutine which_atom(E, Atoms, EMFP_tot, i_at) ! find which atom from the compound electron scatters on
   real(8), intent(in) :: E ! [eV] for which kinetic energy we are looking up the MFP
   type(At_data), dimension(:), intent(in) :: Atoms	! all kinds of atoms of the compound (matter%Atoms)
   real(8), intent(in) :: EMFP_tot ! [1/A] total elastic mean free path
   integer, intent(out) :: i_at ! number of atomic species which an electron is scattering on in this event
   integer :: i
   real(8) :: RN, Lsum, Lsampled, EMFP
   
   if (size(Atoms) > 1) then ! there are more than 1, lets sample it:
      call random_number(RN)
      Lsampled = RN*EMFP_tot ! [1/A]
      Lsum = 0.0d0 ! summed MFP, to start counting
      do i = 1, size(Atoms) ! for each atomic species
         call Mean_free_path(E, Atoms(i)%El_EMFP, EMFP, inversed = .true.) ! finds inverse mean free path in [1/A]
         Lsum = Lsum + EMFP ! sum them up until we find the sampled one, [1/A] elastic MFP
!          print*, 'which_atom', i, Lsum, RN, EMFP, EMFP_tot
         if (Lsum >= Lsampled) exit
      enddo
      if (i > size(Atoms)) i = size(Atoms) ! just in case of numerical precision troubles
   else ! there is only one, no need to sample:
      i = 1
   endif
   i_at = i ! this is the kind of atom that electron scatters on now
end subroutine which_atom


subroutine which_shell(E, Atoms, matter, which_mfp, KOA, SHL)
   real(8), intent(in) :: E ! [eV] for which energy we are looking up the MFP
   type(At_data), dimension(:), intent(in) :: Atoms	! all kinds of atoms of the compound (matter%Atoms)
   type(Solid), intent(in) :: matter              ! all material parameters
   integer, intent(in) :: which_mfp ! 0=photons, 1=electrons
   integer, intent(out) :: KOA, SHL ! kind of atom; number of shell
   !===============================
   real(8) invL, Ltot, Lcur, RN
   integer i, j, N

   ! Total mean free path:
   select case (which_mfp)
   case (0) ! photon
      call Find_in_array_monoton(matter%Ph_MFP_tot%E, E, N) ! module "Little_subroutines"
      Ltot = matter%Ph_MFP_tot%L(N) ! inverse mean free path [1/A]
   case (1) ! electron
      call Find_in_array_monoton(matter%El_MFP_tot%E, E, N) ! module "Little_subroutines"
      Ltot = matter%El_MFP_tot%L(N) ! inverse mean free path [1/A]
   end select

   invL = 0.0d0 ! start with it
   call random_number(RN)

   AT_NUM:do i = 1, size(Atoms) ! for each kind of atoms:
      SH_NUM:do j = 1, Atoms(i)%sh ! for all shells:
         INF_L:if (Ltot .GE. 1d-20) then ! only do for the case when it's not infinite
            select case (which_mfp)
            case (0) ! photon
               POSIBL:if (E > matter%Atoms(i)%Ip(j)) then ! ionization of this shell is possible: 
                  ! Find the closest values in the array of mean free paths
                  call Find_in_array_monoton(Atoms(i)%Ph_MFP(j)%E, E, N) ! module "Little_subroutines"
                  ! Interpolate to the actual value:
                  if (N > 1) then
                     if (Atoms(i)%Ph_MFP(j)%L(N-1) < 1.0d-15) then
                        call linear_interpolation(Atoms(i)%Ph_MFP(j)%E, Atoms(i)%Ph_MFP(j)%L, E, Lcur, N, x0=matter%Atoms(i)%Ip(j), y0=0.0d0, replac=.true.) ! module "Little_subroutines" 
                     else
                        call linear_interpolation(Atoms(i)%Ph_MFP(j)%E, Atoms(i)%Ph_MFP(j)%L, E, Lcur, N, x0=matter%Atoms(i)%Ip(j), y0=0.0d0) ! module "Little_subroutines" 
                     endif
                  else 
                     Lcur = 0.0d0
                  endif
               else POSIBL
                  Lcur = 0.0d0
               endif POSIBL
               invL = invL + Lcur/Ltot  ! inverse the  inverse mean free path [A]
!                print*, 'which_shell photon', E, Atoms(i)%Ip(j), Atoms(i)%Ph_MFP(j)%L(N), N
!                if (N > 1) print*, Atoms(i)%Ph_MFP(j)%L(N-1), 'LAST'
!                print*, Lcur, Ltot, invL, RN
            case (1) ! electron
               POSIBL2:if (E > matter%Atoms(i)%Ip(j)) then ! ionization of this shell is possible: 
                  ! Find the closest values in the array of mean free paths
                  call Find_in_array_monoton(Atoms(i)%El_MFP(j)%E, E, N) ! module "Little_subroutines"
                  ! Interpolate to the actual value:
                  if ((E > matter%Atoms(i)%Ip(j)) .and. (N > 1)) then
                     if (Atoms(i)%Ph_MFP(j)%L(N-1) < 1.0d-15) then
!                         print*, 'which_shell 1'
                        call linear_interpolation(Atoms(i)%El_MFP(j)%E, Atoms(i)%El_MFP(j)%L, E, Lcur, N, x0=matter%Atoms(i)%Ip(j), y0=0.0d0, replac=.true.) ! module "Little_subroutines"
                     else
                        call linear_interpolation(Atoms(i)%El_MFP(j)%E, Atoms(i)%El_MFP(j)%L, E, Lcur, N, x0=matter%Atoms(i)%Ip(j), y0=0.0d0) ! module "Little_subroutines" 
                     endif
                     !invL = invL + Atoms(i)%El_MFP(j)%L(N)/Ltot  ! inverse the inverse mean free path [A]
                  else
                     Lcur = 0.0d0
                     !invL = invL + Atoms(i)%El_MFP(j)%L(N)/Ltot  ! inverse the inverse mean free path [A]
                  endif
               else POSIBL2
                  Lcur = 0.0d0
               endif POSIBL2
               invL = invL + Lcur/Ltot  ! inverse the inverse mean free path [A]

!                call linear_interpolation(Atoms(i)%El_MFP(j)%E, Atoms(i)%El_MFP(j)%L, E, Lcur, N, x0=Atoms(i)%Ip(j), y0=0.0d0) ! module "Little_subroutines"
               !write(*,'(a,i1,es,es,f,f)') 'shell ', j, Atoms(i)%El_MFP(j)%L(N), Ltot, Atoms(i)%El_MFP(j)%L(N)/Ltot, RN
            end select
            if (invL .GE. RN) then
               KOA = i
               SHL = j
               exit AT_NUM
            else ! no shell of an atom could be ionizaed, assume VB:
               KOA = 1
               SHL = Atoms(1)%sh
            endif
         else INF_L ! no shell of an atom could be ionized for such energy, use VB:
            KOA = 1
            SHL = Atoms(1)%sh
         endif INF_L
      enddo SH_NUM
   enddo AT_NUM
end subroutine which_shell


!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Electrons
subroutine get_MFPs(Scell, NSC, matter, laser, numpar, TeeV, Err)
   type(Super_cell), dimension(:), intent(in) :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! supercell index
   type(Solid), intent(inout) :: matter ! parameters of the material
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar  ! all numerical parameters
   real(8), intent(in) :: TeeV   ! [eV] electronic tempreature
   type(Error_handling), intent(inout) :: Err	! error save
   !===============================================
   integer N_grid, N, FN, INFO
   character(300) :: File_name, Error_descript
   character(10) :: chtemp, chtemp2
   logical :: redo, file_opened, file_exists, read_well
   integer i, j, k, Nshl, my_id, Reason, count_lines, N_Te_points
   real(8) :: Ele, L, dEdx, Omega, ksum, fsum, Te_temp
   integer OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
   
   ! Set grid for MFPs:
   call get_grid_4CS(N_grid, maxval(laser(:)%hw), matter%Atoms(1)%El_MFP(1)%E, matter) ! below

   if (.not.allocated(matter%Atoms(1)%El_MFP(1)%L)) allocate(matter%Atoms(1)%El_MFP(1)%L(N_grid))

   ! Temperature dependence (for the valence/conduction band only):
   if (.not. allocated(matter%Atoms(1)%El_MFP_vs_T)) then
      N_Te_points = 51  ! points in Te, by 1000K step
      allocate(matter%Atoms(1)%El_MFP_vs_T(N_Te_points))
      do i = 1, N_Te_points
         allocate(matter%Atoms(1)%El_MFP_vs_T(i)%L(N_grid))
         allocate(matter%Atoms(1)%El_MFP_vs_T(i)%E(N_grid))
      enddo
   endif

   ! Total inelastic mean free paths:
   if (.not.allocated(matter%El_MFP_tot%L)) allocate(matter%El_MFP_tot%L(N_grid))
   if (.not.allocated(matter%El_MFP_tot%E)) allocate(matter%El_MFP_tot%E(N_grid))
   matter%El_MFP_tot%E = matter%Atoms(1)%El_MFP(1)%E
   matter%El_MFP_tot%L = 0.0d0
   ! Elastic mean free paths:
   if (.not.allocated(matter%Atoms(1)%El_EMFP%E)) allocate(matter%Atoms(1)%El_EMFP%E(N_grid))
   if (.not.allocated(matter%Atoms(1)%El_EMFP%L)) allocate(matter%Atoms(1)%El_EMFP%L(N_grid))
   if (.not.allocated(matter%El_EMFP_tot%L)) allocate(matter%El_EMFP_tot%L(N_grid))
   if (.not.allocated(matter%El_EMFP_tot%E)) allocate(matter%El_EMFP_tot%E(N_grid))
   matter%Atoms(1)%El_EMFP%E = matter%Atoms(1)%El_MFP(1)%E ! use the same grid as inelastic
   matter%El_EMFP_tot%E = matter%Atoms(1)%El_EMFP%E
   matter%El_EMFP_tot%L = 0.0d0
   

   ATOMS:do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      SHELLS:do j = 1, Nshl ! for all shells of this atom
         !redo = .false. ! may be there is no need to recalculate MFPs
         redo = numpar%redo_MFP  ! user defined, if there is a need to recalculate MFPs

         ! Check if CDF coefficients are set:
         if (.not. allocated(matter%Atoms(i)%CDF)) then
            allocate(matter%Atoms(i)%CDF(Nshl))
         endif
         ! And if not, use single-pole approximation:
         if (.not. allocated(matter%Atoms(i)%CDF(j)%A) ) then
            call set_single_pole_CDF(Scell, NSC, matter, numpar, i, j)  ! below
         endif

         if ((i .NE. 1) .or. (j .NE. 1)) then
            if (.not.allocated(matter%Atoms(i)%El_MFP(j)%E)) allocate(matter%Atoms(i)%El_MFP(j)%E(N_grid))
            matter%Atoms(i)%El_MFP(j)%E = matter%Atoms(1)%El_MFP(1)%E
            if (.not.allocated(matter%Atoms(i)%El_MFP(j)%L)) allocate(matter%Atoms(i)%El_MFP(j)%L(N_grid))
         endif
         write(chtemp,'(i6)') INT(matter%Atoms(i)%Ip(j))
         if (j > 1) then    ! check if it's a degenerate level:
            if (INT(matter%Atoms(i)%Ip(j-1)) == INT(matter%Atoms(i)%Ip(j))) then
               ! artificially shift it a little bit to make it not exactly degenerate:
               write(chtemp,'(i6)') INT(matter%Atoms(i)%Ip(j) - 0.5d0)
            endif
         endif
         
         select case (matter%Atoms(i)%TOCS(j)) ! which inelastic cross section to use (BEB vs CDF):
         case (1) ! CDF
            !print*, 'get_MFPs-1:', matter%Atoms(i)%TOCS(j), trim(adjustl(numpar%At_base))
            select case (trim(adjustl(numpar%At_base)))
            case ('CDF') ! cdf from file
               write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
                  trim(adjustl(matter%Atoms(i)%Name))//'_CDF_Electron_IMFP_Ip='//trim(adjustl(chtemp))//'eV.txt'
            case ('CDF_sp', 'CDF:EPICS') ! single-pole cdf
               write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
                  trim(adjustl(matter%Atoms(i)%Name))//'_CDFsp_Electron_IMFP_Ip='//trim(adjustl(chtemp))//'eV.txt'
            end select
         case default ! BEB
            !print*, 'get_MFPs-0:', matter%Atoms(i)%TOCS(j), trim(adjustl(numpar%At_base))
            write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
            trim(adjustl(matter%Atoms(i)%Name))//'_BEB_Electron_IMFP_Ip='//trim(adjustl(chtemp))//'eV.txt'
         end select
         FN = 112+i*Nshl+j

         inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if this file already exists
         if (file_exists) then ! IMFPs are already there
            open(UNIT=FN, FILE = trim(adjustl(File_name)))
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then ! there is no such file, create a new one
               open(UNIT=FN, FILE = trim(adjustl(File_name)))
               inquire(file=trim(adjustl(File_name)),opened=file_opened)
               redo = .true. ! no data, need to recalculate the MFPs
               print*, 'Could not open file with MFP => recalculating MFP'
            else
               call Count_lines_in_file(FN, N) ! get how many lines in the file
               !print*, 'N:', N, N_grid
               if (N .NE. N_grid) then ! replace file
                  redo = .true. ! not enough data, need to recalculate the MFPs
                  print*, 'Energy grid mismatch in MFP file => recalculating MFP'
               endif
            endif
         else
            open(UNIT=FN, FILE = trim(adjustl(File_name)))
            redo = .true. ! no data, need to recalculate the MFPs
            if (numpar%verbose) print*, 'No file with MFP present => calculating MFP'
         endif

9900     if (redo) then ! recalculate the MFPs:
            write(0,'(a,f8.1)') ' Calculating electron IMFP for '//trim(adjustl(matter%Atoms(i)%Name))//', Ip=', matter%Atoms(i)%Ip(j)
            !!$omp PARALLEL private(k,Ele,L,dEdx)
            !!$omp do
            do k = 1, N_grid
               Ele = matter%Atoms(i)%El_MFP(j)%E(k) ! [eV] energy
               call TotIMFP(Ele, matter, TeeV, i, j, L, dEdx=dEdx)
               matter%Atoms(i)%El_MFP(j)%L(k) = L ! [A] MFP
               write(FN,'(f25.16,es25.16,es25.16)') Ele, L, dEdx
               call print_progress('Progress:',k,N_grid)    ! module "Little_subroutines"
               !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
               !write(*,*) 'Thread #', my_id, k
               !print*, 'MFP:', k, Ele, L
            enddo
            !!$omp end do
            !!$omp end parallel
            print*, 'Electron IMFPs are saved into file:', trim(adjustl(File_name))
         else ! Just read MFPs from the file:
            !print*, 'Reding electron IMFPs from file: ', trim(adjustl(File_name))
            count_lines = 0
            do k = 1, N_grid
               read(FN,*,IOSTAT=Reason) Ele, L, dEdx ! read the line
               call read_file(Reason, count_lines, read_well)
               if (.not. read_well) then
                  redo = .true. ! no data, need to recalculate the MFPs
                  goto 9900 ! if couldn't read the file, just recalculate it then
               endif
               matter%Atoms(i)%El_MFP(j)%E(k) = Ele ! [eV] energy
               matter%Atoms(i)%El_MFP(j)%L(k) = L   ! [A] MFP
            enddo
         endif
         inquire(file=trim(adjustl(File_name)),opened=file_opened)
         if (file_opened) close(FN)

         ! And check the sum rules, if needed:
         if (redo .or. numpar%verbose) then
            select case (matter%Atoms(i)%TOCS(j)) ! which inelastic cross section to use (BEB vs CDF):
            case (1) ! CDF
               Omega = w_plasma(1d6*matter%At_dens) ! function below, atomic-density plasma frequency [1/s^2]
               if ( (i == 1) .and. (j == Nshl) ) then ! the valence band
                  call sumrules(matter%Atoms(i)%CDF(j)%A, matter%Atoms(i)%CDF(j)%E0, matter%Atoms(i)%CDF(j)%G, &
                                                ksum, fsum, Scell(NSC)%E_gap, Omega) ! below
               else ! core shell
                  call sumrules(matter%Atoms(i)%CDF(j)%A, matter%Atoms(i)%CDF(j)%E0, matter%Atoms(i)%CDF(j)%G, &
                                                ksum, fsum, matter%Atoms(i)%Ip(j), Omega) ! below
               endif
               write(*,'(a,f7.2,e12.2)') 'Sum rules for '//trim(adjustl(matter%Atoms(i)%Name))//', Ip='// &
                                                trim(adjustl(chtemp))//'eV are (k,f):', ksum, fsum !, sqrt(Omega)
            end select
         endif

      enddo SHELLS
      
      ! Get also elastic MFPs for each atomic species:
      call Elastic_MFP(i, numpar, matter) ! see below
      
   enddo ATOMS

   ! Printout plasmon frequency:
   Omega = w_plasma(1d6*matter%At_dens* Scell(NSC)%Ne/dble(Scell(NSC)%Na) ) ! function below, plasma frequency^2 [1/s^2]
   matter%hw_plasma = g_h/g_e*sqrt(Omega) ! [eV]
   write(*,'(a,es12.2,a,f12.3,a)') 'Estimated plasmon frequency: ', sqrt(Omega), ' [1/s],', g_h/g_e*sqrt(Omega), ' [eV]'


   ! Total inelastic mean free path:
   write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Name))//'_Total_Electron_IMFP.txt'
   FN = FN + 1
   open(UNIT=FN, FILE = trim(adjustl(File_name)))
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl ! for all shells of this atom
         matter%El_MFP_tot%L(:) = matter%El_MFP_tot%L(:) + 1.0d0/matter%Atoms(i)%El_MFP(j)%L(:) ! [1/A] inverse MFP
      enddo
   enddo

   do i = 1, N_grid
      if (matter%El_MFP_tot%L(i) .LE. 0.0d0) then
         matter%El_MFP_tot%L(i) = 1d30
      else
         matter%El_MFP_tot%L(i) = 1.0d0/matter%El_MFP_tot%L(i) ! [A]
      endif
      write(FN,'(f25.16,es25.16)') matter%El_MFP_tot%E(i), matter%El_MFP_tot%L(i)
   enddo


   !TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
   ! Get the mean free paths vs Te:
   Nshl = size(matter%Atoms(1)%Ip)
   select case (matter%Atoms(1)%TOCS(Nshl)) ! Valence band and CDF only
   case (1) ! CDF
      !!$omp PARALLEL private(i, Te_temp)
      !!$omp do schedule(dynamic)
      do i = 1, N_Te_points   ! for all electronic temperature points
         Te_temp = dble((i-1)*1000)*g_kb_EV ! electronic temperature [eV]
         call IMFP_vs_Te_files(matter, laser, numpar, Te_temp, i) ! below
      enddo ! i = 1, N_Te_points
      !!$omp end do
      !!$omp end parallel
   endselect

   !iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
   ! Make the inverse MFPs, since that's how we use them in the MC code:
   matter%El_MFP_tot%L = 1.0d0/matter%El_MFP_tot%L ! [1/A]
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl ! for all shells of this atom
         matter%Atoms(i)%El_MFP(j)%L(:) = 1.0d0/matter%Atoms(i)%El_MFP(j)%L(:) ! [1/A] inverse MFP
      enddo
   enddo
   !iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) close(FN)
   
   ! And also elastic total MFP:
   if (.not.numpar%do_elastic_MC) then ! if elastic scattering is excluded
      ! Set nearly zero inverse meen free paths, meaning infinite mean free path:
      do i = 1, size(matter%Atoms) ! for all atoms
         matter%Atoms(i)%El_EMFP%L(:) = 1.0d-30	! [1/A] inverse MFP
      enddo
   endif
   call Get_total_EMFP(numpar, matter) ! see below

9999 continue
end subroutine get_MFPs



subroutine set_single_pole_CDF(Scell, NSC, matter, numpar, i, j)  ! only for VB/CB
   type(Super_cell), dimension(:), intent(in) :: Scell ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! supercell index
   type(Solid), intent(inout) :: matter ! parameters of the material
   type(Numerics_param), intent(in) :: numpar  ! all numerical parameters
   integer, intent(in) :: i, j   ! index of atom and shell
   !----------------
   real(8) :: Omega, NVB, ksum, fsum, Nat, contrib
   integer :: Nshl

   Nshl = size(matter%Atoms(1)%Ip)  ! index of the valence band

   if (.not. allocated(matter%Atoms(i)%N_CDF)) allocate(matter%Atoms(i)%N_CDF(matter%Atoms(i)%sh)) ! allocate number of electrons

   if (.not.allocated(matter%Atoms(i)%CDF(j)%A)) then
      matter%Atoms(i)%N_CDF(j) = 1  ! set single CDF, coefficients to be determined
      allocate(matter%Atoms(i)%CDF(j)%A(matter%Atoms(i)%N_CDF(j)))
      allocate(matter%Atoms(i)%CDF(j)%E0(matter%Atoms(i)%N_CDF(j)))
      allocate(matter%Atoms(i)%CDF(j)%G(matter%Atoms(i)%N_CDF(j)))
   endif

   ! Check if the flag for CS needs to be changed:
   select case (trim(adjustl(numpar%At_base)))
   case('CDF', 'CDF_sp')
      matter%Atoms(i)%TOCS(j) = 1      ! mark it as CDF cross-section
      matter%Atoms(i)%TOCSph(j) = 1    ! mark it as CDF cross-section
   case('CDF:EADL', 'CDF:EPICS')
      matter%Atoms(i)%TOCS(j) = 1      ! mark it as CDF cross-section
      matter%Atoms(i)%TOCSph(j) = 0    ! mark it as EPICS cross-section
   case default   ! EADL, BEB
      ! do not change the flag
   end select

   if ( (i == 1) .and. (j == Nshl) ) then ! the valence band
      NVB = dble(Scell(NSC)%Ne) / dble(Scell(NSC)%Na) ! valence electrons per atom

      ! Set them according to the single-pole approximation:
      Omega = w_plasma(1d6*matter%At_dens*NVB) ! function below, plasma frequency [1/s]

      matter%Atoms(i)%CDF(j)%E0(1) = sqrt((g_h/g_e)*(g_h/g_e) * Omega)  ! [eV]
      ! Gamma set equal to E0 without effective mass (empirical approximation):
      matter%Atoms(i)%CDF(j)%G(1) = matter%Atoms(i)%CDF(j)%E0(1)
      ! A is set vie normalization (sum rule):
      matter%Atoms(i)%CDF(j)%A(1) = 1.0d0   ! just to get sum rule to renormalize below
      ! Get sum rule:
      Omega = w_plasma(1d6*matter%At_dens)   ! below
      call sumrules(matter%Atoms(i)%CDF(j)%A, matter%Atoms(i)%CDF(j)%E0, matter%Atoms(i)%CDF(j)%G, ksum, fsum, Scell(NSC)%E_gap, Omega) ! below

      matter%Atoms(i)%CDF(j)%A(1) = NVB/ksum

   else ! core shell
      NVB = matter%Atoms(i)%Ne_shell(j) ! electrons in this shell
      Nat = SUM(matter%Atoms(:)%percentage)  ! number of different kinds of atoms in the compound
      contrib = matter%Atoms(i)%percentage/Nat  ! contribution of the atoms into the compound
      !print*, i, matter%Atoms(i)%Name, contrib

      ! Set them according to the single-pole approximation:
      Omega = w_plasma(1d6*matter%At_dens*NVB*contrib) ! function below, plasma frequency [1/s]

      matter%Atoms(i)%CDF(j)%E0(1) = matter%Atoms(i)%Ip(j)+10.0d0  ! [eV] -- approximation
      ! Gamma set equal to E0 without effective mass (empirical approximation):
      matter%Atoms(i)%CDF(j)%G(1) = matter%Atoms(i)%CDF(j)%E0(1)
      ! A is set vie normalization (sum rule):
      matter%Atoms(i)%CDF(j)%A(1) = 1.0d0   ! just to get sum rule to renormalize below
      ! Get sum rule:
      Omega = w_plasma(1d6*matter%At_dens)   ! below
      call sumrules(matter%Atoms(i)%CDF(j)%A, matter%Atoms(i)%CDF(j)%E0, matter%Atoms(i)%CDF(j)%G, ksum, fsum, matter%Atoms(i)%Ip(j), Omega) ! below

      matter%Atoms(i)%CDF(j)%A(1) = NVB/ksum
   endif !( (i == 1) .and. (j == Nshl) )
end subroutine set_single_pole_CDF



subroutine IMFP_vs_Te_files(matter, laser, numpar, Te, N_Te)
   type(Solid), intent(inout) :: matter ! parameters of the material
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar  ! all numerical parameters
   real(8), intent(in) :: Te  ! [eV] electronic tempreature
   integer, intent(in) :: N_Te   ! index
   !-------------------
   logical :: redo, file_exists, file_opened, read_well
   character(10) :: chtemp, ch_Te
   character(200) :: File_name
   integer :: Nshl, FN, N, N_grid, k, Reason, count_lines
   real(8) :: Ele, L, dEdx


   redo = .false. ! may be there is no need to recalculate MFPs
   Nshl = size(matter%Atoms(1)%Ip)
   N_grid = size(matter%Atoms(1)%El_MFP_vs_T(1)%L)
   matter%Atoms(1)%El_MFP_vs_T(N_Te)%E(:) = matter%Atoms(1)%El_MFP(Nshl)%E(:) ! copy energy grid [eV]

   write(chtemp,'(i6)') INT(matter%Atoms(1)%Ip(Nshl))
   write(ch_Te,'(i6)') INT(Te*g_kb)

   select case (trim(adjustl(numpar%At_base)))
   case ('CDF') ! cdf from file
      write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
      trim(adjustl(matter%Atoms(1)%Name))//'_CDF_Electron_IMFP_Ip='//trim(adjustl(chtemp))//'eV_'//trim(adjustl(ch_Te))//'K.txt'
   case ('CDF_sp', 'CDF:EPICS') ! single-pole cdf
      write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, &
      trim(adjustl(matter%Atoms(1)%Name))//'_CDFsp_Electron_IMFP_Ip='//trim(adjustl(chtemp))//'eV_'//trim(adjustl(ch_Te))//'K.txt'
   end select
   FN = 312+N_Te   ! file number

   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if this file already exists
   if (file_exists) then ! IMFPs are already there
      open(UNIT=FN, FILE = trim(adjustl(File_name)))
      inquire(file=trim(adjustl(File_name)),opened=file_opened)
      if (.not.file_opened) then ! there is no such file, create a new one
         open(UNIT=FN, FILE = trim(adjustl(File_name)))
         inquire(file=trim(adjustl(File_name)),opened=file_opened)
         redo = .true. ! no data, need to recalculate the MFPs
         print*, 'Could not open file with MFP => recalculating MFP'
      else
         call Count_lines_in_file(FN, N) ! get how many lines in the file
         !print*, 'N:', N, N_grid
         if (N .NE. N_grid) then ! replace file
            redo = .true. ! not enough data, need to recalculate the MFPs
            print*, 'Energy grid mismatch in MFP file => recalculating MFP'
         endif
      endif
   else
      open(UNIT=FN, FILE = trim(adjustl(File_name)))
      redo = .true. ! no data, need to recalculate the MFPs
      if (numpar%verbose) print*, 'No file with MFP present => calculating MFP'
   endif

9901 if (redo) then ! recalculate the MFPs:
      write(0,'(a,f8.1)') ' Calculating electron IMFP for valence band at Te=', Te*g_kb
      do k = 1, N_grid
         Ele = matter%Atoms(1)%El_MFP(Nshl)%E(k) ! [eV] energy grid
         call TotIMFP(Ele, matter, Te, 1, Nshl, L, dEdx=dEdx)   ! below
         matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(k) = L ! [A] MFP
         write(FN,'(f25.16,es25.16,es25.16)') Ele, L, dEdx
         call print_progress('Progress:',k,N_grid)    ! module "Little_subroutines"
         !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
         !write(*,*) 'Thread #', my_id, k
         !print*, 'MFP:', k, Ele, L
      enddo
      print*, 'Electron IMFPs are saved into file:', trim(adjustl(File_name))

   else ! Just read MFPs from the file:
      count_lines = 0 ! to start with
      do k = 1, N_grid
         read(FN,*,IOSTAT=Reason) Ele, L, dEdx ! read the line
         call read_file(Reason, count_lines, read_well)
         if (.not. read_well) then
            redo = .true. ! no data, need to recalculate the MFPs
            goto 9901 ! if couldn't read the file, just recalculate it then
         endif
         matter%Atoms(1)%El_MFP_vs_T(N_Te)%E(k) = Ele ! [eV] energy
         matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(k) = L   ! [A] MFP
      enddo
   endif
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) close(FN)

   matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) = 1.0d0/matter%Atoms(1)%El_MFP_vs_T(N_Te)%L(:) ! [1/A] inverse MFP
end subroutine IMFP_vs_Te_files



subroutine Electron_energy_transfer_inelastic(matter, TeeV, Ele, Nat, Nshl, mfps, dE_out)
    type(solid), intent(in) :: matter	! materil parameters
    real(8), intent(in) :: Ele  ! electron energy [eV]
    real(8), intent(in) :: TeeV  ! electronic tempreature [eV]
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    type(MFP), intent(in) :: mfps	! electron mean free paths for each shell [1/A]
    real(8), intent(out) :: dE_out   ! the transferred energy [eV]
    real(8) :: L_tot    ! [A] total mean free path
    integer i, j, n, coun
    real(8) :: Emin, Emax, E_cur, E, dE, dL, Ltot1, Ltot0, ddEdx, a, b, temp1, temp2, RN, L_need, L_cur, Sigma_cur, Tfact
    real(8) :: E_low, E_high

    if (Ele .LT. matter%Atoms(Nat)%Ip(Nshl)) then
       print*, 'Attention! In subroutine Electron_energy_transfer_inelastic:'
       print*, 'Transferred energy is lower than the ionization potential!'
       write(*,'(f25.16, f25.16, i2, i2)') Ele, matter%Atoms(Nat)%Ip(Nshl), Nat, Nshl
       write(*,'(f25.16,$)') matter%Atoms(Nat)%Ip(:)
       write(*,'(a)') ''
       pause 'Electron_energy_transfer_inelastic'
    endif

    call Mean_free_path(Ele, mfps, L_tot) ! [A] total mean free path
    call random_number(RN)
    L_need = L_tot/RN   ! [A] we need to reach
    Emin = matter%Atoms(Nat)%Ip(Nshl) ! [eV] ionization potential of the shell is minimum possible transferred energy
    if (Emin .LE. 1.0d-3) Emin = 1.0d-3 ! for metals there is no band gap

    L_cur = 1.0d10 ! to start
    select case (matter%Atoms(Nat)%TOCS(Nshl)) ! which inelastic cross section to use (BEB vs CDF):
    case (1) ! CDF cross section
       Emax = (Ele + Emin)/2.0d0 ! [eV] maximum energy transfer, accounting for equality of electrons
       !n = 10*(MAX(INT(Emin),10))    ! number of integration steps ! OLD
       !dE = (Emax - Emin)/(real(n)) ! differential of transferred momentum [kg*m/s]
       n = 100
       ! Define the interval of integration where the peak are (requires fined grid):
       E_low = max( minval(matter%Atoms(Nat)%CDF(Nshl)%E0 - 5.0d0*matter%Atoms(Nat)%CDF(Nshl)%G) , Emin )
       E_high = min( maxval(matter%Atoms(Nat)%CDF(Nshl)%E0 + 5.0d0*matter%Atoms(Nat)%CDF(Nshl)%G) , Emax )

       i = 1       ! to start integration
       E = Emin    ! to start integration
       Ltot1 = 0.0d0
       Ltot0 = 0.0d0
       call Diff_cross_section(Ele, E, matter, Nat, Nshl, Ltot0)
       do while (L_cur .GT. L_need)
          !dE = (1.0d0/(E+1.0d0) + E)/real(n)
          !dE = define_dE(-1, n, E, E0_min=E_low, E0_max=E_high)  ! below
          dE = define_dE(0, n, E, E0_min=E_low, E0_max=E_high)  ! below

          ! Temperature factor:
          Tfact = temperature_factor(dE, TeeV)    ! below

          ! If it's Simpson integration:
          a =  E + dE/2.0d0
          call Diff_cross_section(Ele, a, matter, Nat, Nshl, dL)
          temp1 = dL
          b = E + dE
          call Diff_cross_section(Ele, b, matter, Nat, Nshl, dL)
          temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)

          ! Include temperature factor:
          temp2 = temp2*Tfact

          Ltot1 = Ltot1 + temp2
          Ltot0 = dL
          L_cur = 1.0d0/Ltot1
          E = E + dE  ! [eV]
          if (E .GE. Emax) exit
       enddo
       if (E .GE. Emax) E = Emax
       dE_out = E  ! energy transfer [eV]
    case default ! BEB cross section
       Emax = (Ele-matter%Atoms(Nat)%Ip(Nshl))/2.0d0
       E_cur = 0.0d0
       coun = 0
       do while (ABS(L_cur-L_need)/L_need .GT. 0.0001d0) ! search the transferred energy by bisection:
         coun = coun + 1	! just count the loops
         Sigma_cur = dSigma_int_BEB(Ele, E_cur, matter%Atoms(Nat)%Ip(Nshl), matter%Atoms(Nat)%Ek(Nshl), matter%Atoms(Nat)%Ne_shell(Nshl))
         temp1 = matter%At_Dens*1d-24*(matter%Atoms(Nat)%percentage)/SUM(matter%Atoms(:)%percentage)
         if (Sigma_cur .LE. 0.0d0) then
            L_cur = 1d30 ! [A]
         else
            L_cur = 1.0d0/(temp1*Sigma_cur) ! IMFP [A]
         endif
         if (L_cur .GT. L_need) then
            Emin = E_cur
         else
            Emax = E_cur
         endif
         E_cur = (Emax+Emin)/2.0d0
         if (coun .GE. 1d3) exit
       enddo
       dE_out = matter%Atoms(Nat)%Ip(Nshl) + E_cur ! energy transfer [eV]
    end select
end subroutine Electron_energy_transfer_inelastic


!subroutine TotIMFP(Ele, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, Mat_DOS_inv, NumPar)
subroutine TotIMFP(Ele, matter, TeeV, Nat, Nshl, Sigma, dEdx)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    !type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    type(Solid), intent(in) :: matter ! parameters of the material
    real(8), intent(in) :: TeeV  ! [eV] electronic tempreature
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma       ! calculated inverse mean free path (cross-section) [1/A],
    real(8), intent(out), optional :: dEdx   ! the energy losses [eV/A]
    !=====================================
    integer :: i, j, n, num
    real(8) :: Emin, Emax, E, dE, dL, Ltot1, Ltot0, ddEdx, a, b, temp1, temp2, Tfact
    real(8) :: E_low, E_high

    Emin = matter%Atoms(Nat)%Ip(Nshl) ! [eV] ionization potential of the shell is minimum possible transferred energy

    if (Emin .LE. 1.0d-2) Emin = 1.0d-2 ! for metals there is no band gap
    Emax = (Ele + Emin)/2.0d0 ! [eV] maximum energy transfer, accounting for equality of electrons

    select case (matter%Atoms(Nat)%TOCS(Nshl)) ! which inelastic cross section to use (BEB vs CDF):
    case (1) ! CDF cross section
       if (Ele < Emin) then ! no scattering possible
         Ltot1 = 0.0d0
         ddEdx = 0.0d0
         Sigma = 1d30 ! [A]
       else ! (Ee < Emin) ! scattering possible
         !n = 10*(MAX(INT(Emin),10))    ! number of integration steps   ! OLD
         !dE = (Emax - Emin)/(real(n)) ! differential of transferred energy [eV]
         n = 100
         ! Define the interval of integration where the peak are (requires fined grid):
         E_low = max( minval(matter%Atoms(Nat)%CDF(Nshl)%E0 - 5.0d0*matter%Atoms(Nat)%CDF(Nshl)%G) , Emin )
         E_high = min( maxval(matter%Atoms(Nat)%CDF(Nshl)%E0 + 5.0d0*matter%Atoms(Nat)%CDF(Nshl)%G) , Emax )

         i = 1       ! to start integration
         E = Emin    ! to start integration
         Ltot1 = 0.0d0
         Ltot0 = 0.0d0
         !call Diff_cross_section(Ele, E, Target_atoms, Nat, Nshl, Ltot0, Mass, Matter, Mat_DOS_temp, NumPar)
         call Diff_cross_section(Ele, E, matter, Nat, Nshl, Ltot0)  ! below
         ddEdx = 0.0d0
         do while (E .LE. Emax) ! integration
            !dE = (1.0d0/(E+1.0d0) + E)/real(n)
            !dE = define_dE(-1, n, E, E0_min=E_low, E0_max=E_high)  ! below
            dE = define_dE(0, n, E, E0_min=E_low, E0_max=E_high)  ! below

            ! If it's Simpson integration:
            a =  E + dE/2.0d0
            ! Temperature factor (for the mean point):
            Tfact = temperature_factor(a, TeeV)    ! below
            !call Diff_cross_section(Ele, a, Target_atoms, Nat, Nshl, dL, Mass, Matter, Mat_DOS_temp, NumPar)
            call Diff_cross_section(Ele, a, matter, Nat, Nshl, dL)  ! below
            temp1 = dL
            b = E + dE
            !call Diff_cross_section(Ele, b, Target_atoms, Nat, Nshl, dL, Mass, Matter, Mat_DOS_temp, NumPar)
            call Diff_cross_section(Ele, b, matter, Nat, Nshl, dL)  ! below
            temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
            temp2 = temp2*Tfact    ! include temprature factor

            Ltot1 = Ltot1 + temp2
            if (present(dEdx)) ddEdx = ddEdx + E*temp2
            Ltot0 = dL
            E = E + dE  ! [eV]
         enddo
         Sigma = 1.0d0/Ltot1 !*dE ! [A]
       endif ! (Ee < Emin)
       if (present(dEdx)) dEdx = ddEdx !*dE ! energy losses [eV/A]
    case default ! BEB cross section
       Sigma = Sigma_BEB(Ele, matter%Atoms(Nat)%Ip(Nshl), matter%Atoms(Nat)%Ek(Nshl), matter%Atoms(Nat)%Ne_shell(Nshl)) ! [A^2] cross section, function below

       temp1 = matter%At_Dens*1d-24*(matter%Atoms(Nat)%percentage)/SUM(matter%Atoms(:)%percentage)
       if (Sigma .LE. 0.0d0) then
          Sigma = 1d30 ! [A]
       else
          Sigma = 1.0d0/(temp1*Sigma) ! IMFP [A]
       endif
       if (present(dEdx)) ddEdx = dSigma_w_int_BEB(Ele, (Ele-1.0d0)/2.0d0, matter%Atoms(Nat)%Ip(Nshl), matter%Atoms(Nat)%Ek(Nshl), matter%Atoms(Nat)%Ne_shell(Nshl)) ! cross section integrated with energy [A^2*eV]
       if (present(dEdx)) dEdx = temp1*ddEdx ! energy losses [eV/A]

       !print*, Sigma_BEB(Ele,matter%Atoms(Nat)%Ip(Nshl),matter%Atoms(Nat)%Ek(Nshl),matter%Atoms(Nat)%Ne_shell(Nshl)), Sigma, temp1, matter%At_Dens, matter%Atoms(Nat)%percentage, SUM(matter%Atoms(:)%percentage)
       !pause 'BEB Sigma'

    end select
end subroutine TotIMFP



pure function define_dE(CS_method, n, E, E0_min, E0_max, dE_min) result(dE)
   real(8) dE
   integer, intent(in) :: CS_method, n ! integration grid index; number of grid points
   real(8), intent(in) :: E
   real(8), intent(in), optional :: E0_min, E0_max, dE_min ! both must be present for NEW grid; dE_min for minimal step optional
   !-----------------------
   real(8) :: dE_min_use

   if (present(dE_min)) then
      dE_min_use = dE_min
   else ! use default value (typically the case for inelastic scsattering)
      dE_min_use = 0.001d0   ! [eV] step smaller than this is not allowed
   endif

   select case (CS_method)
   case (-1)   ! Old default
      dE = (1.0d0/(E+1.0d0) + E)/dble(n)

   case default ! new integration grid
      if ( (present(E0_min)) .and. (present(E0_max)) ) then
         if ( (E > E0_min) .and. (E < E0_max) ) then   ! fine grid in the interval
            dE = (E0_max - E0_min)/dble(n)
         elseif (E > E0_max) then ! high-energy - too far from the peak, no need for fine resolution
            dE = E/dble(n)
         else ! low energy - too far from the peak, no need for fine resolution (usually unused, but anyway...)
            dE = (E-E0_min)/dble(n) ! step not smaller than this
         endif

      else ! use old grid
         dE = (1.0d0/(E+1.0d0) + E)/dble(n)  ! start with old, if there are no parameters for the new
      endif
   end select

   dE = max(dE, dE_min_use)  ! step not smaller than this
end function define_dE


!subroutine Diff_cross_section(Ele, hw, Target_atoms, Nat, Nshl, Diff_IMFP, Mass, Matter, Mat_DOS, NumPar)
subroutine Diff_cross_section(Ele, hw, matter, Nat, Nshl, Diff_IMFP)
    real(8), intent(in), target :: Ele  ! SHI energy [eV]
    REAL(8), INTENT(in), target :: hw   ! transferred energy [eV]
    !type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    type(Solid), intent(in) :: matter ! parameters of the material
    integer, intent(in) :: Nat, Nshl    ! number of atom, and number of shell    
    real(8), intent(out) :: Diff_IMFP   ! differential inverse mean free path 1/lambda(Ele,dE)
    !real(8), intent(in) :: Mass ! mass of incomming particle (electron)
    !===========================
    integer i, n
    real(8), pointer :: Ee, dE
    real(8) dLs, qmin, qmax, hq, ddq, dq, Ime, dLs0, dL, hq0, dq_save
    real(8) a, b, x, temp1, temp2
    Ee => Ele        ! energy [eV]
    dE => hw         ! transferred energy [eV]

    !qmin = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee) - sqrt((Ee - dE))) ! min transferred momentum [kg*m/s]
    !qmax = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee) + sqrt((Ee - dE))) ! max transferred momentum [kg*m/s]
    qmin = sqrt(2.0d0*g_me)/g_h*(sqrt(Ee) - sqrt(abs(Ee - dE))) ! min transferred momentum [kg*m/s]
    qmax = sqrt(2.0d0*g_me)/g_h*(sqrt(Ee) + sqrt(abs(Ee - dE))) ! max transferred momentum [kg*m/s]

    dLs = 0.0d0 ! starting integration, mean free path per energy [A/eV]^(-1)
    hq = qmin    ! transient transferred momentum for integration [kg*m/s]
    n = 100
    dq = (qmax - qmin)/real(n) ! differential of transferred momentum [kg*m/s]
    dLs0 = 0.0d0
    do while (hq .LT. qmax) ! no matter how many points, go till the end
        dq = hq/real(n)
        ! If it's Simpson integration:
        a = hq + dq/2.0d0
        !call Imewq(hw, a, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        call Imewq(matter%Atoms(Nat)%CDF, hw, a, Nshl, ImE)
        temp1 = ImE
        b = hq + dq
        !call Imewq(hw, b, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        call Imewq(matter%Atoms(Nat)%CDF, hw, b, Nshl, ImE)
        dL = ImE
        dLs = dLs + dq/6.0d0*(dLs0 + 4.0d0*temp1 + dL)/hq
        dLs0 = dL
        hq = hq + dq
    enddo
    Diff_IMFP = 1.0d0/(g_Pi*g_a0*Ele)*dLs
    nullify(Ee)
    nullify(dE)
end subroutine Diff_cross_section


pure function temperature_factor(W, T) result (F)
   real(8) :: F
   real(8), intent(in) :: W ! [eV] transferred energy
   real(8), intent(in) :: T  ! [eV] target temperature
   real(8) :: eps
   eps = 1.0d-4
   if ((T < eps) .or. (W < 1.0d-10)) then    ! zero temperature or no energy
      F = 1.0d0
   else ! non-zero tempreature
      F = 1.0d0/(1.0d0 - exp(-W/T))
   endif
end function temperature_factor


!subroutine Imewq(matter, hw, hq, shl, ImE) ! constructs Im(-1/e(w,q)) as a sum of Drude-like functions
subroutine Imewq(CDF, hw, hq, shl, ImE) ! constructs Im(-1/e(w,q)) as a sum of Drude-like functions
    !type(solid), intent(in) :: matter	! materil parameters
    type(Ritchi), dimension(:), target :: CDF	! coefficients of CDF
    REAL(8), INTENT(in) :: hw    ! transferred energy [eV]
    REAL(8), INTENT(in) :: hq    ! transferred momentum [kg*m/s]
    integer, intent(in) :: shl   ! number of shell which is being ionized
    REAL(8), INTENT(out) :: ImE  ! Im(-1/e(w,q))
    real(8), pointer :: A, Gamma, E
    real(8) sumf
    integer i, N
    ImE = 0.0d0
    !N = matter%N_CDF(shl) ! number of CDF functions for this shell
    N = size(CDF(shl)%A) ! number of CDF functions

    do i = 1, N !Nff_esh(Nosh)
        !A = matter%CDF(shl)%A(i)
        !E = matter%CDF(shl)%E0(i)
        !Gamma = matter%CDF(shl)%G(i)
        A => CDF(shl)%A(i)
        E => CDF(shl)%E0(i)
        Gamma => CDF(shl)%G(i)
        sumf = Diel_func(A,E,Gamma, hw, hq) ! Diel_func - function, see below
        ImE = ImE + sumf
    enddo
    nullify(A, Gamma, E)
end subroutine Imewq

function Diel_func(A,E,Gamma,dE,dq) ! fit functions in Ritchi algorithm
    real(8) A, E, Gamma, dE, dq  ! parameters and variable
    real(8) Diel_func ! function itself
    real(8) E0, dE2, E02, dE2E02 ! temporary parameters
    E0 = E + g_h*g_h*dq*dq/(2.0d0*g_me)
    dE2 = dE*dE
    E02 = E0*E0
    dE2E02 = dE2 - E02
    !Diel_func = A*Gamma*dE/((dE*dE - E0*E0)*(dE*dE - E0*E0) + Gamma*Gamma*dE*dE)
    Diel_func = A*Gamma*dE/(dE2E02*dE2E02 + Gamma*Gamma*dE2)
end function Diel_func


function w_plasma(At_dens, Mass)
   real(8) w_plasma ! Squared plasma frequency [1/s]^2
   real(8), intent(in) :: At_dens   ! atomic density [1/m^3]
   real(8), intent(in), optional :: Mass ! effective mass [kg]
   !--------------------

   !w_plasma = (4.0d0*g_Pi*At_dens*g_e*g_e/(4.0d0*g_Pi*g_e0*g_me))

    if (present(Mass)) then  ! use user-provided mass:
      w_plasma = At_dens*g_e*g_e/(g_e0*Mass)
   else ! assume free electron
      w_plasma = At_dens*g_e*g_e/(g_e0*g_me)
   endif
end function w_plasma


subroutine sumrules(Afit, Efit, Gfit, ksum, fsum, x_min, Omega)
    real(8), dimension(:),  INTENT(in), target :: Afit ! A coefficient
    real(8), dimension(:),  INTENT(in), target :: Efit ! E0 coefficient
    real(8), dimension(:),  INTENT(in), target :: Gfit ! Gamma coefficient
    real(8), intent (out) :: ksum, fsum
    real(8), intent(in) ::  x_min, Omega ! ionization potneital [eV] and plasma frequency [1/s]
    real(8), pointer ::  A1, E1, Gamma1
    real(8) ne, f
    integer j
    ne = 0.0d0
    f = 0.0d0
    do j = 1,size(Afit)
        A1 => Afit(j)
        E1 => Efit(j)
        Gamma1 => Gfit(j)
        ne = ne + Int_Ritchi_x(A1,E1,Gamma1,1d10) - Int_Ritchi_x(A1,E1,Gamma1,x_min)
        f = f + Int_Ritchi_p_x(A1,E1,Gamma1,1d10) - Int_Ritchi_p_x(A1,E1,Gamma1,x_min)
    enddo
    ksum = 2.0d0*g_e*g_e/(g_Pi*Omega*g_h*g_h)*ne
    fsum = f*2.0d0/g_Pi
    nullify(A1, E1, Gamma1)
end subroutine sumrules


function Int_Ritchi(A,E,Gamma,x) ! integral of the Ritchi function
    real(8) A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi ! function itself
    real(8) S
    complex(8) Sc
    Sc = (sqrt((2.0d0*E)*(2.0d0*E) - Gamma*Gamma))
    S = REAL(Sc)
    Int_Ritchi = A/S*ATAN( (2.0d0*(x*x-E*E) + Gamma*Gamma)/(Gamma*S) )
end function Int_Ritchi


function Int_Ritchi_x(A,E,Gamma,x) ! integral of the Ritchi*x (k-sum rule)
    real(8) A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi_x ! function itself
    real(8) S, sq2, G, s_plus, s_minus,B
    complex(8) Sc, Gc, s_plus_c, s_minus_c, Bc, Ic, arg,arg2, oneI
    sq2 = sqrt(2.0d0)
    oneI = cmplx(0.0d0,1.0d0)
    if ((-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma) .GE. 0.0d0) then
        Sc = cmplx(sqrt(-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma),0.0d0)
    else
        Sc = cmplx(0.0d0, sqrt(ABS(-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma)))
    endif
    Gc = cmplx(Gamma*Gamma - 2.0d0*E*E,0.0d0)
    s_plus_c = sq2*sqrt(Gc + Gamma*Sc)
    s_minus_c = sq2*sqrt(Gc - Gamma*Sc)
    Bc=Gc/(Gamma*Sc)
    !Int_Ritchi_x = real(A*Gamma*( ATAN(2.0e0*x/s_minus)/s_minus*(1.0e0-B) + ATAN(2.0e0*x/s_plus)/s_plus*(1.0e0+B) ))
    arg = 2.0d0*x/(s_minus_c)
    arg2 = 2.0d0*x/(s_plus_c)
    !Ic = (A*Gamma*( ATAN2(aimag(arg),real(arg))/s_minus_c*(1.0e0-Bc) + ATAN2(aimag(arg2),real(arg2))/s_plus_c*(1.0e0+Bc) ))
    Ic = (A*Gamma*(0.5d0*oneI*(log(1.0d0-oneI*arg)-log(1.0d0+oneI*arg))/s_minus_c*(1.0d0-Bc) + (0.5d0*oneI*(log(1.0d0-oneI*arg2)-log(1.0d0+oneI*arg2)))/s_plus_c*(1.0d0+Bc) ))
    Int_Ritchi_x = real(Ic)
end function Int_Ritchi_x


function Int_Ritchi_p_x(A,E,Gamma,x) ! integral of the Ritchi/x (ff-sum rule)
    real(8) A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi_p_x ! function itself
    real(8) S, sq2, G, s_plus, s_minus
    complex(8) Sc, Gc, s_plus_c, s_minus_c, oneI, In_c, arg, arg2
    sq2 = sqrt(2.0d0)
    oneI = cmplx(0.0d0,1.0d0)
    if ((-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma) .GE. 0.0d0) then
        Sc = cmplx(sqrt(-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma),0.0d0)
    else
        Sc = cmplx(0.0d0, sqrt(ABS(-(2.0d0*E)*(2.0d0*E) + Gamma*Gamma)))
    endif
    Gc = cmplx(Gamma*Gamma - 2.0d0*E*E,0.0d0)
    s_plus_c = sqrt(Gc + Gamma*Sc)
    s_minus_c = sqrt(Gc - Gamma*Sc)
    !Int_Ritchi_p_x = sq2*A/S*( ATAN(sq2*x/s_minus)/s_minus - ATAN(sq2*x/s_plus)/s_plus )
    arg = sq2*x/s_minus_c
    arg2 = sq2*x/s_plus_c
    In_c = sq2*A/Sc*( (0.5d0*oneI*(log(1.0d0-oneI*arg)-log(1.0d0+oneI*arg)))/s_minus_c - (0.5d0*oneI*(log(1.0d0-oneI*arg2)-log(1.0d0+oneI*arg2)))/s_plus_c )
    Int_Ritchi_p_x = real(In_c)
end function Int_Ritchi_p_x


subroutine get_grid_4CS(N, Emax_given, grid_array, matter)
   integer, intent(out) :: N ! number of grid points
   real(8), intent(in) :: Emax_given ! [eV] up to which energy we need electron cross-sections
   real(8), intent(out), dimension(:), allocatable :: grid_array ! array of these grid points
   type(Solid), intent(in) :: matter ! parameters of the material: ionization potentials to set special points
   !-----------------------------
   real(8), dimension(:), allocatable :: special_point
   real(8) :: Emin, Emax, E_sp_eps
   integer :: i, NP

   ! Initial definitions:
   Emin = 0.1d0   ! [eV] we start with this minimum
   Emax = 50.0d3  ! [eV] defaul value, may be changed below
   if (Emax_given .GT. Emax) Emax = Emax_given ! [eV] maximum energy for cross section (or mean free path)
   E_sp_eps = 1.0d-3 ! how close a grid point should be around a special point

   ! Get the special points associated with the ionization potentials of all shells:
   call define_special_points(matter, special_point)  ! below

   ! Count how many points, to allocate the grid:
   call go_thru_grid(Emin, Emax, E_sp_eps, special_point, Ngrid=N)   ! below

   ! Save the grid:
   allocate(grid_array(N), source = 0.0d0)
   call go_thru_grid(Emin, Emax, E_sp_eps, special_point, array=grid_array)   ! below

!    do i = 1, size(grid_array)
!        print*, i, grid_array(i)
!    enddo
!    pause 'get_grid_4CS'

   if (allocated(special_point)) deallocate(special_point)
end subroutine get_grid_4CS


subroutine go_thru_grid(Emin, Emax, E_sp_eps, special_point, Ngrid, array)
   real(8), intent(in) :: Emin, Emax, E_sp_eps   ! grid start and end; precision around a special point
   real(8), dimension(:), intent(in) :: special_point
   integer, intent(inout), optional :: Ngrid ! number of grid points
   real(8), dimension(:), intent(inout), optional :: array ! save grid
   !--------------
   integer :: SP_count, N
   real(8) :: E_cur, dE, dE_min
   logical :: point_is_here

   SP_count = 1
   N = 0
   dE_min = 0.1d0
   E_cur = Emin - dE_min  ! start from min
   do while (E_cur < Emax)
      N = N + 1   ! count points
      if (E_cur < 1.0d0-dE_min) then
         dE = dE_min
      else if (E_cur < 100.0d0) then
         dE = 1.0d0
      else
         dE = 10.0d0**(find_order_of_number(E_cur)-2) ! module "Little_subroutines"
      endif
      E_cur = E_cur + dE

      ! save grid points:
      if (present(array)) then
         if (N <= size(array)) then
            array(N) = E_cur
         else
            print*, 'Mismatch #1 of array size in go_thru_grid:', N, size(array)
            return
         endif
      endif

      ! Check if the special points are still there:
      point_is_here = .false. ! by default
      if (SP_count <= size(special_point)) then ! there may be a specila point:
         if ( (E_cur >= special_point(SP_count)) .and. &
              ( ((E_cur-dE) < special_point(SP_count)) .or. (N == 1) ) )  then ! special point inside interval
            SP_count = SP_count + 1 ! this special point is done, do the next one
            point_is_here = .true.
         endif
      endif

      ! Save grid point:
      if (point_is_here) then
         if (present(array)) then
            if (N+2 <= size(array)) then
               ! add two points:
               ! 1) below the special point:
               if (special_point(SP_count-1) - E_sp_eps > 0.0d0) then
                  N = N + 1
                  array(N) = special_point(SP_count-1) - E_sp_eps
               endif
               ! 2) above the special point:
               N = N + 1
               array(N) = special_point(SP_count-1) + E_sp_eps
            else
               print*, 'Mismatch #2 of array size in go_thru_grid', N, size(array)
               return
            endif
         else
            ! add extra points around the special point
            if (special_point(SP_count-1) - E_sp_eps < 0.0d0) then
               N = N + 1
            else ! 2 extra points, around the chosen one
               N = N + 2
            endif
         endif
      endif

   enddo ! while (E_cur .LT. Emax)

   ! prepare output:
   if (present(Ngrid)) Ngrid = N ! save grid size

   if (present(array)) then   ! make sure the array is sorted increasing
      call sort_array(array)  ! module "Algebra_tools"
   endif
end subroutine go_thru_grid



subroutine define_special_points(matter, special_point)
   type(Solid), intent(in) :: matter ! parameters of the material: ionization potentials to set special points
   real(8), dimension(:), allocatable :: special_point
   !-----------------------------
   integer :: Nsiz, i, k, sh, coun

   ! Count how many Ip's are there:
   Nsiz = 0
   do i = 1, size(matter%Atoms)  ! for all elements
      Nsiz = Nsiz + size(matter%Atoms(i)%Ip)
   enddo

   ! allocate and set special_points:
   allocate(special_point(Nsiz), source = 0.0d0)

   ! Copy the special points (ionization potentials):
   coun = 0
   do i = 1, size(matter%Atoms)  ! for all elements
      sh = size(matter%Atoms(i)%Ip)
      do k = 1, sh
         coun = coun + 1
         special_point(coun) = matter%Atoms(i)%Ip(k)
      enddo
   enddo

   ! Sort the array increasing:
   call sort_array(special_point)   ! module "Algebra_tools"

   ! Exclude copies of points, if any:
   call exclude_doubles(special_point) ! module "Little_subroutines"

!    do i = 1, size(special_point)  ! for all elements
!       print*, i, special_point(i)
!    enddo
end subroutine define_special_points




!BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB
! Ionization of an atom: BEB cross-section
! Eq.(57) from [Y.K.Kim, M.E.Rudd, Phys.Rev.A 50 (1994) 3954]
function Sigma_BEB(T, B_in, U, N)
   real(8) Sigma_BEB ! cross-section [A^2]
   real(8), intent(in) :: T   ! kinetic energy of incident electron [eV]
   real(8), intent(in) :: B_in   ! bindning energy of atomic electron [eV]
   real(8), intent(in) :: U   ! mean kinetic energy of electron in sub-shell [eV]
   real(8), intent(in) :: N   ! ocupation number of electrons in this shell
   real(8) t0, u0, S, B, eps

   eps = 1.0d-3   ! precision [eV]
   B = max(eps, B_in)  ! to exclude too small Ip or Egap

   if ((T .LE. B) .or. (abs(B) < 1.0d-6)) then
      Sigma_BEB = 0.0d0 ! [A^2] cross section for energies lower than the Ip
   else
      S = 4.0d0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      t0 = T/B ! energy normalized to the Rydberg constant ! Eq.(4)
      u0 = U/B
      if (t0 <= 0.0d0) then
         !print*, 'Sigma_BEB', t0
         Sigma_BEB = 0.0d0
      else
         Sigma_BEB = S/(t0+u0+1.0d0)*(log(t0)*0.5d0*(1.0d0-1.0d0/(t0*t0)) + (1.0d0-1.0d0/t0) -log(t0)/(t0+1.0d0)) ! [A^2]
      endif
   endif
end function Sigma_BEB


function dSigma_int_BEB(T, w, B_in, U, N)
   real(8), intent(in) :: w 	! [eV] transferred energy
   real(8), intent(in) :: T 	! [eV] kinetic energy of the electron
   real(8), intent(in) ::  B_in	! bindning energy of atomic electron [eV]
   real(8), intent(in) ::  U	! mean kinetic energy of electron in the sub-shell [eV]
   real(8), intent(in) ::  N	! ocupation number of electrons in this shell
   real(8) :: dSigma_int_BEB    ! differential cross-section

   real(8) t0, u0, w0, S, Sigma, dSigma0, B, eps

   eps = 1.0d-3   ! precision [eV]
   B = max(eps, B_in)  ! to exclude too small Ip or Egap

   S = 4.0d0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
   t0 = T/B    ! energy normalized to the Rydberg constant ! Eq.(4)
   u0 = U/B    ! kinetic energy normalized
   w0 = w/B    ! transferred energy normalized

   dSigma0 = dSigma_dw_int(S, t0, u0, 0.0d0) ! function see below
   dSigma_int_BEB = dSigma_dw_int(S, t0, u0, w0) - dSigma0
end function dSigma_int_BEB


function dSigma_dw_int(S, t0, u0, w0)
   real(8) t0, w0, u0, S, dSigma_dw_int
   dSigma_dw_int = S/(t0+u0+1.0d0) * (-(log(w0+1.0d0)-log(abs(t0-w0)))/(t0+1.0d0) + &
                  (1.0d0/(t0-w0)-1.0d0/(w0+1.0d0)) + &
                  log(t0)*0.5d0*(1.0d0/((t0-w0)*(t0-w0)) - &
                  1.0d0/((w0+1.0d0)*(w0+1.0d0))))
end function dSigma_dw_int


function dSigma_w_int_BEB(T, w, B_in, U, N)
   real(8), intent(in) :: w 	! [eV] transferred energy
   real(8), intent(in) :: T 	! [eV] kinetic energy of the electron
   real(8), intent(in) ::  B_in	! bindning energy of atomic electron [eV]
   real(8), intent(in) ::  U	! mean kinetic energy of electron in the sub-shell [eV]
   real(8), intent(in) ::  N	! ocupation number of electrons in this shell
   real(8) :: dSigma_w_int_BEB ! differential cross-section
   real(8) t0, u0, w0, S, Sigma, dSigma0, B, eps

   eps = 1.0d-3   ! precision [eV]
   B = max(eps, B_in)  ! to exclude too small Ip or Egap

   S = 4.0e0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
   t0 = T/B    ! energy normalized to the Rydberg constant ! Eq.(4)
   u0 = U/B    ! kinetic energy normalized
   w0 = w/B    ! transferred energy normalized

   dSigma0 = dSigma_dw_w_int(S, t0, u0, 0.0d0) ! function see below
   dSigma_w_int_BEB = B*(dSigma_dw_w_int(S, t0, u0, w0) - dSigma0)
end function dSigma_w_int_BEB


function dSigma_dw_w_int(S, t0, u0, w0)
   real(8) t0, w0, u0, S, dSigma_dw_w_int
   real(8) A, B, C, tw, logwt, overwt, overw1, overt1
   tw = 1.0d0/(t0+u0+1.0d0)
   logwt = log((w0+1.0d0)*abs(w0-t0))
   overwt = 1.0d0/(w0-t0)
   overw1 = 1.0d0/(w0+1.0d0)
   overt1 = 1.0d0/(t0+1.0d0)
   A = tw*overt1*overt1*(t0*t0*log(abs(w0-t0)) + t0*logwt + log(w0+1.0d0))
   B = tw*( -t0*overwt + overw1 + logwt)
   C = tw*log(t0)*(0.5d0*(overw1*overw1 + t0*overwt*overwt) + overwt - overw1)
   dSigma_dw_w_int = S*(A+B+C)
end function dSigma_dw_w_int
!BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB



!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Elastic electron scattering.
! Definition of elastic scattering we use: collisions in which
! only kinetic energy is exchanged without changing internal energy.
! In our case it means no ionization, which would change the system state.

subroutine Get_total_EMFP(numpar, matter) ! from partial EMFP get the total one
   type(Numerics_param), intent(in) :: numpar  ! all numerical parameters
   type(Solid), intent(inout) :: matter ! parameters of the material
   logical :: file_opened
   character(300) :: File_name
   integer :: i, FN
   ! Total elastic mean free path:
   write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Name))//'_Total_Electron_EMFP.txt'
   
   FN = 9695
   open(UNIT=FN, FILE = trim(adjustl(File_name)))
   do i = 1, size(matter%Atoms) ! for all atoms
      matter%El_EMFP_tot%L(:) = matter%El_EMFP_tot%L(:) + matter%Atoms(i)%El_EMFP%L(:) ! [1/A] inverse MFP
   enddo

   do i = 1, size( matter%El_EMFP_tot%E) ! save into the file EMFP
      write(FN,'(f25.16,es25.16)') matter%El_EMFP_tot%E(i), 1.0d0/matter%El_EMFP_tot%L(i) ! [eV], [A]
   enddo

   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) close(FN)
end subroutine Get_total_EMFP


subroutine Elastic_MFP(i_at, numpar, matter) ! elastic mean free path for scattering on the given atomic species
   integer, intent(in) :: i_at  ! number of atomic species in the array : matter%Atoms(i_at)
   type(Numerics_param), intent(in) :: numpar  ! all numerical parameters
   type(Solid), target, intent(inout) :: matter ! parameters of the material
   !-----------------------------------
   logical :: redo, file_opened, file_exists, read_well
   character(300) :: File_name
   integer :: N_grid, FN, k, count_lines, N, Reason
   real(8) :: Ele, L, sigma_el, Contrib
   real(8), pointer :: Zat
   
   Zat => matter%Atoms(i_at)%Z ! atomic number of an atom of the media
   
   N_grid = size(matter%Atoms(1)%El_EMFP%E) ! size of MFP arrays
   
   redo = .false. ! may be there is no need to recalculate MFPs
   
   if (i_at /= 1) then ! if the arrays are not set yet, do that
      if (.not.allocated(matter%Atoms(i_at)%El_EMFP%E)) allocate(matter%Atoms(i_at)%El_EMFP%E(N_grid))
      matter%Atoms(i_at)%El_EMFP%E = matter%Atoms(1)%El_EMFP%E
      if (.not.allocated(matter%Atoms(i_at)%El_EMFP%L)) allocate(matter%Atoms(i_at)%El_EMFP%L(N_grid))
   endif
   
   ! Set name of the file and its location (in the INPUT_FILES directory):
   write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Atoms(i_at)%Name))//'_Mott_Electron_EMFP.txt'
   FN = 9696      
   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if this file already exists
   if (file_exists) then ! EMFPs are already there
      open(UNIT=FN, FILE = trim(adjustl(File_name)))
      inquire(file=trim(adjustl(File_name)),opened=file_opened)
      if (.not.file_opened) then ! there is no such file, create a new one
         open(UNIT=FN, FILE = trim(adjustl(File_name)))
         inquire(file=trim(adjustl(File_name)),opened=file_opened)
         redo = .true. ! no data, need to recalculate the MFPs
         print*, 'Could not open file with MFP => recalculating MFP'
      else
         call Count_lines_in_file(FN, N) ! get how many lines in the file, module "Dealing_with_files"
         !If something is wrong in the file, number of lines is not the same as what we need:
         if (N .NE. N_grid) then ! replace file
            redo = .true. ! not enough data, need to recalculate the EMFPs
            print*, 'Energy grid mismatch in MFP file => recalculating MFP'
         endif
      endif
   else ! there is no file with EMFPs:
      open(UNIT=FN, FILE = trim(adjustl(File_name)))
      redo = .true. ! no data, need to recalculate the EMFPs
      if (numpar%verbose) print*, 'No file with MFP present => calculating MFP'
   endif
   
   ! If file exists and opened well:
   if (.not.redo) then ! just read the file
      count_lines = 0
      do k = 1, N_grid
         read(FN,*,IOSTAT=Reason) Ele, L ! read the line
         call read_file(Reason, count_lines, read_well) ! module "Dealing_with_files"
         if (.not. read_well) then ! file was not read properly
            redo = .true. ! no data, need to recalculate the MFPs
         endif
         matter%Atoms(i_at)%El_EMFP%E(k) = Ele ! [eV] energy
         matter%Atoms(i_at)%El_EMFP%L(k) = L   ! [A] MFP
      enddo
   endif   

   ! In case file doesn't exist or couldn't be read well:            
   if (redo) then ! recalculate the MFPs:
      write(0,'(a)') ' Calculating electron EMFP for '//trim(adjustl(matter%Atoms(i_at)%Name))//' atom'
      do k = 1, N_grid
         Ele = matter%Atoms(i_at)%El_EMFP%E(k) ! [eV] electron energy
         ! Cross section (Mott's):
         !call Atomic_elastic_sigma(Zat, Ele, sigma_el) ! [A^2] OLD
         sigma_el = Mott_total_CS(Ele, Zat)   ! [A^2] below (relativistic)

         Contrib = matter%At_Dens*1d-24*(matter%Atoms(i_at)%percentage)/SUM(matter%Atoms(:)%percentage)
         if (sigma_el <= 0.0d0) then ! no scattering
            L = 1d30 ! [A] infinite MFP
         else ! there is scattering
            L = 1.0d0/(Contrib*sigma_el) ! [A] EMFP
         endif
         matter%Atoms(i_at)%El_EMFP%L(k) = L ! [A] MFP
         write(FN,'(f25.16,es25.16)') Ele, L
         call print_progress('Progress:',k,N_grid)  ! module "Little_subroutines"
      enddo
      print*, 'Electron EMFPs are saved into file:', trim(adjustl(File_name))
   endif
   
   ! Inverse EMFP, as will be used in the program:
   matter%Atoms(i_at)%El_EMFP%L(:) = 1.0d0/matter%Atoms(i_at)%El_EMFP%L(:) ! [1/A] inverse MFP
   
   nullify(Zat) ! free the pointers
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) close(FN)
   !PAUSE 'Elastic_MFP'
end subroutine Elastic_MFP




! Mott's cross section according to [3]
pure function Mott_total_CS(Ee, Z, mass, mu_max_in) result(sigma)
   real(8) sigma	! [A^2] cross section
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   real(8), intent(in) :: Z	! Ion atomic number
   real(8), intent(in), optional :: mass    ! [me] incident particle mass
   real(8), intent(in), optional :: mu_max_in  ! mu=cos(theta), integration limit
   real(8) :: nu, beta, v, dbleZ, beta2, me, mu_max
   if (present(mass)) then
      me = mass*g_me
   else
      me = g_me
   endif
   if (present(mu_max_in)) then
      if (mu_max_in > 1.0d0) then ! total CS
         mu_max = 1.0d0
      else if (mu_max_in < -1.0d0) then ! zero
         mu_max = -1.0d0
      else  ! partially integrated CS, used for definition of the transfered energy
         mu_max = mu_max_in
      endif
   else ! total CS
      mu_max = 1.0d0
   endif
   v = velosity_from_kinetic_energy(Ee, me, afs=.false.)    ! [m/s] below
   if (v < 1.0d-6) then ! immobile particle does not scatter
      sigma = 0.0d0
   else
      beta = beta_factor(v)   ! module "Relativity"
      beta2 = beta*beta
      dbleZ = dble(Z)
      nu = screening_parameter(beta, dbleZ, Ee, me)   ! below
      sigma = g_Pi*g_r0*g_r0*dbleZ*(dbleZ + 1.0d0) * (mu_max + 1.0d0) / ((2.0d0*nu + 1.0d0 - mu_max)*(nu + 1.0d0)) * &
               (1.0d0 - beta2)/(beta2*beta2)	! [3] int of Eq.(39)
   endif
end function Mott_total_CS



subroutine NRG_transfer_elastic_atomic(Mat, Zat, Ee, dE) ! energy transfer in an elastic scattering on an atom:
   real(8), intent(in) :: Mat   ! mass of an atom of the media [kg]
   real(8), intent(in) :: Zat   ! atomic number of an atom of the media
   real(8), intent(in) :: Ee    ! [eV] incident electron energy
   real(8), intent(out) :: dE   ! [eV] transferred energy

   real(8) :: theta

   ! Sample polar angle according to the differential cross section
   call  get_electron_elastic_polar_angle(Ee, Zat, theta) ! below

   ! Change electron energy accordingly:
   dE = transfered_E_from_theta(Ee, theta, g_me, Mat) ! below
end subroutine NRG_transfer_elastic_atomic



function transfered_E_from_theta(Ekin, theta, M_in, mt) result(dE)
   real(8) dE   ! [eV]  according to Eq.(A.25) in 2015-edition of [1]
   real(8), intent(in) :: Ekin  ! [eV] kinetic energy of the incident particle
   real(8), intent(in) :: theta ! scattering angle
   real(8), intent(in) :: M_in, mt  ! incoming particle and target particle masses
   real(8) :: cos_theta, cos_theta2, sin_theta2, mc2, Mct2, Emc, E2mc, W1, W2, EmcMc
   mc2 = rest_energy(M_in)	! mc^2 [eV],  module "Relativity"
   Mct2 = rest_energy(mt)	! Mc^2 [eV],  module "Relativity"
   cos_theta = cos(theta)
   cos_theta2 = cos_theta*cos_theta
   sin_theta2 = 1.0d0 - cos_theta2
   Emc = Ekin + mc2
   E2mc = Ekin + 2.0d0*mc2
   EmcMc = Emc + Mct2
   W1 = Emc*sin_theta2 + Mct2 - cos_theta * sqrt( Mct2*Mct2 - mc2*mc2*sin_theta2 )
   W2 = Ekin*E2mc / ( EmcMc*EmcMc - Ekin*E2mc*cos_theta2 )
   dE = W1*W2
end function transfered_E_from_theta


subroutine get_electron_elastic_polar_angle(Ekin, Zat, theta)
   real(8), intent(in) :: Ekin  ! [eV] photon energy
   real(8), intent(in) :: Zat   ! element number
   real(8), intent(out) :: theta   ! polar angle of photoelectron emission
   real(8) :: RN
   ! Sample the angle:
   call random_number(RN)
   ! The scattering angle is:
   theta = Mott_sample_mu(Ekin, Zat, RN)  ! below
end subroutine get_electron_elastic_polar_angle


! Solution of diff.Mott's cross section:
pure function Mott_sample_mu(Ee, Z, RN, mass) result(theta)
   real(8) theta	! deflection angle
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   real(8), intent(in) :: Z	! Ion atomic number
   real(8), intent(in) :: RN    ! random number [0,1]
   real(8), intent(in), optional :: mass    ! in units of [me]
   real(8) :: nu, beta, v, beta2, me, mu
   if (present(mass)) then
      me = mass*g_me
   else
      me = g_me
   endif
   v = velosity_from_kinetic_energy(Ee, me, afs=.false.)    ! [m/s] below
   if (v < 1.0d-6) then ! immobile particle does not scatter
      theta = 0.0d0
   else
      beta = beta_factor(v)   ! below
      beta2 = beta*beta
      nu = screening_parameter(beta, Z, Ee, me)   ! below
      mu = (RN*(2.0d0*nu + 1.0d0) - nu) / (RN + nu)
      theta = ACOS(mu)
   endif
end function Mott_sample_mu


pure function screening_parameter(beta, Z, Ee, me) result(nu)
   real(8) nu	! screening parameter ! [4] Page 160, Eq. (7.8)
   real(8), intent(in) :: beta	! relativistic beta
   real(8), intent(in) :: Z	! atomic number
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   real(8), intent(in), optional :: me  ! [kg] mass of the particle
   real(8) :: beta2, tau, Erest
   beta2 = beta*beta
   if (present(me)) then
      Erest = rest_energy(me)	! module "Relativity"
   else
      Erest = rest_energy(g_me)	! module "Relativity"
   endif
   tau = Ee / Erest
   nu = 1.7d-5*Z**m_two_third*(1.0d0 - beta2)/beta2 * ( 1.13d0 + 3.76d0 * g_alpha*g_alpha/beta2 * Z*Z * sqrt(tau/(1.0d0 + tau)) )
end function screening_parameter



pure function velosity_from_kinetic_energy(Ekin, M0, afs) result(v)
   real(8) v        ! velosity  [A/fs], unless afs=.true. then [m/s]
   real(8), intent(in) :: Ekin	! kinetic energy [eV]
   real(8), intent(in) :: M0	! rest mass [kg]
   logical, intent(in), optional :: afs ! output velosity in [A/fs] (true), or in [m/s] (false)
   real(8) :: fact, Erest
   if (M0 < 1.0d-10*g_me) then   ! assume massless particle
      v = g_cvel    ! [m/s]
   else
      if (Ekin < 1.0d-15) then
         v = 0.0d0
      else
         Erest = rest_energy(M0) ! [eV] rest energy
         fact = Ekin/Erest + 1.0d0
         v = g_cvel*sqrt(1.0d0 - 1.0d0/(fact*fact))   ! [m/s]
      endif
   endif
   if (present(afs)) then
      if (afs) v = v * g_ms2Afs ! [m/s] -> [A/fs]
   else ! by default, use [A/fs]
      v = v * g_ms2Afs ! [m/s] -> [A/fs]
   endif
end function velosity_from_kinetic_energy


pure function beta_factor(v, afs_in) result(beta)
   real(8) beta
   real(8), intent(in) :: v	! [m/s] velosity
   logical, intent(in), optional :: afs_in ! v provided in [A/fs] if true, [m/s] otherwise
   real(8) :: conv
   if (present(afs_in)) then
      if (afs_in) then
         conv = g_Afs2ms ! [A/fs] -> [m/s]
      else
         conv = 1.0d0 ! v is in [m/s] by default
      endif
   else
      conv = 1.0d0  ! v is in [m/s] by default
   endif
   beta = v/g_cvel * conv
end function beta_factor


pure function rest_energy(M0) result(E)
   real(8) E	! [eV] total energy
   real(8), intent(in) :: M0	! [kg] rest mass of the particle
   E = M0*g_cvel*g_cvel/g_e	! [eV]
end function rest_energy



! Obsolete Mott's subroutines:
subroutine Atomic_elastic_sigma_OLD(Zat, Ee, sigma_el) ! total cross-section of elastic scattering on an atom:
   real(8), intent(in) :: Zat   ! atomic number of an atom of the media
   real(8), intent(in) :: Ee    ! [eV] incident electron energy
   REAL(8), INTENT(out) :: sigma_el   ! cross-section of elastic scteering [A^2]

   real(8) nc, pc, mec2e, Zat137, RyEe, beta2

   mec2e = g_me*g_cvel*g_cvel/g_e   ! a parameter enterring eq. below
   Zat137 = Zat/137.0d0   ! a parameter enterring eq. below
   RyEe = g_Ry/Ee         ! a parameter enterring eq. below
   beta2 = 2.0d0*Ee/(mec2e) ! v/c, a parameter enterring eq. below:
   ! Modified Moliere's screening parameter as given by Eq.(7.8) in
   ! "Monte Carlo Transport of Electrons and Photons" edited by T.M. Jenkins, W.R. Nelson, A. Rindi:
   pc = 1.7d-5*(Zat**(2.0d0/3.0d0))*(1.0d0-beta2)/beta2
   nc = pc*(1.13d0 + 3.76d0*(Zat137*Zat137)/beta2*sqrt(Ee/(Ee+mec2e)))
   ! Mott's cross-section with screening:
   sigma_el = g_Pi*g_a0*g_a0*Zat*(Zat+1.0d0)/(nc*(nc+1.0d0))*RyEe*RyEe ! [A^2]
end subroutine Atomic_elastic_sigma_OLD

subroutine NRG_transfer_elastic_atomic_OLD(Mat, Zat, Ee, dE) ! energy transfer in an elastic scattering on an atom:
   REAL(8), intent(in) :: Mat   ! mass of an atom of the media [kg]
   real(8), intent(in) :: Zat   ! atomic number of an atom of the media
   real(8), intent(in) :: Ee    ! [eV] incident electron energy
   REAL(8), INTENT(out) :: dE   ! [eV] transferred energy

   real(8) nc, pc, mec2e, Zat137, RyEe, Masses, beta2, RN

   call random_number(RN)

   mec2e = g_me*g_cvel*g_cvel/g_e   ! a parameter enterring eq. below
   Zat137 = Zat/137.0d0   ! a parameter enterring eq. below
   RyEe = g_Ry/Ee         ! a parameter enterring eq. below
   Masses = Mat + g_me    ! a parameter enterring eq. below
   beta2 = 2.0d0*Ee/(mec2e) ! v/c, a parameter enterring eq. below:
   ! Modified Moliere's screening parameter as given by Eq.(7.8) in
   ! "Monte Carlo Transport of Electrons and Photons" edited by T.M. Jenkins, W.R. Nelson, A. Rindi:
   pc = 1.7d-5*(Zat**(2.0d0/3.0d0))*(1.0d0-beta2)/beta2
   nc = pc*(1.13d0 + 3.76d0*(Zat137*Zat137)/beta2*sqrt(Ee/(Ee+mec2e)))
   dE = 4.0d0*Ee*g_me*Mat/(Masses*Masses)*(nc*(1.0d0-RN)/(nc+RN)) ! [eV]
end subroutine NRG_transfer_elastic_atomic_OLD






!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
!Photons:

subroutine get_photon_attenuation(matter, laser, numpar, Err)
   type(Solid), intent(inout) :: matter ! parameters of the material
   type(Pulse), dimension(:), intent(in) :: laser ! Laser pulse parameters
   type(Numerics_param), intent(in) :: numpar  ! all numerical parameters
   type(Error_handling), intent(inout) :: Err	! error save
   !===============================================
   real(8) Ele, L, temp1
   integer Nat, Nshl, i, j, k, N_grid, FN, FN2, N, Reason, count_lines, INFO, Z, indVB
   character(200) :: File_name, Folder_name, File_name_EPDL, Error_descript
   character(10) :: chtemp
   logical :: file_exist, file_exists, file_opened, redo, read_well
   integer, dimension(:), allocatable :: Shl_dsgnr
   real(8), dimension(:,:), allocatable :: Phot_abs_CS_tot
   real(8), dimension(:,:,:), allocatable :: Phot_abs_CS_shl
   real(8), dimension(:), allocatable :: Phot_abs_CS_VB
   
   Nat = size(matter%Atoms) ! number of atoms
   indVB = size(matter%Atoms(1)%Ip) ! the last shell of the first element is VB

   if (.not.allocated(matter%Atoms(1)%Ph_MFP(1)%E)) then
      call get_grid_4CS(N_grid, maxval(laser(:)%hw), matter%Atoms(1)%Ph_MFP(1)%E, matter) ! above
   else
      N_grid = size(matter%Atoms(1)%Ph_MFP(1)%E)
   endif
   if (.not.allocated(matter%Atoms(1)%Ph_MFP(1)%L)) allocate(matter%Atoms(1)%Ph_MFP(1)%L(N_grid))
   ! Total mean free paths:
   if (.not.allocated(matter%Ph_MFP_tot%L)) allocate(matter%Ph_MFP_tot%L(N_grid))
   if (.not.allocated(matter%Ph_MFP_tot%E)) allocate(matter%Ph_MFP_tot%E(N_grid))
   matter%Ph_MFP_tot%E = matter%Atoms(1)%Ph_MFP(1)%E
   matter%Ph_MFP_tot%L = 0.0d0

!    ATOMS:do i = 1, Nat ! for all atoms
   ATOMS:do i = Nat, 1, -1   ! for all atoms, counting backwards (to reach the VB the last)
      ! Before we start, make sure VB variables are allocated:
      if (.not.allocated(matter%Atoms(1)%Ph_MFP(indVB)%E)) then
         allocate(matter%Atoms(1)%Ph_MFP(indVB)%E(N_grid))
         matter%Atoms(1)%Ph_MFP(indVB)%E = matter%Atoms(1)%Ph_MFP(1)%E
      endif
      if (.not.allocated(matter%Atoms(1)%Ph_MFP(indVB)%L)) then
         allocate(matter%Atoms(1)%Ph_MFP(indVB)%L(N_grid))
         matter%Atoms(1)%Ph_MFP(indVB)%L = 0.0d0  ! just to start
      endif
      !------------------------------------------------------------------------------------
      ! Use new subroutines to deal with EPDL (faster versions):
      if (size(matter%Atoms(i)%TOCSph) > 0) then  ! if this element have more than VB
       select case (matter%Atoms(i)%TOCSph(1)) ! which inelastic cross section to use (BEB vs CDF)
       case (1) ! CDF
       case default ! BEB  
         Folder_name = trim(adjustl(numpar%input_path))//'Atomic_parameters'
         !File_name_EPDL = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(m_EPDL_file))
         File_name_EPDL = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(numpar%EPDL_file))
         inquire(file=trim(adjustl(File_name_EPDL)),exist=file_exist)
         if (.not.file_exist) then
            Error_descript = 'File '//trim(adjustl(File_name_EPDL))//' does not exist, the program terminates'
            call Save_error_details(Err, 1, Error_descript)
            print*, trim(adjustl(Error_descript))
            goto 9999
         endif
         inquire(file=trim(adjustl(File_name_EPDL)),opened=file_opened)
         if (.not. file_opened) then
            FN2 = 338
            open(FN2, FILE = trim(adjustl(File_name_EPDL)), status = 'old', action='read') ! EPDL
         else
            inquire(file=trim(adjustl(File_name_EPDL)),number=FN2) ! if it's already opened, just read from it
         endif
         Z = matter%Atoms(i)%Z ! atomic number
         if (allocated(Shl_dsgnr)) deallocate(Shl_dsgnr)
         allocate(Shl_dsgnr(size(matter%Atoms(i)%Shl_dsgnr_atomic)))
         Shl_dsgnr = matter%Atoms(i)%Shl_dsgnr_atomic
         
         if (allocated(Phot_abs_CS_tot)) deallocate(Phot_abs_CS_tot)
         if (allocated(Phot_abs_CS_shl)) deallocate(Phot_abs_CS_shl)
         ! Read cross sections for all photon energies at once:
         call Read_EPDL_data(numpar%path_sep, FN2, File_name_EPDL, INFO, Z, Shl_dsgnr, Phot_abs_CS_tot, Phot_abs_CS_shl) ! module "Dealing_with_EADL"
         rewind(FN2)    ! next element should start from the beginning of the file
!          do j = 1, size(Phot_abs_CS_shl,3)
!             print*, j, Phot_abs_CS_shl(1,1,j), Phot_abs_CS_shl(1,2,j)
!          enddo
!          pause 'Phot_abs_CS_shl'
!          
         Phot_abs_CS_tot(1,:) = Phot_abs_CS_tot(1,:)*1.0d6	! Convert [MeV] -> [eV]
         Phot_abs_CS_tot(2,:) = Phot_abs_CS_tot(2,:)*1.0d-8	! Convert [b] -> [A^2]
         Phot_abs_CS_shl(:,1,:) = Phot_abs_CS_shl(:,1,:)*1.0d6	! Convert [MeV] -> [eV]
         Phot_abs_CS_shl(:,2,:) = Phot_abs_CS_shl(:,2,:)*1.0d-8	! Convert [b] -> [A^2]
       end select
      endif ! (size(matter%Atoms(i)%TOCSph) > 0)
      !------------------------------------------------------------------------------------
   
      Nshl = size(matter%Atoms(i)%Ip)
      SHELLS:do j = 1, Nshl ! for all shells of this atom
         !redo = .false.
         redo = numpar%redo_MFP  ! user defined if MFP need to be recalculated


         if ((i .NE. 1) .or. (j .NE. 1)) then
            if (.not.allocated(matter%Atoms(i)%Ph_MFP(j)%E)) allocate(matter%Atoms(i)%Ph_MFP(j)%E(N_grid))
            matter%Atoms(i)%Ph_MFP(j)%E = matter%Atoms(1)%Ph_MFP(1)%E
            if (.not.allocated(matter%Atoms(i)%Ph_MFP(j)%L)) then
               allocate(matter%Atoms(i)%Ph_MFP(j)%L(N_grid))
               matter%Atoms(i)%Ph_MFP(j)%L = 0.0d0  ! just to start
            endif
         endif

         write(chtemp,'(i6)') INT(matter%Atoms(i)%Ip(j))
         if (j > 1) then    ! check if it's a degenerate level:
            if (INT(matter%Atoms(i)%Ip(j-1)) == INT(matter%Atoms(i)%Ip(j))) then
               write(chtemp,'(i6)') INT(matter%Atoms(i)%Ip(j) - 0.5)    ! artificially shift it a little bit to make it not exactly degenerate
            endif
         endif
         
         select case (matter%Atoms(i)%TOCSph(j)) ! which inelastic cross section to use (BEB vs CDF):
         case (1) ! CDF
            write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Atoms(i)%Name))//'_CDF_Photon_IMFP_Ip='//trim(adjustl(chtemp))//'eV.txt'
         case default ! BEB
            write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Atoms(i)%Name))//'_EPDL_Photon_IMFP_Ip='//trim(adjustl(chtemp))//'eV.txt'
         end select
         FN = 112+i*Nshl+j

         inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if this file already exists
         if (file_exists) then ! IMFPs are already there
            open(UNIT=FN, FILE = trim(adjustl(File_name)))
            inquire(file=trim(adjustl(File_name)),opened=file_opened)
            if (.not.file_opened) then ! there is no such file, create a new one
               open(UNIT=FN, FILE = trim(adjustl(File_name)))
               inquire(file=trim(adjustl(File_name)),opened=file_opened)
               redo = .true. ! no data, need to recalculate the MFPs
               print*, 'Could not open file with MFP => recalculating MFP'
            else
               call Count_lines_in_file(FN, N) ! get how many lines in the file
               if (N .NE. N_grid) then ! replace file
                  redo = .true. ! not enough data, need to recalculate the MFPs
                  print*, 'Energy grid mismatch in MFP file => recalculating MFP'
               endif
            endif
         else
            open(UNIT=FN, FILE = trim(adjustl(File_name)))
            redo = .true. ! no data, need to recalculate the MFPs
            if (numpar%verbose) print*, 'No file with MFP present => calculating MFP'
         endif

9900   if (redo) then ! recalculate the MFPs:
            write(0,'(a,f8.1)') ' Obtaining photon IMFP for '//trim(adjustl(matter%Atoms(i)%Name))//', Ip=', matter%Atoms(i)%Ip(j)
            
            select case (matter%Atoms(i)%TOCSph(j)) ! which inelastic cross section to use (BEB vs CDF):
            case (1) ! CDF
               ! Printout the calculated attenuations:
               do k = 1, N_grid ! for all grid-points
                 Ele = matter%Atoms(i)%Ph_MFP(j)%E(k)
                 call Tot_Phot_IMFP(Ele, matter, numpar, i, j, L, Err=Err)
                 matter%Atoms(i)%Ph_MFP(j)%L(k) = L ! [A] MFP
                 write(FN,'(f25.16,es25.16)') Ele, L
                 call print_progress('Progress:',k,N_grid)    ! module "Little_subroutines"
               enddo
               print*, 'Photon IMFPs are saved into file:', trim(adjustl(File_name))
            case default ! BEB
               if (.not.allocated(Phot_abs_CS_VB)) then
                  allocate(Phot_abs_CS_VB(size(matter%Atoms(1)%Ph_MFP(1)%E)))
               endif

               ! Check VB:
               if (j == Nshl) then  ! this is the outermost level:
                  Phot_abs_CS_VB = 0.0d0
                  if (i == 1) then ! there is VB explicitly
                     do k = j, size(Shl_dsgnr)  ! sum over all atomic shells that create VB:
                        call interpolate_data_on_grid(Phot_abs_CS_shl(k,1,:), Phot_abs_CS_shl(k,2,:), matter%Atoms(i)%Ph_MFP(j)%E, Phot_abs_CS_VB) ! module
                        temp1 = matter%At_Dens*1d-24*(matter%Atoms(i)%percentage)/SUM(matter%Atoms(:)%percentage) ! n_a
                        where(Phot_abs_CS_VB(:) > 1.0d-16)
                           Phot_abs_CS_VB(:) = 1.0d0/(temp1*Phot_abs_CS_VB(:)) ! IMFP [A]
                        elsewhere
                           Phot_abs_CS_VB(:) = 1.0d30 ! IMFP [A]
                        endwhere
                        where(matter%Atoms(i)%Ph_MFP(j)%L(:) > 1.0d-16)
                           matter%Atoms(i)%Ph_MFP(j)%L(:) = 1.0d0/ ( 1.0d0/matter%Atoms(i)%Ph_MFP(j)%L(:) + 1.0d0/Phot_abs_CS_VB(:)) ! IMFP [A]
                        elsewhere
                           matter%Atoms(i)%Ph_MFP(j)%L(:) = Phot_abs_CS_VB(:) ! IMFP [A]
                        endwhere
                     enddo ! k
                     
                      ! If it is the last shell of last atom, we know all VB now, printout the data:
!                      if ( (i == Nat) .and. (j == Nshl) ) then
                        write(chtemp,'(i6)') INT(matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)))
                        if (j > 1) then    ! check if it's a degenerate level:
                           if (INT(matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)-1)) == INT(matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)))) then
                              write(chtemp,'(i6)') INT(matter%Atoms(1)%Ip(size(matter%Atoms(1)%Ip)) - 0.5)    ! artificially shift it a little bit to make it not exactly degenerate
                           endif
                       endif
                       write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Atoms(1)%Name))//'_EPDL_Photon_IMFP_Ip='//trim(adjustl(chtemp))//'eV.txt'
                       open(UNIT=FN, FILE = trim(adjustl(File_name)))
                       ! Printout the calculated attenuations:
                       do k = 1, N_grid ! for all grid-points
                          Ele = matter%Atoms(i)%Ph_MFP(j)%E(k)
                          write(FN,'(f25.16,es25.16)') Ele, matter%Atoms(i)%Ph_MFP(j)%L(k)
                       enddo
                       print*, 'Photon IMFPs are saved into file:', trim(adjustl(File_name))
                       inquire(file=trim(adjustl(File_name)),opened=file_opened)
                       if (file_opened) close(FN)
!                      endif ! ( (i == Nat) .and. (j == Nshl) ) then
                  else ! (i == 1)
                     ! First, take care of the core level:
                     call interpolate_data_on_grid(Phot_abs_CS_shl(j,1,:), Phot_abs_CS_shl(j,2,:), matter%Atoms(i)%Ph_MFP(j)%E, matter%Atoms(i)%Ph_MFP(j)%L) ! module "Little_subroutines"
                     temp1 = matter%At_Dens*1d-24*(matter%Atoms(i)%percentage)/SUM(matter%Atoms(:)%percentage)
                     where(matter%Atoms(i)%Ph_MFP(j)%L(:) > 1.0d-16)
                        matter%Atoms(i)%Ph_MFP(j)%L(:) = 1.0d0/(temp1*matter%Atoms(i)%Ph_MFP(j)%L(:)) ! IMFP [A]
                     elsewhere
                        matter%Atoms(i)%Ph_MFP(j)%L(:) = 1d30    ! IMFP [A]
                     endwhere
                     ! Printout the calculated attenuations:
                     do k = 1, N_grid ! for all grid-points
                        Ele = matter%Atoms(i)%Ph_MFP(j)%E(k)
                        write(FN,'(f25.16,es25.16)') Ele, matter%Atoms(i)%Ph_MFP(j)%L(k)
                     enddo
                     print*, 'Photon IMFPs are saved into file:', trim(adjustl(File_name))
                     inquire(file=trim(adjustl(File_name)),opened=file_opened)
                     if (file_opened) close(FN)
         
                     ! And then add VB contributions:
                     do k = j+1, size(Shl_dsgnr)  ! sum over all atomic shells that create VB:
                        call interpolate_data_on_grid(Phot_abs_CS_shl(k,1,:), Phot_abs_CS_shl(k,2,:), matter%Atoms(1)%Ph_MFP(indVB)%E, Phot_abs_CS_VB) ! module
!                         temp1 = matter%At_Dens*1d-24*(matter%Atoms(i)%percentage)/SUM(matter%Atoms(:)%percentage) ! n_a
                        where(Phot_abs_CS_VB(:) > 1.0d-16)
                           Phot_abs_CS_VB(:) = 1.0d0/(temp1*Phot_abs_CS_VB(:)) ! IMFP [A]
                        elsewhere
                           Phot_abs_CS_VB(:) = 1.0d30 ! IMFP [A]
                        endwhere
                        where(matter%Atoms(1)%Ph_MFP(indVB)%L(:) > 1.0d-16)
                           matter%Atoms(1)%Ph_MFP(indVB)%L(:) = 1.0d0/ ( 1.0d0/matter%Atoms(1)%Ph_MFP(indVB)%L(:) + 1.0d0/Phot_abs_CS_VB(:)) ! IMFP [A]
                        elsewhere
                           matter%Atoms(1)%Ph_MFP(indVB)%L(:) = Phot_abs_CS_VB(:) ! IMFP [A]
                        endwhere
                     enddo ! k
                     
!                      do k = 1, size(matter%Atoms(1)%Ph_MFP(indVB)%L)
!                         print*, matter%Atoms(1)%Ph_MFP(indVB)%E(k), matter%Atoms(1)%Ph_MFP(indVB)%L(k)
!                      enddo
!                      pause 'matter%Atoms(1)%Ph_MFP(indVB)%L'
  
                  endif  ! (i == 1)
               else    ! core levels
                  call interpolate_data_on_grid(Phot_abs_CS_shl(j,1,:), Phot_abs_CS_shl(j,2,:), matter%Atoms(i)%Ph_MFP(j)%E, matter%Atoms(i)%Ph_MFP(j)%L) ! module "Little_subroutines"
                  
                  temp1 = matter%At_Dens*1d-24*(matter%Atoms(i)%percentage)/SUM(matter%Atoms(:)%percentage)
                  where(matter%Atoms(i)%Ph_MFP(j)%L(:) > 1.0d-16)
                     matter%Atoms(i)%Ph_MFP(j)%L(:) = 1.0d0/(temp1*matter%Atoms(i)%Ph_MFP(j)%L(:)) ! IMFP [A]
                  elsewhere
                     matter%Atoms(i)%Ph_MFP(j)%L(:) = 1d30    ! IMFP [A]
                  endwhere
                  ! Printout the calculated attenuations:
                  do k = 1, N_grid ! for all grid-points
                    Ele = matter%Atoms(i)%Ph_MFP(j)%E(k)
!                   call Tot_Phot_IMFP(Ele, matter, numpar, i, j, L, Err=Err)
!                   matter%Atoms(i)%Ph_MFP(j)%L(k) = L ! [A] MFP
                    write(FN,'(f25.16,es25.16)') Ele, matter%Atoms(i)%Ph_MFP(j)%L(k)
!                     print*, Ele, matter%Atoms(i)%Ph_MFP(j)%L(k)
!                     call print_progress('Progress:',k,N_grid)    ! module "Little_subroutines"
                  enddo
!                   pause 'Ele, matter%Atoms(i)%Ph_MFP(j)%L(k) TEST 1'
                  print*, 'Photon IMFPs are saved into file:', trim(adjustl(File_name))
               endif ! if (i == 1) then ! there is VB explicitly
            
            end select
         else ! Just read MFPs from the file:
            !print*, 'Reding photon IMFPs from file: ', trim(adjustl(File_name))
            count_lines = 0
            do k = 1, N_grid
               read(FN,*,IOSTAT=Reason)	Ele, L ! read the line
               call read_file(Reason, count_lines, read_well)
               if (.not. read_well) then
                  redo = .true. ! no data, need to recalculate the MFPs
                  goto 9900 ! if couldn't read the file, just recalculate it then
               endif
               matter%Atoms(i)%Ph_MFP(j)%E(k) = Ele ! [eV] energy
               matter%Atoms(i)%Ph_MFP(j)%L(k) = L   ! [A] MFP
            enddo
         endif
         inquire(file=trim(adjustl(File_name)),opened=file_opened)
         if (file_opened) close(FN)
      enddo SHELLS
   enddo ATOMS

   ! Total photon mean free path:
   write(File_name,'(a,a,a)') trim(adjustl(numpar%input_path)), trim(adjustl(matter%Name))//numpar%path_sep, trim(adjustl(matter%Name))//'_Total_photon_IMFP.txt'
   FN = FN + 1
   open(UNIT=FN, FILE = trim(adjustl(File_name)))
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl ! for all shells of this atom
         matter%Ph_MFP_tot%L(:) = matter%Ph_MFP_tot%L(:) + 1.0d0/matter%Atoms(i)%Ph_MFP(j)%L(:) ! [1/A] inverse MFP
      enddo
   enddo
   do i = 1, size(matter%Ph_MFP_tot%L)
      if (matter%Ph_MFP_tot%L(i) > 1.d-24) then
         matter%Ph_MFP_tot%L(i) = 1.0d0/matter%Ph_MFP_tot%L(i) ! [A]
      else ! infinity:
         matter%Ph_MFP_tot%L(i) = 1.0d30 ! [A] 
      endif
   enddo
   do k = 1, N_grid
      write(FN,'(f25.16,es25.16)') matter%Ph_MFP_tot%E(k), matter%Ph_MFP_tot%L(k)
   enddo

   !iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
   ! Make the inverse MFPs, since that's how we use them in MC:
   matter%Ph_MFP_tot%L = 1.0d0/matter%Ph_MFP_tot%L ! [1/A]
   do i = 1, size(matter%Atoms) ! for all atoms
      Nshl = size(matter%Atoms(i)%Ip)
      do j = 1, Nshl ! for all shells of this atom
         matter%Atoms(i)%Ph_MFP(j)%L(:) = 1.0d0/matter%Atoms(i)%Ph_MFP(j)%L(:) ! [1/A] inverse MFP
      enddo
   enddo
   !iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (file_opened) close(FN)
9999 continue
   if (allocated(Shl_dsgnr)) deallocate(Shl_dsgnr)
   if (allocated(Phot_abs_CS_tot)) deallocate(Phot_abs_CS_tot)
   if (allocated(Phot_abs_CS_shl)) deallocate(Phot_abs_CS_shl)
   if (allocated(Phot_abs_CS_VB)) deallocate(Phot_abs_CS_VB)
end subroutine get_photon_attenuation


subroutine Tot_Phot_IMFP(Ele, matter, numpar, Nat, Nshl, Sigma, Err, dEdx)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Solid), intent(inout) :: matter ! parameters of the material
    type(Numerics_param), intent(in) :: numpar 	! all numerical parameters
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma   ! IMFP [A]
    real(8), intent(out), optional :: dEdx   ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
    type(Error_handling), intent(inout) :: Err ! error save
    !=========================================
    integer i, j, n
    real(8) temp1, ImE, Sigma1
    character(200) :: Folder_name, File_name

    select case (matter%Atoms(Nat)%TOCSph(Nshl)) ! which inelastic cross section to use (EPDL97 vs CDF):
    case (1) ! CDF cross section
       if (Ele .LT. matter%Atoms(Nat)%Ip(Nshl)) then ! no ionization possible, IMFP = infinity
          Sigma = 1d30      ! [A] IMFP
          if (present(dEdx)) dEdx = Ele/Sigma  ! energy losses [eV/A]
       else
          !call Imewq(Ele, 0.0d0, Target_atoms, Nat, Nshl, ImE, Matter, NumPar=NumPar, photon=.true.) ! constructs full Im(-1/e(w,q)) as a sum of Drude-like functions
          call Imewq(matter%Atoms(Nat)%CDF, Ele, 0.0d0, Nshl, ImE) ! constructs Im(-1/e(w,q)) as a sum of Drude-like functions
          Sigma = g_cvel*g_h/(ImE*Ele*g_e)*1d10 ! [A] IMFP
          if (present(dEdx)) dEdx = Ele/Sigma    ! energy losses [eV/A]
       endif
    case default ! extract the cross section from EPDL97:
       if (matter%Atoms(Nat)%Shl_dsgnr(Nshl) .GE. 63) then ! it's VB
          Sigma1 = 0.0d0
          do i = 1, size(matter%Atoms) ! for all atoms, find outermost shells that form the VB:
             if (i == 1) then ! VB
                j = matter%Atoms(i)%Shl_dsgnr(matter%Atoms(i)%sh-1)
                call next_designator(n, j) ! find the designator for the VB (next shell after last one) ! module "Dealing_with_EADL"
                call get_photon_cross_section_EPDL(Ele, matter, numpar, i, Sigma, Err, Shl_dsgntr=n)    ! below
             else
                j = matter%Atoms(i)%Shl_dsgnr(matter%Atoms(i)%sh)
                call next_designator(n, j) ! find the designator for the VB (next shell after last one)  ! module "Dealing_with_EADL"
                call get_photon_cross_section_EPDL(Ele, matter, numpar, i, Sigma, Err, Shl_dsgntr=n)    ! below
                !call get_photon_cross_section_EPDL(Ele, matter, numpar, i, Sigma, Err, Nshl=Nshl)
             endif
             Sigma1 = Sigma1 + Sigma*(matter%Atoms(i)%percentage)/SUM(matter%Atoms(:)%percentage)
          enddo
          temp1 = matter%At_Dens*1d-24
          if ((temp1 .LE. 0.0d0) .or. (Sigma1 .LE. 0.0d0)) then
             Sigma = 1d30
          else
             Sigma = 1.0d0/(temp1*Sigma1) ! IMFP [A]
          endif
       else
          call get_photon_cross_section_EPDL(Ele, matter, numpar, Nat, Sigma, Err, Nshl=Nshl)   ! below
          temp1 = matter%At_Dens*1d-24*(matter%Atoms(Nat)%percentage)/SUM(matter%Atoms(:)%percentage)
          if ((temp1 .LE. 0.0d0) .or. (Sigma .LE. 0.0d0)) then
             Sigma = 1d30
          else
             Sigma = 1.0d0/(temp1*Sigma) ! IMFP [A]
          endif
       endif
    end select

!pause 'Tot_Phot_IMFP'
end subroutine Tot_Phot_IMFP


subroutine get_photon_cross_section_EPDL(Ele, matter, numpar, Nat, Sigma, Err, Shl_dsgntr, Nshl)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Solid), intent(inout) :: matter ! parameters of the material
    type(Numerics_param), intent(in) :: numpar ! all numerical parameters
    integer, intent(in) :: Nat ! number of atom
    real(8), intent(out) :: Sigma ! calculated inverse mean free path (out of cross section) [1/A]
    type(Error_handling), intent(inout) :: Err ! error save
    integer, intent(in), optional :: Nshl, Shl_dsgntr    ! shell number, shell designator
    !==========================================
    integer Z
    integer FN2
    character(100) File_name, Folder_name, Error_descript
    logical :: file_exist, File_opened
    !Folder_name = 'INPUT_EADL'  ! here we keep databases
    !File_name = trim(adjustl(Folder_name))//'/epdl97.all'

    Folder_name = trim(adjustl(numpar%input_path))//'Atomic_parameters'
    !File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//'epdl97.all'
    !File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(m_EPDL_file))
    File_name = trim(adjustl(Folder_name))//trim(adjustl(numpar%path_sep))//trim(adjustl(numpar%EPDL_file))
    inquire(file=trim(adjustl(File_name)),exist=file_exist)
    if (.not.file_exist) then
       Error_descript = 'File '//trim(adjustl(File_name))//' does not exist, the program terminates'
       call Save_error_details(Err, 1, Error_descript)
       print*, trim(adjustl(Error_descript))
       goto 9999
    endif
    inquire(file=trim(adjustl(File_name)),opened=File_opened)
    if (.not. File_opened) then
       FN2 = 338
       open(FN2, FILE = trim(adjustl(File_name)), status = 'old', action='read') ! EPDL97
    else
       inquire(file=trim(adjustl(File_name)),number=FN2) ! if it's already opened, just read from it
    endif

    Z = matter%Atoms(Nat)%Z ! atomic number
    shl_sign:if (present(Shl_dsgntr)) then
       if (Shl_dsgntr .GE. 62) then ! VB
          if (present(Nshl)) then
             call READ_EPDL_TYPE_FILE_real(FN2, File_name, Z, 73, 0, 91, Sigma, Eph=Ele, Shl_design=Shl_dsgntr, Last_dsgtr=matter%Atoms(Nat)%Shl_dsgnr(Nshl-1)) ! by subshells, module "Dealing_with_EADL"
                                         !(FN2, File_name, Z_needed, C_needed, I_needed, S_needed, Photoabs_sigma, Eph, Shl_design, Last_dsgtr, Photoabs_sigma_tot)
          else
             call READ_EPDL_TYPE_FILE_real(FN2, File_name, Z, 73, 0, 91, Sigma, Eph=Ele, Shl_design=Shl_dsgntr) ! by subshells, module "Dealing_with_EADL"
          endif
       else ! correct shell designator for the shell forming VB has already been set in the parent subroutine:
          call READ_EPDL_TYPE_FILE_real(FN2, File_name, Z, 73, 0, 91, Sigma, Eph=Ele, Shl_design=Shl_dsgntr) ! by subshells, module "Dealing_with_EADL"
       endif
    else shl_sign ! Nshl is present
       if (present(Nshl)) then ! by subshells
          call READ_EPDL_TYPE_FILE_real(FN2, File_name, Z, 73, 0, 91, Sigma, Eph=Ele, Shl_design=matter%Atoms(Nat)%Shl_dsgnr(Nshl)) ! module "Dealing_with_EADL"
       else
          print*, 'Niether Nshl nor Shl_dsgntr is present in get_photon_cross_section_EPDL'
          pause 'What to do?'
       endif
    endif shl_sign

9999 continue
end subroutine get_photon_cross_section_EPDL




END MODULE MC_cross_sections
