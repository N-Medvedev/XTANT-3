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
! This module contains subroutines for dealing with the Bond Order Potential TB parameters
! containing numerical sets of radial Hamiltonian and Overlap integrals
! Details are from
! [1]  https://arxiv.org/pdf/1909.04561.pdf

MODULE Dealing_with_BOP
use Universal_Constants   ! let it use universal constants
use Objects   ! since it uses derived types, it must know about them from module 'Objects'
use Dealing_with_files, only : read_file, close_file

implicit none
PRIVATE

! Modular parameters:
character(10), parameter :: m_BOP_directory = 'BOP_data'
character(10), parameter :: m_BOP_file = 'models.bx'
character(20), parameter :: m_BOP_model = 'non-orthogonal_TB'
character(50), parameter :: m_BOP_separator = '/-------------------------------------------------'
character(5), parameter ::  m_model_txt = 'model'  ! BOP marker of model type
character(4), parameter ::  m_bond_txt = 'bond'   ! BOP marker of atoms 1 and 2
character(7), parameter ::  m_basis_txt = 'valence'

character(30), parameter :: m_dimer_bond_length = 'INPUT_dimer_bond_lengths.dat'
character(14), parameter :: m_repulsive = '_repulsive.txt'

integer, parameter :: m_N_BOP_rep_grid = 400  ! number of grid points for BOP potential


public :: m_BOP_directory, m_BOP_file, read_BOP_parameters, idnetify_basis_size_BOP, &
          read_BOP_repulsive, check_if_repulsion_exists, m_repulsive, m_N_BOP_rep_grid


 contains



subroutine read_BOP_repulsive(TB_Repuls, i, j, Folder_name, path_sep, Name1, Name2, error_message)
   type(TB_Rep_BOP), dimension(:,:), intent(inout) ::  TB_Repuls    ! parameters of the repulsive potential
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   character(*), intent(in) :: Folder_name    ! directory where to find BOP parameters
   character(1), intent(in) :: path_sep
   character(*), intent(in) :: Name1, Name2 ! element names
   character(*), intent(inout) :: error_message
   !-------------------
   character(300) :: File_name
   logical :: file_exists, read_well
   integer :: FN_BL, Reason, count_lines, k

   error_message = ''    ! to start with

   ! File with repulsive BOP potential:
   File_name = trim(adjustl(Folder_name))//path_sep// &
      trim(adjustl(Name1))//'_'//trim(adjustl(Name2))//trim(adjustl(m_repulsive))   ! file with repulsive BOP parameters

   ! Check if file with repulsive part of BOP parameterization exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exists)

   if (.not.file_exists) then ! check reversed elements order:
      File_name = trim(adjustl(Folder_name))//path_sep// &
         trim(adjustl(Name2))//'_'//trim(adjustl(Name1))//trim(adjustl(m_repulsive))   ! file with repulsive BOP parameters
      ! Check if file with repulsive part of BOP parameterization exists:
      inquire(file=trim(adjustl(File_name)),exist=file_exists)
   endif

   if (file_exists) then
      FN_BL=112
      count_lines = 1
      open(UNIT=FN_BL, FILE = trim(adjustl(File_name)), status = 'old', action='read')
      read_well = .true.    ! to start with
      if (.not.allocated( TB_Repuls(i,j)%R)) allocate(TB_Repuls(i,j)%R(m_N_BOP_rep_grid))
      if (.not.allocated( TB_Repuls(i,j)%V_rep)) allocate(TB_Repuls(i,j)%V_rep(m_N_BOP_rep_grid))
      ! Read from the file:
      do k = 1, m_N_BOP_rep_grid    ! for all grid points
         read(FN_BL,*,IOSTAT=Reason) TB_Repuls(i,j)%R(k), TB_Repuls(i,j)%V_rep(k)
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
         if (.not. read_well) then
            error_message = trim(adjustl(error_message))//': '//trim(adjustl(File_name))
            goto 2013 ! exit
         endif
      enddo
      ! And the lower triangle of the matrix:
      if (j /= i) then
         if (.not.allocated( TB_Repuls(j,i)%R)) allocate(TB_Repuls(j,i)%R(m_N_BOP_rep_grid))
         if (.not.allocated( TB_Repuls(j,i)%V_rep)) allocate(TB_Repuls(j,i)%V_rep(m_N_BOP_rep_grid))
         TB_Repuls(j,i)%R(:) = TB_Repuls(i,j)%R(:)
         TB_Repuls(j,i)%V_rep(:) = TB_Repuls(i,j)%V_rep(:)
      endif
   else
      error_message = 'File with repulsive BOP parameters: '//trim(adjustl(File_name))//' not found'
   endif

2013 continue
   call close_file('close', FN=FN_BL) ! module "Dealing_with_files"
end subroutine read_BOP_repulsive


subroutine check_if_repulsion_exists(Elem1, Name1, Elem2, Name2, Folder_name, path_sep, file_exists, data_exists, bond_length, error_message)
   integer, intent(in) :: Elem1, Elem2  ! atomic numbers of the two elements we need the parameters for
   character(*), intent(in) :: Name1, Name2 ! element names
   character(*), intent(in) :: Folder_name    ! directory where to find BOP parameters
   character(1), intent(in) :: path_sep
   logical, intent(out) :: file_exists, data_exists
   real(8), intent(out) :: bond_length   ! [A] bond length for dimer
   character(*), intent(inout) :: error_message
   !--------------------------
   real(8), dimension(121,121) :: Bond_lengths  ! array with all bond lengths accross the Periodic Table
   real(8) :: eps
   integer :: FN_BL, Reason, count_lines
   character(300) :: File_name
   logical :: file_exist, read_well

   file_exists = .false. ! to start with
   data_exists = .false. ! to start with
   bond_length = 0.0d0   ! to start with
   error_message = ''    ! to start with
   eps = 1.0d-6 ! presision

   File_name = trim(adjustl(Folder_name))//path_sep// &
      trim(adjustl(Name1))//'_'//trim(adjustl(Name2))//trim(adjustl(m_repulsive))   ! file with repulsive BOP parameters

   ! Check if file with repulsive part of BOP parameterization exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exists)

   ! If there is no file, can new parameterization be created?
   ! For that, we ned the data on the equilibrium bond lengths, so check if they exist:
   if (.not.file_exists) then
      File_name = trim(adjustl(Folder_name))//path_sep//trim(adjustl(m_dimer_bond_length))   ! file with bond lengths
      inquire(file=trim(adjustl(File_name)),exist=file_exist)   ! make sure the file is there
      if (file_exist) then ! there is file with data, check if the required elements are present:
         FN_BL=111
         open(UNIT=FN_BL, FILE = trim(adjustl(File_name)), status = 'old', action='read')
         read(FN_BL,*,IOSTAT=Reason) Bond_lengths   ! read all array at once
         count_lines = 1
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
         if (read_well) then
            bond_length = Bond_lengths(Elem1, Elem2)
            if (bond_length > eps) then ! non-zero value for bond length, we can work with that
               data_exists = .true. ! we found the data we need for construction of the repulsive potential from BOP
            else ! check reversed order of elements:
               bond_length = Bond_lengths(Elem2, Elem1)
               if (bond_length > eps) then ! non-zero value for bond length, we can work with that
                  data_exists = .true. ! we found the data we need for construction of the repulsive potential from BOP
               else ! check reversed order of elements:
                  error_message = 'File '//trim(adjustl(File_name))//' does not have bond length for '//trim(adjustl(Name1))//'-'//trim(adjustl(Name2))
               endif
            endif
         endif
         ! Clean up:
         call close_file('close', FN=FN_BL) ! module "Dealing_with_files"
      else
         error_message = 'File with bond length: '//trim(adjustl(File_name))//' not found'
      endif ! (file_exist)
   endif ! (.not.file_exists)
end subroutine check_if_repulsion_exists


subroutine read_BOP_parameters(FN_BOP, Atom_1, Atom_2, TB_Hamil, KOA1, KOA2, error_message)
   integer, intent(in) :: FN_BOP     ! number of the file with BOP parameters (must be already open)
   character(*), intent(in) :: Atom_1, Atom_2       ! elements
   type (TB_H_BOP), dimension(:,:), intent(inout) :: TB_Hamil   ! All parameters of the Hamiltonian of TB for atom pair 1 and 2
   integer, intent(in) :: KOA1, KOA2    ! indices of atoms in the array
   character(*), intent(inout) :: error_message
   !-------------------------------------------
   integer :: count_lines, Reason, i, Vr_ind
   character(25) :: temp_txt, temp_txt2, temp_txt3, temp_txt4
!    character(25) :: model_txt, bond_txt, basis_txt
   logical :: read_well, found_atoms, correct_order
   
   ! To start with:
   temp_txt = ''
   temp_txt2 = ''
   temp_txt3 = ''
   temp_txt4 = ''
   count_lines = 0
   
!    print*, '=== BOP ==='
   
   MDL:do while (trim(adjustl(temp_txt3)) /= trim(adjustl(m_BOP_model)))
      read(FN_BOP,*,IOSTAT=Reason) temp_txt
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) goto 2012 ! exit
!       print*, trim(adjustl(temp_txt)), (trim(adjustl(temp_txt)) == trim(adjustl(model_txt)))
      
      if (trim(adjustl(temp_txt)) == trim(adjustl(m_model_txt))) then ! found the block containing the model parameters we need
         count_lines = count_lines - 1
         backspace(FN_BOP)  ! get back one line and read the model type
         
         read(FN_BOP,*,IOSTAT=Reason) temp_txt, temp_txt2, temp_txt3    ! read the model name
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
         if (.not.read_well) goto 2012 ! exit
!           print*, trim(adjustl(temp_txt3)), (trim(adjustl(temp_txt3)) == trim(adjustl(m_BOP_model)))
         
         if ( (trim(adjustl(temp_txt3)) == trim(adjustl(m_BOP_model))) ) then   ! found the right model, read parameters
            temp_txt = ''   ! reset to start reading new things
            
            ! Now, find the right elements:
            BND:do while (trim(adjustl(temp_txt)) /= trim(adjustl(m_bond_txt)))
               read(FN_BOP,*,IOSTAT=Reason) temp_txt
               call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
               if (.not.read_well) goto 2012 ! exit
!                 print*, 'One: ', trim(adjustl(temp_txt))
               
               ! If we reached another model without finding what we need, exit it:
               if (trim(adjustl(temp_txt)) == trim(adjustl(m_model_txt))) then
                  error_message = 'File with BOP TB parameterization does not contain required elements: '//trim(adjustl(Atom_1))//' and '//trim(adjustl(Atom_2))
                  print*, error_message
                  goto 2012 ! exit
               endif
               
               if (trim(adjustl(temp_txt)) == trim(adjustl(m_bond_txt))) then ! found the line containing elements
                  count_lines = count_lines - 1
                  backspace(FN_BOP)  ! get back one line and read the model type
                  read(FN_BOP,*,IOSTAT=Reason) temp_txt, temp_txt2, temp_txt3, temp_txt4    ! read the model name
                  call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
                  if (.not.read_well) goto 2012 ! exit
!                    print*, 'Two: ', trim(adjustl(temp_txt3)), trim(adjustl(temp_txt4)), trim(adjustl(Atom_1)), trim(adjustl(Atom_2))
                  
                  ! Check if found the correct pair of atoms (could be in any order):
                  if ( (trim(adjustl(temp_txt3)) == trim(adjustl(Atom_1))) .and. (trim(adjustl(temp_txt4)) == trim(adjustl(Atom_2)))  ) then
                     found_atoms = .true.
                     correct_order = .true.
                  elseif ( (trim(adjustl(temp_txt3)) == trim(adjustl(Atom_2))) .and. (trim(adjustl(temp_txt4)) == trim(adjustl(Atom_1)))  ) then
                     found_atoms = .true.
                     correct_order = .false.
                  else
                     found_atoms = .false.
                  endif
                  
                  ! In case the correct atomic pair is found:
                  if (found_atoms) then
                     
                     !print*, 'FOUND ATOMS:', trim(adjustl(temp_txt3)), ' ', trim(adjustl(temp_txt4)) 
                  
                     ! Set the starting values:
                     TB_Hamil(KOA1, KOA2)%H_ci = 0.0d0
                     TB_Hamil(KOA1, KOA2)%H_li = 0.0d0
                     TB_Hamil(KOA1, KOA2)%H_ni = 0.0d0
                     TB_Hamil(KOA1, KOA2)%E_ci = 80.0d0 ! default shift to high values
                     TB_Hamil(KOA1, KOA2)%E_li = 0.0d0
                     TB_Hamil(KOA1, KOA2)%E_ni = 0.0d0
                     TB_Hamil(KOA1, KOA2)%S_ci = 0.0d0
                     TB_Hamil(KOA1, KOA2)%S_li = 0.0d0
                     TB_Hamil(KOA1, KOA2)%S_ni = 0.0d0
                     if (KOA2 > KOA1) then  ! the back-order of the atoms is set here too (for not identical atoms)
                        TB_Hamil(KOA2, KOA1)%H_ci = 0.0d0
                        TB_Hamil(KOA2, KOA1)%H_li = 0.0d0
                        TB_Hamil(KOA2, KOA1)%H_ni = 0.0d0
                        TB_Hamil(KOA2, KOA1)%E_ci = 80.0d0 ! default shift to high values
                        TB_Hamil(KOA2, KOA1)%E_li = 0.0d0
                        TB_Hamil(KOA2, KOA1)%E_ni = 0.0d0
                        TB_Hamil(KOA2, KOA1)%S_ci = 0.0d0
                        TB_Hamil(KOA2, KOA1)%S_li = 0.0d0
                        TB_Hamil(KOA2, KOA1)%S_ni = 0.0d0
                     endif
                     
                     ! Now, read the given values from the file:
                     temp_txt = ' '  ! to restart reading
                     do while (trim(adjustl(temp_txt)) /= trim(adjustl(m_bond_txt)))
                        
                        read(FN_BOP,'(a4)',IOSTAT=Reason) temp_txt
                        call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
                        !error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
                        if (.not.read_well) goto 2012 ! exit
                        count_lines = count_lines - 1
                        backspace(FN_BOP)  ! get back one line and read it well again
                        
                        if (trim(adjustl(temp_txt)) == '/ E_') then ! it is a line with on-site parameters
                           read(FN_BOP,'(2X,a4)',IOSTAT=Reason) temp_txt
                        else ! it is any other line
                           read(FN_BOP,*,IOSTAT=Reason) temp_txt
                        endif
                        call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
                        !error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
                        if (.not.read_well) goto 2012 ! exit
                        
                        ! Interpret the marker:
                        select case (trim(adjustl(temp_txt)))
                        ! Hamiltonian parameters:
                        case ('ssSigma')  !Vr(1) = (s s sigma)
                           Vr_ind = 1
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('spSigma')  !Vr(2) = (s p sigma)
                           Vr_ind = 2
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('psSigma')  !Vr(3) = (p s sigma) INVERSE
                           Vr_ind = 3
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('sdSigma')  !Vr(4) = (s d sigma)
                           Vr_ind = 4
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('dsSigma')  !Vr(5) = (d s sigma) 
                           Vr_ind = 5
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('ppSigma')  !Vr(6) = (p p sigma)
                           Vr_ind = 6
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('ppPi')  !Vr(7) = (p p pi)
                           Vr_ind = 7
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('pdSigma')  !Vr(8) = (p d sigma)
                           Vr_ind = 8
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('dpSigma')  !Vr(9) = (d p sigma) INVERSE
                           Vr_ind = 9
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('pdPi')  !Vr(10) = (p d pi)
                           Vr_ind = 10
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('dpPi')  !Vr(11) = (d p pi) INVERSE
                           Vr_ind = 11
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('ddSigma')  !Vr(12) = (d d sigma)
                           Vr_ind = 12
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('ddPi')  !Vr(13) = (d d pi)
                           Vr_ind = 13
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('ddDelta')  !Vr(14) = (d d delta)
                           Vr_ind = 14
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        ! On-site parameters:
                        case ('E_s')  ! s
                           Vr_ind = 1
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_p0')  ! p sigma
                           Vr_ind = 2
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_p1')  ! p pi
                           Vr_ind = 3
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_d0')  ! d sigma
                           Vr_ind = 4
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_d1')  ! d pi
                           Vr_ind = 5
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_d2')  ! d delta
                           Vr_ind = 6
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                           
                        ! Overlap parameters:
                        case ('ssSigmaOverlap')  !Vr(1) = (s s sigma)
                           Vr_ind = 1
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('spSigmaOverlap')  !Vr(2) = (s p sigma)
                           Vr_ind = 2
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('psSigmaOverlap')  !Vr(2) = (p s sigma) INVERSE
                           Vr_ind = 3
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('sdSigmaOverlap')  !Vr(3) = (s d sigma)
                           Vr_ind = 4
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('dsSigmaOverlap')  !Vr(3) = (d s sigma)
                           Vr_ind = 5
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('ppSigmaOverlap')  !Vr(4) = (p p sigma)
                           Vr_ind = 6
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('ppPiOverlap')  !Vr(5) = (p p pi)
                           Vr_ind = 7
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('pdSigmaOverlap')  !Vr(6) = (p d sigma)
                           Vr_ind = 8
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('dpSigmaOverlap')  !Vr(6) = (d p sigma) INVERSE
                           Vr_ind = 9
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('pdPiOverlap')  !Vr(7) = (p d pi)
                           Vr_ind = 10
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('dpPiOverlap')  !Vr(7) = (d p pi) INVERSE
                           Vr_ind = 11
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('ddSigmaOverlap')  !Vr(8) = (d d sigma)
                           Vr_ind = 12
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('ddPiOverlap')  !Vr(9) = (d d pi)
                           Vr_ind = 13
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case ('ddDeltaOverlap')  !Vr(10) = (d d delta)
                           Vr_ind = 14
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        case (m_bond_txt)  ! new block of elements
!                            print*, 'Read everything, exit the cycle'
                           exit MDL    ! exit the cycle
                        endselect ! (trim(adjustl(temp_txt)))
                        
                     enddo
                  else
                     temp_txt = ''  ! to continue searching
                  endif ! ( (trim(adjustl(temp_txt3)) == trim(adjustl(Atom_1))) .and. (trim(adjustl(temp_txt4)) == trim(adjustl(Atom_2)))  )
               endif ! (trim(adjustl(temp_txt)) == trim(adjustl(bond_txt)))
            enddo BND
         endif  ! ( (trim(adjustl(temp_txt3)) == trim(adjustl(m_BOP_model))) )
      endif !   (trim(adjustl(temp_txt)) == trim(adjustl(m_model_txt)))
   enddo MDL
   
!    ! Test BOP parameteres reading:
!     print*,  'Test BOP parameteres reading: ', KOA1, KOA2,  Atom_1,' ', Atom_2
!     print*,  'H:', TB_Hamil(KOA1, KOA2)%H_ci, TB_Hamil(KOA1, KOA2)%H_li, TB_Hamil(KOA1, KOA2)%H_ni
!     print*,  'E 1-2:', TB_Hamil(KOA1, KOA2)%E_ci, 'li:', TB_Hamil(KOA1, KOA2)%E_li, 'ni:', TB_Hamil(KOA1, KOA2)%E_ni
!     print*,  'E 2-1:', TB_Hamil(KOA2, KOA1)%E_ci, 'li:', TB_Hamil(KOA2, KOA1)%E_li, 'ni:', TB_Hamil(KOA2, KOA1)%E_ni
!     print*,  'S:', TB_Hamil(KOA1, KOA2)%S_ci, TB_Hamil(KOA1, KOA2)%S_li, TB_Hamil(KOA1, KOA2)%S_ni
   
2012    continue
!    pause 'TEST read_BOP_parameters'
end subroutine read_BOP_parameters




subroutine read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, ord)
   integer, intent(in) :: FN_BOP    ! file to read from
   integer, intent(inout) :: count_lines
   logical, intent(inout) :: read_well
   character(*), intent(inout) :: error_message
   type (TB_H_BOP), dimension(:,:), intent(inout) :: TB_Hamil   ! All parameters of the Hamiltonian of TB for atom pair 1 and 2
   integer, intent(in) :: KOA1, KOA2, Vr_ind    ! indices of atoms in the array, and index of the shell
   logical, intent(in) :: ord  ! are atoms in correct oder? If present, meaning we are reading on-site energies, where order is important
   !------------------------------
   integer :: Reason, atom_ind, atom_given, i, j
   character(20) :: temp_txt, temp_txt2, temp_txt3, temp_txt4, temp_txt5
   real(8) :: eps
   
   atom_ind = -1    ! to start with
   
   count_lines = count_lines - 1
   backspace(FN_BOP)    ! read the same line again, but this time read all the parameters
   
   ! Read parameters in the BOP format:
   read(FN_BOP,'(2X, a4, a5, i2)',IOSTAT=Reason)  temp_txt, temp_txt2, atom_ind    ! e.g.   / E_s atom 0        =  sum7exp
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not.read_well) then
      error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
      print*,  temp_txt, temp_txt2, atom_ind, Reason
      print*, 'Error BOP3: ', error_message
   endif
   count_lines = count_lines - 1
   backspace(FN_BOP)  ! to read it again

   if (ord) then ! correct order of atoms:
      i = KOA1
      j = KOA2
   else ! inverse order
      i = KOA2
      j = KOA1
   endif

   if (atom_ind == 0) then    ! parameters of the first atom
      read(FN_BOP, '(2x)', advance='no')  ! to exclude the slash, which is interpreted by FORTRAN as termination statement
      read(FN_BOP, *, IOSTAT=Reason) temp_txt, temp_txt2, atom_ind, temp_txt4, temp_txt5, &  ! e.g.   / E_s atom 0        =  sum7exp
                            TB_Hamil(i,j)%E_ci(Vr_ind,1), TB_Hamil(i,j)%E_li(Vr_ind,1), TB_Hamil(i,j)%E_ni(Vr_ind,1), &
                            TB_Hamil(i,j)%E_ci(Vr_ind,2), TB_Hamil(i,j)%E_li(Vr_ind,2), TB_Hamil(i,j)%E_ni(Vr_ind,2), &
                            TB_Hamil(i,j)%E_ci(Vr_ind,3), TB_Hamil(i,j)%E_li(Vr_ind,3), TB_Hamil(i,j)%E_ni(Vr_ind,3), &
                            TB_Hamil(i,j)%E_ci(Vr_ind,4), TB_Hamil(i,j)%E_li(Vr_ind,4), TB_Hamil(i,j)%E_ni(Vr_ind,4), &
                            TB_Hamil(i,j)%E_ci(Vr_ind,5), TB_Hamil(i,j)%E_li(Vr_ind,5), TB_Hamil(i,j)%E_ni(Vr_ind,5), &
                            TB_Hamil(i,j)%E_ci(Vr_ind,6), TB_Hamil(i,j)%E_li(Vr_ind,6), TB_Hamil(i,j)%E_ni(Vr_ind,6), &
                            TB_Hamil(i,j)%E_ci(Vr_ind,7), TB_Hamil(i,j)%E_li(Vr_ind,7), TB_Hamil(i,j)%E_ni(Vr_ind,7)
      ! Replace zeros with shifted-away energies:
!       eps = 1.0d-12
!       where (ABS(TB_Hamil(i,j)%E_ci(Vr_ind,:)) < eps)
!          TB_Hamil(i,j)%E_ci(Vr_ind,:) = 80.0d0
!       endwhere
   else ! 1
      read(FN_BOP, '(2x)', advance='no')  ! to exclude the slash, which is interpreted by FORTRAN as termination statement
      read(FN_BOP, *, IOSTAT=Reason) temp_txt, temp_txt2, atom_ind, temp_txt4, temp_txt5, &  ! e.g.   / E_s atom 1        =  sum7exp
                            TB_Hamil(j,i)%E_ci(Vr_ind,1), TB_Hamil(j,i)%E_li(Vr_ind,1), TB_Hamil(j,i)%E_ni(Vr_ind,1), &
                            TB_Hamil(j,i)%E_ci(Vr_ind,2), TB_Hamil(j,i)%E_li(Vr_ind,2), TB_Hamil(j,i)%E_ni(Vr_ind,2), &
                            TB_Hamil(j,i)%E_ci(Vr_ind,3), TB_Hamil(j,i)%E_li(Vr_ind,3), TB_Hamil(j,i)%E_ni(Vr_ind,3), &
                            TB_Hamil(j,i)%E_ci(Vr_ind,4), TB_Hamil(j,i)%E_li(Vr_ind,4), TB_Hamil(j,i)%E_ni(Vr_ind,4), &
                            TB_Hamil(j,i)%E_ci(Vr_ind,5), TB_Hamil(j,i)%E_li(Vr_ind,5), TB_Hamil(j,i)%E_ni(Vr_ind,5), &
                            TB_Hamil(j,i)%E_ci(Vr_ind,6), TB_Hamil(j,i)%E_li(Vr_ind,6), TB_Hamil(j,i)%E_ni(Vr_ind,6), &
                            TB_Hamil(j,i)%E_ci(Vr_ind,7), TB_Hamil(j,i)%E_li(Vr_ind,7), TB_Hamil(j,i)%E_ni(Vr_ind,7)
                            
!       ! Replace zeros with shifted-away energies:
!       eps = 1.0d-12
!       where (ABS(TB_Hamil(j,i)%E_ci(Vr_ind,:)) < eps)
!          TB_Hamil(j,i)%E_ci(Vr_ind,:) = 80.0d0
!       endwhere
   endif
   
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not.read_well) then
      error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
      print*, 'Error BOP4: ', error_message
      print*, temp_txt, temp_txt2, atom_ind, temp_txt4, temp_txt5
   endif

end subroutine read_line_with_BOP_onsite




subroutine read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, ci, li, ni, ord, atom0, invers)
   integer, intent(in) :: FN_BOP    ! file to read from
   integer, intent(inout) :: count_lines
   logical, intent(inout) :: read_well
   character(*), intent(inout) :: error_message
   real(8), dimension(:), intent(inout) :: ci, li, ni
   logical, intent(in), optional :: ord  ! are atoms in correct oder? If present, meaning we are reading on-site energies, where order is important
   integer, intent(in), optional :: atom0
   logical, intent(in), optional :: invers  ! to change the sign of ci or not?
   !------------------------------
   integer :: Reason, atom_ind, atom_given
   character(20) :: temp_txt, temp_txt2, temp_txt3, temp_txt4, temp_txt5
   real(8) :: eps
   
   atom_ind = -1    ! to start with
   
   count_lines = count_lines - 1
   backspace(FN_BOP)    ! read the same line again, but this time read all the parameters
   
   ! Read parameters in the BOP format:
   if (present(ord)) then ! format of line with on-site energy parameters:
      read(FN_BOP,'(2X, a4, a5, i2)',IOSTAT=Reason)  temp_txt, temp_txt2, atom_ind    ! e.g.   / E_s atom 0        =  sum7exp
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) then
         error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
         print*,  temp_txt, temp_txt2, atom_ind, Reason
         print*, 'Error BOP1: ', error_message
      endif
      count_lines = count_lines - 1
      backspace(FN_BOP)  ! to read it again

      if (present(atom0)) then
         atom_given = atom0
      else
         atom_given = 0    ! be default
      endif
   
      if (atom_ind == atom_given) then    ! parameters of the first atom
         read(FN_BOP, '(2x)', advance='no')  ! to exclude the slash, which is interpreted by FORTRAN as termination statement
         read(FN_BOP, *, IOSTAT=Reason) temp_txt, temp_txt2, atom_ind, temp_txt4, temp_txt5, &  ! e.g.   / E_s atom 0        =  sum7exp
                               ci(1), li(1), ni(1), &
                               ci(2), li(2), ni(2), &
                               ci(3), li(3), ni(3), &
                               ci(4), li(4), ni(4), &
                               ci(5), li(5), ni(5), &
                               ci(6), li(6), ni(6), &
                               ci(7), li(7), ni(7)
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
         if (.not.read_well) then
            error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
            print*, 'Error BOP2: ', error_message
            print*, temp_txt, temp_txt2, atom_ind, temp_txt4, temp_txt5
         endif
      else
         read(FN_BOP,*,IOSTAT=Reason)
      endif
      
      ! Replace zeros with shifted-away energies:
!       eps = 1.0d-12
!       where (ABS(ci(:)) < eps)
!          ci(:) = 80.0d0
!       endwhere
   
   else ! format of line with H or S parameters:
      read(FN_BOP,*,IOSTAT=Reason)  temp_txt, temp_txt2, temp_txt3, &  ! e.g.   ssSigma        =  sum7exp
                               ci(1), li(1), ni(1), &
                               ci(2), li(2), ni(2), &
                               ci(3), li(3), ni(3), &
                               ci(4), li(4), ni(4), &
                               ci(5), li(5), ni(5), &
                               ci(6), li(6), ni(6), &
                               ci(7), li(7), ni(7)
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) then
         error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
         print*, 'Error BOP2: ', error_message
         print*, temp_txt, temp_txt2, atom_ind, temp_txt4, temp_txt5
      endif
      ! Change the sign to convert into the K-S format we use in this code:
      !if (present(invers)) then  ! correct
      if (present(invers)) then
!          ci(:) = -ci(:)
         if (invers) ci(:) = -ci(:)
      endif
   endif
end subroutine read_line_with_BOP_params


pure  subroutine idnetify_basis_size_BOP(TB_Hamil, Nsiz)
   integer, intent(out) :: Nsiz  ! index: 0=s, 1=sp3, 2=sp3d5
   type(TB_H_BOP), intent(in) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   ! Check if d-band is used for any element:
   if ( minval(TB_Hamil%E_ci(4:6,:)) > 79.0d0 ) then  ! The basis is smaller than sp3d5
      if ( minval(TB_Hamil%E_ci(2:3,:)) > 79.0d0 ) then  ! The basis is smaller than sp3
         Nsiz = 0   ! s
      else
         Nsiz = 1   ! sp3
      endif
   else
      Nsiz = 2  ! sp3d5
   endif
end subroutine idnetify_basis_size_BOP




!-------------------------------------------------------------------------------------------------------------------
! Obsolete subroutines:


subroutine read_BOP_parameters_OLD(FN_BOP, Atom_1, Atom_2, TB_Hamil, KOA1, KOA2, error_message)
   integer, intent(in) :: FN_BOP     ! number of the file with BOP parameters (must be already open)
   character(*), intent(in) :: Atom_1, Atom_2       ! elements
   type (TB_H_BOP), dimension(:,:), intent(inout) :: TB_Hamil   ! All parameters of the Hamiltonian of TB for atom pair 1 and 2
   integer, intent(in) :: KOA1, KOA2    ! indices of atoms in the array
   character(*), intent(inout) :: error_message
   !-------------------------------------------
   integer :: count_lines, Reason, i, Vr_ind
   character(25) :: temp_txt, temp_txt2, temp_txt3, temp_txt4
!    character(25) :: model_txt, bond_txt, basis_txt
   logical :: read_well, found_atoms, correct_order
   
   ! To start with:
   temp_txt = ''
   temp_txt2 = ''
   temp_txt3 = ''
   temp_txt4 = ''
   count_lines = 0
   
!    print*, '=== BOP ==='
   
   MDL:do while (trim(adjustl(temp_txt3)) /= trim(adjustl(m_BOP_model)))
      read(FN_BOP,*,IOSTAT=Reason) temp_txt
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) goto 2012 ! exit
!       print*, trim(adjustl(temp_txt)), (trim(adjustl(temp_txt)) == trim(adjustl(model_txt)))
      
      if (trim(adjustl(temp_txt)) == trim(adjustl(m_model_txt))) then ! found the block containing the model parameters we need
         count_lines = count_lines - 1
         backspace(FN_BOP)  ! get back one line and read the model type
         
         read(FN_BOP,*,IOSTAT=Reason) temp_txt, temp_txt2, temp_txt3    ! read the model name
         call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
         if (.not.read_well) goto 2012 ! exit
!           print*, trim(adjustl(temp_txt3)), (trim(adjustl(temp_txt3)) == trim(adjustl(m_BOP_model)))
         
         if ( (trim(adjustl(temp_txt3)) == trim(adjustl(m_BOP_model))) ) then   ! found the right model, read parameters
            temp_txt = ''   ! reset to start reading new things
            
            ! Now, find the right elements:
            BND:do while (trim(adjustl(temp_txt)) /= trim(adjustl(m_bond_txt)))
               read(FN_BOP,*,IOSTAT=Reason) temp_txt
               call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
               if (.not.read_well) goto 2012 ! exit
!                 print*, 'One: ', trim(adjustl(temp_txt))
               
               ! If we reached another model without finding what we need, exit it:
               if (trim(adjustl(temp_txt)) == trim(adjustl(m_model_txt))) then
                  error_message = 'File with BOP TB parameterization does not contain required elements: '//trim(adjustl(Atom_1))//' and '//trim(adjustl(Atom_2))
                  print*, error_message
                  goto 2012 ! exit
               endif
               
               if (trim(adjustl(temp_txt)) == trim(adjustl(m_bond_txt))) then ! found the line containing elements
                  count_lines = count_lines - 1
                  backspace(FN_BOP)  ! get back one line and read the model type
                  read(FN_BOP,*,IOSTAT=Reason) temp_txt, temp_txt2, temp_txt3, temp_txt4    ! read the model name
                  call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
                  if (.not.read_well) goto 2012 ! exit
!                    print*, 'Two: ', trim(adjustl(temp_txt3)), trim(adjustl(temp_txt4)), trim(adjustl(Atom_1)), trim(adjustl(Atom_2))
                  
                  ! Check if found the correct pair of atoms (could be in any order):
                  if ( (trim(adjustl(temp_txt3)) == trim(adjustl(Atom_1))) .and. (trim(adjustl(temp_txt4)) == trim(adjustl(Atom_2)))  ) then
                     found_atoms = .true.
                     correct_order = .true.
                  elseif ( (trim(adjustl(temp_txt3)) == trim(adjustl(Atom_2))) .and. (trim(adjustl(temp_txt4)) == trim(adjustl(Atom_1)))  ) then
                     found_atoms = .true.
                     correct_order = .false.
                  else
                     found_atoms = .false.
                  endif
                  
                  ! In case the correct atomic pair is found:
                  if (found_atoms) then
                     
                     !print*, 'FOUND ATOMS:', trim(adjustl(temp_txt3)), ' ', trim(adjustl(temp_txt4)) 
                  
                     ! Set the starting values:
                     TB_Hamil(KOA1, KOA2)%H_ci = 0.0d0
                     TB_Hamil(KOA1, KOA2)%H_li = 0.0d0
                     TB_Hamil(KOA1, KOA2)%H_ni = 0.0d0
                     TB_Hamil(KOA1, KOA2)%E_ci = 80.0d0 ! default shift to high values
                     TB_Hamil(KOA1, KOA2)%E_li = 0.0d0
                     TB_Hamil(KOA1, KOA2)%E_ni = 0.0d0
                     TB_Hamil(KOA1, KOA2)%S_ci = 0.0d0
                     TB_Hamil(KOA1, KOA2)%S_li = 0.0d0
                     TB_Hamil(KOA1, KOA2)%S_ni = 0.0d0
                     if (KOA2 > KOA1) then  ! the back-order of the atoms is set here too (for not identical atoms)
                        TB_Hamil(KOA2, KOA1)%H_ci = 0.0d0
                        TB_Hamil(KOA2, KOA1)%H_li = 0.0d0
                        TB_Hamil(KOA2, KOA1)%H_ni = 0.0d0
                        TB_Hamil(KOA2, KOA1)%E_ci = 80.0d0 ! default shift to high values
                        TB_Hamil(KOA2, KOA1)%E_li = 0.0d0
                        TB_Hamil(KOA2, KOA1)%E_ni = 0.0d0
                        TB_Hamil(KOA2, KOA1)%S_ci = 0.0d0
                        TB_Hamil(KOA2, KOA1)%S_li = 0.0d0
                        TB_Hamil(KOA2, KOA1)%S_ni = 0.0d0
                     endif
                     
                     ! Now, read the given values from the file:
                     temp_txt = ' '  ! to restart reading
                     do while (trim(adjustl(temp_txt)) /= trim(adjustl(m_bond_txt)))
                        
                        read(FN_BOP,'(a4)',IOSTAT=Reason) temp_txt
                        call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
                        !error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
                        if (.not.read_well) goto 2012 ! exit
                        count_lines = count_lines - 1
                        backspace(FN_BOP)  ! get back one line and read it well again
                        
                        if (trim(adjustl(temp_txt)) == '/ E_') then ! it is a line with on-site parameters
                           read(FN_BOP,'(2X,a4)',IOSTAT=Reason) temp_txt
                        else ! it is any other line
                           read(FN_BOP,*,IOSTAT=Reason) temp_txt
                        endif
                        call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
                        !error_message = trim(adjustl(error_message))//'; file with BOP TB parameterization'
                        if (.not.read_well) goto 2012 ! exit
                        
!                         print*, 'Three: ', trim(adjustl(temp_txt))
                  
                        ! Interpret the marker:
                        select case (trim(adjustl(temp_txt)))
                        ! Hamiltonian parameters:
                        case ('ssSigma')  !Vr(1) = (s s sigma)
                           Vr_ind = 1
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                        case ('spSigma')  !Vr(2) = (s p sigma)
                        !case ('psSigma')  !Vr(3) = (p s sigma) 
                           Vr_ind = 2
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                          
                        case ('psSigma')  !Vr(3) = (p s sigma) INVERSE
                        !case ('spSigma')  !Vr(2) = (s p sigma)
                           Vr_ind = 2
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:), invers=.true.)   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:), invers=.true.)   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('sdSigma')  !Vr(3) = (s d sigma)
                           Vr_ind = 3
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                          
                        case ('dsSigma')  !Vr(3) = (d s sigma) 
                           Vr_ind = 3
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('ppSigma')  !Vr(4) = (p p sigma)
                           Vr_ind = 4
                           !Vr_ind = 5   ! test
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                           
                        case ('ppPi')  !Vr(5) = (p p pi)
                           Vr_ind = 5
                           !Vr_ind = 4   ! test
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:)
                           endif
                           
                        case ('pdSigma')  !Vr(6) = (p d sigma)
                           Vr_ind = 6
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                           
                        case ('dpSigma')  !Vr(6) = (d p sigma) INVERSE
                           Vr_ind = 6
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:), invers=.true.)   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:), invers=.true.)   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('pdPi')  !Vr(7) = (p d pi)
                           Vr_ind = 7
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                          
                        case ('dpPi')  !Vr(7) = (d p pi) INVERSE
                           Vr_ind = 7
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:), invers=.true.)   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:), invers=.true.)   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('ddSigma')  !Vr(8) = (d d sigma)
                           Vr_ind = 8
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                          
                        case ('ddPi')  !Vr(9) = (d d pi)
                           Vr_ind = 9
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                          
                        case ('ddDelta')  !Vr(10) = (d d delta)
                           Vr_ind = 10
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%H_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%H_ni(Vr_ind,:))   ! below
                          else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%H_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%H_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                           
                        ! On-site parameters:
                        case ('E_s')  ! s
                           Vr_ind = 1
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_p0')  ! p sigma
                           Vr_ind = 2
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_p1')  ! p pi
                           Vr_ind = 3
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_d0')  ! d sigma
                           Vr_ind = 4
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_d1')  ! d pi
                           Vr_ind = 5
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                        case ('E_d2')  ! d delta
                           Vr_ind = 6
                           call read_line_with_BOP_onsite(FN_BOP, count_lines, read_well, error_message, TB_Hamil, KOA1, KOA2, Vr_ind, correct_order) ! below
                           if (.not.read_well) goto 2012 ! exit
                           
                        ! Overlap parameters:
                        case ('ssSigmaOverlap')  !Vr(1) = (s s sigma)
                           Vr_ind = 1
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case ('spSigmaOverlap')  !Vr(2) = (s p sigma)
                        !case ('psSigmaOverlap')  !Vr(2) = (p s sigma) INVERSE
                           Vr_ind = 2
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                           
                        case ('psSigmaOverlap')  !Vr(2) = (p s sigma) INVERSE
                        !case ('spSigmaOverlap')  !Vr(2) = (s p sigma)
                           Vr_ind = 2
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:), invers=.true.)   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:), invers=.true.)   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('sdSigmaOverlap')  !Vr(3) = (s d sigma)
                           Vr_ind = 3
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                          else
                             call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:))   ! below
                          endif
                          if (.not.read_well) goto 2012 ! exit
                           
                        case ('dsSigmaOverlap')  !Vr(3) = (d s sigma)
                           Vr_ind = 3
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:))   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('ppSigmaOverlap')  !Vr(4) = (p p sigma)
                           Vr_ind = 4
                           !Vr_ind = 5   ! test
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                        
                        case ('ppPiOverlap')  !Vr(5) = (p p pi)
                           Vr_ind = 5
                           !Vr_ind = 4   ! test
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case ('pdSigmaOverlap')  !Vr(6) = (p d sigma)
                           Vr_ind = 6
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case ('dpSigmaOverlap')  !Vr(6) = (d p sigma) INVERSE
                           Vr_ind = 6
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:), invers=.true.)   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:), invers=.true.)   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('pdPiOverlap')  !Vr(7) = (p d pi)
                           Vr_ind = 7
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case ('dpPiOverlap')  !Vr(7) = (d p pi) INVERSE
                           Vr_ind = 7
                           if (correct_order) then
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:), TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:), invers=.true.)   ! below
                           else
                              call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                        TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:), invers=.true.)   ! below
                           endif
                           if (.not.read_well) goto 2012 ! exit
                           
                        case ('ddSigmaOverlap')  !Vr(8) = (d d sigma)
                           Vr_ind = 8
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case ('ddPiOverlap')  !Vr(9) = (d d pi)
                           Vr_ind = 9
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case ('ddDeltaOverlap')  !Vr(10) = (d d delta)
                           Vr_ind = 10
                           call  read_line_with_BOP_params(FN_BOP, count_lines, read_well, error_message, &
                                            TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:), TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:))   ! below
                           if (.not.read_well) goto 2012 ! exit
                           if (KOA2 > KOA1) then
                              TB_Hamil(KOA2, KOA1)%S_ci(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ci(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_li(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_li(Vr_ind,:)
                              TB_Hamil(KOA2, KOA1)%S_ni(Vr_ind,:) = TB_Hamil(KOA1, KOA2)%S_ni(Vr_ind,:)
                           endif
                           
                        case (m_bond_txt)  ! new block of elements
!                            print*, 'Read everything, exit the cycle'
                           exit MDL    ! exit the cycle
                        endselect ! (trim(adjustl(temp_txt)))
                        
                     enddo
                  else
                     temp_txt = ''  ! to continue searching
                  endif ! ( (trim(adjustl(temp_txt3)) == trim(adjustl(Atom_1))) .and. (trim(adjustl(temp_txt4)) == trim(adjustl(Atom_2)))  )
               endif ! (trim(adjustl(temp_txt)) == trim(adjustl(bond_txt)))
            enddo BND
         endif  ! ( (trim(adjustl(temp_txt3)) == trim(adjustl(m_BOP_model))) )
      endif !   (trim(adjustl(temp_txt)) == trim(adjustl(m_model_txt)))
   enddo MDL
   
!    ! Test BOP parameteres reading:
!     print*,  'Test BOP parameteres reading: ', KOA1, KOA2,  Atom_1,' ', Atom_2
!     print*,  'H:', TB_Hamil(KOA1, KOA2)%H_ci, TB_Hamil(KOA1, KOA2)%H_li, TB_Hamil(KOA1, KOA2)%H_ni
!     print*,  'E 1-2:', TB_Hamil(KOA1, KOA2)%E_ci, 'li:', TB_Hamil(KOA1, KOA2)%E_li, 'ni:', TB_Hamil(KOA1, KOA2)%E_ni
!     print*,  'E 2-1:', TB_Hamil(KOA2, KOA1)%E_ci, 'li:', TB_Hamil(KOA2, KOA1)%E_li, 'ni:', TB_Hamil(KOA2, KOA1)%E_ni
!     print*,  'S:', TB_Hamil(KOA1, KOA2)%S_ci, TB_Hamil(KOA1, KOA2)%S_li, TB_Hamil(KOA1, KOA2)%S_ni
   
2012    continue
!    pause 'TEST read_BOP_parameters'
end subroutine read_BOP_parameters_OLD





END MODULE Dealing_with_BOP
