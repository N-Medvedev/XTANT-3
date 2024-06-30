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
! This module contains subroutines for dealing with the xTB input files (UNFINISHED!)
! Available from https://github.com/grimme-lab/xtb
! PArameters pAlpha1 are copied from xTB file "slater.f90"

MODULE Dealing_with_xTB
use Universal_Constants   ! let it use universal constants
use Objects   ! since it uses derived types, it must know about them from module 'Objects'
use Dealing_with_files, only : read_file, close_file
use BS_Basis_sets, only : set_xTB_AO

implicit none
PRIVATE

! Modular parameters:
character(10), parameter :: m_xTB_directory = 'xTB'
character(20), parameter :: m_xTB_file = 'param_gfn0-xtb.txt'

public :: m_xTB_directory, read_xTB_parameters, identify_basis_size_xTB, identify_AOs_xTB


 contains

subroutine read_xTB_parameters(xTB_Folder, i,j, TB_Hamil, TB_Rep, numpar, matter, Error_descript, INFO)
   character(*), intent(in) :: xTB_Folder   ! address of the xTB directory
   integer, intent(in) :: i, j  ! numbers of pair of elements for which we read the data
   type(TB_H_xTB), dimension(:,:), intent(inout) ::  TB_Hamil    ! parameters of the Hamiltonian of TB
   type(TB_Rep_xTB), dimension(:,:), intent(inout) :: TB_Rep    ! repulsive potential
   type(Numerics_param), intent(inout) :: numpar 	! all numerical parameters
   type(Solid), intent(in) :: matter	! all material parameters
   character(*), intent(out) :: Error_descript	! error save
   integer, intent(out) :: INFO	! error description
   !------------------------------------------------------
   character(200) :: File_name
   character(10) :: temp_ch, temp_ch2
   integer count_lines, Reason, i_cur, Z_cur, FN_GFN0, ToA, N_basis_siz
   logical file_exist, file_opened, read_well
   INFO = 0
   ! In xTB, we only need the parameters for element, not pair-wise
   if (i /= j) goto 3426   ! exit subroutine if nothing else to do

   ! xTB parameterization does not exist for Z>86, check it:
   if (matter%Atoms(i)%Z > 86) then
      write(temp_ch, '(i3)') matter%Atoms(i)%Z
      Error_descript = 'Element '//trim(adjustl(temp_ch))//' is not supported in xTB, the program terminates'
      INFO = 5
      goto 3426
   endif

   ! If it is an element, continue
   count_lines = 1
   
   File_name = trim(adjustl(xTB_Folder))//trim(adjustl(m_xTB_file))     ! file with the parameters

   ! Check if such xTB parameterization exists:
   inquire(file=trim(adjustl(File_name)),exist=file_exist)
   if (.not.file_exist) then
      Error_descript = 'File '//trim(adjustl(File_name))//' not found, the program terminates'
      INFO = 1
      goto 3426
   endif
   FN_GFN0=111
   open(UNIT=FN_GFN0, FILE = trim(adjustl(File_name)), status = 'old', action='read')
   inquire(file=trim(adjustl(File_name)),opened=file_opened)
   if (.not.file_opened) then
      Error_descript = 'File '//trim(adjustl(File_name))//' could not be opened, the program terminates'
      INFO = 2
      goto 3426
   endif


   ! Read parameters to find the line we need:
   i_cur = 0
   RF:do
      i_cur = i_cur + 1
      read(FN_GFN0,*,IOSTAT=Reason) temp_ch
      call read_file(Reason, count_lines, read_well)
      if (.not. read_well) then
         write(Error_descript,'(a,i6,a)') 'Could not read line ', count_lines, ' in the file '//trim(adjustl(File_name))
         INFO = 3
         goto 3426
      endif

      if (trim(adjustl(temp_ch(1:3))) == '$Z=') then    ! found the line with descriptor of new element
         backspace(FN_GFN0)
         read(FN_GFN0,'(3X,i2)',IOSTAT=Reason) Z_cur   ! read the element number
         if (.not. read_well) then
            write(Error_descript,'(a,i6,a)') 'Could not read line ', count_lines, ' in the file '//trim(adjustl(File_name))
            INFO = 3
            goto 3426
         endif

         if (Z_cur == matter%Atoms(i)%Z) then   ! correct element found
            read(FN_GFN0,*,IOSTAT=Reason) temp_ch2  ! read orbitals configuration
            call read_file(Reason, count_lines, read_well)
            if (.not. read_well) then
               write(Error_descript,'(a,i6,a)') 'Could not read line ', count_lines, ' in the file '//trim(adjustl(File_name))
               INFO = 3
               goto 3426
            endif
            TB_Hamil(i,j)%AO_names = trim(adjustl(temp_ch2(4:)))    ! save orbitals descriptor
            exit RF ! found what we needed

         endif ! Z_cur == matter%Atoms(i)%Z
      endif ! (trim(adjustl(temp_ch)) == '$Z=')
   enddo RF



3426 continue
   ! Close files that have been read through:
   call close_file('close', FN=FN_GFN0) ! module "Dealing_with_files"
end subroutine read_xTB_parameters



pure subroutine identify_basis_size_xTB(TB_Hamil, Nsiz)
   integer, intent(out) :: Nsiz  ! index: Cartesian: 0=s; 1=ss*; 2=sp3; 3=sp3s*; 4=sp3d6; 5=sp3d6s*
   type(TB_H_xTB), dimension(:,:), intent(in) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !----------------------------
   integer :: Nkoa, i, lngth, j
   logical :: S_present, SS_present, P_present, D_present
   ! defaults to start:
   Nsiz = 0
   S_present    = .false.
   SS_present   = .false.
   P_present    = .false.
   D_present    = .false.

   Nkoa = size(TB_Hamil,1)  !number of types of atoms
   do i = 1, Nkoa ! check all types of atoms
      lngth = LEN(trim(adjustl(TB_Hamil(i,i)%AO_names)))   ! defines how many orbitals
      ! Check three possible orbitals:
      ! 1) first orbital can be s, p or d:
      if ( trim(adjustl(TB_Hamil(i,i)%AO_names(2:2))) == 's' ) then
         S_present = .true.
         if ( trim(adjustl(TB_Hamil(i,i)%AO_names(4:4))) == 's' ) then
            SS_present = .true.
            cycle  ! it is hydrogen, continue with other elements
         endif
      endif
      if ( trim(adjustl(TB_Hamil(i,i)%AO_names(2:2))) == 'p' ) P_present = .true.
      if ( trim(adjustl(TB_Hamil(i,i)%AO_names(2:2))) == 'd' ) D_present = .true.

      ! 2) second orbital can be s, p or d:
      if ( trim(adjustl(TB_Hamil(i,i)%AO_names(4:4))) == 's' ) S_present = .true.
      if ( trim(adjustl(TB_Hamil(i,i)%AO_names(4:4))) == 'p' ) P_present = .true.
      if ( trim(adjustl(TB_Hamil(i,i)%AO_names(4:4))) == 'd' ) D_present = .true.

      ! 3) third orbital can be s, p or d:
      if (lngth > 4) then
         if ( trim(adjustl(TB_Hamil(i,i)%AO_names(6:6))) == 's' ) S_present = .true.
         if ( trim(adjustl(TB_Hamil(i,i)%AO_names(6:6))) == 'p' ) P_present = .true.
         if ( trim(adjustl(TB_Hamil(i,i)%AO_names(6:6))) == 'd' ) D_present = .true.
      endif
   enddo

   ! Now we know which orbitals present, identify the corresponding basis set index:
   if (S_present .and. P_present .and. D_present .and. SS_present) then ! sp3d6s*
      Nsiz = 5   ! sp3d6s*
   elseif (S_present .and. P_present .and. D_present) then ! sp3d6
      Nsiz = 4   ! sp3d6
   elseif (S_present .and. P_present .and. SS_present) then ! sp3s*
      Nsiz = 3   ! sp3s*
   elseif (S_present .and. P_present) then ! sp3
      Nsiz = 2   ! sp3
   elseif (S_present .and. SS_present) then ! ss*
      Nsiz = 1   ! ss*
   elseif (S_present) then ! s
      Nsiz = 0   ! s
   endif
end subroutine identify_basis_size_xTB


pure function identify_Cartesian_BS_block(xTB_BS_ind) result(BS_size)
   integer :: BS_size
   integer, intent(in) :: xTB_BS_ind    ! index of the basis set in xTB
   ! Size of basic overlap matrix block for Cartesian basis set:
   select case (xTB_BS_ind)
   case (0) ! s
      BS_size = 1
   case (1) ! ss*
      BS_size = 2
   case (2) ! sp3
      BS_size = 4
   case (3) ! sp3s*
      BS_size = 5
   case (4) ! sp3d6
      BS_size = 10
   case (5) ! sp3d6s*
      BS_size = 11
   endselect
end function identify_Cartesian_BS_block



subroutine identify_AOs_xTB(TB_Hamil)
   type(TB_H_xTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   !----------------------------
   integer :: Nkoa, i, lngth, j

   Nkoa = size(TB_Hamil,1)  !number of types of atoms
   do i = 1, Nkoa ! check all types of atoms
      lngth = LEN(trim(adjustl(TB_Hamil(i,i)%AO_names)))   ! defines how many orbitals
      if (.not.allocated(TB_Hamil(i,i)%STO)) allocate (TB_Hamil(i,i)%STO(lngth/2))

      ! Three possible types of orbitals:
      ! 1) first orbital can be s, p or d:
      read(TB_Hamil(i,i)%AO_names(1:1),'(i1)') TB_Hamil(i,i)%STO(1)%n   ! principal quantum number
      TB_Hamil(i,i)%STO(1)%MO_type = trim(adjustl(TB_Hamil(i,i)%AO_names(2:2)))  ! type of orbital
      ! Set its angular momenta:
      call set_xTB_AO(TB_Hamil(i,i)%STO(1)%MO_type, TB_Hamil(i,i)%STO(1)) ! module "BS_Basis_sets"

      ! 2) second orbital can be s, p or d:
      read(TB_Hamil(i,i)%AO_names(3:3),'(i1)') TB_Hamil(i,i)%STO(2)%n ! principal quantum number
      TB_Hamil(i,i)%STO(2)%MO_type = trim(adjustl(TB_Hamil(i,i)%AO_names(4:4)))  ! type of orbital
      ! Set its angular momenta:
      call set_xTB_AO(TB_Hamil(i,i)%STO(2)%MO_type, TB_Hamil(i,i)%STO(2)) ! module "BS_Basis_sets"

      ! 3) third orbital can be s, p or d:
      if (lngth > 4) then
         read(TB_Hamil(i,i)%AO_names(5:5),'(i1)') TB_Hamil(i,i)%STO(3)%n
         TB_Hamil(i,i)%STO(3)%MO_type = trim(adjustl(TB_Hamil(i,i)%AO_names(6:6)))
         call set_xTB_AO(TB_Hamil(i,i)%STO(3)%MO_type, TB_Hamil(i,i)%STO(3)) ! module "BS_Basis_sets"
      endif
   enddo
end subroutine identify_AOs_xTB



END MODULE Dealing_with_xTB
