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
! This module contains subroutines for dealing with the Koster-Slater files (format .skf)
! containing numerical sets of radial Hamiltonian and Overlap integrals: http://www.dftb.org/parameters/

MODULE Dealing_with_DFTB
use Universal_Constants   ! let it use universal constants
use Objects   ! since it uses derived types, it must know about them from module 'Objects'
use Dealing_with_files, only : read_file

implicit none
PRIVATE

! Modular parameters:
character(10) :: m_DFTB_directory
character(20) :: m_DFTB_norep_directory

parameter (m_DFTB_directory = 'DFTB')
parameter (m_DFTB_norep_directory = 'DFTB_no_repulsion')

public :: read_skf_file, same_or_different_atom_types, m_DFTB_directory, construct_skf_filename, idnetify_basis_size, &
m_DFTB_norep_directory, read_skf_file_no_rep

 contains


subroutine read_skf_file(FN, TB_Hamil, TB_Rep, ToA, error_message, rep_pot)
   integer, intent(in) :: FN    ! skf file to read from (must be already opened)
   type(TB_H_DFTB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   type(TB_Rep_DFTB), intent(inout) :: TB_Rep   ! repulsive potential
   integer, intent(in) :: ToA    ! types of atoms (0=the same or 1=different)
   character(100), intent(inout) :: error_message
   logical, intent(in), optional :: rep_pot   ! flag for no repulsion file
   !-----------------------
   real(8) :: dr, temp_r, rd_param(20), temp2, temp3, temp4, temp8
   integer :: Reason, count_lines, Ngrid, i, temp_i, N_basis
   logical :: read_well, rep_present
   character(100) :: temp_ch
   error_message = ''   ! no error at the start
   N_basis = 10 ! number of overlap functions for sp3d5
   count_lines = 0
   temp_ch = ''   ! to start with
   read(FN,*,IOSTAT=Reason) temp_ch
   if (temp_ch(1:1) == '@') then    ! it is parameterization with f orbital, currently not suported!
      error_message = 'DFTB parameterization includes f orbital, not supported in XTANT!'
      goto 2012 ! exit
   endif
   backspace(FN)    ! to read it again
   
   read(FN,'(a)',IOSTAT=Reason) temp_ch
   dr = 0.0d0
   Ngrid = 0   ! to start with
   !read(FN,*,IOSTAT=Reason) dr, Ngrid   ! radial grid step [Bohr], size of the grid
   read(temp_ch,*,IOSTAT=Reason) dr, Ngrid   ! radial grid step [Bohr], size of the grid
   if (Reason /= 0) then   ! try to skip the first character, in case there is something wrong with it
      print*, 'Problem in read_skf_file: first line is unreadable: ', trim(adjustl(temp_ch))
      ! Try to skip possible functional characters befor the number:
      read(temp_ch(4:LEN(temp_ch)),*,IOSTAT=Reason) dr, Ngrid
   endif

   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not.read_well) goto 2012 ! exit
   Ngrid = Ngrid - 1
   dr = dr*g_au2A ! [Bohr] -> [A]
   ! Knowing the size, allocate the arrays:
   if (.not.allocated(TB_Hamil%Rr)) allocate(TB_Hamil%Rr(Ngrid))    ! Grid [A]
   if (.not.allocated(TB_Hamil%Vr)) allocate(TB_Hamil%Vr(Ngrid,N_basis))    ! Hamiltonian radial function [eV]
   if (.not.allocated(TB_Hamil%Sr)) allocate(TB_Hamil%Sr(Ngrid,N_basis))    ! Overlap radial function [eV]
   
   ! Read on-site energies:
   if (ToA == 0) then   ! the same element, on-site energies:
      read(FN,*,IOSTAT=Reason) TB_Hamil%Ed, TB_Hamil%Ep, TB_Hamil%Es, temp_r, &
                               TB_Hamil%Ud, TB_Hamil%Up, TB_Hamil%Us  ! [a.u.] on-site energies; Hubbard coefficients
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) goto 2012 ! exit
      TB_Hamil%Ed = TB_Hamil%Ed * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Ep = TB_Hamil%Ep * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Es = TB_Hamil%Es * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Ud = TB_Hamil%Ud * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Up = TB_Hamil%Up * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Us = TB_Hamil%Us * g_au2ev  ! [Hartree] -> [eV]
   else ! no on-site energies for non-like elements
      TB_Hamil%Ed = 0.0d0
      TB_Hamil%Ep = 0.0d0
      TB_Hamil%Es = 0.0d0
      TB_Hamil%Us = 0.0d0
      TB_Hamil%Up = 0.0d0
      TB_Hamil%Ud = 0.0d0
   endif

   ! Read parameters of repulsive potential in the polinomial form:
   read(FN,*,IOSTAT=Reason) temp_r, TB_Rep%c(:), TB_Rep%rcut   ! [a.u.] coefficients for repulsive potential in polinomial form
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not.read_well) goto 2012 ! exit
   ! Convert in the units used in this code:
   TB_Rep%c(:) = TB_Rep%c(:) * g_au2ev  ! [Hartree] -> [eV]
   TB_Rep%rcut = TB_Rep%rcut * g_au2A  ! [Bohr] -> [A]
   temp2 = g_A2au * g_A2au
   TB_Rep%c(1) = TB_Rep%c(1) * temp2   ! [eV/Bohr^2] -> [eV/A^2]
   temp3 = temp2 * g_A2au 
   TB_Rep%c(2) = TB_Rep%c(2) * temp3   ! [eV/Bohr^3] -> [eV/A^3]
   temp4 = temp2 * temp2
   TB_Rep%c(3) = TB_Rep%c(3) * temp4   ! [eV/Bohr^4] -> [eV/A^4]
   TB_Rep%c(4) = TB_Rep%c(4) * temp2 * temp3   ! [eV/Bohr^5] -> [eV/A^5]
   TB_Rep%c(5) = TB_Rep%c(5) * temp3 * temp3   ! [eV/Bohr^6] -> [eV/A^6]
   TB_Rep%c(6) = TB_Rep%c(6) * temp4 * temp3   ! [eV/Bohr^7] -> [eV/A^7]
   temp8 = temp4 * temp4
   TB_Rep%c(7) = TB_Rep%c(7) * temp4 * temp4   ! [eV/Bohr^8] -> [eV/A^8]
   TB_Rep%c(8) = TB_Rep%c(8) * temp8 * g_A2au   ! [eV/Bohr^9] -> [eV/A^9]

   ! Now we know the grid step, can set the radial grid and read the parametres from the file:
   TB_Hamil%Rr(1) = dr  ! [A]
   do i = 1, Ngrid  ! read the data for all grid points
      if (i > 1) TB_Hamil%Rr(i) = TB_Hamil%Rr(i - 1) + dr  ! [A] set grid points
      read(FN,*, IOSTAT=Reason) rd_param(:)  ! read the parameters line by line
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not. read_well) goto 2012 ! exit
      ! Now distribute the read data into corresponding variables:
      call distribute_parameters_DFTB(TB_Hamil, rd_param, i)   ! below
   enddo
   ! Convert in the units used in this code:
   TB_Hamil%Vr = TB_Hamil%Vr * g_au2ev  ! [Hartree] -> [eV]
   
   ! Read the precalculated repulsive potential in Spline form:
   call read_DFTB_spline(FN, count_lines, TB_Rep, error_message)   ! below

2012 continue
end subroutine read_skf_file



subroutine read_skf_file_no_rep(FN, TB_Hamil, TB_Rep, ToA, error_message, rep_pot)
   integer, intent(in) :: FN    ! skf file to read from (must be already opened)
   type(TB_H_DFTB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   type(TB_Rep_DFTB_no), intent(inout) :: TB_Rep   ! repulsive potential
   integer, intent(in) :: ToA    ! types of atoms (0=the same or 1=different)
   character(100), intent(inout) :: error_message
   logical, intent(in), optional :: rep_pot   ! flag for no repulsion file
   !-----------------------
   real(8) :: dr, temp_r, rd_param(20), temp2, temp3, temp4, temp8
   real(8), dimension(8) :: c   ! [eV] c2, . . . , c9 are the polynomial coefficients
   real(8) :: rcut  ! [A] cutoff radius of the repulsive interaction (can be zero, if the repulsive is described by splines)
   integer :: Reason, count_lines, Ngrid, i, temp_i, N_basis
   logical :: read_well, rep_present
   character(100) :: temp_ch
   error_message = ''   ! no error at the start
   N_basis = 10 ! number of overlap functions for sp3d5
   count_lines = 0
   read(FN,*,IOSTAT=Reason) temp_ch
   if (temp_ch(1:1) == '@') then    ! it is parameterization with f orbital, currently not suported!
      error_message = 'DFTB parameterization includes f orbital, not supported in XTANT!'
      goto 2012 ! exit
   endif
   backspace(FN)    ! to read it again

   read(FN,*,IOSTAT=Reason) dr, Ngrid   ! radial grid step [Bohr], size of the grid
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not.read_well) goto 2012 ! exit
   Ngrid = Ngrid - 1
   dr = dr*g_au2A ! [Bohr] -> [A]
   ! Knowing the size, allocate the arrays:
   if (.not.allocated(TB_Hamil%Rr)) allocate(TB_Hamil%Rr(Ngrid))    ! Grid [A]
   if (.not.allocated(TB_Hamil%Vr)) allocate(TB_Hamil%Vr(Ngrid,N_basis))    ! Hamiltonian radial function [eV]
   if (.not.allocated(TB_Hamil%Sr)) allocate(TB_Hamil%Sr(Ngrid,N_basis))    ! Overlap radial function [-]

   ! Read on-site energies:
   if (ToA == 0) then   ! the same element, on-site energies:
      read(FN,*,IOSTAT=Reason) TB_Hamil%Ed, TB_Hamil%Ep, TB_Hamil%Es, temp_r, &
                               TB_Hamil%Ud, TB_Hamil%Up, TB_Hamil%Us  ! [a.u.] on-site energies; Hubbard coefficients
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not.read_well) goto 2012 ! exit
      TB_Hamil%Ed = TB_Hamil%Ed * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Ep = TB_Hamil%Ep * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Es = TB_Hamil%Es * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Ud = TB_Hamil%Ud * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Up = TB_Hamil%Up * g_au2ev  ! [Hartree] -> [eV]
      TB_Hamil%Us = TB_Hamil%Us * g_au2ev  ! [Hartree] -> [eV]
   else ! no on-site energies for non-like elements
      TB_Hamil%Ed = 0.0d0
      TB_Hamil%Ep = 0.0d0
      TB_Hamil%Es = 0.0d0
      TB_Hamil%Us = 0.0d0
      TB_Hamil%Up = 0.0d0
      TB_Hamil%Ud = 0.0d0
   endif

   ! Read parameters of repulsive potential in the polinomial form (zeros in the line present):
   read(FN,*,IOSTAT=Reason) temp_r, c(:), rcut   ! [a.u.] coefficients for repulsive potential in polinomial form
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not.read_well) goto 2012 ! exit

   ! Now we know the grid step, can set the radial grid and read the parametres from the file:
   TB_Hamil%Rr(1) = dr  ! [A]
   do i = 1, Ngrid  ! read the data for all grid points
      if (i > 1) TB_Hamil%Rr(i) = TB_Hamil%Rr(i - 1) + dr  ! [A] set grid points
      read(FN,*, IOSTAT=Reason) rd_param(:)  ! read the parameters line by line
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not. read_well) goto 2012 ! exit
      ! Now distribute the read data into corresponding variables:
      call distribute_parameters_DFTB(TB_Hamil, rd_param, i)   ! below
   enddo
   ! Convert in the units used in this code:
   TB_Hamil%Vr = TB_Hamil%Vr * g_au2ev  ! [Hartree] -> [eV]

2012 continue
end subroutine read_skf_file_no_rep


subroutine read_DFTB_spline(FN, count_lines, TB_Rep, error_message)
   integer, intent(in) :: FN    ! skf file to read from (must be already opened)
   integer, intent(inout) :: count_lines    ! to read from this line on
   type(TB_Rep_DFTB), intent(inout) :: TB_Rep   ! repulsive potential
   character(100), intent(inout) :: error_message
   !--------------------------------
   real(8) :: temp_r, temp2, temp3
   character(6) :: temp_ch
   integer :: i, Reason, N_spline_grid
   logical :: read_well, found_line
   found_line = .false.
   do while (.not.found_line) 
      read(FN,*, IOSTAT=Reason) temp_ch  ! read the parameters line by line
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not. read_well) goto 2013 ! exit
      found_line = (trim(adjustl(temp_ch)) == 'Spline')
   enddo
   
   ! the number of (subsequent) intervals being described by various cubic splines and cutoff the cutoff of the repulsive interaction:
   read(FN,*, IOSTAT=Reason) N_spline_grid, TB_Rep%rcut_spline
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not. read_well) goto 2013 ! exit
   ! Convert into units used in the code:
   TB_Rep%rcut_spline = TB_Rep%rcut_spline * g_au2A  ! [Bohr] -> [A]
   
   ! Now we know how many intervals are there, allocate the arrays for splines:
   if (.not.allocated(TB_Rep%R)) allocate(TB_Rep%R(N_spline_grid))
   if (.not.allocated(TB_Rep%V_rep)) allocate(TB_Rep%V_rep(N_spline_grid,6))
   TB_Rep%R = 0.0d0 ! just to start
   TB_Rep%V_rep = 0.0d0 ! just to start
   
   ! coefficients of the exponential for too short distances
   read(FN,*, IOSTAT=Reason) TB_Rep%a(:)
   call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
   if (.not. read_well) goto 2013 ! exit
   ! Convert into units used in the code:
   TB_Rep%a(1) = TB_Rep%a(1) * g_A2au  ! 1/[Bohr] -> 1/[A]
   TB_Rep%a(3) = TB_Rep%a(3) * g_au2ev  ! [Hartree] -> [eV]

   ! spline parameters:
   do i = 1, N_spline_grid
      if (i < N_spline_grid) then
         !read(FN,*, IOSTAT=Reason) TB_Rep%R(i), temp_r, TB_Rep%V_rep(i,1:4)
         read(FN,*, IOSTAT=Reason) TB_Rep%R(i), temp_r, TB_Rep%V_rep(i,1), TB_Rep%V_rep(i,2), TB_Rep%V_rep(i,3), TB_Rep%V_rep(i,4)
      else
         !read(FN,*, IOSTAT=Reason) TB_Rep%R(i), temp_r, TB_Rep%V_rep(i,1:6)
         read(FN,*, IOSTAT=Reason) TB_Rep%R(i), temp_r, TB_Rep%V_rep(i,1), TB_Rep%V_rep(i,2), TB_Rep%V_rep(i,3), &
                                    TB_Rep%V_rep(i,4), TB_Rep%V_rep(i,5), TB_Rep%V_rep(i,6)
      endif
      call read_file(Reason, count_lines, read_well, error_message)   ! module "Dealing_with_files"
      if (.not. read_well) goto 2013 ! exit
      !print*, i, TB_Rep%R(i), TB_Rep%R(i)*g_au2A, dble(INT(TB_Rep%R(i)*1.0d6))*1.0d-6
   enddo
   
   ! Convert into units used in the code:
   ! TB_Rep%R = TB_Rep%R * g_au2A  ! [Bohr] -> [A] ORIGINAL
   TB_Rep%R = dble(INT(TB_Rep%R(:)*1.0d6))*1.0d-6 * g_au2A  ! [Bohr] -> [A] Rounded up to 6th digit
   TB_Rep%V_rep = TB_Rep%V_rep * g_au2ev  ! [Hartree] -> [eV]
   TB_Rep%V_rep(:,2) = TB_Rep%V_rep(:,2) * g_A2au   ! [eV/Bohr] -> [eV/A]
   temp2 = g_A2au * g_A2au
   TB_Rep%V_rep(:,3) = TB_Rep%V_rep(:,3) * temp2    ! [eV/Bohr^2] -> [eV/A^2]
   temp3 = temp2 * g_A2au
   TB_Rep%V_rep(:,4) = TB_Rep%V_rep(:,4) * temp3    ! [eV/Bohr^3] -> [eV/A^3]
   TB_Rep%V_rep(:,5) = TB_Rep%V_rep(:,5) * temp2 * temp2    ! [eV/Bohr^4] -> [eV/A^4]
   TB_Rep%V_rep(:,6) = TB_Rep%V_rep(:,6) * temp3 * temp3    ! [eV/Bohr^5] -> [eV/A^5]

!    pause 'read_DFTB_spline'
2013 continue
end subroutine read_DFTB_spline



pure subroutine idnetify_basis_size(TB_Hamil, Nsiz)
   integer, intent(out) :: Nsiz  ! index: 0=s, 1=sp3, 2=sp3d5
   type(TB_H_DFTB), dimension(:,:), intent(inout) ::  TB_Hamil ! parameters of the Hamiltonian of TB
   integer :: Nd, Np, Ns
   integer :: i,j 
   Ns = 0   ! to start with
   Np = 1   ! to start with
   Nd = 1   ! to start with
   ! Check if d-band is used for any element:   
   if ( maxval(ABS(TB_Hamil(:,:)%Ed)) < 1.0d-8) then
      Nd = 0
      ! Check if p band is used for any element:
      if ( maxval(ABS(TB_Hamil(:,:)%Ep)) < 1.0d-8) Np = 0
   endif

   ! Set basis set size:
   Nsiz = Ns + Np + Nd
   ! And shift "empty" energy levels away, not to interfer with real calculations
   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,2)
!          print*, i, j, TB_Hamil(i,j)%Ed
         if (ABS(TB_Hamil(i,j)%Ed) < 1.0d-8) TB_Hamil(i,j)%Ed = 80.0d0 
      enddo
   enddo
   do i = 1, size(TB_Hamil,1)
      do j = 1, size(TB_Hamil,2)
         if (ABS(TB_Hamil(i,j)%Ep) < 1.0d-8) TB_Hamil(i,j)%Ep = 100.0d0
      enddo
   enddo
end subroutine idnetify_basis_size



pure subroutine distribute_parameters_DFTB(TB_Hamil, rd_param, ir) ! according to the file "slakoformat.pdf" from www.dftb.org
   type(TB_H_DFTB), intent(inout) :: TB_Hamil    ! Hamiltonian and Overlap data
   real(8), dimension(:), intent(in) :: rd_param    ! parameters read from a file
   integer, intent(in) :: ir    ! index of the radial grid step
   !    1       2       3       4       5       6       7       8       9   10
   ! Hdd0 Hdd1 Hdd2 Hpd0 Hpd1 Hpp0 Hpp1 Hsd0 Hsp0 Hss0 
   ! Sdd0 Sdd1 Sdd2 Spd0 Spd1 Spp0 Spp1 Ssd0 Ssp0 Sss0
   TB_Hamil%Vr(ir,1) = rd_param(10) ! V(1) = (s s sigma)
   TB_Hamil%Vr(ir,2) = rd_param(9) ! V(2) = (s p sigma)
   TB_Hamil%Vr(ir,3) = rd_param(8) ! V(3) = (s d sigma)
   TB_Hamil%Vr(ir,4) = rd_param(6) ! V(4) = (p p sigma)
   TB_Hamil%Vr(ir,5) = rd_param(7) ! V(5) = (p p pi)
   TB_Hamil%Vr(ir,6) = rd_param(4) ! V(6) = (p d sigma)
   TB_Hamil%Vr(ir,7) = rd_param(5) ! V(7) = (p d pi)
   TB_Hamil%Vr(ir,8) = rd_param(1) ! V(8) = (d d sigma)
   TB_Hamil%Vr(ir,9) = rd_param(2) ! V(9) = (d d pi)
   TB_Hamil%Vr(ir,10) = rd_param(3) ! V(10) = (d d delta)
   TB_Hamil%Sr(ir,1) = rd_param(20) ! S(1) = (s s sigma)
   TB_Hamil%Sr(ir,2) = rd_param(19) ! S(2) = (s p sigma)
   TB_Hamil%Sr(ir,3) = rd_param(18) ! S(3) = (s d sigma)
   TB_Hamil%Sr(ir,4) = rd_param(16) ! S(4) = (p p sigma)
   TB_Hamil%Sr(ir,5) = rd_param(17) ! S(5) = (p p pi)
   TB_Hamil%Sr(ir,6) = rd_param(14) ! S(6) = (p d sigma)
   TB_Hamil%Sr(ir,7) = rd_param(15) ! S(7) = (p d pi)
   TB_Hamil%Sr(ir,8) = rd_param(11) ! S(8) = (d d sigma)
   TB_Hamil%Sr(ir,9) = rd_param(12) ! S(9) = (d d pi)
   TB_Hamil%Sr(ir,10) = rd_param(13) ! S(10) = (d d delta)
end subroutine distribute_parameters_DFTB



pure function same_or_different_atom_types(Name1, Name2) result(ToA)
   character(*), intent(in) :: Name1, Name2 ! two names of atoms / chemical elements
   integer :: ToA   ! 0 = same, 1= different
   if ( trim(adjustl(Name1)) == trim(adjustl(Name2)) ) then
      ToA = 0
   else
      ToA = 1
   endif
end function same_or_different_atom_types



pure subroutine construct_skf_filename(Element1, Element2, Filename, add_text)
   character(*), intent(in) :: Element1, Element2   ! elements names for which we require parameterization
   character(*), intent(inout) :: Filename  ! corresponding filename with skd file
   character(*), intent(in), optional :: add_text
   !------------------
   if (present(add_text)) then
      write(Filename,'(a,a,a,a)') trim(adjustl(Element1)), '-', trim(adjustl(Element2))//trim(adjustl(add_text)), '.skf'
   else
      write(Filename,'(a,a,a,a)') trim(adjustl(Element1)), '-', trim(adjustl(Element2)), '.skf'
   endif
end subroutine construct_skf_filename
 

END MODULE Dealing_with_DFTB
