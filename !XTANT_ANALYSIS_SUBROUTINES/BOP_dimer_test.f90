PROGRAM BOP_dimer_test
! Test of the tight binding bond parameterization from
! [1] https://arxiv.org/abs/1909.04561
! [2] https://doi.org/10.1016/j.cpc.2018.08.013 
!<===========================================
! Compilation under Windows:
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp  /Qmkl=parallel /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs BOP_dimer_test.f90 -o BOP_dimer_test.exe /link /stack:9999999999 
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays  /Qmkl=parallel BOP_dimer_test.f90 -o BOP_dimer_test.exe /link /stack:9999999999 
!<===========================================
USE IFLPORT
USE OMP_LIB

implicit none

! Parameters of BOP functions (test case of silicon dimer, copied from the file "models.bx"):
! Orthagonal parameters:
real(8), dimension(3), parameter :: & ! maximum 3 exp terms in the current data
Hamilt_ssSigma_ci = (/ -1.2700421114e+02, 1.6146546857e+03, 0.0e0 /), &
Hamilt_ssSigma_li = (/ 1.6954484550e+00, 4.3002825501e+00, 0.0e0 /), &
Hamilt_ssSigma_ni = (/ 1.0000000000e+00, 8.2747683124e-01, 0.0e0 /), &
Hamilt_spSigma_ci = (/ 5.5331123426e+01, -2.1748480288e+05, 0.0e0 /), &
Hamilt_spSigma_li = (/ 1.1638760790e+00, 9.2045509104e+00, 0.0e0 /), &
Hamilt_spSigma_ni = (/ 1.0000000000e+00, 4.2689518014e-01, 0.0e0 /), &
Hamilt_psSigma_ci = (/ -5.5331123426e+01, 2.1748480288e+05, 0.0e0 /), &
Hamilt_psSigma_li = (/ 1.1638760790e+00, 9.2045509104e+00, 0.0e0 /), &
Hamilt_psSigma_ni = (/ 1.0000000000e+00, 4.2689518014e-01, 0.0e0 /), &
Hamilt_ppSigma_ci = (/ 4.9199757780e+01, -5.0943616881e+01, -9.1169375549e+00 /), &
Hamilt_ppSigma_li = (/ 1.0611981009e+00, 1.4908570274e+00, 3.7747307847e-01 /), &
Hamilt_ppSigma_ni = (/ 1.0000000000e+00, 1.0836132492e+00, 3.5626986052e+00 /), &
Hamilt_ppPi_ci = (/ -2.3827392938e+01, 5.3036278288e-01, 0.0e0 /), &
Hamilt_ppPi_li = (/ 1.2511004637e+00, 9.4503161756e-02, 0.0e0 /), &
Hamilt_ppPi_ni = (/ 1.0000000000e+00, 6.4675261893e+00, 0.0e0 /), &
Onsite_s_ci = (/ -10.5973551802, 11.7459035206, -4.8747408590 /), &
Onsite_s_li = (/ 0.0000000000, 0.8546945228, 0.3548522019 /), &
Onsite_s_ni = (/ 1.0000000000, 1.4254599330, 4.3920307652 /), &
Onsite_p0_ci = (/ -3.8993235302, 0.4548483102, -4.3488116614 /), &
Onsite_p0_li = (/ 0.0000000000, 0.0539715288, 0.2881897180 /), &
Onsite_p0_ni = (/ 1.0000000000, 2.8746370880, 3.5586989375 /), &
Onsite_p1_ci = (/ -3.8987155243, -5.3807428689, -0.5607656054 /), &
Onsite_p1_li = (/ 0.0000000000, 1.4346988314, 0.2414858703 /), &
Onsite_p1_ni = (/ 1.0000000000, 1.1308996985, 4.4733434982/)

! Parameters to user:
real(8), parameter :: Dist_start=1.5e0, Dist_end=8.0e0  ! [A] distances to create plots
real(8), dimension(3) :: R_at1, R_at2   ! [A]
real(8) :: Dist ! [A] distance between them
real(8) :: angles   ! cosine directions between the atoms
real(8) :: Dist_step    ! [A] step to print out data
real(8) :: Hopping_ssSigma, Hopping_spSigma, Hopping_psSigma, Hopping_ppSigma, Hopping_ppPi ! Radial hopping inegrals [eV]
real(8) :: Onsite_s, Onsite_p0, Onsite_p1   ! Radial onsite functions [A]

! Variables to be used:
integer :: i, i_step
integer, parameter :: FN1=1111, FN2=1112, FN3=1113

real(8) :: l, m, n  ! directional cosines
real(8), dimension(8,8) :: Total_H  ! total hamiltonian
real(8), dimension(4,4) :: Hij_onsite, Hij_offdiag, Hij_offdiag2  ! block matrices to construct total Hamiltonian
real(8), dimension(4,4) :: Hij_onsite_rotated1, Hij_onsite_rotated2    ! block matrices to construct total Hamiltonian after rotation
real(8), dimension(4,4) :: Hij_rotated1, Hij_rotated2    ! block matrices to construct total Hamiltonian after rotation
real(8), dimension(4,4) :: Rot_matrix1, Rot_matrix2    ! rotation matrix for sp3 basis for first and second atom
real(8), dimension(8) :: Eigenval_diag, Eigenval_SK

! Output file names:
character(50), parameter :: File_out_radial='OUT_radial_functions.dat', &   ! Radial functions, Eqs. (11-12) [1]
File_out_eigenvalues_diag='OUT_eigenvalues_direct_diag.dat', &  ! Direct diagonalization with rotation matrices [2]
File_out_eigenvalues_SK='OUT_eigenvalues_SK.dat'    ! Slater-Koster method

character(200) :: Error_descript

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Error_descript = '' ! to start with, no error

! Create/Open output files:
open(UNIT=FN1, FILE = trim(adjustl(File_out_radial)) )
open(UNIT=FN2, FILE = trim(adjustl(File_out_eigenvalues_diag)) )
open(UNIT=FN3, FILE = trim(adjustl(File_out_eigenvalues_SK)) )
! Set description of the parameters:
write(FN1, '(a)') 'Distance Onsite_s    Onsite_p0   Onsite_p1   Hopping_ssSigma Hopping_spSigma Hopping_psSigma Hopping_ppSigma Hopping_ppPi'

! For the interatomic distances in the required interval:
Dist_step = (Dist_end - Dist_start)/100.0e0 ! [A] step
Dist = Dist_start   ! starting interatomic distance [A]
do i = 1, 101   ! to increase interatomic distance

   !---------------------------------------
   ! 0) Set atomic coordinates:
   R_at1 = 0.0d0    ! first atom is at the origin
   R_at2 = (/ Dist, 0.0d0, 0.0d0 /)    ! second atom is shifted along X
   ! And directional cosines:
   l = R_at2(1)/Dist     ! cos x
   m = R_at2(2)/Dist   ! cos y
   n = R_at2(3)/Dist    ! cos z

   !---------------------------------------
   ! 1) Construct the radial functions for the matrix hamiltonian:
   
   ! Create hopping radial functions:
   Hopping_ssSigma = BOP_radial_function(Hamilt_ssSigma_ci, Hamilt_ssSigma_li, Hamilt_ssSigma_ni, Dist)   ! see function below
   Hopping_spSigma = BOP_radial_function(Hamilt_spSigma_ci, Hamilt_spSigma_li, Hamilt_spSigma_ni, Dist)   ! see function below
   Hopping_psSigma = BOP_radial_function(Hamilt_psSigma_ci, Hamilt_psSigma_li, Hamilt_psSigma_ni, Dist)   ! see function below
   Hopping_ppSigma = BOP_radial_function(Hamilt_ppSigma_ci, Hamilt_ppSigma_li, Hamilt_ppSigma_ni, Dist)   ! see function below
   Hopping_ppPi  = BOP_radial_function(Hamilt_ppPi_ci, Hamilt_ppPi_li, Hamilt_ppPi_ni, Dist)   ! see function below

   ! Create on-site radial functions:
   Onsite_s = BOP_radial_function(Onsite_s_ci, Onsite_s_li, Onsite_s_ni, Dist)   ! see function below
   Onsite_p0 = BOP_radial_function(Onsite_p0_ci, Onsite_p0_li, Onsite_p0_ni, Dist)   ! see function below
   Onsite_p1 = BOP_radial_function(Onsite_p1_ci, Onsite_p1_li, Onsite_p1_ni, Dist)   ! see function below
   
   ! Printout parameters of the radial functions (to make sure, we construct radial functions correctly):
   write(FN1, '(es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16,es25.16)') Dist, Onsite_s, Onsite_p0, Onsite_p1, &
         Hopping_ssSigma, Hopping_spSigma, Hopping_psSigma, Hopping_ppSigma, Hopping_ppPi   ! tested, all agree with Fig. 2 from [1]
   
   !---------------------------------------
   ! 2) Construct the hamiltonian matrix, including angular part:

   ! 2a) Via rotation matrix according to [2]:
   ! Create rotation matrix:
   call set_rotation_matrix(l, m, n, Rot_matrix1)  ! subroutine below
   call set_rotation_matrix(-l, -m, -n, Rot_matrix2)  ! subroutine below

   ! Create Hamiltonian matrix:
   Hij_onsite = RESHAPE( (/ Onsite_s, 0.0d0, 0.0d0, 0.0d0,  &
                                            0.0d0, Onsite_p0, 0.0d0, 0.0d0,  &
                                            0.0d0, 0.0d0, Onsite_p1, 0.0d0,  &
                                            0.0d0, 0.0d0, 0.0d0,  Onsite_p1  /), (/4,4/) ) ! Unrotated onsite elements

   Hij_offdiag = RESHAPE( (/ Hopping_ssSigma, Hopping_spSigma, 0.0d0, 0.0d0,  &
                        Hopping_psSigma, Hopping_ppSigma, 0.0d0, 0.0d0,  &
                        0.0d0, 0.0d0, Hopping_ppPi, 0.0d0,  &
                        0.0d0, 0.0d0, 0.0d0,  Hopping_ppPi  /), (/4,4/) )   ! Unrotated off-diagonal (hopping integrals) elements

   ! Rotate them via rotation matrices:
   Hij_onsite_rotated1 = MATMUL (TRANSPOSE(Rot_matrix1) , MATMUL (Hij_onsite , (Rot_matrix1)))   ! completed rotation (atom 1)
   Hij_onsite_rotated2 = MATMUL (TRANSPOSE(Rot_matrix2) , MATMUL (Hij_onsite , (Rot_matrix2)))   ! completed rotation (atom 2)
   Hij_rotated1 = MATMUL (TRANSPOSE(Rot_matrix1) , MATMUL (Hij_offdiag , (Rot_matrix1)))   ! completed rotation (atom 1)
   Hij_rotated2 = MATMUL (TRANSPOSE(Rot_matrix2) , MATMUL (Hij_offdiag , (Rot_matrix2)))   ! completed rotation (atom 2)
   
   ! Construct the total (rotated) hamiltonian matrix:
   Total_H(1:4,1:4) = Hij_onsite_rotated1
   Total_H(1:4,5:8) = Hij_rotated1
   Total_H(5:8,1:4) = Hij_rotated2
   Total_H(5:8,5:8) = Hij_onsite_rotated2

   ! Diagonalize the Hamiltonian to get eigenvalues:
   call r_diagonalize(Total_H, Eigenval_diag, Error_descript, check_M=.true.)  ! subroutine below
   if (LEN(trim(adjustl(Error_descript))) > 0) then
      Error_descript = 'Subroutine r_diagonalize: '//trim(adjustl(Error_descript))
      print*, trim(adjustl(Error_descript))
   endif
   ! Printout parameters of the radial functions (to make sure, we construct radial functions correctly):
   write(FN2,'(f16.6,f16.6,f16.6,f16.6,f16.6,f16.6,f16.6,f16.6,f16.6)') Dist, Eigenval_diag(1:8)
   
   !---------------------------------------
   ! 2b) Via Slater-Koster approximation:
   Hij_onsite = RESHAPE( (/ Onsite_s, 0.0d0, 0.0d0, 0.0d0,  &
                        0.0d0, Onsite_p0 * l*l + (1.0d0-l*l)*Onsite_p1, 0.0d0, 0.0d0,  &
                        0.0d0, 0.0d0, Onsite_p0 * m*m + (1.0d0-m*m)*Onsite_p1, 0.0d0,  &
                        0.0d0, 0.0d0, 0.0d0,  Onsite_p0 * n*n + (1.0d0-n*n)*Onsite_p1 /), (/4,4/) )  ! Onsite elements

   Hij_offdiag = RESHAPE( (/ Hopping_ssSigma, Hopping_spSigma*l,  Hopping_spSigma*m,  Hopping_spSigma*n,  & ! Off-diagonal (hopping integrals) elements
                    Hopping_psSigma*l, Hopping_ppSigma*l*l + (1.0d0-l*l)*Hopping_ppPi, l*m*(Hopping_ppSigma-Hopping_ppPi),  l*n*(Hopping_ppSigma-Hopping_ppPi),  &
                    Hopping_psSigma*m, l*m*(Hopping_ppSigma-Hopping_ppPi),  Hopping_ppSigma*m*m + (1.0d0-m*m)*Hopping_ppPi,  m*n*(Hopping_ppSigma-Hopping_ppPi),  &
                    Hopping_psSigma*n,  l*n*(Hopping_ppSigma-Hopping_ppPi), m*n*(Hopping_ppSigma-Hopping_ppPi),  Hopping_ppSigma*n*n + (1.0d0-n*n)*Hopping_ppPi /), (/4,4/) )

   ! Construct the total hamiltonian matrix:
   Total_H(1:4,1:4) = Hij_onsite
   Total_H(1:4,5:8) = Hij_offdiag
   Total_H(5:8,1:4) = TRANSPOSE(Hij_offdiag)
   Total_H(5:8,5:8) = Hij_onsite
    ! Diagonalize the Hamiltonian to get eigenvalues:
   call r_diagonalize(Total_H, Eigenval_diag, Error_descript, check_M=.true.)  ! subroutine below
   if (LEN(trim(adjustl(Error_descript))) > 0) then
      Error_descript = 'Subroutine r_diagonalize: '//trim(adjustl(Error_descript))
      print*, trim(adjustl(Error_descript))
   endif
   ! Printout parameters of the radial functions (to make sure, we construct radial functions correctly):
   write(FN3,'(f16.6,f16.6,f16.6,f16.6,f16.6,f16.6,f16.6,f16.6,f16.6)') Dist, Eigenval_diag(1:8)
   
   !---------------------------------------
   Dist = Dist + Dist_step  ! increase the interatomic distance
enddo


! Close output files:
close(FN1)
close(FN2)
close(FN3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains  ! Subroutines and functions



! wrapper for LAPACK diagonalization of real symmetric matrix:
subroutine r_diagonalize(M, Ev, Error_descript, print_Ei, check_M, use_DSYEV) 
   real(8), dimension(:,:), intent(inout) :: M	! matrix
   real(8), dimension(:), intent(out) :: Ev	! eigenvalues
   character(*), intent(inout) :: Error_descript	! error description
   logical, intent(in), optional :: print_Ei ! print out eigenvalues or not
   logical, intent(in), optional :: check_M  ! check the diagonalized matrix
   logical, intent(in), optional :: use_DSYEV	! do not use divide-and-conquare algorithm
   integer :: LIWORK, LWORK, N, INFO
   real(8), dimension(:), allocatable :: LAPWORK
   INTEGER, dimension(:), allocatable :: IWORK
   !real(8), dimension(size(M,1),size(M,2)) :: M_save	! matrix
   real(8), dimension(:,:), allocatable :: M_save	! matrix
   integer :: i_countin, FN, i, j

   N = size(M,1)
   allocate(M_save(N,N))
   LIWORK = 30 + 5*N
   LWORK = 10 + 6*N + 2*N*N
   
   allocate(LAPWORK(LWORK))
   allocate(IWORK(LIWORK))

   !$OMP WORKSHARE
   M_save = M ! save matrix before diagonalization just in case
   !$OMP END WORKSHARE

   if (.not.present(use_DSYEV)) then ! diagonalize real symmetric matrix
      call dsyevd('V','U', N, M, N, Ev, LAPWORK, LWORK, IWORK, LIWORK, INFO) ! LAPACK
   else
      call dsyev('V', 'U', N, M, N, Ev, LAPWORK, LWORK, INFO)   ! LAPACK
   endif
   
   if (INFO .NE. 0) then ! if divide-n-conquare LAPACK diagonalization procidure failed, try regular one:
      !$OMP WORKSHARE
      M = M_save
      !$OMP END WORKSHARE
      call dsyev('V', 'U', N, M, N, Ev, LAPWORK, LWORK, INFO)
      if (INFO .NE. 0) then ! if LAPACK diagonalization procidure failed:
         Error_descript = 'Module Algebra_tools: real matrix diagonalization failed!'
         print*, trim(adjustl(Error_descript))
      endif
   endif
   
   if (present(check_M)) then
      call check_Ha(M_save, M, Ev)   ! to make sure diagonalization went well
   endif

   if (present(print_Ei)) then
      if (print_Ei) then
         do i_countin = 1, size(Ev)
            write(*,'(a,i4,a,f12.5)') 'Eigenvalue #', i_countin, ' is ', Ev(i_countin)
         enddo
      endif
   endif
end subroutine r_diagonalize


! check of eigenvalues and eigenvectors produced by LAPACK diagonalization are really eigenvalues of the hamiltonian:
subroutine check_Ha(Mat, Eigenvec, Eigenval) 
   real(8), dimension(:,:), intent(in) :: Mat	! matrix
   real(8), dimension(:,:), intent(in) :: Eigenvec	! eigenvectors of this matrix
   real(8), dimension(:), intent(in) :: Eigenval	! eigenvalues of this matrix
   real(8), dimension(size(Mat,1), size(Mat,2)) :: M_temp
   real(8) Ev, Evec(size(Eigenval))
   integer i, j, k, M

   M = size(Eigenval)
   Ev = 0.0d0
   do i = 1,M
      do k = 1,M
          Evec(k) = SUM(Mat(:,k)*Eigenvec(:,i))
      enddo
      Ev = SUM(Evec(:)*Eigenvec(:,i))
      if (ABS(Ev - Eigenval(i)) .GT. 1.0d-10*min(ABS(Ev),ABS(Eigenval(i)))) then ! diagonalization went wrong:
         print*, 'Diagonalization error: Ev:', i, Ev, Eigenval(i)
      endif
   enddo
end subroutine check_Ha


! Construct rotation matrix written in terms of direction cosines:
subroutine set_rotation_matrix(l, m, n, Rot_mat)
   real(8), intent(in) :: l, m, n    ! directional cosines to contruct rotation matrix
   real(8), dimension(4,4), intent(out) :: Rot_mat  ! sp3 rotation matrix
   real(8) sqrt_l_m
   sqrt_l_m = sqrt(l*l+m*m)
   Rot_mat(1,1) = 1.0d0
   Rot_mat(1,2) = 0.0d0
   Rot_mat(1,3) = 0.0d0
   Rot_mat(1,4) = 0.0d0
   Rot_mat(2,1) = 0.0d0
   Rot_mat(2,2) = n
   Rot_mat(2,3) = -m
   Rot_mat(2,4) = l
   Rot_mat(3,1) = 0.0d0
   Rot_mat(3,2) = sqrt(1.0d0 - n*n)
   if (ABS(sqrt_l_m) > 1.0d-9) then ! not along Z
      Rot_mat(3,3) = m*n/sqrt_l_m
      Rot_mat(3,4) = -l*n/sqrt_l_m
      Rot_mat(4,3) = l/sqrt_l_m
      Rot_mat(4,4) = m/sqrt_l_m
   else ! along Z
      Rot_mat(3,3) = n
      Rot_mat(3,4) = -n
      Rot_mat(4,3) = 1.0d0
      Rot_mat(4,4) = 1.0d0
   endif
   Rot_mat(4,1) = 0.0d0
   Rot_mat(4,2) = 0.0d0
end subroutine set_rotation_matrix


! Radial function as a sum of exponentials, Eqs.(11) and (12) [1]:
function BOP_radial_function(ci, lambda, ni, R) result(Rfunc)
   real(8) Rfunc
   real(8), dimension(3), intent(in) :: ci, lambda, ni
   real(8), intent(in) :: R
   real(8) :: eps, arg
   integer :: i
   eps = 1.0d-10     ! precision
   Rfunc = 0.0d0    ! to start summation
   do i = 1, 3    ! sum up BOP fitting functions (maximum 3 terms in the current data for Si)
      if ( abs(ci(i)) < eps ) exit    ! exclude zeros
      if ( (abs(R) < eps) .or. (abs(lambda(i)) < eps ) ) then ! no exponential needed
         Rfunc = Rfunc + ci(i)
      else
         arg = lambda(i) * R**ni(i)
         Rfunc = Rfunc + ci(i) * exp(-arg)
      endif
   enddo
end function BOP_radial_function


END PROGRAM BOP_dimer_test
