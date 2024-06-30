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
! This module contains all defined objects

MODULE Objects
 implicit none


!==============================================
! Different TB parametrizations will have different functional dependencies,
! meaning it require different data types.
! All of them must be defined here, and later appropriate ones must be chosen.
! Repulsive potential TB parametrization
type :: TB_repulsive ! parent type
   character(25) :: Param ! name of parametrization
end type
! Hamiltonian TB parametrization
type :: TB_Hamiltonian ! parent type
   character(25) :: Param ! name of parametrization
end type
! van der Waals within TB parametrization
type :: TB_vdW ! parent type
   character(25) :: Param ! name of parametrization
   real(8) :: d0_cut, dd_cut  ! [A] [A] cutoff radius and smoothing distance
end type
! Coulomb potential:
type :: TB_Coulomb ! parent type
   character(25) :: Param ! name of parametrization
end type
! Exponential wall within TB parametrization
type :: TB_Exp_wall ! parent type
   character(25) :: Param ! name of parametrization
end type



!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
! Primitive gaussian parameters:
type Primitive_gaussian
   real(8) :: alpha ! alpha {1/r^2} => [1/au^2]
   real(8) :: C     ! coefficients before exponents
   real(8) :: Norm  ! normalization factor
endtype Primitive_gaussian

! Basis set type GTO:
type Basis_set      ! Gaussian basis set parameteres
   integer :: n     ! principal quantum number
   character(1) :: MO_type  ! type of orbital: S, P, D ...
   type(Primitive_gaussian), dimension(:), allocatable :: GTO	! contracted Gaussian Type Orbital, defined above
   integer, dimension(3) :: AM  ! Cartesian angular momentum quantum number in X,Y,Z directions; Spherical: n, l, m
   real(8) :: zeta  ! zeta for STO adjustment: [Szabo, Ostlund, "Modern Quantum Chemistry", p. 185]
   integer :: index ! index different exponents (not to repeat the same calculations for different orbitals with the same exponents)
endtype Basis_set


! Basis set type STO:
type AngMom ! angular momenta for STO
   integer, dimension(3) :: AM  ! Cartesian angular momentum quantum number in X,Y,Z directions: l, m, n
endtype AngMom


type Basis_set_STO      ! Slatter-type basis set parameteres
   integer :: n     ! principal quantum number
   character(1) :: MO_type  ! type of orbital: S, P, D ...
   type(Primitive_gaussian), dimension(:), allocatable :: GTO	! contracted Gaussian Type Orbital, defined above
   type(AngMom), dimension(:), allocatable :: AM_set  ! Set of angular momenta defined by AO type (1 for s, 3 for p, 6 for d)
   real(8) :: zeta  ! zeta for STO adjustment: [Szabo, Ostlund, "Modern Quantum Chemistry", p. 185]
endtype Basis_set_STO


!+111111111111111111111111111111111111111111111
! 1) TB parameters by Pettifor (or similar):
! a) repulsive part:
type, EXTENDS (TB_repulsive) :: TB_Rep_Pettifor	! repulsive potential coefficients:
   real(8) :: E0_TB   ! [eV] constant shift of the Hamiltonain, H. Jeschke PhD thesis, p. 142
   real(8), dimension(5) :: a0	! coefficients of the polinomial form for repulsive potential
   real(8) :: phi0, d0, m, mc, dc, d1, dm
   real(8), dimension(4) :: c0
end type TB_Rep_Pettifor

type, EXTENDS (TB_Rep_Pettifor) :: TB_Rep_Fu	! repulsive potential coefficients:
   ! [C.-C. Fu and M. Weissmann, Phys. Rev. B 60 (1999) 2762]
   real(8) :: C_a	! coefficients for anions and cations
end type TB_Rep_Fu

type, EXTENDS (TB_repulsive) :: TB_Rep_Molteni	! repulsive potential coefficients:
   ! [C. Molteni et al. J. Phys.: Condens. Matter 6 (1994) 5243]
   integer :: NP, m ! Np = number of potential: 1=exp (Table 3), 2=rational (Table 4), 3=Allen's [Graves, Allen, PRB 58, 13627 (1998)]
   real(8) :: phi1, phi2, r0, alpha
   real(8) :: rcut, d
   real(8) :: a, b, c	! Allen's potential
end type TB_Rep_Molteni

type, EXTENDS (TB_repulsive) :: TB_Rep_NRL	! repulsive potential coefficients:
   ! D.A. Papaconstantopoulos and M.J. Mehl, J. Phys.: Condens. Matter 15 (2003) R413–R440
   ! No repulsive part in this parameterization
end type TB_Rep_NRL

type, EXTENDS (TB_repulsive) :: TB_Rep_BOP	! repulsive potential coefficients:
   ! https://arxiv.org/pdf/1909.04561.pdf
   ! Since BOP model does not provide repulsive terms, this potential
   ! is reconstructed from ZBL repulsive potential:
   ! https://en.wikipedia.org/wiki/Stopping_power_(particle_radiation)#Repulsive_interatomic_potentials
   real(8), dimension(:), allocatable :: R  ! [A] distance
   real(8), dimension(:), allocatable :: V_rep  ! [eV] parameterized repulsive potential
   ! unfinished, do not use!
end type TB_Rep_BOP


type, EXTENDS (TB_repulsive) :: TB_Rep_DFTB	! repulsive potential coefficients:
   ! www.dftb.org
   character(200) :: param_name  ! name of parameterization used
   integer :: ToP   ! type of parameterization: 0=polinomial, 1=spline
   ! Polinomial coefficients:
   real(8), dimension(8) :: c   ! [eV] c2, . . . , c9 are the polynomial coefficients
   real(8) :: rcut  ! [A] cutoff radius of the repulsive interaction (can be zero, if the repulsive is described by splines)
   ! Spline:
   real(8) :: rcut_spline   ! [A] cutoff radius of the repulsive interaction in case of spline
   real(8), dimension(3) :: a   ! parameters for spline: exponential at short distances
   real(8), dimension(:), allocatable :: R  ! [A] distance
   real(8), dimension(:,:), allocatable :: V_rep  ! [eV] parameterized repulsive potential for spline
end type TB_Rep_DFTB


type, EXTENDS (TB_repulsive) :: TB_Rep_DFTB_no  ! repulsive potential coefficients:
   ! Testing construction of the repulsive potential:
   character(200) :: param_name  ! name of parameterization used
!    integer :: ToP   ! type of parameterization: 0=polinomial, 1=spline
!    ! In case there is no repulsive terms, this potential is reconstructed from ZBL repulsive potential:
!    ! https://en.wikipedia.org/wiki/Stopping_power_(particle_radiation)#Repulsive_interatomic_potentials
!    real(8), dimension(:), allocatable :: R  ! [A] distance
!    real(8), dimension(:), allocatable :: V_rep  ! [eV] parameterized repulsive potential
end type TB_Rep_DFTB_no


type, EXTENDS (TB_repulsive) :: TB_Rep_3TB	! repulsive potential coefficients:
   ! https://arxiv.org/pdf/2112.11585.pdf
   ! No repulsive part in this parameterization
end type TB_Rep_3TB


type, EXTENDS (TB_repulsive) :: TB_Rep_xTB  ! Repulsive part of the xTB Hamiltonian
   ! https://github.com/grimme-lab/xtb
   ! unfinished, do not use!
end type TB_Rep_xTB


! b) Hamiltonian:
type, EXTENDS (TB_Hamiltonian) :: TB_H_Pettifor ! hamiltonian coefficients:
   real(8) :: Es, Ep
   real(8), dimension(4) :: nc, rc, c0, c1, c2, c3, V0
   real(8) :: r0, n, r1, rm
end type TB_H_Pettifor

type, EXTENDS (TB_H_Pettifor) :: TB_H_Fu ! hamiltonian coefficients:
   ! [C.-C. Fu and M. Weissmann, Phys. Rev. B 60 (1999) 2762]
   real(8) :: C_a	! coefficients for anions and cations
end type TB_H_Fu

type, EXTENDS (TB_Hamiltonian) :: TB_H_Molteni ! hamiltonian coefficients:
   ! [C. Molteni et al. J. Phys.: Condens. Matter 6 (1994) 5243]
   real(8) :: Es, Ep, Esa
   real(8), dimension(5) :: V0
   real(8) :: nc, rc
   real(8) :: r0, n
   real(8) :: rcut, d
end type TB_H_Molteni

type, EXTENDS (TB_Hamiltonian) :: TB_H_NRL ! hamiltonian coefficients:
   ! [D.A. Papaconstantopoulos and M.J. Mehl, J. Phys.: Condens. Matter 15 (2003) R413–R440]
   integer :: ind_split, ind_overlap 	! which parameterization format
   real(8) :: lambd	! exp. parameter in rho function, Eq.(15)
   real(8) :: Rc, lden	! parameters in F(R), Eq.(16)
   real(8), dimension(4) :: al, bl, cl, dl	! set of parameters in Eq.(17)
   ! Hopping integral parameters:
   real(8), dimension(10) :: ellm, fllm, gllm, hllm	! set of params in Eq.(19)
   ! Overlap matrix parameters:
   real(8), dimension(10) :: pllm, qllm, rllm, sllm	! set of params in Eq.(20)
end type TB_H_NRL

type, EXTENDS (TB_Hamiltonian) :: TB_H_DFTB ! hamiltonian coefficients:
   ! http://www.dftb.org
   character(200) :: param_name  ! name of parameterization used (or a path to it)
   real(8) :: rcut, d  ! cut-off radius [A] and smoothing distance for Fermi-like cut-off [A]
   real(8) :: Ed, Ep, Es    ! Ed, Ep and Es are the on-site energies for the angular momenta d, p and s for the given atom
   real(8) :: Ud, Up, Us    ! Hubbard U for the angular momenta d, p and s for the given atom
   real(8), dimension(:), allocatable :: Rr ! Radial grid [A]
   real(8), dimension(:,:), allocatable :: Vr ! Hopping integrals [eV]
   real(8), dimension(:,:), allocatable :: Sr ! Overlaps
   ! Reminder:
   !Vr(ir,1) = (s s sigma)
   !Vr(ir,2) = (s p sigma)
   !Vr(ir,3) = (s d sigma)
   !Vr(ir,4) = (p p sigma)
   !Vr(ir,5) = (p p pi)
   !Vr(ir,6) = (p d sigma)
   !Vr(ir,7) = (p d pi)
   !Vr(ir,8) = (d d sigma)
   !Vr(ir,9) = (d d pi)
   !Vr(ir,10) = (d d delta)
end type TB_H_DFTB


type, EXTENDS (TB_Hamiltonian) :: TB_H_3TB ! hamiltonian coefficients:
   ! [1] https://arxiv.org/pdf/2112.11585.pdf
   real(8) :: rcut, d  ! cut-off radius [A] and smoothing distance for Fermi-like cut-off [A]
   real(8) :: rc       ! rescaling coefficient for the distance entering inside Laguerres
   logical :: include_3body   ! include or exclude 3-body parts
   logical :: nullify_diag_cf   ! exclude diagonal part of crystal field
   ! Onsite parameters:
   real(8) :: Ed, Ep, Es      ! Ed, Ep and Es are the on-site energies for the angular momenta d, p and s for the given atom
   real(8), dimension(3,4)   :: Hhavg     ! coefficients before Leguerres in 2-body onsite, Eq.(19) in [1]
   real(8), dimension(3,3,4) :: Hhcf      ! coefficients before Leguerres in crystal-field onsite, Eq.(20) in [1]
   real(8), dimension(3,4)   :: Hh3bdy    ! coefficients before Leguerres in 3-body onsite, Eq.(20) in [1]
   ! Overlap parameters:
   real(8), dimension(10,5)    :: Vrfx    ! coefficients before Leguerres in 2-body Hamiltonian, Eq.(12) in [1]
   real(8), dimension(10,6)    :: Srfx    ! coefficients before Leguerres in 2-body overlap, Eq.(12) in [1]
   real(8), dimension(3,3,3,4) :: V3bdy   ! coefficients before Leguerres in 3-body H, Eq.(14) in [1]
   ! Reminder:
   ! V3bdy (atoms index; shell 1; shell 2; Laguerre parameter index)
   ! "atoms indices" are defined in the function find_3bdy_ind, module "Dealing_with_3TB"
   ! Reminder:
   !Vr(1,5) = (s s sigma)
   !Vr(2,5) = (s p sigma)
   !Vr(3,5) = (s d sigma)
   !Vr(4,5) = (p p sigma)
   !Vr(5,5) = (p p pi)
   !Vr(6,5) = (p d sigma)
   !Vr(7,5) = (p d pi)
   !Vr(8,5) = (d d sigma)
   !Vr(9,5) = (d d pi)
   !Vr(10,5) = (d d delta)
end type TB_H_3TB


type, EXTENDS (TB_Hamiltonian) :: TB_H_BOP ! hamiltonian coefficients:
   ! According to [1] https://arxiv.org/pdf/1909.04561.pdf
   real(8), dimension(14,7) :: H_ci, H_li, H_ni    ! Hamiltonian parameters of Eq.(11)
   real(8), dimension(6,7) :: E_ci, E_li, E_ni      ! On-site parameters of Eq.(12)
   real(8), dimension(14,7) :: S_ci, S_li, S_ni    ! Overlap parameters of Eq.(12)
   real(8) :: rcut, dcut  ! cut-off radius [A] and smoothing distance for Fermi-lik cut-off [A], Eq.(19)
   ! Reminder:
   !Vr(1) = (s s sigma)
   !Vr(2) = (s p sigma)         
   !Vr(3) = (p s sigma)
   !Vr(4) = (s d sigma)         
   !Vr(5) = (d s sigma)
   !Vr(6) = (p p sigma)
   !Vr(7) = (p p pi)
   !Vr(8) = (p d sigma)         
   !Vr(9) = (d p sigma)
   !Vr(10) = (p d pi)               
   !Vr(11) = (d p pi)
   !Vr(12) = (d d sigma)         
   !Vr(13) = (d d pi)
   !Vr(14) = (d d delta)
   ! On-site:
   !E_s(1)  ! s
   !E_s(2)  ! p sigma
   !E_s(3)  ! p pi
   !E_s(4)  ! d sigma
   !E_s(5)  ! d pi
   !E_s(6)  ! d delta
end type TB_H_BOP


type, EXTENDS (TB_Hamiltonian) :: TB_H_xTB ! hamiltonian coefficients:
   ! [1] R.F. Stewart, The Journal of Chemical Physics 52, 431 (1970); doi: 10.1063/1.1672702
   character(200) :: param_name  ! name of parameterization used
   character(10) :: AO_names    ! orbitals descriptor
   real(8) :: rcut, d  ! cut-off radius [A] and smoothing distance for Fermi-lik cut-off [A]
   integer :: Nprim     ! number of primitive gaussian-type orbitals (GTO) in Slater-type orb. (STO) (can be from 1 to 6 [1])
   type(Basis_set_STO), dimension(:), allocatable :: STO ! STO parameteres (via GTO)
   ! unfinished, do not use!
end type TB_H_xTB




! c) van der Waals energy contribution:
type, EXTENDS (TB_vdW) :: TB_vdW_Girifalco
   ! An idea based on this work: [http://www.me.umn.edu/~dtraian/tony-thesis.pdf]
   ! Based on Lenard-Jones potential, van der Waals force can be constructed,
   ! but needs to be cut-off at large distances
   ! and, most importantly, at small distances where covalent bonds dominate
   real(8) :: C12	! [eV*A^12] Lenard-Jones C12
   real(8) :: C6	! [eV*A^6] Lenard-Jones C6
   real(8) :: r_L	! [A] cut-off radius at large distances
   real(8) :: d_L	! [A] cut-off range at large distances
   real(8) :: r_S	! [A] cut-off radius at small distances
   real(8) :: d_S	! [A] cut-off range at small distances
   real(8) :: r_LJ	! [A] shift of coordinate to reproduce equilibrium distance
   real(8) :: dm 	! [A] radius where to switch to polinomial
   real(8) :: d_cut 	! [A] cut-off radius
   real(8) :: a		! fitting polinomial coefficients: a*x^3+b*x^2+c*x+d
   real(8) :: b 
   real(8) :: c
   real(8) :: d
   real(8) :: dsm 	! [A] radius where to switch to polinomial at small distances
   real(8) :: ds_cut 	! [A] cut-off radius on small distances
   real(8) :: as	! fitting polinomial coefficients: a*x^5+b*x^4+c*x^3+d*x^2+e*x+f
   real(8) :: bs 
   real(8) :: cs
   real(8) :: ds
   real(8) :: es
   real(8) :: fs
end type TB_vdW_Girifalco

type, EXTENDS (TB_vdW) :: TB_vdW_Dumitrica ! UNFINISHED
   ! [A. Carlson and T. Dumitrica, Nanotechnology 18 (2007) 065706]
   real(8) :: C6	! [eV*A^6]
   real(8) :: alpha	! [1/A]
end type TB_vdW_Dumitrica

type, EXTENDS (TB_vdW) :: TB_vdW_LJ_cut
   ! Mie potential, generalized LJ in n-exp form: V=eps*( (r0/r)^(2*n) + (r0/r)^n )
   ! reducing to LJ for n=6
   real(8) :: eps    ! [eV*A^12] Lennard-Jones C12
   real(8) :: r0     ! [eV*A^6] Lennard-Jones C6
   real(8) :: n      ! power
   real(8) :: d0_short, dd_short ! [A] [A] cutoff radius and width at short distances
end type TB_vdW_LJ_cut


type, EXTENDS (TB_vdW) :: TB_vdW_ILJ_cut
   ! Improved Lennard-Jones: https://www.mdpi.com/1420-3049/26/13/3906
   ! V=eps*( m/(n-m)*(r0/r)^(n) - n/(n-m)*(r0/r)^m )
   ! reducing to LJ for n=12, m=6
   real(8) :: eps    ! [eV*A^12] Lennard-Jones C12
   real(8) :: r0     ! [eV*A^6] Lennard-Jones C6
   real(8) :: n, m   ! powers
   real(8) :: d0_short, dd_short ! [A] [A] cutoff radius and width at short distances
end type TB_vdW_ILJ_cut



! d) Coulomb energy contribution:
type, EXTENDS (TB_Coulomb) :: TB_Coulomb_cut
   ! Coulomb with smooth cut-off at lond distances:
   integer :: ind    ! 0=softly truncated; 1=Wolf's truncated; 2=Ewalds method
   real(8) :: k		! [N*A^2] coupling constant of Coulomb field (e/(4*Pi*e0))
   real(8) :: dm	! [A] cut-off radius: f_cut = 1/(1+exp((r-dm)/dd))
   real(8) :: dd	! [A] cut-off smoothing distance
   real(8) :: alpha  ! Wolf's truncation parameter
end type TB_Coulomb_cut



!+111111111111111111111111111111111111111111111
! Short-range repulsion in various functional formats:

! Enveloping Fermi-like function: f_cut = 1/(1+exp((r-d0)/dd))
type :: Fermi_cutoff
   real(8) :: d0  ! [A] cut-off radius
   real(8) :: dd  ! [A] cut-off smoothing distance
end type Fermi_cutoff

! Enveloping inverse Fermi-like function: f_cut' = 1 - 1/(1+exp((r-d0)/dd))
type :: Fermi_cutoff_inv
   logical :: use_it ! flag to use this function or not
   real(8) :: d0  ! [A] cut-off radius
   real(8) :: dd  ! [A] cut-off smoothing distance
end type Fermi_cutoff_inv

! Exponential of 1/r: C*exp(1/(r-r0))
type :: Rep_inv_exp
   logical :: use_it  ! flag to use this function or not
   real(8) :: C   ! [eV] energy of the "wall"
   real(8) :: r0  ! [A] "wall" position
   !logical :: use_short_cutoff ! flag to use short-range cutoff or not
end type Rep_inv_exp

! Exponential of r: Phi*exp(-(r-r0)/a)
type :: Rep_exp
   logical :: use_it ! flag to use this function or not
   real(8) :: Phi    ! [eV] energy of the "wall"
   real(8) :: r0     ! [A] "wall" position
   real(8) :: a   ! [A] width
   !logical :: use_short_cutoff ! flag to use short-range cutoff or not
end type Rep_exp

! Power of r: Phi*(r/r0)^m
type :: Rep_pow
   logical :: use_it ! flag to use this function or not
   real(8) :: Phi    ! [eV] energy of the "wall"
   real(8) :: r0     ! [A] position
   real(8) :: m      ! [-] power
   !logical :: use_short_cutoff ! flag to use short-range cutoff or not
end type Rep_pow

! ZBL potential - see separate module
type :: Rep_ZBL
   logical :: use_it
end type Rep_ZBL

! Tabulated potential read from file:
type :: Rep_tab
   logical :: use_it ! flag to use this function or not
   real(8), dimension(:), allocatable :: R      ! [A] array of interatomic distances
   real(8), dimension(:), allocatable :: E      ! [eV/atom] array of repulsive energies
   logical :: use_spline
   ! natural cubic spline parameters: a + b*x + c*x^2 + d*x^3
   real(8), dimension(:), allocatable :: a, b, c, d
end type Rep_tab


! Combined functions together:
type, EXTENDS (TB_Exp_wall) :: TB_Short_Rep
   ! 1) Enveloping Fermi function (cut off at large distances):
   type(Fermi_cutoff) :: f_cut
   ! 1.5) Enveloping inverse Fermi function (cut off at short distances):
   type(Fermi_cutoff_inv) :: f_cut_inv
   ! 2) Inverse exponential:
   type(Rep_inv_exp) :: f_inv_exp
   ! 3) Exponential:
   type(Rep_exp) :: f_exp
   ! 4) Power (array to be able to construct arbitrary polynomial):
   type(Rep_pow), dimension(:), allocatable  :: f_pow
   ! 5) ZBL:
   type(Rep_ZBL) :: f_ZBL
   ! 6) Tabulated potential:
   type(Rep_tab) :: f_tab
end type TB_Short_Rep


! e) Exponential wall at short distances:
type, EXTENDS (TB_Exp_wall) :: TB_Exp_wall_simple
! Exponential wall with smooth cut-off at lond distances:
! C*exp(1/(r-r0)) * fcut
! f_cut = 1/(1+exp((r-dm)/dd))
   real(8) :: C	! [eV] energy of the "wall"
   real(8) :: r0	! [A] "wall" position
   real(8) :: d0	! [A] cut-off radius: f_cut = 1/(1+exp((r-d0)/dd))
   real(8) :: dd	! [A] cut-off smoothing distance
end type TB_Exp_wall_simple




!+111111111111111111111111111111111111111111111
type, EXTENDS (TB_Coulomb) :: Cutie
   character(10) :: Test 
end type Cutie
!==============================================

type Energies
   real(8) :: At_pot       ! [eV] potential energy of atoms
   real(8) :: E_rep        ! [eV] repulsive part of the potential energy
   real(8) :: E_vdW        ! [eV] van der Waals potential energy
   real(8) :: E_coul       ! [eV] Coulomb potential energy (electron emission contribution)
   real(8) :: E_coul_scc   ! [eV] Coulomb potential energy (self-consistent charge contribution)
   real(8) :: E_expwall    ! [eV] Exponential wall potential energy
   real(8) :: At_kin       ! [eV] kinetic energy of atoms
   real(8) :: E_tot        ! [eV] total electron energy in the super-cell
   real(8) :: El_low       ! [eV] total energy of electrons (low-energy domain)
   real(8) :: El_high      ! [eV] total energy of electrons (high-energy domain)
   real(8) :: Eh_tot       ! [eV] total energy of holes
   real(8) :: Total        ! [eV] total energy within the super-cell
   real(8) :: E_supce      ! [eV] total energy of the supercell
   real(8) :: E_glob       ! [eV] total energy (atoms + electrons) in the super-cell
   real(8) :: E_high_heating ! [eV] energy transfer from high-energy electrons to atoms
end type Energies

type Forces
   real(8), dimension(3) :: rep   ! repulsive parts
   real(8), dimension(3) :: att   ! attractive parts
   real(8), dimension(3) :: total ! attractive + repulsive parts, normalized properly
end type Forces

type Supce_force
   real(8), dimension(3,3) :: rep	! repulsive parts for the supercell
   real(8), dimension(3,3) :: att	! attractive parts for the supercell
   real(8), dimension(3,3) :: total	! attractive + repulsive parts for the supercell, normalized properly
   real(8), dimension(3,3) :: total0	! attractive + repulsive parts for the supercell on the last timestep 
end type Supce_force


type Drude
   real(8) :: ReEps, ReEps0 ! real part of the dielectric constant
   real(8) :: ImEps, ImEps0 ! imaginary part of the dielectric constant
   !complex(8) :: Eps_xx, Eps_yy, Eps_zz	! diagonal components of the complex dielectric tensor
   !complex(8) :: Eps_xy, Eps_xz, Eps_yx, Eps_yz, Eps_zx, Eps_zy	! off-diagonal components of the complex dielectric tensor
   complex :: Eps_xx, Eps_yy, Eps_zz	! diagonal components of the complex dielectric tensor
   complex :: Eps_xy, Eps_xz, Eps_yx, Eps_yz, Eps_zx, Eps_zy	! off-diagonal components of the complex dielectric tensor
   real(8) :: n, k ! optical coefficients
   real(8) :: R, T, A ! reflectivity, transmittion, absorbtion
   real(8) :: dc_cond	! dc-conductivity []
   real(8) :: w	  ! frequency of the prob-pulse [1/s]
   real(8) :: l	  ! wavelength of the prob-pulse [nm]
   real(8) :: tau ! FWHM duration of probe pulse [fs], to convolve with
   real(8) :: me_eff, mh_eff	! [kg] effective mass of CB electron and VB hole
   real(8) :: tau_e, tau_h	! [fs] mean scattering times of electrons and holes
   real(8) :: teta	! angle of prob-pulse [radians]
   real(8) :: dd 	! [nm] experimental layer thickness
   logical :: all_w	! calculate for all hw or not
   logical :: KK	! calculate real part of CDF via Kramers-Kronig relation or not
   real(8) :: E_min, E_max, dE ! min, max, and step for the array of Eps_hw (given in input file)
   real(8), dimension(:,:), allocatable :: Eps_hw ! array of all eps vs hw
   ! Reminder:
!    Eps_hw(1,i) = hw     ! energy [eV]
!    Eps_hw(2,i) = Re_eps ! real part of CDF
!    Eps_hw(3,i) = Im_eps ! imaginary part of CDF
!    Eps_hw(4,i) = LF  ! loss function
!    Eps_hw(5,i) = R   ! reflectivity
!    Eps_hw(6,i) = T   ! transmission
!    Eps_hw(7,i) = A   ! absorption
!    Eps_hw(8,i) = n   ! optical n
!    Eps_hw(9,i) = k   ! optical k
!    Eps_hw(10,i) = k   ! dc-conductivity
!    Eps_hw(11,i) = Re_E_xx
!    Eps_hw(12,i) = Im_E_xx
!    Eps_hw(13,i) = Re_E_yy
!    Eps_hw(14,i) = Im_E_yy
!    Eps_hw(15,i) = Re_E_zz
!    Eps_hw(16,i) = Im_E_zz
end type Drude


!==============================================
! Atoms as objects:
type Atom
   integer :: KOA		! kind of atom in the compound (according to type At_data)
   real(8), dimension(3) :: R	! [A] coordinates (x,y,z)
   real(8), dimension(3) :: S	! [a.u.] relative coordinates (sx, sy, sz)
   real(8), dimension(3) :: V	! [A/fs] velocities (Vx, Vy, Vz)
   real(8), dimension(3) :: SV	! [a.u.] relative velocities (SVx, SVy, SVz)
   real(8), dimension(3) :: A	! [A^2/fs] accelerations
   ! Forces acting on atoms:
   type(Forces) :: forces	! all interatomic forces
   ! on the last time-step:
   real(8), dimension(3) :: R0	! [A] coordinates
   real(8), dimension(3) :: S0	! [a.u.] relative coordinates
   real(8), dimension(3) :: V0	! [A/fs] velocities
   real(8), dimension(3) :: SV0	! [a.u.] relative velocities
   real(8), dimension(3) :: A0	! [A^2/fs] accelerations
   ! Kinetic energy:
   real(8) :: Ekin	! [eV]
   ! Potential energy:
   real(8) :: Epot	! [eV]
   ! On the first time-step, equilibrium positions of atoms to get their mean displacements:
   real(8), dimension(3) :: R_eq	! [A] coordinates (x,y,z)
   real(8), dimension(3) :: S_eq	! [a.u.] relative coordinates (sx, sy, sz)
   ! For Martyna algorithm, fictitious parameters:
   real(8), dimension(3) :: A_tild, v_F, v_J, A_tild0, v_F0, v_J0
   ! For linear scaling TB, subcell numbers this atom belongs to:
   integer :: Nx_subcel, Ny_subcel, Nz_subcel
   ! Additional data:
   real(8) :: q   ! charge (Mulliken)
   real(8), dimension(3) :: accel, accel0  ! acceleration for calculation of configurational temperature
end type Atom

type :: MC_atoms ! to treat holes in each shell of each atom
   real(8), dimension(:), allocatable :: Noh	! total number of holes in all shells of this atom
end type MC_atoms


type :: Orbital_resolved_data ! to count electrons in orbitals separately
   real(8), dimension(:), allocatable :: Ne  ! number of electrons in each orbital
   real(8), dimension(:), allocatable :: Ee  ! energy of electrons in each orbital
   real(8), dimension(:,:), allocatable :: fe  ! electron distribution function over TB levels in each orbital
end type Orbital_resolved_data


!==============================================
! Subcells for linear scaling TB (NOT READY):
type Sub_cell
   real(8), dimension(:), allocatable :: fe ! low-energy electron distribution
   real(8) :: mu	! [eV] electron chemical potential
   real(8) :: TeeV, TaeV ! [eV] electrons and atoms temperatures
   real(8), dimension(:), allocatable :: Ta_sub    ! [K] temperatures of different sublattices
   real(8) :: Ne_low	! current number of low-energy electrons of VB and CB
   type(Energies) :: nrg		! [eV] energies in the super-cell
   real(8), dimension(:,:,:), allocatable :: DOS_weights ! to identify and separate different bands
   real(8) :: E_gap	! [eV] current band gap
   real(8), dimension(:,:), allocatable :: Ha	! hamiltonian matrix
   real(8), dimension(:,:), allocatable :: Ha0	! hamiltonian matrix on the last step
   real(8), dimension(:,:), allocatable :: H_non	! nondiagonalized Hamiltonian
   real(8), dimension(:,:), allocatable :: H_non0	! nondiagonalized Hamiltonian at last step
   real(8), dimension(:,:), allocatable :: Sij	! Overlap matrix for non-orthogonal TB
   real(8), dimension(:,:), allocatable :: Hij	! Non-orthogonal TB Hamiltonian
   real(8), dimension(:,:), allocatable :: Hij_sol	! Eigenvectors of non-orthogonal TB Hamiltonian
   real(8), dimension(:), allocatable :: Ei	! energy levels, eigenvalues of the hamiltonian matrix
   real(8), dimension(:), allocatable :: Ei0	! energy levels, eigenvalues of the hamiltonian matrix on the last step
   real(8), dimension(:,:), allocatable :: Aij	! coefficients used for forces in TB
   ! Complex Hamiltonian for multiple k-points:
   complex, dimension(:,:,:,:,:), allocatable :: CHa	! Complex hamiltonian matrix for each (kx, ky, kz) points
   complex, dimension(:,:,:,:,:), allocatable :: CHa0	! Complex hamiltonian matrix on the last step
end type Sub_cell


!==============================================
type Displacement_analysis
   character(100) :: mask_name   ! name of the mask
   real(8) :: MSD_power          ! power of the mean displacement for this particular analysis
   real(8) :: mean_disp          ! [A^MSD_power] mean displacement of all atoms
   real(8), dimension(:), allocatable :: mean_disp_sort     ! [A^MSD_power] mean displacement of atoms by sort
   logical, dimension(3) :: print_r       ! along which exis the user requirested the analysis
   real(8), dimension(3) :: mean_disp_r   ! [A^MSD_power] mean displacements along X, Y, Z axes
   real(8), dimension(:,:), allocatable :: mean_disp_r_sort   ! [A^MSD_power] mean displacements along X, Y, Z axes by sort
   logical, dimension(:), allocatable :: Atomic_mask  ! atomic mask to be used
   ! Definition of spatial-section mask:
   ! Currently, only 2 sections are supported, connected by a logical operator 'and' or 'or':
   logical :: logical_and, logical_or  ! which one is used
   integer, dimension(2) :: axis_ind   ! index of the axis for the section
   real(8), dimension(2,3) :: r_start, r_end   ! section starting and ending points
end type Displacement_analysis



!==============================================
! Supercell as object:
type Super_cell
   ! Sub-cells for linear scaling TB, if used:
   type(Sub_cell), dimension(:,:,:), allocatable :: Subcell ! Subcells along 3 axes: X, Y, Z (currently unused)
   ! Data for the entire simulation box:
   ! How many electrons and atoms in the super-cell:
   integer :: Na, Ne	! number of atoms, electrons, in the super-cell
   real(8) :: Ne_low	! current number of low-energy electrons of VB and CB
   ! Results of MC run:
   real(8) :: Ne_high	! current number of high-energy electrons of CB
   real(8) :: Ne_emit	! current number of emitted electrons (above the work function)
   real(8) :: Nh	! current number of deep-shell holes
   type(MC_atoms), dimension(:), allocatable :: MChole ! number of holes in each shell of each atom
   real(8) :: Nph	! current number of absorbed photons
   real(8) :: Q		! mean unballanced charge per atom: <Ne_emit/Na>
   ! Other data:
   real(8) :: Ne_CB	! current number of electrons in CB only (above the band gap)
   real(8) :: Te, Ta, Tconf, Tconf2   ! [K] temperature of electrons, kinetic temperature of atoms, (configurational Ta)
   real(8) :: Ta_var(8)     ! [K] various definition of temperature: kinetic, configurational, etc.
   real(8) :: Ta_r_var(6)   ! [K] dimension-resolved temperatures along X,Y,Z
   real(8) :: Ta_conf_run_average(10)     ! [K] laset 10 points of Ta_config for running average
   real(8) :: Fv  ! <F*v> force-velocity correlator
   real(8), dimension(:), allocatable :: Ta_sub    ! [K] temperatures of different sublattices
   real(8) :: TeeV, TaeV ! [eV] electrons and atoms temperatures
   real(8) :: Pressure	! [Pa] pressure in the atomic system
   real(8), dimension(3,3) :: Stress	! [Pa] stress tensor in the atomic system
   real(8) :: Pot_Pressure  ! [Pa] potential contribution to pressure (without kinteic energy of atoms)
   real(8), dimension(3,3) :: Pot_Stress	! [Pa] potential contribution to stress tensor (without kinetic energies)
   real(8) :: MSD	! [A^2] mean square displacements average over all atoms
   real(8), dimension(:), allocatable :: MSDP	! [A^2] mean square displacements for atoms of different sorts
   type(Displacement_analysis), dimension(:), allocatable :: Displ ! [A] mean displacements for atoms masked
   real(8) :: mu	! [eV] electron chemical potential
   real(8), dimension(:), allocatable :: fe ! low-energy electron distribution
   real(8), dimension(:), allocatable :: fe_eq ! equivalent Fermi electron distribution
   real(8) :: Se, Se_eq  ! electron entropy [K/eV], and equivalent equilibrium entropy
   ! Atomic distribution function:
   real(8), dimension(:), allocatable :: Ea_grid, Ea_grid_out   ! energy grid for atomic distribution
   real(8), dimension(:), allocatable :: fa, fa_eq, fa_out, fa_eq_out ! atomic distribution and equivalent Maxwell distribution
   real(8), dimension(:), allocatable :: Ea_pot_grid, Ea_pot_grid_out   ! energy grid for atomic distribution
   real(8), dimension(:), allocatable :: fa_pot, fa_eq_pot, fa_pot_out, fa_eq_pot_out ! atomic distribution of potential energies
   real(8), dimension(:), allocatable :: Ea_tot_grid, Ea_tot_grid_out   ! total-energy grid for atomic distribution
   real(8), dimension(:), allocatable :: fa_tot, fa_tot_out    ! atomic distribution of total energies
   real(8) :: Sa, Sa_eq, Sa_eq_num  ! atomic entropy [K/eV], and equivalent equilibrium entropy (analytical & numerically calculated)
   real(8) :: Sa_conf, Sa_tot   ! configurational and total atomic entropy [K/eV]
   real(8) :: Pot_distr_E_shift  ! [eV] shift of the distribution
   ! Separate for VB and CB, if needed:
   real(8) :: Ne_low_CB, Ne_low_VB, El_low_CB, El_low_VB ! number and energy of electrons in VB and CB
   real(8) :: Te_VB, Te_CB ! [eV] electron temperatures for conduction and valence bands
   real(8) :: mu_VB, mu_CB ! [eV] electron chemical potentials for conduction and valence bands
   real(8), dimension(:), allocatable :: fe_eq_CB  ! equivalent Fermi electron distribution in CB
   real(8), dimension(:), allocatable :: fe_eq_VB  ! equivalent Fermi electron distribution in VB
   real(8) :: Se_VB, Se_eq_VB, Se_CB, Se_eq_CB  ! equivalent electron entropy in VB and CB [K/eV]
   ! Orbital-specified parameters:
   type(Orbital_resolved_data), dimension(:), allocatable :: Orb_data ! electron data in different orbitals of each kind of atom
   ! grid for electron distribution; distribution on this grid; high-energy electrons distribution on grid:
   real(8), dimension(:), allocatable :: E_fe_grid, fe_on_grid, fe_high_on_grid  ! electron spectrum on grid (fe*DOS)
   real(8), dimension(:), allocatable :: fe_norm_on_grid, fe_norm_high_on_grid   ! electron distribution on grid (fe)
   real(8), dimension(:), allocatable :: I_ij ! electron-ion collision integral [1/s]
   real(8), dimension(:), allocatable :: Norm_WF ! Normalization of wave functions
   real(8) :: Ce  ! electron heat capacity [J/(m^3 K)]
   real(8), dimension(:), allocatable :: Ce_part   ! band-resolved electron heat capacity [J/(m^3 K)]
   real(8), dimension(:), allocatable :: Ce_i   ! energy level resolved electron heat capacity [J/(m^3 K)]
   real(8) :: kappa_e  ! electron heat conductivity [W/(m K)]
   real(8), dimension(:), allocatable :: kappa_e_vs_Te  ! total electron heat conductivity [W/(m K)] vs Te [K]
   real(8), dimension(:), allocatable :: kappa_ee_vs_Te  ! e-e contribution to electron heat conductivity [W/(m K)] vs Te [K]
   real(8), dimension(:), allocatable :: kappa_Te_grid  ! Te grid for electron heat conductivity [K]
   real(8), dimension(:), allocatable :: kappa_mu_grid  ! mu vs Te (in electron heat conductivity) [eV]
   real(8), dimension(:), allocatable :: kappa_Ce_grid  ! Ce vs Te (in electron heat conductivity) [J/(m^3 K)]
   real(8), dimension(:), allocatable :: kappa_e_part   ! band-resolved electron heat conductivity [W/(m K)]
   real(8), dimension(:,:), allocatable :: Mij ! matrix element for electron-ion coupling or dynamic kappa
   ! Atoms:
   type(Atom), dimension(:), allocatable :: MDAtoms ! all atoms in MD
   type(Energies) :: nrg		! [eV] energies in the super-cell
   ! Hamiltonians, wave-functions and related parameters:
   real(8), dimension(:,:), allocatable :: Ha	! hamiltonian matrix
   real(8), dimension(:,:), allocatable :: Ha0	! hamiltonian matrix on the last step
   real(8), dimension(:,:), allocatable :: H_non	! nondiagonalized Hamiltonian
   real(8), dimension(:,:), allocatable :: H_non0	! nondiagonalized Hamiltonian at last step
!    real(8), dimension(:,:), allocatable :: Vij	! Radial part of the hamiltonian matrix
   real(8), dimension(:,:), allocatable :: Sij	! Overlap matrix for non-orthogonal TB
   real(8), dimension(:,:), allocatable :: Hij	! Non-orthogonal TB Hamiltonian
   real(8), dimension(:,:), allocatable :: Hij_sol ! Eigenvectors of non-orthogonal TB Hamiltonian
   real(8), dimension(:), allocatable :: eigen_S   ! eigenvalues of the hamiltonian overlap matrix
   real(8), dimension(:), allocatable :: Ei	    ! energy levels, eigenvalues of the hamiltonian matrix
   real(8), dimension(:), allocatable :: Ei0	! energy levels, eigenvalues of the hamiltonian matrix on the last step
   real(8), dimension(:), allocatable :: Ei_scc_part  ! eigenvalues of the non-SCC part of the hamiltonian
   real(8), dimension(:,:), allocatable :: Aij	! coefficients used for forces in TB
   ! Complex Hamiltonian for multiple k-points:
   complex, dimension(:,:,:,:,:), allocatable :: CHa	! Complex hamiltonian matrix for each (kx, ky, kz) points
   complex, dimension(:,:,:,:,:), allocatable :: CHa0	! Complex hamiltonian matrix on the last step
   ! Parameters of TB:
   class(TB_repulsive), allocatable, dimension(:,:)   :: TB_Repuls   ! parameters of the repulsive part of TB (shape to be defined)
   class(TB_Hamiltonian), allocatable, dimension(:,:) :: TB_Hamil ! parameters of the Hamiltonian of TB (shape to be defined)
   class(TB_vdW), allocatable, dimension(:,:) :: TB_Waals         ! parameters of the van der Waals for TB (shape to be defined)
   class(TB_Coulomb), allocatable, dimension(:,:) :: TB_Coul      ! parameters of the Coulomb together with TB (shape to be defined)
   class(TB_Exp_wall), allocatable, dimension(:,:) :: TB_Expwall  ! parameters of the exponential wall potential (shape to be defined)
   ! Forces:
   type(Supce_force) :: SCforce ! forces acting on the supercell
   ! Super-cell size and velocities (used within Parrinello-Rahman):
   real(8), dimension(3,3) :: supce 	! [A] length of super-cell
   real(8), dimension(3,3) :: supce0 	! [A] length of super-cell on the previous time-step
   real(8), dimension(3,3) :: supce_t	! [A] transposed super-cell
   real(8), dimension(3,3) :: GG	    ! [A^2] super-cell^2, g-matrix
   real(8), dimension(3,3) :: Vsupce    ! Derivatives of Super-cell vectors (velocities)
   real(8), dimension(3,3) :: Vsupce0   ! Derivatives of Super-cell vectors (velocities) on last time-step
   real(8), dimension(3,3) :: k_supce 	! [1/A] length of the reciprocal super-cell
   real(8), dimension(3,3) :: supce_eq 	! [A] equilibrium lengths of super-cell
   real(8) :: V      ! super-cell volume [A^3]
   ! For atmoic calculations, lists of nearest neighbors:
   integer, dimension(:,:), allocatable :: Near_neighbor_list  	! list of nearest neighbors for all atoms
   real(8), dimension(:,:,:), allocatable :: Near_neighbor_dist	! distances to this atoms [A]
   real(8), dimension(:,:,:), allocatable :: Near_neighbor_dist_s	! distances to this atoms, relative
   integer, dimension(:), allocatable :: Near_neighbor_size  	! how many nearest neighbours there are for each atom
   integer, dimension(:), allocatable :: Near_neighbors_user  	! how many nearest neighbours within radius defined by the user
   ! For dielectric function:
   type(Drude) :: eps	! epsylon, dielectric function and its parameters
   real(8), dimension(:,:), allocatable :: PRRx	! matrix element for dielectric function (along x)
   real(8), dimension(:,:), allocatable :: PRRy	! matrix element for dielectric function (along y)
   real(8), dimension(:,:), allocatable :: PRRz	! matrix element for dielectric function (along z)
   complex, dimension(:,:), allocatable :: cPRRx ! matrix element for dielectric function (along x)
   complex, dimension(:,:), allocatable :: cPRRy ! matrix element for dielectric function (along y)
   complex, dimension(:,:), allocatable :: cPRRz ! matrix element for dielectric function (along z)
   !Energy parameters:
   real(8) :: E_gap	! [eV] current band gap
   integer :: N_Egap	! HOMO level to define bandgap
   real(8) :: E_bottom	! [eV] current bottom of the conduction band
   real(8) :: E_top	! [eV] current top of the conduction band
   real(8) :: E_VB_bottom	! [eV] current bottom of the valence band
   real(8) :: E_VB_top	! [eV] current bottom of the valence band
   real(8) :: G_ei	! [W/(m^3 K)] electron-ion coupling parameter
   real(8), dimension(:,:), allocatable :: G_ei_partial	! [W/(m^3 K)] partial electron-ion coupling parameter per all orbitals pairwise
   real(8), dimension(:,:), allocatable :: DOS	! DOS
   real(8), dimension(:,:,:), allocatable :: partial_DOS	! partial DOS made of different orbitals
   ! Testmode additional data (usually not needed):
   real(8) :: V_CoM(3)     ! center of mass velosity
   real(8) :: I_tot(3,3)   ! moment of inertia tensor
   real(8) :: F_tot(3)     ! total force in the supercell
end type Super_cell


!==============================================
! Material and laser parameters:
type MFP
   real(8), dimension(:), allocatable :: E	! [eV] energy
   real(8), dimension(:), allocatable :: L	! [A] mean free path
end type MFP

type Ritchi
   real(8), dimension(:), allocatable :: A, E0, G	! coefficients of Ritchi CDF
end type Ritchi

type At_data
   character(3) Name	! Chemical element name
   real(8) :: Ma	! [kg] atomic mass
   real(8) :: Z		! atomic number
   integer :: sh	! number of shells of the element
   real(8) :: percentage ! contribution to the compound
   integer :: NVB	! number valence electrons according to periodic table
   real(8) :: Hubbard_U ! [eV] Hubbard U for SCC
   type(Basis_set), dimension(:), allocatable :: Cart_Basis  ! Cartesian GTO basis set functions for this element (if xTB is used)
   type(Basis_set), dimension(:), allocatable :: Spher_Basis  ! Spherical (pure) GTO basis set functions for this element (if xTB is used)
   real(8) :: mulliken_Ne   ! electron population according to Mulliken analysis
   real(8) :: mulliken_q   ! charge according to Mulliken analysis
   integer, dimension(:), allocatable :: Shl_dsgnr  ! EADL shell designator
   integer, dimension(:), allocatable :: Shl_dsgnr_atomic  ! EADL shell designator for atomic shells without VB
   character(11), dimension(:), allocatable :: Shell_name   ! names of the shells
   real(8), dimension(:), allocatable :: Ip   ! [eV] ionization potentials for all shells
   real(8), dimension(:), allocatable :: Ek   ! [eV] mean kinetic energy of all shells
   real(8), dimension(:), allocatable :: Ne_shell  ! number of electrons in each shell
   !logical, dimension(:), allocatable :: valent    ! flag if this shell is valent or core-shell
   real(8), dimension(:), allocatable :: Auger  ! [fs] Auger-decay times for all shells
   integer, dimension(:), allocatable :: TOCS   ! type of electron scattering cross-section used for each shell
   integer, dimension(:), allocatable :: N_CDF  ! number of CDF-functions in each shell
   integer, dimension(:), allocatable :: TOCSph ! type of photon absorption cross-section used for each shell
   type(Ritchi), dimension(:), allocatable :: CDF  ! coefficients of CDF
   real(8), dimension(:), allocatable :: Nh_shell  ! current number of deep-shell holes in each shell
   type(MFP), dimension(:), allocatable :: El_MFP  ! electron inelastic mean free paths for each shell (inversed [1/A])
   type(MFP), dimension(:), allocatable :: Ph_MFP  ! photon mean free paths for each shell (inversed [1/A])
   type(MFP), dimension(:), allocatable :: El_MFP_vs_T ! electron MFP for (inversed [1/A]) for different Te
   type(MFP) :: El_EMFP ! electron elastic mean free paths (inversed [1/A])
end type At_data

type Solid
   character(100) Name	! material name
   character(100) Chem	! material composition (chemical formulae)
   integer :: N_KAO	! number of different kinds of atoms in this compound
   type(At_data), dimension(:), allocatable :: Atoms	! all kinds of atoms of the compound
   type(MFP) :: El_MFP_tot  ! Total electron inelastic mean free paths (inversed [1/A])
   type(MFP) :: El_EMFP_tot ! Total electron elastic mean free paths (inversed [1/A])
   type(MFP) :: Ph_MFP_tot  ! Total photon mean free paths for each shell (inversed [1/A])
   real(8), dimension(:,:), allocatable :: PCF	! pair correlation function (if required)
   integer :: cell_x, cell_y, cell_z	! number of unit-cells in x, y and z directions
   real(8) :: W_PR	! [kg] Parinello_Rahman super-cell mass
   real(8) :: p_ext	! [Pa] external pressure applied
   real(8) :: dens	! [g/cm^3] density of the material to be used in MC in case of BEB-cross-sections
   real(8) :: At_dens   ! [1/cm^3] atomic density of the material to be used in MC in case of BEB-cross-sections
   real(8) :: hw_plasma ! [eV] plasmon energy
   real(8) :: T_bath	! [K] Bath temeprature for Berendsen thermostat for atoms
   real(8) :: T_bath_e	! [K] Bath temeprature for Berendsen thermostat for electrons
   real(8) :: tau_bath	! [fs] time constant of cooling via thermostat for atoms
   real(8) :: tau_bath_e	! [fs] time constant of cooling via thermostat for electrons
end type Solid


type Pulse
   integer :: KOP    ! kind of pulse: 0 = flat-top, 1 = Gaussian, 2 = SASE
   real(8) :: hw     ! [eV] photon energy
   real(8) :: FWHM_hw   ! [eV] distribution of photon energy spectrum (assumed gaussian)
   real(8) :: t      ! [fs] pulse duration
   real(8) :: t0     ! [fs] pulse maximum position
   real(8) :: F_in   ! [J/cm^2]  incoming fluence
   real(8) :: F      ! [eV/atom] absorbed dose
   real(8) :: Fabs   ! [eV] total absorbed energy per simulation box
   real(8) :: Nph    ! number of absorbed photons
end type Pulse


! Atomic parameters provided by the user, to overwrite the default values
type User_overwrite_data
   logical :: do_overwrite
   character(3), dimension(:), allocatable :: name ! element name
   real(8), dimension(:), allocatable :: mass   ! [a.m.u.] atomic mass
   real(8), dimension(:,:), allocatable :: Ip   ! [eV] ionization potentials for all shells
   real(8), dimension(:,:), allocatable :: Ek   ! [eV] mean kinetic energy of all shells
   real(8), dimension(:,:), allocatable :: Ne_shell   ! number of electron in each shell
   real(8), dimension(:,:), allocatable :: auger   ! [fs] Auger decay times
endtype User_overwrite_data


! MPI parameters
type :: Used_MPI_parameters
    ! Processes characteristics:
    integer :: process_rank     ! index of the process
    integer :: size_of_cluster  ! total number of processes
    integer :: ierror           ! error handler
    integer :: INFO             ! info to use in subroutines calls
    !----------------
    character(10) :: rank_ch    ! rank process as character (for printout)
    real(8) :: Wt0, Wt1   ! wall time defined by MPI
end type Used_MPI_parameters



type Numerics_param
   ! MPI parameters:
   type(Used_MPI_parameters) :: MPI_param
   ! Subcell parameters for linear scaling TB (not ready):
   integer :: lin_scal    ! use linear scaling TB (1), or not (0)
   integer, dimension(3) :: N_subcels   ! number of subcells along each axis: X, Y, Z
   real(8), dimension(:), allocatable :: Subcell_coord_sx, Subcell_coord_sy, Subcell_coord_sz
   ! Other parameters:
   integer :: which_input  ! number of input file used (for using more then one sequentially)
   logical :: verbose, nonverbose
   logical :: redo_MFP     ! flag to recalculate mean free paths, if needed
   logical :: print_MFP    ! flag to printout mean free paths
   ! Electronic distribution function parameters:
   logical :: fe_input_exists ! flag to use the distribution from a file
   character(100) :: fe_filename ! file name with user-provided initial electronic distribution
   real(8), dimension(:), allocatable :: fe_input  ! initial distribution function
   real(8), dimension(:), allocatable :: high_DOS  ! DOS for high-energy electron distribution, if required
   integer :: fe_aver_num  ! number of time-steps over which to average the distribution on the grid
   real(8) :: tau_fe ! [fs] characteristic relaxation time of ALL electrons (used for the relaxation time approximation)
   real(8) :: tau_fe_CB ! [fs] characteristic relaxation time of CB electrons (used for the relaxation time approximation)
   real(8) :: tau_fe_VB ! [fs] characteristic relaxation time of VB electrons (used for the relaxation time approximation)
   logical :: do_partial_thermal ! flag for band-resolved thermalization
   ! MD:
   real(8) :: dt	      ! [fs] time-step for MD
   real(8) :: halfdt      ! dt/2, often used
   real(8) :: dtsqare	  ! dt*dt/2, often used
   real(8) :: dt3	      ! dt^3/6, often used
   real(8) :: dt4	      ! dt^4/48, often used
   integer :: MD_algo     ! 0=Verlet (2d order); 1=Yoshida (4th order); 2=Martyna(4th order)
   real(8), dimension(:), allocatable :: dt_MD_reset_grid   ! grid, when to change the MD timestep
   real(8), dimension(:), allocatable :: dt_MD_grid         ! grid, which MD timestep to use
   integer :: i_dt        ! which timestep from the array "dt_MD_grid" to use now
   character(100) :: MD_step_grid_file   ! filename with MD time step grid
   !logical, dimension(:), allocatable :: Atomic_masks ! user-defined masks for atomic analysis
   logical :: print_Ta ! flag for various atomic temperature definitions
   integer :: ind_starting_V  ! index to set starting velocity distribution: 1=linear; 2=Maxwellian
   logical :: vel_from_file   ! index to mark whether velocities were read from a file or set
   !-----------------
   ! EPICS data: EADL, EPDL, EEDL:
   character(200) :: EADL_file, EPDL_file, EEDL_file
   !-----------------
   real(8), dimension(:), allocatable :: At_bath_reset_grid   ! grid, when to change the Atomic bath parameters
   real(8), dimension(:), allocatable :: At_bath_grid_Ta         ! Atomic bath temperatures array [K]
   real(8), dimension(:), allocatable :: At_bath_grid_tau        ! Atomic bath characteristic times array [fs]
   integer :: i_At_bath_dt        ! which timestep from the array to use now
   character(100) :: At_bath_step_grid_file   ! filename with Atomic bath parameters
   real(8), dimension(:), allocatable :: El_bath_reset_grid   ! grid, when to change the Electronic bath parameters
   real(8), dimension(:), allocatable :: El_bath_grid_Ta         ! Electronic bath temperatures array [K]
   real(8), dimension(:), allocatable :: El_bath_grid_tau        ! Electronic bath characteristic times array [fs]
   integer :: i_El_bath_dt        ! which timestep from the array to use now
   character(100) :: El_bath_step_grid_file   ! filename with Electronic bath parameters
   !-----------------
   real(8) :: t_start	  ! [fs] starting time of simulation
   real(8) :: t_total	  ! [fs] total time of simulation
   real(8) :: t_Te_Ee	  ! time when we switch from Te=const, to Ee=const [fs] / negative value, when not using it
   real(8) :: t_NA	      ! time when we switch on nonadiabatic terms [fs]
   real(8) :: dt_cooling  ! time-step how often to cool down the atoms [fs]
   logical :: p_const	  ! P=const, otherwise V=const for Parinello-Rahman MD simulations
   logical :: Nonadiabat  ! true=included / false=excluded nonadiabatic coupling
   integer :: NA_kind     ! non-adiabatic model kind
   integer :: ind_at_distr ! model for atomic distribution in the nonadiabatic coupling
   logical :: scc         ! include self-consistent charge corrections in TB?
   integer :: scc_gam_ind ! index for the model for gamma in scc term
   real(8) :: scc_mix     ! mixin factor for scc calculations
   integer :: basis_size_ind  ! index for the size of the basis set used: s=0; sp3=1; sp3d5=2; sp3s*=3; sp3d5s*=4;
   integer :: N_basis_size    ! basis set size (orbitals per atom)
   logical :: Transport		 ! true=included / false=excluded ; for atoms
   logical :: Transport_e		 ! true=included / false=excluded ; for electrons
   real(8) :: E_cut	    ! [eV] cut-off energy, separating low-energy-electrons from high-energy-electrons
   logical :: E_cut_dynamic ! include evolution of E_cut due to changes of the band struture
   real(8) :: E_work	! [eV] cut-off energy (work function) above which electrons are emitted from the material
   real(8) :: E_Coulomb	! [eV] Coulomb energy attracting electrons back to the material
   integer :: NMC	! number of MC iterations
   integer :: NOMP	! Number of threads for openmp
   real(8) :: dt_save 	! save data into files every 'dt_save' [fs]
   real(8) :: fctr	! multiplication factor often used in MD
   real(8) :: acc_window	! [eV] acceptance window for nonadiabatic coupling
   real(8) :: degeneracy_eV	! [eV] window to exclude quasidegenerate levels in nonadiabatic coupling
   real(8) :: M2_scaling    ! scaling factor for electron-ion coupling matrix element
   real(8) :: at_cool_start, at_cool_dt	! [fs], cooling of atoms: when to start, and how often to do
   real(8) :: Smear_DOS	! smearing to use for DOS calculations
   real(8) :: NN_radius ! user-defined radius to count nearest neighbors (for printout only, not for MD calculations!)
   integer :: MSD_power ! power of mean displacement to print out (set integer N: <u^N>-<u0^N>)
   integer :: DOS_splitting ! model to split DOS
   logical, dimension(:,:,:), allocatable :: mask_DOS ! to identify and separate different bands
   integer :: weigthed_DOS  ! use weights on energy levels accourding to DOS_weights
   real(8), dimension(:,:,:), allocatable :: DOS_weights ! to identify and separate different bands
   character(200) :: input_path	        ! input folder address
   character(200) :: output_path	! output folder address
   character(1) :: path_sep   ! path separator
   character(10) :: At_base   ! where to take atomic data from (EADL, CDF, XATOM...)
   real(8) :: user_defined_E_gap ! [eV] to start with, user-defined gap
   ! Setting supercell for biomolecules, embedding in water:
   logical :: embed_water  ! flag for embedding in water
   integer :: N_water_mol  ! how many water molecules to use
   ! numbers of files:
   integer :: FN_temperatures, FN_energies, FN_atoms_R, FN_atoms_S, FN_supercell, FN_electron_properties, FN_numbers, FN_all_w
   integer :: FN_deep_holes, FN_Ei, FN_fe, FN_PCF, FN_optics, FN_parameters, FN_communication, FN_cif, FN_pressure, FN_DOS
   integer :: FN_coupling, FN_neighbors, FN_Ce, FN_kappa, FN_kappa_dyn, FN_Se, FN_fe_on_grid, FN_Te, FN_mu, FN_orb_resolved
   integer :: FN_fa, FN_Sa, FN_Ta, FN_fa_pot, FN_Ta_part, FN_fa_tot, FN_testmode
   integer, dimension(:), allocatable :: FN_displacements
   integer :: MOD_TIME ! time when the communication.txt file was last modified
   integer :: drude_ray, optic_model
   integer :: el_ion_scheme
   integer :: ixm, iym, izm	! number of k-points in each direction for eps-calculations
   real(8), dimension(:,:), allocatable :: k_grid	! for the case of user-provided grid for k-space (for CDF and DOS calculations)
   logical :: r_periodic(3)	! periodic boundaries in each of the three spatial dimensions
   ! Different output, what to save:
   logical :: save_Ei, save_fe, save_PCF, save_XYZ, do_drude, do_cool, do_atoms, change_size, allow_rotate
   logical :: save_fe_grid, save_fe_orb, save_fa, save_testmode
   ! Reminder: codes of save_XYZ_extra indices: (1) atomic mass; (2) atomic charge; (3) kinetic energy
   logical :: save_XYZ_extra(3)  ! additional properties of atoms to print (or not)
   logical :: do_elastic_MC, do_path_coordinate, do_kappa, do_DOS, do_kappa_dyn
   logical :: save_CIF, save_pressure, save_DOS, save_raw, save_NN, save_CDF
   integer :: Mulliken_model
   integer :: ind_fig_extention, change_size_step
   real(8) :: change_size_max, change_size_min
   real(8) :: kappa_Te_min, kappa_Te_max, kappa_dTe
   integer :: kappa_model
   character(200) :: input_CDF_file
   character(4) :: fig_extention
   ! BOP parameters creation:
   logical :: create_BOP_repulse
   character(200) :: BOP_Folder_name, Filename_communication
   real(8) :: BOP_bond_length   ! [A]
   ! Potential DOS power parameter (used for testing):
   real(8) :: power_b
   ! Initial cell data file:
   character(200) :: Cell_filename
   ! General flags:
   integer :: save_files_used ! to mark for the output printing out, whether SAVE-files are used or unit cell
   logical :: numpar_in_input ! numpar were already read from the single file, no need for separate file reading
   character(200) :: output_extra_name, output_name_add ! user-defined additional text in the output folder name
end type Numerics_param
!==============================================



!==============================================
! types used in Monte Carlo:
! One photon as an object, to make an array of them:
type :: Particle ! basic class for all particles
   real(8) E ! energy [eV]
   real(8) ti ! time of impact
end type Particle

type, EXTENDS (Particle) :: Photon ! photon as an object, contains the following info:
end type Photon

type, EXTENDS (Particle) :: Electron ! electron as an object contains the same into as photon + more:
   real(8) :: t ! current time [fs]
   integer :: colls ! number of collisions that the electron made
end type Electron

type, EXTENDS (Particle) :: Hole ! hole as an object, contains the same info as electron + more:
   real(8) :: t ! current time [fs]
   integer :: KOA ! kind of atom
   integer :: Sh  ! shell
end type Hole

type MC_data	! all MC arrays storred here
   type(Photon), dimension(:), allocatable :: photons
   type(Electron), dimension(:), allocatable :: electrons
   type(Hole), dimension(:), allocatable :: holes
   integer :: noe, noh_tot ! numbers of existing electron and holes
   integer :: noe_emit     ! number of emitted electrons
end type MC_data


!==============================================
! Error handling as an object:
type Error_handling
   logical :: Err          ! indicates that an error occured
   logical :: Warn         ! indicates that an warning (not an error) occured
   logical :: Stopsignal   ! flag to stop calculations, even if there was no error
   integer :: Err_Num      ! assign a number to an error
   integer :: File_Num     ! number of the file with error log
   character(300) Err_descript	! describes more details about the error
end type
!==============================================


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 contains
! How to write the log about an error:
subroutine Save_error_details(Err_name, Err_num, Err_data, Err_data2, empty, Warn)
   class(Error_handling) :: Err_name
   integer, intent(in) :: Err_num
   character(*), intent(in) :: Err_data
   character(*), intent(in), optional :: Err_data2
   logical, intent(in), optional :: empty, Warn
   !-----------------------
   integer :: FN
   logical :: no_err

   if (present(empty)) then ! maybe just initialize with no error
      no_err = empty
   else  ! by default, describe error
      no_err = .false.
   endif

   if (no_err) then  ! no error, just initialize
      Err_name%Err = .false.
      if (present(Warn)) then ! warning present
         Err_name%Warn = Warn
      else
         Err_name%Warn = .false.
      endif
      Err_name%Stopsignal = .false.
      Err_name%Err_descript = '' ! start with an empty string
      Err_name%File_Num = 99     ! default file number for error message
   else ! there is an error, describe it
      FN = Err_name%File_Num
      Err_name%Err = .true.
      Err_name%Stopsignal = .false.
      Err_name%Err_Num = Err_num
      Err_name%Err_descript = Err_data
      ! Main error info:
      write(FN, '(a,i2,1x,a)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript))
      ! Additional indo:
      if (present(Err_data2)) then
         write(FN, '(a)') Err_data2
      endif
   endif ! (no_err)

!  Reminder:
!    Error #1: file not found
!    Error #2: file could not be opened
!    Error #3: file could not be read on the line number given
!    Error #4: some problem with databases (EADL, EPDL97, periodic table file)
!    Error #5: inconsistent TB parametrization (only the same type of parametrization is allowed for all within compound)
!    Error #6: diagonalization subroutine with LAPACK failed (uses MKL library)
!    Error #7: some errors in low-energy electrons (probably in tempereature or chem.potential calculation)
!    Error #8: error in optical coefficients (probably in complex Hamiltonian)
!    Error #9: error in conversion between fluence and dose
end subroutine Save_error_details

END MODULE Objects
