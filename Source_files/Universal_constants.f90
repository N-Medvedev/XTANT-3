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
! This module contains universal constants
! The values are taken from: http://physics.nist.gov/constants

MODULE Universal_constants
 implicit none

! Constants:
real(8) :: g_Pi, g_2Pi, g_half_Pi, g_exp, g_e, g_me, g_cvel, g_Mp, g_h, g_kb, g_kb_EV, g_kb_J, g_e0, g_ke, g_mu0, &
g_Ry, g_a0, g_v0, g_alpha, g_P_atm, g_e_m, g_h_MeVs, g_e_ESU, g_me_MeV, g_u_MeV, g_amu, g_re, g_lambda_e, g_SIGMA_0, &
g_MU_B_MeV_T, g_E_M_e, g_E_M_P, g_G, g_g_Earth, g_N_A, g_V_MOLAR, g_LAMBDAT, g_SIGMA_SB, g_G_F, g_M_W, g_M_Z0, &
g_G_S, g_AUENERGY, g_AUACTION, g_AUTIME, g_AUFORCE, g_AUVELOCITY, g_AUMOMENTUM, g_AUEFIELD, g_AUEDIPOLE, &
g_AUMFLUX, g_AUMDIPOLE, g_ASTRONOMICALUNIT, g_NA, g_ms2Afs, g_Afs2ms, g_r0, g_Eulers_gamma, g_sqrt_Pi

complex :: g_CI

! Conversion coefficients:
real(8) :: g_au2A, g_A2au, g_au2ev, g_ev2au, g_au2am, g_am2au, g_au2fs, g_fs2au, g_ev2kc, g_kc2ev, g_au2kc, g_kc2au, g_au2ic, &
g_ic2au, g_ev2ic, g_ic2ev, g_ev2kj, g_kj2ev, g_in2cm, g_ft2m, g_yd2m, g_cm2in, g_m2ft, g_m2yd, g_ev2Ry, g_rad2deg, g_deg2rad

!UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
! Universal constants:
parameter (g_Pi 	= 3.1415926535897932384626433832795d0)	! Pi
parameter (g_2Pi   = 2.0d0*g_Pi)        ! 2*Pi
parameter (g_half_Pi = 0.5d0*g_Pi)      ! Pi/2
parameter (g_sqrt_Pi = sqrt(g_Pi))      ! sqrt(Pi)
parameter (g_exp 	= dexp(1.0d0))      ! e
parameter (g_Eulers_gamma = 0.5772156649015328606d0) ! Euler's gamma
parameter (g_e 		= 1.602176634d-19)	! Electron charge	[Coulomb]
parameter (g_e_ESU	= 4.8032068d-10)	! Electron charge magnitude	[esu]
parameter (g_me	= 9.1093821545d-31)	    ! Electron mass	[kg]
parameter (g_me_MeV	= 0.51099906d0)		! Electron mass	[MeV/c^2]
parameter (g_u_MeV	= 931.49432d0)		! unified atomic mass unit	[MeV/c^2]
parameter (g_amu  = 1.6605390666050e-27)    ! atomic mass unit (Dalton) [kg]
parameter (g_cvel	= 299792458.0d0)	! Light velocity	[m/sec]
parameter (g_Mp	= 1836.1526724780d0*g_me)	! Proton mass		[kg]
parameter (g_h		= 1.05457162853d-34)	! Plank constant	[J*sec]
parameter (g_h_MeVs	= 6.582122d-22)		! Planck constant, reduced	[MeV*s]
parameter (g_kb		= 11604.0d0)		! Boltzmann constant	[K/eV]
parameter (g_kb_EV	= 1.0d0/g_kb)		! Boltzmann constant	[eV/K]
parameter (g_kb_J	= 1.38064852d-23)	! Boltzmann constant	[J/K]
parameter (g_e0		= 8.854187817620d-12)	! Electrical constant	[F/m]
parameter (g_ke     = 1.0d0/(4.0d0*g_Pi*g_e0))  ! Coulomb constant
parameter (g_mu0	= 1.2566370614359d-6)	! Magnetic constant	[H*A^-2]
parameter (g_Ry	= 13.6056981d0)		! Rydberg constant	[eV]
parameter (g_alpha	= g_e*g_e/(g_h*g_cvel*4.0d0*g_Pi*g_e0))	! Fine structure constant // 0.0072973530796448
parameter (g_a0		= 0.5291772085936d0)	! Bohr radius		[A]
parameter (g_r0		= g_alpha*g_alpha*g_a0)	! classical electron radius [A]
parameter (g_re		= 2.81794092d-15)	! classical electron radius	[m]
parameter (g_v0		= sqrt(2.0d0*g_Ry*g_e/g_me))	! Bohr velocity	[m/s]
parameter (g_P_atm	= 101325.0d0)		! Atmospheric pressure	[Pa]
parameter (g_e_m	= g_e/g_me)		! Ratio of electron charge to its mass	[Coulomb/kg]
parameter (g_lambda_e	= 3.86159323d-13)	! electron Compton wavelength	[m]
parameter (g_SIGMA_0	= 0.66524616d0)		! Thomson cross section	[barn]
parameter (g_MU_B_MeV_T	= 5.78838263d-11)	! Bohr magneton		[MeV/T]
parameter (g_E_M_e	= 1758819620d0)		! electron cyclotron freq./field	[C/kg (rad/sT)]
parameter (g_E_M_P	= 95788309d0)		! proton cyclotron freq./field	[C/kg (rad/sT)]
parameter (g_G		= 6.67259d-11)		! gravitational constant	[m^3/kgs^2]
parameter (g_g_Earth	= 9.80665d0)		! standard grav. accel., sea level	[m/s^2]
parameter (g_NA		= 6.0221367d+23)	! Avogadro constant	[1/mole]
parameter (g_V_MOLAR	= 0.0224141d0)		! molar volume, ideal gas at STP	[m^3/mole]
parameter (g_LAMBDAT	= 0.002897756d0)	! Wien displacement law constant	[m K]
parameter (g_SIGMA_SB	= 5.67051d-08)		! Stefan-Boltzmann constant	[W/m^2K^4]
parameter (g_G_F	= 1.16639d-05)		! Fermi coupling constant	[GeV^{-2}]
parameter (g_M_W	= 80.22d0)		! W boson mass	[GeV/c^2]
parameter (g_M_Z0	= 91.187d0)		! Z_0 boson mass	[GeV/c^2]
parameter (g_G_S	= 0.117d0)		! strong coupling constant	[at M_Z]
parameter (g_AUENERGY	= 4.3597482d-18)	! 1 a.u. energy in SI, hartree energy
parameter (g_AUACTION	= 1.05457266d-34)	! 1 a.u. action in SI, Planck constant/(2*PI)
parameter (g_AUTIME	= 2.4188843341d-17)	! 1 a.u. time in SI
parameter (g_AUFORCE	= 8.2387295d-08)	! 1 a.u. force in  SI
parameter (g_AUVELOCITY	= 2187691.42d0)		! 1 a.u. velocity in SI
parameter (g_AUMOMENTUM	= 1.9928534d-34)	! 1 a.u. momentum in SI
parameter (g_AUEFIELD	= 514220820d0)		! 1 a.u. el. field in SI
parameter (g_AUEDIPOLE	= 8.4783579d-30)	! 1 a.u. el. dipole in SI
parameter (g_AUMFLUX	= 235051.808d0)		! 1 a.u. magn. flux in SI
parameter (g_AUMDIPOLE	= 1.85480308d-23)	! 1 a.u. magn. flux in SI
parameter (g_ASTRONOMICALUNIT	= 14959787d0)	! 1 AE length in SI
parameter (g_CI		= dcmplx(0.0d0,1.0d0) )	! complex unity
! Conversion between units:
parameter (g_au2A	= 0.529177249d0)	! [a.u.] -> [Angstroem]
parameter (g_A2au	= 1.8897259885789d0)	! [Angstroem] -> [a.u.]
parameter (g_au2ev	= 27.211396131788d0)	! [a.u.] -> [eV]
parameter (g_ev2au	= 0.036749308824762d0)	! [eV] -> [a.u.]
parameter (g_ev2Ry	= 0.07349861764)		! [eV] -> [Ry]
parameter (g_au2am	= 0.00054857989586762d0)	! [a.u.] -> [atomar mass units]
parameter (g_am2au	= 1822.8885300626d0)	! [atomar mass units] -> [a.u.]
parameter (g_au2fs	= 0.024188843341d0)	! [a.u.] -> [fs]
parameter (g_fs2au	= 41.341373206755d0)	! [fs] -> [a.u.]
parameter (g_ev2kc	= 23.05d0)		! [eV] -> [kcal/mol]
parameter (g_kc2ev	= 0.043364476254697d0)	! [kcal/mol] -> [eV]
parameter (g_au2kc	= 627.50431878767d0)	! [a.u.] -> [kcal/mol]
parameter (g_kc2au	= 0.0015936145299079d0)	! [kcal/mol] -> [a.u.]
parameter (g_au2ic	= 219474.7d0)		! [au] -> [1/cm]
parameter (g_ic2au	= 4.5563338279993d-06)	! [1/cm] -> [au]
parameter (g_ev2ic	= 8065.54353d0)		! [eV] -> [1/cm]
parameter (g_ic2ev	= 0.0001239842047d0)	! [1/cm] -> [eV]
parameter (g_ev2kj	= 96.485d0)		! [eV] -> [kJoule/mol]
parameter (g_kj2ev	= 0.0103643d0)		! [kJoule/mol] -> [eV]
parameter (g_in2cm	= 2.54d0)			! [inch] -> [cm]
parameter (g_ft2m	= 0.3048d0)		! [feet] -> [m]
parameter (g_yd2m	= 0.9144d0)		! [yard] -> [m]
parameter (g_cm2in	= 0.39370078740157d0)	! [cm]-> [inch]
parameter (g_m2ft	= 3.2808398950131d0)	! [m] ->  [feet]
parameter (g_m2yd	= 1.0936132983377d0)	! [m] -> [yard]
parameter (g_ms2Afs = 1.0d-5)                   ! [m/s] -> [A/fs]
parameter (g_Afs2ms = 1.0d5)                    ! [A/fs] -> [m/s]
parameter (g_rad2deg = 180.0d0/g_Pi)    ! [rad] -> [deg]
parameter (g_deg2rad = g_Pi/180.0d0)    ! [deg] -> [rad]
!UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

END MODULE Universal_constants
