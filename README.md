# XTANT-3
 X-ray-induced Thermal And Nonthermal Transitions
 
It is a hybrid code aimed at modeling femtosecond X-ray-induced effects in matter. The code combines the following methods into one model with feed-backs:

 a) Monte Carlo (MC) method for modeling X-ray-induced electron kinetics
 
 b) tight binding molecular dynamics (TBMD) for tracing atomic responce to modification of the interatomic potential due to electronic excitation
 
 c) Boltzmann equation for low-energy electrons: relaxation-time approximation (RTA) for electron-electron scattering; Boltzmann collision integrals (BCI) for nonadiabatic electron-ion (electron-phonon) coupling
 
 d) random phase approximation (RPA) for calculation of the change of the optical properties of the material

## Disclamer

_This code is work in progress, anything might change without a notice, bugfixes and patches are expected!_

Although we endeavour to ensure that the code XTANT and results delivered are correct, no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of this code or its parts or any results produced with it, or from any action or decision taken as a result of using this code or any related material.

This code is distributed _as is_ for non-commercial peaceful purposes only, such as research and education. It is __explicitly prohibited__ to use the code, its parts, its results or any related material for military-related or other than peaceful purposes.

By using this code or its materials, you agree with these terms and conditions. 

## How to cite

The use of the code is at your own risk. Should you chose to use it, an appropriate citation is mandatory. The manual and the code can be cited as:

* N. Medvedev “XTANT-3: X-ray-induced Thermal And Nonthermal Transitions in Matter: theory, numerical details, user manual” (2023) https://github.com/N-Medvedev/XTANT-3

Journal citations may be used as follows: 

* N. Medvedev, V. Tkachenko, V. Lipp, Z. Li, B. Ziaja, Various damage mechanisms in carbon and silicon materials under femtosecond x-ray irradiation, 4open. 1 (2018) 3. https://doi.org/10.1051/fopen/2018003

Should you use electron-phonon coupling in the calculations, the following citation should be included in addition to the abovementioned one:

* N. Medvedev, I. Milov, Electron-phonon coupling in metals at high electronic temperatures, Phys. Rev. B. 102 (2020) 064302. https://doi.org/10.1103/PhysRevB.102.064302 

In a publication, at least the following parameters should be mentioned for reproducibility of the results:
Material, its initial structure, the number of atoms in the supercell, the initial conditions (atomic and electronic temperatures), an ensemble used, a type of boundary conditions, a scheme of low-energy electron thermalization (and value of thermalization type is nonequilibrium scheme is used), a type of cross sections in Monte Carlo simulation, a type of tight binding parametrization, whether the electron emission was included or not (if yes, whether Coulomb potential for atoms was accounted for and what model for electron emission was used), whether an additional short-range repulsive potential was used (and what type), time step of MD simulation, parameters of the incoming laser pulse (its photon energy, deposited dose, duration).
