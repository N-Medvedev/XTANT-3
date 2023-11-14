![XTANT3_logo_small_new](https://github.com/N-Medvedev/XTANT-3/assets/104917286/4258b260-88ed-42a3-bfc5-5e228591ce27)

[![DOI](https://zenodo.org/badge/490215542.svg)](https://zenodo.org/badge/latestdoi/490215542)

`X-ray-induced Thermal And Nonthermal Transitions`

It is a hybrid code aimed at modeling femtosecond X-ray-induced effects in matter. The code combines the following methods into one model with feed-backs:

> a) Monte Carlo (MC) method for modeling X-ray-induced electron kinetics
 
> b) Boltzmann equation for low-energy electrons: relaxation-time approximation (RTA) for electron-electron scattering; Boltzmann collision integrals (BCI) for nonadiabatic electron-ion (electron-phonon) coupling
 
> c) Transferable tight binding (TB) for tracing electronic structure and interatomic forces
 
> d) Molecular dynamics (MD) for tracing atomic response to modification of the interatomic potential due to electronic excitation
 
> e) Kubo-Greenwood (or Random phase approximation, RPA) for calculation of the optical properties and electronic heat conductivity of the material

## Disclaimer

_This code is work in progress, anything might change without a notice, bugfixes and patches are expected!_

Although we endeavour to ensure that the code XTANT and results delivered are correct, no warranty is given as to its accuracy (for details, see GPL-3.0 license). This code was developed for non-commercial peaceful purposes only, such as research and education.

## Tables with calculated electron-phonon coupling

Two-tempreature model parameters: the tables with the results with calculated electron-ion (electron-phonon) coupling, and electronic heat capacity, at high electronic temepratures can be found here:
https://github.com/N-Medvedev/XTANT-3_coupling_data

## How to cite

The use of the code is at your own risk. Should you chose to use it, please cite the code and/or the manual:
* N. Medvedev (2023). XTANT-3 [Computer software]. https://doi.org/10.5281/zenodo.8392569

* N. Medvedev “_XTANT-3: X-ray-induced Thermal And Nonthermal Transitions in Matter: theory, numerical details, user manual_” (2023) https://doi.org/10.48550/arXiv.2307.03953

Journal citations may be used as follows: 

* N. Medvedev, V. Tkachenko, V. Lipp, Z. Li, B. Ziaja, "_Various damage mechanisms in carbon and silicon materials under femtosecond x-ray irradiation"_, 4open. 1 (2018) 3. https://doi.org/10.1051/fopen/2018003

Should you use electron-phonon coupling in the calculations, the following citation should be included in addition to the above-mentioned one:

* N. Medvedev, I. Milov, "_Electron-phonon coupling in metals at high electronic temperatures_", Phys. Rev. B. 102 (2020) 064302. https://doi.org/10.1103/PhysRevB.102.064302 

In a publication, at least the following parameters should be mentioned for reproducibility of the results:
Material, its initial structure, the number of atoms in the supercell, the initial conditions (atomic and electronic temperatures), an ensemble used, a type of boundary conditions, a scheme of low-energy electron thermalization (and value of thermalization type is nonequilibrium scheme is used), a type of cross sections in Monte Carlo simulation, a type of tight binding parametrization, whether the electron emission was included or not (if yes, whether Coulomb potential for atoms was accounted for and what model for electron emission was used), whether an additional short-range repulsive potential was used (and what type), time step of MD simulation, parameters of the incoming laser pulse (its photon energy, deposited dose, duration).
