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
! This module contains subroutines to deal with TB hamiltonian in the GFN-xTB formalism:
! https://github.com/grimme-lab/xtb
! UNFINISHED, UNUSED!


MODULE TB_xTB

use Universal_constants
use TB_Koster_Slater
use Objects
use Little_subroutines, only : Fermi_function, d_Fermi_function
use Electron_tools, only : find_band_gap
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, test_orthogonalization_c, Loewdin_Orthogonalization, &
                  Loewdin_Orthogonalization_c
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity
use Atomic_tools, only : Reciproc_rel_to_abs


implicit none
PRIVATE

public :: Construct_Vij_xTB, construct_TB_H_xTB, get_Erep_s_xTB, identify_xTB_orbitals_per_atom

!-----------------------------------
! Global GNF0-xTB parameters:

integer, parameter :: m_maxElem = 86    ! number of elements from Periodic Table

! Shell coefficients:
real(8), parameter :: m_kss = 2.0d0
real(8), parameter :: m_kpp = 2.4868d0
real(8), parameter :: m_kdd = 2.27d0

! Electronegativity shell coefficients:
real(8), parameter :: m_ENsh_ss = 0.6d0
real(8), parameter :: m_ENsh_pp = -0.1d0
real(8), parameter :: m_ENsh_dd = -0.2d0

real(8), parameter :: m_bEN = 4.0d0

! Short range bond:
real(8), parameter :: m_k_srb = -0.0129d0
real(8), parameter :: m_eta_srb = 3.4847d0
real(8), parameter :: m_gscal = 0.5097d0
real(8), parameter :: m_c1 = -1.78069d0 ! In the xTB code, it has +
real(8), parameter :: m_c2 = 2.11d0

! Dispersion:
real(8), parameter :: m_a1 = 0.8d0
real(8), parameter :: m_a2 = 4.6d0
real(8), parameter :: m_s8 = 2.85d0


!-----------------------------------
! Repulsive parameters:
real(8), parameter :: m_rep_alpha(1:m_maxElem) = [&
      & 2.1885472d0, 2.2714498d0, 0.6634645d0, 0.9267640d0, 1.1164621d0, &
      & 1.2680750d0, 1.6211038d0, 2.1037547d0, 2.2062651d0, 1.9166982d0, &
      & 0.8129781d0, 0.8408742d0, 0.8361156d0, 0.8859465d0, 1.0684151d0, &
      & 1.1882871d0, 1.4429448d0, 1.1993811d0, 0.5700050d0, 0.8345430d0, &
      & 0.6840185d0, 0.7915733d0, 1.0676223d0, 0.9216746d0, 1.1151815d0, &
      & 1.1883881d0, 1.1895339d0, 1.2692713d0, 1.1734165d0, 1.0018764d0, &
      & 1.1597304d0, 1.1708353d0, 1.2085038d0, 1.1161800d0, 1.3193094d0, &
      & 0.7670615d0, 0.6171015d0, 0.8421909d0, 0.6513468d0, 0.6906528d0, &
      & 0.8705783d0, 0.9711021d0, 1.0252504d0, 0.9847071d0, 1.0559061d0, &
      & 1.0645317d0, 0.9139636d0, 0.9095541d0, 0.9965441d0, 1.0676257d0, &
      & 1.0759855d0, 0.8659486d0, 0.9301733d0, 0.8139884d0, 0.5842740d0, &
      & 0.8070627d0, 0.6961124d0, 0.7599095d0, 0.7667071d0, 0.7735047d0, &
      & 0.7803023d0, 0.7870999d0, 0.7938975d0, 0.8006951d0, 0.8074927d0, &
      & 0.8142903d0, 0.8210879d0, 0.8278855d0, 0.8346831d0, 0.8414808d0, &
      & 0.8482784d0, 0.8803684d0, 0.9915500d0, 0.9875716d0, 1.1535600d0, &
      & 1.1418384d0, 1.1434832d0, 1.1783705d0, 1.0591477d0, 0.9794378d0, &
      & 1.2439938d0, 1.0437958d0, 1.1391049d0, 0.9115474d0, 0.9157573d0, &
      & 0.8137168d0]

real(8), parameter :: m_rep_Zeff(1:m_maxElem) = [&
      & 1.2455414d0,  1.3440060d0,  1.1710492d0,  2.9064151d0,  4.4020866d0, &
      & 4.3101011d0,  4.5460146d0,  4.7850603d0,  7.3393960d0,  4.2503997d0, &
      &10.5220970d0,  7.7916659d0, 11.3886282d0, 13.9495563d0, 16.7912135d0, &
      &13.3874290d0, 13.9700526d0, 14.4971987d0, 13.8061512d0, 13.9719788d0, &
      &10.9127447d0, 13.4067871d0, 16.7322903d0, 21.8192969d0, 22.8754319d0, &
      &25.2196212d0, 26.9753662d0, 27.2652026d0, 26.2195102d0, 14.3840374d0, &
      &25.4102208d0, 43.7565690d0, 34.9344472d0, 22.8724870d0, 34.2378269d0, &
      &15.1027639d0, 39.1086736d0, 32.7340796d0, 18.6398784d0, 22.6163764d0, &
      &27.6545601d0, 37.8625561d0, 40.9844265d0, 30.0686254d0, 35.5737255d0, &
      &28.4443233d0, 25.9740558d0, 28.8257081d0, 53.9657064d0, 88.0203443d0, &
      &82.7978295d0, 39.3120212d0, 49.7072042d0, 45.1199137d0, 55.2536842d0, &
      &50.0381164d0, 48.0939804d0, 46.1827790d0, 46.0844595d0, 45.9861400d0, &
      &45.8878205d0, 45.7895010d0, 45.6911815d0, 45.5928620d0, 45.4945424d0, &
      &45.3962229d0, 45.2979034d0, 45.1995839d0, 45.1012644d0, 45.0029449d0, &
      &44.9046254d0, 41.1538255d0, 46.6524574d0, 53.4995959d0, 73.8197012d0, &
      &59.6567627d0, 50.0720023d0, 49.4064531d0, 44.5201114d0, 39.7677937d0, &
      &58.8051943d0,103.0123579d0, 85.5566053d0, 70.6036525d0, 82.8260761d0, &
      &68.9676875d0]


!-----------------------------------
! GAO functions:

! Number of functions
integer, parameter :: m_nf = 15

! Exponents from first row Table I-V.
real(8), parameter :: m_Alpha1(m_nf) = [ &
      & 2.709498091e-1, & ! 1s
      & 1.012151084e-1, & ! 2s
      & 5.296881757e-2, & ! 3s
      & 3.264600274e-2, & ! 4s
      & 2.216912938e-2, & ! 5s
      & 1.759666885e-1, & ! 2p
      & 9.113614253e-2, & ! 3p
      & 5.578350235e-2, & ! 4p
      & 3.769845216e-2, & ! 5p
      & 1.302270363e-1, & ! 3d
      & 7.941656339e-2, & ! 4d
      & 5.352200793e-2, & ! 5d
      & 1.033434062e-1, & ! 4f
      & 6.952785407e-2, & ! 5f
      & 8.565417784e-2]   ! 5g


! Exponents from second row Table I-V.
real(8), parameter :: m_Alpha2(2, m_nf) = reshape([ &
      & 8.518186635e-1, 1.516232927e-1, & ! 1s
      & 1.292278611e-1, 4.908584205e-2, & ! 2s
      & 6.694095822e-1, 5.837135094e-2, & ! 3s
      & 2.441785453e-1, 4.051097664e-2, & ! 4s
      & 1.213425654e-1, 3.133152144e-2, & ! 5s
      & 4.323908358e-1, 1.069139065e-1, & ! 2p
      & 1.458620964e-1, 5.664210742e-2, & ! 3p
      & 6.190052680e-2, 2.648418407e-2, & ! 4p
      & 2.691294191e-1, 3.980805011e-2, & ! 5p
      & 2.777427345e-1, 8.336507714e-2, & ! 3d
      & 1.330958892e-1, 5.272119659e-2, & ! 4d
      & 6.906014388e-2, 3.399457777e-2, & ! 5d
      & 2.006693538e-1, 6.865384900e-2, & ! 4f
      & 1.156094555e-1, 4.778940916e-2, & ! 5f
      & 1.554531559e-1, 5.854079811e-2],& ! 5g
      & shape(m_Alpha2))

! Coefficients from second row Table I-V.
real(8), parameter :: m_Coeff2(2, m_nf) = reshape([ &
      & 4.301284983e-1, 6.789135305e-1, & ! 1s
      & 7.470867124e-1, 2.855980556e-1, & ! 2s
      &-1.529645716e-1, 1.051370110e+0, & ! 3s
      &-3.046656896e-1, 1.146877294e+0, & ! 4s
      &-5.114756049e-1, 1.307377277e+0, & ! 5s
      & 4.522627513e-1, 6.713122642e-1, & ! 2p
      & 5.349653114e-1, 5.299607212e-1, & ! 3p
      & 8.743116767e-1, 1.513640107e-1, & ! 4p
      &-1.034227010e-1, 1.033376378e+0, & ! 5p
      & 4.666137923e-1, 6.644706516e-1, & ! 3d
      & 4.932764167e-1, 5.918727866e-1, & ! 4d
      & 6.539405185e-1, 3.948945302e-1, & ! 5d
      & 4.769346276e-1, 6.587383976e-1, & ! 4f
      & 4.856637346e-1, 6.125980914e-1, & ! 5f
      & 4.848298074e-1, 6.539381621e-1],& ! 5g
      & shape(m_Coeff2))


! Exponents from third row Table I-V.
real(8), parameter :: m_Alpha3(3, m_nf) = reshape([ &
      & 2.227660584e+0, 4.057711562e-1, 1.098175104e-1, & ! 1s
      & 2.581578398e+0, 1.567622104e-1, 6.018332272e-2, & ! 2s
      & 5.641487709e-1, 6.924421391e-2, 3.269529097e-2, & ! 3s
      & 2.267938753e-1, 4.448178019e-2, 2.195294664e-2, & ! 4s
      & 1.080198458e-1, 4.408119382e-2, 2.610811810e-2, & ! 5s
      & 9.192379002e-1, 2.359194503e-1, 8.009805746e-2, & ! 2p
      & 2.692880368e+0, 1.489359592e-1, 5.739585040e-2, & ! 3p
      & 4.859692220e-1, 7.430216918e-2, 3.653340923e-2, & ! 4p
      & 2.127482317e-1, 4.729648620e-2, 2.604865324e-2, & ! 5p
      & 5.229112225e-1, 1.639595876e-1, 6.386630021e-2, & ! 3d
      & 1.777717219e-1, 8.040647350e-2, 3.949855551e-2, & ! 4d
      & 4.913352950e-1, 7.329090601e-2, 3.594209290e-2, & ! 5d
      & 3.483826963e-1, 1.249380537e-1, 5.349995725e-2, & ! 4f
      & 1.649233885e-1, 7.487066646e-2, 3.135787219e-2, & ! 5f
      & 2.545432122e-1, 1.006544376e-1, 4.624463922e-2],& ! 5g
      & shape(m_Alpha3))

! Coefficients from third row Table I-V.
real(8), parameter :: m_Coeff3(3, m_nf) = reshape([ &
      & 1.543289673e-1, 5.353281423e-1, 4.446345422e-1, & ! 1s
      &-5.994474934e-2, 5.960385398e-1, 4.581786291e-1, & ! 2s
      &-1.782577972e-1, 8.612761663e-1, 2.261841969e-1, & ! 3s
      &-3.349048323e-1, 1.056744667e+0, 1.256661680e-1, & ! 4s
      &-6.617401158e-1, 7.467595004e-1, 7.146490945e-1, & ! 5s
      & 1.623948553e-1, 5.661708862e-1, 4.223071752e-1, & ! 2p
      &-1.061945788e-2, 5.218564264e-1, 5.450015143e-1, & ! 3p
      &-6.147823411e-2, 6.604172234e-1, 3.932639495e-1, & ! 4p
      &-1.389529695e-1, 8.076691064e-1, 2.726029342e-1, & ! 5p
      & 1.686596060e-1, 5.847984817e-1, 4.056779523e-1, & ! 3d
      & 2.308552718e-1, 6.042409177e-1, 2.595768926e-1, & ! 4d
      &-2.010175008e-2, 5.899310608e-1, 4.658445960e-1, & ! 5d
      & 1.737856685e-1, 5.973380628e-1, 3.929395614e-1, & ! 4f
      & 1.909729355e-1, 6.146060459e-1, 3.059611271e-1, & ! 5f
      & 1.780980905e-1, 6.063757846e-1, 3.828552923e-1],& ! 5g
      & shape(m_Coeff3))


! Exponents from forth row Table I-V.
real(8), parameter :: m_Alpha4(4, m_nf) = reshape([ &
      & 5.216844534e+0, 9.546182760e-1, & ! 1s
      & 2.652034102e-1, 8.801862774e-2, &
      & 1.161525551e+1, 2.000243111e+0, & ! 2s
      & 1.607280687e-1, 6.125744532e-2, &
      & 1.513265591e+0, 4.262497508e-1, & ! 3s
      & 7.643320863e-2, 3.760545063e-2, &
      & 3.242212833e-1, 1.663217177e-1, & ! 4s
      & 5.081097451e-2, 2.829066600e-2, &
      & 8.602284252e-1, 1.189050200e-1, & ! 5s
      & 3.446076176e-2, 1.974798796e-2, &
      & 1.798260992e+0, 4.662622228e-1, & ! 2p
      & 1.643718620e-1, 6.543927065e-2, &
      & 1.853180239e+0, 1.915075719e-1, & ! 3p
      & 8.655487938e-2, 4.184253862e-2, &
      & 1.492607880e+0, 4.327619272e-1, & ! 4p
      & 7.553156064e-2, 3.706272183e-2, &
      & 3.962838833e-1, 1.838858552e-1, & ! 5p
      & 4.943555157e-2, 2.750222273e-2, &
      & 9.185846715e-1, 2.920461109e-1, & ! 3d
      & 1.187568890e-1, 5.286755896e-2, &
      & 1.995825422e+0, 1.823461280e-1, & ! 4d
      & 8.197240896e-2, 4.000634951e-2, &
      & 4.230617826e-1, 8.293863702e-2, & ! 5d
      & 4.590326388e-2, 2.628744797e-2, &
      & 5.691670217e-1, 2.074585819e-1, & ! 4f
      & 9.298346885e-2, 4.473508853e-2, &
      & 2.017831152e-1, 1.001952178e-1, & ! 5f
      & 5.441006630e-2, 3.037569283e-2, &
      & 3.945205573e-1, 1.588100623e-1, & ! 5g
      & 7.646521729e-1, 3.898703611e-2],&
      & shape(m_Alpha4))

! Coefficients from forth row Table I-V.
real(8), parameter :: m_Coeff4(4, m_nf) = reshape([ &
      & 5.675242080e-2, 2.601413550e-1, & ! 1s
      & 5.328461143e-1, 2.916254405e-1, &
      &-1.198411747e-2,-5.472052539e-2, & ! 2s
      & 5.805587176e-1, 4.770079976e-1, &
      &-3.295496352e-2,-1.724516959e-1, & ! 3s
      & 7.518511194e-1, 3.589627317e-1, &
      &-1.120682822e-1,-2.845426863e-1, & ! 4s
      & 8.909873788e-1, 3.517811205e-1, &
      & 1.103657561e-2,-5.606519023e-1, & ! 5s
      & 1.179429987e+0, 1.734974376e-1, &
      & 5.713170255e-2, 2.857455515e-1, & ! 2p
      & 5.517873105e-1, 2.632314924e-1, &
      &-1.434249391e-2, 2.755177589e-1, & ! 3p
      & 5.846750879e-1, 2.144986514e-1, &
      &-6.035216774e-3,-6.013310874e-2, & ! 4p
      & 6.451518200e-1, 4.117923820e-1, &
      &-1.801459207e-2,-1.360777372e-1, & ! 5p
      & 7.533973719e-1, 3.409304859e-1, &
      & 5.799057705e-2, 3.045581349e-1, & ! 3d
      & 5.601358038e-1, 2.432423313e-1, &
      &-2.816702620e-3, 2.177095871e-1, & ! 4d
      & 6.058047348e-1, 2.717811257e-1, &
      &-2.421626009e-2, 3.937644956e-1, & ! 5d
      & 5.489520286e-1, 1.190436963e-1, &
      & 5.902730589e-2, 3.191828952e-1, & ! 4f
      & 5.639423893e-1, 2.284796537e-1, &
      & 9.174268830e-2, 4.023496947e-1, & ! 5f
      & 4.937432100e-1, 1.254001522e-1, &
      & 6.010484250e-2, 3.309738329e-1, & ! 5g
      & 5.655207585e-1, 2.171122608e-1],&
      & shape(m_Coeff4))


! Exponents from fifth row Table I-V.
real(8), parameter :: m_Alpha5(5, m_nf) = reshape([ &
      & 1.130563696e+1, 2.071728178e+0, 5.786484833e-1, & ! 1s
      & 1.975724573e-1, 7.445271746e-2, &
      & 8.984956862e+0, 1.673710636e+0, 1.944726668e-1, & ! 2s
      & 8.806345634e-2, 4.249068522e-2, &
      & 4.275877914e+0, 1.132409433e+0, 4.016256968e-1, & ! 3s
      & 7.732370620e-2, 3.800708627e-2, &
      & 2.980263783e+0, 3.792228833e-1, 1.789717224e-1, & ! 4s
      & 5.002110360e-2, 2.789361681e-2, &
      & 7.403763257e-1, 1.367990863e-1, 9.135301779e-2, & ! 5s
      & 3.726907315e-2, 2.241490836e-2, &
      & 3.320386533e+0, 8.643257633e-1, 3.079819284e-1, & ! 2p
      & 1.273309895e-1, 5.606243164e-2, &
      & 6.466803859e+0, 1.555914802e+0, 1.955925255e-1, & ! 3p
      & 8.809647701e-2, 4.234835707e-2, &
      & 1.091977298e+0, 3.719985051e-1, 8.590019352e-2, & ! 4p
      & 4.786503860e-2, 2.730479990e-2, &
      & 3.422168934e-1, 1.665099900e-1, 5.443732013e-2, & ! 5p
      & 3.367775277e-2, 2.091949042e-2, &
      & 1.539033958e+0, 4.922090297e-1, 2.029756928e-1, & ! 3d
      & 9.424112917e-2, 4.569058269e-2, &
      & 1.522122079e+0, 2.173041823e-1, 1.084876577e-1, & ! 4d
      & 5.836797641e-2, 3.206682246e-2, &
      & 9.702946470e-1, 3.603270196e-1, 8.668717752e-2, & ! 5d
      & 4.833708379e-2, 2.751899341e-2, &
      & 8.925960415e-1, 3.277589120e-1, 1.492869962e-1, & ! 4f
      & 7.506099109e-2, 3.892475795e-2, &
      & 1.670735676e+0, 2.072477219e-1, 1.024709357e-1, & ! 5f
      & 5.531913898e-2, 3.072866652e-2, &
      & 5.895429375e-1, 2.393343780e-1, 1.172646904e-1, & ! 5g
      & 6.254074479e-2, 3.411243214e-2],&
      & shape(m_Alpha5))

! Coefficients from fifth row Table I-V.
real(8), parameter :: m_Coeff5(5, m_nf) = reshape([ &
      & 2.214055312e-2, 1.135411520e-1, 3.318161484e-1, & ! 1s
      & 4.825700713e-1, 1.935721966e-1, &
      &-1.596349096e-2,-5.685884883e-2, 3.698265599e-1, & ! 2s
      & 5.480512593e-1, 1.472634893e-1, &
      &-3.920358850e-3,-4.168430506e-2,-1.637440990e-1, & ! 3s
      & 7.419373723e-1, 3.724364929e-1, &
      & 1.513948997e-3,-7.316801518e-2,-3.143703799e-1, & ! 4s
      & 9.032615169e-1, 3.294210848e-1, &
      & 1.375523371e-2,-3.097344179e-1,-3.199192259e-1, & ! 5s
      & 1.084547038e+0, 3.345288361e-1, &
      & 2.079051117e-2, 1.235472099e-1, 3.667738986e-1, & ! 2p
      & 4.834930290e-1, 1.653444074e-1, &
      &-2.329023747e-3,-1.357395221e-2, 2.632185383e-1, & ! 3p
      & 5.880427024e-1, 2.242794445e-1, &
      &-1.143929558e-2,-6.322651538e-2, 4.398907721e-1, & ! 4p
      & 5.245859166e-1, 1.017072253e-1, &
      &-3.113958289e-2,-1.374007017e-1, 5.573881018e-1, & ! 5p
      & 4.855428100e-1, 6.605423564e-2, &
      & 2.020869128e-2, 1.321157923e-1, 3.911240346e-1, & ! 3d
      & 4.779609701e-1, 1.463662294e-1, &
      &-3.673711876e-3, 1.167122499e-1, 4.216476416e-1, & ! 4d
      & 4.547673415e-1, 1.037803318e-1, &
      &-3.231527611e-3,-2.434931372e-2, 3.440817054e-1, & ! 5d
      & 5.693674316e-1, 1.511340183e-1, &
      & 1.999839052e-2, 1.395427440e-1, 4.091508237e-1, & ! 4f
      & 4.708252119e-1, 1.328082566e-1, &
      &-7.301193568e-4, 8.414991343e-2, 3.923683153e-1, & ! 5f
      & 5.040033146e-1, 1.328979300e-1, &
      & 1.998085812e-2, 1.460384050e-1, 4.230565459e-1, & ! 5g
      & 4.635699665e-1, 1.226411691e-1],&
      & shape(m_Coeff5))


! Exponents from sixth row Table I-V.
real(8), parameter :: m_Alpha6(6, m_nf) = reshape([ &
      & 2.310303149e+1, 4.235915534e+0, 1.185056519e+0, & ! 1s
      & 4.070988982e-1, 1.580884151e-1, 6.510953954e-2, &
      & 2.768496241e+1, 5.077140627e+0, 1.426786050e+0, & ! 2s
      & 2.040335729e-1, 9.260298399e-2, 4.416183978e-2, &
      & 3.273031938e+0, 9.200611311e-1, 3.593349765e-1, & ! 3s
      & 8.636686991e-2, 4.797373812e-2, 2.724741144e-2, &
      !& 3.232838646e+0, 3.605788802e-1, 1.717905487e-1, & ! 4s
      !& 5.277666487e-2, 3.163400284e-2, 1.874093091e-2, &
      & 1.365346e+00,   4.393213e-01,   1.877069e-01,   & ! 4s (old)
      & 9.360270e-02,   5.052263e-02,   2.809354e-02,   &
      & 1.410128298e+0, 5.077878915e-1, 1.847926858e-1, & ! 5s
      & 1.061070594e-1, 3.669584901e-2, 2.213558430e-2, &
      & 5.868285913e+0, 1.530329631e+0, 5.475665231e-1, & ! 2p
      & 2.288932733e-1, 1.046655969e-1, 4.948220127e-2, &
      & 5.077973607e+0, 1.340786940e+0, 2.248434849e-1, & ! 3p
      & 1.131741848e-1, 6.076408893e-2, 3.315424265e-2, &
      !& 2.389722618e+0, 7.960947826e-1, 3.415541380e-1, & ! 4p
      !& 8.847434525e-2, 4.958248334e-2, 2.816929784e-2, &
      & 1.365346e+00,   4.393213e-01,   1.877069e-01,   & ! 4p (old)
      & 9.360270e-02,   5.052263e-02,   2.809354e-02,   &
      & 3.778623374e+0, 3.499121109e-1, 1.683175469e-1, & ! 5p
      & 5.404070736e-2, 3.328911801e-2, 2.063815019e-2, &
      & 2.488296923e+0, 7.981487853e-1, 3.311327490e-1, & ! 3d
      & 1.559114463e-1, 7.817734732e-2, 4.058484363e-2, &
      & 4.634239420e+0, 1.341648295e+0, 2.209593028e-1, & ! 4d
      & 1.101467943e-1, 5.904190370e-2, 3.232628887e-2, &
      & 8.820520428e-1, 3.410838409e-1, 9.204308840e-2, & ! 5d
      & 5.472831774e-2, 3.391202830e-2, 2.108227374e-2, &
      & 1.357718039e+0, 5.004907278e-1, 2.296565064e-1, & ! 4f
      & 1.173146814e-1, 6.350097171e-2, 3.474556673e-2, &
      & 1.334096840e+0, 2.372312347e-1, 1.269485144e-1, & ! 5f
      & 7.290318381e-2, 4.351355997e-2, 2.598071843e-2, &
      & 8.574668996e-1, 3.497184772e-1, 1.727917060e-1, & ! 5g
      & 9.373643151e-2, 5.340032759e-2, 3.051364464e-2],&
      & shape(m_Alpha6))

! Coefficients from sixth row Table I-V.
real(8), parameter :: m_Coeff6(6, m_nf) = reshape([ &
      & 9.163596280e-3, 4.936149294e-2, 1.685383049e-1, & ! 1s
      & 3.705627997e-1, 4.164915298e-1, 1.303340841e-1, &
      &-4.151277819e-3,-2.067024148e-2,-5.150303337e-2, & ! 2s
      & 3.346271174e-1, 5.621061301e-1, 1.712994697e-1, &
      &-6.775596947e-3,-5.639325779e-2,-1.587856086e-1, & ! 3s
      & 5.534527651e-1, 5.015351020e-1, 7.223633674e-2, &
      !& 1.374817488e-3,-8.666390043e-2,-3.130627309e-1, & ! 4s
      !& 7.812787397e-1, 4.389247988e-1, 2.487178756e-2, &
      & 3.775056e-03,  -5.585965e-02,  -3.192946e-01,   & ! 4s (old)
      &-2.764780e-02,   9.049199e-01,   3.406258e-01,   &
      & 2.695439582e-3, 1.850157487e-2,-9.588628125e-2, & ! 5s
      &-5.200673560e-1, 1.087619490e+0, 3.103964343e-1, &
      & 7.924233646e-3, 5.144104825e-2, 1.898400060e-1, & ! 2p
      & 4.049863191e-1, 4.012362861e-1, 1.051855189e-1, &
      &-3.329929840e-3,-1.419488340e-2, 1.639395770e-1, & ! 3p
      & 4.485358256e-1, 3.908813050e-1, 7.411456232e-2, &
      !&-1.665913575e-3,-1.657464971e-2,-5.958513378e-2, & ! 4p
      !& 4.053115554e-1, 5.433958189e-1, 1.204970491e-1, &
      &-7.052075e-03,  -5.259505e-02,  -3.773450e-02,   & ! 4p (old)
      & 3.874773e-01,   5.791672e-01,   1.221817e-01,   &
      & 1.163246387e-4,-2.920771322e-2,-1.381051233e-1, & ! 5p
      & 5.706134877e-1, 4.768808140e-1, 6.021665516e-2, &
      & 7.283828112e-3, 5.386799363e-2, 2.072139149e-1, & ! 3d
      & 4.266269092e-1, 3.843100204e-1, 8.902827546e-2, &
      &-4.749842876e-4,-3.566777891e-3, 1.108670481e-1, & ! 4d
      & 4.159646930e-1, 4.621672517e-1, 1.081250196e-1, &
      &-4.097311019e-3,-2.508271857e-2, 2.648458555e-1, & ! 5d
      & 5.097437054e-1, 2.654483467e-1, 2.623132212e-2, &
      & 6.930234381e-3, 5.634263745e-2, 2.217065797e-1, & ! 4f
      & 4.411388883e-1, 3.688112625e-1, 7.787514504e-2, &
      &-9.486751531e-4, 4.624275998e-2, 2.373699784e-1, & ! 5f
      & 4.589112231e-1, 3.205010548e-1, 5.077063693e-2, &
      & 6.729778096e-3, 5.874145170e-2, 2.339955227e-1, & ! 5g
      & 4.512983737e-1, 3.552053926e-1, 6.974153145e-2],&
      & shape(m_Coeff6))

real(8), parameter :: m_Alpha6s(6) = [&
      & 5.800292686e-1, 2.718262251e-1, 7.938523262e-2, & ! 6s
      & 4.975088254e-2, 2.983643556e-2, 1.886067216e-2]
real(8), parameter :: m_Coeff6s(6) = [&
      & 4.554359511e-3, 5.286443143e-2,-7.561016358e-1, & ! 6s
      &-2.269803820e-1, 1.332494651e+0, 3.622518293e-1]

real(8), parameter :: m_Alpha6p(6) = [&
      & 6.696537714e-1, 1.395089793e-1, 8.163894960e-2, & ! 6p
      & 4.586329272e-2, 2.961305556e-2, 1.882221321e-2]
real(8), parameter :: m_Coeff6p(6) = [&
      & 2.782723680e-3,-1.282887780e-1,-2.266255943e-1, & ! 6p
      & 4.682259383e-1, 6.752048848e-1, 1.091534212e-1]



 contains


subroutine Construct_Vij_xTB(numpar, TB, Scell, NSC, M_Vij, M_dVij, M_SVij, M_dSVij)
   type(Numerics_param), intent(in), target :: numpar 	! all numerical parameters
   type(TB_H_xTB), dimension(:,:), intent(in), target :: TB	! All parameters of the Hamiltonian of TB
   type(Super_cell), dimension(:), intent(inout), target :: Scell  ! supercell with all the atoms as one object
   integer, intent(in) :: NSC ! number of supercell
   real(8), dimension(:,:,:), allocatable :: M_Vij	! matrix of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dVij	! matrix of derivatives of Overlap functions for Hamiltonian for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_SVij	! matrix of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   real(8), dimension(:,:,:), allocatable :: M_dSVij	! matrix of derivatives of Overlap functions for S-matrix for all pairs of atoms, all orbitals
   !----------------------------------
   real(8) :: x, y, z, r, sx, sy, sz, r1, fcut, d_fcut, Fermi, dFermi
   integer :: i, j, atom_2, ki, N_bs, ihop
   real(8), pointer :: rm
   integer, pointer :: nat, m, KOA1, KOA2

   nat => Scell(NSC)%Na	! number of atoms in the supercell
   ! number of hopping integrals for this basis set in xTB:
   N_bs = identify_xTB_basis_size(numpar%basis_size_ind)  ! below

   if (.not.allocated(M_Vij)) allocate(M_Vij(nat,nat,N_bs))	! each pair of atoms, all  V functions
   if (.not.allocated(M_dVij)) allocate(M_dVij(nat,nat,N_bs))	! each pair of atoms, all  dV functions
   if (.not.allocated(M_SVij)) allocate(M_SVij(nat,nat,N_bs))	! each pair of atoms, all S functions
   if (.not.allocated(M_dSVij)) allocate(M_dSVij(nat,nat,N_bs))	! each pair of atoms, all dS functions


   M_Vij = 0.0d0
   M_dVij = 0.0d0
   M_SVij = 0.0d0
   M_dSVij = 0.0d0


   ! Construct matrix of all the radial functions for each pair of atoms:
!$omp PARALLEL
!$omp do private(j, m, atom_2, i, KOA1, KOA2, r, ihop, Fermi, dFermi)
   AT1:do j = 1,nat	! all atoms
      m => Scell(NSC)%Near_neighbor_size(j)
      AT2:do atom_2 = 1,m ! do only for atoms close to that one
         i = Scell(NSC)%Near_neighbor_list(j,atom_2) ! this is the list of such close atoms

         KOA2 => Scell(NSC)%MDatoms(j)%KOA  ! Correct order, checked by cohesive energy minimum
         KOA1 => Scell(NSC)%MDatoms(i)%KOA
         r = Scell(NSC)%Near_neighbor_dist(j,atom_2,4) ! at this distance, R [A]

         ! All radial functions for Overlap and Hamiltonian:
!          call xTB_overlap(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 1, M_SVij(j,i,1), M_dVij(j,i,1)) ! (s s sigma)

!          M_Vij(j,i,1) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 1)   ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
!             M_Vij(j,i,2) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 2)   ! (s p sigma)
!             M_SVij(j,i,2) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 2) ! (s p sigma)
!             M_Vij(j,i,3) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 4)   ! (p p sigma)
!             M_SVij(j,i,3) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 4) ! (p p sigma)
!             M_Vij(j,i,4) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 5)   ! (p p pi)
!             M_SVij(j,i,4) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 5) ! (p p pi)
         case default    ! sp3d5
            do ihop = 2, 10
!                M_Vij(j,i,ihop) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, ihop)
!                M_SVij(j,i,ihop) = DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, ihop)
            enddo
         endselect

         KOA1 => Scell(NSC)%MDatoms(i)%KOA    ! testing
         KOA2 => Scell(NSC)%MDatoms(j)%KOA
         ! All derivatives of the radial functions and the radial functions for Overlap matrix:
!          M_dVij(j,i,1) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 1)   ! (s s sigma)
!          M_dSVij(j,i,1) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 1) ! (s s sigma)
         select case (numpar%basis_size_ind)
         case (1)    ! sp3
!             M_dVij(j,i,2) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 2)   ! (s p sigma)
!             M_dSVij(j,i,2) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 2) ! (s p sigma)
!             M_dVij(j,i,3) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 4)   ! (p p sigma)
!             M_dSVij(j,i,3) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 4) ! (p p sigma)
!             M_dVij(j,i,4) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, 5)   ! (p p pi)
!             M_dSVij(j,i,4) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, 5) ! (p p pi)
         case default    ! sp3d5
            do ihop = 2, 10
!                M_dVij(j,i,ihop) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Vr, ihop)
!                M_dSVij(j,i,ihop) = d_DFTB_radial_function(r, TB(KOA1,KOA2)%Rr, TB(KOA1,KOA2)%Sr, ihop)
            enddo
         endselect

         ! Add Fermi-like smoothing approach to zero:
         Fermi = Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         dFermi = d_Fermi_function(TB(KOA1,KOA2)%rcut, TB(KOA1,KOA2)%d, r)   ! module "Little_subroutines"
         M_dVij(j,i,:) = M_dVij(j,i,:)*Fermi + M_Vij(j,i,:)*dFermi
         M_dSVij(j,i,:) = M_dSVij(j,i,:)*Fermi + M_SVij(j,i,:)*dFermi
         M_Vij(j,i,:) = M_Vij(j,i,:)*Fermi
         M_SVij(j,i,:) = M_SVij(j,i,:)*Fermi

!          ! Testing:
!          M_dVij = 0.0d0
!          M_dSVij = 0.0d0
!          M_Vij = 0.0d0
!          M_SVij = 0.d00

      enddo AT2
   enddo AT1
!$omp end do
   nullify(m, KOA1, KOA2)	! clean up for each thread
!$omp END PARALLEL

   nullify(rm, nat, m, KOA1, KOA2)	! clean up at the end
end subroutine Construct_Vij_xTB


subroutine construct_TB_H_xTB()
end subroutine construct_TB_H_xTB



subroutine get_Erep_s_xTB()
end subroutine get_Erep_s_xTB



pure function identify_xTB_basis_size(ind) result (N)
   integer, intent(in) :: ind   ! index of the basis set type
   integer :: N ! number of radial overlap functions
   select case (ind)
   case (0)    ! s
      N = 1
   case (1)    ! s s*
      N = 3
   case (2)    ! sp3
      N = 4
   case (3)    ! sp3s*
      N = 7
   case (4)    ! sp3d5
      N = 10
   case (5)    ! sp3d5s*
      N = 14
   endselect
end function identify_xTB_basis_size

pure function identify_xTB_orbitals_per_atom(ind) result(N)
  integer, intent(in) :: ind   ! index of the basis set type
   integer :: N ! number of orbitals per atom
   select case (ind)
   case (0)    ! s
      N = 1
   case (1)    ! s s*
      N = 2
   case (2)    ! sp3
      N = 4
   case (3)    ! sp3s*
      N = 5
   case (4)    ! sp3d6
      N = 10
   case (5)    ! sp3d6s*
      N = 11
   endselect
end function identify_xTB_orbitals_per_atom


END MODULE TB_xTB
