! 000000000000000000000000000000000000000000000000000000000000
! This file is part of XTANT
!
! Copyright (C) 2022 Nikita Medvedev
!
! XTANT is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Although we endeavour to ensure that the code XTANT and results delivered are correct,
! no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions.
! We shall not be liable for any damage arising from the use of this code or its parts
! or any results produced with it, or from any action or decision taken
! as a result of using this code or any related material.
!
! This code is distributed as is for non-commercial peaceful purposes only,
! such as research and education. It is explicitly prohibited to use the code,
! its parts, its results or any related material for military-related and other than peaceful purposes.
!
! By using this code or its materials, you agree with these terms and conditions.
!
! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with 3TB hamiltonian: https://github.com/usnistgov/ThreeBodyTB.jl

MODULE TB_3TB

use Universal_constants
use TB_Koster_Slater
use Objects
use Little_subroutines, only : linear_interpolation, Fermi_function, d_Fermi_function, Find_in_array_monoton
use Electron_tools, only : find_band_gap
use TB_NRL, only : test_nonorthogonal_solution, test_orthogonalization_r, test_orthogonalization_c, Loewdin_Orthogonalization, Loewdin_Orthogonalization_c
use Algebra_tools, only : mkl_matrix_mult, sym_diagonalize, Reciproc, check_hermiticity
use Atomic_tools, only : Reciproc_rel_to_abs


implicit none

 contains




END MODULE TB_3TB
