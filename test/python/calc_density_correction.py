# Copyright (c) 2013 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
# Copyright (c) 2013 Centre national de la recherche scientifique (CNRS)
# Copyright (c) 2019-2024 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https:#www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: A. Hampel

import shutil
from triqs.gf import Gf, MeshImFreq
from triqs_dft_tools.sumk_dft import SumkDFT

# Simple run w/o error test for all DFT codes to write the density correction to a file
# Comparison with real self energy against refrence data should be done in the future

# define mesh for all calculations
beta = 40
mesh = MeshImFreq(beta, statistic='Fermion', n_iw=1024)

# Wien2k test
sumk = SumkDFT(hdf_file='SrVO3.ref.h5', mesh=mesh)

Sigma_iw = [
    sumk.block_structure.create_gf(ish=iineq, gf_function=Gf, space='solver', mesh=sumk.mesh) for iineq in range(sumk.n_inequiv_shells)
]

for iineq in range(sumk.n_inequiv_shells):
    Sigma_iw[iineq] << 0.1 + 0.0j

sumk.set_Sigma(Sigma_iw)

deltaN, dens = sumk.calc_density_correction(dm_type='wien2k')

####################################
# Elk
sumk = SumkDFT(hdf_file='elk/elk_convert/elk_convert.ref.h5', mesh=mesh)

Sigma_iw = [
    sumk.block_structure.create_gf(ish=iineq, gf_function=Gf, space='solver', mesh=sumk.mesh) for iineq in range(sumk.n_inequiv_shells)
]

for iineq in range(sumk.n_inequiv_shells):
    Sigma_iw[iineq] << 0.1 + 0.0j

sumk.set_Sigma(Sigma_iw)

deltaN, dens = sumk.calc_density_correction(dm_type='elk')

####################################
# Vasp
sumk = SumkDFT(hdf_file='plovasp/converter/lunio3.ref.h5', mesh=mesh)

Sigma_iw = [
    sumk.block_structure.create_gf(ish=iineq, gf_function=Gf, space='solver', mesh=sumk.mesh) for iineq in range(sumk.n_inequiv_shells)
]

for iineq in range(sumk.n_inequiv_shells):
    Sigma_iw[iineq] << 0.1 + 0.0j

sumk.set_Sigma(Sigma_iw)

deltaN, dens, en_corr = sumk.calc_density_correction(dm_type='vasp')

####################################
# QE
shutil.copy('w90_convert/SrVO3_col_blochbasis.ref.h5', 'w90_convert/SrVO3_col_blochbasis.test.h5')
sumk = SumkDFT(hdf_file='w90_convert/SrVO3_col_blochbasis.test.h5', mesh=mesh)

Sigma_iw = [
    sumk.block_structure.create_gf(ish=iineq, gf_function=Gf, space='solver', mesh=sumk.mesh) for iineq in range(sumk.n_inequiv_shells)
]

for iineq in range(sumk.n_inequiv_shells):
    Sigma_iw[iineq] << 0.1 + 0.0j

sumk.set_Sigma(Sigma_iw)
# little hack to speed up the calculation
sumk.n_k = 1

deltaN, dens, en_corr = sumk.calc_density_correction(dm_type='qe', filename='dump.h5')

