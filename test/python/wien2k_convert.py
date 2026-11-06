################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from h5 import *
from triqs.utility.comparison_tests import *
from triqs.utility.h5diff import h5diff 
import triqs.utility.mpi as mpi

from triqs_dft_tools.converters import Wien2kConverter
import os
import tempfile
import numpy as np

Converter = Wien2kConverter(filename='SrVO3')
Converter.hdf_file = 'wien2k_convert.out.h5'
Converter.convert_dft_input()

Converter.convert_parproj_input()

if mpi.is_master_node():
    h5diff('wien2k_convert.out.h5','wien2k_convert.ref.h5') 

    with tempfile.TemporaryDirectory() as tmp:
        filename = os.path.join(tmp, 'magnetic_band')
        converter = Wien2kConverter(filename=filename)
        with HDFArchive(converter.hdf_file, 'w') as archive:
            archive['dft_input'] = dict(
                SP=1, SO=0, n_shells=1, n_corr_shells=1,
                shells=[dict(dim=1)], corr_shells=[dict(dim=1)], energy_unit=2.0)

        # outband.f writes all real spin blocks before all imaginary blocks.
        with open(converter.band_file, 'w') as output:
            output.write('''1
2 1
1 2 3
4 5 6
0.1 0.2 0.3
2
11 12 21
31 32 41
13 14 22
33 34 42
1 1 GAMMA
''')
        converter.convert_bands_input()

        expected = np.zeros((1, 2, 1, 2, 1, 2), complex)
        expected[0, 0, 0, :, 0, :] = [[11 + 31j, 12 + 32j],
                                      [13 + 33j, 14 + 34j]]
        expected[0, 1, 0, :, 0, 0] = [21 + 41j, 22 + 42j]
        with HDFArchive(converter.hdf_file, 'r') as archive:
            band = archive[converter.bands_subgrp]
            np.testing.assert_array_equal(band['n_orbitals'], [[2, 1]])
            np.testing.assert_array_equal(band['n_parproj'], [2])
            np.testing.assert_array_equal(band['proj_mat_all'], expected)
