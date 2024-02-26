
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

from h5 import HDFArchive
from triqs_dft_tools.sumk_dft_tools import SumkDFTTools
import triqs.utility.mpi as mpi
from triqs.utility.comparison_tests import *
from triqs.utility.h5diff import h5diff
import numpy as np
SK = SumkDFTTools(hdf_file = 'SrVO3.ref.h5')

dm = SK.density_matrix(method = 'using_gf', transform_to_solver_blocks=False, with_Sigma=False)
dm_pc = SK.partial_charges(with_Sigma=False, with_dc=False)
dm_pi = SK.density_matrix_using_point_integration()

for key, value in dm[0].items():
    assert (np.allclose(value, dm_pi[0][key], atol=1e-6, rtol=1e-6))

with HDFArchive('sumkdft_basic.out.h5','w') as ar:
    ar['dm'] = dm
    ar['dm_pc'] = dm_pc

if mpi.is_master_node():
    h5diff('sumkdft_basic.out.h5','sumkdft_basic.ref.h5')
