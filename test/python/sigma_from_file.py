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
from triqs.gf import *
from triqs.gf.tools import *
from triqs_dft_tools.sumk_dft_tools import *
from triqs.utility.comparison_tests import *
import numpy as np

# Read self energy from hdf file
ar = HDFArchive('SrVO3_Sigma.h5','r')
Sigma_hdf = ar['dmft_transp_input']['Sigma_w']

# Save self energy to txt files
for name, s in Sigma_hdf:
    mesh = np.array([p for p in s.mesh]).reshape(-1,1).real
    re_data = s.data.real.reshape((s.data.shape[0],-1))
    im_data = s.data.imag.reshape((s.data.shape[0],-1))

    mesh_a_data = np.hstack((mesh,re_data,im_data))
    np.savetxt('Sigma_' + name + '.dat', mesh_a_data)

# Read self energy from txt files
SK = SumkDFTTools(hdf_file =  'SrVO3.ref.h5', use_dft_blocks = True)

# the order in the orig SrVO3 file is not assured, hence order it here
a_list = sorted([a for a,al in SK.gf_struct_solver[0].items()])
g_list = [read_gf_from_txt([['Sigma_' + a + '.dat']], a)  for a in a_list]
Sigma_txt = BlockGf(name_list = a_list, block_list = g_list, make_copies=False)

SK.set_Sigma([Sigma_txt])

SK.hdf_file = 'sigma_from_file.out.h5'
SK.save(['Sigma_imp_w'])

assert_block_gfs_are_close(Sigma_txt, Sigma_hdf)
