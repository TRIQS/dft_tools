
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

from pytriqs.archive import *
from pytriqs.applications.dft.sumk_lda_tools import SumkLDATools


SK = SumkLDATools(hdf_file = 'SrVO3.h5')

dm = SK.density_matrix(method = 'using_gf', beta = 40)
dm_pc = SK.partial_charges(40)

ar = HDFArchive('sumklda_basic.output.h5','w')
ar['dm'] = dm
ar['dm_pc'] = dm_pc
del ar
