
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

"""
SOC SrVO3 tests: one correlated site and uncorrelated orbitals, one impurity
"""

from triqs_dft_tools.converters import Wannier90Converter
from triqs.utility import h5diff
from triqs.utility import mpi

subfolder = 'w90_convert/'

# bloch_basis False
seedname = subfolder+'SrVO3_soc'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'.out.h5',
                               rot_mat_type='wannier', bloch_basis=False)
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'.out.h5', seedname+'.ref.h5')
