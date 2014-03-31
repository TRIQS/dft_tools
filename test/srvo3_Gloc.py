

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
from pytriqs.applications.dft.sumk_lda import *
from pytriqs.applications.dft.converters.wien2k_converter import *

#=====================================================
#Basic input parameters:
LDAFilename = 'SrVO3'
U = 4.0
J = 0.6
Beta = 40
DC_type = 1                      # DC type: 0 FLL, 1 Held, 2 AMF
useBlocs = True                  # use bloc structure from LDA input
useMatrix = False                # True: Slater parameters, False: Kanamori parameters U+2J, U, U-J
use_spinflip = False             # use the full rotational invariant interaction?
#=====================================================

#U=U-2*J

HDFfilename = LDAFilename+'.h5'

# Init the SumK class
SK=SumkLDA(hdf_file='SrVO3.h5',use_lda_blocks=True)


Norb = SK.corr_shells[0][3]
l = SK.corr_shells[0][2]


from pytriqs.applications.dft.solver_multiband import *

S=SolverMultiBand(beta=Beta,n_orb=Norb,gf_struct=SK.gf_struct_solver[0],map=SK.map[0])

SK.put_Sigma([S.Sigma])
Gloc=SK.extract_G_loc()

ar = HDFArchive('srvo3_Gloc.output.h5','w')
ar['Gloc'] = Gloc[0]
del ar
