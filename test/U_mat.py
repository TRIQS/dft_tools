
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
import numpy
#from pytriqs.applications.dft.U_matrix import Umatrix
from U_matrix import Umatrix

U = Umatrix(U_interact = 2.0, J_hund = 0.5, l=2)

T = numpy.zeros([5,5],numpy.complex_)
sqtwo = 1.0/numpy.sqrt(2.0)
T[0,0] = 1j*sqtwo
T[0,4] = -1j*sqtwo
T[1,1] = -1j*sqtwo
T[1,3] = -1j*sqtwo
T[2,2] = 1.0
T[3,1] = -sqtwo
T[3,3] = sqtwo
T[4,0] = sqtwo
T[4,4] = sqtwo

U(T=T)

U.reduce_matrix()

ar = HDFArchive('U_mat.output.h5')
ar['U'] = U.U
ar['Up'] = U.Up
ar['Ufull'] = U.Ufull
del ar

