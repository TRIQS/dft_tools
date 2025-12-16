
##########################################################################
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
##########################################################################
"""
Backward-compatible converters module - re-exports from triqs_dftkit
"""

from triqs_dftkit.wien2k import Converter as Wien2kConverter
from triqs_dftkit.hk import Converter as HkConverter
from triqs_dftkit.vasp import Converter as VaspConverter
from triqs_dftkit.wannier90 import Converter as Wannier90Converter
from triqs_dftkit.elk import Converter as ElkConverter

# Re-export plovasp and elktools submodules
from triqs_dftkit.vasp import plovasp
from triqs_dftkit.elk import elktools

__all__ = ['Wien2kConverter', 'HkConverter', 'Wannier90Converter',
           'VaspConverter', 'ElkConverter', 'plovasp', 'elktools']
