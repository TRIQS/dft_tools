################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
# Copyright (c) 2022-2023 Simons Foundation
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
# Authors: M. Aichhorn, S. Beck, A. Hampel, L. Pourovskii, V. Vildosola
################################################################################

from numpy import *
from h5 import HDFArchive
from triqs_dft_tools.converters.wien2k import *
from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.sumk_dft_tools import *
from triqs_dft_tools.sumk_dft_transport import transport_distribution, init_spectroscopy, conductivity_and_seebeck, write_output_to_hdf
from triqs.utility.comparison_tests import *
from triqs.utility import h5diff

beta = 40

Converter = Wien2kConverter(filename='SrVO3', repacking=True)
Converter.convert_dft_input()
Converter.convert_transport_input()

with HDFArchive('SrVO3_Sigma_transport.h5', 'r') as ar:
    Sigma = ar['dmft_transp_input']['Sigma_w']
    SK = SumkDFTTools(hdf_file='SrVO3.ref.h5', mesh=Sigma.mesh, use_dft_blocks=True)
    SK.set_Sigma([Sigma])
    SK.chemical_potential = ar['dmft_transp_input']['chemical_potential']
    SK.dc_imp = ar['dmft_transp_input']['dc_imp']

SK = init_spectroscopy(SK, code='wien2k')
Gamma_w, omega, Om_mesh = transport_distribution(SK, directions=['xx'], broadening=0.0, energy_window=[-0.3,0.3],
                                                 Om_mesh=[0.00, 0.02], beta=beta, with_Sigma=True, code='wien2k')

# SK.save(['Gamma_w','Om_meshr','omega','directions'])
# SK.load(['Gamma_w','Om_meshr','omega','directions'])
optic_cond, seebeck, kappa = conductivity_and_seebeck(Gamma_w, omega, Om_mesh, SK.SP, ['xx'], beta=beta)
output_dict = {'seebeck': seebeck, 'optic_cond': optic_cond, 'kappa': kappa}

# comparison of the output transport data
if mpi.is_master_node():
    write_output_to_hdf(SK, output_dict, 'transp_output')
    out = HDFArchive('SrVO3.ref.h5','r')
    ref = HDFArchive('srvo3_transp.ref.h5', 'r')
    h5diff.compare('', out['transp_output'], ref['transp_output'], 0, 1e-8)