#!/bin/python3
import numpy as np

import triqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs.gf import MeshReFreq

from triqs_dft_tools.sumk_dft import SumkDFT
from triqs_dft_tools.sumk_dft_transport import transport_distribution, init_spectroscopy, conductivity_and_seebeck, write_output_to_hdf
from triqs_dft_tools.converters.wannier90 import Wannier90Converter

kx = 5
seedname = 'sro'
w90_params = {'seedname': 'w90_optics/'+seedname, 'nk_optics': [kx, kx, kx]}
n_wf = 3

# ----------------- set-up sum_k -----------------------

Converter = Wannier90Converter(seedname=w90_params['seedname'], rot_mat_type='wannier', bloch_basis=False)
Converter.convert_dft_input()

# add group 'dft_transp_input' and some other stuff in 'dft_misc_input' needed for transport
h5_archive = ''.join([w90_params['seedname'], '.h5'])
if mpi.is_master_node():
    with HDFArchive(h5_archive, 'a') as ar:
        band_window = np.array([(1, n_wf) for k in range(ar['dft_input']['kpts'].shape[0])])
        band_window_optics = np.array([(1, n_wf) for k in range(ar['dft_input']['kpts'].shape[0])])
        for group_name, group_prop in [('band_window', [band_window]), ('rot_symmetries', [np.eye(3)]), ('n_symmetries', 1)]:
            ar['dft_misc_input'].create_group(group_name)
            ar['dft_misc_input'][group_name] = group_prop
        # dft_transp_input
        ar.create_group('dft_transp_input')
        ar['dft_transp_input'].create_group('band_window_optics')
        ar['dft_transp_input']['band_window_optics'] = [band_window_optics]

# create SumK object
window = (-4., 4.)
n_w = 1001
beta = 40
mu = 11.369
w_mesh = MeshReFreq(window=window, n_w=n_w)
sum_k = SumkDFT(h5_archive, use_dft_blocks=False, mesh=w_mesh)

# ----------------- create Sigma / set into sum_k -----------------------

sigma_freq = [sum_k.block_structure.create_gf(ish=iineq, mesh=w_mesh) for iineq in range(sum_k.n_inequiv_shells)]

# Sigma const
eta_sigma = 0.01
for block, gf in sigma_freq[0]:
    gf << - eta_sigma * 1j

sum_k.set_mu(mu)
sum_k.set_Sigma(sigma_freq)

# ----------------- compute optical properties -----------------------

Om_mesh = [0.0, 0.2]
eta_lattice = 0.001
directions = ['xx']

sum_k, cell_volume = init_spectroscopy(sum_k, code='wannier90', w90_params=w90_params)

Gamma_w, omega, Om_mesh = transport_distribution(sum_k, directions=directions,
                                                 Om_mesh=Om_mesh, energy_window=window,
                                                 with_Sigma=True, broadening=eta_lattice,
                                                 beta=beta, code='wannier90', cell_volume=cell_volume)

if mpi.is_master_node():
    optic_cond, seebeck, kappa = conductivity_and_seebeck(Gamma_w, omega, Om_mesh, sum_k.SP, directions, beta=beta, method=None)
