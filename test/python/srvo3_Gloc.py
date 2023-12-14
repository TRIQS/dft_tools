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
from triqs.gf import MeshImFreq, MeshDLRImFreq
from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.converters.wien2k import *
from triqs.operators.util import set_operator_structure
from triqs.utility.comparison_tests import *
from triqs.utility.h5diff import h5diff

import time

# Basic input parameters
beta = 40
n_iw = 1025

# classic full Matsubara mesh
mpi.report(f"{'#'*12}\nregular Matsubara mesh test\n")

# Init the SumK class (reference data with n_iw=1025)
iw_mesh = MeshImFreq(n_iw=n_iw,beta=beta, statistic='Fermion')
SK=SumkDFT(hdf_file='SrVO3.ref.h5',mesh=iw_mesh,use_dft_blocks=True)

num_orbitals = SK.corr_shells[0]['dim']
l = SK.corr_shells[0]['l']
spin_names = ['down','up']
orb_hybridized = False

gf_struct = set_operator_structure(spin_names,num_orbitals,orb_hybridized)
glist = [ Gf(target_shape=(bl_size,bl_size),mesh=iw_mesh) for bl, bl_size in gf_struct]
Sigma_iw = BlockGf(name_list = [bl for bl, bl_size in gf_struct], block_list = glist, make_copies = False)

SK.set_Sigma([Sigma_iw])

if mpi.is_master_node():
    start_time = time.time()

Gloc = SK.extract_G_loc()

if mpi.is_master_node():
    mpi.report(f'extract_G_loc time: {(time.time()-start_time)*1000:.1f} msec')

if mpi.is_master_node():
    with HDFArchive('srvo3_Gloc.out.h5','w') as ar:
        ar['Gloc'] = Gloc[0]

if mpi.is_master_node():
    h5diff("srvo3_Gloc.out.h5","srvo3_Gloc.ref.h5")

mpi.report(f"{'#'*12}\n")


# DLR Matsubara mesh
mpi.report(f"{'#'*12}\nDLR Matsubara mesh test\n")

dlr_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=10, eps=1e-10)
SK=SumkDFT(hdf_file='SrVO3.ref.h5',mesh=dlr_mesh,use_dft_blocks=True)

glist_dlr = [ Gf(target_shape=(bl_size,bl_size),mesh=dlr_mesh) for bl, bl_size in gf_struct]
Sigma_dlr = BlockGf(name_list = [bl for bl, bl_size in gf_struct], block_list = glist_dlr, make_copies = False)
SK.set_Sigma([Sigma_dlr])

if mpi.is_master_node():
    start_time = time.time()

Gloc_dlr_iw = SK.extract_G_loc()

if mpi.is_master_node():
    mpi.report(f'extract_G_loc time: {(time.time()-start_time)*1000:.1f} msec')

    with HDFArchive('srvo3_Gloc.out.h5','a') as ar:
        ar['Gloc_dlr'] = make_gf_imfreq(make_gf_dlr(Gloc_dlr_iw[0]),n_iw=n_iw)
    # get full Giw and compare
    Gloc_iw_full = make_gf_imfreq(make_gf_dlr(Gloc_dlr_iw[0]),n_iw=n_iw)
    assert_block_gfs_are_close(Gloc[0], Gloc_iw_full)

mpi.report(f"{'#'*12}\n")
