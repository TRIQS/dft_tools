
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


import os
import sys
from triqs_dft_tools.converters import Wannier90Converter
from triqs_dft_tools import SumkDFT
from triqs.utility import h5diff
from triqs.utility import mpi
from h5 import HDFArchive

def custom_h5diff(filename1, filename2):
    """
    Compares the dft_input section of two h5 archives, omitting the rot_mat.
    This is useful if there are degenerate eigenvalues where the rot_mat ends
    up depending on the diagonalization algorithm.
    """
    with HDFArchive(filename1, 'r') as out, HDFArchive(filename2, 'r') as ref:
        key = 'dft_input'
        a = out[key]
        b = ref[key]
        if sorted(list(a.keys())) != sorted(list(b.keys())):
            h5diff.failures.append('Two archive groups \'%s\' with different keys \n %s \n vs\n %s'%(key,list(a.keys()), list(b.keys())))
        for k in a.keys():
            # We cannot compare the rot_mats here because of degeneracies from TR symmetry
            if k == 'rot_mat':
                continue
            h5diff.compare(key + '/'+ k, a[k], b[k], 2, 1e-6)
    if h5diff.failures:
        print('-'*50, file=sys.stderr )
        print('-'*20 + '  FAILED  ' +  '-'*20, file=sys.stderr)
        print('-'*50, file=sys.stderr)
        for x in h5diff.failures:
            print(x, file=sys.stderr)
            print('-'*50, file=sys.stderr)
        raise RuntimeError('FAILED')

subfolder = 'w90_convert/'

# --------------------- LaVO3 tests, four correlated sites, one impurity ---------------------

# test rot_mat_type='hloc_diag'
seedname = subfolder+'LaVO3-Pbnm'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'_hloc_diag.out.h5',
                               rot_mat_type='hloc_diag')
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_hloc_diag.out.h5', seedname+'_hloc_diag.ref.h5')

# test rot_mat_type='wannier'
seedname = subfolder+'LaVO3-Pnma'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'_wannier.out.h5',
                               rot_mat_type='wannier')
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_wannier.out.h5', seedname+'_wannier.ref.h5')


# test add_lambda
seedname = subfolder+'LaVO3-Pnma'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'_lambda.out.h5',
                               rot_mat_type='wannier', add_lambda=(.2, .2, .2))
converter.convert_dft_input()

if mpi.is_master_node():
    custom_h5diff(seedname+'_lambda.out.h5', seedname+'_lambda.ref.h5')

# --------------------- Collinear SrVO3 tests, two correlated sites (V t2g and V eg)
# and uncorrelated orbitals (O p), two impurities ---------------------

# test svo, bloch_basis False
seedname = subfolder+'SrVO3_col'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'_wannierbasis.out.h5',
                               rot_mat_type='hloc_diag', bloch_basis=False)
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_wannierbasis.out.h5', seedname+'_wannierbasis.ref.h5')

# test svo, bloch_basis True
if mpi.is_master_node():
    os.rename(seedname + '.OUTCAR', subfolder + 'OUTCAR')
    os.rename(seedname + '.LOCPROJ', subfolder + 'LOCPROJ')
mpi.barrier()
try:
    converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'_blochbasis.out.h5',
                                   rot_mat_type='hloc_diag', bloch_basis=True)
    converter.convert_dft_input()
    fermi_energy = converter.fermi_energy
finally:
    if mpi.is_master_node():
        os.rename(subfolder + 'OUTCAR', seedname + '.OUTCAR')
        os.rename(subfolder + 'LOCPROJ', seedname + '.LOCPROJ')

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_blochbasis.out.h5', seedname+'_blochbasis.ref.h5')

# Compares effective atomic levels between wannier and bloch basis
sum_k_wannier = SumkDFT(seedname+'_wannierbasis.out.h5')
eff_levels_wannier = sum_k_wannier.eff_atomic_levels()
sum_k_bloch = SumkDFT(seedname+'_blochbasis.out.h5')
sum_k_bloch.set_mu(-fermi_energy)
eff_levels_bloch = sum_k_bloch.eff_atomic_levels()
h5diff.compare('', eff_levels_wannier, eff_levels_bloch, 0, 1e-5)
if h5diff.failures :
    print('-'*50, file=sys.stderr)
    print('-'*20 + '  FAILED  ' +  '-'*20, file=sys.stderr)
    print('-'*50, file=sys.stderr)
    for x in h5diff.failures:
        print(x, file=sys.stderr)
        print('-'*50, file=sys.stderr)
    raise RuntimeError('FAILED')

# --------------------- Collinear SrVO3 tests, one correlated site and uncorrelated orbitals, one impurity ---------------------
# test SVO, bloch_basis False and add_lambda
seedname = subfolder+'SrVO3_t2g_col'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'_lambda.out.h5',
                               rot_mat_type='wannier', bloch_basis=False, add_lambda=(.3, .3, .3))
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_lambda.out.h5', seedname+'_lambda.ref.h5')

# --------------------- SOC SrVO3 tests, one correlated site and uncorrelated orbitals, one impurity ---------------------
seedname = subfolder+'SrVO3_soc'
converter = Wannier90Converter(seedname=seedname, hdf_filename=seedname+'.out.h5',
                               rot_mat_type='wannier', bloch_basis=False)
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'.out.h5', seedname+'.ref.h5')
