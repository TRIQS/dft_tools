##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2023 by A. Carta
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



import numpy as np
from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.util import compute_DC_from_density
from triqs.gf import  *
from h5 import HDFArchive
from triqs.operators.util import *
import triqs.utility.mpi as mpi


Uval = 5
Jval = 0.3

method_dict = {
        "FLL":{
            "numbering_convention":0,
            "new_convention":"cFLL"
            },
        "AMF":{
            "numbering_convention":2,
            "new_convention":"cAMF"
            },
        "Held":{
            "numbering_convention":1,
            "new_convention":"cHeld"
            },
        }




def test_dc(SK_compat, SK_new, method, method_dict, dens, Uval, Jval, filename):

    dc_no = method_dict[method]["numbering_convention"]
    dc_string = method_dict[method]["new_convention"]
    
    mpi.report("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    mpi.report(f"\n Testing interface {method} \n")
    mpi.report("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

    mpi.report("\nTesting legacy compatibility layer:\n")
    SK_compat.calc_dc(dens_mat=dens[0], U_interact=Uval, J_hund=Jval, use_dc_formula=dc_no )
    mpi.report("Up DC matrix:")
    mpi.report(SK_compat.dc_imp[0]['up'])
    mpi.report(f"Double counting energy = {SK_compat.dc_energ} ")

    mpi.report("\nTesting new dc interface:\n")
    SK_new.calc_dc(dens_mat=dens[0], U_interact=Uval, J_hund=Jval, use_dc_formula=dc_string)
    mpi.report("Up DC matrix:")
    mpi.report(SK_new.dc_imp[0]['up'])
    mpi.report(f"Double counting energy = {SK_new.dc_energ} ")
    
    # Load previously computed DC values from h5 archive 
    R = HDFArchive(f'./{filename}', 'r')
    dc_comp = R[f'DC_{method}_benchmark']['dc_imp']
    en_comp = R[f'DC_{method}_benchmark']['dc_energ']
    del R

    mpi.report(f"\nAsserting comparison for method: {method}")
    assert np.allclose(SK_compat.dc_imp[0]['up'], dc_comp, atol=1e-12), f"Assertion failed comparing legacy Vdc to reference, method: {method} "
    assert np.allclose(SK_compat.dc_energ, en_comp, atol=1e-12), f"Assertion failed comparing legacy energy to reference, method {method} "
    assert np.allclose(SK_new.dc_imp[0]['up'], dc_comp, atol=1e-12), f"Assertion failed comparing Vdc to reference, method: {method} "
    assert np.allclose(SK_new.dc_energ, en_comp, atol=1e-12), f"Assertion failed comparing energy to reference, method: {method} "
    mpi.report("Comparison with stored DC values successfull!\n")
    



# %% 5 orbitals testing

mpi.report("\n############################################")
mpi.report("############################################")
mpi.report(f"\n \n Starting tests for 5 orbitals \n \n")
mpi.report("############################################")
mpi.report("############################################")
dft_filename = "./NiO.ref"
use_blocks = False
SK_compat = SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=use_blocks)
SK_compat.set_mu(13.9)
SK_new    = SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=use_blocks)
SK_new.set_mu(13.9)


icrsh = 0
dens = SK_compat.density_matrix()

with np.printoptions(precision=5):
    for key in dens[0].keys():
        mpi.report(f"{key} channel")
        mpi.report(dens[0][key].real)

N_up = np.trace(dens[0]['up'].real)
N_down = np.trace(dens[0]['down'].real)
N_tot = N_up + N_down

mpi.report(f"{N_up=} ,{N_down=}, {N_tot=}\n")

for method in ["FLL", "AMF", "Held"]: 
    test_dc(SK_compat, SK_new, method, method_dict, dens, Uval, Jval, filename = f"{dft_filename}.h5")

    #in case implementation changes, to write new testing data into archive
    #R = HDFArchive('./NiO.ref.h5', 'a')
    #R.create_group(f'DC_{method}_benchmark')
    #R[f'DC_{method}_benchmark']['dc_imp']= SK_new.dc_imp[0]['up']
    #R[f'DC_{method}_benchmark']['dc_energ']= SK_new.dc_energ
    #del R


# 3 orbital testing

mpi.report("############################################")
mpi.report("############################################")
mpi.report(f"\n \n Starting tests for 3 orbitals \n \n")
mpi.report("############################################")
mpi.report("############################################")
dft_filename = "./SrVO3.ref"
use_blocks = False
SK_compat = SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=use_blocks)
SK_new    = SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=use_blocks)

icrsh = 0
dens = SK_compat.density_matrix()


with np.printoptions(precision=5):
    for key in dens[0].keys():
        mpi.report(f"{key} channel")
        mpi.report(dens[0][key].real)

N_up = np.trace(dens[0]['up'].real)
N_down = np.trace(dens[0]['down'].real)
N_tot = N_up + N_down

mpi.report(f"{N_up=} ,{N_down=}, {N_tot=}\n")

Uval = 5
Jval = 0.3

for method in ["FLL", "AMF", "Held"]: 
    test_dc(SK_compat, SK_new, method, method_dict, dens, Uval, Jval, filename = f"{dft_filename}.h5" )
    
    #in case implementation changes, to write new testing data into archive
    #R = HDFArchive(f'./{dft_filename}.h5', 'a')
    #R.create_group(f'DC_{method}_benchmark')
    #R[f'DC_{method}_benchmark']['dc_imp']= SK_new.dc_imp[0]['up']
    #R[f'DC_{method}_benchmark']['dc_energ']= SK_new.dc_energ
    #del R
