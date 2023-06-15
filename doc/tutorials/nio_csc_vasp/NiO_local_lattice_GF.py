from itertools import *
import numpy as np
import triqs.utility.mpi as mpi
from h5 import *
from triqs.gf import *
from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.sumk_dft_tools import *
from triqs.operators.util.hamiltonians import *
from triqs.operators.util.U_matrix import *
from triqs_cthyb import *
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

filename = 'vasp'
beta = 5.0
mesh = MeshImFreq(beta=beta, S='Fermion', n_iw=1000)

SK = SumkDFT(hdf_file = filename+'.h5', use_dft_blocks = False, mesh=mesh)


# We analyze the block structure of the Hamiltonian
Sigma = SK.block_structure.create_gf(mesh=mesh)

SK.put_Sigma([Sigma])


# Setup CTQMC Solver
n_orb = SK.corr_shells[0]['dim']
spin_names = ['up','down']

# Print some information on the master node
iteration_offset = 0
Sigma_iw = 0
block_structure = None
if mpi.is_master_node():
    ar = HDFArchive(filename+'.h5','a')
    if not 'DMFT_results' in ar: ar.create_group('DMFT_results')
    if not 'Iterations' in ar['DMFT_results']: ar['DMFT_results'].create_group('Iterations')
    if 'DMFT_input' in ar:
        block_structure = ar['DMFT_input']['sumk_block_structure']
    if 'iteration_count' in ar['DMFT_results']:
        iteration_offset = ar['DMFT_results']['iteration_count']+1
        print(('offset',iteration_offset))
        Sigma_iw = ar['DMFT_results']['Iterations']['Sigma_it'+str(iteration_offset-1)]
        SK.dc_imp = ar['DMFT_results']['Iterations']['dc_imp'+str(iteration_offset-1)]
        SK.dc_energ = ar['DMFT_results']['Iterations']['dc_energ'+str(iteration_offset-1)]
        SK.chemical_potential = ar['DMFT_results']['Iterations']['chemical_potential'+str(iteration_offset-1)].real

block_structure = mpi.bcast(block_structure)
iteration_offset = mpi.bcast(iteration_offset)
Sigma_iw = mpi.bcast(Sigma_iw)
SK.dc_imp = mpi.bcast(SK.dc_imp)
SK.dc_energ = mpi.bcast(SK.dc_energ)
SK.chemical_potential = mpi.bcast(SK.chemical_potential)

if block_structure:
    SK.block_structure = block_structure
else:
    G = SK.extract_G_loc()
    SK.analyse_block_structure_from_gf(G, threshold = 1e-3)

SK.put_Sigma(Sigma_imp = [Sigma_iw])

ikarray = numpy.array(list(range(SK.n_k)))

gf_csc = Gf(mesh=SK.mesh, target_shape=(SK.proj_mat_csc.shape[2],SK.proj_mat_csc.shape[2]))
G_latt_orb = BlockGf(name_list=['up','down'], block_list=[gf_csc, gf_csc], make_copies=True)

for ik in mpi.slice_array(ikarray):
    G_latt_KS = SK.lattice_gf(ik=ik)*SK.bz_weights[ik]
    for bname, gf in G_latt_orb:
        gf += SK.downfold(ik, 0, bname, G_latt_KS[bname], gf, shells='csc', ir=None)

G_latt_orb << mpi.all_reduce(G_latt_orb)

mpi.barrier()

if mpi.is_master_node():
    ar['DMFT_results']['Iterations']['G_latt_orb_it'+str(iteration_offset-1)] = G_latt_orb
if mpi.is_master_node(): del ar
