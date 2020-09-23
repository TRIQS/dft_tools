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

filename = 'nio'
SK = SumkDFT(hdf_file = filename+'.h5', use_dft_blocks = False)

beta = 5.0

# We analyze the block structure of the Hamiltonian
Sigma = SK.block_structure.create_gf(beta=beta)

SK.put_Sigma([Sigma])
G = SK.extract_G_loc()
SK.analyse_block_structure_from_gf(G, threshold = 1e-3)


# Setup CTQMC Solver
n_orb = SK.corr_shells[0]['dim']
spin_names = ['up','down']
orb_names = [i for i in range(0,n_orb)]


# Print some information on the master node
iteration_offset = 0
Sigma_iw = 0
if mpi.is_master_node():
    ar = HDFArchive(filename+'.h5','a')
    if not 'DMFT_results' in ar: ar.create_group('DMFT_results')
    if not 'Iterations' in ar['DMFT_results']: ar['DMFT_results'].create_group('Iterations')
    if 'iteration_count' in ar['DMFT_results']: 
        iteration_offset = ar['DMFT_results']['iteration_count']+1
        print(('offset',iteration_offset))
        Sigma_iw = ar['DMFT_results']['Iterations']['Sigma_it'+str(iteration_offset-1)]
        SK.dc_imp = ar['DMFT_results']['Iterations']['dc_imp'+str(iteration_offset-1)]
        SK.dc_energ = ar['DMFT_results']['Iterations']['dc_energ'+str(iteration_offset-1)]
        SK.chemical_potential = ar['DMFT_results']['Iterations']['chemical_potential'+str(iteration_offset-1)].real

iteration_offset = mpi.bcast(iteration_offset)
Sigma_iw = mpi.bcast(Sigma_iw)
SK.dc_imp = mpi.bcast(SK.dc_imp)
SK.dc_energ = mpi.bcast(SK.dc_energ)
SK.chemical_potential = mpi.bcast(SK.chemical_potential)


SK.put_Sigma(Sigma_imp = [Sigma_iw])

ikarray = numpy.array(list(range(SK.n_k)))

# set up the orbitally resolved local lattice greens function:
n_orbs = SK.proj_mat_csc.shape[2]
spn = SK.spin_block_names[SK.SO]
mesh = Sigma_iw.mesh
block_structure = [list(range(n_orbs)) for sp in spn]
gf_struct = [(spn[isp], block_structure[isp])
         for isp in range(SK.n_spin_blocks[SK.SO])]
block_ind_list = [block for block, inner in gf_struct]
glist = lambda: [GfImFreq(indices=inner, mesh=mesh)
    for block, inner in gf_struct]

G_latt_orb = BlockGf(name_list=block_ind_list,
             block_list=glist(), make_copies=False)
G_latt_orb.zero()

for ik in mpi.slice_array(ikarray):
    G_latt_KS = SK.lattice_gf(ik=ik, beta=beta)
    G_latt_KS *= SK.bz_weights[ik]
    for bname, gf in G_latt_orb:
        add_g_ik = gf.copy()
        add_g_ik.zero()
        add_g_ik << SK.downfold(ik, 0, bname, G_latt_KS[bname], gf, shells='csc', ir=None)
        gf << gf + add_g_ik
        
G_latt_orb << mpi.all_reduce(
                mpi.world, G_latt_orb, lambda x, y: x + y)

mpi.barrier()

if mpi.is_master_node(): 
    ar['DMFT_results']['Iterations']['G_latt_orb_it'+str(iteration_offset-1)] = G_latt_orb
if mpi.is_master_node(): del ar
