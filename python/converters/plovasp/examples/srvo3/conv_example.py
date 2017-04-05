
import numpy as np
from pytriqs.gf import *
#from sumk_dft import SumkDFT
from sumk_dft_tools import SumkDFTTools
from converters.vasp_converter import VaspConverter

np.set_printoptions(suppress=True)

def density_matrix_and_overlap(sk):
    glist = [GfImFreq(indices=inds,beta=1000.0) for bl, inds in sk.gf_struct_sumk[0]]
    sigma = BlockGf(name_list=[bl for bl, inds in sk.gf_struct_sumk[0]], block_list=glist)
    sigma.zero()
    sk.put_Sigma(Sigma_imp=[sigma])

    print "Overlap matrix:"
    dm_blocks = sk.check_projectors()
    print dm_blocks[0].real

    sk.calc_mu(precision=0.001)

    print
    print "Denisty matrix:"
    dm_blocks = sk.density_matrix(method='using_gf', beta=1000.0)
    ntot = 0.0
    for bl in dm_blocks[0]:
        print
        print bl
        print dm_blocks[0][bl].real
        ntot += dm_blocks[0][bl].real.trace()

    print
    print "  Impurity density:", ntot

def density_of_states(sk):
    glist = [GfReFreq(indices=inds, window=(-10.0, 10.0), n_points=2000) for bl, inds in sk.gf_struct_sumk[0]]
    sigma = BlockGf(name_list=[bl for bl, inds in sk.gf_struct_sumk[0]], block_list=glist)
    sigma.zero()
    sk.put_Sigma(Sigma_imp=[sigma])

    print "  Evaluating DOS..."
    sk.dos_wannier_basis(broadening=0.03, with_dc=False)

if __name__ == '__main__':
    conv = VaspConverter('vasp')
    conv.convert_dft_input()

    sk = SumkDFTTools(hdf_file='vasp.h5')

    density_matrix_and_overlap(sk)
#    density_of_states(sk)
