import numpy as np
from triqs_dft_tools.converters.wannier90 import Wannier90Converter
from triqs_dft_tools import SumkDFT

subfolder = 'w90_convert/'
seedname = subfolder+'LaVO3-Pnma'

conv1 = Wannier90Converter(seedname=seedname)
conv1.convert_dft_input()
SK1 = SumkDFT(hdf_file=seedname+'.h5')

conv2 = Wannier90Converter(seedname=seedname+'_ef')
conv2.convert_dft_input()
SK2 = SumkDFT(hdf_file=seedname+'_ef.h5')

for ik in range(SK1.n_k):
    assert np.allclose(SK1.hopping[ik,0] - conv2.fermi_energy*np.identity(SK1.n_orbitals[ik][0]), SK2.hopping[ik,0], atol=1e-12, rtol=0)
