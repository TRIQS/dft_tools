import unittest
import numpy as np
import sys
sys.path.insert(1, '../python/converters/')
from triqs_dft_tools.converters.wannier90 import Wannier90Converter
from triqs_dft_tools import SumkDFT

class test_w90_conv(unittest.TestCase):

    def test_hopping(self):

        conv1 = Wannier90Converter(seedname='LaVO3-Pnma')
        conv1.convert_dft_input()
        SK1 = SumkDFT(hdf_file='LaVO3-Pnma.h5')

        conv2 = Wannier90Converter(seedname='LaVO3-Pnma_ef')
        conv2.convert_dft_input()
        SK2 = SumkDFT(hdf_file='LaVO3-Pnma_ef.h5')

        for ik in range(SK1.n_k):
            self.assertTrue(np.all(SK1.hopping[ik,0] - conv2.fermi_energy*np.identity(SK1.n_orbitals[ik][0]) - SK2.hopping[ik,0] < 1e-12))

if __name__ == '__main__':
    unittest.main()

