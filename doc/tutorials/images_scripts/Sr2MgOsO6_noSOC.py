# Conversion:
from triqs_dft_tools.converters.wien2k_converter import *
Converter = Wien2kConverter(filename = "Sr2MgOsO6_noSOC")
Converter.convert_dft_input()

import numpy
numpy.set_printoptions(precision=3, suppress=True)

# Set up SK class:
from triqs_dft_tools.sumk_dft import *
SK = SumkDFT(hdf_file='Sr2MgOsO6.h5',use_dft_blocks=True)

eal = SK.eff_atomic_levels()

mat = SK.calculate_diagonalization_matrix(prop_to_be_diagonal='eal')




