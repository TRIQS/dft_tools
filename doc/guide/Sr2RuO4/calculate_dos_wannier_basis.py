from triqs_dft_tools.converters.wien2k import Wien2kConverter
from triqs_dft_tools import SumkDFTTools

filename = 'Sr2RuO4'

conv = Wien2kConverter(filename = filename,hdf_filename=filename+'.h5')
conv.convert_dft_input()

SK = SumkDFTTools(filename+'.h5')
mesh = (-10.0,10.0,500)
SK.dos_wannier_basis(broadening=(mesh[1]-mesh[0])/float(mesh[2]),
                     mesh=mesh,
                     save_to_file=True,
                     with_Sigma=False,
                     with_dc=False)
