import os
from h5 import *
from triqs.utility.comparison_tests import *
from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi

from triqs_dft_tools.converters import ElkConverter
from triqs_dft_tools.sumk_dft_tools import *
#get current working directory path
cwd = format(os.getcwd())
#location of test directory
testdir = cwd+'/occ_test'
#change to test directory
os.chdir(testdir)

Converter = ElkConverter(filename='SrVO3', repacking=True)
Converter.hdf_file = 'elk_occ_convert.out.h5'
Converter.convert_dft_input()

SK = SumkDFTTools(hdf_file='elk_occ_convert.out.h5', use_dft_blocks=True)
SK.occupations(with_Sigma=False, with_dc=False)

omin = -1.0
omax = 1.0
oN = 3
mesh = MeshReFreq(omin,omax,oN)
dos_occ = SK.density_of_states(broadening=0.01, mesh=mesh, with_Sigma=False, with_dc=False, dosocc=True, save_to_file=False)

if mpi.is_master_node():

    with HDFArchive('elk_occ_convert.ref.h5', 'a') as ar:
        ar['dos_occ'] = dos_occ
        ar['dos_mesh'] = [omin,omax,oN]
    with HDFArchive('elk_occ_convert.out.h5', 'a') as ar:
        ar['dos_occ'] = dos_occ
        ar['dos_mesh'] = [omin,omax,oN]

if mpi.is_master_node():
    h5diff('elk_occ_convert.out.h5','elk_occ_convert.ref.h5')

#return to cwd
os.chdir(cwd)
