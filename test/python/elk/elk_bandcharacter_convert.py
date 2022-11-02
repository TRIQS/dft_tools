import os
from h5 import *
from triqs.utility.comparison_tests import *
from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi

from triqs_dft_tools.converters import ElkConverter
#get current working directory path
cwd = format(os.getcwd())
#location of test directory
testdir = cwd+'/elk_bandcharacter_convert'
#change to test directory
os.chdir(testdir)

Converter = ElkConverter(filename='SrVO3', repacking=True)
Converter.hdf_file = 'elk_bc_convert.out.h5'
Converter.convert_dft_input()
Converter.dft_band_characters()

if mpi.is_master_node():
    h5diff('elk_bc_convert.out.h5','elk_bc_convert.ref.h5')

#return to cwd
os.chdir(cwd)
