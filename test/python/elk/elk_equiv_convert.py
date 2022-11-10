import os
from h5 import *
from triqs.utility.comparison_tests import *
from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi

from triqs_dft_tools.converters import ElkConverter
#get current working directory path
cwd = format(os.getcwd())
#location of test directory
testdir = cwd+'/elk_equiv_convert'
#change to test directory
os.chdir(testdir)

Converter = ElkConverter(filename='Ni3Al', repacking=True)
Converter.hdf_file = 'elk_equiv_convert.out.h5'
Converter.convert_dft_input()

if mpi.is_master_node():
    h5diff('elk_equiv_convert.out.h5','elk_equiv_convert.ref.h5')

#return to cwd
os.chdir(cwd)
