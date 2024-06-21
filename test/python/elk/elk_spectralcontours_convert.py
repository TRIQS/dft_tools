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
testdir = cwd+'/elk_spectralcontours_convert'
#change to test directory
os.chdir(testdir)

#default k-mesh
Converter = ElkConverter(filename='SrVO3', repacking=True)
Converter.hdf_file = 'elk_spectralcontours_convert.out.h5'
Converter.convert_dft_input()
Converter.convert_contours_input()

omin = -1.0
omax = 1.0
oN = 3
mesh = MeshReFreq(omin,omax,oN)
SK = SumkDFTTools(hdf_file='elk_spectralcontours_convert.out.h5', use_dft_blocks=True)
fs_elk = SK.spectral_contours(broadening=0.01, mesh=mesh, with_Sigma=False, with_dc=False, FS=True, proj_type='wann', save_to_file=False)
omega_elk = SK.spectral_contours(broadening=0.01, mesh=mesh, with_Sigma=False, with_dc=False, FS=False, proj_type='wann', save_to_file=False)
omega_range_elk = SK.spectral_contours(broadening=0.01, mesh=mesh, plot_range=(-0.5,2), with_Sigma=False, with_dc=False, FS=False, proj_type='wann', save_to_file=False)

#user specified k-mesh - has to be same as used in elk.in
Converter = ElkConverter(filename='SrVO3', repacking=True)
Converter.hdf_file = 'elk_spectralcontours_convert.out.h5'
ngrid=np.array([10,10,1],np.int_)
kgrid=np.array([[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],np.float64)
Converter.convert_contours_input(kgrid=kgrid,ngrid=ngrid)
SK2 = SumkDFTTools(hdf_file='elk_spectralcontours_convert.out.h5', use_dft_blocks=True)
fs_elk_user = SK2.spectral_contours(broadening=0.01, mesh=mesh, with_Sigma=False, with_dc=False, FS=True, proj_type='wann', save_to_file=False)

if mpi.is_master_node():

    #with HDFArchive('elk_spectralcontours_convert.ref.h5', 'a') as ar:
    #    ar['fs_elk'] = fs_elk
    #    ar['fs_elk_user'] = fs_elk_user
    #    ar['omega_elk'] = omega_elk
    #    ar['omega_range_elk'] = omega_range_elk
    #    ar['mesh'] = [omin,omax,oN]
    with HDFArchive('elk_spectralcontours_convert.out.h5', 'a') as ar:
        ar['fs_elk'] = fs_elk
        ar['fs_elk_user'] = fs_elk_user
        ar['omega_elk'] = omega_elk
        ar['omega_range_elk'] = omega_range_elk
        ar['mesh'] = [omin,omax,oN]


if mpi.is_master_node():
    h5diff('elk_spectralcontours_convert.out.h5','elk_spectralcontours_convert.ref.h5')

#return to cwd
os.chdir(cwd)
