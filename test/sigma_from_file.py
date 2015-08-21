from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft_tools import *
from pytriqs.applications.dft.build_sigma_from_txt import *
import numpy as np

# Read self energy from hdf file
ar = HDFArchive('SrVO3_Sigma.h5','r')
Sigma_hdf = ar['dmft_transp_input']['Sigma_w']

# Save self energy to txt files
for name, s in Sigma_hdf:
    mesh = np.array([p for p in s.mesh]).reshape(-1,1).real
    re_data = s.data.real.reshape((s.data.shape[0],-1))
    im_data = s.data.imag.reshape((s.data.shape[0],-1))

    mesh_a_data = np.hstack((mesh,re_data,im_data))
    np.savetxt('Sigma_' + name + '.dat', mesh_a_data)

# Read self energy from txt files
SK = SumkDFTTools(hdf_file =  'SrVO3.h5', use_dft_blocks = True)
Sigma_txt = constr_Sigma_real_axis(SK, 'Sigma', hdf=False, n_om=101, orb=0)
SK.put_Sigma(Sigma_imp = [Sigma_txt])

SK.hdf_file = 'sigma_from_file.output.h5'
SK.save(['Sigma_imp_w'])


if ((Sigma_txt - Sigma_hdf).real < 1e-6) & ((Sigma_txt - Sigma_hdf).imag < 1e-6):
        print 'Conversion: HDF -> TRIQS -> TXT -> TRIQS successful!'
