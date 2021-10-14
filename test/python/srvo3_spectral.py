from triqs_dft_tools.sumk_dft_tools import *
from triqs.utility.h5diff import h5diff 
from h5 import HDFArchive

beta = 40

SK = SumkDFTTools(hdf_file='SrVO3_spectral.h5', use_dft_blocks=True)

if mpi.is_master_node():
    with HDFArchive('SrVO3_Sigma.h5', 'a') as ar:
        Sigma = ar['dmft_transp_input']['Sigma_w']
        SK.chemical_potential = ar['dmft_transp_input']['chemical_potential']
        SK.dc_imp = ar['dmft_transp_input']['dc_imp']

Sigma = mpi.bcast(Sigma)
SK.chemical_potential = mpi.bcast(SK.chemical_potential)
SK.dc_imp = mpi.bcast(SK.dc_imp)
SK.set_Sigma([Sigma])

dos_wannier = SK.dos_wannier_basis(broadening=0.01, with_Sigma=True, with_dc=True, save_to_file=False)
dos_parproj = SK.dos_parproj_basis(broadening=0.01, with_Sigma=True, with_dc=True, save_to_file=False)
spaghetti = SK.spaghettis(broadening=0.01, plot_shift=0.0, plot_range=(-1,1), ishell=None, save_to_file=False)

if mpi.is_master_node():

    # with HDFArchive('srvo3_spectral.ref.h5', 'a') as ar:
    #     ar['dos_wannier'] = dos_wannier 
    #     ar['dos_parproj'] = dos_parproj
    #     ar['spaghetti'] = spaghetti

    with HDFArchive('srvo3_spectral.out.h5', 'a') as ar:
        ar['dos_wannier'] = dos_wannier 
        ar['dos_parproj'] = dos_parproj
        ar['spaghetti'] = spaghetti

    h5diff('srvo3_spectral.out.h5', 'srvo3_spectral.ref.h5')
