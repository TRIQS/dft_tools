import numpy as np

from h5 import HDFArchive
from triqs.gf import BlockGf, Gf
from triqs.utility.comparison_tests import assert_arrays_are_close

from triqs_dft_tools.sumk_dft_tools import SumkDFTTools


def assert_spectral_data_is_valid(Akw, pAkw, pAkw_orb, n_k, n_om, dim):
    for sp in ('up', 'down'):
        assert Akw[sp].shape == (n_k, n_om)
        assert pAkw[0][sp].shape == (n_k, n_om)
        assert pAkw_orb[0][sp].shape == (n_k, n_om, dim, dim)

        assert np.all(np.isfinite(Akw[sp]))
        assert np.all(np.isfinite(pAkw[0][sp]))
        assert np.all(np.isfinite(pAkw_orb[0][sp]))
        assert np.min(Akw[sp]) > -1.0e-12
        assert np.min(pAkw[0][sp]) > -1.0e-12

        assert_arrays_are_close(
            pAkw[0][sp],
            pAkw_orb[0][sp].trace(axis1=2, axis2=3),
            precision=1.0e-12,
        )


def make_diagonal_sigma_from_srvo3():
    with HDFArchive('SrVO3_Sigma.h5', 'r') as ar:
        sigma_srvo3 = ar['dmft_transp_input']['Sigma_w']

    block_list = [
        Gf(mesh=sigma_srvo3.mesh, target_shape=(3, 3)),
        Gf(mesh=sigma_srvo3.mesh, target_shape=(3, 3)),
    ]
    sigma = BlockGf(name_list=['up', 'down'], block_list=block_list, make_copies=False)
    sigma.zero()

    for sp in ('up', 'down'):
        for orb in range(3):
            sigma[sp].data[:, orb, orb] = sigma_srvo3[f'{sp}_{orb}'].data[:, 0, 0]

    return sigma


Sigma = make_diagonal_sigma_from_srvo3()
n_om = len(Sigma.mesh)

SK = SumkDFTTools(hdf_file='svo_vasp_bands.test.h5', mesh=Sigma.mesh)

Akw_vasp, pAkw_vasp, pAkw_orb_vasp = SK.spaghettis(
    broadening=0.05,
    mesh=Sigma.mesh,
    with_Sigma=False,
    with_dc=False,
    proj_type='vasp',
    save_to_file=False,
)

Akw_wann, pAkw_wann, pAkw_orb_wann = SK.spaghettis(
    broadening=0.05,
    mesh=Sigma.mesh,
    with_Sigma=False,
    with_dc=False,
    proj_type='wann',
    save_to_file=False,
)

assert_spectral_data_is_valid(Akw_vasp, pAkw_vasp, pAkw_orb_vasp, SK.n_k, n_om, 3)

for sp in ('up', 'down'):
    assert_arrays_are_close(Akw_vasp[sp], Akw_wann[sp], precision=1.0e-12)
    assert_arrays_are_close(pAkw_vasp[0][sp], pAkw_wann[0][sp], precision=1.0e-12)
    assert_arrays_are_close(pAkw_orb_vasp[0][sp], pAkw_orb_wann[0][sp], precision=1.0e-12)

SK.set_Sigma([Sigma], transform_to_sumk_blocks=False)

Akw_sigma, pAkw_sigma, pAkw_orb_sigma = SK.spaghettis(
    with_Sigma=True,
    with_dc=False,
    proj_type='vasp',
    save_to_file=False,
)

assert_spectral_data_is_valid(Akw_sigma, pAkw_sigma, pAkw_orb_sigma, SK.n_k, n_om, 3)

for sp in ('up', 'down'):
    assert np.max(np.abs(Akw_sigma[sp] - Akw_vasp[sp])) > 1.0e-6
    assert np.max(np.abs(pAkw_sigma[0][sp] - pAkw_vasp[0][sp])) > 1.0e-6
