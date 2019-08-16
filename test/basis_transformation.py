from pytriqs.utility.comparison_tests import *
from triqs_dft_tools.sumk_dft import *
import numpy as np

def is_diagonal_matrix(M):
    return abs(np.sum(M-np.diag(np.diagonal(M)))) < 1e-10

def call_diagonalize(SK):
    SK.block_structure.transformation = None
    t_sumk_eal = SK.calculate_diagonalization_matrix(prop_to_be_diagonal='eal', calc_in_solver_blocks=False, write_to_blockstructure = True)
    SK.block_structure.transformation = None
    t_solver_eal = SK.calculate_diagonalization_matrix(prop_to_be_diagonal='eal', calc_in_solver_blocks=True, write_to_blockstructure = True)
    SK.block_structure.transformation = None
    t_solver_dm = SK.calculate_diagonalization_matrix(prop_to_be_diagonal='dm', calc_in_solver_blocks=False, write_to_blockstructure = True)
    SK.block_structure.transformation = None
    t_sumk_dm = SK.calculate_diagonalization_matrix(prop_to_be_diagonal='dm', calc_in_solver_blocks=True, write_to_blockstructure = True)
    SK.block_structure.transformation = None
    return t_sumk_eal, t_solver_eal, t_sumk_dm, t_solver_dm

SK = SumkDFT(hdf_file = 'SrVO3.h5', use_dft_blocks=True)

# only eal and dm are allowed
SK.block_structure.transformation = None
assert not SK.calculate_diagonalization_matrix(prop_to_be_diagonal='test')

# check for shell index
assert not SK.calculate_diagonalization_matrix(shells = [15])

# calling the function twice leads to block_structure.transformation already being set
SK.calculate_diagonalization_matrix()
assert not SK.calculate_diagonalization_matrix()
SK.block_structure.transformation = None

# Check writing to block_structure
SK.calculate_diagonalization_matrix(write_to_blockstructure=False)
assert SK.block_structure.transformation is None
SK.block_structure.transformation = None
SK.calculate_diagonalization_matrix(write_to_blockstructure=True)
assert SK.block_structure.transformation is not None
SK.block_structure.transformation = None

t_sumk_eal, t_solver_eal, t_sumk_dm, t_solver_dm = call_diagonalize(SK)

# All matrices should be identities
for orb in range(SK.n_corr_shells):
    for block in t_solver_eal[orb]:
        assert_arrays_are_close(t_sumk_eal[orb][block],np.identity(3), precision=1e-6)
        assert_arrays_are_close(t_sumk_dm[orb][block],np.identity(3), precision=1e-6)
        assert_arrays_are_close(t_solver_eal[orb][block],np.identity(3), precision=1e-6)
        assert_arrays_are_close(t_solver_dm[orb][block],np.identity(3), precision=1e-6)

SK = SumkDFT(hdf_file = 'w90_convert.ref.h5', use_dft_blocks=True)

t_sumk_eal, t_solver_eal, t_sumk_dm, t_solver_dm = call_diagonalize(SK)

# In this example solver and sumk should be the same
for orb in range(SK.n_corr_shells):
    for block in t_solver_eal[orb]:
        assert_arrays_are_close(t_sumk_eal[orb][block],t_solver_eal[orb][block], precision=1e-6)
        assert_arrays_are_close(t_sumk_dm[orb][block],t_solver_dm[orb][block], precision=1e-6)

# Check if transformations make eal and dm really diagonal
eal = SK.eff_atomic_levels()[0]
for e in eal:
    assert is_diagonal_matrix(np.dot(np.dot(t_solver_eal[0][e], eal[e].conj().T),t_solver_eal[0][e].conj().T))

dm = SK.density_matrix(method='using_point_integration')
for dmi in dm:
    for e in dmi:
        assert is_diagonal_matrix(np.dot(np.dot(t_solver_dm[0][e], dmi[e].conj().T),t_solver_dm[0][e].conj().T))
