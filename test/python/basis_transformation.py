from triqs.utility.comparison_tests import *
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

SK = SumkDFT(hdf_file = 'SrVO3.ref.h5', use_dft_blocks=True)

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


# Test convert_operator
SK = SumkDFT(hdf_file = 'SrVO3.ref.h5', use_dft_blocks=True)
BS = SK.block_structure
from triqs.operators.util import h_int_slater, U_matrix, t2g_submatrix, transform_U_matrix


U3x3 = t2g_submatrix(U_matrix(2, U_int=2, J_hund=0.2, basis='spheric'))

BS.transformation = [{'up':np.eye(3), 'down': np.eye(3)}]
H0 = h_int_slater(spin_names=['up','down'], orb_names=range(3), U_matrix=U3x3, off_diag=False)
H1 = h_int_slater(spin_names=['up','down'], orb_names=range(3), U_matrix=U3x3, off_diag=True)
assert( H0 == BS.convert_operator(H1) )

# Trafo Matrix switching index 1 & 2
BS.transformation = [{'up':np.array([[1,0,0],[0,0,1],[0,1,0]]), 'down': np.array([[1,0,0],[0,0,1],[0,1,0]])}]
H2 = BS.convert_operator(h_int_slater(spin_names=['up','down'], orb_names=[0,2,1], U_matrix=U3x3, off_diag=True))
assert( H0 == H2 )

BS.transformation = [{'up':np.array([[1,0,0],[0,1/np.sqrt(2),1/np.sqrt(2)],[0,1/np.sqrt(2),-1/np.sqrt(2)]]), 'down': np.array([[1,0,0],[0,1/np.sqrt(2),1/np.sqrt(2)],[0,1/np.sqrt(2),-1/np.sqrt(2)]])}]
H3 = BS.convert_operator(h_int_slater(spin_names=['up','down'], orb_names=[0,1,2], U_matrix=U3x3, off_diag=True))
for op in H3:
    for c_op in op[0]:
        assert(BS.gf_struct_solver_dict[0][c_op[1][0]][c_op[1][1]] is not None) # This crashes with a key error if the operator structure is not the solver structure

U_trafod = transform_U_matrix(U3x3, BS.transformation[0]['up'].conjugate()) # The notorious .conjugate()
H4 = h_int_slater(spin_names=['up','down'], orb_names=range(3), U_matrix=U_trafod, map_operator_structure=BS.sumk_to_solver[0])
assert( H4  == H3 ) # check that convert_operator does the same as transform_U_matrix
