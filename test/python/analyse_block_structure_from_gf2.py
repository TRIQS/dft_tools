from triqs.gf import *
from triqs_dft_tools.sumk_dft import SumkDFT
import numpy as np
from triqs.utility.comparison_tests import assert_block_gfs_are_close

# here we test the SK.analyse_block_structure_from_gf function
# with GfReFreq, GfReTime


# helper function to get random Hermitian matrix
def get_random_hermitian(dim):
    herm = np.random.rand(dim,dim)+1.0j*np.random.rand(dim,dim)
    herm = herm + herm.conjugate().transpose()
    return herm

# helper function to get random unitary matrix
def get_random_transformation(dim):
    herm = get_random_hermitian(dim)
    T = expm(1.0j*herm)
    return T

# construct a random block-diagonal Hloc
Hloc = np.zeros((10,10), dtype=np.complex_)
# the Hloc of the first three 2x2 blocks is equal
Hloc0 = get_random_hermitian(2)
Hloc[:2,:2] = Hloc0
Hloc[2:4,2:4] = Hloc0
Hloc[4:6,4:6] = Hloc0
# the Hloc of the last two 2x2 blocks is equal
Hloc1 = get_random_hermitian(2)
Hloc[6:8,6:8] = Hloc1
Hloc[8:,8:] = Hloc1
# construct the hybridization delta
# this is equal for all 2x2 blocks
V = get_random_hermitian(2) # the hopping elements from impurity to bath
b1 = np.random.rand() # the bath energy of the first bath level
b2 = np.random.rand() # the bath energy of the second bath level
delta = GfReFreq(window=(-10,10), indices=list(range(2)), n_points=1001)
delta[0,0] << (V[0,0]*V[0,0].conjugate()*inverse(Omega-b1)+V[0,1]*V[0,1].conjugate()*inverse(Omega-b2+0.02j))/2.0
delta[0,1] << (V[0,0]*V[1,0].conjugate()*inverse(Omega-b1)+V[0,1]*V[1,1].conjugate()*inverse(Omega-b2+0.02j))/2.0
delta[1,0] << (V[1,0]*V[0,0].conjugate()*inverse(Omega-b1)+V[1,1]*V[0,1].conjugate()*inverse(Omega-b2+0.02j))/2.0
delta[1,1] << (V[1,0]*V[1,0].conjugate()*inverse(Omega-b1)+V[1,1]*V[1,1].conjugate()*inverse(Omega-b2+0.02j))/2.0
# construct G
G = BlockGf(name_block_generator=[('ud',GfReFreq(window=(-10,10), indices=list(range(10)), n_points=1001))], make_copies=False)
for i in range(0,10,2):
    G['ud'][i:i+2,i:i+2] << inverse(Omega-delta+0.02j)
G['ud'] << inverse(inverse(G['ud']) - Hloc)


SK = SumkDFT(hdf_file = 'SrIrO3_rot.h5', use_dft_blocks=False)
G_new = SK.analyse_block_structure_from_gf([G])
G_new_symm = G_new[0].copy()
SK.symm_deg_gf(G_new_symm, 0)
assert_block_gfs_are_close(G_new[0], G_new_symm)


assert SK.gf_struct_sumk == [[('ud', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])], [('ud', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])]],\
    "wrong gf_struct_sumk"
for i in range(5):
    assert 'ud_{}'.format(i) in SK.gf_struct_solver[0], "missing block"
    assert SK.gf_struct_solver[0]['ud_{}'.format(i)] == list(range(2)), "wrong block size"
for i in range(10):
    assert SK.sumk_to_solver[0]['ud',i] == ('ud_{}'.format(i//2), i%2), "wrong mapping"

assert len(SK.deg_shells[0]) == 2, "wrong number of equivalent groups found"
assert sorted([len(d) for d in SK.deg_shells[0]]) == [2,3], "wrong number of members in the equivalent groups found"
for d in SK.deg_shells[0]:
    if len(d)==2:
        assert 'ud_3' in d, "shell ud_3 missing"
        assert 'ud_4' in d, "shell ud_4 missing"
    if len(d)==3:
        assert 'ud_0' in d, "shell ud_0 missing"
        assert 'ud_1' in d, "shell ud_1 missing"
        assert 'ud_2' in d, "shell ud_2 missing"



def get_delta_from_mesh(mesh):
    w0 = None
    for w in mesh:
        if w0 is None:
            w0 = w
        else:
            return w-w0

Gt = BlockGf(name_block_generator = [(name,
            GfReTime(window=(-np.pi*(len(block.mesh)-1) / (len(block.mesh)*get_delta_from_mesh(block.mesh)), np.pi*(len(block.mesh)-1) / (len(block.mesh)*get_delta_from_mesh(block.mesh))),
                n_points=len(block.mesh),
                indices=block.indices)) for name, block in G], make_copies=False)

known_moments = np.zeros((2,10,10), dtype=np.complex)
known_moments[1,:] = np.eye(10)
Gt['ud'].set_from_fourier(G['ud'], known_moments)

G_new = SK.analyse_block_structure_from_gf([Gt])
G_new_symm = G_new[0].copy()
SK.symm_deg_gf(G_new_symm, 0)
assert_block_gfs_are_close(G_new[0], G_new_symm)

assert SK.gf_struct_sumk == [[('ud', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])], [('ud', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])]],\
    "wrong gf_struct_sumk"
for i in range(5):
    assert 'ud_{}'.format(i) in SK.gf_struct_solver[0], "missing block"
    assert SK.gf_struct_solver[0]['ud_{}'.format(i)] == list(range(2)), "wrong block size"
for i in range(10):
    assert SK.sumk_to_solver[0]['ud',i] == ('ud_{}'.format(i//2), i%2), "wrong mapping"

assert len(SK.deg_shells[0]) == 2, "wrong number of equivalent groups found"
assert sorted([len(d) for d in SK.deg_shells[0]]) == [2,3], "wrong number of members in the equivalent groups found"
for d in SK.deg_shells[0]:
    if len(d)==2:
        assert 'ud_3' in d, "shell ud_3 missing"
        assert 'ud_4' in d, "shell ud_4 missing"
    if len(d)==3:
        assert 'ud_0' in d, "shell ud_0 missing"
        assert 'ud_1' in d, "shell ud_1 missing"
        assert 'ud_2' in d, "shell ud_2 missing"
