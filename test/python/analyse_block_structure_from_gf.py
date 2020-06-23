from triqs.gf import *
from triqs_dft_tools.sumk_dft import SumkDFT
from scipy.linalg import expm
import numpy as np
from triqs.utility.comparison_tests import assert_gfs_are_close, assert_arrays_are_close, assert_block_gfs_are_close
from h5 import *
import itertools

# The full test checks all different possible combinations of conjugated
# blocks. This takes a few minutes. For a quick test, just checking one
# random value suffices.
# (this parameter affects the second test)
full_test = False

#######################################################################
# First test                                                          #
# where we check the analyse_block_structure_from_gf function         #
# for the SrIrO3_rot.h5 file                                          #
#######################################################################

beta = 40
SK = SumkDFT(hdf_file = 'SrIrO3_rot.h5')
Sigma = SK.block_structure.create_gf(beta=beta)
SK.put_Sigma([Sigma])
G = SK.extract_G_loc()

# the original block structure
block_structure1 = SK.block_structure.copy()
G_new = SK.analyse_block_structure_from_gf(G)

# the new block structure
block_structure2 = SK.block_structure.copy()

with HDFArchive('analyse_block_structure_from_gf.out.h5','w') as ar:
    ar['bs1'] = block_structure1
    ar['bs2'] = block_structure2

# check whether the block structure is the same as in the reference
with HDFArchive('analyse_block_structure_from_gf.out.h5','r') as ar,\
     HDFArchive('analyse_block_structure_from_gf.ref.h5','r') as ar2:
    assert ar['bs1'] == ar2['bs1'], 'bs1 not equal'
    a1 = ar['bs2']
    a2 = ar2['bs2']
    assert a1==block_structure2, "writing/reading block structure incorrect"
    # we set the deg_shells to None because the transformation matrices
    # have a phase freedom and will, therefore, not be equal in general
    a1.deg_shells = None
    a2.deg_shells = None
    assert a1==a2, 'bs2 not equal'

# check if deg shells are correct
assert len(SK.deg_shells[0])==1, "wrong number of equivalent groups"

# check if the Green's functions that are found to be equal in the
# routine are indeed equal
for d in SK.deg_shells[0]:
    assert len(d)==2, "wrong number of shells in equivalent group"
    # the convention is that for every degenerate shell, the transformation
    # matrix v and the conjugate bool is saved
    # then,
    # maybe_conjugate1( v1^dagger G1 v1 ) = maybe_conjugate2( v2^dagger G2 v2 )
    # therefore, to test, we calculate
    #   maybe_conjugate( v^dagger G v )
    # for all degenerate shells and check that they are all equal
    normalized_gfs = []
    for key in d:
        normalized_gf = G_new[0][key].copy()
        normalized_gf.from_L_G_R(d[key][0].conjugate().transpose(), G_new[0][key], d[key][0])
        if d[key][1]:
            normalized_gf << normalized_gf.transpose()
        normalized_gfs.append(normalized_gf)
    for i in range(len(normalized_gfs)):
        for j in range(i+1,len(normalized_gfs)):
            assert_arrays_are_close(normalized_gfs[i].data, normalized_gfs[j].data, 1.e-5)

#######################################################################
# Second test                                                         #
# where a Green's function is constructed from a random model         #
# and the analyse_block_structure_from_gf function is tested for that #
# model                                                               #
#######################################################################

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

# we will conjugate the Green's function blocks according to the entries
# of conjugate_values
# for each of the 5 blocks that will be constructed, there is an entry
# True or False that says whether it will be conjugated
if full_test:
    # in the full test we check all combinations
    conjugate_values = list(itertools.product([False, True], repeat=5))
else:
    # in the quick test we check a random combination
    conjugate_values = [np.random.rand(5)>0.5]

for conjugate in conjugate_values:
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
    delta = G[0]['ud'][:2,:2].copy()
    delta[0,0] << (V[0,0]*V[0,0].conjugate()*inverse(Omega-b1)+V[0,1]*V[0,1].conjugate()*inverse(Omega-b2))/2.0
    delta[0,1] << (V[0,0]*V[1,0].conjugate()*inverse(Omega-b1)+V[0,1]*V[1,1].conjugate()*inverse(Omega-b2))/2.0
    delta[1,0] << (V[1,0]*V[0,0].conjugate()*inverse(Omega-b1)+V[1,1]*V[0,1].conjugate()*inverse(Omega-b2))/2.0
    delta[1,1] << (V[1,0]*V[1,0].conjugate()*inverse(Omega-b1)+V[1,1]*V[1,1].conjugate()*inverse(Omega-b2))/2.0
    # construct G
    G[0].zero()
    for i in range(0,10,2):
        G[0]['ud'][i:i+2,i:i+2] << inverse(Omega-delta)
    G[0]['ud'] << inverse(inverse(G[0]['ud']) - Hloc)

    # for testing symm_deg_gf below, we need this
    # we construct it so that for every group of degenerate blocks of G[0], the
    # mean of the blocks of G_noisy is equal to G[0]
    G_noisy = G[0].copy()
    noise1 = np.random.randn(*delta.target_shape)
    G_noisy['ud'][:2,:2].data[:,:,:] += noise1
    G_noisy['ud'][2:4,2:4].data[:,:,:] -= noise1/2.0
    G_noisy['ud'][4:6,4:6].data[:,:,:] -= noise1/2.0
    noise2 = np.random.randn(*delta.target_shape)
    G_noisy['ud'][6:8,6:8].data[:,:,:] += noise2
    G_noisy['ud'][8:,8:].data[:,:,:] -= noise2

    # for testing backward-compatibility in symm_deg_gf, we need the
    # un-transformed Green's functions
    G_pre_transform = G[0].copy()
    G_noisy_pre_transform = G_noisy.copy()

    # transform each block using a random transformation matrix
    for i in range(0,10,2):
        T = get_random_transformation(2)
        G[0]['ud'][i:i+2,i:i+2].from_L_G_R(T, G[0]['ud'][i:i+2,i:i+2], T.conjugate().transpose())
        G_noisy['ud'][i:i+2,i:i+2].from_L_G_R(T, G_noisy['ud'][i:i+2,i:i+2], T.conjugate().transpose())
        # if that block shall be conjugated, go ahead and do it
        if conjugate[i//2]:
            G[0]['ud'][i:i+2,i:i+2] << G[0]['ud'][i:i+2,i:i+2].transpose()
            G_noisy['ud'][i:i+2,i:i+2] << G_noisy['ud'][i:i+2,i:i+2].transpose()

    # analyse the block structure
    G_new = SK.analyse_block_structure_from_gf(G, 1.e-7)

    # transform G_noisy etc. to the new block structure
    G_noisy = SK.block_structure.convert_gf(G_noisy, block_structure1, beta = G_noisy.mesh.beta, space_from='sumk')
    G_pre_transform = SK.block_structure.convert_gf(G_pre_transform, block_structure1, beta = G_noisy.mesh.beta, space_from='sumk')
    G_noisy_pre_transform = SK.block_structure.convert_gf(G_noisy_pre_transform, block_structure1, beta = G_noisy.mesh.beta, space_from='sumk')

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

        # the convention is that for every degenerate shell, the transformation
        # matrix v and the conjugate bool is saved
        # then,
        # maybe_conjugate1( v1^dagger G1 v1 ) = maybe_conjugate2( v2^dagger G2 v2 )
        # therefore, to test, we calculate
        #   maybe_conjugate( v^dagger G v )
        # for all degenerate shells and check that they are all equal
        normalized_gfs = []
        for key in d:
            normalized_gf = G_new[0][key].copy()
            normalized_gf.from_L_G_R(d[key][0].conjugate().transpose(), G_new[0][key], d[key][0])
            if d[key][1]:
                normalized_gf << normalized_gf.transpose()
            normalized_gfs.append(normalized_gf)
        for i in range(len(normalized_gfs)):
            for j in range(i+1,len(normalized_gfs)):
                # here, we use a threshold that is 1 order of magnitude less strict
                # because of numerics
                assert_gfs_are_close(normalized_gfs[i], normalized_gfs[j], 1.e-6)

    # now we check symm_deg_gf
    # symmetrizing the GF has is has to leave it unchanged
    G_new_symm = G_new[0].copy()
    SK.symm_deg_gf(G_new_symm, 0)
    assert_block_gfs_are_close(G_new[0], G_new_symm, 1.e-6)

    # symmetrizing the noisy GF, which was carefully constructed,
    # has to give the same result as G_new[0]
    SK.symm_deg_gf(G_noisy, 0)
    assert_block_gfs_are_close(G_new[0], G_noisy, 1.e-6)

    # check backward compatibility of symm_deg_gf
    # first, construct the old format of the deg shells
    for ish in range(len(SK.deg_shells)):
        for gr in range(len(SK.deg_shells[ish])):
            SK.deg_shells[ish][gr] = list(SK.deg_shells[ish][gr].keys())

    # symmetrizing the GF as is has to leave it unchanged
    G_new_symm << G_pre_transform
    SK.symm_deg_gf(G_new_symm, 0)
    assert_block_gfs_are_close(G_new_symm, G_pre_transform, 1.e-6)

    # symmetrizing the noisy GF pre transform, which was carefully constructed,
    # has to give the same result as G_pre_transform
    SK.symm_deg_gf(G_noisy_pre_transform, 0)
    assert_block_gfs_are_close(G_noisy_pre_transform, G_pre_transform, 1.e-6)
