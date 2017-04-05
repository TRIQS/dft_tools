from pytriqs.applications.dft.sumk_dft import *
from pytriqs.utility.h5diff import h5diff
from pytriqs.gf import *
from pytriqs.utility.comparison_tests import assert_block_gfs_are_close
from pytriqs.applications.dft import BlockStructure

SK = SumkDFT('blockstructure.in.h5',use_dft_blocks=True)

original_bs = SK.block_structure

# check pick_gf_struct_solver
pick1 = original_bs.copy()
pick1.pick_gf_struct_solver([{'up_0': [1], 'up_1': [0], 'down_1': [0]}])

# check loading a block_structure from file
SK.block_structure = SK.load(['block_structure'],'mod')[0]
assert SK.block_structure == pick1, 'loading SK block structure from file failed'

# check SumkDFT backward compatibility
sk_pick1 = BlockStructure(gf_struct_sumk = SK.gf_struct_sumk,
                          gf_struct_solver = SK.gf_struct_solver,
                          solver_to_sumk = SK.solver_to_sumk,
                          sumk_to_solver = SK.sumk_to_solver,
                          solver_to_sumk_block = SK.solver_to_sumk_block)
assert sk_pick1 == pick1, 'constructing block structure from SumkDFT properties failed'

# check pick_gf_struct_sumk
pick2 = original_bs.copy()
pick2.pick_gf_struct_sumk([{'up': [1, 2], 'down': [0,1]}])

# check map_gf_struct_solver
mapping = [{ ('down_0', 0):('down', 0),
             ('down_0', 1):('down', 2),
             ('down_1', 0):('down', 1),
             ('up_0', 0)  :('down_1', 0),
             ('up_0', 1)  :('up_0', 0) }]
map1 = original_bs.copy()
map1.map_gf_struct_solver(mapping)

# check create_gf
G1 = original_bs.create_gf(beta=40,n_points=3)
i = 1
for block,gf in G1:
    gf << SemiCircular(i)
    i+=1

# check approximate_as_diagonal
offd = original_bs.copy()
offd.approximate_as_diagonal()

# check map_gf_struct_solver
G2 = map1.convert_gf(G1,original_bs,beta=40,n_points=3,show_warnings=False)

# check full_structure
full = BlockStructure.full_structure([{'up_0': [0, 1], 'up_1': [0], 'down_1': [0], 'down_0': [0, 1]}],None)

# check __eq__
assert full==full, 'equality not correct (equal structures not equal)'
assert pick1==pick1, 'equality not correct (equal structures not equal)'
assert pick1!=pick2, 'equality not correct (different structures not different)'
assert original_bs!=offd, 'equality not correct (different structures not different)'

if mpi.is_master_node():
    with HDFArchive('blockstructure.out.h5','w') as ar:
        ar['original_bs'] = original_bs
        ar['pick1'] = pick1
        ar['pick2'] = pick2
        ar['map1'] = map1
        ar['offd'] = offd
        ar['G1'] = G1
        ar['G2'] = G2
        ar['full'] = full

    # cannot use h5diff because BlockStructure testing is not implemented
    # there (and seems difficult to implement because it would mix triqs
    # and dft_tools)
    with HDFArchive('blockstructure.out.h5','r') as ar,\
         HDFArchive('blockstructure.ref.h5','r') as ar2:
            for k in ar2:
                if isinstance(ar[k],BlockGf):
                    assert_block_gfs_are_close(ar[k],ar2[k],1.e-6)
                else:
                    assert ar[k]==ar2[k], '{} not equal'.format(k)
