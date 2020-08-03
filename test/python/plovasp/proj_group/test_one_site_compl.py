
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import VaspData
from triqs_dft_tools.converters.plovasp.elstruct import ElectronicStructure
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters
from triqs_dft_tools.converters.plovasp.proj_shell import ProjectorShell
from triqs_dft_tools.converters.plovasp.proj_group import ProjectorGroup
from h5 import HDFArchive
import mytest

################################################################################
#
# TestProjectorGroup
#
################################################################################
class TestProjectorGroupCompl(mytest.MyTestCase):
    """
    Class:

    ProjectorGroupCompl(sh_pars, proj_raw)

    Scenarios:
    - **test** that unequal number of bands at different k-points gives error
    - **test** that COMLEMENT=TRUE gives orthonormal projectors
    """
    def setUp(self):
        conf_file = _rpath + 'example.cfg'
        self.pars = ConfigParameters(conf_file)
        self.pars.parse_input()
        vasp_data = VaspData(_rpath + 'one_site/')
        self.el_struct = ElectronicStructure(vasp_data)

        efermi = self.el_struct.efermi
        self.eigvals = self.el_struct.eigvals - efermi

        struct = self.el_struct.structure
        kmesh = self.el_struct.kmesh

        self.proj_sh = ProjectorShell(self.pars.shells[0], vasp_data.plocar.plo, vasp_data.plocar.proj_params, kmesh, struct, 0)


    def test_num_bands(self):
        self.pars.groups[0]['complement'] = True
        err_mess = "At each band the same number"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.proj_gr = ProjectorGroup(self.pars.groups[0], [self.proj_sh], self.eigvals)

    def test_compl(self):
        self.pars.groups[0]['complement'] = True
        self.pars.groups[0]['bands'] = [10, 25]

        self.proj_gr = ProjectorGroup(self.pars.groups[0], [self.proj_sh], self.eigvals)

        self.proj_gr.orthogonalize()
        self.proj_gr.calc_complement(self.eigvals)

        temp = self.proj_gr.normion
        self.proj_gr.normion = False
        block_maps, ndim = self.proj_gr.get_block_matrix_map()
        self.proj_gr.normion = temp

        _, ns, nk, _, _ = self.proj_gr.shells[0].proj_win.shape

# Note that 'ns' and 'nk' are the same for all shells
        for isp in range(ns):
            for ik in range(nk):
                print(('ik',ik))
                bmin = self.proj_gr.ib_win[ik, isp, 0]
                bmax = self.proj_gr.ib_win[ik, isp, 1]+1

                nb = bmax - bmin
                p_mat = np.zeros((ndim, nb), dtype=np.complex128)
                #print(bmin,bmax,nb)
# Combine all projectors of the group to one block projector
                for bl_map in block_maps:
                    p_mat[:, :] = 0.0j  # !!! Clean-up from the last k-point and block!
                    for ibl, block in enumerate(bl_map):
                        i1, i2 = block['bmat_range']
                        ish, ion = block['shell_ion']
                        nlm = i2 - i1 + 1
                        shell = self.proj_gr.shells[ish]
                        p_mat[i1:i2, :nb] = shell.proj_win[ion, isp, ik, :nlm, :nb]

                overlap_L = np.dot(p_mat.conjugate().transpose(),p_mat)
                overlap_N = np.dot(p_mat,p_mat.conjugate().transpose())

                assert np.all(np.abs(np.eye(overlap_N.shape[0]) - overlap_N) < 1e-13)
                assert np.all(np.abs(np.eye(overlap_L.shape[0]) - overlap_L) < 1e-13)


