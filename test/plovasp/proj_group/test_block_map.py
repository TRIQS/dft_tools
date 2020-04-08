
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters
from triqs_dft_tools.converters.plovasp.proj_shell import ProjectorShell
from triqs_dft_tools.converters.plovasp.proj_group import ProjectorGroup
import mytest

################################################################################
#
# TestBlockMap
#
################################################################################
class TestBlockMap(mytest.MyTestCase):
    """
    Function:

    def ProjectorGroup.get_block_matrix_map()

    Scenarios:
    - **test** block matrix for NORMION = False
    - **test** block matrix for NORMION = True
    """
    def setUp(self):
# Mock data
        self.mock_eigvals = np.zeros((1, 11, 1))

        nproj = 16
        self.mock_plo = np.zeros((nproj, 1, 1, 11), dtype=np.complex128)
        self.mock_proj_params = [{} for i in range(nproj)]
        ip = 0
# Mock d-sites
        for isite in range(2):
            for im in range(5):
                self.mock_proj_params[ip]['label'] = 'd-orb'
                self.mock_proj_params[ip]['isite'] = isite + 1
                self.mock_proj_params[ip]['l'] = 2
                self.mock_proj_params[ip]['m'] = im
                ip += 1
# Mock p-sites
        for isite in range(2, 4):
            for im in range(3):
                self.mock_proj_params[ip]['label'] = 'p-orb'
                self.mock_proj_params[ip]['isite'] = isite + 1
                self.mock_proj_params[ip]['l'] = 1
                self.mock_proj_params[ip]['m'] = im
                ip += 1
# Mock k-mesh
        self.mock_kmesh = {'kpoints': np.zeros((1, 3))}
# Mock structure
        self.mock_struct = {'qcoords': np.zeros((4, 3))}

# Scenario 1
    def test_normion_false(self):
        conf_file = _rpath + 'block_matrix.cfg'
        self.pars = ConfigParameters(conf_file)
        self.pars.parse_input()

        shells = []
        for sh_par in self.pars.shells:
            shells.append(ProjectorShell(sh_par, self.mock_plo, self.mock_proj_params, self.mock_kmesh, self.mock_struct, 0))

        proj_gr = ProjectorGroup(self.pars.groups[0], shells, self.mock_eigvals)

        proj_gr.normion = False
        block_maps, ndim = proj_gr.get_block_matrix_map()

        ndim_exp = 16
        block_maps_exp = [[{'bmat_range': (0, 5), 'shell_ion': (0, 0)},
                           {'bmat_range': (5, 10), 'shell_ion': (0, 1)},
                           {'bmat_range': (10, 13), 'shell_ion': (1, 0)},
                           {'bmat_range': (13, 16), 'shell_ion': (1, 1)}]]

        self.assertEqual(ndim, ndim_exp)
        self.assertEqual(block_maps, block_maps_exp)

# Scenario 2
    def test_normion_true(self):
        conf_file = _rpath + 'block_matrix.cfg'
        self.pars = ConfigParameters(conf_file)
        self.pars.parse_input()

        shells = []
        for sh_par in self.pars.shells:
            shells.append(ProjectorShell(sh_par, self.mock_plo, self.mock_proj_params, self.mock_kmesh, self.mock_struct, 0))

        proj_gr = ProjectorGroup(self.pars.groups[0], shells, self.mock_eigvals)

        proj_gr.normion = True
        block_maps, ndim = proj_gr.get_block_matrix_map()

        ndim_exp = 5
        block_maps_exp = [[{'bmat_range': (0, 5), 'shell_ion': (0, 0)}],
                          [{'bmat_range': (0, 5), 'shell_ion': (0, 1)}],
                          [{'bmat_range': (0, 3), 'shell_ion': (1, 0)}],
                          [{'bmat_range': (0, 3), 'shell_ion': (1, 1)}]]

        self.assertEqual(ndim, ndim_exp)
        self.assertEqual(block_maps, block_maps_exp)

