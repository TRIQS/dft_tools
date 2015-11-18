
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import numpy as np
import vaspio
import elstruct
from inpconf import ConfigParameters
from proj_shell import ProjectorShell
from proj_group import ProjectorGroup
import mytest

################################################################################
#
# TestProjectorGroup
#
################################################################################
class TestProjectorGroup(mytest.MyTestCase):
    """
    Class:

    ProjectorGroup(sh_pars, proj_raw)

    Scenarios:
    - **test** that orthogonalization is correct
    """
    def setUp(self):
        conf_file = _rpath + 'example.cfg'
        self.pars = ConfigParameters(conf_file)
        self.pars.parse_input()
        vasp_data = vaspio.VaspData(_rpath + 'one_site/')
        self.el_struct = elstruct.ElectronicStructure(vasp_data)

        efermi = vasp_data.doscar.efermi
        self.eigvals = vasp_data.eigenval.eigs - efermi
        ferw = vasp_data.eigenval.ferw

        self.proj_sh = ProjectorShell(self.pars.shells[0], vasp_data.plocar.plo, vasp_data.plocar.proj_params, 0)
        self.proj_gr = ProjectorGroup(self.pars.groups[0], [self.proj_sh], self.eigvals, ferw)

# Scenario 1
    def test_ortho(self):
        self.proj_gr.orthogonalize()

        dens_mat, overl = self.proj_sh.density_matrix(self.el_struct)

        testout = _rpath + 'projortho.out.test'
        with open(testout, 'wt') as f:
            f.write("density matrix: %s\n"%(dens_mat))
            f.write("overlap matrix: %s\n"%(overl))

        self.assertEqual(overl, np.eye(5))

        expected_file = _rpath + 'projortho.out'
        self.assertFileEqual(testout, expected_file)


