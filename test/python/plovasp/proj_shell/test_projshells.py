
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import VaspData
from triqs_dft_tools.converters.plovasp.elstruct import ElectronicStructure
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters
from triqs_dft_tools.converters.plovasp.proj_shell import ProjectorShell
from triqs_dft_tools.converters.plovasp.proj_group import ProjectorGroup
import mytest

################################################################################
#
# TestProjectorShell
#
################################################################################
class TestProjectorShell(mytest.MyTestCase):
    """
    Class:

    ProjectorShell(sh_pars, proj_raw)

    Scenarios:
    - **if** a correct input is given **compare** output files
    - **if** a correct input is given **compare** density matrices
    """
    def setUp(self):
        """
        """
        conf_file = _rpath + 'example.cfg'
        self.pars = ConfigParameters(conf_file)
        self.pars.parse_input()
        vasp_data = VaspData(_rpath + 'one_site/')
        self.el_struct = ElectronicStructure(vasp_data)

#        efermi = vasp_data.doscar.efermi
#        eigvals = vasp_data.eigenval.eigs - efermi
        efermi = self.el_struct.efermi
        eigvals = self.el_struct.eigvals - efermi
        emin, emax = self.pars.groups[0]['ewindow']
        struct = self.el_struct.structure
        kmesh = self.el_struct.kmesh

        self.proj_sh = ProjectorShell(self.pars.shells[0], vasp_data.plocar.plo, vasp_data.plocar.proj_params, kmesh, struct, 0)
        self.proj_gr = ProjectorGroup(self.pars.groups[0], [self.proj_sh], eigvals)

# Scenario 1
    def test_example(self):
        testout = _rpath + 'projshells.out.test'
        nion, ns, nk, nlm, nbtot = self.proj_sh.proj_win.shape
        with open(testout, 'wt') as f:
            f.write("pars: %s\n"%(self.pars.shells[0]))
            for ion in range(nion):
                for isp in range(ns):
                    for ik in range(nk):
                        ib1 = self.proj_sh.ib_win[ik, 0, 0]
                        ib2 = self.proj_sh.ib_win[ik, 0, 1]
                        f.write("%i  %i\n"%(ib1, ib2))
                        for ib in range(ib2 - ib1 + 1):
                            for ilm in range(nlm):
                                p = self.proj_sh.proj_win[ion, isp, ik, ilm, ib]
                                f.write("%5i  %f  %f\n"%(ilm+1, p.real, p.imag))

        expected_file = _rpath + 'projshells.out'
        self.assertFileEqual(testout, expected_file)

# Scenario 2
    def test_dens_mat(self):
        dens_mat, overl = self.proj_sh.density_matrix(self.el_struct)
        testout = _rpath + 'densmat.out.test'
        with open(testout, 'wt') as f:
            f.write("density matrix: %s\n"%(dens_mat))
            f.write("overlap matrix: %s\n"%(overl))

        expected_file = _rpath + 'densmat.out'
        self.assertFileEqual(testout, expected_file)

