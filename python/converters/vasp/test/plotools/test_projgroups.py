
import numpy as np
import vaspio
from inpconf import ConfigParameters
from plotools import ProjectorShell, ProjectorGroup
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
    - test output for a correct input
    - test the output of 'orthogonalization()' (sanity check)
    """
    def setUp(self):
        conf_file = 'example.cfg'
        self.pars = ConfigParameters(conf_file)
        self.pars.parse_input()
        self.vasp_data = vaspio.VaspData('./')

        efermi = self.vasp_data.doscar.efermi
        eigvals = self.vasp_data.eigenval.eigs - efermi

        self.shells = [ProjectorShell(self.pars.shells[0], self.vasp_data.plocar.plo)]
        self.proj_gr = ProjectorGroup(self.pars.groups[0], self.shells, eigvals)

# Scenario 1
    def test_example(self):
#        proj_sh.select_projectors(ib_win, nb_min, nb_max)
#
        testout = 'projgroups.out.test'
        nion, ns, nk, nlm, nbtot = self.proj_gr.shells[0].proj_win.shape
        with open(testout, 'wt') as f:
            f.write("pars: %s\n"%(self.pars.groups[0]))
            for ion in xrange(nion):
                for isp in xrange(ns):
                    for ik in xrange(nk):
                        ib1 = self.proj_gr.ib_win[ik, 0, 0]
                        ib2 = self.proj_gr.ib_win[ik, 0, 1]
                        f.write("%i  %i\n"%(ib1, ib2))
                        for ib in xrange(ib2 - self.proj_gr.nb_min + 1):
                            for ilm in xrange(nlm):
                                p = self.proj_gr.shells[0].proj_win[ion, isp, ik, ilm, ib]
                                f.write("%5i  %s\n"%(ilm+1, p))
 
# Scenario 2
    def test_ortho(self):
        self.proj_gr.orthogonalize()

        testout = 'projortho.out.test'
        nion, ns, nk, nlm, nbtot = self.proj_gr.shells[0].proj_win.shape
        with open(testout, 'wt') as f:
            f.write("pars: %s\n"%(self.pars.groups[0]))
            for ion in xrange(nion):
                for isp in xrange(ns):
                    for ik in xrange(nk):
                        ib1 = self.proj_gr.ib_win[ik, 0, 0]
                        ib2 = self.proj_gr.ib_win[ik, 0, 1]
                        f.write("%i  %i\n"%(ib1, ib2))
                        for ib in xrange(ib2 - self.proj_gr.nb_min + 1):
                            for ilm in xrange(nlm):
                                p = self.proj_gr.shells[0].proj_win[ion, isp, ik, ilm, ib]
                                f.write("%5i  %s\n"%(ilm+1, p))

        expected_file = 'projortho.out'
        self.assertFileEqual(testout, expected_file)

