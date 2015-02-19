
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
    - compare output for a correct input
    """
# Scenario 1
    def test_example(self):
        conf_file = 'example.cfg'
        pars = ConfigParameters(conf_file)
        pars.parse_input()
        print pars.groups
        vasp_data = vaspio.VaspData('./')

        efermi = vasp_data.doscar.efermi
        eigvals = vasp_data.eigenval.eigs - efermi

        shells = [ProjectorShell(pars.shells[0], vasp_data.plocar.plo)]
        proj_gr = ProjectorGroup(pars.groups[0], shells, eigvals)

#        proj_sh.select_projectors(ib_win, nb_min, nb_max)
#
        testout = 'projgroups.out.test'
        nion, ns, nk, nbtot, nlm = proj_gr.shells[0].proj_win.shape
        with open(testout, 'wt') as f:
            f.write("pars: %s\n"%(pars.groups[0]))
            for ion in xrange(nion):
                for isp in xrange(ns):
                    for ik in xrange(nk):
                        ib1 = proj_gr.ib_win[ik, 0, 0]
                        ib2 = proj_gr.ib_win[ik, 0, 1]
                        f.write("%i  %i\n"%(ib1, ib2))
                        for ib in xrange(ib2 - proj_gr.nb_min + 1):
                            for ilm in xrange(nlm):
                                p = proj_gr.shells[0].proj_win[ion, isp, ik, ib, ilm]
                                f.write("%5i  %s\n"%(ilm+1, p))

        expected_file = 'projgroups.out'
        self.assertFileEqual(testout, expected_file)
 
