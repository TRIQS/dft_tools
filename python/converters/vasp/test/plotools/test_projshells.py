
import numpy as np
import vaspio
from inpconf import ConfigParameters
from plotools import select_bands, ProjectorShell
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
    - compare output for a correct input
    """
# Scenario 1
    def test_example(self):
        conf_file = 'example.cfg'
        pars = ConfigParameters(conf_file)
        pars.parse_input()
        vasp_data = vaspio.VaspData('./')

        efermi = vasp_data.doscar.efermi
        eigvals = vasp_data.eigenval.eigs - efermi
        emin = pars.groups[0]['emin']
        emax = pars.groups[0]['emax']
        ib_win, nb_min, nb_max = select_bands(eigvals, emin, emax)

        proj_sh = ProjectorShell(pars.shells[0], vasp_data.plocar.plo)

        proj_sh.select_projectors(ib_win, nb_min, nb_max)

        testout = 'projshells.out.test'
        nion, ns, nk, nlm, nbtot = proj_sh.proj_win.shape
        with open(testout, 'wt') as f:
            f.write("pars: %s\n"%(pars.shells[0]))
            for ion in xrange(nion):
                for isp in xrange(ns):
                    for ik in xrange(nk):
                        ib1 = ib_win[ik, 0, 0]
                        ib2 = ib_win[ik, 0, 1]
                        f.write("%i  %i\n"%(ib1, ib2))
                        for ib in xrange(ib2 - nb_min + 1):
                            for ilm in xrange(nlm):
                                p = proj_sh.proj_win[ion, isp, ik, ilm, ib]
                                f.write("%5i  %s\n"%(ilm+1, p))

        expected_file = 'projshells.out'
        self.assertFileEqual(testout, expected_file)
 
