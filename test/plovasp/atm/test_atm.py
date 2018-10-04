
import os

import numpy as np
from triqs_dft_tools.converters.plovasp.atm import dos_tetra_weights_3d
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
    - **if** a correct input is given **compare** output arrays
    """
# Scenario 1
    def test_example(self):
        eigs = np.array([-1.5, -1.309017, -1.0, -0.5])
        en = -0.55
        itt = np.array([[1, 0, 1, 2, 3]]).T

        res = dos_tetra_weights_3d(eigs, en, itt)[:, 0]

        r_should = np.zeros(4)
        r_should[0] = 0.000309016992226;
        r_should[1] = 0.000381966005939;
        r_should[2] = 0.000618033984453;
        r_should[3] = 0.017232002550965;

        self.assertEqual(res, r_should)


