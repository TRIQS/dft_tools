
import numpy as np
import vaspio
from inpconf import ConfigParameters
from plotools import select_bands
import mytest

################################################################################
#
# TestSelectBands
#
################################################################################
class TestSelectBands(mytest.MyTestCase):
    """
    Function:

    def select_bands(eigvals, emin, emax)

    Scenarios:
    - compare output for a correct input
    - **if** emin > max(eigvals) **raise** Exception
    - **if** emax > min(eigvals) **raise** Exception
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

        nb_min_exp = 3
        nb_max_exp = 8 
        ib_win_exp = np.array([[[3, 8]], [[3, 8]], [[3, 7]], [[3, 7]]])

        self.assertEqual(nb_min, nb_min_exp)
        self.assertEqual(nb_max, nb_max_exp)
        self.assertEqual(ib_win, ib_win_exp)

# Scenario 2
    def test_emin_too_large(self):
        conf_file = 'example.cfg'
        pars = ConfigParameters(conf_file)
        pars.parse_input()
        vasp_data = vaspio.VaspData('./')

        efermi = vasp_data.doscar.efermi
        eigvals = vasp_data.eigenval.eigs - efermi
        emin = 20.0
        emax = 25.0
        with self.assertRaisesRegexp(Exception, "Energy window does not overlap"):
            ib_win, nb_min, nb_max = select_bands(eigvals, emin, emax)

# Scenario 3
    def test_emax_too_small(self):
        conf_file = 'example.cfg'
        pars = ConfigParameters(conf_file)
        pars.parse_input()
        vasp_data = vaspio.VaspData('./')

        efermi = vasp_data.doscar.efermi
        eigvals = vasp_data.eigenval.eigs - efermi
        emin = -50.0
        emax = -55.0
        with self.assertRaisesRegexp(Exception, "Energy window does not overlap"):
            ib_win, nb_min, nb_max = select_bands(eigvals, emin, emax)


