r"""
Tests of 'parse_shells()' defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters

################################################################################
#
# TestParseShells
#
################################################################################
class TestParseShells(arraytest.ArrayTestCase):
    """
    Function:

    def parse_shells(self)

    Scenarios:

    - **if** config-file does not contain a valid [Shell] section
      **raise** AssertionError
    - **if** a [Shell] section does not contain a valid index
      **raise** ValueError
    - **if** a [Shell] section does not contain all required parameters
      **raise** Exception
    - **if** two correct [Shell] sections are defined
      **return** a dictionary of shell parameters
    - **if** two correct [Shell] sections (one has CORR=False are defined
      **return** a dictionary of shell parameters
    """
# Scenario 1
    def test_no_shell(self):
        conf_pars = ConfigParameters(_rpath + 'parse_shells_1.cfg')
        err_mess = "No projected shells"
        with self.assertRaisesRegex(AssertionError, err_mess):
            conf_pars.parse_shells()

# Scenario 2
    def test_bad_indices(self):
        conf_pars = ConfigParameters(_rpath + 'parse_shells_2.cfg')
        err_mess = "Failed to extract shell indices"
        with self.assertRaisesRegex(ValueError, err_mess):
            conf_pars.parse_shells()

# Scenario 3
    def test_sh_required(self):
        conf_pars = ConfigParameters(_rpath + 'parse_shells_3.cfg')
        err_mess = "Required parameter"
        with self.assertRaisesRegex(Exception, err_mess):
            conf_pars.parse_shells()

# Scenario 4
    def test_two_shells(self):
        conf_pars = ConfigParameters(_rpath + 'parse_shells_4.cfg')
        conf_pars.parse_shells()
        res = conf_pars.shells
        expected = [{'user_index': 1, 'lshell': 2, 'ions': {'nion': 4, 'ion_list': [[4],[5],[6],[7]]},'corr': True,'ion_sort':None},
                    {'user_index': 2, 'lshell': 1, 'ions': {'nion': 4, 'ion_list': [[0],[1],[2],[3]]},
                        'tmatrix': np.array([[ 0.,  1.,  0.], [ 1.,  0.,  0.], [ 0.,  0.,  1.]]),'corr': True,'ion_sort':None}]
# ...lousy way to test equality of two dictionaries containing numpy arrays
        self.assertEqual(len(res), len(expected))

        arr = res[0].pop('ions')
        arr_exp = expected[0].pop('ions')
        self.assertDictEqual(arr, arr_exp)

        arr = res[1].pop('ions')
        arr_exp = expected[1].pop('ions')
        self.assertDictEqual(arr, arr_exp)

        arr = res[1].pop('tmatrix')
        arr_exp = expected[1].pop('tmatrix')
        self.assertEqual(arr, arr_exp)

        self.assertListEqual(res, expected)

# Scenario 5
    def test_two_shells_with_corr_false(self):
        conf_pars = ConfigParameters(_rpath + 'parse_shells_5.cfg')
        conf_pars.parse_shells()
        res = conf_pars.shells
        expected = [{'user_index': 1, 'lshell': 2, 'ions': {'nion': 4, 'ion_list': [[4],[5],[6],[7]]},'corr': True,'ion_sort':None},
                    {'user_index': 2, 'lshell': 1, 'ions': {'nion': 4, 'ion_list': [[0],[1],[2],[3]]},'corr': False,'ion_sort':None}]
        self.assertEqual(len(res), len(expected))

        arr = res[0].pop('ions')
        arr_exp = expected[0].pop('ions')
        self.assertDictEqual(arr, arr_exp)

        arr = res[1].pop('ions')
        arr_exp = expected[1].pop('ions')
        self.assertDictEqual(arr, arr_exp)

        self.assertListEqual(res, expected)

if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)

