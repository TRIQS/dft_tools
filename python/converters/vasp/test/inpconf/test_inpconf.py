r"""
Test suite for module `inpconf.py`.
"""

import unittest
import numpy as np
from inpconf import ConfigParameters
import ConfigParser

class TestSpecialParsers(unittest.TestCase):
    """
    Tests of special parsers.
    """
    def __init__(self, *args, **kwargs):
        """
        Initializes a custom equality function for comparing numpy arrays.
        """
        super(TestSpecialParsers, self).__init__(*args, **kwargs)
        self.addTypeEqualityFunc(np.ndarray, self.is_arrays_equal)

    def is_arrays_equal(self, arr1, arr2, msg=None):
        """
        Raises self.failureException is arrays arr1 and arr2
        are not equal.
        """
        if not np.allclose(arr1, arr2):
            raise self.failureException(msg)

    def setUp(self):
        """
        """
        pass

################################################################################
#
# test_parse_string_logical()
#
################################################################################
    def test_parse_string_logical(self):
        """
        Function:

        def parse_string_logical(self, par_str)

        Scenarios:

        - **if** par_str == 'True' **return** True
        - **if** par_str == 'False' **return** False
        - **if** par_str == '0' **raise** assertion
        """
        conf_pars = ConfigParameters('test1.cfg')

# Scenario 1
        res = conf_pars.parse_string_logical('True')
        self.assertEqual(res, True)

# Scenario 2
        res = conf_pars.parse_string_logical('False')
        self.assertEqual(res, False)

# Scenario 3
        with self.assertRaises(AssertionError):
            conf_pars.parse_string_logical('0')

################################################################################
#
# test_parse_string_ion_list()
#
################################################################################
    def test_parse_string_ion_list(self):
        """
        Function:

        def parse_string_ion_list(self, par_str)

        Scenarios:

        - **if** par_str == '5 6 7 8' **return** array([4, 5, 6, 7])
        - **if** par_str == 'Ni' **raise** NotImplementedError
        - **if** par_str == '0 1' **raise** AssertionError
        - **if** par_str == '5..8' **return** array([4, 5, 6, 7])
        - **if** par_str == '8..5' **raise** AssertionError
        """
        conf_pars = ConfigParameters('test1.cfg')

# Scenario 1
        expected = np.array([4, 5, 6, 7])
        res = conf_pars.parse_string_ion_list('5 6 7 8')
        self.assertEqual(res, expected)

# Scenario 2
        with self.assertRaises(NotImplementedError):
            conf_pars.parse_string_ion_list('Ni')

# Scenario 3
        with self.assertRaises(AssertionError):
            conf_pars.parse_string_ion_list('0 1')

# Scenario 4
        res = conf_pars.parse_string_ion_list('5..8')
        self.assertEqual(res, expected)

# Scenario 5
        err_mess = "First index of the range"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            conf_pars.parse_string_ion_list('8..5')


################################################################################
#
# test_parse_string_tmatrix()
#
################################################################################
    def test_parse_string_tmatrix(self):
        """
        Function:

        def parse_string_tmatrix(self, par_str)

        Parses a matrix defined as a set of rows in the conf-file.

        Scenarios:

        - **if** number of columns is not the same **raise** AssertionError
        - **if** complex matrix is read and the number of columns is odd
          **raise** AssertionError
        - **if** a correct matrix is given **return** an array

        """
        conf_pars = ConfigParameters('test1.cfg')
# Scenario 1
        par_str = "1.0 0.0\n1.0"
        err_mess = "Number of columns"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            conf_pars.parse_string_tmatrix(par_str, real=True)
 
# Scenario 2
        par_str = "1.0 0.0 2.0 1.0 0.0\n0.0 1.0 2.0 3.0 -1.0"
        err_mess = "Complex matrix must"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            conf_pars.parse_string_tmatrix(par_str, real=False)

# Scenario 3
        par_str = "1.0 0.0 2.0 -3.0\n0.0 1.0 -1.0 1.0"
        res = conf_pars.parse_string_tmatrix(par_str, real=False)
        expected = np.array([[1.0, 2.0 - 3.0j], [1.0j, -1.0 + 1.0j]])
        self.assertEqual(res, expected)


################################################################################
#
# test_parse_parameter_set()
#
################################################################################
    def test_parse_parameter_set(self):
        """
        Function:

        def parse_parameter_set(self, section, param_set, excpetion=False)

        Scenarios:

        - **if** config-file section [Shell 1] contains 'LSHELL = 2' **and**
          'lshell' and 'ions' are in `param_set` **return** a dictionary {'lshell': 2}

        - **if** config-file section [Shell 1] contains 'LSHELL = 2' **and**
          'lshell' and 'ions' are in `param_set` and
          exception=True **raise** Exception
        """
        conf_pars = ConfigParameters('test1.cfg')
        param_set = conf_pars.sh_required  # contains 'lshell' and 'ions'

# Scenario 1
        res = conf_pars.parse_parameter_set('Shell 1', param_set)
        expected = {'lshell': 2}
        self.assertDictEqual(res, expected)

# Scenario 2
        section = 'Shell 1'
        err_mess = "Required parameter" # .* in section [%s]"%(section)
        with self.assertRaisesRegexp(Exception, err_mess):
            conf_pars.parse_parameter_set(section, param_set, exception=True)
        
################################################################################
#
# test_parse_shells()
#
################################################################################
    def test_parse_shells(self):
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
        """
# Scenario 1
        conf_pars = ConfigParameters('test2.cfg')
        err_mess = "No projected shells"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            conf_pars.parse_shells()

# Scenario 2
        conf_pars = ConfigParameters('test3.cfg')
        err_mess = "Failed to extract shell indices"
        with self.assertRaisesRegexp(ValueError, err_mess):
            conf_pars.parse_shells()

# Scenario 3
        conf_pars = ConfigParameters('test4.cfg')
        err_mess = "Required parameter"
        with self.assertRaisesRegexp(Exception, err_mess):
            conf_pars.parse_shells()

# Scenario 4
        conf_pars = ConfigParameters('test5.cfg')
        conf_pars.parse_shells()
        res = conf_pars.shells
        expected = {1: {'lshell': 2, 'ion_list': np.array([4, 5, 6, 7])},
                    2: {'lshell': 1, 'ion_list': np.array([0, 1, 2, 3]), 
                        'tmatrix': np.array([[ 0.,  1.,  0.], [ 1.,  0.,  0.], [ 0.,  0.,  1.]])}}
        self.assertSetEqual(set(res.keys()), set(expected.keys()))

        arr = res[1].pop('ion_list')
        arr_exp = expected[1].pop('ion_list')
        self.assertEqual(arr, arr_exp)

        arr = res[2].pop('ion_list')
        arr_exp = expected[2].pop('ion_list')
        self.assertEqual(arr, arr_exp)

        arr = res[2].pop('tmatrix')
        arr_exp = expected[2].pop('tmatrix')
        self.assertEqual(arr, arr_exp)

        self.assertDictEqual(res, expected)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSpecialParsers)
#    unittest.TextTestRunner(verbosity=2, buffer=False).run(suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)

