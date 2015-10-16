r"""
Tests of special parseres defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from inpconf import ConfigParameters

################################################################################
#
# TestParseStringLogical
#
################################################################################
class TestParseStringLogical(arraytest.ArrayTestCase):
    """
    Function:

    def parse_string_logical(self, par_str)

    Scenarios:

    - **if** par_str == 'True' **return** True
    - **if** par_str == 'False' **return** False
    - **if** par_str == '0' **raise** assertion
    """
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_true(self):
        res = self.cpars.parse_string_logical('True')
        self.assertEqual(res, True)

# Scenario 2
    def test_false(self):
        res = self.cpars.parse_string_logical('False')
        self.assertEqual(res, False)

# Scenario 3
    def test_incorrect(self):
        with self.assertRaises(AssertionError):
            self.cpars.parse_string_logical('0')

################################################################################
#
# TestParseStringIonList
#
################################################################################
class TestParseStringIonList(arraytest.ArrayTestCase):
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
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_simple_list(self):
        expected = np.array([4, 5, 6, 7])
        res = self.cpars.parse_string_ion_list('5 6 7 8')
        self.assertEqual(res, expected)

# Scenario 2
    def test_atomic_symbol(self):
        with self.assertRaises(NotImplementedError):
            self.cpars.parse_string_ion_list('Ni')

# Scenario 3
    def test_out_of_bounds(self):
        err_mess = "Lowest ion index is"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            self.cpars.parse_string_ion_list('0 1')

# Scenario 4
    def test_list_range(self):
        expected = np.array([4, 5, 6, 7])
        res = self.cpars.parse_string_ion_list('5..8')
        self.assertEqual(res, expected)

# Scenario 5
    def test_range_wrong_order(self):
        err_mess = "First index of the range"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            self.cpars.parse_string_ion_list('8..5')


################################################################################
#
# TestParseStringTmatrix
#
################################################################################
class TestParseStringTmatrix(arraytest.ArrayTestCase):
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
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_number_of_columns(self):
        par_str = "1.0 0.0\n1.0"
        err_mess = "Number of columns"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            self.cpars.parse_string_tmatrix(par_str, real=True)
 
# Scenario 2
    def test_complex_matrix_odd(self):
        par_str = "1.0 0.0 2.0 1.0 0.0\n0.0 1.0 2.0 3.0 -1.0"
        err_mess = "Complex matrix must"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            self.cpars.parse_string_tmatrix(par_str, real=False)

# Scenario 3
    def test_complex_matrix(self):
        par_str = "1.0 0.0 2.0 -3.0\n0.0 1.0 -1.0 1.0"
        res = self.cpars.parse_string_tmatrix(par_str, real=False)
        expected = np.array([[1.0, 2.0 - 3.0j], [1.0j, -1.0 + 1.0j]])
        self.assertEqual(res, expected)


