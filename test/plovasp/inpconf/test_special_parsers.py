r"""
Tests of special parseres defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters

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

    - **if** par_str == '5 6 7 8' **return** 'ions' with ion_list = [[4], [5], [6], [7]]
    - **if** par_str == 'Ni' **raise** NotImplementedError
    - **if** par_str == '0 1' **raise** AssertionError
    - **if** par_str == '5..8' **return** 'ions' with ion_list = [[4], [5], [6], [7]]
    - **if** par_str == '8..5' **raise** AssertionError
    - **if** par_str == '[5, 8] [6 7]' **return** 'ions' with ion_list = [[4, 7], [5, 6]]
    """
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_simple_list(self):
        expected = {'nion': 4, 'ion_list': [[4], [5], [6], [7]]}
        res = self.cpars.parse_string_ion_list('5 6 7 8')
        self.assertDictEqual(res, expected)

# Scenario 2
    def test_atomic_symbol(self):
        with self.assertRaises(NotImplementedError):
            self.cpars.parse_string_ion_list('Ni')

# Scenario 3
    def test_out_of_bounds(self):
        err_mess = "Lowest ion index is"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_string_ion_list('0 1')

# Scenario 4
    def test_list_range(self):
        expected = {'nion': 4, 'ion_list': [[4], [5], [6], [7]]}
        res = self.cpars.parse_string_ion_list('5..8')
        self.assertDictEqual(res, expected)

# Scenario 5
    def test_range_wrong_order(self):
        err_mess = "First index of the range"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_string_ion_list('8..5')

# Scenario 6
    def test_eq_classes(self):
        expected = {'nion': 4, 'ion_list': [[4, 7], [5, 6]]}
        res = self.cpars.parse_string_ion_list('[5, 8] [6 7]')
        self.assertDictEqual(res, expected)


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
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_string_tmatrix(par_str, real=True)

# Scenario 2
    def test_complex_matrix_odd(self):
        par_str = "1.0 0.0 2.0 1.0 0.0\n0.0 1.0 2.0 3.0 -1.0"
        err_mess = "Complex matrix must"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_string_tmatrix(par_str, real=False)

# Scenario 3
    def test_complex_matrix(self):
        par_str = "1.0 0.0 2.0 -3.0\n0.0 1.0 -1.0 1.0"
        res = self.cpars.parse_string_tmatrix(par_str, real=False)
        expected = np.array([[1.0, 2.0 - 3.0j], [1.0j, -1.0 + 1.0j]])
        self.assertEqual(res, expected)


################################################################################
#
# TestParseEnergyWindow
#
################################################################################
class TestParseEnergyWindow(arraytest.ArrayTestCase):
    """
    Function:

    def parse_energy_window(self, par_str)

    Scenarios:

    - **if** par_str == '-1.5 3.0' **return** (-1.5, 3.0)
    - **if** par_str == '3.0 -1.5' **raise** AssertionError
    - **if** par_str == '1.0' **raise** AssertionError
    - **if** par_str == 'aaa' **raise** ValueError
    - **if** par_str == '1.5 3.0 2.0' **raise** AssertionError
    """
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_correct_range(self):
        expected = (-1.5, 3.0)
        res = self.cpars.parse_energy_window('-1.5 3.0')
        self.assertEqual(res, expected)

# Scenario 2
    def test_wrong_range(self):
        err_mess = "The first float in EWINDOW"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_energy_window('3.0 -1.5')

# Scenario 3
    def test_one_float(self):
        err_mess = "EWINDOW must be specified"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_energy_window('1.0')

# Scenario 4
    def test_wrong_string(self):
        with self.assertRaises(ValueError):
            self.cpars.parse_energy_window('aaa')

# Scenario 5
    def test_three_floats(self):
        err_mess = "EWINDOW must be specified"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_energy_window('1.5 3.0 2.0')

################################################################################
#
# TestParseBandWindow
#
################################################################################
class TestParseBandWindow(arraytest.ArrayTestCase):
    """
    Function:

    def parse_band_window(self, par_str)

    Scenarios:

    - **if** par_str == '1 10' **return** (1, 10)
    - **if** par_str == '3.0 -1.5' **raise** AssertionError
    - **if** par_str == '1.0' **raise** AssertionError
    - **if** par_str == 'aaa' **raise** ValueError
    - **if** par_str == '1.5 3.0 2.0' **raise** AssertionError
    """
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_correct_range(self):
        expected = (1, 10)
        res = self.cpars.parse_band_window('1 10')
        self.assertEqual(res, expected)

# Scenario 2
    def test_wrong_range(self):
        err_mess = "The first int in BANDS"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_band_window('10 1')

# Scenario 3
    def test_one_float(self):
        err_mess = "BANDS must be specified"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_band_window('1')

# Scenario 4
    def test_wrong_string(self):
        with self.assertRaises(ValueError):
            self.cpars.parse_band_window('aaa')

# Scenario 5
    def test_three_ints(self):
        err_mess = "BANDS must be specified"
        with self.assertRaisesRegex(AssertionError, err_mess):
            self.cpars.parse_band_window('1 2 3')

################################################################################
#
# TestParseFileTmatrix
#
################################################################################
class TestParseFileTmatrix(arraytest.ArrayTestCase):
    """
    Function:

    def parse_file_tmatrix(self, par_str)

    Scenarios:

    - **if** file is correct **return** array()
    - **if** file is incorrect **raise** (-1.5, 3.0)
    """
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_correct_file(self):
        expected = np.array(
[[ -2.52000000e-04, -0.00000000e+00, -4.27145000e-02,  3.00000000e-07, -9.99087300e-01],
 [ -4.13570000e-03, -2.00000000e-07, -9.99078700e-01, -1.00000000e-07,  4.27152000e-02],
 [ -3.80200000e-04,  0.00000000e+00,  6.04452000e-02, -1.00000000e-07, -9.98171400e-01],
 [ -5.14500000e-04, -0.00000000e+00, -9.98171400e-01,  0.00000000e+00, -6.04450000e-02]])

        res = self.cpars.parse_file_tmatrix(_rpath + 'tmatrix_file.dat')
        self.assertEqual(res, expected)

# Scenario 2
    def test_wrong_file(self):
        with self.assertRaises(ValueError):
            self.cpars.parse_file_tmatrix(_rpath + 'test1.cfg')

################################################################################
#
# TestParseStringDosmesh
#
################################################################################
class TestParseStringDosmesh(arraytest.ArrayTestCase):
    """
    Function:

    def parse_string_dosmesh(self, par_str)

    Scenarios:

    - **if** par_str == '-8.0 4.0 101' **return** dictionary
    - **if** par_str == '101' **return** dictionary
    - **if** par_str == '-8.0 101' **raise** ValueError
    - **if** par_str == '8.0' **raise** ValueError
    """
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_range_npoints(self):
        expected = {'n_points': 101, 'emin': -8.0, 'emax': 4.0}
        res = self.cpars.parse_string_dosmesh('-8.0 4.0 101')
        self.assertEqual(res, expected)

# Scenario 2
    def test_only_npoints(self):
        expected_npoints = 101
        res = self.cpars.parse_string_dosmesh('101')
        self.assertTrue(np.isnan(res['emin']))
        self.assertTrue(np.isnan(res['emax']))
        self.assertEqual(res['n_points'], expected_npoints)

# Scenario 3
    def test_two_numbers(self):
        err_mess = "DOSMESH must be either"
        with self.assertRaisesRegex(ValueError, err_mess):
            self.cpars.parse_string_dosmesh('-8.0 101')

# Scenario 4
    def test_wrong_input(self):
        with self.assertRaises(ValueError):
            self.cpars.parse_string_dosmesh('8.0')


if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)

