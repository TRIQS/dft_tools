r"""
Tests of 'parse_input()' defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from inpconf import ConfigParameters

################################################################################
#
# TestParseInput
#
################################################################################
class TestParseInput(arraytest.ArrayTestCase):
    """
    Function:

    def parse_input(self)

    Scenarios:

    - **if** no [Group] section exists and more than one [Shell] section
      is given **raise** AssertionError
    - **if** no [Group] section exists but the single [Shell] section
      does not contain required group information **raise** KeyError
    - **if** a shell referenced in a group does not exist
      **raise** Exception
    - **if** not all defined shells are referenced in the groups
      **raise** Exception
    - **if** all sections are parsed error-free check that the output
      is correct
    - correct example with a single shell and no explicit groups
    """
# Scenario 1
    def test_no_group(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_1.cfg')
        err_mess = "At least one group"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            conf_pars.parse_input()

# Scenario 2
    def test_gr_required(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_2.cfg')
        err_mess = "One \[Shell\] section is"
        with self.assertRaisesRegexp(KeyError, err_mess):
            conf_pars.parse_input()

# Scenario 3
    def test_no_shell(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_3.cfg')
        err_mess = "Shell 3 referenced in"
        with self.assertRaisesRegexp(Exception, err_mess):
            conf_pars.parse_input()

# Scenario 4
    def test_shell_outside_groups(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_4.cfg')
        err_mess = "Some shells are not inside"
        with self.assertRaisesRegexp(AssertionError, err_mess):
            conf_pars.parse_input()

# Scenario 5
    def test_example(self):
        conf_pars = ConfigParameters(_rpath + 'example.cfg')
        conf_pars.parse_input()
#        with open('parse_input.output.test', 'wt') as f:
#            f.write("Shells:\n")
#            f.write(conf_pars.shells.__repr__() + '\n\n')
#            f.write("Groups:\n")
#            f.write(conf_pars.groups.__repr__() + '\n')
        res = "Shells:\n"
        res += conf_pars.shells.__repr__() + '\n\n'
        res += "Groups:\n"
        res += conf_pars.groups.__repr__()

        expected = r"""Shells:
[{'ion_list': array([4, 5, 6, 7]), 'user_index': 1, 'lshell': 2}, {'tmatrix': array([[ 0.,  1.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.]]), 'ion_list': array([0, 1, 2, 3]), 'user_index': 2, 'lshell': 1}, {'ion_list': array([0, 1, 2, 3]), 'user_index': 3, 'lshell': 3}]

Groups:
[{'normalize': True, 'index': 1, 'ewindow': (-7.6, 3.0), 'normion': True, 'shells': [0, 1]}, {'normalize': True, 'index': 2, 'ewindow': (-1.6, 2.0), 'normion': True, 'shells': [2]}]"""

        self.assertEqual(res, expected)

# Scenario 6
    def test_example_no_groups(self):
        conf_pars = ConfigParameters(_rpath + 'example_nogroup.cfg')
        conf_pars.parse_input()
#        with open('parse_input.output.test', 'wt') as f:
#            f.write("Shells:\n")
#            f.write(conf_pars.shells.__repr__() + '\n\n')
#            f.write("Groups:\n")
#            f.write(conf_pars.groups.__repr__() + '\n')
        res = "Shells:\n"
        res += conf_pars.shells.__repr__() + '\n\n'
        res += "Groups:\n"
        res += conf_pars.groups.__repr__()

        expected = r"""Shells:
[{'ion_list': array([4, 5, 6, 7]), 'user_index': 1, 'lshell': 2}]

Groups:
[{'normalize': True, 'index': '1', 'ewindow': (-7.6, 3.0), 'shells': [0], 'normion': True}]"""

        self.assertEqual(res, expected)


