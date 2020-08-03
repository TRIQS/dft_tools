r"""
Tests of 'parse_input()' defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters

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
        with self.assertRaisesRegex(AssertionError, err_mess):
            conf_pars.parse_input()

# Scenario 2
    def test_gr_required(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_2.cfg')
        err_mess = "One \[Shell\] section is"
        with self.assertRaisesRegex(KeyError, err_mess):
            conf_pars.parse_input()

# Scenario 3
    def test_no_shell(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_3.cfg')
        err_mess = "Shell 3 referenced in"
        with self.assertRaisesRegex(Exception, err_mess):
            conf_pars.parse_input()

# Scenario 4
    def test_shell_outside_groups(self):
        conf_pars = ConfigParameters(_rpath + 'input_test_4.cfg')
        err_mess = "Some shells are not inside"
        with self.assertRaisesRegex(AssertionError, err_mess):
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
        res = res.replace(" ","") # Remove spaces for comparison

        expected = r"""Shells:
[{'user_index':1,'ions':{'ion_list':[[4],[5],[6],[7]],'nion':4},'lshell':2,'ion_sort':None,'corr':True},{'user_index':2,'ions':{'ion_list':[[0],[1],[2],[3]],'nion':4},'lshell':1,'tmatrix':array([[0.,1.,0.],
[1.,0.,0.],
[0.,0.,1.]]),'ion_sort':None,'corr':True},{'user_index':3,'ions':{'ion_list':[[0],[1],[2],[3]],'nion':4},'lshell':3,'ion_sort':None,'corr':True}]

Groups:
[{'index':1,'shells':[0,1],'ewindow':(-7.6,3.0),'normalize':True,'normion':True,'complement':False},{'index':2,'shells':[2],'ewindow':(-1.6,2.0),'normalize':True,'normion':True,'complement':False}]"""

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
        res = res.replace(" ","") # Remove spaces for comparison

        expected = r"""Shells:
[{'user_index':1,'ions':{'ion_list':[[4],[5],[6],[7]],'nion':4},'lshell':2,'ion_sort':None,'corr':True}]

Groups:
[{'index':'1','ewindow':(-7.6,3.0),'normalize':True,'normion':True,'complement':False,'shells':[0]}]"""

        self.assertEqual(res, expected)


if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)

