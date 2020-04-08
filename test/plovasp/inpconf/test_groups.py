r"""
Tests of 'parse_groups()' defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters

################################################################################
#
# TestParseGroups
#
################################################################################
class TestParseGroups(arraytest.ArrayTestCase):
    """
    Function:

    def parse_groups(self)

    Scenarios:

    - **if** a [Group] section does not contain all required parameters
      **raise** Exception
    - **if** a correct group section is defined **return** a list of dictionaries
    """
# Scenario 1
    def test_gr_required(self):
        conf_pars = ConfigParameters(_rpath + 'parse_groups_1.cfg')
        err_mess = "Required parameter"
        with self.assertRaisesRegex(Exception, err_mess):
            conf_pars.parse_groups()

# Scenario 2
    def test_example(self):
        conf_pars = ConfigParameters(_rpath + 'example.cfg')
        conf_pars.parse_groups()
        res = conf_pars.groups
        expected = [{'index': 1, 'shells': [1, 2], 'ewindow': (-7.6, 3.0),
                     'normalize': True, 'normion': True,'complement': False},
                    {'index': 2, 'shells': [3], 'ewindow': (-1.6, 2.0),
                     'normalize': True, 'normion': True,'complement': False}]
        print(res)
        print(expected)
        self.assertListEqual(res, expected)



