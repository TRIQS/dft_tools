r"""
Tests of 'parse_groups()' defined in ConfigParameters class
"""
import arraytest
import numpy as np
from inpconf import ConfigParameters

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
        conf_pars = ConfigParameters('parse_groups_1.cfg')
        err_mess = "Required parameter"
        with self.assertRaisesRegexp(Exception, err_mess):
            conf_pars.parse_groups()

# Scenario 2
    def test_example(self):
        conf_pars = ConfigParameters('example.cfg')
        conf_pars.parse_groups()
        res = conf_pars.groups
        expected = [{'index': 1, 'shells': [1, 2], 'emin': -7.6, 'emax': 3.0},
                    {'index': 2, 'shells': [3], 'emin': -1.6, 'emax': 2.0}]
        self.assertListEqual(res, expected)



