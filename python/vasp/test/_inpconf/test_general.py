r"""
Tests of 'parse_general()' defined in ConfigParameters class
"""
import arraytest
import numpy as np
from inpconf import ConfigParameters

################################################################################
#
# TestParseGeneral
#
################################################################################
class TestParseGeneral(arraytest.ArrayTestCase):
    """
    Function:

    def parse_general(self)

    Scenarios:

    - **if** a correct [General] section is defined **return** a dictionary
    """
# Scenario 1
    def test_example(self):
        conf_pars = ConfigParameters('example.cfg')
        conf_pars.parse_general()
        res = conf_pars.general
        expected = {'basename': 'test_base', 'efermi': 0.1}
        self.assertDictEqual(res, expected)



