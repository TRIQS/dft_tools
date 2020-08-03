r"""
Tests of 'parse_general()' defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters

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
        conf_pars = ConfigParameters(_rpath + 'example.cfg')
        conf_pars.parse_general()
        res = conf_pars.general
        expected = {'basename': 'test_base', 'efermi': 0.1,
                    'dosmesh': {'n_points': 101, 'emin': -8.0, 'emax': 4.0},
                    'hk' : False}
        self.assertDictEqual(res, expected)



