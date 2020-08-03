r"""
Tests of 'parse_parameter_set()' defined in ConfigParameters class
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
import numpy as np
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters

################################################################################
#
# TestParseParameterSet
#
################################################################################
class TestParseParameterSet(arraytest.ArrayTestCase):
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
    def setUp(self):
        """
        """
# Dummy ConfigParameters object
        self.cpars = ConfigParameters(_rpath + 'test1.cfg')

# Scenario 1
    def test_sh_required(self):
        param_set = self.cpars.sh_required  # contains 'lshell' and 'ions'
        res = self.cpars.parse_parameter_set('Shell 1', param_set)
        expected = {'lshell': 2}
        self.assertDictEqual(res, expected)

# Scenario 2
    def test_sh_required_exception(self):
        section = 'Shell 1'
        param_set = self.cpars.sh_required  # contains 'lshell' and 'ions'
        err_mess = "Required parameter" # .* in section [%s]"%(section)
        with self.assertRaisesRegex(Exception, err_mess):
            self.cpars.parse_parameter_set(section, param_set, exception=True)

