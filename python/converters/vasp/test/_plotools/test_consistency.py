
import vaspio
from inpconf import ConfigParameters
import mytest

################################################################################
#
# TestDataConsistency
#
################################################################################
class TestDataConsistency(mytest.MyTestCase):
    """
    Function:

    def read_plocar(filename)

    Scenarios:

    - **if** file PLOCAR does not exist **raise** IOError
    - **if** PLOCAR is truncated **raise** IOError
    - **if** the precision flag is not 4 or 8 **raise** ValueError
    - **if** PLOCAR with prec=8 is read **compare** the output
    - **if** PLOCAR with prec=4 is read **compare** the output
    """
# Scenario 1
    def test_example(self):
        conf_file = 'example.cfg'
        pars = ConfigParameters(conf_file)
        pars.parse_input()
        vasp_data = vaspio.VaspData('./')

        print pars.shells
        print pars.groups

