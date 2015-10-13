
import vaspio
from inpconf import ConfigParameters
from plotools import check_data_consistency
from elstruct import ElectronicStructure
import mytest

################################################################################
#
# TestDataConsistency
#
################################################################################
class TestDataConsistency(mytest.MyTestCase):
    """
    Function:

    def check_data_consistency(pars, el_struct)

    Scenarios:

    - **if** a shell contains ions of different types **raise** AssertionError
    """
# Scenario 1
    def test_shell_ion_types(self):
        pass
#        conf_file = 'wrong_shell.cfg'
#        pars = ConfigParameters(conf_file)
#        pars.parse_input()
#        vasp_data = vaspio.VaspData('./', read_all=False)
#        vasp_data.poscar.from_file('./', poscar_filename='POSCAR.complex')
#        el_strct = ElectronicStructure(vasp_data)
#
#        print pars.shells
#        print vasp_data.poscar.type_of_ion
#
#        err_mess = "Each projected shell must"
#        with self.assertRaisesRegexp(Exception, err_mess):
#            check_data_consistency(pars, el_struct)

