r"""
Tests for class 'Poscar' from module 'vaspio'
"""
import mytest
import numpy as np
from vaspio import Poscar

################################################################################
#
# TestPoscar
#
################################################################################
class TestPoscar(mytest.MyTestCase):
    """
    Function:

    def Poscar.from_file(vasp_dir, poscar_filename)

    Scenarios:
    - correct POSCAR file

    """
# Scenario 1
    def test_example(self):
        filename = 'POSCAR.example'
        poscar = Poscar()
        poscar.from_file(vasp_dir='./', poscar_filename=filename)

        testout = 'POSCAR.example.out.test'
        with open(testout, 'w') as f:
            writeline = lambda s: f.write(s + '\n')
            writeprop = lambda pname: writeline("%s = %s"%(pname, poscar.__dict__[pname]))

            writeprop('nq')
            writeprop('ntypes')
            writeprop('nions')
            writeprop('el_names')
            writeline("a_brav:\n%s"%(poscar.a_brav))
            writeline("q_types:\n%s"%(poscar.q_types))

        expected = 'POSCAR.example.out'
        self.assertFileEqual(testout, expected)

