r"""
Tests for class 'Poscar' from module 'vaspio'
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import mytest
import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import Poscar

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
    - check 'type_of_ion' array for a complex POSCAR file

    """
# Scenario 1
    def test_example(self):
        filename = 'POSCAR.example'
        poscar = Poscar()
        poscar.from_file(vasp_dir=_rpath, poscar_filename=filename)

        testout = _rpath + 'POSCAR.example.out.test'
        with open(testout, 'w') as f:
            writeline = lambda s: f.write(s + '\n')
            writeprop = lambda pname: writeline("%s = %s"%(pname, poscar.__dict__[pname]))

            writeprop('nq')
            writeprop('ntypes')
            writeprop('nions')
            writeprop('el_names')
            writeline("a_brav:\n%s"%(poscar.a_brav))
            writeline("q_types:\n%s"%(poscar.q_types))

        expected = _rpath + 'POSCAR.example.out'
        self.assertFileEqual(testout, expected)

# Scenario 2
    def test_type_of_ion(self):
        filename = 'POSCAR.complex'
        poscar = Poscar()
        poscar.from_file(vasp_dir=_rpath, poscar_filename=filename)

        test_types = 4 * [0] + 4 * [1] + 12 * [2]
        self.assertListEqual(test_types, poscar.type_of_ion)

