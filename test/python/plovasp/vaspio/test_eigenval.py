r"""
Tests for class 'Eigneval' from module 'vaspio'
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import mytest
import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import Eigenval

################################################################################
#
# TestEigenval
#
################################################################################
class TestEigenval(mytest.MyTestCase):
    """
    Function:

    def Eigenval.from_file(vasp_dir, eig_filename)

    Scenarios:
    - correct EIGENVAL file
    - wrong EIGENVAL file from old versions of VASP

    """
# Scenario 1
    def test_example(self):
        filename = 'EIGENVAL.example'
        eigenval = Eigenval()
        eigenval.from_file(vasp_dir=_rpath, eig_filename=filename)

        testout = _rpath + 'EIGENVAL.example.out.test'
        with open(testout, 'w') as f:
            writeline = lambda s: f.write(s + '\n')
            writeprop = lambda pname: writeline("%s = %s"%(pname, eigenval.__dict__[pname]))

            writeprop('nq')
            writeprop('ispin')
            writeprop('nelect')
            writeprop('nktot')
            writeprop('nband')
            writeline("kpts:\n%s"%(eigenval.kpts))
            writeline("kwghts:\n%s"%(eigenval.kwghts))
            writeline("eigs:\n%s"%(eigenval.eigs))
            writeline("ferw:\n%s"%(eigenval.ferw))

        expected = _rpath + 'EIGENVAL.example.out'
        self.assertFileEqual(testout, expected)

# Scenario 2
    def test_bad_example(self):
        filename = 'EIGENVAL.wrong'
        eigenval = Eigenval()

        err_mess = "EIGENVAL file is incorrect"
        with self.assertRaisesRegex(AssertionError, err_mess):
            eigenval.from_file(vasp_dir=_rpath, eig_filename=filename)

