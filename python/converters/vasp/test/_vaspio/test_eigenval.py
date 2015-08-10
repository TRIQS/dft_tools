r"""
Tests for class 'Eigneval' from module 'vaspio'
"""
import mytest
import numpy as np
from vaspio import Eigenval

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

    """
# Scenario 1
    def test_example(self):
        filename = 'EIGENVAL.example'
        eigenval = Eigenval()
        eigenval.from_file(vasp_dir='./', eig_filename=filename)

        testout = 'EIGENVAL.example.out.test'
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

        expected = 'EIGENVAL.example.out'
        self.assertFileEqual(testout, expected)

