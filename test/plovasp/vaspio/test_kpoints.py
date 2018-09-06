r"""
Tests for class 'Ibzkpt' from module 'vaspio'
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import mytest
import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import Kpoints

################################################################################
#
# TestIbzkpt
#
################################################################################
class TestIbzkpt(mytest.MyTestCase):
    """
    Function:

    def read_plocar(filename)

    Scenarios:
    - full IBZKPT file with tetrahedra
    - partial IBZKPT file with k-points only

    """
# Scenario 1
    def test_example(self):
        ibz_file = 'IBZKPT.example'
        kpoints = Kpoints()
        kpoints.from_file(vasp_dir=_rpath, ibz_filename=ibz_file)

        testout = _rpath + 'IBZKPT.example.out.test'
        with open(testout, 'w') as f:
            writeline = lambda s: f.write(s + '\n')
            writeline("nktot = %s"%(kpoints.nktot))
            writeline("ntet = %s"%(kpoints.ntet))
            writeline("volt = %s"%(kpoints.volt))
            writeline("kpts:\n%s"%(kpoints.kpts))
            writeline("tets:\n%s"%(kpoints.itet))

        expected = _rpath + 'IBZKPT.example.out'
        self.assertFileEqual(testout, expected)

# Scenario 2
    def test_notet(self):
        ibz_file = 'IBZKPT.notet'
        kpoints = Kpoints()
        kpoints.from_file(vasp_dir=_rpath, ibz_filename=ibz_file)

        testout = _rpath + 'IBZKPT.notet.out.test'
        with open(testout, 'w') as f:
            writeline = lambda s: f.write(s + '\n')
            writeline("nktot = %s"%(kpoints.nktot))
            writeline("kpts:\n%s"%(kpoints.kpts))

        expected = _rpath + 'IBZKPT.notet.out'
        self.assertFileEqual(testout, expected)


