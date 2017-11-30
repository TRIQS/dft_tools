r"""
Tests for class 'Doscar' from module 'vaspio'
"""
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import mytest
import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import Doscar

################################################################################
#
# TestDoscar
#
################################################################################
class TestDoscar(mytest.MyTestCase):
    """
    Function:

    def Doscar.from_file(vasp_dir, dos_filename)

    Scenarios:
    - correct DOSCAR file

    """
# Scenario 1
    def test_example(self):
        filename = 'DOSCAR.example'
        doscar = Doscar()
        doscar.from_file(vasp_dir=_rpath, dos_filename=filename)

        test_efermi = doscar.efermi
        expected = 5.84395237
        self.assertAlmostEqual(test_efermi, expected)

