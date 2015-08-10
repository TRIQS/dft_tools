r"""
Tests for class 'Doscar' from module 'vaspio'
"""
import mytest
import numpy as np
from vaspio import Doscar

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
        doscar.from_file(vasp_dir='./', dos_filename=filename)

        test_efermi = doscar.efermi
        expected = 5.84395237
        self.assertAlmostEqual(test_efermi, expected)

