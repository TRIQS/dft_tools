
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from triqs_dft_tools.converters.plovasp.converter import generate_and_output_as_text
from triqs_dft_tools.converters import VaspConverter
import mytest

################################################################################
#
# TestConverterLuNiO3
#
################################################################################
class TestConverterLuNiO3(mytest.MyTestCase):
    """
    Function:

    def generate_and_output_as_text(pars, el_struct)
        and
    VaspConverter

    Scenarios:

    - Parse the config file and produce a correct h5-file for DFTTools.
    """
# Scenario 1
    def test_convert_lunio3(self):
        generate_and_output_as_text(_rpath + 'lunio3.cfg', _rpath + 'lunio3/')

        test_file = _rpath + 'lunio3.test.h5'
        converter = VaspConverter(filename=_rpath + 'lunio3',
                                  hdf_filename=test_file)

        converter.convert_dft_input()

        expected_file = _rpath + 'lunio3.ref.h5'
        self.assertH5FileEqual(test_file, expected_file)

if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
