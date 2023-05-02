
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from triqs_dft_tools.converters.plovasp.converter import generate_and_output_as_text
from triqs_dft_tools.converters import VaspConverter
import mytest

################################################################################
#
# TestConverterNiO
#
################################################################################
class TestConverterNiO(mytest.MyTestCase):
    """
    Function:

    def generate_and_output_as_text(pars, el_struct)
        and
    VaspConverter

    Scenarios:

    - Parse the config file and produce a correct h5-file for DFTTools.
      In this case we also test that two correlated shells are converted correctly.
    """
# Scenario 1
    def test_convert_nio(self):
        generate_and_output_as_text(_rpath + 'nio.cfg', _rpath + 'nio/')

        test_file = _rpath + 'nio.test.h5'
        converter = VaspConverter(filename=_rpath + 'nio',
                                  hdf_filename=test_file)

        converter.convert_dft_input()

        expected_file = _rpath + 'nio.ref.h5'
        self.assertH5FileEqual(test_file, expected_file)

if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
