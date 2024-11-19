
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from triqs_dft_tools.converters.plovasp.converter import generate_and_output_as_text
from triqs_dft_tools.converters import VaspConverter
import mytest

################################################################################
#
# TestConverterOneSite
#
################################################################################
class TestConverterSVO(mytest.MyTestCase):
    """
    Function:

    def generate_and_output_as_text(pars, el_struct)
        and
    VaspConverter

    Scenarios:

    - Parse config file and produce a correct converted h5-file
    """
# Scenario 1
    def test_convert_svo(self):
        generate_and_output_as_text(_rpath + 'svo.cfg', _rpath + 'svo/')

        test_file = _rpath + 'svo.test.h5'
        converter = VaspConverter(filename=_rpath + 'svo',
                                  hdf_filename=test_file)

        converter.convert_dft_input()

        expected_file = _rpath + 'svo.ref.h5'
        self.assertH5FileEqual(test_file, expected_file)

if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
