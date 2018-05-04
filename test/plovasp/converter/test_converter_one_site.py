
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from pytriqs.applications.dft.converters.plovasp.converter import generate_and_output_as_text
from pytriqs.applications.dft.converters import VaspConverter
import mytest

################################################################################
#
# TestConverterOneSite
#
################################################################################
class TestConverterOneSite(mytest.MyTestCase):
    """
    Function:

    def generate_and_output_as_text(pars, el_struct)
        and
    VaspConverter

    Scenarios:

    - Parse config file and produce a correct converted h5-file
    """
# Scenario 1
    def test_convert_one_site(self):
        generate_and_output_as_text(_rpath + 'example.cfg', _rpath + 'one_site/')

        test_file = _rpath + 'pg_output.test.h5'
        converter = VaspConverter(filename=_rpath + 'vasp', 
                                  hdf_filename=test_file)

        converter.convert_dft_input()

        expected_file = _rpath + 'pg_output.out.h5'
        self.assertH5FileEqual(test_file, expected_file)

if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
