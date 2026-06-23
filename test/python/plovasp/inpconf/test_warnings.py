r"""
Tests that misplaced or unknown tags in the config-file trigger a warning
(see TRIQS/dft_tools#293).
"""
import os
import io
import tempfile
import contextlib
import unittest

import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

import arraytest
from triqs_dft_tools.converters.plovasp.inpconf import ConfigParameters


def _parse_capturing_output(cfg_text):
    """Writes `cfg_text` to a temporary file, parses it, and returns
    (ConfigParameters, captured_stdout)."""
    with tempfile.NamedTemporaryFile('w', suffix='.cfg', delete=False) as fh:
        fh.write(cfg_text)
        fname = fh.name
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            conf_pars = ConfigParameters(fname, verbosity=0)
            conf_pars.parse_input()
        return conf_pars, buf.getvalue()
    finally:
        os.remove(fname)


################################################################################
#
# TestMisplacedTags
#
################################################################################
class TestMisplacedTags(arraytest.ArrayTestCase):
    """
    Function:

    def check_for_unknown_parameters(self, section, known_pars)

    Scenarios:

    - **if** a group-tag is placed in [General] **warn** and ignore it
    - **if** an unknown tag is placed anywhere **warn** that it is unrecognized
    - **if** the config is clean **emit** no warning
    """

# Scenario 1: COMPLEMENT (a [Group] tag) misplaced in [General]
    def test_misplaced_complement(self):
        cfg = ("[General]\nHK = True\nCOMPLEMENT = True\n\n"
               "[Group 1]\nSHELLS = 1\nEWINDOW = -9 2\n\n"
               "[Shell 1]\nLSHELL = 2\nIONS = 1\n")
        conf_pars, out = _parse_capturing_output(cfg)
# The tag is ignored, so complement keeps its default value (False)
        self.assertFalse(conf_pars.groups[0]['complement'])
# A warning is issued that points to the correct section
        self.assertIn('WARNING', out)
        self.assertIn('complement', out)
        self.assertIn('[General]', out)
        self.assertIn('[Group]', out)

# Scenario 2: an unrecognized keyword (typo)
    def test_unknown_keyword(self):
        cfg = ("[General]\nHK = True\n\n"
               "[Group 1]\nSHELLS = 1\nEWINDOW = -9 2\nNORMALISE = True\n\n"
               "[Shell 1]\nLSHELL = 2\nIONS = 1\n")
        _, out = _parse_capturing_output(cfg)
        self.assertIn('WARNING', out)
        self.assertIn('normalise', out)
        self.assertIn('not a recognized keyword', out)

# Scenario 3: a clean config produces no warning
    def test_clean_config(self):
        cfg = ("[General]\nHK = True\n\n"
               "[Group 1]\nSHELLS = 1\nEWINDOW = -9 2\nCOMPLEMENT = True\n\n"
               "[Shell 1]\nLSHELL = 2\nIONS = 1\n")
        conf_pars, out = _parse_capturing_output(cfg)
        self.assertTrue(conf_pars.groups[0]['complement'])
        self.assertNotIn('WARNING', out)


if __name__ == '__main__':
    unittest.main()
