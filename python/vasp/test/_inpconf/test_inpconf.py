r"""
Test suite for module `inpconf.py`.
"""
import unittest

if __name__ == '__main__':
#    suite = unittest.TestLoader().discover('./')
    suite = unittest.TestLoader().discover('./', pattern='test_shells.py')
    unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)

