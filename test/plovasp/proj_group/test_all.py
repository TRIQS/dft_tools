r"""
Test suite for module `plotools`.
"""
import unittest

if __name__ == '__main__':
#    suite = unittest.TestLoader().discover('./')
    suite = unittest.TestLoader().discover('./', pattern='test_one*')
    unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)
#    unittest.TextTestRunner(verbosity=2, buffer=False).run(suite)

