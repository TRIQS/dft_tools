r"""
Searches and runs all available test suites.
"""
import unittest

if __name__ == '__main__':
    suite = unittest.TestLoader().discover('./')
    unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)
#    unittest.TextTestRunner(verbosity=2, buffer=False).run(suite)

