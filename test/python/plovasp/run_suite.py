r"""
Searches and runs all available test suites.
"""
import unittest
import sys

if __name__ == '__main__':
    if len(sys.argv) == 1:
        suite = unittest.TestLoader().discover('./')
    else:
        suite = unittest.TestLoader().discover(sys.argv[1] + '/')

#    def list_tests(suite):
#        for test in suite:
#            if isinstance(test, unittest.TestSuite):
#                list_tests(test)
#            elif isinstance(test, unittest.TestCase):
##                print test.__class__.__bases__
#                tmp = test.__str__().split()
#                test_method, test_class = tmp[0], tmp[1]
#                test_class = test_class.strip('()')
#                print test_class + '.' + test_method

#    list_tests(suite)
    results = unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)
#    unittest.TextTestRunner(verbosity=2, buffer=False).run(suite)
    if results.wasSuccessful():
        raise SystemExit(0)
    else:
        print("Failed tests:")
        for failure in results.failures:
            print(failure[0].__str__())
        raise SystemExit(1)

