r"""
Module defining a custom TestCase with extra functionality.
"""

import unittest
import numpy as np
import difflib

class MyTestCase(unittest.TestCase):
    """
    Custom TestCase class supporting additional equality checks:
    - numpy array equality
    - file equality
    """
    def __init__(self, *args, **kwargs):
        """
        Initializes a custom equality function for comparing numpy arrays.
        """
        super(MyTestCase, self).__init__(*args, **kwargs)
        np.set_printoptions(suppress=True)
        self.addTypeEqualityFunc(np.ndarray, self.is_arrays_equal)

    def is_arrays_equal(self, arr1, arr2, msg=None):
        """
        Raises self.failureException is arrays arr1 and arr2
        are not equal.
        """
        if not np.allclose(arr1, arr2):
            raise self.failureException(msg)

    def assertFileEqual(self, file1, file2):
        """
        Compares two files using difflib.
        Empty lines are ignored.
        Files are assumed to be relatively small;
        the data is truncated for files larger than MAX_SIZE bytes.
        """
        MAX_SIZE = 100000
        with open(file1, 'r') as f1:
            str1 = f1.read(MAX_SIZE)
        with open(file2, 'r') as f2:
            str2 = f2.read(MAX_SIZE)
#
# Make a diff
#
# Remove empty lines
        lstr1 = [s for s in str1.splitlines(True) if s.strip() != '']
        lstr2 = [s for s in str2.splitlines(True) if s.strip() != '']
# diff
        delta = difflib.unified_diff(lstr1, lstr2)
# combine delta's to a string
        diff = ''.join(delta)
# if 'diff' is non-empty, files are different
        if diff:
            return self.fail("Files '%s' and '%s' differ"%(file1, file2))


