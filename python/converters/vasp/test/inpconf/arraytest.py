r"""
Module defining a custom TestCase with extra functionality.
"""

import unittest
import numpy as np

class ArrayTestCase(unittest.TestCase):
    """
    Custom TestCase class supporting array equality.
    """
    def __init__(self, *args, **kwargs):
        """
        Initializes a custom equality function for comparing numpy arrays.
        """
        super(ArrayTestCase, self).__init__(*args, **kwargs)
        self.addTypeEqualityFunc(np.ndarray, self.is_arrays_equal)

    def is_arrays_equal(self, arr1, arr2, msg=None):
        """
        Raises self.failureException is arrays arr1 and arr2
        are not equal.
        """
        if not np.allclose(arr1, arr2):
            raise self.failureException(msg)


