#!/usr/bin/env python

import unittest

from app4triqs import chain

class test_chain(unittest.TestCase):

    def test_chain(self):

        i = 111
        j = 222
        ij = chain(i,j)
        self.assertEqual(ij, 111222)

if __name__ == '__main__':
    unittest.main()
