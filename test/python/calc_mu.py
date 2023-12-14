##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2023 by A. Hampel
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

import unittest
import numpy as np

from h5 import *
from triqs.gf import MeshImFreq

from triqs_dft_tools.sumk_dft import *


class test_solver(unittest.TestCase):

    def setUp(self):
        self.iw_mesh = MeshImFreq(beta=40, statistic='Fermion', n_iw=300)
        self.dlr_mesh = MeshDLRImFreq(beta=40, statistic='Fermion', w_max=10, eps=1e-10)
        self.w_mesh = MeshReFreq(n_w=1001, window=(-3,3))
        # magic reference value for the Wien2k SVO t2g example
        self.ref_mu = 0.281
        self.ref_mu_real = 0.215

    def test_dichotomy(self):
        sumk = SumkDFT('SrVO3.ref.h5', mesh=self.iw_mesh)
        mu = sumk.calc_mu(method='dichotomy', precision=0.001, delta=0.1)
        self.assertTrue(abs(self.ref_mu - mu) < 0.01)

    def test_dichotomy_dlr(self):
        sumk = SumkDFT('SrVO3.ref.h5', mesh=self.dlr_mesh)
        mu = sumk.calc_mu(method='dichotomy', precision=0.001, delta=0.1)
        self.assertTrue(abs(self.ref_mu - mu) < 0.01)

    def test_dichotomy_real(self):
        sumk = SumkDFT('SrVO3.ref.h5', mesh=self.w_mesh)
        mu = sumk.calc_mu(method='dichotomy', precision=0.001, delta=0.1, broadening = 0.01, beta=1000)
        self.assertTrue(abs(self.ref_mu_real - mu) < 0.001)

    def test_brent(self):
        sumk = SumkDFT('SrVO3.ref.h5', mesh=self.iw_mesh)
        mu = sumk.calc_mu(method='brent', precision=0.001, delta=0.1)
        self.assertTrue(abs(self.ref_mu - mu) < 0.01)

    # Newton seems to be quite unstable here. Not tested right now!
    # def test_newton(self):
    #     sumk = SumkDFT('SrVO3.ref.h5', mesh = self.iw_mesh)
    #     mu = sumk.calc_mu(method='newton', precision= 0.1, delta=0.01)
    #     self.assertTrue(abs(self.ref_mu - mu) < 0.01)


if __name__ == '__main__':
    unittest.main()
