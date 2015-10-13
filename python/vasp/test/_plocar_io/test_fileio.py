r"""
Tests of 'read_plocar()' from module 'plocar_io.c_plocar_io'
"""
import mytest
import numpy as np
from plocar_io.c_plocar_io import read_plocar

################################################################################
#
# TestFileIO
#
################################################################################
class TestFileIO(mytest.MyTestCase):
    """
    Function:

    def read_plocar(filename)

    Scenarios:

    - **if** file PLOCAR does not exist **raise** IOError
    - **if** PLOCAR is truncated **raise** IOError
    - **if** the precision flag is not 4 or 8 **raise** ValueError
    - **if** PLOCAR with prec=8 is read **compare** the output
    - **if** PLOCAR with prec=4 is read **compare** the output
    """
# Scenario 1
    def test_no_plocar(self):
        err_mess = "Error opening xPLOCAR"
        with self.assertRaisesRegexp(IOError, err_mess):
            read_plocar('xPLOCAR')

# Scenario 2
    def test_end_of_file(self):
        err_mess = "End-of-file reading"
        with self.assertRaisesRegexp(IOError, err_mess):
            read_plocar('PLOCAR.trunc')

# Scenario 3
    def test_wrong_prec(self):
        err_mess = "only 'prec = 4, 8' are supported"
        with self.assertRaisesRegexp(ValueError, err_mess):
            read_plocar('PLOCAR.noprec')

# Scenario 4
    def test_plocar_prec8(self):
        pars, plo, ferw = read_plocar('PLOCAR.prec8')
        nion, ns, nk, nb, nlm = plo.shape

        test_file = 'PLOCAR.prec8.out.test' 
        with open(test_file, 'wt') as f:
            f.write(" nlm =%5i\n"%(nlm))
            ion = 1
            isp = 1
            for ik in xrange(nk):
                for ib in xrange(nb):
                    f.write("%5i%5i%5i%5i%10.5f\n"%(ion, isp, ik+1, ib+1, ferw[0, 0, ik, ib]))
                    for ilm in xrange(nlm):
                        p = plo[0, 0, ik, ib, ilm]
                        f.write("%5i%15.7f%15.7f\n"%(ilm+1, p.real, p.imag))
            
        expected_file = 'PLOCAR.prec8.out'
        self.assertFileEqual(test_file, expected_file)

# Scenario 5
    def test_plocar_example(self):
        pars, plo, ferw = read_plocar('PLOCAR.example')
        nion, ns, nk, nb, nlm = plo.shape

        self.assertEqual(pars['nion'], nion)
        self.assertEqual(pars['ns'], ns)
        self.assertEqual(pars['nk'], nk)
        self.assertEqual(pars['nb'], nb)

        test_file = 'PLOCAR.example.out.test' 
        with open(test_file, 'wt') as f:
            f.write("pars: %s\n"%(pars))
            for ion in xrange(nion):
                for isp in xrange(ns):
                    for ik in xrange(nk):
                        for ib in xrange(nb):
                            f.write("%5i%5i%5i%5i  %s\n"%(ion+1, isp+1, ik+1, ib+1,
                                    ferw[ion, isp, ik, ib]))
                            for ilm in xrange(nlm):
                                p = plo[ion, isp, ik, ib, ilm]
                                f.write("%5i  %s\n"%(ilm+1, p))
            
        expected_file = 'PLOCAR.example.out'
        self.assertFileEqual(test_file, expected_file)




