from pytriqs.applications.dft.sumk_dft import *
from pytriqs.applications.dft.converters import Wien2kConverter
from pytriqs.gf import *
from pytriqs.archive import *
import pytriqs.utility.mpi as mpi
import numpy
import copy


class TransBasis:
    """
    Computates rotations into a new basis, using the condition that a given property is diagonal in the new basis.
    """

    def __init__(self, SK=None, hdf_datafile=None):
        """
        Initialization of the class. There are two ways to do so:

        - existing SumkLDA class :  when you have an existing SumkLDA instance
        - from hdf5 archive :  when you want to use data from hdf5 archive

        Giving the class instance overrides giving the string for the hdf5 archive.

        Parameters
        ----------
        SK : class SumkLDA, optional
             Existing instance of SumkLDA class.
        hdf5_datafile : string, optional
                        Name of hdf5 archive to be used.

        """

        if SK is None:
            # build our own SK instance
            if hdf_datafile is None:
                mpi.report("trans_basis: give SK instance or HDF filename!")
                return 0

            Converter = Wien2kConverter(filename=hdf_datafile, repacking=False)
            Converter.convert_dft_input()
            del Converter

            self.SK = SumkDFT(hdf_file=hdf_datafile +
                              '.h5', use_dft_blocks=False)
        else:
            self.SK = SK

        self.T = copy.deepcopy(self.SK.T[0])
        self.w = numpy.identity(SK.corr_shells[0]['dim'])

    def calculate_diagonalisation_matrix(self, prop_to_be_diagonal='eal'):
        """
        Calculates the diagonalisation matrix w, and stores it as member of the class.

        Parameters
        ----------
        prop_to_be_diagonal : string, optional
                              Defines the property to be diagonalized.

                              - 'eal' : local hamiltonian (i.e. crystal field)
                              - 'dm' : local density matrix

        Returns
        -------
        wsqr : double
               Measure for the degree of rotation done by the diagonalisation. wsqr=1 means no rotation.

        """

        if prop_to_be_diagonal == 'eal':
            prop = self.SK.eff_atomic_levels()[0]
        elif prop_to_be_diagonal == 'dm':
            prop = self.SK.density_matrix(method='using_point_integration')[0]
        else:
            mpi.report(
                "trans_basis: not a valid quantitiy to be diagonal. Choices are 'eal' or 'dm'.")
            return 0

        if self.SK.SO == 0:
            self.eig, self.w = numpy.linalg.eigh(prop['up'])
            # calculate new Transformation matrix
            self.T = numpy.dot(self.T.transpose().conjugate(),
                               self.w).conjugate().transpose()
        else:
            self.eig, self.w = numpy.linalg.eigh(prop['ud'])
            # calculate new Transformation matrix
            self.T = numpy.dot(self.T.transpose().conjugate(),
                               self.w).conjugate().transpose()

        # measure for the 'unity' of the transformation:
        wsqr = sum(abs(self.w.diagonal())**2) / self.w.diagonal().size
        return wsqr

    def rotate_gf(self, gf_to_rot):
        """
        Uses the diagonalisation matrix w to rotate a given GF into the new basis.

        Parameters
        ----------
        gf_to_rot : BlockGf
                    Green's function block to rotate.

        Returns
        -------
        gfreturn : BlockGf
                   Green's function rotated into the new basis.
        """

        # build a full GF
        gfrotated = BlockGf(name_block_generator=[(block, GfImFreq(
            indices=inner, mesh=gf_to_rot.mesh)) for block, inner in self.SK.gf_struct_sumk[0]], make_copies=False)

        # transform the CTQMC blocks to the full matrix:
        # ish is the index of the inequivalent shell corresponding to icrsh
        ish = self.SK.corr_to_inequiv[0]
        for block, inner in self.gf_struct_solver[ish].iteritems():
            for ind1 in inner:
                for ind2 in inner:
                    gfrotated[self.SK.solver_to_sumk_block[ish][block]][
                        ind1, ind2] << gf_to_rot[block][ind1, ind2]

        # Rotate using the matrix w
        for bname, gf in gfrotated:
            gfrotated[bname].from_L_G_R(
                self.w.transpose().conjugate(), gfrotated[bname], self.w)

        gfreturn = gf_to_rot.copy()
        # Put back into CTQMC basis:
        for block, inner in self.gf_struct_solver[ish].iteritems():
            for ind1 in inner:
                for ind2 in inner:
                    gfreturn[block][ind1, ind2] << gfrotated[
                        self.SK.solver_to_sumk_block[0][block]][ind1, ind2]

        return gfreturn

    def write_trans_file(self, filename):
        """
        Writes the new transformation T into a file readable by dmftproj. By that, the requested quantity is
        diagonal already at input.

        Parameters
        ----------
        filename : string
                   Name of the file where the transformation is stored.
        """

        f = open(filename, 'w')
        Tnew = self.T.conjugate()
        dim = self.SK.corr_shells[0]['dim']

        if self.SK.SO == 0:

            for i in range(dim):
                st = ''
                for k in range(dim):
                    st += " %9.6f" % (Tnew[i, k].real)
                    st += " %9.6f" % (Tnew[i, k].imag)
                for k in range(2 * dim):
                    st += " 0.0"

                if i < (dim - 1):
                    f.write("%s\n" % (st))
                else:
                    st1 = st.replace(' ', '*', 1)
                    f.write("%s\n" % (st1))

            for i in range(dim):
                st = ''
                for k in range(2 * dim):
                    st += " 0.0"
                for k in range(dim):
                    st += " %9.6f" % (Tnew[i, k].real)
                    st += " %9.6f" % (Tnew[i, k].imag)

                if i < (dim - 1):
                    f.write("%s\n" % (st))
                else:
                    st1 = st.replace(' ', '*', 1)
                    f.write("%s\n" % (st1))

        else:

            for i in range(dim):
                st = ''
                for k in range(dim):
                    st += " %9.6f" % (Tnew[i, k].real)
                    st += " %9.6f" % (Tnew[i, k].imag)

                if i < (dim - 1):
                    f.write("%s\n" % (st))
                else:
                    st1 = st.replace(' ', '*', 1)
                    f.write("%s\n" % (st1))

        f.close()
