from pytriqs.applications.dft.sumk_lda import *
from pytriqs.applications.dft.converters import Wien2kConverter
from pytriqs.gf.local.block_gf import BlockGf
from pytriqs.gf.local.gf_imfreq import GfImFreq
from pytriqs.archive import *
import pytriqs.utility.mpi as mpi
import numpy
import copy

class TransBasis:
    '''Computates rotations into a new basis in order to make certain quantities diagonal.'''

    def __init__(self, SK=None, hdf_datafile=None):
        '''Inits the class by reading the input.'''

        if SK is None:
            # build our own SK instance
            if hdf_datafile is None:
                mpi.report("trans_basis: give SK instance or HDF filename!")
                return 0

            Converter = Wien2kConverter(filename=hdf_datafile,repacking=False)
            Converter.convert_dmft_input()
            del Converter

            self.SK = SumkLDA(hdf_file=hdf_datafile+'.h5',use_lda_blocks=False)
        else:
            self.SK = SK

        self.T = copy.deepcopy(self.SK.T[0])
        self.w = numpy.identity(SK.corr_shells[0][3])


    def __call__(self, prop_to_be_diagonal = 'eal'):
        '''Calculates the diagonalisation.'''

        if prop_to_be_diagonal == 'eal':
            prop = self.SK.eff_atomic_levels()[0]
        elif prop_to_be_diagonal == 'dm':
            prop = self.SK.simple_point_dens_mat()[0]
        else:
            mpi.report("trans_basis: not a valid quantitiy to be diagonal. Choices are 'eal' or 'dm'.")
            return 0

        if self.SK.SO == 0:
            self.eig,self.w = numpy.linalg.eigh(prop['up'])
            # calculate new Transformation matrix
            self.T = numpy.dot(self.T.transpose().conjugate(),self.w).conjugate().transpose()
        else:
            self.eig,self.w = numpy.linalg.eigh(prop['ud'])
            # calculate new Transformation matrix
            self.T = numpy.dot(self.T.transpose().conjugate(),self.w).conjugate().transpose()

        # measure for the 'unity' of the transformation:
        wsqr = sum(abs(self.w.diagonal())**2)/self.w.diagonal().size
        return wsqr


    def rotate_gf(self,gf_to_rot):
        '''Rotates a given GF into the new basis.'''

        # build a full GF
        gfrotated = BlockGf( name_block_generator = [ (block,GfImFreq(indices = inner, mesh = gf_to_rot.mesh)) for block,inner in self.SK.gf_struct_sumk[0] ], make_copies = False)

        # transform the CTQMC blocks to the full matrix:
        ish = self.SK.corr_to_inequiv[0]    # ish is the index of the inequivalent shell corresponding to icrsh
        for block, inner in self.gf_struct_solver[ish].iteritems():
            for ind1 in inner:
                for ind2 in inner:
                    gfrotated[self.SK.solver_to_sumk_block[ish][block]][ind1,ind2] << gf_to_rot[block][ind1,ind2]

        # Rotate using the matrix w
        for bname,gf in gfrotated:
            gfrotated[bname].from_L_G_R(self.w.transpose().conjugate(),gfrotated[bname],self.w)

        gfreturn = gf_to_rot.copy()
        # Put back into CTQMC basis:
        for block, inner in self.gf_struct_solver[ish].iteritems():
            for ind1 in inner:
                for ind2 in inner:
                    gfreturn[block][ind1,ind2] << gfrotated[self.SK.solver_to_sumk_block[0][block]][ind1,ind2]

        return gfreturn


    def write_trans_file(self, filename):
        '''Writes the new transformation into a file readable by dmftproj.'''

        f = open(filename,'w')
        Tnew = self.T.conjugate()
        dim = self.SK.corr_shells[0][3]

        if self.SK.SO == 0:

           for i in range(dim):
               st = ''
               for k in range(dim):
                   st += " %9.6f"%(Tnew[i,k].real)
                   st += " %9.6f"%(Tnew[i,k].imag)
               for k in range(2*dim):
                   st += " 0.0"

               if i < (dim-1):
                   f.write("%s\n"%(st))
               else:
                   st1 = st.replace(' ','*',1)
                   f.write("%s\n"%(st1))

           for i in range(dim):
               st = ''
               for k in range(2*dim):
                   st += " 0.0"
               for k in range(dim):
                   st += " %9.6f"%(Tnew[i,k].real)
                   st += " %9.6f"%(Tnew[i,k].imag)

               if i < (dim-1):
                   f.write("%s\n"%(st))
               else:
                   st1 = st.replace(' ','*',1)
                   f.write("%s\n"%(st1))

        else:

            for i in range(dim):
               st = ''
               for k in range(dim):
                   st += " %9.6f"%(Tnew[i,k].real)
                   st += " %9.6f"%(Tnew[i,k].imag)

               if i < (dim-1):
                   f.write("%s\n"%(st))
               else:
                   st1 = st.replace(' ','*',1)
                   f.write("%s\n"%(st1))

        f.close()
