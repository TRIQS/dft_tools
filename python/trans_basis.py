from pytriqs.applications.dft.sumk_lda import *
from pytriqs.applications.dft.converters import Wien2kConverter
from pytriqs.gf.local.block_gf import BlockGf
from pytriqs.gf.local.gf_imfreq import GfImFreq
import numpy
from pytriqs.archive import *
import copy
import pytriqs.utility.mpi as mpi

class TransBasis:
    '''Computates rotations into a new basis in order to make certain quantities diagonal.'''


    def __init__(self, SK=None, hdf_datafile=None):
        '''Inits the class by reading the input.'''

        if SK is None:
            # build our own SK instance
            if hdf_datafile is None:
                mpi.report("Give SK instance or HDF filename!")
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

        if (prop_to_be_diagonal=='eal'):
            eal = self.SK.eff_atomic_levels()[0]
        elif (prop_to_be_diagonal=='dm'):
            eal = self.SK.simple_point_dens_mat()[0]
        else:
            mpi.report("Not a valid quantitiy to be diagonal! Choices are 'eal' or 'dm'")
            return 0

        if (self.SK.SO==0):
            self.eig,self.w = numpy.linalg.eigh(eal['up'])

            # now calculate new Transformation matrix
            self.T = numpy.dot(self.T.transpose().conjugate(),self.w).conjugate().transpose()


            #return numpy.dot(self.w.transpose().conjugate(),numpy.dot(eal['up'],self.w))

        else:

            self.eig,self.w = numpy.linalg.eigh(eal['ud'])

            # now calculate new Transformation matrix
            self.T = numpy.dot(self.T.transpose().conjugate(),self.w).conjugate().transpose()


            #MPI.report("SO not implemented yet!")
            #return 0

        # measure for the 'unity' of the transformation:
        wsqr = sum(abs(self.w.diagonal())**2)/self.w.diagonal().size
        return wsqr


    def rotate_gf(self,gf_to_rot):
        '''Rotates a given GF into the new basis.'''

        # build a full GF
        gfrotated = BlockGf( name_block_generator = [ (block,GfImFreq(indices = inner, mesh = gf_to_rot.mesh)) for block,inner in self.SK.gf_struct_sumk[0] ], make_copies = False)

        # transform the CTQMC blocks to the full matrix:
        s = self.SK.shellmap[0]    # s is the index of the inequivalent shell corresponding to icrsh
        for block, inner in self.gf_struct_solver[s].iteritems():
            for i in range(len(inner)):
                for j in range(len(inner)):
                    ind1 = inner[i]
                    ind2 = inner[j]
                    gfrotated[self.SK.solver_to_sumk_block[s][block]][ind1,ind2] << gf_to_rot[block][ind1,ind2]

        # Rotate using the matrix w
        for bname,gf in gfrotated:
            gfrotated[bname].from_L_G_R(self.w.transpose().conjugate(),gfrotated[bname],self.w)

        gfreturn = gf_to_rot.copy()
        # Put back into CTQMC basis:
        for block, inner in self.gf_struct_solver[s].iteritems():
            for i in range(len(inner)):
                for j in range(len(inner)):
                    ind1 = inner[i]
                    ind2 = inner[j]
                    gfreturn[block][ind1,ind2] << gfrotated[self.SK.solver_to_sumk_block[0][block]][ind1,ind2]

        return gfreturn


    def write_trans_file(self, filename):
        '''Writes the new transformation into a file readable by dmftproj.'''

        f=open(filename,'w')

        Tnew = self.T.conjugate()
        N = self.SK.corr_shells[0][3]

        if (self.SK.SO==0):

           for i in range(N):
               st = ''
               for k in range(N):
                   st += " %9.6f"%(Tnew[i,k].real)
                   st += " %9.6f"%(Tnew[i,k].imag)
               for k in range(2*N):
                   st += " 0.0"

               if (i<(N-1)):
                   f.write("%s\n"%(st))
               else:
                   st1=st.replace(' ','*',1)
                   f.write("%s\n"%(st1))


           for i in range(N):
               st = ''
               for k in range(2*N):
                   st += " 0.0"
               for k in range(N):
                   st += " %9.6f"%(Tnew[i,k].real)
                   st += " %9.6f"%(Tnew[i,k].imag)

               if (i<(N-1)):
                   f.write("%s\n"%(st))
               else:
                   st1=st.replace(' ','*',1)
                   f.write("%s\n"%(st1))

        else:

            for i in range(N):
               st = ''
               for k in range(N):
                   st += " %9.6f"%(Tnew[i,k].real)
                   st += " %9.6f"%(Tnew[i,k].imag)

               if (i<(N-1)):
                   f.write("%s\n"%(st))
               else:
                   st1=st.replace(' ','*',1)
                   f.write("%s\n"%(st1))
            #MPI.report("SO not implemented!")

        f.close()
