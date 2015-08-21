import numpy
import string
from pytriqs.gf.local import *

def read_fortran_file (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exist."%filename
    for line in open(filename,'r') :
        for x in line.replace('D','E').split() :
            yield string.atof(x)

def constr_Sigma_real_axis(self, filename, hdf=True, hdf_dataset='SigmaReFreq',n_om=0,orb=0, tol_mesh=1e-6):
    """Uses Data from files to construct Sigma (or GF)  on the real axis."""

    if not hdf: # then read sigma from text files

        # first get the mesh out of any one of the files:
        bl = self.gf_struct_solver[orb].items()[0][0] # block name
        ol = self.gf_struct_solver[orb].items()[0][1] # list of orbital indices
        if (len(ol)==1): # if blocks are of size one
            Fname = filename+'_'+bl+'.dat'
        else:
            print 'TEST'
            Fname = filename+'_'+bl+'/'+str(ol[0])+'_'+str(ol[0])+'.dat'

        R = read_fortran_file(Fname)
        mesh = numpy.zeros([n_om],numpy.float_)
        try:
            for i in xrange(n_om):
                mesh[i] = R.next()
                sk = R.next()
                sk = R.next()

        except StopIteration : # a more explicit error if the file is corrupted.
            raise "SumkDFT.read_Sigma_ME : reading mesh failed!"
        R.close()

        # check whether the mesh is uniform
        bin = (mesh[n_om-1]-mesh[0])/(n_om-1)
        for i in xrange(n_om):
            assert abs(i*bin+mesh[0]-mesh[i]) < tol_mesh, 'constr_Sigma_ME: real-axis mesh is non-uniform!'

        # construct Sigma
        a_list = [a for a,al in self.gf_struct_solver[orb].iteritems()]
        glist = lambda : [ GfReFreq(indices = al, window=(mesh[0],mesh[n_om-1]),n_points=n_om) for a,al in self.gf_struct_solver[orb].iteritems()]
        SigmaME = BlockGf(name_list = a_list, block_list = glist(),make_copies=False)

        #read Sigma
        for i,g in SigmaME:
            mesh=[w for w in g.mesh]
            for iL in g.indices:
                for iR in g.indices:
                    if (len(g.indices) == 1):
                        Fname = filename+'_%s'%(i)+'.dat'
                    else:
                        Fname = 'SigmaME_'+'%s'%(i)+'_%s'%(iL)+'_%s'%(iR)+'.dat'
                    R = read_fortran_file(Fname)
                    try:
                        for iom in xrange(n_om):
                            sk = R.next()
                            rsig = R.next()
                            isig = R.next()
                            g.data[iom,iL,iR]=rsig+1j*isig
                    except StopIteration : # a more explicit error if the file is corrupted.
                        raise "SumkDFT.read_Sigma_ME : reading Sigma from file failed!"
                    R.close()


    else: # read sigma from hdf

        omega_min=0.0
        omega_max=0.0
        n_om=0
        if (mpi.is_master_node()):
            ar = HDFArchive(filename)
            SigmaME = ar[hdf_dataset]
            del ar

        #OTHER SOLUTION FIXME
        #else:
        #    SigmaME=0
        #SigmaME = mpi.broadcast..

            # we need some parameters to construct Sigma on other nodes
            omega_min=SigmaME.mesh.omega_min
            omega_max=SigmaME.mesh.omega_max
            n_om=len(SigmaME.mesh)
        omega_min=mpi.bcast(omega_min)
        omega_max=mpi.bcast(omega_max)
        n_om=mpi.bcast(n_om)
        mpi.barrier()
        # construct Sigma on other nodes
        if (not mpi.is_master_node()):
            a_list = [a for a,al in self.gf_struct_solver[orb].iteritems()]
            glist = lambda : [ GfReFreq(indices = al, window=(omega_min,omega_max),n_points=n_om) for a,al in self.gf_struct_solver[orb].iteritems()]
            SigmaME = BlockGf(name_list = a_list, block_list = glist(),make_copies=False)
        # pass SigmaME to other nodes
        SigmaME = mpi.bcast(SigmaME)
        mpi.barrier()

    SigmaME.note='ReFreq'

    return SigmaME
