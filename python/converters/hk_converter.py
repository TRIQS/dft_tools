
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn
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
################################################################################

from types import *
import numpy
from pytriqs.archive import *
import pytriqs.utility.mpi as mpi
import string
from math import sqrt


def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').replace('(',' ').replace(')',' ').replace(',',' ').split() : 
	    yield string.atof(x)



class HkConverter:
    """
    Conversion from general H(k) file to an hdf5 file, that can be used as input for the SumK_LDA class.
    """

    def __init__(self, hk_file, hdf_file, lda_subgrp = 'SumK_LDA', symm_subgrp = 'SymmCorr', repacking = False):
        """
        Init of the class. Variable Filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so
        on. 
        """

        assert type(nmto_file)==StringType,"LDA_file must be a filename"
        self.hdf_file = hdf_file
        self.lda_file = hk_file
        #self.Symm_file = Filename+'.symqmc'
        #self.Parproj_file = Filename+'.parproj'
        #self.Symmpar_file = Filename+'.sympar'
        #self.Band_file = Filename+'.outband'
        self.lda_subgrp = lda_subgrp
        self.symm_subgrp = symm_subgrp

        # Checks if h5 file is there and repacks it if wanted:
        import os.path
        if (os.path.exists(self.hdf_file) and repacking):
            self.__repack()
        
        

    def convert_dmft_input(self, only_upper_triangle = True, weights_in_file = False):
        """
        Reads the input files, and stores the data in the HDFfile
        """
        
                   
        if not (mpi.is_master_node()): return # do it only on master:
        mpi.report("Reading input from %s..."%self.lda_file)

        # Read and write only on Master!!!
        # R is a generator : each R.Next() will return the next number in the file
        R = Read_Fortran_File(self.lda_file)
        try:
            energy_unit = 1.0                              # the energy conversion factor is 1.0, we assume eV in files
            n_k = int(R.next())                            # read the number of k points
            k_dep_projection = 0                          
            SP = 0                                        # no spin-polarision
            SO = 0                                        # no spin-orbit 
            charge_below = 0.0                            # total charge below energy window is set to 0
            density_required = R.next()                   # density required, for setting the chemical potential
            symm_op = 0                                   # No symmetry groups for the k-sum

            # the information on the non-correlated shells is needed for defining dimension of matrices:
            n_shells = int(R.next())                      # number of shells considered in the Wanniers
                                                          # corresponds to index R in formulas
            # now read the information about the shells:
            shells = [ [ int(R.next()) for i in range(4) ] for icrsh in range(n_shells) ]    # reads iatom, sort, l, dim

            n_corr_shells = int(R.next())                 # number of corr. shells (e.g. Fe d, Ce f) in the unit cell, 
                                                          # corresponds to index R in formulas
            # now read the information about the shells:
            corr_shells = [ [ int(R.next()) for i in range(6) ] for icrsh in range(n_corr_shells) ]    # reads iatom, sort, l, dim, SO flag, irep

            self.inequiv_shells(corr_shells)              # determine the number of inequivalent correlated shells, has to be known for further reading...


            use_rotations = 0
            rot_mat = [numpy.identity(corr_shells[icrsh][3],numpy.complex_) for icrsh in xrange(n_corr_shells)]
            rot_mat_time_inv = [0 for i in range(n_corr_shells)]
            
            # Representative representations are read from file
            n_reps = [1 for i in range(self.n_inequiv_corr_shells)]
            dim_reps = [0 for i in range(self.n_inequiv_corr_shells)]
            
            for icrsh in range(self.n_inequiv_corr_shells):
                n_reps[icrsh] = int(R.next())   # number of representatives ("subsets"), e.g. t2g and eg
                dim_reps[icrsh] = [int(R.next()) for i in range(n_reps[icrsh])]   # dimensions of the subsets
            
            # The transformation matrix:
            # it is of dimension 2l+1, it is taken to be standard d (as in Wien2k)
            T = []
            for icrsh in range(self.n_inequiv_corr_shells):
            #for ish in xrange(self.N_inequiv_corr_shells):
                ll = 2*corr_shells[self.invshellmap[icrsh]][2]+1
                lmax = ll * (corr_shells[self.invshellmap[icrsh]][4] + 1)
                T.append(numpy.zeros([lmax,lmax],numpy.complex_))
                
                T[icrsh] = numpy.array([[0.0, 0.0, 1.0, 0.0, 0.0],
                                       [1.0/sqrt(2.0), 0.0, 0.0, 0.0, 1.0/sqrt(2.0)],
                                       [-1.0/sqrt(2.0), 0.0, 0.0, 0.0, 1.0/sqrt(2.0)],
                                       [0.0, 1.0/sqrt(2.0), 0.0, -1.0/sqrt(2.0), 0.0],
                                       [0.0, 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0), 0.0]])

    
            # Spin blocks to be read:
            n_spin_blocks = SP + 1 - SO   # number of spins to read for Norbs and Ham, NOT Projectors
                 
        
            # define the number of N_Orbitals for all k points: it is the number of total bands and independent of k!
            n_orb = sum([ shells[ish][3] for ish in range(n_shells)])
            #n_orbitals = [ [n_orb for isp in range(n_spin_blocks)] for ik in xrange(n_k)]
            n_orbitals = numpy.ones([n_k,n_spin_blocs],numpy.int) * n_orb
            #print N_Orbitals

            # Initialise the projectors:
            #proj_mat = [ [ [numpy.zeros([corr_shells[icrsh][3], n_orbitals[ik][isp]], numpy.complex_) 
            #                for icrsh in range (n_corr_shells)] 
            #               for isp in range(n_spin_blocks)] 
            #             for ik in range(n_k) ]
            proj_mat = numpy.zeros([n_k,n_spin_blocs,n_corr_shells,max(numpy.array(corr_shells)[:,3]),max(n_orbitals)],numpy.complex_)


            # Read the projectors from the file:
            for ik in xrange(n_k):
                for icrsh in range(n_corr_shells):
                    # calculate the offset:
                    offset = 0
                    no = 0
                    for i in range(n_shells):
                        if (no==0):
                            if ((shells[i][0]==corr_shells[icrsh][0]) and (shells[i][1]==corr_shells[icrsh][1])):
                                no = corr_shells[icrsh][3]
                            else:
                                offset += shells[i][3]

                    proj_mat[ik,isp,icrsh,0:no,offset:offset+no] = numpy.identity(no)
                    
                      
          
            # now define the arrays for weights and hopping ...
            bz_weights = numpy.ones([n_k],numpy.float_)/ float(n_k)  # w(k_index),  default normalisation 
            #hopping = [ [numpy.zeros([n_orbitals[ik][isp],n_orbitals[ik][isp]],numpy.complex_) 
            #             for isp in range(n_spin_blocks)] for ik in xrange(n_k) ]
            hopping = numpy.zeros([n_k,n_spin_blocs,max(n_orbitals),max(n_orbitals)],numpy.complex_)

            if (weights_in_file):
                # weights in the file
                for ik in xrange(n_k) : bz_weights[ik] = R.next()
                
            # if the sum over spins is in the weights, take it out again!!
            sm = sum(bz_weights)
            bz_weights[:] /= sm 


            # Grab the H
            for isp in range(n_spin_blocks):
                for ik in xrange(n_k) :
                    no = n_orbitals[ik][isp]
                    for i in xrange(no):
                        if (only_upper_triangle):
                            ii=i
                        else:
                            ii = 0
                        for j in xrange(ii,no):
                            hopping[ik,isp,i,j] = R.next()
                            hopping[ik,isp,i,j] += R.next() * 1j
                            
                            if ((only_upper_triangle)and(i!=j)): hopping[ik,isp,j,i] = hopping[ik,isp,i,j].conjugate()
            
            #keep some things that we need for reading parproj:
            self.n_shells = n_shells
            self.shells = shells
            self.n_corr_shells = n_corr_shells
            self.corr_shells = corr_shells
            self.n_spin_blocks = n_spin_blocks
            self.n_orbitals = n_orbitals
            self.n_k = n_k
            self.SO = SO
            self.SP = SP
            self.energy_unit = energy_unit
        except StopIteration : # a more explicit error if the file is corrupted.
            raise "SumK_LDA : reading file HMLT_file failed!"

        R.close()
        
        #print Proj_Mat[0]

        #-----------------------------------------
        # Store the input into HDF5:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.lda_subgrp in ar): ar.create_group(self.lda_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created.
        # If it exists, the data is overwritten!!!
        
        ar[self.lda_subgrp]['energy_unit'] = energy_unit
        ar[self.lda_subgrp]['n_k'] = n_k
        ar[self.lda_subgrp]['k_dep_projection'] = k_dep_projection
        ar[self.lda_subgrp]['SP'] = SP
        ar[self.lda_subgrp]['SO'] = SO
        ar[self.lda_subgrp]['charge_below'] = charge_below
        ar[self.lda_subgrp]['density_required'] = density_required
        ar[self.lda_subgrp]['symm_op'] = symm_op
        ar[self.lda_subgrp]['n_shells'] = n_shells
        ar[self.lda_subgrp]['shells'] = shells
        ar[self.lda_subgrp]['n_corr_shells'] = n_corr_shells
        ar[self.lda_subgrp]['corr_shells'] = corr_shells
        ar[self.lda_subgrp]['use_rotations'] = use_rotations
        ar[self.lda_subgrp]['rot_mat'] = rot_mat
        ar[self.lda_subgrp]['rot_mat_time_inv'] = rot_mat_time_inv
        ar[self.lda_subgrp]['n_reps'] = n_reps
        ar[self.lda_subgrp]['dim_reps'] = dim_reps
        ar[self.lda_subgrp]['T'] = T
        ar[self.lda_subgrp]['n_orbitals'] = n_orbitals
        ar[self.lda_subgrp]['proj_mat'] = proj_mat
        ar[self.lda_subgrp]['bz_weights'] = bz_weights
        ar[self.lda_subgrp]['hopping'] = hopping
        
        del ar
             
       

        

    def __repack(self):
        """Calls the h5repack routine, in order to reduce the file size of the hdf5 archive.
           Should only be used BEFORE the first invokation of HDF_Archive in the program, otherwise
           the hdf5 linking is broken!!!"""

        import subprocess

        if not (mpi.is_master_node()): return

        mpi.report("Repacking the file %s"%self.hdf_file)

        retcode = subprocess.call(["h5repack","-i%s"%self.hdf_file, "-otemphgfrt.h5"])
        if (retcode!=0):
            mpi.report("h5repack failed!")
        else:
            subprocess.call(["mv","-f","temphgfrt.h5","%s"%self.hdf_file])
            


    def inequiv_shells(self,lst):
        """
        The number of inequivalent shells is calculated from lst, and a mapping is given as
        map(i_corr_shells) = i_inequiv_corr_shells
        invmap(i_inequiv_corr_shells) = i_corr_shells
        in order to put the Self energies to all equivalent shells, and for extracting Gloc
        """

        tmp = []
        self.shellmap = [0 for i in range(len(lst))]
        self.invshellmap = [0]
        self.n_inequiv_corr_shells = 1
        tmp.append( lst[0][1:3] )
        
        if (len(lst)>1):
            for i in range(len(lst)-1):
               
                fnd = False
                for j in range(self.n_inequiv_corr_shells):
                    if (tmp[j]==lst[i+1][1:3]):
                        fnd = True
                        self.shellmap[i+1] = j
                if (fnd==False):
                    self.shellmap[i+1] = self.n_inequiv_corr_shells
                    self.n_inequiv_corr_shells += 1
                    tmp.append( lst[i+1][1:3] )
                    self.invshellmap.append(i+1)
                                
