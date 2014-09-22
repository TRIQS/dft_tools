
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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


def read_fortran_file (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exist."%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
	    yield string.atof(x)



class Wien2kConverter:
    """
    Conversion from Wien2k output to an hdf5 file that can be used as input for the SumkLDA class.
    """

    def __init__(self, filename, lda_subgrp = 'SumK_LDA', symm_subgrp = 'SymmCorr', repacking = False):
        """
        Init of the class. Variable filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so on. 
        """

        assert type(filename)==StringType,"LDA_file must be a filename"
        self.hdf_file = filename+'.h5'
        self.lda_file = filename+'.ctqmcout'
        self.symm_file = filename+'.symqmc'
        self.parproj_file = filename+'.parproj'
        self.symmpar_file = filename+'.sympar'
        self.band_file = filename+'.outband'
        self.lda_subgrp = lda_subgrp
        self.symm_subgrp = symm_subgrp

        # Checks if h5 file is there and repacks it if wanted:
        import os.path
        if (os.path.exists(self.hdf_file) and repacking):
            self.__repack()
        
        

    def convert_dmft_input(self):
        """
        Reads the input files, and stores the data in the HDFfile
        """
        
                   
        # Read and write only on the master node
        if not (mpi.is_master_node()): return
        mpi.report("Reading input from %s..."%self.lda_file)

        # R is a generator : each R.Next() will return the next number in the file
        R = read_fortran_file(self.lda_file)
        try:
            energy_unit = R.next()                         # read the energy convertion factor
            n_k = int(R.next())                            # read the number of k points
            k_dep_projection = 1                          
            SP = int(R.next())                            # flag for spin-polarised calculation
            SO = int(R.next())                            # flag for spin-orbit calculation
            charge_below = R.next()                       # total charge below energy window
            density_required = R.next()                   # total density required, for setting the chemical potential
            symm_op = 1                                   # Use symmetry groups for the k-sum

            # the information on the non-correlated shells is not important here, maybe skip:
            n_shells = int(R.next())                      # number of shells (e.g. Fe d, As p, O p) in the unit cell, 
                                                               # corresponds to index R in formulas
            shells = [ [ int(R.next()) for i in range(4) ] for icrsh in range(n_shells) ]    # reads iatom, sort, l, dim

            n_corr_shells = int(R.next())                 # number of corr. shells (e.g. Fe d, Ce f) in the unit cell, 
                                                          # corresponds to index R in formulas
            # now read the information about the shells:
            corr_shells = [ [ int(R.next()) for i in range(6) ] for icrsh in range(n_corr_shells) ]    # reads iatom, sort, l, dim, SO flag, irep

            self.inequiv_shells(corr_shells)              # determine the number of inequivalent correlated shells, has to be known for further reading...

            use_rotations = 1
            rot_mat = [numpy.identity(corr_shells[icrsh][3],numpy.complex_) for icrsh in xrange(n_corr_shells)]
           
            # read the matrices
            rot_mat_time_inv = [0 for i in range(n_corr_shells)]

            for icrsh in xrange(n_corr_shells):
                for i in xrange(corr_shells[icrsh][3]):    # read real part:
                    for j in xrange(corr_shells[icrsh][3]):
                        rot_mat[icrsh][i,j] = R.next()
                for i in xrange(corr_shells[icrsh][3]):    # read imaginary part:
                    for j in xrange(corr_shells[icrsh][3]):
                        rot_mat[icrsh][i,j] += 1j * R.next()

                if (SP==1):             # read time inversion flag:
                    rot_mat_time_inv[icrsh] = int(R.next())
                    
                  
            
            # Read here the info for the transformation of the basis:
            n_reps = [1 for i in range(self.n_inequiv_corr_shells)]
            dim_reps = [0 for i in range(self.n_inequiv_corr_shells)]
            T = []
            for icrsh in range(self.n_inequiv_corr_shells):
                n_reps[icrsh] = int(R.next())   # number of representatives ("subsets"), e.g. t2g and eg
                dim_reps[icrsh] = [int(R.next()) for i in range(n_reps[icrsh])]   # dimensions of the subsets
            
                # The transformation matrix:
                # is of dimension 2l+1 without SO, and 2*(2l+1) with SO!
                ll = 2*corr_shells[self.invshellmap[icrsh]][2]+1
                lmax = ll * (corr_shells[self.invshellmap[icrsh]][4] + 1)
                T.append(numpy.zeros([lmax,lmax],numpy.complex_))
                
                # now read it from file:
                for i in xrange(lmax):
                    for j in xrange(lmax):
                        T[icrsh][i,j] = R.next()
                for i in xrange(lmax):
                    for j in xrange(lmax):
                        T[icrsh][i,j] += 1j * R.next()

    
            # Spin blocks to be read:
            n_spin_blocs = SP + 1 - SO   
                 
        
            # read the list of n_orbitals for all k points
            n_orbitals = numpy.zeros([n_k,n_spin_blocs],numpy.int)
            for isp in range(n_spin_blocs):
                for ik in xrange(n_k):
                    n_orbitals[ik,isp] = int(R.next())
            

            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k,n_spin_blocs,n_corr_shells,max(numpy.array(corr_shells)[:,3]),max(n_orbitals)],numpy.complex_)


            # Read the projectors from the file:
            for ik in xrange(n_k):
                for icrsh in range(n_corr_shells):
                    no = corr_shells[icrsh][3]
                    # first Real part for BOTH spins, due to conventions in dmftproj:
                    for isp in range(n_spin_blocs):
                        for i in xrange(no):
                            for j in xrange(n_orbitals[ik][isp]):
                                proj_mat[ik,isp,icrsh,i,j] = R.next()
                    # now Imag part:
                    for isp in range(n_spin_blocs):
                        for i in xrange(no):
                            for j in xrange(n_orbitals[ik][isp]):
                                proj_mat[ik,isp,icrsh,i,j] += 1j * R.next()

          
            # now define the arrays for weights and hopping ...
            bz_weights = numpy.ones([n_k],numpy.float_)/ float(n_k)  # w(k_index),  default normalisation 
            hopping = numpy.zeros([n_k,n_spin_blocs,max(n_orbitals),max(n_orbitals)],numpy.complex_)

            # weights in the file
            for ik in xrange(n_k) : bz_weights[ik] = R.next()         
                
            # if the sum over spins is in the weights, take it out again!!
            sm = sum(bz_weights)
            bz_weights[:] /= sm 

            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian!!!!
            for isp in range(n_spin_blocs):
                for ik in xrange(n_k) :
                    no = n_orbitals[ik,isp]
                    for i in xrange(no):
                        hopping[ik,isp,i,i] = R.next() * energy_unit
            
            # keep some things that we need for reading parproj:
            self.n_shells = n_shells
            self.shells = shells
            self.n_corr_shells = n_corr_shells
            self.corr_shells = corr_shells
            self.n_spin_blocs = n_spin_blocs
            self.n_orbitals = n_orbitals
            self.n_k = n_k
            self.SO = SO
            self.SP = SP
            self.energy_unit = energy_unit
        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Wien2k_converter : reading file lda_file failed!"

        R.close()
        
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
              
        # Symmetries are used, 
        # Now do the symmetries for correlated orbitals:
        self.read_symmetry_input(orbits=corr_shells,symm_file=self.symm_file,symm_subgrp=self.symm_subgrp,SO=SO,SP=SP)


    def convert_parproj_input(self, par_proj_subgrp='SumK_LDA_ParProj', symm_par_subgrp='SymmPar'):
        """
        Reads the input for the partial charges projectors from case.parproj, and stores it in the symm_par_subgrp
        group in the HDF5.
        """

        if not (mpi.is_master_node()): return

        self.par_proj_subgrp = par_proj_subgrp
        self.symm_par_subgrp = symm_par_subgrp

        mpi.report("Reading parproj input from %s..."%self.parproj_file)

        dens_mat_below = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.n_shells)] 
                           for isp in range(self.n_spin_blocs) ]

        R = read_fortran_file(self.parproj_file)

        n_parproj = [int(R.next()) for i in range(self.n_shells)]
        n_parproj = numpy.array(n_parproj)
                
        # Initialise P, here a double list of matrices:
        proj_mat_pc = numpy.zeros([self.n_k,self.n_spin_blocs,self.n_shells,max(n_parproj),max(numpy.array(self.shells)[:,3]),max(self.n_orbitals)],numpy.complex_)
        
        rot_mat_all = [numpy.identity(self.shells[ish][3],numpy.complex_) for ish in xrange(self.n_shells)]
        rot_mat_all_time_inv = [0 for i in range(self.n_shells)]

        for ish in range(self.n_shells):
            # read first the projectors for this orbital:
            for ik in xrange(self.n_k):
                for ir in range(n_parproj[ish]):
                    for isp in range(self.n_spin_blocs):
                                    
                        for i in xrange(self.shells[ish][3]):    # read real part:
                            for j in xrange(self.n_orbitals[ik][isp]):
                                proj_mat_pc[ik,isp,ish,ir,i,j] = R.next()
                            
                    for isp in range(self.n_spin_blocs):
                        for i in xrange(self.shells[ish][3]):    # read imaginary part:
                            for j in xrange(self.n_orbitals[ik][isp]):
                                proj_mat_pc[ik,isp,ish,ir,i,j] += 1j * R.next()
                                        
                    
            # now read the Density Matrix for this orbital below the energy window:
            for isp in range(self.n_spin_blocs):
                for i in xrange(self.shells[ish][3]):    # read real part:
                    for j in xrange(self.shells[ish][3]):
                        dens_mat_below[isp][ish][i,j] = R.next()
            for isp in range(self.n_spin_blocs):
                for i in xrange(self.shells[ish][3]):    # read imaginary part:
                    for j in xrange(self.shells[ish][3]):
                        dens_mat_below[isp][ish][i,j] += 1j * R.next()
                if (self.SP==0): dens_mat_below[isp][ish] /= 2.0

            # Global -> local rotation matrix for this shell:
            for i in xrange(self.shells[ish][3]):    # read real part:
                for j in xrange(self.shells[ish][3]):
                    rot_mat_all[ish][i,j] = R.next()
            for i in xrange(self.shells[ish][3]):    # read imaginary part:
                for j in xrange(self.shells[ish][3]):
                    rot_mat_all[ish][i,j] += 1j * R.next()
                    
            #print Dens_Mat_below[0][ish],Dens_Mat_below[1][ish]
            
            if (self.SP):
                rot_mat_all_time_inv[ish] = int(R.next())

        R.close()

        #-----------------------------------------
        # Store the input into HDF5:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.par_proj_subgrp in ar): ar.create_group(self.par_proj_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created.
        # If it exists, the data is overwritten!!!
        thingstowrite = ['dens_mat_below','n_parproj','proj_mat_pc','rot_mat_all','rot_mat_all_time_inv']
        for it in thingstowrite: exec "ar['%s']['%s'] = %s"%(self.par_proj_subgrp,it,it)
        del ar

        # Symmetries are used, 
        # Now do the symmetries for all orbitals:
        self.read_symmetry_input(orbits=self.shells,symm_file=self.symmpar_file,symm_subgrp=self.symm_par_subgrp,SO=self.SO,SP=self.SP)


    def convert_bands_input(self, bands_subgrp = 'SumK_LDA_Bands'):
        """
        Converts the input for momentum resolved spectral functions, and stores it in bands_subgrp in the
        HDF5.
        """

        if not (mpi.is_master_node()): return

        self.bands_subgrp = bands_subgrp
        mpi.report("Reading bands input from %s..."%self.band_file)

        R = read_fortran_file(self.band_file)
        try:
            n_k = int(R.next())

            # read the list of n_orbitals for all k points
            n_orbitals = numpy.zeros([n_k,self.n_spin_blocs],numpy.int)
            for isp in range(self.n_spin_blocs):
                for ik in xrange(n_k):
                    n_orbitals[ik,isp] = int(R.next())

            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k,self.n_spin_blocs,self.n_corr_shells,max(numpy.array(self.corr_shells)[:,3]),max(n_orbitals)],numpy.complex_)

            # Read the projectors from the file:
            for ik in xrange(n_k):
                for icrsh in range(self.n_corr_shells):
                    no = self.corr_shells[icrsh][3]
                    # first Real part for BOTH spins, due to conventions in dmftproj:
                    for isp in range(self.n_spin_blocs):
                        for i in xrange(no):
                            for j in xrange(n_orbitals[ik,isp]):
                                proj_mat[ik,isp,icrsh,i,j] = R.next()
                    # now Imag part:
                    for isp in range(self.n_spin_blocs):
                        for i in xrange(no):
                            for j in xrange(n_orbitals[ik,isp]):
                                proj_mat[ik,isp,icrsh,i,j] += 1j * R.next()

            hopping = numpy.zeros([n_k,self.n_spin_blocs,max(n_orbitals),max(n_orbitals)],numpy.complex_)
         	    
            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian!!!!
            for isp in range(self.n_spin_blocs):
                for ik in xrange(n_k) :
                    no = n_orbitals[ik,isp]
                    for i in xrange(no):
                        hopping[ik,isp,i,i] = R.next() * self.energy_unit

            # now read the partial projectors:
            n_parproj = [int(R.next()) for i in range(self.n_shells)]
            n_parproj = numpy.array(n_parproj)
            
            # Initialise P, here a double list of matrices:
            proj_mat_pc = numpy.zeros([n_k,self.n_spin_blocs,self.n_shells,max(n_parproj),max(numpy.array(self.shells)[:,3]),max(n_orbitals)],numpy.complex_)


            for ish in range(self.n_shells):
               
                for ik in xrange(n_k):
                    for ir in range(n_parproj[ish]):
                        for isp in range(self.n_spin_blocs):
                                    
                            for i in xrange(self.shells[ish][3]):    # read real part:
                                for j in xrange(n_orbitals[ik,isp]):
                                    proj_mat_pc[ik,isp,ish,ir,i,j] = R.next()
                            
                            for i in xrange(self.shells[ish][3]):    # read imaginary part:
                                for j in xrange(n_orbitals[ik,isp]):
                                    proj_mat_pc[ik,isp,ish,ir,i,j] += 1j * R.next()

        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Wien2k_converter : reading file band_file failed!"

        R.close()
        # reading done!

        #-----------------------------------------
        # Store the input into HDF5:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.bands_subgrp in ar): ar.create_group(self.bands_subgrp) 

        # The subgroup containing the data. If it does not exist, it is created.
        # If it exists, the data is overwritten!!!
        thingstowrite = ['n_k','n_orbitals','proj_mat','hopping','n_parproj','proj_mat_pc']
        for it in thingstowrite: exec "ar['%s']['%s'] = %s"%(self.bands_subgrp,it,it)
        del ar
   




    def read_symmetry_input(self, orbits, symm_file, symm_subgrp, SO, SP):
        """
        Reads input for the symmetrisations from symm_file, which is case.sympar or case.symqmc.
        """

        if not (mpi.is_master_node()): return

        mpi.report("Reading symmetry input from %s..."%symm_file)

        n_orbits = len(orbits)
        R=read_fortran_file(symm_file)

        try:
            n_s = int(R.next())           # Number of symmetry operations
            n_atoms = int(R.next())       # number of atoms involved
            perm = [ [int(R.next()) for i in xrange(n_atoms)] for j in xrange(n_s) ]    # list of permutations of the atoms
            if SP: 
                time_inv = [ int(R.next()) for j in xrange(n_s) ]           # timeinversion for SO xoupling
            else:
                time_inv = [ 0 for j in xrange(n_s) ] 

            # Now read matrices:
            mat = []  
            for in_s in xrange(n_s):
                
                mat.append( [ numpy.zeros([orbits[orb][3], orbits[orb][3]],numpy.complex_) for orb in xrange(n_orbits) ] )
                for orb in range(n_orbits):
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat[in_s][orb][i,j] = R.next()            # real part
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat[in_s][orb][i,j] += 1j * R.next()      # imaginary part

            # determine the inequivalent shells:
            #SHOULD BE FINALLY REMOVED, PUT IT FOR ALL ORBITALS!!!!! (PS: FIXME?)
            #self.inequiv_shells(orbits)
            mat_tinv = [numpy.identity(orbits[orb][3],numpy.complex_)
                        for orb in range(n_orbits)]

            if ((SO==0) and (SP==0)):
                # here we need an additional time inversion operation, so read it:
                for orb in range(n_orbits):
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat_tinv[orb][i,j] = R.next()            # real part
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat_tinv[orb][i,j] += 1j * R.next()      # imaginary part
                


        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Wien2k_converter : reading file symm_file failed!"
        
        R.close()

        # Save it to the HDF:
        ar=HDFArchive(self.hdf_file,'a')
        if not (symm_subgrp in ar): ar.create_group(symm_subgrp)
        thingstowrite = ['n_s','n_atoms','perm','orbits','SO','SP','time_inv','mat','mat_tinv']
        for it in thingstowrite: exec "ar['%s']['%s'] = %s"%(symm_subgrp,it,it)
        del ar
        
        

    def __repack(self):
        """Calls the h5repack routine, in order to reduce the file size of the hdf5 archive.
           Should only be used BEFORE the first invokation of HDFArchive in the program, otherwise
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
                                
