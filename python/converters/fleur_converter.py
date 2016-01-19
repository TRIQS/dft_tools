
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2014 by M. Betzinger, P. Seth
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
from converter_tools import *

class FleurConverter(ConverterTools):
    """
    Conversion from Fleur output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, filename, hdf_filename = None,
                       dft_subgrp = 'dft_input', parproj_subgrp='dft_parproj_input',
                       repacking = False):
        """
        Init of the class. Variable filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so on. 
        """

        assert type(filename)==StringType, "Please provide the DFT files' base name as a string."
        if hdf_filename is None: hdf_filename = filename
        self.hdf_file = hdf_filename+'.h5'
        #self.dft_file = filename+'.ctqmcout'
        #self.parproj_file = filename+'.parproj'
        self.dft_file = 'ctqmcout'
        self.parproj_file = 'parproj'
        self.dft_subgrp = dft_subgrp
        self.parproj_subgrp = parproj_subgrp
        self.fortran_to_replace = {'D':'E'}

        # Checks if h5 file is there and repacks it if wanted:
        import os.path
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)
        

    def convert_dft_input(self):
        """
        Reads the input files, and stores the data in the HDFfile
        """
        
        # Read and write only on the master node
        if not (mpi.is_master_node()): return
        mpi.report("Reading input from %s..."%self.dft_file)

        # R is a generator : each R.Next() will return the next number in the file
        R = ConverterTools.read_fortran_file(self,self.dft_file,self.fortran_to_replace)
        try:
            energy_unit = R.next()                        # read the energy convertion factor
            n_k = int(R.next())                           # read the number of k points
            k_dep_projection = 1                          
            SP = int(R.next())                            # flag for spin-polarised calculation
            SO = int(R.next())                            # flag for spin-orbit calculation
            charge_below = R.next()                       # total charge below energy window
            density_required = R.next()                   # total density required, for setting the chemical potential
            symm_op = 0                                   # Use symmetry groups for the k-sum

            # the information on the non-correlated shells is not important here, maybe skip:
            n_shells = int(R.next())                      # number of shells (e.g. Fe d, As p, O p) in the unit cell, 
                                                          # corresponds to index R in formulas
            # now read the information about the shells (atom, sort, l, dim):
            shell_entries = ['atom', 'sort', 'l', 'dim']
            shells = [ {name: int(val) for name, val in zip(shell_entries, R)} for ish in range(n_shells) ]

            n_corr_shells = int(R.next())                 # number of corr. shells (e.g. Fe d, Ce f) in the unit cell, 
                                                          # corresponds to index R in formulas
            # now read the information about the shells (atom, sort, l, dim, SO flag, irep):
            corr_shell_entries = ['atom', 'sort', 'l', 'dim', 'SO', 'irep']
            corr_shells = [ {name: int(val) for name, val in zip(corr_shell_entries, R)} for icrsh in range(n_corr_shells) ]

            # determine the number of inequivalent correlated shells and maps, needed for further reading
            n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = ConverterTools.det_shell_equivalence(self,corr_shells)

            use_rotations = 1
            rot_mat = [numpy.identity(corr_shells[icrsh]['dim'],numpy.complex_) for icrsh in range(n_corr_shells)]
           
            # read the matrices
            rot_mat_time_inv = [0 for i in range(n_corr_shells)]

            for icrsh in range(n_corr_shells):
                for i in range(corr_shells[icrsh]['dim']):    # read real part:
                    for j in range(corr_shells[icrsh]['dim']):
                        rot_mat[icrsh][i,j] = R.next()
                for i in range(corr_shells[icrsh]['dim']):    # read imaginary part:
                    for j in range(corr_shells[icrsh]['dim']):
                        rot_mat[icrsh][i,j] += 1j * R.next()

                if (SP==1):             # read time inversion flag:
                    rot_mat_time_inv[icrsh] = int(R.next())
                    
            # Read here the info for the transformation of the basis:
            n_reps = [1 for i in range(n_inequiv_shells)]
            dim_reps = [0 for i in range(n_inequiv_shells)]
            T = []
            for ish in range(n_inequiv_shells):
                n_reps[ish] = int(R.next())   # number of representatives ("subsets"), e.g. t2g and eg
                dim_reps[ish] = [int(R.next()) for i in range(n_reps[ish])]   # dimensions of the subsets
            
                # The transformation matrix:
                # is of dimension 2l+1 without SO, and 2*(2l+1) with SO!
                ll = 2*corr_shells[inequiv_to_corr[ish]]['l']+1
                lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
                T.append(numpy.zeros([lmax,lmax],numpy.complex_))
                
                # now read it from file:
                for i in range(lmax):
                    for j in range(lmax):
                        T[ish][i,j] = R.next()
                for i in range(lmax):
                    for j in range(lmax):
                        T[ish][i,j] += 1j * R.next()
    
            # Spin blocks to be read:
            n_spin_blocs = SP + 1 - SO   
                 
            # read the list of n_orbitals for all k points
            n_orbitals = numpy.zeros([n_k,n_spin_blocs],numpy.int)
            for isp in range(n_spin_blocs):
                for ik in range(n_k):
                    n_orbitals[ik,isp] = int(R.next())
            
            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k,n_spin_blocs,n_corr_shells,max([crsh['dim'] for crsh in corr_shells]),max(n_orbitals)],numpy.complex_)

            # Read the projectors from the file:
            for ik in range(n_k):
                for icrsh in range(n_corr_shells):
                    n_orb = corr_shells[icrsh]['dim']
                    # first Real part for BOTH spins, due to conventions in dmftproj:
                    for isp in range(n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik][isp]):
                                proj_mat[ik,isp,icrsh,i,j] = R.next()
                    # now Imag part:
                    for isp in range(n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik][isp]):
                                proj_mat[ik,isp,icrsh,i,j] += 1j * R.next()
          
            # now define the arrays for weights and hopping ...
            bz_weights = numpy.ones([n_k],numpy.float_)/ float(n_k)  # w(k_index),  default normalisation 
            hopping = numpy.zeros([n_k,n_spin_blocs,max(n_orbitals),max(n_orbitals)],numpy.complex_)

            # weights in the file
            for ik in range(n_k) : bz_weights[ik] = R.next()         
                
            # if the sum over spins is in the weights, take it out again!!
            sm = sum(bz_weights)
            bz_weights[:] /= sm 

            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian -- convention for Fleur
            for isp in range(n_spin_blocs):
                for ik in range(n_k) :
                    n_orb = n_orbitals[ik,isp]
                    for i in range(n_orb):
                        hopping[ik,isp,i,i] = R.next() * energy_unit
            
            # keep some things that we need for reading parproj:
            things_to_set = ['n_shells','shells','n_corr_shells','corr_shells','n_spin_blocs','n_orbitals','n_k','SO','SP','energy_unit'] 
            for it in things_to_set: setattr(self,it,locals()[it])
        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Fleur_converter : reading file %s failed!"%filename

        R.close()
        # Reading done!
        
        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.dft_subgrp in ar): ar.create_group(self.dft_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created. If it exists, the data is overwritten!
        things_to_save = ['energy_unit','n_k','k_dep_projection','SP','SO','charge_below','density_required',
                          'symm_op','n_shells','shells','n_corr_shells','corr_shells','use_rotations','rot_mat',
                          'rot_mat_time_inv','n_reps','dim_reps','T','n_orbitals','proj_mat','bz_weights','hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
        for it in things_to_save: ar[self.dft_subgrp][it] = locals()[it]
        del ar


    def convert_parproj_input(self):
        """
        Reads the input for the partial charges projectors from case.parproj, and stores it in the symmpar_subgrp
        group in the HDF5.
        """

        if not (mpi.is_master_node()): return
        mpi.report("Reading parproj input from %s..."%self.parproj_file)

        dens_mat_below = [ [numpy.zeros([self.shells[ish]['dim'],self.shells[ish]['dim']],numpy.complex_) for ish in range(self.n_shells)] 
                           for isp in range(self.n_spin_blocs) ]

        R = ConverterTools.read_fortran_file(self,self.parproj_file,self.fortran_to_replace)

        n_parproj = [int(R.next()) for i in range(self.n_shells)]
        n_parproj = numpy.array(n_parproj)
                
        # Initialise P, here a double list of matrices:
        proj_mat_pc = numpy.zeros([self.n_k,self.n_spin_blocs,self.n_shells,max(n_parproj),max([sh['dim'] for sh in self.shells]),max(self.n_orbitals)],numpy.complex_)
        
        rot_mat_all = [numpy.identity(self.shells[ish]['dim'],numpy.complex_) for ish in range(self.n_shells)]
        rot_mat_all_time_inv = [0 for i in range(self.n_shells)]

        for ish in range(self.n_shells):
            # read first the projectors for this orbital:
            for ik in range(self.n_k):
                for ir in range(n_parproj[ish]):

                    for isp in range(self.n_spin_blocs):
                        for i in range(self.shells[ish]['dim']):    # read real part:
                            for j in range(self.n_orbitals[ik][isp]):
                                proj_mat_pc[ik,isp,ish,ir,i,j] = R.next()
                            
                    for isp in range(self.n_spin_blocs):
                        for i in range(self.shells[ish]['dim']):    # read imaginary part:
                            for j in range(self.n_orbitals[ik][isp]):
                                proj_mat_pc[ik,isp,ish,ir,i,j] += 1j * R.next()
                                        
                    
            # now read the Density Matrix for this orbital below the energy window:
            for isp in range(self.n_spin_blocs):
                for i in range(self.shells[ish]['dim']):    # read real part:
                    for j in range(self.shells[ish]['dim']):
                        dens_mat_below[isp][ish][i,j] = R.next()
            for isp in range(self.n_spin_blocs):
                for i in range(self.shells[ish]['dim']):    # read imaginary part:
                    for j in range(self.shells[ish]['dim']):
                        dens_mat_below[isp][ish][i,j] += 1j * R.next()
                if (self.SP==0): dens_mat_below[isp][ish] /= 2.0

            # Global -> local rotation matrix for this shell:
            for i in range(self.shells[ish]['dim']):    # read real part:
                for j in range(self.shells[ish]['dim']):
                    rot_mat_all[ish][i,j] = R.next()
            for i in range(self.shells[ish]['dim']):    # read imaginary part:
                for j in range(self.shells[ish]['dim']):
                    rot_mat_all[ish][i,j] += 1j * R.next()
                    
            if (self.SP):
                rot_mat_all_time_inv[ish] = int(R.next())

        R.close()
        # Reading done!

        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.parproj_subgrp in ar): ar.create_group(self.parproj_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created. If it exists, the data is overwritten!
        things_to_save = ['dens_mat_below','n_parproj','proj_mat_pc','rot_mat_all','rot_mat_all_time_inv']
        for it in things_to_save: ar[self.parproj_subgrp][it] = locals()[it]
        del ar
