
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
from converter_tools import *
import os.path

class Wien2kConverter(ConverterTools):
    """
    Conversion from Wien2k output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, filename, hdf_filename = None,
                       dft_subgrp = 'dft_input', symmcorr_subgrp = 'dft_symmcorr_input',
                       parproj_subgrp='dft_parproj_input', symmpar_subgrp='dft_symmpar_input',
                       bands_subgrp = 'dft_bands_input', transp_subgrp = 'dft_transp_input', repacking = False):
        """
        Init of the class. Variable filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so on. 
        """

        assert type(filename)==StringType, "Please provide the DFT files' base name as a string."
        if hdf_filename is None: hdf_filename = filename
        self.hdf_file = hdf_filename+'.h5'
        self.dft_file = filename+'.ctqmcout'
        self.symmcorr_file = filename+'.symqmc'
        self.parproj_file = filename+'.parproj'
        self.symmpar_file = filename+'.sympar'
        self.band_file = filename+'.outband'
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.parproj_subgrp = parproj_subgrp
        self.symmpar_subgrp = symmpar_subgrp
        self.bands_subgrp = bands_subgrp
        self.transp_subgrp = transp_subgrp
        self.fortran_to_replace = {'D':'E'}
        self.vel_file = filename+'.pmat'
        self.outputs_file = filename+'.outputs'
        self.struct_file = filename+'.struct'
        self.oubwin_file = filename+'.oubwin'

        # Checks if h5 file is there and repacks it if wanted:
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
            symm_op = 1                                   # Use symmetry groups for the k-sum

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
            # we use now the convention of a DIAGONAL Hamiltonian -- convention for Wien2K.
            for isp in range(n_spin_blocs):
                for ik in range(n_k) :
                    n_orb = n_orbitals[ik,isp]
                    for i in range(n_orb):
                        hopping[ik,isp,i,i] = R.next() * energy_unit
            
            # keep some things that we need for reading parproj:
            things_to_set = ['n_shells','shells','n_corr_shells','corr_shells','n_spin_blocs','n_orbitals','n_k','SO','SP','energy_unit'] 
            for it in things_to_set: setattr(self,it,locals()[it])
        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Wien2k_converter : reading file %s failed!"%filename

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

        # Symmetries are used, so now convert symmetry information for *correlated* orbitals:
        self.convert_symmetry_input(orbits=corr_shells,symm_file=self.symmcorr_file,symm_subgrp=self.symmcorr_subgrp,SO=self.SO,SP=self.SP)


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

        # Symmetries are used, so now convert symmetry information for *all* orbitals:
        self.convert_symmetry_input(orbits=self.shells,symm_file=self.symmpar_file,symm_subgrp=self.symmpar_subgrp,SO=self.SO,SP=self.SP)


    def convert_bands_input(self):
        """
        Converts the input for momentum resolved spectral functions, and stores it in bands_subgrp in the
        HDF5.
        """

        if not (mpi.is_master_node()): return
        mpi.report("Reading bands input from %s..."%self.band_file)

        R = ConverterTools.read_fortran_file(self,self.band_file,self.fortran_to_replace)
        try:
            n_k = int(R.next())

            # read the list of n_orbitals for all k points
            n_orbitals = numpy.zeros([n_k,self.n_spin_blocs],numpy.int)
            for isp in range(self.n_spin_blocs):
                for ik in range(n_k):
                    n_orbitals[ik,isp] = int(R.next())

            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k,self.n_spin_blocs,self.n_corr_shells,max([crsh['dim'] for crsh in corr_shells]),max(n_orbitals)],numpy.complex_)

            # Read the projectors from the file:
            for ik in range(n_k):
                for icrsh in range(self.n_corr_shells):
                    n_orb = self.corr_shells[icrsh]['dim']
                    # first Real part for BOTH spins, due to conventions in dmftproj:
                    for isp in range(self.n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik,isp]):
                                proj_mat[ik,isp,icrsh,i,j] = R.next()
                    # now Imag part:
                    for isp in range(self.n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik,isp]):
                                proj_mat[ik,isp,icrsh,i,j] += 1j * R.next()

            hopping = numpy.zeros([n_k,self.n_spin_blocs,max(n_orbitals),max(n_orbitals)],numpy.complex_)
         	    
            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian!!!!
            for isp in range(self.n_spin_blocs):
                for ik in range(n_k) :
                    n_orb = n_orbitals[ik,isp]
                    for i in range(n_orb):
                        hopping[ik,isp,i,i] = R.next() * self.energy_unit

            # now read the partial projectors:
            n_parproj = [int(R.next()) for i in range(self.n_shells)]
            n_parproj = numpy.array(n_parproj)
            
            # Initialise P, here a double list of matrices:
            proj_mat_pc = numpy.zeros([n_k,self.n_spin_blocs,self.n_shells,max(n_parproj),max([sh['dim'] for sh in self.shells]),max(n_orbitals)],numpy.complex_)

            for ish in range(self.n_shells):
                for ik in range(n_k):
                    for ir in range(n_parproj[ish]):
                        for isp in range(self.n_spin_blocs):
                                    
                            for i in range(self.shells[ish]['dim']):    # read real part:
                                for j in range(n_orbitals[ik,isp]):
                                    proj_mat_pc[ik,isp,ish,ir,i,j] = R.next()
                            
                            for i in range(self.shells[ish]['dim']):    # read imaginary part:
                                for j in range(n_orbitals[ik,isp]):
                                    proj_mat_pc[ik,isp,ish,ir,i,j] += 1j * R.next()

        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Wien2k_converter : reading file band_file failed!"

        R.close()
        # Reading done!

        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.bands_subgrp in ar): ar.create_group(self.bands_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created. If it exists, the data is overwritten!
        things_to_save = ['n_k','n_orbitals','proj_mat','hopping','n_parproj','proj_mat_pc']
        for it in things_to_save: ar[self.bands_subgrp][it] = locals()[it]
        del ar


    def convert_transport_input(self, spinbl=['']):
        """ 
        Reads the input files necessary for transport calculations
        and stores the data in the HDFfile
        """

        #Read and write files only on the master node
        if not (mpi.is_master_node()): return
        
        # Check if SP, SO and n_k are already in h5
        ar = HDFArchive(self.hdf_file, 'a')
        if not (self.dft_subgrp in ar): raise IOError, "No %s subgroup in hdf file found! Call convert_dmft_input first." %self.dft_subgrp
        SP = ar[self.dft_subgrp]['SP']
        SO = ar[self.dft_subgrp]['SO']
        n_k = ar[self.dft_subgrp]['n_k']
        del ar

        # Read relevant data from .pmat file
        ############################################
       
        vk = [[] for sp in spinbl]
        kp = [[] for sp in spinbl]
        bandwin_opt = []
        
        for isp, sp in enumerate(spinbl):
            bandwin_opt_isp = []
            if not (os.path.exists(self.vel_file + sp)) : raise IOError, "File %s does not exist" %self.vel_file+sp
            print "Reading input from %s..."%self.vel_file+sp

            with open(self.vel_file + sp) as f:
                while 1:
                    try:
                        s = f.readline()
                        if (s == ''):
                            break
                    except:
                        break
                    try:
                       [k, nu1, nu2] = [int(x) for x in s.strip().split()]
                       bandwin_opt_isp.append((nu1,nu2))
                       dim = nu2 - nu1 +1
                       v_xyz = numpy.zeros((dim,dim,3), dtype = complex)
                       temp = f.readline().strip().split()
                       kp[isp].append(numpy.array([float(t) for t in temp[0:3]]))
                       for nu_i in xrange(dim):
                           for nu_j in xrange(nu_i, dim):
                               for i in xrange(3):
                                   s = f.readline().strip("\n ()").split(',')
                                   v_xyz[nu_i][nu_j][i] = float(s[0]) + float(s[1])*1j
                                   if (nu_i != nu_j):
                                        v_xyz[nu_j][nu_i][i] = v_xyz[nu_i][nu_j][i].conjugate()

                       vk[isp].append(v_xyz)
        
                    except IOError:
                        raise "Wien2k_converter : reading file %s failed" %self.vel_file
            bandwin_opt.append(numpy.array(bandwin_opt_isp))

        print "Read in %s file done!" %self.vel_file


        # Read relevant data from .struct file
        ############################################
        if not (os.path.exists(self.struct_file)) : raise IOError, "File %s does not exist" %self.struct_file
        print "Reading input from %s..."%self.struct_file
        
        with open(self.struct_file) as f:
            try:
                f.readline() #title
                temp = f.readline() #lattice
                #latticetype = temp[0:10].split()[0]
                latticetype = temp.split()[0]
                f.readline()
                temp = f.readline().strip().split() # lattice constants
                latticeconstants = numpy.array([float(t) for t in temp[0:3]])
                latticeangles = numpy.array([float(t) for t in temp[3:6]])
                latticeangles *= numpy.pi/180.0
                print 'Lattice: ', latticetype
                print 'Lattice constants: ', latticeconstants
                print 'Lattice angles: ', latticeangles

            except IOError:
                raise "Wien2k_converter : reading file %s failed" %self.struct_file

        print "Read in %s file done!" %self.struct_file


        # Read relevant data from .outputs file
        ############################################
        if not (os.path.exists(self.outputs_file)) : raise IOError, "File %s does not exist" %self.outputs_file
        print "Reading input from %s..."%self.outputs_file
        
        symmcartesian = []
        taucartesian = []

        with open(self.outputs_file) as f:
            try:
                while 1:
                    temp = f.readline().strip(' ').split()
                    if (temp[0] =='PGBSYM:'):
                         nsymm = int(temp[-1])
                         break
                for i in range(nsymm):
                    while 1:
                        temp = f.readline().strip().split()
                        if (temp[0] == 'Symmetry'):
                            break

                    # read cartesian symmetries
                    symmt = numpy.zeros((3, 3), dtype = float)
                    taut = numpy.zeros(3, dtype = float)
                    for ir in range(3):
                        temp = f.readline().strip().split()
                        for ic in range(3):
                            symmt[ir, ic] = float(temp[ic])
                    temp = f.readline().strip().split()
                    for ir in range(3):
                        taut[ir] = float(temp[ir])

                    symmcartesian.append(symmt)
                    taucartesian.append(taut)
            except IOError:
                 raise "Wien2k_converter: reading file %s failed" %self.outputs_file
     
        print "Read in %s file done!" %self.outputs_file


        # Read relevant data from .oubwin/up/down files
        ############################################
        
        if (SP == 0 or SO == 1):        
            files = [self.oubwin_file]
        elif SP == 1:
            files = [self.oubwin_file+'up', self.oubwin_file+'dn']
        else: # SO and SP can't both be 1
            assert 0, "Reding oubwin error! Check SP and SO!"
        bandwin = [numpy.zeros((n_k, 2), dtype=int) for isp in range(SP + 1 - SO)]

        for isp, f in enumerate(files):
            if not os.path.exists(f): raise IOError, "File %s does not exist" %f
            print "Reading input from %s..."%f
            R = ConverterTools.read_fortran_file(self, f, self.fortran_to_replace)
            assert int(R.next()) == n_k, "Wien2k_converter: Number of k-points is inconsistent in oubwin file!"
            assert int(R.next()) == SO, "Wien2k_converter: SO is inconsistent in oubwin file!"

            for ik in xrange(n_k):
                R.next()
                bandwin[isp][ik,0] = R.next()
                bandwin[isp][ik,1] = R.next()
                R.next()

        print "Read in %s files done!" %self.oubwin_file

        # Put data to HDF5 file
        ar = HDFArchive(self.hdf_file, 'a')
        if not (self.transp_subgrp in ar): ar.create_group(self.transp_subgrp)
        # The subgroup containing the data. If it does not exist, it is created. If it exists, the data is overwritten!!!
        things_to_save = ['bandwin_opt', 'kp', 'vk', 'latticetype', 'latticeconstants', 'latticeangles', 'nsymm', 'symmcartesian',
                        'taucartesian', 'bandwin']
        for it in things_to_save: ar[self.transp_subgrp][it] = locals()[it]
        del ar

    def convert_symmetry_input(self, orbits, symm_file, symm_subgrp, SO, SP):
        """
        Reads input for the symmetrisations from symm_file, which is case.sympar or case.symqmc.
        """

        if not (mpi.is_master_node()): return
        mpi.report("Reading symmetry input from %s..."%symm_file)

        n_orbits = len(orbits)

        R = ConverterTools.read_fortran_file(self,symm_file,self.fortran_to_replace)

        try:
            n_symm = int(R.next())           # Number of symmetry operations
            n_atoms = int(R.next())       # number of atoms involved
            perm = [ [int(R.next()) for i in range(n_atoms)] for j in range(n_symm) ]    # list of permutations of the atoms
            if SP: 
                time_inv = [ int(R.next()) for j in range(n_symm) ]           # time inversion for SO coupling
            else:
                time_inv = [ 0 for j in range(n_symm) ]

            # Now read matrices:
            mat = []  
            for i_symm in range(n_symm):
                
                mat.append( [ numpy.zeros([orbits[orb]['dim'], orbits[orb]['dim']],numpy.complex_) for orb in range(n_orbits) ] )
                for orb in range(n_orbits):
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            mat[i_symm][orb][i,j] = R.next()            # real part
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            mat[i_symm][orb][i,j] += 1j * R.next()      # imaginary part

            mat_tinv = [numpy.identity(orbits[orb]['dim'],numpy.complex_)
                        for orb in range(n_orbits)]

            if ((SO==0) and (SP==0)):
                # here we need an additional time inversion operation, so read it:
                for orb in range(n_orbits):
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            mat_tinv[orb][i,j] = R.next()            # real part
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            mat_tinv[orb][i,j] += 1j * R.next()      # imaginary part
                


        except StopIteration : # a more explicit error if the file is corrupted.
            raise "Wien2k_converter : reading file symm_file failed!"
        
        R.close()
        # Reading done!

        # Save it to the HDF:
        ar=HDFArchive(self.hdf_file,'a')
        if not (symm_subgrp in ar): ar.create_group(symm_subgrp)
        things_to_save = ['n_symm','n_atoms','perm','orbits','SO','SP','time_inv','mat','mat_tinv']
        for it in things_to_save: ar[symm_subgrp][it] = locals()[it]
        del ar
