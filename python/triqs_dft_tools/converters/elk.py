
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2019 by A. D. N. James, A. Hampel and M. Aichhorn
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
"""
Elk converter
"""

from types import *
import numpy
from .converter_tools import *
import os.path
from h5 import *
from locale import atof
from triqs_dft_tools.converters.elktools import readElkfiles as read_Elk
from triqs_dft_tools.converters.elktools import ElkConverterTools as Elk_tools

class ElkConverter(ConverterTools,Elk_tools,read_Elk):
    """
    Conversion from Elk output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, filename, hdf_filename=None,
                 dft_subgrp='dft_input', symmcorr_subgrp='dft_symmcorr_input',
                 bc_subgrp='dft_bandchar_input', symmpar_subgrp='dft_symmpar_input',
                 bands_subgrp='dft_bands_input', misc_subgrp='dft_misc_input',
                 transp_subgrp='dft_transp_input',cont_subgrp='dft_contours_input',
                 repacking=False):
        """
        Initialise the class.

        Parameters
        ----------
        filename : string
                   Base name of DFT files.
        hdf_filename : string, optional
                       Name of hdf5 archive to be created.
        dft_subgrp : string, optional
                     Name of subgroup storing necessary DFT data.
        symmcorr_subgrp : string, optional
                          Name of subgroup storing correlated-shell symmetry data.
        parproj_subgrp : string, optional
                         Name of subgroup storing partial projector data.
        symmpar_subgrp : string, optional
                         Name of subgroup storing partial-projector symmetry data.
        bands_subgrp : string, optional
                       Name of subgroup storing band data.
        misc_subgrp : string, optional
                      Name of subgroup storing miscellaneous DFT data.
        transp_subgrp : string, optional
                        Name of subgroup storing transport data.
        repacking : boolean, optional
                    Does the hdf5 archive need to be repacked to save space?

        """

        assert isinstance(filename, str), "ElkConverter: Please provide the DFT files' base name as a string."
        if hdf_filename is None:
            hdf_filename = filename + '.h5'
        self.hdf_file = hdf_filename
        self.dft_file = 'PROJ.OUT'
        self.band_file = 'BAND.OUT'
        self.eval_file = 'EIGVAL.OUT'
        self.efermi_file = 'EFERMI.OUT'
        self.kp_file = 'KPOINTS.OUT'
        self.geom_file='GEOMETRY.OUT'
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.bc_subgrp = bc_subgrp
        self.symmpar_subgrp = symmpar_subgrp
        self.bands_subgrp = bands_subgrp
        self.misc_subgrp = misc_subgrp
        self.transp_subgrp = transp_subgrp
        self.cont_subgrp = cont_subgrp
        self.fortran_to_replace = {'D': 'E'}

        # Checks if h5 file is there and repacks it if wanted:
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)

    def check_dens(self,n_k,nstsv,occ,bz_weights,n_spin_blocs,band_window,SO):
      """
      Check the charge density below the correlated energy window and up to the Fermi level
      """

      density_required=0.0E-7
      charge_below=0.0E-7
      #calculate the valence charge and charge below lower energy window bound
      #Elk does not use the tetrahedron method when calculating these charges
      for ik in range(n_k):
        for ist in range(nstsv):
          #calculate the charge over all the bands
          density_required+=occ[ik][ist]*bz_weights[ik]
        for isp in range(n_spin_blocs):
          #Convert occ list from elk to two index format for spins
          jst=int((isp)*nstsv*0.5)
          #Take lowest index in band_window for SO system
          if(SO==0):
            nst=band_window[isp][ik, 0]-1
          else:
            band=[band_window[0][ik, 0],band_window[1][ik, 0]]
            nst=min(band)-1
          #calculate the charge below energy window
          for ist in range(jst,nst):
            charge_below+=occ[ik][ist]*bz_weights[ik]
      #return charges
      return(density_required,charge_below)

    def rotsym(self,n_shells,shells,n_symm,ind,basis,T,mat):
      """
      Rotates the symmetry matrices into basis defined by the T unitary matrix
      the outputted projectors are rotated to the irreducible representation
      and then reduced in size to the orbitals used to construct the projectors.
      """

      for ish in range(n_shells):
           #check that the T matrix is not the Identity (i.e. not using spherical
           #harmonics).
           if(basis[ish]!=0):
           #put mat into temporary matrix
             temp=mat
           #index range of lm values used to create the Wannier projectors
             min_ind=numpy.min(ind[ish][:])
             max_ind=numpy.max(ind[ish][:])+1
           #dimension of lm values used to construct the projectors
             dim=shells[ish]['dim']
           #loop over all symmetries
             for isym in range(n_symm):
           #rotate symmetry matrix into basis defined by T
               mat[isym][ish]=numpy.matmul(T[ish],mat[isym][ish])
               mat[isym][ish]=numpy.matmul(mat[isym][ish],T[ish].conjugate().transpose())
           #put desired subset of transformed symmetry matrix into temp matrix for symmetry isym
               for id in range(len(ind[ish])):
                 i=ind[ish][id]
                 for jd in range(len(ind[ish][:])):
                   j=ind[ish][jd]
                   temp[isym][ish][id,jd]=mat[isym][ish][i,j]
           #put temp matrix into mat
             mat=temp
           #reduce size of lm arrays in mat lm dim
             for isym in range(n_symm):
               dim=shells[ish]['dim']
               mat[isym][ish]=mat[isym][ish][:dim,:dim]
      return mat

    def update_so_quatities(self,n_shells,shells,n_corr_shells,corr_shells,n_inequiv_shells,dim_reps,n_k,n_symm,n_orbitals,proj_mat,T,su2,mat,sym=True):
        """
        Changes the array sizes and elements for arrays used in spin-orbit coupled calculations.
        """

        #change dim for each shell
        for ish in range(n_shells):
          shells[ish]['dim'] *= 2
        for ish in range(n_corr_shells):
          corr_shells[ish]['dim'] *= 2
        for ish in range(n_inequiv_shells):
          dim_reps[ish]=[2*dim_reps[ish][i] for i in range(len(dim_reps[ish]))]
        #Make temporary array of original n_orbitals
        n_orbitals_orig=n_orbitals
        #Make SO n_orbitals array
        #loop over k-points
        for ik in range(n_k):
          #new orbital array
          n_orbitals[ik,0]=max(n_orbitals[ik,:])
        #reduce array size
        n_orbitals=n_orbitals[:,:1]
        #Resize proj_mat, mat, T
        #make temporary projector array
        proj_mat_tmp = numpy.zeros([n_k, 1, n_corr_shells, max([crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], complex)
        for ish in range(n_corr_shells):
          #update proj_mat
          for ik in range(n_k):
            #extra array elements in "dim" dimension
            size=int(0.5*corr_shells[ish]['dim'])
            #put each spinor into tmp array and ensure elements are assigned correctly in case of change of max(n_orbitals)
            proj_mat_tmp[ik][0][ish][0:size][0:n_orbitals_orig[ik,0]]=proj_mat[ik][0][ish][0:size][0:n_orbitals_orig[ik,0]]
            #put other spinor projectors into extra "dim" elements
            proj_mat_tmp[ik][0][ish][size:2*size][0:n_orbitals_orig[ik,1]]=proj_mat[ik][1][ish][0:size][0:n_orbitals_orig[ik,1]]
          #update T
          #extra array elements in each dimension
          size=2*corr_shells[ish]['l']+1
          #extend the arrays
          T[ish]=numpy.lib.pad(T[ish],((0,size),(0,size)),'constant',constant_values=(0.0))
          #make block diagonal
          T[ish][size:2*size,size:2*size]=T[ish][0:size,0:size]
          #update the symmetries arrays if needed
          if(sym):
            #update mat - This includes the spin SU(2) matrix for spin-coupled calculations
            #size of each quadrant in the lm symmetry array.
            size=int(0.5*corr_shells[ish]['dim'])
            #temporary spin block array for SU(2) spin operations on mat
            spinmat = numpy.zeros([size,2,size,2],complex)
            for isym in range(n_symm):
              #expand size of array
              mat[isym][ish]=numpy.lib.pad(mat[isym][ish],((0,size),(0,size)),'constant',constant_values=(0.0))
              #make arraye block diagonal
              mat[isym][ish][size:2*size,size:2*size]=mat[isym][ish][0:size,0:size]
              #apply SU(2) spin matrices to lm symmetries
              #put mat into array of spin blocks
              for i1,i2 in numpy.ndindex(2,2):
                spinmat[0:size,i1,0:size,i2] = mat[isym][ish][i1*size:(i1+1)*size,i2*size:(i2+1)*size]
              #apply the SU(2) spin matrices
              for ilm,jlm in numpy.ndindex(size,size):
                spinmat[ilm,:,jlm,:] = numpy.dot(su2[isym][:,:],spinmat[ilm,:,jlm,:])
              #put spinmat into back in mat format
              for i1,i2 in numpy.ndindex(2,2):
                mat[isym][ish][i1*size:(i1+1)*size,i2*size:(i2+1)*size] = spinmat[0:size,i1,0:size,i2]
        #assign arrays and delete temporary arrays
        del proj_mat
        proj_mat = proj_mat_tmp
        del proj_mat_tmp, n_orbitals_orig
        return shells,corr_shells,dim_reps,n_orbitals,proj_mat,T,mat

    def sort_dft_eigvalues(self,n_spin_blocs,SO,n_k,n_orbitals,band_window,en,energy_unit):
        """
        Rearranges the energy eigenvalue arrays into TRIQS format
        """

        hopping = numpy.zeros([n_k, n_spin_blocs, numpy.max(n_orbitals), numpy.max(n_orbitals)], complex)
        #loop over spin
        for isp in range(n_spin_blocs):
          #loop over k-points
          for ik in range(n_k):
            #loop over bands
            for ist in range(0,n_orbitals[ik, isp]):
              #converter index for spin polarised Elk indices and take SO into consideration
              if(SO==0):
                jst=int(band_window[isp][ik, 0]-1+ist)
              else:
                band=[band_window[0][ik, 0],band_window[1][ik, 0]]
                jst=int(min(band)-1+ist)
              #correlated window energies
              hopping[ik,isp,ist,ist]=en[ik][jst]*energy_unit
        return hopping

    def convert_dft_input(self):
        """
        Reads the appropriate files and stores the data for the

        - dft_subgrp
        - symmcorr_subgrp
        - misc_subgrp

        in the hdf5 archive.

        """
        # Read and write only on the master node
        if not (mpi.is_master_node()):
            return
        filext='.OUT'
        dft_file='PROJ'+filext
        mpi.report("Reading %s" % dft_file)
        #Energy conversion - Elk uses Hartrees
        energy_unit = 27.2113850560                        # Elk uses hartrees
        #The projectors change size per k-point
        k_dep_projection = 1
        #Symmetries are used
        symm_op = 1                                   # Use symmetry groups for the k-sum
        shells=[]
        #read information about projectors calculated in the Elk calculation
        [gen_info,n_corr_shells,n_inequiv_shells,corr_to_inequiv,inequiv_to_corr,corr_shells,n_reps,dim_reps,ind,basis,T] = read_Elk.read_proj(self,dft_file)
        #get info for HDF5 file from gen_info
        n_k=gen_info['n_k']
        SP=gen_info['spinpol']-1
        #Elk uses spinor wavefunctions. Therefore these two spinor wavefunctions have spin-orbit coupling incorporated in them. Here we read in the spinors
        n_spin_blocs = SP + 1
        SO=gen_info['SO']
        n_atoms=gen_info['natm']
        #Elk only calculates Wannier projectors (no theta projectors generated):
        n_shells=n_corr_shells
        for ish in range(n_shells):
           shells.append(corr_shells[ish].copy())
           #remove last 2 entries from corr_shlls
           del shells[ish]['SO']
           del shells[ish]['irep']
           shells[ish]['dim'] = 2*shells[ish]['l']+1
        #read eigenvalues calculated in the Elk calculation
        mpi.report("Reading %s and EFERMI.OUT" % self.eval_file)
        [en,occ,nstsv]=read_Elk.read_eig(self)
        #read projectors calculated in the Elk calculation
        proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max([crsh['dim'] for crsh in corr_shells]), nstsv], complex)
        mpi.report("Reading projector(s)")
        for ish in range(n_corr_shells):
          [n_orbitals,band_window,rep,proj_mat]=read_Elk.read_projector(self,corr_shells,n_spin_blocs,ish,proj_mat,ind,T,basis,filext)
        #read kpoints calculated in the Elk calculation
        mpi.report("Reading %s" % self.kp_file)
        [bz_weights,vkl]=read_Elk.read_kpoints(self)
        #symmetry matrix
        mpi.report("Reading GEOMETRY.OUT")
        #read in atom positions, the symmetry operators (in lattice coordinates) and lattice vectors
        [ns, na, atpos]=read_Elk.read_geometry(self)
        #Read symmetry files
        mpi.report("Reading SYMCRYS.OUT")
        [n_symm,spinmat,symmat,tr] = read_Elk.readsym(self)
        mpi.report("Reading LATTICE.OUT")
        [amat,amatinv,bmat,bmatinv,cell_vol] = read_Elk.readlat(self)
        #calculating atom permutations
        perm = Elk_tools.gen_perm(self,n_symm,ns,na,n_atoms,symmat,tr,atpos)
        #determine the cartesian lattice symmetries and the spin axis rotations
        #required for the spinors (for SO for now)
        su2 = []
        symmatc=[]
        for isym in range(n_symm):
          #convert the lattice symmetry matrices into cartesian coordinates
          tmp = numpy.matmul(amat,symmat[isym])
          symmatc.append(numpy.matmul(tmp,amatinv))
          #convert the spin symmetry matrices into cartesian coordinates
          spinmatc = numpy.matmul(amat,spinmat[isym])
          spinmatc = numpy.matmul(spinmatc,amatinv)
          #calculate the rotation angle and spin axis vector in cartesian coordinates
          [v,th] = self.rotaxang(spinmatc[:,:])
          #calculate the SU(2) matrix from the angle and spin axis vector
          su2.append(self.axangsu2(v,th))
        del tmp

        #calculating the symmetries in complex harmonics
        mat = Elk_tools.symlat_to_complex_harmonics(self,n_symm,n_corr_shells,symmatc,corr_shells)
        mat = self.rotsym(n_corr_shells,corr_shells,n_symm,ind,basis,T,mat)

        #The reading is done. Some variables may need to change for TRIQS compatibility.
        #Alter size of some of the arrays if spin orbit coupling is enabled.
        #For SO in Elk, the eigenvalues and eigenvector band indices are in asscending order w.r.t energy
        if(SO==1):
          [shells,corr_shells,dim_reps,n_orbitals,proj_mat,T,mat]=self.update_so_quatities(n_shells,shells,n_corr_shells,corr_shells,n_inequiv_shells,dim_reps,n_k,n_symm,n_orbitals,proj_mat,T,su2,mat)
          #reduce n_spin_blocs
          n_spin_blocs = SP + 1 - SO

        #put the energy eigenvalues arrays in TRIQS format
        hopping = self.sort_dft_eigvalues(n_spin_blocs,SO,n_k,n_orbitals,band_window,en,energy_unit)

        #Elk does not use global to local matrix rotation (Rotloc) as is done in Wien2k. However, the projectors
        #require a symmetry matrix to rotate from jatom to iatom. Below finds the non inversion
        #symmetric matrices which were used in calculating the projectors
        use_rotations = 1
        rot_mat = [numpy.identity(corr_shells[icrsh]['dim'], complex) for icrsh in range(n_corr_shells)]
        for icrsh in range(n_corr_shells):
          #return inequivalent index  
          incrsh = corr_to_inequiv[icrsh]
          #return first inequivalent corr_shell index
          jcrsh = inequiv_to_corr[incrsh]
          #want to rotate atom to first inequivalent atom in list
          iatom = corr_shells[jcrsh]['atom']
          for isym in range(n_symm):
            jatom=perm[isym][corr_shells[icrsh]['atom']-1]
            #determinant determines if crystal symmetry matrix has inversion symmetry (=-1)
            det = numpy.linalg.det(symmat[isym][:,:])
            if((jatom==iatom)&(det>0.0)):
              #local rotation which rotates equivalent atom into its local coordinate system
              #(inverse of the symmetry operator applied to the projectors in Elk)
              rot_mat[icrsh][:,:]=mat[isym][icrsh][:,:].conjugate().transpose()
              #used first desired symmetry in crystal symmetry list
              break
        # Elk does not currently use time inversion symmetry
        rot_mat_time_inv = [0 for i in range(n_corr_shells)]

        #Check that the charge of all the bands and below the correlated window have been calculated correctly
        [density_required, charge_below] = self.check_dens(n_k,nstsv,occ,bz_weights,n_spin_blocs,band_window,SO)
        #calculate the required charge (density_required) to remain charge neutral
        mpi.report("The total charge of the system = %f" %density_required)
        mpi.report("The charge below the correlated window = %f" %charge_below)
        mpi.report("The charge within the correlated window = %f" %(density_required - charge_below))

        #Elk interface does not calculate theta projectors, hence orbits are the same as Wannier projectors
        orbits=[]
        #remove the spatom index to avoid errors in the symmetry routines
        for ish in range(n_corr_shells):
           #remove "spatom"
           del corr_shells[ish]['spatom']
           orbits.append(corr_shells[ish].copy())
        for ish in range(n_shells):
           #remove "spatom"
           del shells[ish]['spatom']
        n_orbits=len(orbits)
        #Note that the T numpy array is defined for all shells.

        #new variable: dft_code - this determines which DFT code the inputs come from.
        #used for certain routines within dft_tools if treating the inputs differently is required.
        dft_code = 'elk'

        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file, 'a')
        if not (self.dft_subgrp in ar):
            ar.create_group(self.dft_subgrp)
        # The subgroup containing the data. If it does not exist, it is
        # created. If it exists, the data is overwritten!
        things_to_save = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                          'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                          'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr', 'dft_code']
        for it in things_to_save:
            ar[self.dft_subgrp][it] = locals()[it]
        del ar
        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file, 'a')
        symm_subgrp=self.symmcorr_subgrp
        #Elk does not use time inversion symmetry
        time_inv = [0 for j in range(n_symm)]
        mat_tinv = [numpy.identity(orbits[orb]['dim'], complex)
                        for orb in range(n_orbits)]
        #Save all the symmetry data
        if not (symm_subgrp in ar):
            ar.create_group(symm_subgrp)
        things_to_save_sym = ['n_symm', 'n_atoms', 'perm',
                          'orbits', 'SO', 'SP', 'time_inv', 'mat', 'mat_tinv']
        for it in things_to_save_sym:
            ar[symm_subgrp][it] = locals()[it]
        del ar
        #Save misc info
        things_to_save_misc = ['band_window','vkl','nstsv']
        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file, 'a')
        if not (self.misc_subgrp in ar):
            ar.create_group(self.misc_subgrp)
        for it in things_to_save_misc:
            ar[self.misc_subgrp][it] = locals()[it]
        del ar
        mpi.report('Converted the Elk ground state data')


    def convert_bands_input(self):
        """
        Reads the appropriate files and stores the data for the bands_subgrp in the hdf5 archive.

        """
        # Read and write only on the master node
        if not (mpi.is_master_node()):
            return
        filext='_WANBAND.OUT'
        dft_file='PROJ'+filext
        mpi.report("Reading %s" % dft_file)
        #Energy conversion - Elk uses Hartrees
        energy_unit = 27.2113850560                        # Elk uses hartrees
        shells=[]
        #read information about projectors calculated in the Elk calculation
        [gen_info,n_corr_shells,n_inequiv_shells,corr_to_inequiv,inequiv_to_corr,corr_shells,n_reps,dim_reps,ind,basis,T] = read_Elk.read_proj(self,dft_file)
        #get info for HDF5 file from gen_info
        n_k=gen_info['n_k']
        SP=gen_info['spinpol']-1
        #Elk uses spinor wavefunctions. Therefore these two spinor wavefunctions have spin-orbit coupling incorporated in them. Here we read in the spinors
        n_spin_blocs = SP + 1
        SO=gen_info['SO']
        #Elk only calculates Wannier projectors (no theta projectors generated):
        n_shells=n_corr_shells
        for ish in range(n_shells):
           shells.append(corr_shells[ish].copy())
           #remove last 2 entries from corr_shlls
           del shells[ish]['SO']
           del shells[ish]['irep']
           shells[ish]['dim'] = 2*shells[ish]['l']+1

        #read in the band eigenvalues
        mpi.report("Reading BAND.OUT")
        en=numpy.loadtxt('BAND.OUT')
        nstsv=int(len(en[:,1])/n_k)
        #convert the en array into a workable format
        entmp = numpy.zeros([n_k,nstsv], complex)
        enj=0
        for ist in range(nstsv):
          for ik in range(n_k):
            entmp[ik,ist]=en[enj,1]
            enj+=1
        del en
        #read projectors
        proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max([crsh['dim'] for crsh in corr_shells]), nstsv], complex)
        mpi.report("Reading projector(s)")
        for ish in range(n_corr_shells):
          [n_orbitals,band_window,rep,proj_mat]=read_Elk.read_projector(self,corr_shells,n_spin_blocs,ish,proj_mat,ind,T,basis,filext)

        #alter arrays for spin-orbit coupling
        if(SO==1):
          mat=[]
          su2=[]
          n_symm=1
          [shells,corr_shells,dim_reps,n_orbitals,proj_mat,T,mat]=self.update_so_quatities(n_shells,shells,n_corr_shells,corr_shells,n_inequiv_shells,dim_reps,n_k,n_symm,n_orbitals,proj_mat,T,su2,mat,sym=False)
          #reduce n_spin_blocs
          n_spin_blocs = SP + 1 - SO

        #put the energy eigenvalues arrays in TRIQS format
        hopping = self.sort_dft_eigvalues(n_spin_blocs,SO,n_k,n_orbitals,band_window,entmp,energy_unit)

        # No partial projectors generated here, so set to 0:
        n_parproj = numpy.array([0])
        proj_mat_all = numpy.array([0])

        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file, 'a')
        if not (self.bands_subgrp in ar):
            ar.create_group(self.bands_subgrp)
        # The subgroup containing the data. If it does not exist, it is
        # created. If it exists, the data is overwritten!
        things_to_save = ['n_k', 'n_orbitals', 'proj_mat',
                          'hopping', 'n_parproj', 'proj_mat_all']
        for it in things_to_save:
            ar[self.bands_subgrp][it] = locals()[it]
        del ar
        mpi.report('Converted the band data')

    def convert_contours_input(self,kgrid=None,ngrid=None):
        r"""
        Reads the appropriate files and stores the data for the cont_subgrp in the hdf5 archive.

        Parameters
        ----------
        kgrid : size (4,3) double numpy array, optional
                Numpy array defining the reciprocal lattice vertices used in the Elk Fermi
                surface calculation. Each row has the following meaning:
                grid3d[0,:] - origin lattice vertex
                grid3d[1,:] - b1 lattice vertex
                grid3d[2,:] - b2 lattice vertex
                grid3d[3,:] - b3 lattice vertex
        ngrid : size (3) integer numpy array, optional
                Numpy array for the number of points along each (b1,b2,b3) lattice vertices

        Note that these inputs relate to the plot3d input of Elk.
        """

        if not (mpi.is_master_node()):
            return
        filext='_FS.OUT'
        dft_file='PROJ'+filext
        mpi.report("Reading %s" % dft_file)
        #read the symmetries and k-points first
        #read kpoints calculated in the Elk FS calculation
        mpi.report("Reading KPOINT_FS.OUT")
        [bz_weights,vkl]=read_Elk.read_kpoints(self,filext=filext)
        n_k=vkl[:,0].size
        #Need lattice symmetries to unfold the irreducible BZ
        #Read symmetry files
        mpi.report("Reading SYMCRYS.OUT")
        [n_symm,spinmat,symlat,tr] = read_Elk.readsym(self)
        #generate full vectors for Fermi surface plotting along with index mapping
        #to irreducible vector set.
        if (ngrid is not None and kgrid is not None):
          mpi.report('Using User defined k-mesh')
        #check variables are in correct format
          if ngrid.size != 3:
            assert 0, "The input numpy ngrid is not the required size of 3!"
          elif ngrid.dtype != int:
            assert 0, "The input numpy ngrid is not an array of integers."
          elif kgrid.shape != (4,3):
            assert 0, "The input numpy kgrid is not the required size of (4x3)!"
        #generate full set of k-points with mapping to reduced set    
          [BZ_vkl, BZ_iknr, BZ_n_k] = Elk_tools.plotpt3d(self,n_k,vkl,n_symm,symlat,kgrid,ngrid)
        elif (ngrid is None and kgrid is None):
          mpi.report('No grid dimension input for Fermi surface.')
          mpi.report('Calculating k-points by folding out irreducible vectors instead if using symmetries.')
          mpi.report('Warning! This may not equate to the same set of vectors used to generate the Fermi surface data.')
          [BZ_vkl, BZ_iknr, BZ_n_k] = Elk_tools.bzfoldout(self,n_k,vkl,n_symm,symlat)
        else:
          assert 0, "Either input both ngrid and kgrid numpy arrays or neither."
        #return all threads apart from master
        if not (mpi.is_master_node()):
            return

        # Read and write the following only on the master thread
        #Energy conversion - Elk uses Hartrees
        energy_unit = 27.2113850560                        # Elk uses hartrees
        shells=[]
        #read information about projectors calculated in the Elk calculation
        [gen_info,n_corr_shells,n_inequiv_shells,corr_to_inequiv,inequiv_to_corr,corr_shells,n_reps,dim_reps,ind,basis,T] = read_Elk.read_proj(self,dft_file)
        #get info for HDF5 file from gen_info
        n_k=gen_info['n_k']
        SP=gen_info['spinpol']-1
        #Elk uses spinor wavefunctions. Therefore these two spinor wavefunctions have spin-orbit coupling incorporated in them. Here we read in the spinors
        n_spin_blocs = SP + 1
        SO=gen_info['SO']
        #Elk only calculates Wannier projectors (no theta projectors generated):
        n_shells=n_corr_shells
        for ish in range(n_shells):
           shells.append(corr_shells[ish].copy())
           #remove last 2 entries from corr_shlls
           del shells[ish]['SO']
           del shells[ish]['irep']
           shells[ish]['dim'] = 2*shells[ish]['l']+1

        #read in the eigenvalues used for the FS calculation
        mpi.report("Reading EIGVAL_FS.OUT and EFERMI.OUT")
        [en,occ,nstsv]=read_Elk.read_eig(self,filext=filext)

        #read projectors
        proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max([crsh['dim'] for crsh in corr_shells]), nstsv], complex)
        mpi.report("Reading projector(s)")
        for ish in range(n_corr_shells):
          [n_orbitals,band_window,rep,proj_mat]=read_Elk.read_projector(self,corr_shells,n_spin_blocs,ish,proj_mat,ind,T,basis,filext)

        mpi.report("Reading LATTICE.OUT")
        [amat,amatinv,bmat,bmatinv,cell_vol] = read_Elk.readlat(self)
        #Put eigenvalues into array of eigenvalues for the correlated window
        #alter arrays for spin-orbit coupling
        if(SO==1):
          mat=[]
          su2=[]
          [shells,corr_shells,dim_reps,n_orbitals,proj_mat,T,mat]=self.update_so_quatities(n_shells,shells,n_corr_shells,corr_shells,n_inequiv_shells,dim_reps,n_k,n_symm,n_orbitals,proj_mat,T,su2,mat,sym=False)
        #reduce n_spin_blocs
          n_spin_blocs = SP + 1 - SO

        #put the energy eigenvalues arrays in TRIQS format
        hopping = self.sort_dft_eigvalues(n_spin_blocs,SO,n_k,n_orbitals,band_window,en,energy_unit)

        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file, 'a')
        if not (self.cont_subgrp in ar):
            ar.create_group(self.cont_subgrp)
        # The subgroup containing the data. If it does not exist, it is
        # created. If it exists, the data is overwritten!
        things_to_save = ['n_k','n_orbitals', 'proj_mat','bmat',
                           'BZ_n_k','BZ_iknr','BZ_vkl','hopping']
        for it in things_to_save:
            ar[self.cont_subgrp][it] = locals()[it]
        del ar
        mpi.report('Converted the Contours data')

# commented out for now - unsure using this produces DFT+DMFT PDOS.
# The data from BC.OUT are the band-resolved diagonal muffin-tin DFT density matrix elements used in Elk to calculate PDOS 
# (the PDOS is calculated from the Trace over the bands indices). Although this is equivalent to using using projectors in DFT and is likely valid for DFT+DMFT, 
# the equivalence needs to be thoroughly checked for DFT+DMFT, but would require theta (or similar) projectors from Elk to do so. 
# code left here just in case.
        
#    def dft_band_characters(self):
#        """
#        Reads in the band-resolved muffin-tin density matrix (band characters) generated in Elk 
#        to be used for PDOS plots.
#        """

#        #determine file extension
#        fileext='.OUT'
#        #read number of k-points and eigenstates
#        things_to_read = ['n_k','n_orbitals']
#        ar = HDFArchive(self.hdf_file, 'r')
#        for it in things_to_read:
#            setattr(self, it, ar[self.dft_subgrp][it])
#        del ar

#        if not (mpi.is_master_node()):
#            return
#        mpi.report("Reading BC%s"%(fileext))

#        # get needed data from hdf file
#        # from general info
#        ar = HDFArchive(self.hdf_file, 'a')
#        things_to_read = ['SP', 'SO']
#        for it in things_to_read:
#           if not hasattr(self, it):
#              setattr(self, it, ar[self.dft_subgrp][it])
#        #from misc info
#        things_to_read = ['nstsv','band_window']
#        for it in things_to_read:
#           if not hasattr(self, it):
#              setattr(self, it, ar[self.misc_subgrp][it])
#        #from sym info
#        things_to_read = ['n_atoms','perm']
#        symm_subgrp=self.symmcorr_subgrp
#        for it in things_to_read:
#           if not hasattr(self, it):
#              setattr(self, it, ar[symm_subgrp][it])
#        del ar

#        #read in band characters
#        [bc,maxlm] = read_Elk.read_bc(self,fileext)
        #note that bc is the band resolved inner product of the Elk muffin-tin wave functions in 
        #a diagonal lm basis (by default). These are used in Elk to calculate the DOS in a diagonal 
        #irreducible lm basis. This band resolved density matrix in the lm basis will be used
        #to project to spectral funtion to get the muffin-tin contributions. There will be an 
        #interstitial contribution (as Elk uses an APW+lo basis) which is the difference between 
        #the total and summed muffin-tin contributions. Also note that the bc array should be 
        #symmetrised within Elk.
        #general variables
#        lmax = int(numpy.sqrt(maxlm)-1)
#        n_spin_blocs = self.SP + 1 - self.SO
#        so = self.SO + 1
#        #get the sort entry which is just the species index for Elk
#        [ns, na, atpos]=read_Elk.read_geometry(self)
#        isrt=0
#        sort=numpy.zeros([self.n_atoms],int)
#        #arrange sort(species) order
#        for i in range(ns):
#          for ia in range(na[i]):
#            sort[isrt]=i
#            isrt+=1
#        #updating n_shells to include all the atoms and l used in the Elk calculation.
#        n_shells = self.n_atoms * (lmax+1)
#        shells = []
#        shell_entries = ['atom', 'sort', 'l', 'dim']
#        for iat in range(self.n_atoms):
#          for l in range(lmax+1):
#            #sort is not known from Elk outputs  
#            tmp = [iat+1, sort[iat]+1, l, so*(2*l+1)]
#            shells.append({name: int(val) for name, val in zip(shell_entries, tmp)})
#        del tmp, ns, na, atpos, isrt, shell_entries
#        #overwrite n_shells and shells
#        things_to_save = ['n_shells', 'shells']
#        ar = HDFArchive(self.hdf_file, 'a')
#        for it in things_to_save:
#          ar[self.dft_subgrp][it] = locals()[it]


#        # Initialise P, here a double list of matrices:
#        band_dens_muffin = numpy.zeros([self.n_k, n_spin_blocs, n_shells, so*(2*lmax+1), numpy.max(self.n_orbitals)], float)
#        for ik in range(self.n_k):
#          for isp in range(n_spin_blocs):
#            ish=0
#            for iat in range(self.n_atoms):
#              for l in range(lmax+1):
#                #variables for putting subset of bc in proj_mat_all  
#                lm_min=l**2
#                lm_max=(l+1)**2
#                nst=self.n_orbitals[ik,isp]
#                ibot=self.band_window[isp][ik, 0]-1
#                itop=ibot+nst
#                dim=l*2+1
                #check use of abs (negative values should be close to 0)
#                band_dens_muffin[ik,isp,ish,0:dim,0:nst] = \
#                      bc[lm_min:lm_max,isp,iat,ibot:itop,ik]
#                if(self.SO==1):
#                  band_dens_muffin[ik,isp,ish,dim:2*dim,0:nst] = \
#                        bc[lm_min:lm_max,1,iat,ibot:itop,ik]
#                ish+=1

#        things_to_save = ['band_dens_muffin']
#        # Save it all to the HDF:
#        with HDFArchive(self.hdf_file, 'a') as ar:
#            if not (self.bc_subgrp in ar):
#                ar.create_group(self.bc_subgrp)
#            # The subgroup containing the data. If it does not exist, it is
#            # created. If it exists, the data is overwritten!
#            things_to_save = ['band_dens_muffin']
#            for it in things_to_save:
#                ar[self.bc_subgrp][it] = locals()[it]


    def convert_transport_input(self):
        """
        Reads the necessary information for transport calculations on:

        - the optical band window and the velocity matrix elements from :file:`case.pmat`

        and stores the data in the hdf5 archive.

        """

        if not (mpi.is_master_node()):
            return

        # get needed data from hdf file
        with HDFArchive(self.hdf_file, 'r') as ar:
            if not (self.dft_subgrp in ar):
                raise IOError("convert_transport_input: No %s subgroup in hdf file found! Call convert_dft_input first." % self.dft_subgrp)
            things_to_read = ['SP', 'SO','n_k','n_orbitals']
            for it in things_to_read:
               if not hasattr(self, it):
                  setattr(self, it, ar[self.dft_subgrp][it])
            #from misc info
            things_to_read = ['band_window','vkl','nstsv']
            for it in things_to_read:
               if not hasattr(self, it):
                  setattr(self, it, ar[self.misc_subgrp][it])

        #unlike in WIEN2k, Elk writes the velocities (momentum) matrix elements for all bands.
        #Therefore, we can use the indices in the n_orbitals array to extract the desired elements.
        #However, the PMAT.OUT file is in Fortran-binary, so the file is read in by python wrappers
        #around the reading fortran code.

        # Read relevant data from PMAT.OUT binary file
        ###########################################
        # band_window_optics: same as Elk converter's band_window, but rearranged to be compatible
        #                     for the transport calculations.
        # velocities_k: velocity (momentum) matrix elements between all bands in band_window_optics
        #               and each k-point.

        #load fortran wrapper module
        import triqs_dft_tools.converters.elktools.elkwrappers.getpmatelk as et
        #elk velocities for all bands
        pmat=numpy.zeros([self.nstsv,self.nstsv,3],dtype=complex)

        n_spin_blocks = self.SP + 1 - self.SO
        #TRIQS' velocities array used in its transport routines
        velocities_k = [[] for isp in range(n_spin_blocks)]
        #TRIQS' band_window array used in its transport routines
        band_window_optics = []


        mpi.report("Reading PMAT.OUT")
        #read velocities for each k-point
        for ik in range(self.n_k):
          #need to use a fortran array for wrapper
          f_vkl = numpy.asfortranarray(self.vkl[ik,:])
          #read the ik velocity using the wrapper
          pmat[:,:,:]=et.getpmatelk(ik+1,self.nstsv,f_vkl)
          #loop over spin
          for isp in range(n_spin_blocks):
            #no. correlated bands at ik
            nu1=self.band_window[isp][ik,0]-1
            nu2=self.band_window[isp][ik,1]-1
            n_bands=nu2-nu1+1
            #put into velocity array (code similar to that in wien.py.
            if n_bands <= 0:
                velocity_xyz = numpy.zeros((1, 1, 3), dtype=complex)
            else:
                velocity_xyz = numpy.zeros(
                    (n_bands, n_bands, 3), dtype=complex)
            #CHECK these lines
            velocity_xyz[:,:,:]=pmat[nu1:nu2+1,nu1:nu2+1,:]
            velocities_k[isp].append(velocity_xyz)

        #rearrange Elk's band_window array into band_window_optics array format
        for isp in range(n_spin_blocks):
          band_window_optics_isp = []
          for ik in range(self.n_k):
            nu1=self.band_window[isp][ik,0]
            nu2=self.band_window[isp][ik,1]
            band_window_optics_isp.append((nu1, nu2))
            n_bands=nu2-nu1+1
          band_window_optics.append(numpy.array(band_window_optics_isp))

        #read in the cell volume from LATTICE.OUT
        mpi.report("Reading LATTICE.OUT")
        [amat,amatinv,bmat,bmatinv,cell_vol] = read_Elk.readlat(self)

        #read in the crystal symmetries
        mpi.report("Reading SYMCRYS.OUT")
        [n_symmetries,spinmat,rot_symmetries,tr] = read_Elk.readsym(self)


        # Put data to HDF5 file
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.transp_subgrp in ar):
                ar.create_group(self.transp_subgrp)
            # The subgroup containing the data. If it does not exist, it is
            # created. If it exists, the data is overwritten!!!
            things_to_save = ['band_window_optics', 'velocities_k']
            for it in things_to_save:
                ar[self.transp_subgrp][it] = locals()[it]
            things_to_save_misc = ['n_symmetries', 'rot_symmetries','cell_vol']
            for it in things_to_save_misc:
              ar[self.misc_subgrp][it] = locals()[it]
        mpi.report("Reading complete!")
