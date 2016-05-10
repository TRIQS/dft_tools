 
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
import pytriqs.utility.dichotomy as dichotomy
from pytriqs.gf.local import *
import pytriqs.utility.mpi as mpi
from pytriqs.archive import *
from symmetry import *
from sets import Set
from itertools import product

class SumkDFT:
    """This class provides a general SumK method for combining ab-initio code and pytriqs."""


    def __init__(self, hdf_file, h_field = 0.0, use_dft_blocks = False, 
	         dft_data = 'dft_input', symmcorr_data = 'dft_symmcorr_input', parproj_data = 'dft_parproj_input', 
                 symmpar_data = 'dft_symmpar_input', bands_data = 'dft_bands_input', transp_data = 'dft_transp_input',
                 misc_data = 'dft_misc_input'):
        r"""
        Initialises the class from data previously stored into an hdf5 archive.

        Parameters
        ----------
        hdf_file : string
                   Name of hdf5 containing the data.
        h_field : scalar, optional
                  The value of magnetic field to add to the DFT Hamiltonian. 
                  The contribution -h_field*sigma is added to diagonal elements of the Hamiltonian.
                  It cannot be used with the spin-orbit coupling on; namely h_field is set to 0 if self.SO=True.
        use_dft_blocks : boolean, optional
                         If True, the local Green's function matrix for each spin is divided into smaller blocks 
                          with the block structure determined from the DFT density matrix of the corresponding correlated shell.
        dft_data : string, optional
                   Name of hdf5 subgroup in which DFT data for projector and lattice Green's function construction are stored.
        symmcorr_data : string, optional
                        Name of hdf5 subgroup in which DFT data on symmetries of correlated shells 
                        (symmetry operations, permutaion matrices etc.) are stored.
        parproj_data : string, optional
                       Name of hdf5 subgroup in which DFT data on non-normalized projectors for non-correlated
                       states (used in the partial density of states calculations) are stored.
        symmpar_data : string, optional
                       Name of hdf5 subgroup in which DFT data on symmetries of the non-normalized projectors
                       are stored.
        bands_data : string, optional
                     Name of hdf5 subgroup in which DFT data necessary for band-structure/k-resolved spectral
                     function calculations (projectors, DFT Hamiltonian for a chosen path in the Brillouin zone etc.)
                     are stored.
        transp_data : string, optional
                      Name of hdf5 subgroup in which DFT data necessary for transport calculations are stored.
        misc_data : string, optional
                    Name of hdf5 subgroup in which miscellaneous DFT data are stored.
        """

        if not type(hdf_file) == StringType:
            mpi.report("Give a string for the hdf5 filename to read the input!")
        else:
            self.hdf_file = hdf_file
            self.dft_data = dft_data
            self.symmcorr_data = symmcorr_data
            self.parproj_data = parproj_data
            self.symmpar_data = symmpar_data
            self.bands_data = bands_data
            self.transp_data = transp_data
            self.misc_data = misc_data
            self.h_field = h_field

            # Read input from HDF:
            things_to_read = ['energy_unit','n_k','k_dep_projection','SP','SO','charge_below','density_required',
                              'symm_op','n_shells','shells','n_corr_shells','corr_shells','use_rotations','rot_mat',
                              'rot_mat_time_inv','n_reps','dim_reps','T','n_orbitals','proj_mat','bz_weights','hopping',
                              'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
            self.subgroup_present, self.value_read = self.read_input_from_hdf(subgrp = self.dft_data, things_to_read = things_to_read)
            if self.symm_op: self.symmcorr = Symmetry(hdf_file,subgroup=self.symmcorr_data)

            if self.SO and (abs(self.h_field) > 0.000001):
                self.h_field = 0.0
                mpi.report("For SO, the external magnetic field is not implemented, setting it to 0!")

            self.spin_block_names = [ ['up','down'], ['ud'] ]
            self.n_spin_blocks = [2,1]
            # Convert spin_block_names to indices -- if spin polarized, differentiate up and down blocks 
            self.spin_names_to_ind = [{}, {}]
            for iso in range(2): # SO = 0 or 1
                for isp in range(self.n_spin_blocks[iso]):
                    self.spin_names_to_ind[iso][self.spin_block_names[iso][isp]] = isp * self.SP

            # GF structure used for the local things in the k sums
            # Most general form allowing for all hybridisation, i.e. largest blocks possible
            self.gf_struct_sumk = [ [ (sp, range( self.corr_shells[icrsh]['dim'])) for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']] ]
                                   for icrsh in range(self.n_corr_shells) ]
            # First set a standard gf_struct solver:
            self.gf_struct_solver = [ dict([ (sp, range(self.corr_shells[self.inequiv_to_corr[ish]]['dim']) )
                                                    for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']] ])
                                      for ish in range(self.n_inequiv_shells) ]
            # Set standard (identity) maps from gf_struct_sumk <-> gf_struct_solver
            self.sumk_to_solver = [ {} for ish in range(self.n_inequiv_shells) ]
            self.solver_to_sumk = [ {} for ish in range(self.n_inequiv_shells) ]
            self.solver_to_sumk_block = [ {} for ish in range(self.n_inequiv_shells) ]
            for ish in range(self.n_inequiv_shells):
                for block,inner_list in self.gf_struct_sumk[self.inequiv_to_corr[ish]]:
                    self.solver_to_sumk_block[ish][block] = block
                    for inner in inner_list:
                        self.sumk_to_solver[ish][(block,inner)] = (block,inner)
                        self.solver_to_sumk[ish][(block,inner)] = (block,inner)
            self.deg_shells = [ [] for ish in range(self.n_inequiv_shells) ] # assume no shells are degenerate

            self.chemical_potential = 0.0 # initialise mu
            self.init_dc() # initialise the double counting

            # Analyse the block structure and determine the smallest gf_struct blocks and maps, if desired
            if use_dft_blocks: self.analyse_block_structure()

################
# hdf5 FUNCTIONS
################

    def read_input_from_hdf(self, subgrp, things_to_read):
        r"""
        Reads data from the HDF file. Prints a warning if a requested dataset is not found.

        Parameters
        ----------
        subgrp : string
                 Name of hdf5 file subgroup from which the data are to be read.
        things_to_read : list of strings
                         List of datasets to be read from the hdf5 file.

        Returns
        -------
        subgroup_present : boolean
                           Is the subgrp is present in hdf5 file?
        value_read : boolean
                     Did the reading of requested datasets succeed?

        """

        value_read = True
        # initialise variables on all nodes to ensure mpi broadcast works at the end
        for it in things_to_read: setattr(self,it,0)
        subgroup_present = 0

        if mpi.is_master_node():
            ar = HDFArchive(self.hdf_file,'r')
            if subgrp in ar:
                subgroup_present = True
                # first read the necessary things:
                for it in things_to_read:
                    if it in ar[subgrp]:
                        setattr(self,it,ar[subgrp][it])
                    else:
                        mpi.report("Loading %s failed!"%it)
                        value_read = False
            else:
                if (len(things_to_read) != 0): mpi.report("Loading failed: No %s subgroup in hdf5!"%subgrp)
                subgroup_present = False
                value_read = False
            del ar
        # now do the broadcasting:
        for it in things_to_read: setattr(self,it,mpi.bcast(getattr(self,it)))
        subgroup_present = mpi.bcast(subgroup_present)
        value_read = mpi.bcast(value_read)

        return subgroup_present, value_read


    def save(self, things_to_save, subgrp='user_data'):
       
        r"""
        Saves data from a list into the HDF file. Prints a warning if a requested data is not found in SumkDFT object.

        Parameters
        ----------
        things_to_save : list of strings
                         List of datasets to be saved into the hdf5 file.
        subgrp : string, optional
                 Name of hdf5 file subgroup in which the data are to be stored.
        """

        if not (mpi.is_master_node()): return # do nothing on nodes
        ar = HDFArchive(self.hdf_file,'a')
        if not subgrp in ar: ar.create_group(subgrp)
        for it in things_to_save: 
            try:
                ar[subgrp][it] = getattr(self,it)
            except:
                mpi.report("%s not found, and so not saved."%it)
        del ar


    def load(self, things_to_load, subgrp='user_data'):
        r"""
        Loads user data from the HDF file. Raises an exeption if a requested dataset is not found.

        Parameters
        ----------
        things_to_read : list of strings
                         List of datasets to be read from the hdf5 file.
        subgrp : string, optional
                 Name of hdf5 file subgroup from which the data are to be read.

        Returns
        -------
        list_to_return : list
                         A list containing data read from hdf5.
        """

        if not (mpi.is_master_node()): return # do nothing on nodes
        ar = HDFArchive(self.hdf_file,'r')
        if not subgrp in ar: mpi.report("Loading %s failed!"%subgrp)
        list_to_return = []
        for it in things_to_load: 
            try:
                list_to_return.append(ar[subgrp][it])
            except:
                raise ValueError, "load: %s not found, and so not loaded."%it 
        del ar
        return list_to_return

################
# CORE FUNCTIONS
################

    def downfold(self,ik,ish,bname,gf_to_downfold,gf_inp,shells='corr',ir=None):
        r"""
        Downfolds a block of the Green's function for a given shell and k-point using the corresponding projector matrices.

        Parameters
        ----------
        ik : integer
             k-point index for which the downfolding is to be done.
        ish : integer
              Shell index of GF to be downfolded.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        bname : string
                Block name of the target block of the lattice Green's function.
        gf_to_downfold : Gf 
                       Block of the Green's function that is to be downfolded.
        gf_inp : Gf 
                 FIXME 
        shells : string, optional

                 - if shells='corr': orthonormalized projectors for correlated shells are used for the downfolding.
                 - if shells='all': non-normalized projectors for all included shells are used for the downfolding.

        ir : integer, optional
             Index of equivalent site in the non-correlated shell 'ish', only used if shells='all'.
             
        Returns
        -------
        gf_downfolded : Gf
                      Downfolded block of the lattice Green's function.
        """
  
        gf_downfolded = gf_inp.copy()
        isp = self.spin_names_to_ind[self.SO][bname]       # get spin index for proj. matrices
        n_orb = self.n_orbitals[ik,isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            projmat = self.proj_mat[ik,isp,ish,0:dim,0:n_orb]
        elif shells == 'all':
            if ir is None: raise ValueError, "downfold: provide ir if treating all shells."
            dim = self.shells[ish]['dim']
            projmat = self.proj_mat_all[ik,isp,ish,ir,0:dim,0:n_orb]
  
        gf_downfolded.from_L_G_R(projmat,gf_to_downfold,projmat.conjugate().transpose())
  
        return gf_downfolded


    def upfold(self,ik,ish,bname,gf_to_upfold,gf_inp,shells='corr',ir=None):
        r"""
        Upfolds a block of the Green's function for a given shell and k-point using the corresponding projector matrices.

        Parameters
        ----------
        ik : integer
             k-point index for which the upfolding is to be done.
        ish : integer
              Shell index of GF to be upfolded.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        bname : string
                Block name of the target block of the lattice Green's function.
        gf_to_upfold : Gf 
                       Block of the Green's function that is to be upfolded.
        gf_inp : Gf 
                 FIXME 
        shells : string, optional

                 - if shells='corr': orthonormalized projectors for correlated shells are used for the upfolding.
                 - if shells='all': non-normalized projectors for all included shells are used for the upfolding.

        ir : integer, optional
             Index of equivalent site in the non-correlated shell 'ish', only used if shells='all'.
             
        Returns
        -------
        gf_upfolded : Gf
                      Upfolded block of the lattice Green's function.
        """
  
        gf_upfolded = gf_inp.copy()
        isp = self.spin_names_to_ind[self.SO][bname]       # get spin index for proj. matrices
        n_orb = self.n_orbitals[ik,isp]
        if shells == 'corr':
            dim = self.corr_shells[ish]['dim']
            projmat = self.proj_mat[ik,isp,ish,0:dim,0:n_orb]
        elif shells == 'all':
            if ir is None: raise ValueError, "upfold: provide ir if treating all shells."
            dim = self.shells[ish]['dim']
            projmat = self.proj_mat_all[ik,isp,ish,ir,0:dim,0:n_orb]
  
        gf_upfolded.from_L_G_R(projmat.conjugate().transpose(),gf_to_upfold,projmat)
  
        return gf_upfolded


    def rotloc(self,ish,gf_to_rotate,direction,shells='corr'):
        r"""
        Rotates a block of the local Green's function from the local frame to the global frame and vice versa.

        Parameters
        ----------
        ish : integer
              Shell index of GF to be rotated.

              - if shells='corr': ish labels all correlated shells (equivalent or not)
              - if shells='all': ish labels only representative (inequivalent) non-correlated shells

        gf_to_rotate : Gf 
                       Block of the Green's function that is to be rotated.
        direction : string
                    The direction of rotation can be either 

                    - 'toLocal' : global -> local transformation,
                    - 'toGlobal' : local -> global transformation.

        shells : string, optional

                 - if shells='corr': the rotation matrix for the correlated shell 'ish' is used,
                 - if shells='all': the rotation matrix for the generic (non-correlated) shell 'ish' is used.

        Returns
        -------
        gf_rotated : Gf
                     Rotated block of the local Green's function.
        """

        assert ((direction == 'toLocal') or (direction == 'toGlobal')),"rotloc: Give direction 'toLocal' or 'toGlobal'."
        gf_rotated = gf_to_rotate.copy()
        if shells == 'corr':
            rot_mat_time_inv = self.rot_mat_time_inv
            rot_mat = self.rot_mat
        elif shells == 'all':
            rot_mat_time_inv = self.rot_mat_all_time_inv
            rot_mat = self.rot_mat_all

        if direction == 'toGlobal':

            if (rot_mat_time_inv[ish] == 1) and self.SO:
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(rot_mat[ish].conjugate(),gf_rotated,rot_mat[ish].transpose())
            else:
                gf_rotated.from_L_G_R(rot_mat[ish],gf_rotated,rot_mat[ish].conjugate().transpose())

        elif direction == 'toLocal':

            if (rot_mat_time_inv[ish] == 1) and self.SO:
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(rot_mat[ish].transpose(),gf_rotated,rot_mat[ish].conjugate())
            else:
                gf_rotated.from_L_G_R(rot_mat[ish].conjugate().transpose(),gf_rotated,rot_mat[ish])

        return gf_rotated


    def lattice_gf(self, ik, mu=None, iw_or_w="iw", beta=40, broadening=None, mesh=None, with_Sigma=True, with_dc=True):
        r"""
        Calculates the lattice Green function for a given k-point from the DFT Hamiltonian and the self energy. 

        Parameters
        ----------
        ik : integer
             k-point index.
        mu : real, optional
             Chemical potential for which the Green's function is to be calculated.
             If not provided, self.chemical_potential is used for mu.
        iw_or_w : string, optional

                  - `iw_or_w` = 'iw' for a imaginary-frequency self-energy
                  - `iw_or_w` = 'w' for a real-frequency self-energy

        beta : real, optional
               Inverse temperature.
        broadening : real, optional
                     Imaginary shift for the axis along which the real-axis GF is calculated.
                     If not provided, broadening will be set to double of the distance between mesh points in 'mesh'.
        mesh : list, optional
               Data defining mesh on which the real-axis GF will be calculated, given in the form
               (om_min,om_max,n_points), where om_min is the minimum omega, om_max is the maximum omega and n_points is the number of points.
        with_Sigma : boolean, optional
                     If True the GF will be calculated with the self-energy stored in self.Sigmaimp_(w/iw), for real/Matsubara GF, respectively. 
                     In this case the mesh is taken from the self.Sigma_imp object.
                     If with_Sigma=True but self.Sigmaimp_(w/iw) is not present, with_Sigma is reset to False.
        with_dc : boolean, optional
                  if True and with_Sigma=True, the dc correction is substracted from the self-energy before it is included into GF.
                     
        Returns
        -------
        G_latt : BlockGf
                 Lattice Green's function.
       
        """
        if mu is None: mu = self.chemical_potential
        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        if (iw_or_w != "iw") and (iw_or_w != "w"): raise ValueError, "lattice_gf: Implemented only for Re/Im frequency functions."
        if not hasattr(self,"Sigma_imp_"+iw_or_w): with_Sigma = False
        if broadening is None:
            if mesh is None:
                broadening = 0.01
            else: # broadening = 2 * \Delta omega, where \Delta omega is the spacing of omega points
                broadening = 2.0 * ( (mesh[1]-mesh[0])/(mesh[2]-1) )

        # Are we including Sigma?
        if with_Sigma:
            Sigma_imp = getattr(self,"Sigma_imp_"+iw_or_w)
            sigma_minus_dc = [s.copy() for s in Sigma_imp]
            if with_dc: sigma_minus_dc = self.add_dc(iw_or_w)
            if iw_or_w == "iw":
                beta = Sigma_imp[0].mesh.beta   # override beta if Sigma_iw is present
                mesh = Sigma_imp[0].mesh
            elif iw_or_w == "w":
                mesh = Sigma_imp[0].mesh
        else:
            if iw_or_w == "iw":
                if beta is None: raise ValueError, "lattice_gf: Give the beta for the lattice GfReFreq."
                mesh = MeshImFreq(beta=beta, S='Fermion', n_max=1025) # Default number of Matsubara frequencies
            elif iw_or_w == "w":
                if mesh is None: raise ValueError, "lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq."
                mesh = MeshReFreq(mesh[0],mesh[1],mesh[2])

        # Check if G_latt is present
        set_up_G_latt = False                       # Assume not
        if not hasattr(self,"G_latt_"+iw_or_w):
            set_up_G_latt = True                    # Need to create G_latt_(i)w
        else:                                       # Check that existing GF is consistent
            G_latt = getattr(self,"G_latt_"+iw_or_w)
            GFsize = [ gf.N1 for bname,gf in G_latt]
            unchangedsize = all( [ self.n_orbitals[ik,ntoi[spn[isp]]]==GFsize[isp] for isp in range(self.n_spin_blocks[self.SO]) ] )
            if not unchangedsize: set_up_G_latt = True
            if (iw_or_w == "iw") and (self.G_latt_iw.mesh.beta != beta): set_up_G_latt = True # additional check for ImFreq

        # Set up G_latt
        if set_up_G_latt:
            block_structure = [ range(self.n_orbitals[ik,ntoi[sp]]) for sp in spn ]
            gf_struct = [ (spn[isp], block_structure[isp]) for isp in range(self.n_spin_blocks[self.SO]) ]
            block_ind_list = [block for block,inner in gf_struct]
            if iw_or_w == "iw":
               glist = lambda : [ GfImFreq(indices=inner,mesh=mesh) for block,inner in gf_struct ]
            elif iw_or_w == "w":
               glist = lambda : [ GfReFreq(indices=inner,mesh=mesh) for block,inner in gf_struct ]
            G_latt = BlockGf(name_list = block_ind_list, block_list = glist(), make_copies = False)
            G_latt.zero()

        if iw_or_w == "iw":
            G_latt << iOmega_n
        elif iw_or_w == "w":
            G_latt << Omega + 1j*broadening

        idmat = [numpy.identity(self.n_orbitals[ik,ntoi[sp]],numpy.complex_) for sp in spn]
        M = copy.deepcopy(idmat)
        for ibl in range(self.n_spin_blocks[self.SO]):
            ind = ntoi[spn[ibl]]
            n_orb = self.n_orbitals[ik,ind]
            M[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb] - (idmat[ibl]*mu) - (idmat[ibl] * self.h_field * (1-2*ibl))
        G_latt -= M

        if with_Sigma:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in G_latt: gf -= self.upfold(ik,icrsh,bname,sigma_minus_dc[icrsh][bname],gf)

        G_latt.invert()
        setattr(self,"G_latt_"+iw_or_w,G_latt)

        return G_latt

    def set_Sigma(self,Sigma_imp):
        self.put_Sigma(Sigma_imp)

    def put_Sigma(self, Sigma_imp):
        r"""
        Inserts the impurity self-energies into the sumk_dft class.

        Parameters
        ----------
        Sigma_imp : list of BlockGf (Green's function) objects
                    List containing impurity self-energy for all inequivalent correlated shells.
                    Self-energies for equivalent shells are then automatically set by this function.
                    The self-energies can be of the real or imaginary-frequency type.
        """

        assert isinstance(Sigma_imp,list), "put_Sigma: Sigma_imp has to be a list of Sigmas for the correlated shells, even if it is of length 1!"
        assert len(Sigma_imp) == self.n_inequiv_shells, "put_Sigma: give exactly one Sigma for each inequivalent corr. shell!"

        # init self.Sigma_imp_(i)w:
        if all(type(gf) == GfImFreq for bname,gf in Sigma_imp[0]):
            # Imaginary frequency Sigma:
            self.Sigma_imp_iw = [ BlockGf( name_block_generator = [ (block,GfImFreq(indices = inner, mesh = Sigma_imp[0].mesh)) 
                                                                    for block,inner in self.gf_struct_sumk[icrsh] ], make_copies = False) 
                                  for icrsh in range(self.n_corr_shells) ]
            SK_Sigma_imp = self.Sigma_imp_iw
        elif all(type(gf) == GfReFreq for bname,gf in Sigma_imp[0]):
            # Real frequency Sigma:
            self.Sigma_imp_w = [ BlockGf( name_block_generator = [ (block,GfReFreq(indices = inner, mesh = Sigma_imp[0].mesh)) 
                                                                   for block,inner in self.gf_struct_sumk[icrsh] ], make_copies = False) 
                                 for icrsh in range(self.n_corr_shells) ]
            SK_Sigma_imp = self.Sigma_imp_w
        else:
            raise ValueError, "put_Sigma: This type of Sigma is not handled."

        # transform the CTQMC blocks to the full matrix:
        for icrsh in range(self.n_corr_shells):
            ish = self.corr_to_inequiv[icrsh]    # ish is the index of the inequivalent shell corresponding to icrsh
            for block,inner in self.gf_struct_solver[ish].iteritems():
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk,ind1_sumk = self.solver_to_sumk[ish][(block,ind1)]
                        block_sumk,ind2_sumk = self.solver_to_sumk[ish][(block,ind2)]
                        SK_Sigma_imp[icrsh][block_sumk][ind1_sumk,ind2_sumk] << Sigma_imp[ish][block][ind1,ind2]

        # rotation from local to global coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in SK_Sigma_imp[icrsh]: gf << self.rotloc(icrsh,gf,direction='toGlobal')

    def extract_G_loc(self, mu=None, with_Sigma=True, with_dc=True):
        r"""
        Extracts the local downfolded Green function by the Brillouin-zone integration of the lattice Green's function.

        Parameters
        ----------
        mu : real, optional
             Input chemical potential. If not provided the value of self.chemical_potential is used as mu.
        with_Sigma : boolean, optional
                     If True then the local GF is calculated with the self-energy self.Sigma_imp.
        with_dc : boolean, optional
                  If True then the double-counting correction is subtracted from the self-energy in calculating the GF.

        Returns
        -------
        G_loc_inequiv : list of BlockGf (Green's function) objects
                        List of the local Green's functions for all inequivalent correlated shells, 
                        rotated into the corresponding local frames.
       
        """

        if mu is None: mu = self.chemical_potential
        G_loc = [ self.Sigma_imp_iw[icrsh].copy() for icrsh in range(self.n_corr_shells) ]   # this list will be returned
        for icrsh in range(self.n_corr_shells): G_loc[icrsh].zero()                          # initialize to zero
        beta = G_loc[0].mesh.beta

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):

            G_latt_iw = self.lattice_gf(ik=ik, mu=mu, iw_or_w="iw", with_Sigma=with_Sigma, with_dc=with_dc, beta=beta)
            G_latt_iw *= self.bz_weights[ik]

            for icrsh in range(self.n_corr_shells):
                tmp = G_loc[icrsh].copy()                  # init temporary storage
                for bname,gf in tmp: tmp[bname] << self.downfold(ik,icrsh,bname,G_latt_iw[bname],gf)
                G_loc[icrsh] += tmp

        # Collect data from mpi
        for icrsh in range(self.n_corr_shells): 
            G_loc[icrsh] << mpi.all_reduce(mpi.world, G_loc[icrsh], lambda x,y : x+y)
        mpi.barrier()

        # G_loc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:
        if self.symm_op != 0: G_loc = self.symmcorr.symmetrize(G_loc)

        # G_loc is rotated to the local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in G_loc[icrsh]: G_loc[icrsh][bname] << self.rotloc(icrsh,gf,direction='toLocal')

        # transform to CTQMC blocks:
        G_loc_inequiv = [ BlockGf( name_block_generator = [ (block,GfImFreq(indices = inner, mesh = G_loc[0].mesh)) for block,inner in self.gf_struct_solver[ish].iteritems() ],
                        make_copies = False) for ish in range(self.n_inequiv_shells)  ]
        for ish in range(self.n_inequiv_shells):
            for block,inner in self.gf_struct_solver[ish].iteritems():
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk,ind1_sumk = self.solver_to_sumk[ish][(block,ind1)]
                        block_sumk,ind2_sumk = self.solver_to_sumk[ish][(block,ind2)]
                        G_loc_inequiv[ish][block][ind1,ind2] << G_loc[self.inequiv_to_corr[ish]][block_sumk][ind1_sumk,ind2_sumk]

        # return only the inequivalent shells:
        return G_loc_inequiv


    def analyse_block_structure(self, threshold = 0.00001, include_shells = None, dm = None):
        r"""
        Determines the block structure of local Green's functions by analysing the structure of 
        the corresponding density matrices. The resulting block structures for correlated shells
        are stored in self.gf_struct_solver list.

        Parameters
        ----------
        threshold : real, optional
                    If the difference between density matrix elements is below threshold,
                    they are considered to be equal.
        include_shells : list of integers, optional
                         List of correlated shells to be analysed.
                         If include_shells is not provided all correlated shells will be analysed.
        dm : list of dict, optional
             List of density matrices from which block stuctures are to be analysed.
             Each density matrix is a dict {block names: 2d numpy arrays}.
             If not provided, dm will be calculated from the DFT Hamiltonian by a simple-point BZ integration.
        """

        self.gf_struct_solver = [ {} for ish in range(self.n_inequiv_shells) ]
        self.sumk_to_solver = [ {} for ish in range(self.n_inequiv_shells) ]
        self.solver_to_sumk = [ {} for ish in range(self.n_inequiv_shells) ]
        self.solver_to_sumk_block = [ {} for ish in range(self.n_inequiv_shells) ]

        if dm is None: dm = self.density_matrix(method = 'using_point_integration')
        dens_mat = [ dm[self.inequiv_to_corr[ish]] for ish in range(self.n_inequiv_shells) ]

        if include_shells is None: include_shells = range(self.n_inequiv_shells)
        for ish in include_shells:


            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                n_orb = self.corr_shells[self.inequiv_to_corr[ish]]['dim']
                dmbool = (abs(dens_mat[ish][sp]) > threshold)  # gives an index list of entries larger that threshold

                # Determine off-diagonal entries in upper triangular part of density matrix
                offdiag = Set([])
                for i in range(n_orb):
                    for j in range(i+1,n_orb):
                        if dmbool[i,j]: offdiag.add((i,j))

                # Determine the number of non-hybridising blocks in the gf
                blocs = [ [i] for i in range(n_orb) ]
                while len(offdiag) != 0:
                    pair = offdiag.pop()
                    for b1,b2 in product(blocs,blocs):
                       if (pair[0] in b1) and (pair[1] in b2):
                           if blocs.index(b1) != blocs.index(b2):     # In separate blocks?
                               b1.extend(blocs.pop(blocs.index(b2)))  # Merge two blocks
                               break                                  # Move on to next pair in offdiag

                # Set the gf_struct for the solver accordingly
                num_blocs = len(blocs)
                for i in range(num_blocs):
                    blocs[i].sort()
                    self.gf_struct_solver[ish].update( [('%s_%s'%(sp,i),range(len(blocs[i])))] )

                # Construct sumk_to_solver taking (sumk_block, sumk_index) --> (solver_block, solver_inner)
                #       and solver_to_sumk taking (solver_block, solver_inner) --> (sumk_block, sumk_index)
                for i in range(num_blocs):
                    for j in range(len(blocs[i])):
                        block_sumk = sp
                        inner_sumk = blocs[i][j]
                        block_solv = '%s_%s'%(sp,i)
                        inner_solv = j
                        self.sumk_to_solver[ish][(block_sumk,inner_sumk)] = (block_solv,inner_solv)
                        self.solver_to_sumk[ish][(block_solv,inner_solv)] = (block_sumk,inner_sumk)
                        self.solver_to_sumk_block[ish][block_solv] = block_sumk

            # Now calculate degeneracies of orbitals
            dm = {}
            for block,inner in self.gf_struct_solver[ish].iteritems():
                # get dm for the blocks:
                dm[block] = numpy.zeros([len(inner),len(inner)],numpy.complex_)
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk,ind1_sumk = self.solver_to_sumk[ish][(block,ind1)]
                        block_sumk,ind2_sumk = self.solver_to_sumk[ish][(block,ind2)]
                        dm[block][ind1,ind2] = dens_mat[ish][block_sumk][ind1_sumk,ind2_sumk]

            for block1 in self.gf_struct_solver[ish].iterkeys():
                for block2 in self.gf_struct_solver[ish].iterkeys(): 
                    if dm[block1].shape == dm[block2].shape:
                        if ( (abs(dm[block1] - dm[block2]) < threshold).all() ) and (block1 != block2):
                            ind1 = -1
                            ind2 = -2
                            # check if it was already there:
                            for n,ind in enumerate(self.deg_shells[ish]):
                                if block1 in ind: ind1 = n
                                if block2 in ind: ind2 = n
                            if (ind1 < 0) and (ind2 >= 0):
                                self.deg_shells[ish][ind2].append(block1)
                            elif (ind1 >= 0) and (ind2 < 0):
                                self.deg_shells[ish][ind1].append(block2)
                            elif (ind1 < 0) and (ind2 < 0):
                                self.deg_shells[ish].append([block1,block2])


    def density_matrix(self, method = 'using_gf', beta = 40.0):
        """Calculate density matrices in one of two ways.

        Parameters
        ----------
        method : string, optional

                 - if 'using_gf': First get lattice gf (g_loc is not set up), then density matrix.
                                  It is useful for Hubbard I, and very quick.
                                  No assumption on the hopping structure is made (ie diagonal or not).
                 - if 'using_point_integration': Only works for diagonal hopping matrix (true in wien2k).

        beta : float, optional
               Inverse temperature.      

        Returns
        -------
        dens_mat : list of dicts
                   Density matrix for each spin in each correlated shell.
        """
        dens_mat = [ {} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']]:
                dens_mat[icrsh][sp] = numpy.zeros([self.corr_shells[icrsh]['dim'],self.corr_shells[icrsh]['dim']], numpy.complex_)

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):

            if method == "using_gf":

                G_latt_iw = self.lattice_gf(ik = ik, mu = self.chemical_potential, iw_or_w = "iw", beta = beta)
                G_latt_iw *= self.bz_weights[ik]
                dm = G_latt_iw.density()
                MMat = [dm[sp] for sp in self.spin_block_names[self.SO]]

            elif method == "using_point_integration":

                ntoi = self.spin_names_to_ind[self.SO]
                spn = self.spin_block_names[self.SO]
                unchangedsize = all( [self.n_orbitals[ik,ntoi[sp]] == self.n_orbitals[0,ntoi[sp]] for sp in spn] )
                if unchangedsize:
                    dim = self.n_orbitals[0,ntoi[sp]]
                else:
                    dim = self.n_orbitals[ik,ntoi[sp]]
                MMat = [numpy.zeros( [dim,dim], numpy.complex_) for sp in spn]

                for isp, sp in enumerate(spn):
                    ind = ntoi[sp]
                    for inu in range(self.n_orbitals[ik,ind]):
                        if (self.hopping[ik,ind,inu,inu] - self.h_field*(1-2*isp)) < 0.0: # only works for diagonal hopping matrix (true in wien2k)
                            MMat[isp][inu,inu] = 1.0
                        else:
                            MMat[isp][inu,inu] = 0.0

            else: raise ValueError, "density_matrix: the method '%s' is not supported."%method

            for icrsh in range(self.n_corr_shells):
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[self.corr_shells[icrsh]['SO']][sp]
                    dim = self.corr_shells[icrsh]['dim']
                    n_orb = self.n_orbitals[ik,ind]
                    projmat = self.proj_mat[ik,ind,icrsh,0:dim,0:n_orb]
                    if method == "using_gf":
                        dens_mat[icrsh][sp] += numpy.dot( numpy.dot(projmat,MMat[isp]),
                                                          projmat.transpose().conjugate() )
                    elif method == "using_point_integration":
                        dens_mat[icrsh][sp] += self.bz_weights[ik] * numpy.dot( numpy.dot(projmat,MMat[isp]) ,
                                                                                projmat.transpose().conjugate() )

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for sp in dens_mat[icrsh]:
                dens_mat[icrsh][sp] = mpi.all_reduce(mpi.world, dens_mat[icrsh][sp], lambda x,y : x+y)
        mpi.barrier()

        if self.symm_op != 0: dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for sp in dens_mat[icrsh]:
                    if self.rot_mat_time_inv[icrsh] == 1: dens_mat[icrsh][sp] = dens_mat[icrsh][sp].conjugate()
                    dens_mat[icrsh][sp] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh][sp]),
                                                    self.rot_mat[icrsh] )

        return dens_mat


    # For simple dft input, get crystal field splittings.
    def eff_atomic_levels(self):
        r"""
        Calculates the effective local Hamiltonian required as an input for
        the Hubbard I Solver.
        The local Hamiltonian (effective atomic levels) is calculated by
        projecting the on-site Bloch Hamiltonian:

        .. math:: H^{loc}_{m m'} = \sum_{k} P_{m \nu}(k) H_{\nu\nu'}(k) P^{*}_{\nu' m'}(k),

        where

        .. math:: H_{\nu\nu'}(k) = [\epsilon_{\nu k} - h_{z} \sigma_{z}] \delta_{\nu\nu'}.

        Parameters
        ----------
        None

        Returns
        -------
        eff_atlevels : gf_struct_solver like
                       Effective local Hamiltonian :math:`H^{loc}_{m m'}` for each correlated shell.

        """

        # define matrices for inequivalent shells:
        eff_atlevels = [ {} for ish in range(self.n_inequiv_shells) ]
        for ish in range(self.n_inequiv_shells):
            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                eff_atlevels[ish][sp] = numpy.identity(self.corr_shells[self.inequiv_to_corr[ish]]['dim'], numpy.complex_)
                eff_atlevels[ish][sp] *= -self.chemical_potential
                eff_atlevels[ish][sp] -= self.dc_imp[self.inequiv_to_corr[ish]][sp]

        # sum over k:
        if not hasattr(self,"Hsumk"):
            # calculate the sum over k. Does not depend on mu, so do it only once:
            self.Hsumk = [ {} for icrsh in range(self.n_corr_shells) ]
            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh]['dim']
                for sp in self.spin_block_names[self.corr_shells[icrsh]['SO']]:
                    self.Hsumk[icrsh][sp] = numpy.zeros([dim,dim],numpy.complex_)
                for isp, sp in enumerate(self.spin_block_names[self.corr_shells[icrsh]['SO']]):
                    ind = self.spin_names_to_ind[self.corr_shells[icrsh]['SO']][sp]
                    for ik in range(self.n_k):
                        n_orb = self.n_orbitals[ik,ind]
                        MMat = numpy.identity(n_orb, numpy.complex_)
                        MMat = self.hopping[ik,ind,0:n_orb,0:n_orb] - (1-2*isp) * self.h_field * MMat
                        projmat = self.proj_mat[ik,ind,icrsh,0:dim,0:n_orb]
                        self.Hsumk[icrsh][sp] += self.bz_weights[ik] * numpy.dot( numpy.dot(projmat,MMat),
                                                                                  projmat.conjugate().transpose() )
            # symmetrisation:
            if self.symm_op != 0: self.Hsumk = self.symmcorr.symmetrize(self.Hsumk)

            # Rotate to local coordinate system:
            if self.use_rotations:
                for icrsh in range(self.n_corr_shells):
                    for sp in self.Hsumk[icrsh]:
                        if self.rot_mat_time_inv[icrsh] == 1: self.Hsumk[icrsh][sp] = self.Hsumk[icrsh][sp].conjugate()
                        self.Hsumk[icrsh][sp] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),self.Hsumk[icrsh][sp]) ,
                                                           self.rot_mat[icrsh] )

        # add to matrix:
        for ish in range(self.n_inequiv_shells):
            for sp in eff_atlevels[ish]:
                eff_atlevels[ish][sp] += self.Hsumk[self.inequiv_to_corr[ish]][sp]


        return eff_atlevels


    def init_dc(self):
        r"""
        Initializes the double counting terms.

        Parameters
        ----------
        None

        """
        self.dc_imp = [ {} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            dim = self.corr_shells[icrsh]['dim']
            spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]
            for sp in spn: self.dc_imp[icrsh][sp] = numpy.zeros([dim,dim],numpy.float_)
        self.dc_energ = [0.0 for icrsh in range(self.n_corr_shells)]


    def set_dc(self,dc_imp,dc_energ):
        r"""
        Sets double counting corrections to given values.
        
        Parameters
        ----------
        dc_imp : gf_struct_sumk like
                 Double-counting self-energy term.
        dc_energ : list of floats
                   Double-counting energy corrections for each correlated shell. 
                 
        """

        self.dc_imp = dc_imp
        self.dc_energ = dc_energ


    def calc_dc(self,dens_mat,orb=0,U_interact=None,J_hund=None,use_dc_formula=0,use_dc_value=None):
        r"""
        Calculates and sets the double counting corrections.

        If 'use_dc_value' is provided the double-counting term is uniformly initialized
        with this constant and 'U_interact' and 'J_hund' are ignored.

        If 'use_dc_value' is None the correction is evaluated according to 
        one of the following formulae:

        * use_dc_formula = 0: fully-localised limit (FLL)
        * use_dc_formula = 1: Held's formula, i.e. mean-field formula for the Kanamori
                              type of the interaction Hamiltonian
        * use_dc_formula = 2: around mean-field (AMF)

        Note that FLL and AMF formulae were derived assuming a full Slater-type interaction
        term and should be thus used accordingly. For the Kanamori-type interaction
        one should use formula 1.

        The double-counting self-energy term is stored in `self.dc_imp` and the energy
        correction in `self.dc_energ`.

        Parameters
        ----------
        dens_mat : gf_struct_solver like
                   Density matrix for the specified correlated shell.
        orb : int, optional
              Index of an inequivalent shell.
        U_interact : float, optional
                     Value of interaction parameter `U`.
        J_hund : float, optional
                 Value of interaction parameter `J`.
        use_dc_formula : int, optional
                         Type of double-counting correction (see description).
        use_dc_value : float, optional
                       Value of the double-counting correction. If specified
                       `U_interact`, `J_hund` and `use_dc_formula` are ignored.

        """

        for icrsh in range(self.n_corr_shells):

            ish = self.corr_to_inequiv[icrsh]    # ish is the index of the inequivalent shell corresponding to icrsh
            if ish != orb: continue # ignore this orbital
            dim = self.corr_shells[icrsh]['dim'] #*(1+self.corr_shells[icrsh]['SO'])
            spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]

            Ncr = { sp: 0.0 for sp in spn }
            for block,inner in self.gf_struct_solver[ish].iteritems():
                bl = self.solver_to_sumk_block[ish][block]
                Ncr[bl] += dens_mat[block].real.trace()
            Ncrtot = sum(Ncr.itervalues())
            for sp in spn: 
                self.dc_imp[icrsh][sp] = numpy.identity(dim,numpy.float_)
                if self.SP == 0: # average the densities if there is no SP:
                    Ncr[sp] = Ncrtot / len(spn)
                elif self.SP == 1 and self.SO == 1: # correction for SO: we have only one block in this case, but in DC we need N/2
                    Ncr[sp] = Ncrtot / 2.0

            if use_dc_value is None:

                if U_interact is None and J_hund is None: raise ValueError, "set_dc: either provide U_interact and J_hund or set use_dc_value to dc value."

                if use_dc_formula == 0: # FLL

                    self.dc_energ[icrsh] = U_interact / 2.0 * Ncrtot * (Ncrtot-1.0)
                    for sp in spn:
                        Uav = U_interact*(Ncrtot-0.5) - J_hund*(Ncr[sp] - 0.5)
                        self.dc_imp[icrsh][sp] *= Uav
                        self.dc_energ[icrsh]  -= J_hund / 2.0 * (Ncr[sp]) * (Ncr[sp]-1.0)
                        mpi.report("DC for shell %(icrsh)i and block %(sp)s = %(Uav)f"%locals())

                elif use_dc_formula == 1: # Held's formula, with U_interact the interorbital onsite interaction

                    self.dc_energ[icrsh] = (U_interact + (dim-1)*(U_interact-2.0*J_hund) + (dim-1)*(U_interact-3.0*J_hund))/(2*dim-1) / 2.0 * Ncrtot * (Ncrtot-1.0)
                    for sp in spn:
                        Uav =(U_interact + (dim-1)*(U_interact-2.0*J_hund) + (dim-1)*(U_interact-3.0*J_hund))/(2*dim-1) * (Ncrtot-0.5)
                        self.dc_imp[icrsh][sp] *= Uav
                        mpi.report("DC for shell %(icrsh)i and block %(sp)s = %(Uav)f"%locals())

                elif use_dc_formula == 2: # AMF

                    self.dc_energ[icrsh] = 0.5 * U_interact * Ncrtot * Ncrtot
                    for sp in spn:
                        Uav = U_interact*(Ncrtot - Ncr[sp]/dim) - J_hund * (Ncr[sp] - Ncr[sp]/dim)
                        self.dc_imp[icrsh][sp] *= Uav
                        self.dc_energ[icrsh] -= (U_interact + (dim-1)*J_hund)/dim * 0.5 * Ncr[sp] * Ncr[sp]
                        mpi.report("DC for shell %(icrsh)i and block %(sp)s = %(Uav)f"%locals())

                mpi.report("DC energy for shell %s = %s"%(icrsh,self.dc_energ[icrsh]))

            else: # use value provided for user to determine dc_energ and dc_imp

                self.dc_energ[icrsh] = use_dc_value * Ncrtot
                for sp in spn: self.dc_imp[icrsh][sp] *= use_dc_value

                mpi.report("DC for shell %(icrsh)i = %(use_dc_value)f"%locals())
                mpi.report("DC energy = %s"%self.dc_energ[icrsh])


    def add_dc(self,iw_or_w="iw"):
        r"""
        Subtracts the double counting term from the impurity self energy.
        
        Parameters
        ----------
        iw_or_w : string, optional

                  - `iw_or_w` = 'iw' for a imaginary-frequency self-energy
                  - `iw_or_w` = 'w' for a real-frequency self-energy

        Returns
        -------
        sigma_minus_dc : gf_struct_sumk like
                         Self-energy with a subtracted double-counting term.

        """

        # Be careful: Sigma_imp is already in the global coordinate system!!
        sigma_minus_dc = [s.copy() for s in getattr(self,"Sigma_imp_"+iw_or_w)]
        for icrsh in range(self.n_corr_shells):
            for bname,gf in sigma_minus_dc[icrsh]:
                # Transform dc_imp to global coordinate system
                dccont = numpy.dot(self.rot_mat[icrsh],numpy.dot(self.dc_imp[icrsh][bname],self.rot_mat[icrsh].conjugate().transpose()))
                sigma_minus_dc[icrsh][bname] -= dccont

        return sigma_minus_dc


    def symm_deg_gf(self,gf_to_symm,orb):
        r"""
        Averages a GF over degenerate shells.

        Degenerate shells of an inequivalent correlated shell are defined by
        `self.deg_shells`. This function enforces corresponding degeneracies
        in the input GF.

        Parameters
        ----------
        gf_to_symm : gf_struct_solver like
                     Input GF.
        orb : int
              Index of an inequivalent shell.

        """

        for degsh in self.deg_shells[orb]:
            ss = gf_to_symm[degsh[0]].copy()
            ss.zero()
            n_deg = len(degsh)
            for bl in degsh: ss += gf_to_symm[bl] / (1.0*n_deg)
            for bl in degsh: gf_to_symm[bl] << ss


    def total_density(self, mu=None, with_Sigma=True, with_dc=True):
        r"""
        Calculates the total charge within the energy window for a given chemical potential. 
        The chemical potential is either given by parameter `mu` or, if it is not specified,
        taken from `self.chemical_potential`.

        The total charge is calculated from the trace of the GF in the Bloch basis.
        By deafult, a full interacting GF is used. To use the non-interacting GF, set
        parameter `with_Sigma = False`.

        The number of bands within the energy windows generally depends on `k`. The trace is
        therefore calculated separately for each `k`-point.

        Since in general n_orbitals depends on k, the calculation is done in the following order:
        ..math:: n_{tot} = \sum_{k} n(k),
          with
        ..math:: n(k) = Tr G_{\nu\nu'}(k, i\omega_{n}).

        The calculation is done in the global coordinate system, if distinction is made between local/global.

        Parameters
        ----------
        mu : float, optional
             Input chemical potential. If not specified, `self.chemical_potential` is used instead.
        with_Sigma : boolean, optional
             If `True` the full interacing GF is evaluated, otherwise the self-energy is not
             included and the charge would correspond to a non-interacting system.
        with_dc : boolean, optional
             Whether or not to subtract the double-counting term from the self-energy.

        Returns
        -------
        dens : float
               Total charge :math:`n_{tot}`.

        """
        if mu is None: mu = self.chemical_potential
        dens = 0.0
        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):
            G_latt_iw = self.lattice_gf(ik=ik, mu=mu, iw_or_w="iw", with_Sigma=with_Sigma, with_dc=with_dc)
            dens += self.bz_weights[ik] * G_latt_iw.total_density()
        # collect data from mpi:
        dens = mpi.all_reduce(mpi.world, dens, lambda x,y : x+y)
        mpi.barrier()

        return dens


    def set_mu(self,mu):
        r"""
        Sets a new chemical potential.

        Parameters
        ----------
        mu : float
             New value of the chemical potential.

        """
        self.chemical_potential = mu


    def calc_mu(self, precision=0.01):
        r"""
        Searches for the chemical potential that gives the DFT total charge.
        A simple bisection method is used.

        Parameters
        ----------
        precision : float, optional
                    A desired precision of the resulting total charge.

        Returns
        -------
        mu : float
             Value of the chemical potential giving the DFT total charge
             within specified precision.

        """
        F = lambda mu : self.total_density(mu=mu)
        density = self.density_required - self.charge_below

        self.chemical_potential = dichotomy.dichotomy(function = F,
                                         x_init = self.chemical_potential, y_value = density,
                                         precision_on_y = precision, delta_x = 0.5, max_loops = 100, 
                                         x_name = "Chemical Potential", y_name = "Total Density",
                                         verbosity = 3)[0]

        return self.chemical_potential
 

    def calc_density_correction(self, filename=None, dm_type='wien2k'):
        r"""
        Calculates the charge density correction and stores it into a file.
        
        The charge density correction is needed for charge-self-consistent DFT+DMFT calculations.
        It represents a density matrix of the interacting system defined in Bloch basis
        and it is calculated from the sum over Matsubara frequecies of the full GF,

        ..math:: N_{\nu\nu'}(k) = \sum_{i\omega_{n}} G_{\nu\nu'}(k, i\omega_{n})
        
        The density matrix for every `k`-point is stored into a file.

        Parameters
        ----------
        filename : string
                   Name of the file to store the charge density correction.

        Returns
        -------
        (deltaN, dens) : tuple
                         Returns a tuple containing the density matrix `deltaN` and
                         the corresponing total charge `dens`.

        """

#        assert type(filename) == StringType, "calc_density_correction: filename has to be a string!"
        assert dm_type in ('vasp', 'wien2k'), "'dm_type' must be either 'vasp' or 'wienk'"

        if filename is None:
            if dm_type == 'wien2k':
                filename = 'dens_mat.dat'
            elif dm_type == 'vasp':
                filename = 'GAMMA'
        
        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        dens = {sp: 0.0 for sp in spn}
        band_en_correction = 0.0

# Fetch Fermi weights and energy window band indices
        if dm_type == 'vasp':
            fermi_weights = 0
            band_window = 0
            if mpi.is_master_node():
                ar = HDFArchive(self.hdf_file,'r')
                fermi_weights = ar['dft_misc_input']['dft_fermi_weights']
                band_window = ar['dft_misc_input']['band_window']
                del ar
            fermi_weights = mpi.bcast(fermi_weights)
            band_window = mpi.bcast(band_window)

# Convert Fermi weights to a density matrix
            dens_mat_dft = {}
            for sp in spn:
                dens_mat_dft[sp] = [fermi_weights[ik, ntoi[sp], :].astype(numpy.complex_) for ik in xrange(self.n_k)]


        # Set up deltaN:
        deltaN = {}
        for sp in spn:
            deltaN[sp] = [numpy.zeros([self.n_orbitals[ik,ntoi[sp]],self.n_orbitals[ik,ntoi[sp]]], numpy.complex_) for ik in range(self.n_k)]

        ikarray = numpy.array(range(self.n_k))
        for ik in mpi.slice_array(ikarray):
            G_latt_iw = self.lattice_gf(ik = ik, mu = self.chemical_potential, iw_or_w = "iw")
            for bname,gf in G_latt_iw:
                deltaN[bname][ik] = G_latt_iw[bname].density()
                dens[bname] += self.bz_weights[ik] * G_latt_iw[bname].total_density()
                if dm_type == 'vasp':
# In 'vasp'-mode subtract the DFT density matrix
                    nb = self.n_orbitals[ik, ntoi[bname]]
                    diag_inds = numpy.diag_indices(nb)
                    deltaN[bname][ik][diag_inds] -= dens_mat_dft[bname][ik][:nb]
                    dens[bname] -= self.bz_weights[ik] * dens_mat_dft[bname][ik].sum().real
                    isp = ntoi[bname]
                    b1, b2 = self.band_window[isp][ik, :2]
                    nb = b2 - b1 + 1
                    assert nb == self.n_orbitals[ik, ntoi[bname]], "Number of bands is inconsistent at ik = %s"%(ik)
                    band_en_correction += numpy.dot(deltaN[bname][ik], self.hopping[ik, isp, :nb, :nb]).trace().real * self.bz_weights[ik]


        # mpi reduce:
        for bname in deltaN:
            for ik in range(self.n_k):
                deltaN[bname][ik] = mpi.all_reduce(mpi.world, deltaN[bname][ik], lambda x,y : x+y)
            dens[bname] = mpi.all_reduce(mpi.world, dens[bname], lambda x,y : x+y)
        mpi.barrier()
        band_en_correction = mpi.all_reduce(mpi.world, band_en_correction, lambda x,y : x+y)

        # now save to file:
        if dm_type == 'wien2k':
            if mpi.is_master_node():
                if self.SP == 0:
                    f = open(filename,'w')
                else:
                    f = open(filename+'up','w')
                    f1 = open(filename+'dn','w')
                # write chemical potential (in Rydberg):
                f.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
                if self.SP != 0: f1.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
                # write beta in rydberg-1
                f.write("%.14f\n"%(G_latt_iw.mesh.beta*self.energy_unit))
                if self.SP != 0: f1.write("%.14f\n"%(G_latt_iw.mesh.beta*self.energy_unit))

                if self.SP == 0: # no spin-polarization

                    for ik in range(self.n_k):
                        f.write("%s\n"%self.n_orbitals[ik,0])
                        for inu in range(self.n_orbitals[ik,0]):
                            for imu in range(self.n_orbitals[ik,0]):
                                valre = (deltaN['up'][ik][inu,imu].real + deltaN['down'][ik][inu,imu].real) / 2.0
                                valim = (deltaN['up'][ik][inu,imu].imag + deltaN['down'][ik][inu,imu].imag) / 2.0
                                f.write("%.14f  %.14f "%(valre,valim))
                            f.write("\n")
                        f.write("\n")
                    f.close()

                elif self.SP == 1: # with spin-polarization

                    # dict of filename: (spin index, block_name)
                    if self.SO == 0: to_write = {f: (0, 'up'), f1: (1, 'down')}
                    if self.SO == 1: to_write = {f: (0, 'ud'), f1: (0, 'ud')}
                    for fout in to_write.iterkeys():
                        isp, sp = to_write[fout]
                        for ik in range(self.n_k):
                            fout.write("%s\n"%self.n_orbitals[ik,isp])
                            for inu in range(self.n_orbitals[ik,isp]):
                                for imu in range(self.n_orbitals[ik,isp]):
                                    fout.write("%.14f  %.14f "%(deltaN[sp][ik][inu,imu].real,deltaN[sp][ik][inu,imu].imag))
                                fout.write("\n")
                            fout.write("\n")
                        fout.close()
        elif dm_type == 'vasp':
            assert self.SP == 0, "Spin-polarized density matrix is not implemented"

            if mpi.is_master_node():
                with open(filename, 'w') as f:
                    f.write(" %i  -1  ! Number of k-points, default number of bands\n"%(self.n_k))
                    for ik in xrange(self.n_k):
                        ib1 = band_window[0][ik, 0]
                        ib2 = band_window[0][ik, 1]
                        f.write(" %i  %i  %i\n"%(ik + 1, ib1, ib2))
                        for inu in xrange(self.n_orbitals[ik, 0]):
                            for imu in xrange(self.n_orbitals[ik, 0]):
                                valre = (deltaN['up'][ik][inu, imu].real + deltaN['down'][ik][inu, imu].real) / 2.0
                                valim = (deltaN['up'][ik][inu, imu].imag + deltaN['down'][ik][inu, imu].imag) / 2.0
                                f.write(" %.14f  %.14f"%(valre, valim))
                            f.write("\n")
        else:
            raise NotImplementedError("Unknown density matrix type: '%s'"%(dm_type))

        res = deltaN, dens

        if dm_type == 'vasp':
            res += (band_en_correction,)

        return res


################
# FIXME LEAVE UNDOCUMENTED
################

    def calc_dc_for_density(self,orb,dc_init,dens_mat,density=None,precision=0.01):
        """Searches for DC in order to fulfill charge neutrality.
           If density is given, then DC is set such that the LOCAL charge of orbital
           orb coincides with the given density."""

        def F(dc):
            self.calc_dc(dens_mat = dens_mat, U_interact = 0, J_hund = 0, orb = orb, use_dc_value = dc)
            if dens_req is None:
                return self.total_density(mu = mu)
            else:
                return self.extract_G_loc()[orb].total_density()

        if density is None: density = self.density_required - self.charge_below

        dc = dichotomy.dichotomy(function = F,
                                    x_init = dc_init, y_value = density,
                                    precision_on_y = precision, delta_x = 0.5,
                                    max_loops = 100, x_name = "Double Counting", y_name= "Total Density",
                                    verbosity = 3)[0]

        return dc


    def check_projectors(self):
        """Calculated the density matrix from projectors (DM = P Pdagger) to check that it is correct and 
           specifically that it matches DFT."""
        dens_mat = [numpy.zeros([self.corr_shells[icrsh]['dim'],self.corr_shells[icrsh]['dim']],numpy.complex_)
                   for icrsh in range(self.n_corr_shells)]

        for ik in range(self.n_k):
            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh]['dim']
                n_orb = self.n_orbitals[ik,0]
                projmat = self.proj_mat[ik,0,icrsh,0:dim,0:n_orb]
                dens_mat[icrsh][:,:] += numpy.dot(projmat, projmat.transpose().conjugate()) * self.bz_weights[ik]

        if self.symm_op != 0: dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                if self.rot_mat_time_inv[icrsh] == 1: dens_mat[icrsh] = dens_mat[icrsh].conjugate()
                dens_mat[icrsh] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh]),
                                            self.rot_mat[icrsh] )

        return dens_mat


    def sorts_of_atoms(self,shells):
        """
        Determine the number of inequivalent sorts.
        """
        sortlst = [ shells[i]['sort'] for i in range(len(shells)) ]
        n_sorts = len(set(sortlst))
        return n_sorts


    def number_of_atoms(self,shells):
        """
        Determine the number of inequivalent atoms.
        """
        atomlst = [ shells[i]['atom'] for i in range(len(shells)) ]
        n_atoms = len(set(atomlst))
        return n_atoms
