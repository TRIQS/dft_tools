 
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

class SumkDFT:
    """This class provides a general SumK method for combining ab-initio code and pytriqs."""


    def __init__(self, hdf_file, h_field = 0.0, use_dft_blocks = False, 
	         dft_data = 'dft_input', symmcorr_data = 'dft_symmcorr_input', parproj_data = 'dft_parproj_input', 
                 symmpar_data = 'dft_symmpar_input', bands_data = 'dft_bands_input', transp_data = 'dft_transp_input',
                 misc_data = 'dft_misc_input'):
        """
        Initialises the class from data previously stored into an HDF5
        """

        if not type(hdf_file) == StringType:
            mpi.report("Give a string for the HDF5 filename to read the input!")
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
# HDF5 FUNCTIONS
################

    def read_input_from_hdf(self, subgrp, things_to_read):
        """
        Reads data from the HDF file
        """

        value_read = True
        # initialise variables on all nodes to ensure mpi broadcast works at the end
        for it in things_to_read: setattr(self,it,0)
        subgroup_present = 0

        if mpi.is_master_node():
            ar = HDFArchive(self.hdf_file,'a')
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
                if (len(things_to_read) != 0): mpi.report("Loading failed: No %s subgroup in HDF5!"%subgrp)
                subgroup_present = False
                value_read = False
            del ar
        # now do the broadcasting:
        for it in things_to_read: setattr(self,it,mpi.bcast(getattr(self,it)))
        subgroup_present = mpi.bcast(subgroup_present)
        value_read = mpi.bcast(value_read)

        return subgroup_present, value_read


    def save(self, things_to_save, subgrp='user_data'):
        """Saves given quantities into the subgroup ('user_data' by default) of the HDF5 archive"""

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
        """Loads given quantities from the subgroup ('user_data' by default) of the HDF5 archive"""

        if not (mpi.is_master_node()): return # do nothing on nodes
        ar = HDFArchive(self.hdf_file,'a')
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
        """Downfolding a block of the Greens function"""
  
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
        """Upfolding a block of the Greens function"""
  
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
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

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
        """Calculates the lattice Green function from the DFT hopping and the self energy at k-point number ik
           and chemical potential mu."""

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
        n_iw = 1025 # Default number of Matsubara frequencies

        # Are we including Sigma?
        if with_Sigma:
            if with_dc: sigma_minus_dc = self.add_dc(iw_or_w)
            Sigma_imp = getattr(self,"Sigma_imp_"+iw_or_w)
            if iw_or_w == "iw":
                beta = Sigma_imp[0].mesh.beta   # override beta if Sigma_iw is present
                n_iw = len(Sigma_imp[0].mesh)
        else:
            if (iw_or_w == "w") and (mesh is None):
                raise ValueError, "lattice_gf: Give the mesh=(om_min,om_max,n_points) for the lattice GfReFreq."

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
               glist = lambda : [ GfImFreq(indices=inner,beta=beta,n_points=n_iw) for block,inner in gf_struct]
            elif iw_or_w == "w":
                 if with_Sigma:
                    glist = lambda : [ GfReFreq(indices=inner,mesh=Sigma_imp[0].mesh) for block,inner in gf_struct]
                 else:
                    glist = lambda : [ GfReFreq(indices=inner,window=(mesh[0],mesh[1]),n_points=mesh[2]) for block,inner in gf_struct]
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


    def put_Sigma(self, Sigma_imp):
        """Sets the impurity self energies for inequivalent atoms into the class, respects the multiplicity of the atoms."""

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
        """
        Extracts the local downfolded Green function at the chemical potential of the class.
        At the end, the local G is rotated from the global coordinate system to the local system.
        if with_Sigma = False: Sigma is not included => non-interacting local GF
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
        """ Determines the Green's function block structure from simple point integration."""

        self.gf_struct_solver = [ {} for ish in range(self.n_inequiv_shells) ]
        self.sumk_to_solver = [ {} for ish in range(self.n_inequiv_shells) ]
        self.solver_to_sumk = [ {} for ish in range(self.n_inequiv_shells) ]
        self.solver_to_sumk_block = [ {} for ish in range(self.n_inequiv_shells) ]

        if dm is None: dm = self.density_matrix(method = 'using_point_integration')
        dens_mat = [ dm[self.inequiv_to_corr[ish]] for ish in range(self.n_inequiv_shells) ]

        if include_shells is None: include_shells = range(self.n_inequiv_shells)
        for ish in include_shells:


            for sp in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]]['SO']]:
                dmbool = (abs(dens_mat[ish][sp]) > threshold)  # gives an index list of entries larger that threshold

                # Determine off-diagonal entries in upper triangular part of density matrix
                offdiag = []
                for i in range(len(dmbool)):
                    for j in range(i+1,len(dmbool)):
                        if dmbool[i,j]: offdiag.append([i,j])

                # Determine the number of non-hybridising blocks in the gf
                num_blocs = len(dmbool)
                blocs = [ [i] for i in range(num_blocs) ]
                for i in range(len(offdiag)):
                    for j in range(len(blocs[offdiag[i][1]])): blocs[offdiag[i][0]].append(blocs[offdiag[i][1]][j])
                    del blocs[offdiag[i][1]]
                    for j in range(i+1,len(offdiag)):
                        if offdiag[j][0] == offdiag[i][1]: offdiag[j][0] = offdiag[i][0]
                        if offdiag[j][1] == offdiag[i][1]: offdiag[j][1] = offdiag[i][0]
                        if offdiag[j][0] > offdiag[i][1]: offdiag[j][0] -= 1
                        if offdiag[j][1] > offdiag[i][1]: offdiag[j][1] -= 1
                        offdiag[j].sort()
                    num_blocs -= 1

                # Set the gf_struct for the solver accordingly
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
        """Calculate density matrices in one of two ways:
            if 'using_gf': First get lattice gf (g_loc is not set up), then density matrix.
                              It is useful for Hubbard I, and very quick.
                              No assumption on the hopping structure is made (ie diagonal or not).
            if 'using_point_integration': Only works for diagonal hopping matrix (true in wien2k).
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
        """Calculates the effective atomic levels needed as input for the Hubbard I Solver."""

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
        """ Initialise the double counting terms to have the correct shape."""

        self.dc_imp = [ {} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            dim = self.corr_shells[icrsh]['dim']
            spn = self.spin_block_names[self.corr_shells[icrsh]['SO']]
            for sp in spn: self.dc_imp[icrsh][sp] = numpy.zeros([dim,dim],numpy.float_)
        self.dc_energ = [0.0 for icrsh in range(self.n_corr_shells)]


    def set_dc(self,dc_imp,dc_energ):
        """Sets double counting terms dc_imp and dc_energ to known values."""

        self.dc_imp = dc_imp
        self.dc_energ = dc_energ


    def calc_dc(self,dens_mat,orb=0,U_interact=None,J_hund=None,use_dc_formula=0,use_dc_value=None):
        """Sets the double counting corrections in the correct form for inequiv orbital orb:
           1) either using U_interact, J_hund and 
           use_dc_formula = 0: fully-localised limit (FLL),
           use_dc_formula = 1: Held's formula,
           use_dc_formula = 2: around mean-field (AMF).
           2) or using a given dc value in use_dc_value.
           Be sure that you are using the correct interaction Hamiltonian!"""

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
        """Substracts the double counting term from the impurity self energy."""

        # Be careful: Sigma_imp is already in the global coordinate system!!
        sigma_minus_dc = [s.copy() for s in getattr(self,"Sigma_imp_"+iw_or_w)]
        for icrsh in range(self.n_corr_shells):
            for bname,gf in sigma_minus_dc[icrsh]:
                # Transform dc_imp to global coordinate system
                dccont = numpy.dot(self.rot_mat[icrsh],numpy.dot(self.dc_imp[icrsh][bname],self.rot_mat[icrsh].conjugate().transpose()))
                sigma_minus_dc[icrsh][bname] -= dccont

        return sigma_minus_dc


    def symm_deg_gf(self,gf_to_symm,orb):
        """Symmetrises a GF for the given degenerate shells self.deg_shells"""

        for degsh in self.deg_shells[orb]:
            ss = gf_to_symm[degsh[0]].copy()
            ss.zero()
            n_deg = len(degsh)
            for bl in degsh: ss += gf_to_symm[bl] / (1.0*n_deg)
            for bl in degsh: gf_to_symm[bl] << ss


    def total_density(self, mu=None, with_Sigma=True, with_dc=True):
        """
        Calculates the total charge for the energy window for a given chemical potential mu. 
        Since in general n_orbitals depends on k, the calculation is done in the following order:
        G_aa'(k,iw) -> n(k) = Tr G_aa'(k,iw) -> sum_k n_k

        The calculation is done in the global coordinate system, if distinction is made between local/global!
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
        """Sets a new chemical potential"""

        self.chemical_potential = mu


    def calc_mu(self, precision = 0.01):
        """
        Searches for mu in order to give the desired charge
        A desired precision can be specified in precision.
        """

        F = lambda mu : self.total_density(mu = mu)
        density = self.density_required - self.charge_below

        self.chemical_potential = dichotomy.dichotomy(function = F,
                                         x_init = self.chemical_potential, y_value = density,
                                         precision_on_y = precision, delta_x = 0.5, max_loops = 100, 
                                         x_name = "Chemical Potential", y_name = "Total Density",
                                         verbosity = 3)[0]

        return self.chemical_potential
 

    def calc_density_correction(self,filename = 'dens_mat.dat'):
        """ Calculates the density correction in order to feed it back to the DFT calculations."""

        assert type(filename) == StringType, "calc_density_correction: filename has to be a string!"

        ntoi = self.spin_names_to_ind[self.SO]
        spn = self.spin_block_names[self.SO]
        dens = {sp: 0.0 for sp in spn}

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

        # mpi reduce:
        for bname in deltaN:
            for ik in range(self.n_k):
                deltaN[bname][ik] = mpi.all_reduce(mpi.world, deltaN[bname][ik], lambda x,y : x+y)
            dens[bname] = mpi.all_reduce(mpi.world, dens[bname], lambda x,y : x+y)
        mpi.barrier()

        # now save to file:
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

        return deltaN, dens

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
