 
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

class SumkLDA:
    """This class provides a general SumK method for combining ab-initio code and pytriqs."""


    def __init__(self, hdf_file, mu = 0.0, h_field = 0.0, use_lda_blocks = False, 
	         lda_data = 'lda_input', symmcorr_data = 'lda_symmcorr_input', parproj_data = 'lda_parproj_input', 
                 symmpar_data = 'lda_symmpar_input', bands_data = 'lda_bands_input', lda_output = 'lda_output'):
        """
        Initialises the class from data previously stored into an HDF5
        """

        if not type(hdf_file) == StringType:
            mpi.report("Give a string for the HDF5 filename to read the input!")
        else:
            self.hdf_file = hdf_file
            self.lda_data = lda_data
            self.symmcorr_data = symmcorr_data
            self.parproj_data = parproj_data
            self.symmpar_data = symmpar_data
            self.bands_data = bands_data
            self.lda_output = lda_output
            self.G_upfold = None
            self.h_field = h_field

            # read input from HDF:
            things_to_read = ['energy_unit','n_k','k_dep_projection','SP','SO','charge_below','density_required',
                              'symm_op','n_shells','shells','n_corr_shells','corr_shells','use_rotations','rot_mat',
                              'rot_mat_time_inv','n_reps','dim_reps','T','n_orbitals','proj_mat','bz_weights','hopping',
                              'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
            self.subgroup_present, self.value_read = self.read_input_from_hdf(subgrp = self.lda_data, things_to_read = things_to_read)

            if self.SO and (abs(self.h_field) > 0.000001):
                self.h_field = 0.0
                mpi.report("For SO, the external magnetic field is not implemented, setting it to 0!")

            self.spin_block_names = [ ['up','down'], ['ud'] ]
            self.n_spin_blocks = [2,1]
            # convert spin_block_names to indices -- if spin polarized, differentiate up and down blocks 
            self.spin_names_to_ind = [{}, {}]
            for iso in range(2): # SO = 0 or 1
                for ibl in range(self.n_spin_blocks[iso]):
                    self.spin_names_to_ind[iso][self.spin_block_names[iso][ibl]] = ibl * self.SP

            # GF structure used for the local things in the k sums
            # Most general form allowing for all hybridisation, i.e. largest blocks possible
            self.gf_struct_sumk = [ [ (b, range( self.corr_shells[icrsh][3])) for b in self.spin_block_names[self.corr_shells[icrsh][4]] ]
                                   for icrsh in range(self.n_corr_shells) ]

            #-----
            # If these quantities are not in HDF, set them up
            optional_things = ['gf_struct_solver','sumk_to_solver','solver_to_sumk','solver_to_sumk_block','chemical_potential','dc_imp','dc_energ','deg_shells']
            self.subgroup_present, self.value_read = self.read_input_from_hdf(subgrp = self.lda_output, things_to_read = [], 
                                                                              optional_things = optional_things)
            if (not self.subgroup_present) or (not self.value_read['gf_struct_solver']):
                # No gf_struct was stored in HDF, so first set a standard one:
                self.gf_struct_solver = [ dict([ (b, range(self.corr_shells[self.inequiv_to_corr[ish]][3]) )
                                                        for b in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]][4]] ])
                                          for ish in range(self.n_inequiv_shells)
                                        ]
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

            if (not self.subgroup_present) or (not self.value_read['dc_imp']):
                self.__init_dc() # initialise the double counting

            if (not self.subgroup_present) or (not self.value_read['chemical_potential']):
                self.chemical_potential = mu

            if (not self.subgroup_present) or (not self.value_read['deg_shells']):
                self.deg_shells = [ [] for ish in range(self.n_inequiv_shells)]
            #-----

            if self.symm_op:
                self.symmcorr = Symmetry(hdf_file,subgroup=self.symmcorr_data)

            # Analyse the block structure and determine the smallest blocks, if desired
            if use_lda_blocks: dm = self.analyse_block_structure()

            # Now save new things to HDF5: 
            # FIXME WHAT HAPPENS TO h_field? INPUT TO __INIT__? ADD TO OPTIONAL_THINGS?
            things_to_save = ['chemical_potential','dc_imp','dc_energ','h_field']
            self.save(things_to_save)

################
# HDF5 FUNCTIONS
################

    def read_input_from_hdf(self, subgrp, things_to_read=[], optional_things=[]):
        """
        Reads data from the HDF file
        """

        value_read = True
        # initialise variables on all nodes to ensure mpi broadcast works at the end
        for it in things_to_read: setattr(self,it,0)
        for it in optional_things: setattr(self,it,0)

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

                if value_read and (len(optional_things) > 0):
                    # if successfully read necessary items, read optional things:
                    value_read = {}
                    for it in optional_things:
                        if it in ar[subgrp]:
                            setattr(self,it,ar[subgrp][it])
                            value_read['%s'%it] = True
                        else:
                            value_read['%s'%it] = False
            else:
                if (len(things_to_read) != 0): mpi.report("Loading failed: No %s subgroup in HDF5!"%subgrp)
                subgroup_present = False
                value_read = False

            del ar

        # now do the broadcasting:
        for it in things_to_read: setattr(self,it,mpi.bcast(getattr(self,it)))
        for it in optional_things: setattr(self,it,mpi.bcast(getattr(self,it)))
        subgroup_present = mpi.bcast(subgroup_present)
        value_read = mpi.bcast(value_read)

        return subgroup_present, value_read


    def save(self,things_to_save):
        """Saves some quantities into an HDF5 archive"""

        if not (mpi.is_master_node()): return # do nothing on nodes
        ar = HDFArchive(self.hdf_file,'a')
        if not self.lda_output in ar: ar.create_group(self.lda_output)
        for it in things_to_save: 
            try:
                ar[self.lda_output][it] = getattr(self,it)
            except:
                mpi.report("%s not found, and so not stored."%it)
        del ar

################
# CORE FUNCTIONS
################

    def downfold(self,ik,icrsh,bname,gf_to_downfold,gf_inp):
        """Downfolding a block of the Greens function"""

        gf_downfolded = gf_inp.copy()
        isp = self.spin_names_to_ind[self.SO][bname]       # get spin index for proj. matrices
        dim = self.corr_shells[icrsh][3]
        n_orb = self.n_orbitals[ik,isp]
        projmat = self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb]

        gf_downfolded.from_L_G_R(projmat,gf_to_downfold,projmat.conjugate().transpose())

        return gf_downfolded


    def upfold(self,ik,icrsh,bname,gf_to_upfold,gf_inp):
        """Upfolding a block of the Greens function"""

        gf_upfolded = gf_inp.copy()
        isp = self.spin_names_to_ind[self.SO][bname]       # get spin index for proj. matrices
        dim = self.corr_shells[icrsh][3]
        n_orb = self.n_orbitals[ik,isp]
        projmat = self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb]

        gf_upfolded.from_L_G_R(projmat.conjugate().transpose(),gf_to_upfold,projmat)

        return gf_upfolded


    def rotloc(self,icrsh,gf_to_rotate,direction):
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

        assert ((direction == 'toLocal') or (direction == 'toGlobal')),"Give direction 'toLocal' or 'toGlobal' in rotloc!"

        gf_rotated = gf_to_rotate.copy()

        if direction == 'toGlobal':

            if (self.rot_mat_time_inv[icrsh] == 1) and self.SO:
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rot_mat[icrsh].conjugate(),gf_rotated,self.rot_mat[icrsh].transpose())
            else:
                gf_rotated.from_L_G_R(self.rot_mat[icrsh],gf_rotated,self.rot_mat[icrsh].conjugate().transpose())

        elif direction == 'toLocal':

            if (self.rot_mat_time_inv[icrsh] == 1) and self.SO:
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rot_mat[icrsh].transpose(),gf_rotated,self.rot_mat[icrsh].conjugate())
            else:
                gf_rotated.from_L_G_R(self.rot_mat[icrsh].conjugate().transpose(),gf_rotated,self.rot_mat[icrsh])

        return gf_rotated


    def lattice_gf_matsubara(self,ik,mu,beta=40,with_Sigma=True):
        """Calculates the lattice Green function from the LDA hopping and the self energy at k-point number ik
           and chemical potential mu."""

        ntoi = self.spin_names_to_ind[self.SO]
        bln = self.spin_block_names[self.SO]

        if not hasattr(self,"Sigma_imp"): with_Sigma = False

        if with_Sigma:
            stmp = self.add_dc()
            beta = self.Sigma_imp[0].mesh.beta        # override beta if Sigma is present

        # Do we need to set up G_upfold?
        set_up_G_upfold = False # assume not
        if self.G_upfold is None: # yes if not G_upfold provided
            set_up_G_upfold = True
        else: # yes if inconsistencies present in existing G_upfold
            GFsize = [ gf.N1 for bname,gf in self.G_upfold]
            unchangedsize = all( [ self.n_orbitals[ik,ntoi[bln[ibl]]] == GFsize[ibl]
                                   for ibl in range(self.n_spin_blocks[self.SO]) ] )
            if (not unchangedsize) or (self.G_upfold.mesh.beta != beta): set_up_G_upfold = True

        # Set up G_upfold
        if set_up_G_upfold:
            block_structure = [ range(self.n_orbitals[ik,ntoi[b]]) for b in bln ]
            gf_struct = [ (bln[ibl], block_structure[ibl]) for ibl in range(self.n_spin_blocks[self.SO]) ]
            block_ind_list = [block for block,inner in gf_struct]
            if with_Sigma:
                glist = lambda : [ GfImFreq(indices = inner, mesh = self.Sigma_imp[0].mesh) for block,inner in gf_struct]
            else:
                glist = lambda : [ GfImFreq(indices = inner, beta = beta) for block,inner in gf_struct]
            self.G_upfold = BlockGf(name_list = block_ind_list, block_list = glist(), make_copies=False)
            self.G_upfold.zero()

        self.G_upfold << iOmega_n
        idmat = [numpy.identity(self.n_orbitals[ik,ntoi[bl]],numpy.complex_) for bl in bln]
        M = copy.deepcopy(idmat)
        for ibl in range(self.n_spin_blocks[self.SO]):
            ind = ntoi[bln[ibl]]
            n_orb = self.n_orbitals[ik,ind]
            M[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb] - (idmat[ibl]*mu) - (idmat[ibl] * self.h_field * (1-2*ibl))
        self.G_upfold -= M

        if with_Sigma:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in self.G_upfold: gf -= self.upfold(ik,icrsh,bname,stmp[icrsh][bname],gf)

        self.G_upfold.invert()

	return self.G_upfold


    def simple_point_dens_mat(self):

        ntoi = self.spin_names_to_ind[self.SO]
        bln = self.spin_block_names[self.SO]

        MMat = [numpy.zeros( [self.n_orbitals[0,ntoi[bl]],self.n_orbitals[0,ntoi[bl]]], numpy.complex_) for bl in bln]

        dens_mat = [ {} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            for bl in self.spin_block_names[self.corr_shells[icrsh][4]]:
                dens_mat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]], numpy.complex_)

        ikarray = numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            unchangedsize = all( [ self.n_orbitals[ik,ntoi[bln[ibl]]] == len(MMat[ibl])
                                   for ibl in range(self.n_spin_blocks[self.SO]) ] )

            if not unchangedsize:
                MMat = [numpy.zeros( [self.n_orbitals[ik,ntoi[bl]],self.n_orbitals[ik,ntoi[bl]]], numpy.complex_) for bl in bln]

            for ibl, bl in enumerate(bln):
                ind = ntoi[bl]
                for inu in range(self.n_orbitals[ik,ind]):
                    if (self.hopping[ik,ind,inu,inu] - self.h_field*(1-2*ibl)) < 0.0: # only works for diagonal hopping matrix (true in wien2k)
                        MMat[ibl][inu,inu] = 1.0
                    else:
                        MMat[ibl][inu,inu] = 0.0

            for icrsh in range(self.n_corr_shells):
                for ibl, bn in enumerate(self.spin_block_names[self.corr_shells[icrsh][4]]):
                    isp = self.spin_names_to_ind[self.corr_shells[icrsh][4]][bn]
                    dim = self.corr_shells[icrsh][3]
                    n_orb = self.n_orbitals[ik,isp]
                    projmat = self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb]
                    dens_mat[icrsh][bn] += self.bz_weights[ik] * numpy.dot( numpy.dot(projmat,MMat[ibl]) ,
                                                                           projmat.transpose().conjugate() )

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for bname in dens_mat[icrsh]:
                dens_mat[icrsh][bname] = mpi.all_reduce(mpi.world, dens_mat[icrsh][bname], lambda x,y : x+y)
        mpi.barrier()

        if self.symm_op != 0: dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bn in dens_mat[icrsh]:
                    if self.rot_mat_time_inv[icrsh] == 1: dens_mat[icrsh][bn] = dens_mat[icrsh][bn].conjugate()
                    dens_mat[icrsh][bn] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh][bn]) ,
                                                    self.rot_mat[icrsh] )

        return dens_mat

    # Calculate upfolded gf, then density matrix. No assumption on the structure made here (ie diagonal or not).
    def density_gf(self,beta):
        """Calculates the density without setting up Gloc. It is useful for Hubbard I, and very quick."""

        dens_mat = [ {} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            for bl in self.spin_block_names[self.corr_shells[icrsh][4]]:
                dens_mat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]], numpy.complex_)

        ikarray = numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            G_upfold = self.lattice_gf_matsubara(ik=ik, beta=beta, mu=self.chemical_potential)
            G_upfold *= self.bz_weights[ik]
            dm = G_upfold.density()
            MMat = [dm[bl] for bl in self.spin_block_names[self.SO]]

            for icrsh in range(self.n_corr_shells):
                for ibl, bn in enumerate(self.spin_block_names[self.corr_shells[icrsh][4]]):
                    isp = self.spin_names_to_ind[self.corr_shells[icrsh][4]][bn]
                    dim = self.corr_shells[icrsh][3]
                    n_orb = self.n_orbitals[ik,isp]
                    projmat = self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb]
                    dens_mat[icrsh][bn] += numpy.dot( numpy.dot(projmat,MMat[ibl]),
                                                      projmat.transpose().conjugate() )

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for bname in dens_mat[icrsh]:
                dens_mat[icrsh][bname] = mpi.all_reduce(mpi.world, dens_mat[icrsh][bname], lambda x,y : x+y)
        mpi.barrier()

        if self.symm_op != 0: dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bn in dens_mat[icrsh]:
                    if self.rot_mat_time_inv[icrsh] == 1: dens_mat[icrsh][bn] = dens_mat[icrsh][bn].conjugate()
                    dens_mat[icrsh][bn] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh][bn]) ,
                                                    self.rot_mat[icrsh] )

        return dens_mat



    def analyse_block_structure(self, threshold = 0.00001, include_shells = None, dm = None):
        """ Determines the Green function block structure from simple point integration."""

        self.gf_struct_solver = [ {} for ish in range(self.n_inequiv_shells) ]
        self.sumk_to_solver = [ {} for ish in range(self.n_inequiv_shells) ]
        self.solver_to_sumk = [ {} for ish in range(self.n_inequiv_shells) ]
        self.solver_to_sumk_block = [ {} for ish in range(self.n_inequiv_shells) ]

        if dm is None: dm = self.simple_point_dens_mat()
        dens_mat = [ dm[self.inequiv_to_corr[ish]] for ish in range(self.n_inequiv_shells) ]

        if include_shells is None: include_shells = range(self.n_inequiv_shells)
        for ish in include_shells:

            block_ind_list = [ block for block,inner in self.gf_struct_sumk[self.inequiv_to_corr[ish]] ]
            for block in block_ind_list:
                dm = dens_mat[ish][block]
                dmbool = (abs(dm) > threshold)          # gives an index list of entries larger that threshold

                # Determine off-diagonal entries in upper triangular part of density matrix
                offdiag = []
                for i in range(len(dmbool)):
                    for j in range(i+1,len(dmbool)):
                        if dmbool[i,j]: offdiag.append([i,j])

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

                for i in range(num_blocs):
                    blocs[i].sort()
                    self.gf_struct_solver[ish].update( [('%s_%s'%(block,i),range(len(blocs[i])))] )

                # Construct sumk_to_solver taking (sumk_block, sumk_index) --> (solver_block, solver_inner)
                #       and solver_to_sumk taking (solver_block, solver_inner) --> (sumk_block, sumk_index)
                for i in range(num_blocs):
                    for j in range(len(blocs[i])):
                        block_sumk = block
                        inner_sumk = blocs[i][j]
                        block_solv = '%s_%s'%(block,i)
                        inner_solv = j
                        self.sumk_to_solver[ish][(block_sumk,inner_sumk)] = (block_solv,inner_solv)
                        self.solver_to_sumk[ish][(block_solv,inner_solv)] = (block_sumk,inner_sumk)
                        self.solver_to_sumk_block[ish][block_solv] = block_sumk

            # now calculate degeneracies of orbitals:
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
                            # check if it was already there:
                            ind1 = -1
                            ind2 = -2
                            for n,ind in enumerate(self.deg_shells[ish]):
                                if block1 in ind: ind1 = n
                                if block2 in ind: ind2 = n
                            if (ind1 < 0) and (ind2 >= 0):
                                self.deg_shells[ish][ind2].append(block1)
                            elif (ind1 >= 0) and (ind2 < 0):
                                self.deg_shells[ish][ind1].append(block2)
                            elif (ind1 < 0) and (ind2 < 0):
                                self.deg_shells[ish].append([block1,block2])

        things_to_save = ['gf_struct_solver','sumk_to_solver','solver_to_sumk','solver_to_sumk_block','deg_shells']
        self.save(things_to_save)

        return dens_mat


    def symm_deg_gf(self,gf_to_symm,orb):
        """Symmetrises a GF for the given degenerate shells self.deg_shells"""

        for degsh in self.deg_shells[orb]:
            #loop over degenerate shells:
            ss = gf_to_symm[degsh[0]].copy()
            ss.zero()
            Ndeg = len(degsh)
            for bl in degsh: ss += gf_to_symm[bl] / (1.0*Ndeg)
            for bl in degsh: gf_to_symm[bl] << ss

# for simple dft input, get crystal field splittings. 
    def eff_atomic_levels(self):
        """Calculates the effective atomic levels needed as input for the Hubbard I Solver."""

        # define matrices for inequivalent shells:
        eff_atlevels = [ {} for ish in range(self.n_inequiv_shells) ]
        for ish in range(self.n_inequiv_shells):
            for bn in self.spin_block_names[self.corr_shells[self.inequiv_to_corr[ish]][4]]:
                eff_atlevels[ish][bn] = numpy.identity(self.corr_shells[self.inequiv_to_corr[ish]][3], numpy.complex_)

        # Chemical Potential:
        for ish in range(self.n_inequiv_shells):
            for ii in eff_atlevels[ish]: eff_atlevels[ish][ii] *= -self.chemical_potential

        # double counting term:
        for ish in range(self.n_inequiv_shells):
            for ii in eff_atlevels[ish]:
                eff_atlevels[ish][ii] -= self.dc_imp[self.inequiv_to_corr[ish]][ii]

        # sum over k:
        if not hasattr(self,"Hsumk"):
            # calculate the sum over k. Does not depend on mu, so do it only once:
            self.Hsumk = [ {} for icrsh in range(self.n_corr_shells) ]
            for icrsh in range(self.n_corr_shells):
                for bn in self.spin_block_names[self.corr_shells[icrsh][4]]:
                    dim = self.corr_shells[icrsh][3]  #*(1+self.corr_shells[icrsh][4])
                    self.Hsumk[icrsh][bn] = numpy.zeros([dim,dim],numpy.complex_)

            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh][3]
                for ibl, bn in enumerate(self.spin_block_names[self.corr_shells[icrsh][4]]):
                    isp = self.spin_names_to_ind[self.corr_shells[icrsh][4]][bn]
                    for ik in range(self.n_k):
                        n_orb = self.n_orbitals[ik,isp]
                        MMat = numpy.identity(n_orb, numpy.complex_)
                        MMat = self.hopping[ik,isp,0:n_orb,0:n_orb] - (1-2*ibl) * self.h_field * MMat
                        projmat = self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb]
                        self.Hsumk[icrsh][bn] += self.bz_weights[ik] * numpy.dot( numpy.dot(projmat,MMat),
                                                                                  projmat.conjugate().transpose() )

            # symmetrisation:
            if self.symm_op != 0: self.Hsumk = self.symmcorr.symmetrize(self.Hsumk)

            # Rotate to local coordinate system:
            if self.use_rotations:
                for icrsh in range(self.n_corr_shells):
                    for bn in self.Hsumk[icrsh]:

                        if self.rot_mat_time_inv[icrsh] == 1: self.Hsumk[icrsh][bn] = self.Hsumk[icrsh][bn].conjugate()
                        self.Hsumk[icrsh][bn] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),self.Hsumk[icrsh][bn]) ,
                                                           self.rot_mat[icrsh] )

        # add to matrix:
        for ish in range(self.n_inequiv_shells):
            for bn in eff_atlevels[ish]:
                eff_atlevels[ish][bn] += self.Hsumk[self.inequiv_to_corr[ish]][bn]


        return eff_atlevels



    def __init_dc(self):

        # construct the density matrix dm_imp and double counting arrays
        self.dc_imp = [ {} for icrsh in range(self.n_corr_shells)]
        for icrsh in range(self.n_corr_shells):
            dim = self.corr_shells[icrsh][3]
            for j in range(len(self.gf_struct_sumk[icrsh])):
                self.dc_imp[icrsh]['%s'%self.gf_struct_sumk[icrsh][j][0]] = numpy.zeros([dim,dim],numpy.float_)
        self.dc_energ = [0.0 for icrsh in range(self.n_corr_shells)]



    def set_dc(self,dens_mat,U_interact,J_hund,orb=0,use_dc_formula=0,use_val=None):
        """Sets the double counting term for inequiv orbital orb:
           use_dc_formula=0: LDA+U FLL double counting,
           use_dc_formula=1: Held's formula,
           use_dc_formula=2: AMF.
           Be sure that you are using the correct interaction Hamiltonian!"""


        for icrsh in range(self.n_corr_shells):

            iorb = self.corr_to_inequiv[icrsh]    # iorb is the index of the inequivalent shell corresponding to icrsh

            if iorb != orb: continue # ignore this orbital

            Ncr = {}
            dim = self.corr_shells[icrsh][3] #*(1+self.corr_shells[icrsh][4])

            for j in range(len(self.gf_struct_sumk[icrsh])):
                self.dc_imp[icrsh]['%s'%self.gf_struct_sumk[icrsh][j][0]] = numpy.identity(dim,numpy.float_)
                blname = self.gf_struct_sumk[icrsh][j][0]
                Ncr[blname] = 0.0

            for block,inner in self.gf_struct_solver[iorb].iteritems():
                bl = self.solver_to_sumk_block[iorb][block]
                Ncr[bl] += dens_mat[block].real.trace()

            Ncrtot = 0.0
            block_ind_list = [block for block,inner in self.gf_struct_sumk[icrsh]]
            for bl in block_ind_list:
                Ncrtot += Ncr[bl]

            # average the densities if there is no SP:
            if self.SP == 0:
                for bl in block_ind_list:
                    Ncr[bl] = Ncrtot / len(block_ind_list)
            # correction for SO: we have only one block in this case, but in DC we need N/2
            elif self.SP == 1 and self.SO == 1:
                for bl in block_ind_list:
                    Ncr[bl] = Ncrtot / 2.0

            if use_val is None:

                if use_dc_formula == 0: # FLL

                    self.dc_energ[icrsh] = U_interact / 2.0 * Ncrtot * (Ncrtot-1.0)
                    for bl in block_ind_list:
                        Uav = U_interact*(Ncrtot-0.5) - J_hund*(Ncr[bl] - 0.5)
                        self.dc_imp[icrsh][bl] *= Uav
                        self.dc_energ[icrsh]  -= J_hund / 2.0 * (Ncr[bl]) * (Ncr[bl]-1.0)
                        mpi.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())

                elif use_dc_formula == 1: # Held's formula, with U_interact the interorbital onsite interaction

                    self.dc_energ[icrsh] = (U_interact + (dim-1)*(U_interact-2.0*J_hund) + (dim-1)*(U_interact-3.0*J_hund))/(2*dim-1) / 2.0 * Ncrtot * (Ncrtot-1.0)
                    for bl in block_ind_list:
                        Uav =(U_interact + (dim-1)*(U_interact-2.0*J_hund) + (dim-1)*(U_interact-3.0*J_hund))/(2*dim-1) * (Ncrtot-0.5)
                        self.dc_imp[icrsh][bl] *= Uav
                        mpi.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())

                elif use_dc_formula == 2: # AMF

                    self.dc_energ[icrsh] = 0.5 * U_interact * Ncrtot * Ncrtot
                    for bl in block_ind_list:
                        Uav = U_interact*(Ncrtot - Ncr[bl]/dim) - J_hund * (Ncr[bl] - Ncr[bl]/dim)
                        self.dc_imp[icrsh][bl] *= Uav
                        self.dc_energ[icrsh] -= (U_interact + (dim-1)*J_hund)/dim * 0.5 * Ncr[bl] * Ncr[bl]
                        mpi.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())

                # output:
                mpi.report("DC energy for shell %s = %s"%(icrsh,self.dc_energ[icrsh]))

            else:

                block_ind_list = [block for block,inner in self.gf_struct_sumk[icrsh]]
                for bl in block_ind_list:
                    self.dc_imp[icrsh][bl] *= use_val

                self.dc_energ[icrsh] = use_val * Ncrtot

                # output:
                mpi.report("DC for shell %(icrsh)i = %(use_val)f"%locals())
                mpi.report("DC energy = %s"%self.dc_energ[icrsh])



    def put_Sigma(self, Sigma_imp):
        """Puts the impurity self energies for inequivalent atoms into the class, respects the multiplicity of the atoms."""

        assert isinstance(Sigma_imp,list), "Sigma_imp has to be a list of Sigmas for the correlated shells, even if it is of length 1!"
        assert len(Sigma_imp) == self.n_inequiv_shells, "give exactly one Sigma for each inequivalent corr. shell!"

        # init self.Sigma_imp:
        if all(type(gf) == GfImFreq for bname,gf in Sigma_imp[0]):
            # Imaginary frequency Sigma:
            self.Sigma_imp = [ BlockGf( name_block_generator = [ (block,GfImFreq(indices = inner, mesh = Sigma_imp[0].mesh)) for block,inner in self.gf_struct_sumk[icrsh] ],
                                  make_copies = False) for icrsh in range(self.n_corr_shells) ]
        elif all(type(gf) == GfReFreq for bname,gf in Sigma_imp[0]):
            # Real frequency Sigma:
            self.Sigma_imp = [ BlockGf( name_block_generator = [ (block,GfReFreq(indices = inner, mesh = Sigma_imp[0].mesh)) for block,inner in self.gf_struct_sumk[icrsh] ],
                                  make_copies = False) for icrsh in range(self.n_corr_shells) ]
        else:
            raise ValueError, "This type of Sigma is not handled."

        # transform the CTQMC blocks to the full matrix:
        for icrsh in range(self.n_corr_shells):
            ish = self.corr_to_inequiv[icrsh]    # ish is the index of the inequivalent shell corresponding to icrsh

            for block,inner in self.gf_struct_solver[ish].iteritems():
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk,ind1_sumk = self.solver_to_sumk[ish][(block,ind1)]
                        block_sumk,ind2_sumk = self.solver_to_sumk[ish][(block,ind2)]
                        self.Sigma_imp[icrsh][block_sumk][ind1_sumk,ind2_sumk] << Sigma_imp[ish][block][ind1,ind2]

        # rotation from local to global coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in self.Sigma_imp[icrsh]: self.Sigma_imp[icrsh][bname] << self.rotloc(icrsh,gf,direction='toGlobal')



    def add_dc(self):
        """Substracts the double counting term from the impurity self energy."""

        # Be careful: Sigma_imp is already in the global coordinate system!!
        sres = [s.copy() for s in self.Sigma_imp]
        for icrsh in range(self.n_corr_shells):
            for bname,gf in sres[icrsh]:
                # Transform dc_imp to global coordinate system
                dccont = numpy.dot(self.rot_mat[icrsh],numpy.dot(self.dc_imp[icrsh][bname],self.rot_mat[icrsh].conjugate().transpose()))
                sres[icrsh][bname] -= dccont

        return sres # list of self energies corrected by DC



    def set_mu(self,mu):
        """Sets a new chemical potential"""

        self.chemical_potential = mu


    def total_density(self, mu):
        """
        Calculates the total charge for the energy window for a given chemical potential mu. 
        Since in general n_orbitals depends on k, the calculation is done in the following order:
        G_aa'(k,iw) -> n(k) = Tr G_aa'(k,iw) -> sum_k n_k

        The calculation is done in the global coordinate system, if distinction is made between local/global!
        """

        dens = 0.0
        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            S = self.lattice_gf_matsubara(ik=ik,mu=mu)
            dens += self.bz_weights[ik] * S.total_density()

        # collect data from mpi:
        dens = mpi.all_reduce(mpi.world, dens, lambda x,y : x+y)
        mpi.barrier()

        return dens


    def find_mu(self, precision = 0.01):
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


    def extract_G_loc(self, mu=None, with_Sigma = True):
        """
        Extracts the local downfolded Green function at the chemical potential of the class.
        At the end, the local G is rotated from the global coordinate system to the local system.
        if with_Sigma = False: Sigma is not included => non-interacting local GF
        """

        if mu is None: mu = self.chemical_potential

        Gloc = [ self.Sigma_imp[icrsh].copy() for icrsh in range(self.n_corr_shells) ]   # this list will be returned
        for icrsh in range(self.n_corr_shells): Gloc[icrsh].zero()                       # initialize to zero
        beta = Gloc[0].mesh.beta

        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            S = self.lattice_gf_matsubara(ik=ik,mu=mu,with_Sigma = with_Sigma, beta = beta)
            S *= self.bz_weights[ik]

            for icrsh in range(self.n_corr_shells):
                tmp = Gloc[icrsh].copy()                  # init temporary storage
                for bname,gf in tmp: tmp[bname] << self.downfold(ik,icrsh,bname,S[bname],gf)
                Gloc[icrsh] += tmp

        #collect data from mpi:
        for icrsh in range(self.n_corr_shells):
            Gloc[icrsh] << mpi.all_reduce(mpi.world, Gloc[icrsh], lambda x,y : x+y)
        mpi.barrier()

        # Gloc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:
        if self.symm_op != 0: Gloc = self.symmcorr.symmetrize(Gloc)

        # Gloc is rotated to the local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                for bname,gf in Gloc[icrsh]: Gloc[icrsh][bname] << self.rotloc(icrsh,gf,direction='toLocal')

        # transform to CTQMC blocks:
        Glocret = [ BlockGf( name_block_generator = [ (block,GfImFreq(indices = inner, mesh = Gloc[0].mesh)) for block,inner in self.gf_struct_solver[ish].iteritems() ],
                        make_copies = False) for ish in range(self.n_inequiv_shells)  ]
        for ish in range(self.n_inequiv_shells):

            for block,inner in self.gf_struct_solver[ish].iteritems():
                for ind1 in inner:
                    for ind2 in inner:
                        block_sumk,ind1_sumk = self.solver_to_sumk[ish][(block,ind1)]
                        block_sumk,ind2_sumk = self.solver_to_sumk[ish][(block,ind2)]
                        Glocret[ish][block][ind1,ind2] << Gloc[self.inequiv_to_corr[ish]][block_sumk][ind1_sumk,ind2_sumk]

        # return only the inequivalent shells:
        return Glocret


    def calc_density_correction(self,filename = 'dens_mat.dat'):
        """ Calculates the density correction in order to feed it back to the DFT calculations."""

        assert type(filename) == StringType, "filename has to be a string!"

        ntoi = self.spin_names_to_ind[self.SO]
        bln = self.spin_block_names[self.SO]

        # Set up deltaN:
        deltaN = {}
        for b in bln:
            deltaN[b] = [ numpy.zeros( [self.n_orbitals[ik,ntoi[b]],self.n_orbitals[ik,ntoi[b]]], numpy.complex_) for ik in range(self.n_k)]

        ikarray=numpy.array(range(self.n_k))

        dens = {}
        for b in bln:
            dens[b] = 0.0

        for ik in mpi.slice_array(ikarray):
            S = self.lattice_gf_matsubara(ik=ik,mu=self.chemical_potential)
            for bname,gf in S:
                deltaN[bname][ik] = S[bname].density()
                dens[bname] += self.bz_weights[ik] * S[bname].total_density()

        #put mpi Barrier:
        for bname in deltaN:
            for ik in range(self.n_k):
                deltaN[bname][ik] = mpi.all_reduce(mpi.world, deltaN[bname][ik], lambda x,y : x+y)
            dens[bname] = mpi.all_reduce(mpi.world, dens[bname], lambda x,y : x+y)
        mpi.barrier()

        # now save to file:
        if mpi.is_master_node():
            if self.SP == 0:
                f=open(filename,'w')
            else:
                f=open(filename+'up','w')
                f1=open(filename+'dn','w')
            # write chemical potential (in Rydberg):
            f.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
            if self.SP != 0: f1.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
            # write beta in ryderg-1
            f.write("%.14f\n"%(S.mesh.beta*self.energy_unit))
            if self.SP != 0: f1.write("%.14f\n"%(S.mesh.beta*self.energy_unit))

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
                    isp, bn = to_write[fout]
                    for ik in range(self.n_k):
                        fout.write("%s\n"%self.n_orbitals[ik,isp])
                        for inu in range(self.n_orbitals[ik,isp]):
                            for imu in range(self.n_orbitals[ik,isp]):
                                fout.write("%.14f  %.14f "%(deltaN[bn][ik][inu,imu].real,deltaN[bn][ik][inu,imu].imag))
                            fout.write("\n")
                        fout.write("\n")
                    fout.close()

        return deltaN, dens


################
# FIXME LEAVE UNDOCUMENTED
################

    # FIXME Merge with find_mu?
    def find_mu_nonint(self, dens_req, orb = None, precision = 0.01):

        def F(mu):
            gnonint = self.extract_G_loc(mu=mu,with_Sigma=False)

            if orb is None:
                dens = 0.0
                for ish in range(self.n_inequiv_shells):
                    dens += gnonint[ish].total_density()
            else:
                dens = gnonint[orb].total_density()

            return dens


        self.chemical_potential = dichotomy.dichotomy(function = F,
                                         x_init = self.chemical_potential, y_value = dens_req,
                                         precision_on_y = precision, delta_x = 0.5, max_loops = 100, 
                                         x_name = "Chemical Potential", y_name = "Total Density",
                                         verbosity = 3)[0]

        return self.chemical_potential


    def find_dc(self,orb,guess,dens_mat,dens_req=None,precision=0.01):
        """Searches for DC in order to fulfill charge neutrality.
           If dens_req is given, then DC is set such that the LOCAL charge of orbital
           orb coincides with dens_req."""

        mu = self.chemical_potential

        def F(dc):
            self.set_dc(dens_mat=dens_mat,U_interact=0,J_hund=0,orb=orb,use_val=dc)
            if dens_req is None:
                return self.total_density(mu=mu)
            else:
                return self.extract_G_loc()[orb].total_density()


        if dens_req is None:
            density = self.density_required - self.charge_below
        else:
            density = dens_req

        dcnew = dichotomy.dichotomy(function = F,
                                    x_init = guess, y_value = density,
                                    precision_on_y = precision, delta_x=0.5,
                                    max_loops = 100, x_name="Double-Counting", y_name= "Total Density",
                                    verbosity = 3)[0]

        return dcnew


    # Check that the density matrix from projectors (DM = P Pdagger) is correct (ie matches DFT)
    def check_projectors(self):

        dens_mat = [numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]],numpy.complex_)
                   for icrsh in range(self.n_corr_shells)]

        for ik in range(self.n_k):
            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh][3]
                n_orb = self.n_orbitals[ik,0]
                projmat = self.proj_mat[ik,0,icrsh,0:dim,0:n_orb]
                dens_mat[icrsh][:,:] += numpy.dot(projmat, projmat.transpose().conjugate()) * self.bz_weights[ik]

        if self.symm_op != 0: dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if self.use_rotations:
            for icrsh in range(self.n_corr_shells):
                if self.rot_mat_time_inv[icrsh] == 1: dens_mat[icrsh] = dens_mat[icrsh].conjugate()
                dens_mat[icrsh] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh]) ,
                                            self.rot_mat[icrsh] )

        return dens_mat


    # Determine the number of equivalent shells
    def sorts_of_atoms(self,lst):
        """
        This routine should determine the number of sorts in the double list lst
        """
        sortlst = [ lst[i][1] for i in range(len(lst)) ]
        sorts = len(set(sortlst))
        return sorts


    # Determine the number of atoms
    def number_of_atoms(self,lst):
        """
        This routine should determine the number of atoms in the double list lst
        """
        atomlst = [ lst[i][0] for i in range(len(lst)) ]
        atoms = len(set(atomlst))
        return atoms
