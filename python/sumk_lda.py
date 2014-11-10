
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
from symmetry import *
import numpy
import pytriqs.utility.dichotomy as dichotomy
from pytriqs.gf.local import *
from pytriqs.archive import *
import pytriqs.utility.mpi as mpi

class SumkLDA:
    """This class provides a general SumK method for combining ab-initio code and pytriqs."""


    def __init__(self, hdf_file, mu = 0.0, h_field = 0.0, use_lda_blocks = False, 
	         lda_data = 'lda_input', symmcorr_data = 'lda_symmcorr_input', parproj_data = 'lda_parproj_input', 
                 symmpar_data = 'lda_symmpar_input', bands_data = 'lda_bands_input', lda_output = 'lda_output'):
        """
        Initialises the class from data previously stored into an HDF5
        """

        if not (type(hdf_file)==StringType):
            mpi.report("Give a string for the HDF5 filename to read the input!")
        else:
            self.hdf_file = hdf_file
            self.lda_data = lda_data
            self.symmcorr_data = symmcorr_data
            self.parproj_data = parproj_data
            self.symmpar_data = symmpar_data
            self.bands_data = bands_data
            self.lda_output = lda_output
            self.block_names = [ ['up','down'], ['ud'] ]
            self.n_spin_blocks_gf = [2,1]
            self.G_upfold = None
            self.h_field = h_field

            # read input from HDF:
            things_to_read = ['energy_unit','n_k','k_dep_projection','SP','SO','charge_below','density_required',
                              'symm_op','n_shells','shells','n_corr_shells','corr_shells','use_rotations','rot_mat',
                              'rot_mat_time_inv','n_reps','dim_reps','T','n_orbitals','proj_mat','bz_weights','hopping']
            self.subgroup_present, self.value_read = self.read_input_from_hdf(subgrp = self.lda_data, things_to_read = things_to_read)

            if (self.SO) and (abs(self.h_field)>0.000001):
                self.h_field=0.0
                mpi.report("For SO, the external magnetic field is not implemented, setting it to 0!")

            # determine the number of inequivalent correlated shells (self.n_inequiv_corr_shells)
            # and related maps (self.shellmap, self.invshellmap)
            self.inequiv_shells(self.corr_shells)

            # field to convert block_names to indices
            self.names_to_ind = [{}, {}]
            for ibl in range(2):
                for inm in range(self.n_spin_blocks_gf[ibl]):
                    self.names_to_ind[ibl][self.block_names[ibl][inm]] = inm * self.SP

            # GF structure used for the local things in the k sums
            # Most general form allowing for all hybridisation, i.e. largest blocks possible
            self.gf_struct_corr = [ [ (al, range( self.corr_shells[i][3])) for al in self.block_names[self.corr_shells[i][4]] ]
                                   for i in xrange(self.n_corr_shells) ]

            #-----
            # If these quantities are not in HDF, set them up and save into HDF
            optional_things = ['gf_struct_solver','map_inv','map','chemical_potential','dc_imp','dc_energ','deg_shells']
            self.subgroup_present, self.value_read = self.read_input_from_hdf(subgrp = self.lda_output, things_to_read = [], 
                                                                              optional_things = optional_things)
            if (not self.subgroup_present) or (not self.value_read['gf_struct_solver']):
                # No gf_struct was stored in HDF, so first set a standard one:
                self.gf_struct_solver = [ dict([ (al, range(self.corr_shells[self.invshellmap[i]][3]) )
                                                        for al in self.block_names[self.corr_shells[self.invshellmap[i]][4]] ])
                                          for i in range(self.n_inequiv_corr_shells)
                                        ]
                self.map = [ {} for i in xrange(self.n_inequiv_corr_shells) ]
                self.map_inv = [ {} for i in xrange(self.n_inequiv_corr_shells) ]
                for i in xrange(self.n_inequiv_corr_shells):
                    for al in self.block_names[self.corr_shells[self.invshellmap[i]][4]]:
                        self.map[i][al] = [al for j in range( self.corr_shells[self.invshellmap[i]][3] ) ]
                        self.map_inv[i][al] = al

            if (not self.subgroup_present) or (not self.value_read['dc_imp']):
                # init the double counting:
                self.__init_dc()

            if (not self.subgroup_present) or (not self.value_read['chemical_potential']):
                self.chemical_potential = mu

            if (not self.subgroup_present) or (not self.value_read['deg_shells']):
                self.deg_shells = [ [] for i in range(self.n_inequiv_corr_shells)]
            #-----

            if self.symm_op:
                #mpi.report("Do the init for symm:")
                self.symmcorr = Symmetry(hdf_file,subgroup=self.symmcorr_data)

            # Analyse the block structure and determine the smallest blocs, if desired
            if (use_lda_blocks): dm=self.analyse_BS()

            # Now save things again to HDF5: 
            # FIXME WHAT HAPPENS TO h_field? INPUT TO __INIT__? ADD TO OPTIONAL_THINGS?
            things_to_save=['chemical_potential','dc_imp','dc_energ','h_field']
            self.save(things_to_save)


    def read_input_from_hdf(self, subgrp, things_to_read=[], optional_things=[]):
        """
        Reads data from the HDF file
        """

        value_read = True
        # init variables on all nodes to ensure mpi broadcast works at the end
        for it in things_to_read: setattr(self,it,0)
        for it in optional_things: setattr(self,it,0)

        if mpi.is_master_node():
            ar=HDFArchive(self.hdf_file,'a')
            if subgrp in ar:
                subgroup_present = True
                # first read the necessary things:
                for it in things_to_read:
                    if it in ar[subgrp]:
                        setattr(self,it,ar[subgrp][it])
                    else:
                        mpi.report("Loading %s failed!"%it)
                        value_read = False

                if (value_read and (len(optional_things)>0)):
                    # if necessary things worked, now read optional things:
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

        return (subgroup_present, value_read)


    def save(self,things_to_save):
        """Saves some quantities into an HDF5 arxiv"""

        if not (mpi.is_master_node()): return # do nothing on nodes
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.lda_output in ar): ar.create_group(self.lda_output)
        for it in things_to_save: 
            try:
                ar[self.lda_output][it] = getattr(self,it)
            except:
                mpi.report("%s not found, and so not stored."%it)
        del ar


    def downfold(self,ik,icrsh,sig,gf_to_downfold,gf_inp):
        """Downfolding a block of the Greens function"""

        gf_downfolded = gf_inp.copy()
        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        dim = self.corr_shells[icrsh][3]
        n_orb = self.n_orbitals[ik,isp]

        gf_downfolded.from_L_G_R(self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb],gf_to_downfold,self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb].conjugate().transpose())  # downfolding G

        return gf_downfolded


    def upfold(self,ik,icrsh,sig,gf_to_upfold,gf_inp):
        """Upfolding a block of the Greens function"""

        gf_upfolded = gf_inp.copy()

        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        dim = self.corr_shells[icrsh][3]
        n_orb = self.n_orbitals[ik,isp]

        gf_upfolded.from_L_G_R(self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb].conjugate().transpose(),gf_to_upfold,self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb])

        return gf_upfolded


    def rotloc(self,icrsh,gf_to_rotate,direction):
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

        assert ((direction=='toLocal')or(direction=='toGlobal')),"Give direction 'toLocal' or 'toGlobal' in rotloc!"

        gf_rotated = gf_to_rotate.copy()
        if (direction=='toGlobal'):

            if ((self.rot_mat_time_inv[icrsh]==1) and (self.SO)):
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rot_mat[icrsh].conjugate(),gf_rotated,self.rot_mat[icrsh].transpose())
            else:
                gf_rotated.from_L_G_R(self.rot_mat[icrsh],gf_rotated,self.rot_mat[icrsh].conjugate().transpose())

        elif (direction=='toLocal'):

            if ((self.rot_mat_time_inv[icrsh]==1)and(self.SO)):
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rot_mat[icrsh].transpose(),gf_rotated,self.rot_mat[icrsh].conjugate())
            else:
                gf_rotated.from_L_G_R(self.rot_mat[icrsh].conjugate().transpose(),gf_rotated,self.rot_mat[icrsh])

        return gf_rotated


    def lattice_gf_matsubara(self,ik,mu,beta=40,with_Sigma=True):
        """Calculates the lattice Green function from the LDA hopping and the self energy at k-point number ik
           and chemical potential mu."""

        ntoi = self.names_to_ind[self.SO]
        bln = self.block_names[self.SO]

        if (not hasattr(self,"Sigma_imp")): with_Sigma=False

        if (with_Sigma):
            stmp = self.add_dc()
            beta = self.Sigma_imp[0].mesh.beta        # override beta if Sigma is present

        # Do we need to set up G_upfold?
        set_up_G_upfold = False # assume not
        if self.G_upfold is None: # yes if not G_upfold provided
            set_up_G_upfold = True
        else: # yes if inconsistencies present in existing G_upfold
            GFsize = [ gf.N1 for sig,gf in self.G_upfold]
            unchangedsize = all( [ self.n_orbitals[ik,ntoi[bln[ib]]]==GFsize[ib]
                                   for ib in range(self.n_spin_blocks_gf[self.SO]) ] )
            if ( (not unchangedsize) or (self.G_upfold.mesh.beta != beta) ): set_up_G_upfold = True

        # Set up G_upfold
        if set_up_G_upfold:
            BS = [ range(self.n_orbitals[ik,ntoi[ib]]) for ib in bln ]
            gf_struct = [ (bln[ib], BS[ib]) for ib in range(self.n_spin_blocks_gf[self.SO]) ]
            a_list = [a for a,al in gf_struct]
            if (with_Sigma):
                glist = lambda : [ GfImFreq(indices = al, mesh = self.Sigma_imp[0].mesh) for a,al in gf_struct]
            else:
                glist = lambda : [ GfImFreq(indices = al, beta = beta) for a,al in gf_struct]
            self.G_upfold = BlockGf(name_list = a_list, block_list = glist(),make_copies=False)
            self.G_upfold.zero()

        self.G_upfold << iOmega_n
        idmat = [numpy.identity(self.n_orbitals[ik,ntoi[bl]],numpy.complex_) for bl in bln]
        M = copy.deepcopy(idmat)
        for ibl in range(self.n_spin_blocks_gf[self.SO]):
            ind = ntoi[bln[ibl]]
            n_orb = self.n_orbitals[ik,ind]
            M[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb] - (idmat[ibl]*mu) - (idmat[ibl] * self.h_field * (1-2*ibl))
        self.G_upfold -= M

        if (with_Sigma):
            for icrsh in xrange(self.n_corr_shells):
                for sig,gf in self.G_upfold: gf -= self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)

        self.G_upfold.invert()

	return self.G_upfold




    def simple_point_dens_mat(self):


        ntoi = self.names_to_ind[self.SO]
        bln = self.block_names[self.SO]

        MMat = [numpy.zeros( [self.n_orbitals[0,ntoi[bl]],self.n_orbitals[0,ntoi[bl]]], numpy.complex_) for bl in bln]

        dens_mat = [ {} for icrsh in xrange(self.n_corr_shells)]
        for icrsh in xrange(self.n_corr_shells):
            for bl in self.block_names[self.corr_shells[icrsh][4]]:
                dens_mat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]], numpy.complex_)

        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            unchangedsize = all( [ self.n_orbitals[ik,ntoi[bln[ib]]]==len(MMat[ib])
                                   for ib in range(self.n_spin_blocks_gf[self.SO]) ] )

            if (not unchangedsize):
                MMat = [numpy.zeros( [self.n_orbitals[ik,ntoi[bl]],self.n_orbitals[ik,ntoi[bl]]], numpy.complex_) for bl in bln]

            for ibl,bl in enumerate(bln):
                ind = ntoi[bl]
                for inu in range(self.n_orbitals[ik,ind]):
                    if ( (self.hopping[ik,ind,inu,inu]-self.h_field*(1-2*ibl)) < 0.0): # ONLY WORKS FOR DIAGONAL HOPPING MATRIX (TRUE IN WIEN2K)
                        MMat[ibl][inu,inu] = 1.0
                    else:
                        MMat[ibl][inu,inu] = 0.0


            for icrsh in range(self.n_corr_shells):
                for ibn,bn in enumerate(self.block_names[self.corr_shells[icrsh][4]]):
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    dim = self.corr_shells[icrsh][3]
                    n_orb = self.n_orbitals[ik,isp]

                    #print ik, bn, isp
                    dens_mat[icrsh][bn] += self.bz_weights[ik] * numpy.dot( numpy.dot(self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb],MMat[ibn]) ,
                                                                           self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb].transpose().conjugate() )

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for sig in dens_mat[icrsh]:
                dens_mat[icrsh][sig] = mpi.all_reduce(mpi.world,dens_mat[icrsh][sig],lambda x,y : x+y)
        mpi.barrier()


        if (self.symm_op!=0): dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.n_corr_shells):
                for bn in dens_mat[icrsh]:
                    if (self.rot_mat_time_inv[icrsh]==1): dens_mat[icrsh][bn] = dens_mat[icrsh][bn].conjugate()
                    dens_mat[icrsh][bn] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh][bn]) ,
                                                    self.rot_mat[icrsh])


        return dens_mat

    # calculate upfolded gf, then density matrix -- no assumptions on structure (ie diagonal or not)
    def density_gf(self,beta):
        """Calculates the density without setting up Gloc. It is useful for Hubbard I, and very fast."""

        dens_mat = [ {} for icrsh in xrange(self.n_corr_shells)]
        for icrsh in xrange(self.n_corr_shells):
            for bl in self.block_names[self.corr_shells[icrsh][4]]:
                dens_mat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]], numpy.complex_)

        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            G_upfold = self.lattice_gf_matsubara(ik=ik, beta=beta, mu=self.chemical_potential)
            G_upfold *= self.bz_weights[ik]
            dm = G_upfold.density()
            MMat = [dm[bl] for bl in self.block_names[self.SO]]

            for icrsh in range(self.n_corr_shells):
                for ibn,bn in enumerate(self.block_names[self.corr_shells[icrsh][4]]):
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    dim = self.corr_shells[icrsh][3]
                    n_orb = self.n_orbitals[ik,isp]
                    #print ik, bn, isp
                    dens_mat[icrsh][bn] += numpy.dot( numpy.dot(self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb],MMat[ibn]),
                                                      self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb].transpose().conjugate() )

        # get data from nodes:
        for icrsh in range(self.n_corr_shells):
            for sig in dens_mat[icrsh]:
                dens_mat[icrsh][sig] = mpi.all_reduce(mpi.world,dens_mat[icrsh][sig],lambda x,y : x+y)
        mpi.barrier()


        if (self.symm_op!=0): dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.n_corr_shells):
                for bn in dens_mat[icrsh]:
                    if (self.rot_mat_time_inv[icrsh]==1): dens_mat[icrsh][bn] = dens_mat[icrsh][bn].conjugate()
                    dens_mat[icrsh][bn] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh][bn]) ,
                                                    self.rot_mat[icrsh] )

        return dens_mat



    def analyse_BS(self, threshold = 0.00001, include_shells = None, dm = None):
        """ Determines the Green function block structure from simple point integration."""

        if dm is None: dm = self.simple_point_dens_mat()

        dens_mat = [dm[self.invshellmap[ish]] for ish in xrange(self.n_inequiv_corr_shells) ]

        if include_shells is None: include_shells=range(self.n_inequiv_corr_shells)
        for ish in include_shells:

            self.gf_struct_solver[ish] = {} 
            gf_struct_temp = []

            a_list = [a for a,al in self.gf_struct_corr[self.invshellmap[ish]] ]
            for a in a_list:
                dm = dens_mat[ish][a]
                dmbool = (abs(dm) > threshold)          # gives an index list of entries larger that threshold

                offdiag = []
                for i in xrange(len(dmbool)):
                    for j in xrange(i,len(dmbool)):
                        if ((dmbool[i,j])&(i!=j)): offdiag.append([i,j])

                NBlocs = len(dmbool)
                blocs = [ [i] for i in range(NBlocs) ]

                for i in range(len(offdiag)):
                    if (offdiag[i][0]!=offdiag[i][1]):
                        for j in range(len(blocs[offdiag[i][1]])): blocs[offdiag[i][0]].append(blocs[offdiag[i][1]][j])
                        del blocs[offdiag[i][1]]
                        for j in range(i+1,len(offdiag)):
                            if (offdiag[j][0]==offdiag[i][1]): offdiag[j][0]=offdiag[i][0]
                            if (offdiag[j][1]==offdiag[i][1]): offdiag[j][1]=offdiag[i][0]
                            if (offdiag[j][0]>offdiag[i][1]): offdiag[j][0] -= 1
                            if (offdiag[j][1]>offdiag[i][1]): offdiag[j][1] -= 1
                            offdiag[j].sort()
                        NBlocs-=1

                for i in range(NBlocs):
                    blocs[i].sort()
                    self.gf_struct_solver[ish].update( [('%s_%s'%(a,i),range(len(blocs[i])))] )
                    gf_struct_temp.append( ('%s_%s'%(a,i),blocs[i]) )


                # map is the mapping of the blocs from the SK blocs to the CTQMC blocs:
                self.map[ish][a] = range(len(dmbool))
                for ibl in range(NBlocs):
                    for j in range(len(blocs[ibl])):
                        self.map[ish][a][blocs[ibl][j]] = '%s_%s'%(a,ibl)
                        self.map_inv[ish]['%s_%s'%(a,ibl)] = a


            # now calculate degeneracies of orbitals:
            dm = {}
            for bl in gf_struct_temp:
                bln = bl[0]
                ind = bl[1]
                # get dm for the blocks:
                dm[bln] = numpy.zeros([len(ind),len(ind)],numpy.complex_)
                for i in range(len(ind)):
                    for j in range(len(ind)):
                        dm[bln][i,j] = dens_mat[ish][self.map_inv[ish][bln]][ind[i],ind[j]]

            for bl in gf_struct_temp:
                for bl2 in gf_struct_temp:
                    if (dm[bl[0]].shape==dm[bl2[0]].shape) :
                        if ( ( (abs(dm[bl[0]]-dm[bl2[0]])<threshold).all() ) and (bl[0]!=bl2[0]) ):
                            # check if it was already there:
                            ind1=-1
                            ind2=-2
                            for n,ind in enumerate(self.deg_shells[ish]):
                                if (bl[0] in ind): ind1=n
                                if (bl2[0] in ind): ind2=n
                            if ((ind1<0)and(ind2>=0)):
                                self.deg_shells[ish][ind2].append(bl[0])
                            elif ((ind1>=0)and(ind2<0)):
                                self.deg_shells[ish][ind1].append(bl2[0])
                            elif ((ind1<0)and(ind2<0)):
                                self.deg_shells[ish].append([bl[0],bl2[0]])

        things_to_save=['gf_struct_solver','map','map_inv','deg_shells']
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

# for simply dft input, get crystal field splittings. 
    def eff_atomic_levels(self):
        """Calculates the effective atomic levels needed as input for the Hubbard I Solver."""

        # define matrices for inequivalent shells:
        eff_atlevels = [ {} for ish in range(self.n_inequiv_corr_shells) ]
        for ish in range(self.n_inequiv_corr_shells):
            for bn in self.block_names[self.corr_shells[self.invshellmap[ish]][4]]:
                eff_atlevels[ish][bn] = numpy.identity(self.corr_shells[self.invshellmap[ish]][3], numpy.complex_)

        # Chemical Potential:
        for ish in xrange(self.n_inequiv_corr_shells):
            for ii in eff_atlevels[ish]: eff_atlevels[ish][ii] *= -self.chemical_potential

        # double counting term:
        #if hasattr(self,"dc_imp"):
        for ish in xrange(self.n_inequiv_corr_shells):
            for ii in eff_atlevels[ish]:
                eff_atlevels[ish][ii] -= self.dc_imp[self.invshellmap[ish]][ii]

        # sum over k:
        if not hasattr(self,"Hsumk"):
            # calculate the sum over k. Does not depend on mu, so do it only once:
            self.Hsumk = [ {} for ish in range(self.n_corr_shells) ]
            for icrsh in range(self.n_corr_shells):
                for bn in self.block_names[self.corr_shells[icrsh][4]]:
                    dim = self.corr_shells[icrsh][3]  #*(1+self.corr_shells[icrsh][4])
                    self.Hsumk[icrsh][bn] = numpy.zeros([dim,dim],numpy.complex_)

            for icrsh in range(self.n_corr_shells):
                dim = self.corr_shells[icrsh][3]
                for ibn, bn in enumerate(self.block_names[self.corr_shells[icrsh][4]]):
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    for ik in xrange(self.n_k):
                        n_orb = self.n_orbitals[ik,isp]
                        MMat = numpy.identity(n_orb, numpy.complex_)
                        MMat = self.hopping[ik,isp,0:n_orb,0:n_orb] - (1-2*ibn) * self.h_field * MMat
                        self.Hsumk[icrsh][bn] += self.bz_weights[ik] * numpy.dot( numpy.dot(self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb],MMat),
                                                                                  self.proj_mat[ik,isp,icrsh,0:dim,0:n_orb].conjugate().transpose() )

            # symmetrisation:
            if (self.symm_op!=0): self.Hsumk = self.symmcorr.symmetrize(self.Hsumk)

            # Rotate to local coordinate system:
            if (self.use_rotations):
                for icrsh in xrange(self.n_corr_shells):
                    for bn in self.Hsumk[icrsh]:

                        if (self.rot_mat_time_inv[icrsh]==1): self.Hsumk[icrsh][bn] = self.Hsumk[icrsh][bn].conjugate()
                        #if (self.corr_shells[icrsh][4]==0): self.Hsumk[icrsh][bn] = self.Hsumk[icrsh][bn].conjugate()

                        self.Hsumk[icrsh][bn] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),self.Hsumk[icrsh][bn]) ,
                                                           self.rot_mat[icrsh] )

        # add to matrix:
        for ish in xrange(self.n_inequiv_corr_shells):
            for bn in eff_atlevels[ish]:
                eff_atlevels[ish][bn] += self.Hsumk[self.invshellmap[ish]][bn]


        return eff_atlevels



    def __init_dc(self):

        # construct the density matrix dm_imp and double counting arrays
        #self.dm_imp = [ {} for i in xrange(self.n_corr_shells)]
        self.dc_imp = [ {} for i in xrange(self.n_corr_shells)]
        for i in xrange(self.n_corr_shells):
            l = self.corr_shells[i][3]
            for j in xrange(len(self.gf_struct_corr[i])):
                self.dc_imp[i]['%s'%self.gf_struct_corr[i][j][0]] = numpy.zeros([l,l],numpy.float_)
        self.dc_energ = [0.0 for i in xrange(self.n_corr_shells)]



    def set_dc(self,dens_mat,U_interact,J_hund,orb=0,use_dc_formula=0,use_val=None):
        """Sets the double counting term for inequiv orbital orb:
           use_dc_formula=0: LDA+U FLL double counting,
           use_dc_formula=1: Held's formula,
           use_dc_formula=2: AMF.
           Be sure that you are using the correct interaction Hamiltonian!"""


        for icrsh in xrange(self.n_corr_shells):

            iorb = self.shellmap[icrsh]    # iorb is the index of the inequivalent shell corresponding to icrsh

            if (iorb==orb):
                # do this orbital
                Ncr = {}
                l = self.corr_shells[icrsh][3] #*(1+self.corr_shells[icrsh][4])

                for j in xrange(len(self.gf_struct_corr[icrsh])):
                    self.dc_imp[icrsh]['%s'%self.gf_struct_corr[icrsh][j][0]] = numpy.identity(l,numpy.float_)
                    blname = self.gf_struct_corr[icrsh][j][0]
                    Ncr[blname] = 0.0

                for a,al in self.gf_struct_solver[iorb].iteritems():
                    bl = self.map_inv[iorb][a]
                    Ncr[bl] += dens_mat[a].real.trace()

                M = self.corr_shells[icrsh][3]

                Ncrtot = 0.0
                a_list = [a for a,al in self.gf_struct_corr[icrsh]]
                for bl in a_list:
                    Ncrtot += Ncr[bl]

                # average the densities if there is no SP:
                if (self.SP==0):
                    for bl in a_list:
                        Ncr[bl] = Ncrtot / len(a_list)
                # correction for SO: we have only one block in this case, but in DC we need N/2
                elif (self.SP==1 and self.SO==1):
                    for bl in a_list:
                        Ncr[bl] = Ncrtot / 2.0

                if (use_val is None):

                    if (use_dc_formula==0): # FLL

                        self.dc_energ[icrsh] = U_interact / 2.0 * Ncrtot * (Ncrtot-1.0)
                        for bl in a_list:
                            Uav = U_interact*(Ncrtot-0.5) - J_hund*(Ncr[bl] - 0.5)
                            self.dc_imp[icrsh][bl] *= Uav
                            self.dc_energ[icrsh]  -= J_hund / 2.0 * (Ncr[bl]) * (Ncr[bl]-1.0)
                            mpi.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())

                    elif (use_dc_formula==1): # Held's formula, with U_interact the interorbital onsite interaction

                        self.dc_energ[icrsh] = (U_interact + (M-1)*(U_interact-2.0*J_hund) + (M-1)*(U_interact-3.0*J_hund))/(2*M-1) / 2.0 * Ncrtot * (Ncrtot-1.0)
                        for bl in a_list:
                            Uav =(U_interact + (M-1)*(U_interact-2.0*J_hund) + (M-1)*(U_interact-3.0*J_hund))/(2*M-1) * (Ncrtot-0.5)
                            self.dc_imp[icrsh][bl] *= Uav
                            mpi.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())

                    elif (use_dc_formula==2): # AMF

                        self.dc_energ[icrsh] = 0.5 * U_interact * Ncrtot * Ncrtot
                        for bl in a_list:
                            Uav = U_interact*(Ncrtot - Ncr[bl]/M) - J_hund * (Ncr[bl] - Ncr[bl]/M)
                            self.dc_imp[icrsh][bl] *= Uav
                            self.dc_energ[icrsh] -= (U_interact + (M-1)*J_hund)/M * 0.5 * Ncr[bl] * Ncr[bl]
                            mpi.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())

                    # output:
                    mpi.report("DC energy for shell %s = %s"%(icrsh,self.dc_energ[icrsh]))

                else:

                    a_list = [a for a,al in self.gf_struct_corr[icrsh]]
                    for bl in a_list:
                        self.dc_imp[icrsh][bl] *= use_val

                    self.dc_energ[icrsh] = use_val * Ncrtot

                    # output:
                    mpi.report("DC for shell %(icrsh)i = %(use_val)f"%locals())
                    mpi.report("DC energy = %s"%self.dc_energ[icrsh])



    def put_Sigma(self, Sigma_imp):
        """Puts the impurity self energies for inequivalent atoms into the class, respects the multiplicity of the atoms."""

        assert isinstance(Sigma_imp,list), "Sigma_imp has to be a list of Sigmas for the correlated shells, even if it is of length 1!"
        assert len(Sigma_imp)==self.n_inequiv_corr_shells, "give exactly one Sigma for each inequivalent corr. shell!"

        # init self.Sigma_imp:
        if all(type(g) == GfImFreq for name,g in Sigma_imp[0]):
            # Imaginary frequency Sigma:
            self.Sigma_imp = [ BlockGf( name_block_generator = [ (a,GfImFreq(indices = al, mesh = Sigma_imp[0].mesh)) for a,al in self.gf_struct_corr[i] ],
                                  make_copies = False) for i in xrange(self.n_corr_shells) ]
        elif all(type(g) == GfReFreq for name,g in Sigma_imp[0]):
            # Real frequency Sigma:
            self.Sigma_imp = [ BlockGf( name_block_generator = [ (a,GfReFreq(indices = al, mesh = Sigma_imp[0].mesh)) for a,al in self.gf_struct_corr[i] ],
                                  make_copies = False) for i in xrange(self.n_corr_shells) ]
        else:
            raise ValueError, "This type of Sigma is not handled."

        # transform the CTQMC blocks to the full matrix:
        for icrsh in xrange(self.n_corr_shells):
            s = self.shellmap[icrsh]    # s is the index of the inequivalent shell corresponding to icrsh

            # setting up the index map:
            map_ind={}
            cnt = {}
            for blname in self.map[s]:
                cnt[blname] = 0

            for a,al in self.gf_struct_solver[s].iteritems():
                blname = self.map_inv[s][a]
                map_ind[a] = range(len(al))
                for i in al:
                    map_ind[a][i] = cnt[blname]
                    cnt[blname]+=1

            for bl, orblist in self.gf_struct_solver[s].iteritems():
                for i in range(len(orblist)):
                    for j in range(len(orblist)):
                        ind1 = orblist[i]	
                        ind2 = orblist[j]	
                        ind1_imp = map_ind[bl][ind1]
                        ind2_imp = map_ind[bl][ind2]
                        self.Sigma_imp[icrsh][self.map_inv[s][bl]][ind1_imp,ind2_imp] << Sigma_imp[s][bl][ind1,ind2]

        # rotation from local to global coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.n_corr_shells):
                for sig,gf in self.Sigma_imp[icrsh]: self.Sigma_imp[icrsh][sig] << self.rotloc(icrsh,gf,direction='toGlobal')



    def add_dc(self):
        """Substracts the double counting term from the impurity self energy."""

        # Be careful: Sigma_imp is already in the global coordinate system!!
        sres = [s.copy() for s in self.Sigma_imp]
        for icrsh in xrange(self.n_corr_shells):
            for bl,gf in sres[icrsh]:
                # Transform dc_imp to global coordinate system
                dccont = numpy.dot(self.rot_mat[icrsh],numpy.dot(self.dc_imp[icrsh][bl],self.rot_mat[icrsh].conjugate().transpose()))
                sres[icrsh][bl] -= dccont

        return sres # list of self energies corrected by DC



    def set_mu(self,mu):
        """Sets a new chemical potential"""

        self.chemical_potential = mu


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



    def total_density(self, mu):
        """
        Calculates the total charge for the energy window for a given mu. Since in general n_orbitals depends on k,
        the calculation is done in the following order:
        G_aa'(k,iw) -> n(k) = Tr G_aa'(k,iw) -> sum_k n_k

        mu: chemical potential

        The calculation is done in the global coordinate system, if distinction is made between local/global!
        """

        dens = 0.0
        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            S = self.lattice_gf_matsubara(ik=ik,mu=mu)
            dens += self.bz_weights[ik] * S.total_density()

        # collect data from mpi:
        dens = mpi.all_reduce(mpi.world,dens,lambda x,y : x+y)
        mpi.barrier()

        return dens


    def find_mu(self, precision = 0.01):
        """Searches for mu in order to give the desired charge
        A desired precision can be specified in precision."""

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
        extracts the local downfolded Green function at the chemical potential of the class.
        At the end, the local G is rotated from the global coordinate system to the local system.
        if with_Sigma = False: Sigma is not included => non-interacting local GF
        """

        if (mu is None): mu = self.chemical_potential

        Gloc = [ self.Sigma_imp[icrsh].copy() for icrsh in xrange(self.n_corr_shells) ]   # this list will be returned
        for icrsh in xrange(self.n_corr_shells): Gloc[icrsh].zero()                # initialize to zero
        beta = Gloc[0].mesh.beta

        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            S = self.lattice_gf_matsubara(ik=ik,mu=mu,with_Sigma = with_Sigma, beta = beta)
            S *= self.bz_weights[ik]


            for icrsh in xrange(self.n_corr_shells):
                tmp = Gloc[icrsh].copy()                  # init temporary storage
                for sig,gf in tmp: tmp[sig] << self.downfold(ik,icrsh,sig,S[sig],gf)
                Gloc[icrsh] += tmp

        #collect data from mpi:
        for icrsh in xrange(self.n_corr_shells):
            Gloc[icrsh] << mpi.all_reduce(mpi.world,Gloc[icrsh],lambda x,y : x+y)
        mpi.barrier()


        # Gloc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:
        if (self.symm_op!=0): Gloc = self.symmcorr.symmetrize(Gloc)

        # Gloc is rotated to the local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.n_corr_shells):
                for sig,gf in Gloc[icrsh]: Gloc[icrsh][sig] << self.rotloc(icrsh,gf,direction='toLocal')

        # transform to CTQMC blocks:
        Glocret = [ BlockGf( name_block_generator = [ (a,GfImFreq(indices = al, mesh = Gloc[0].mesh)) for a,al in self.gf_struct_solver[i].iteritems() ],
                        make_copies = False) for i in xrange(self.n_inequiv_corr_shells)  ]
        for ish in xrange(self.n_inequiv_corr_shells):

            # setting up the index map:
            map_ind={}
            cnt = {}
            for blname in self.map[ish]:
                cnt[blname] = 0

            for a,al in self.gf_struct_solver[ish].iteritems():
                blname = self.map_inv[ish][a]
                map_ind[a] = range(len(al))
                for i in al:
                    map_ind[a][i] = cnt[blname]
                    cnt[blname]+=1

            for bl, orblist in self.gf_struct_solver[ish].iteritems():
                for i in range(len(orblist)):
                    for j in range(len(orblist)):
                        ind1 = orblist[i]	
                        ind2 = orblist[j]	
                        ind1_imp = map_ind[bl][ind1]
                        ind2_imp = map_ind[bl][ind2]
                        Glocret[ish][bl][ind1,ind2] << Gloc[self.invshellmap[ish]][self.map_inv[ish][bl]][ind1_imp,ind2_imp]


        # return only the inequivalent shells:
        return Glocret


    def calc_density_correction(self,filename = 'dens_mat.dat'):
        """ Calculates the density correction in order to feed it back to the DFT calculations."""


        assert (type(filename)==StringType), "filename has to be a string!"

        ntoi = self.names_to_ind[self.SO]
        bln = self.block_names[self.SO]

        # Set up deltaN:
        deltaN = {}
        for ib in bln:
            deltaN[ib] = [ numpy.zeros( [self.n_orbitals[ik,ntoi[ib]],self.n_orbitals[ik,ntoi[ib]]], numpy.complex_) for ik in range(self.n_k)]

        ikarray=numpy.array(range(self.n_k))

        dens = {}
        for ib in bln:
            dens[ib] = 0.0

        for ik in mpi.slice_array(ikarray):

            S = self.lattice_gf_matsubara(ik=ik,mu=self.chemical_potential)
            for sig,g in S:
                deltaN[sig][ik] = S[sig].density()
                dens[sig] += self.bz_weights[ik] * S[sig].total_density()



        #put mpi Barrier:
        for sig in deltaN:
            for ik in range(self.n_k):
                deltaN[sig][ik] = mpi.all_reduce(mpi.world,deltaN[sig][ik],lambda x,y : x+y)
            dens[sig] = mpi.all_reduce(mpi.world,dens[sig],lambda x,y : x+y)
        mpi.barrier()


        # now save to file:
        if (mpi.is_master_node()):
            if (self.SP==0):
                f=open(filename,'w')
            else:
                f=open(filename+'up','w')
                f1=open(filename+'dn','w')
            # write chemical potential (in Rydberg):
            f.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
            if (self.SP!=0): f1.write("%.14f\n"%(self.chemical_potential/self.energy_unit))
            # write beta in ryderg-1
            f.write("%.14f\n"%(S.mesh.beta*self.energy_unit))
            if (self.SP!=0): f1.write("%.14f\n"%(S.mesh.beta*self.energy_unit))
            if (self.SP==0):
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
            elif ((self.SP==1)and(self.SO==0)):
                for ik in range(self.n_k):
                    f.write("%s\n"%self.n_orbitals[ik,0])
                    for inu in range(self.n_orbitals[ik,0]):
                        for imu in range(self.n_orbitals[ik,0]):
                            f.write("%.14f  %.14f "%(deltaN['up'][ik][inu,imu].real,deltaN['up'][ik][inu,imu].imag))
                        f.write("\n")
                    f.write("\n")
                f.close()
                for ik in range(self.n_k):
                    f1.write("%s\n"%self.n_orbitals[ik,1])
                    for inu in range(self.n_orbitals[ik,1]):
                        for imu in range(self.n_orbitals[ik,1]):
                            f1.write("%.14f  %.14f "%(deltaN['down'][ik][inu,imu].real,deltaN['down'][ik][inu,imu].imag))
                        f1.write("\n")
                    f1.write("\n")
                f1.close()
            else:
                for ik in range(self.n_k):
                    f.write("%s\n"%self.n_orbitals[ik,0])
                    for inu in range(self.n_orbitals[ik,0]):
                        for imu in range(self.n_orbitals[ik,0]):
                            f.write("%.14f  %.14f "%(deltaN['ud'][ik][inu,imu].real,deltaN['ud'][ik][inu,imu].imag))
                        f.write("\n")
                    f.write("\n")
                f.close()
                for ik in range(self.n_k):
                    f1.write("%s\n"%self.n_orbitals[ik,0])
                    for inu in range(self.n_orbitals[ik,0]):
                        for imu in range(self.n_orbitals[ik,0]):
                            f1.write("%.14f  %.14f "%(deltaN['ud'][ik][inu,imu].real,deltaN['ud'][ik][inu,imu].imag))
                        f1.write("\n")
                    f1.write("\n")
                f1.close()


        return deltaN, dens



################
# FIXME LEAVE UNDOCUMENTED
################

    def find_mu_nonint(self, dens_req, orb = None, beta = 40, precision = 0.01):

        def F(mu):
            gnonint = self.extract_G_loc(mu=mu,with_Sigma=False)

            if (orb is None):
                dens = 0.0
                for ish in range(self.n_inequiv_corr_shells):
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
            if (dens_req is None):
                return self.total_density(mu=mu)
            else:
                return self.extract_G_loc()[orb].total_density()


        if (dens_req is None):
            density = self.density_required - self.charge_below
        else:
            density = dens_req

        dcnew = dichotomy.dichotomy(function = F,
                                    x_init = guess, y_value = density,
                                    precision_on_y = precision, delta_x=0.5,
                                    max_loops = 100, x_name="Double-Counting", y_name= "Total Density",
                                    verbosity = 3)[0]

        return dcnew


# FIXME Check that dens matrix from projectors (DM=PPdagger) is correct (ie matches DFT)
    def check_projectors(self):

        dens_mat = [numpy.zeros([self.corr_shells[ish][3],self.corr_shells[ish][3]],numpy.complex_)
                   for ish in range(self.n_corr_shells)]

        for ik in range(self.n_k):

            for ish in range(self.n_corr_shells):
                dim = self.corr_shells[ish][3]
                n_orb = self.n_orbitals[ik,0]
                dens_mat[ish][:,:] += numpy.dot(self.proj_mat[ik,0,ish,0:dim,0:n_orb],self.proj_mat[ik,0,ish,0:dim,0:n_orb].transpose().conjugate()) * self.bz_weights[ik]

        if (self.symm_op!=0): dens_mat = self.symmcorr.symmetrize(dens_mat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.n_corr_shells):
                if (self.rot_mat_time_inv[icrsh]==1): dens_mat[icrsh] = dens_mat[icrsh].conjugate()
                dens_mat[icrsh] = numpy.dot( numpy.dot(self.rot_mat[icrsh].conjugate().transpose(),dens_mat[icrsh]) ,
                                            self.rot_mat[icrsh] )


        return dens_mat


# FIXME DETERMINE EQUIVALENCY OF SHELLS
    def sorts_of_atoms(self,lst):
        """
        This routine should determine the number of sorts in the double list lst
        """
        sortlst = [ lst[i][1] for i in xrange(len(lst)) ]
        sortlst.sort()
        sorts = 1
        for i in xrange(len(sortlst)-1):
            if sortlst[i+1]>sortlst[i]: sorts += 1

        return sorts


    def number_of_atoms(self,lst):
        """
        This routine should determine the number of atoms in the double list lst
        """
        atomlst = [ lst[i][0] for i in xrange(len(lst)) ]
        atomlst.sort()
        atoms = 1
        for i in xrange(len(atomlst)-1):
            if atomlst[i+1]>atomlst[i]: atoms += 1

        return atoms
