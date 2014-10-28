
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
from pytriqs.operators import *
import pytriqs.utility.mpi as mpi
from datetime import datetime
from symmetry import *
from sumk_lda import SumkLDA
import string

def read_fortran_file (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exist."%filename
    for line in open(filename,'r') :
        for x in line.replace('D','E').split() :
            yield string.atof(x)


class SumkLDATools(SumkLDA):
    """Extends the SumkLDA class with some tools for analysing the data."""


    def __init__(self, hdf_file, mu = 0.0, h_field = 0.0, use_lda_blocks = False, lda_data = 'SumK_LDA', symm_corr_data = 'SymmCorr',
                 par_proj_data = 'SumK_LDA_ParProj', symm_par_data = 'SymmPar', bands_data = 'SumK_LDA_Bands'):

        self.G_upfold_refreq = None
        SumkLDA.__init__(self,hdf_file=hdf_file,mu=mu,h_field=h_field,use_lda_blocks=use_lda_blocks,lda_data=lda_data,
                          symm_corr_data=symm_corr_data,par_proj_data=par_proj_data,symm_par_data=symm_par_data,
                          bands_data=bands_data)


    def downfold_pc(self,ik,ir,ish,sig,gf_to_downfold,gf_inp):
        """Downfolding a block of the Greens function"""

        gf_downfolded = gf_inp.copy()
        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        dim = self.shells[ish][3]
        n_orb = self.n_orbitals[ik,isp]
        L=self.proj_mat_pc[ik,isp,ish,ir,0:dim,0:n_orb]
        R=self.proj_mat_pc[ik,isp,ish,ir,0:dim,0:n_orb].conjugate().transpose()
        gf_downfolded.from_L_G_R(L,gf_to_downfold,R)

        return gf_downfolded


    def rotloc_all(self,ish,gf_to_rotate,direction):
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

        assert ((direction=='toLocal')or(direction=='toGlobal')),"Give direction 'toLocal' or 'toGlobal' in rotloc!"


        gf_rotated = gf_to_rotate.copy()
        if (direction=='toGlobal'):
            if ((self.rot_mat_all_time_inv[ish]==1) and (self.SO)):
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rot_mat_all[ish].conjugate(),gf_rotated,self.rot_mat_all[ish].transpose())
            else:
                gf_rotated.from_L_G_R(self.rot_mat_all[ish],gf_rotated,self.rot_mat_all[ish].conjugate().transpose())

        elif (direction=='toLocal'):
            if ((self.rot_mat_all_time_inv[ish]==1)and(self.SO)):
                gf_rotated << gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rot_mat_all[ish].transpose(),gf_rotated,self.rot_mat_all[ish].conjugate())
            else:
                gf_rotated.from_L_G_R(self.rot_mat_all[ish].conjugate().transpose(),gf_rotated,self.rot_mat_all[ish])


        return gf_rotated


    def lattice_gf_realfreq(self, ik, mu, broadening, mesh=None, with_Sigma=True):
        """Calculates the lattice Green function on the real frequency axis. If self energy is
           present and with_Sigma=True, the mesh is taken from Sigma. Otherwise, the mesh has to be given."""

        ntoi = self.names_to_ind[self.SO]
        bln = self.block_names[self.SO]

        if (not hasattr(self,"Sigma_imp")): with_Sigma=False
        if (with_Sigma):
            assert all(type(g) == GfReFreq for name,g in self.Sigma_imp[0]), "Real frequency Sigma needed for lattice_gf_realfreq!"
            stmp = self.add_dc()
        else:
            assert (not (mesh is None)),"Without Sigma, give the mesh=(om_min,om_max,n_points) for lattice_gf_realfreq!"

        if (self.G_upfold_refreq is None):
            # first setting up of G_upfold_refreq
            BS = [ range(self.n_orbitals[ik,ntoi[ib]]) for ib in bln ]
            gf_struct = [ (bln[ib], BS[ib]) for ib in range(self.n_spin_blocks_gf[self.SO]) ]
            a_list = [a for a,al in gf_struct]
            if (with_Sigma):
                glist = lambda : [ GfReFreq(indices = al, mesh =self.Sigma_imp[0].mesh) for a,al in gf_struct]
            else:
                glist = lambda : [ GfReFreq(indices = al, window=(mesh[0],mesh[1]),n_points=mesh[2]) for a,al in gf_struct]
            self.G_upfold_refreq = BlockGf(name_list = a_list, block_list = glist(),make_copies=False)
            self.G_upfold_refreq.zero()

        GFsize = [ gf.N1 for sig,gf in self.G_upfold_refreq]
        unchangedsize = all( [ self.n_orbitals[ik,ntoi[bln[ib]]]==GFsize[ib]
                               for ib in range(self.n_spin_blocks_gf[self.SO]) ] )

        if (not unchangedsize):
            BS = [ range(self.n_orbitals[ik,ntoi[ib]]) for ib in bln ]
            gf_struct = [ (bln[ib], BS[ib]) for ib in range(self.n_spin_blocks_gf[self.SO]) ]
            a_list = [a for a,al in gf_struct]
            if (with_Sigma):
                glist = lambda : [ GfReFreq(indices = al, mesh =self.Sigma_imp[0].mesh) for a,al in gf_struct]
            else:
                glist = lambda : [ GfReFreq(indices = al, window=(mesh[0],mesh[1]),n_points=mesh[2]) for a,al in gf_struct]
            self.G_upfold_refreq = BlockGf(name_list = a_list, block_list = glist(),make_copies=False)
            self.G_upfold_refreq.zero()

        idmat = [numpy.identity(self.n_orbitals[ik,ntoi[bl]],numpy.complex_) for bl in bln]

        self.G_upfold_refreq << Omega + 1j*broadening
        M = copy.deepcopy(idmat)
        for ibl in range(self.n_spin_blocks_gf[self.SO]):
            ind = ntoi[bln[ibl]]
            n_orb = self.n_orbitals[ik,ind]
            M[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb] - (idmat[ibl]*mu) - (idmat[ibl] * self.h_field * (1-2*ibl))
        self.G_upfold_refreq -= M

        if (with_Sigma):
            tmp = self.G_upfold_refreq.copy()    # init temporary storage
            for icrsh in xrange(self.n_corr_shells):
                for sig,gf in tmp: tmp[sig] << self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)
                self.G_upfold_refreq -= tmp      # adding to the upfolded GF

        self.G_upfold_refreq.invert()

        return self.G_upfold_refreq



    def check_input_dos(self, om_min, om_max, n_om, beta=10, broadening=0.01):


        delta_om = (om_max-om_min)/(n_om-1)
        om_mesh = numpy.zeros([n_om],numpy.float_)
        for i in range(n_om): om_mesh[i] = om_min + delta_om * i

        DOS = {}
        for bn in self.block_names[self.SO]:
            DOS[bn] = numpy.zeros([n_om],numpy.float_)

        DOSproj     = [ {} for icrsh in range(self.n_inequiv_corr_shells) ]
        DOSproj_orb = [ {} for icrsh in range(self.n_inequiv_corr_shells) ]
        for icrsh in range(self.n_inequiv_corr_shells):
            for bn in self.block_names[self.corr_shells[self.invshellmap[icrsh]][4]]:
                dl = self.corr_shells[self.invshellmap[icrsh]][3]
                DOSproj[icrsh][bn] = numpy.zeros([n_om],numpy.float_)
                DOSproj_orb[icrsh][bn] = numpy.zeros([n_om,dl,dl],numpy.float_)

        # init:
        Gloc = []
        for icrsh in range(self.n_corr_shells):
            b_list = [a for a,al in self.gf_struct_corr[icrsh]]
            #glist = lambda : [ GfReFreq(indices = al, beta = beta, mesh_array = mesh) for a,al in self.gf_struct_corr[icrsh]]
            glist = lambda : [ GfReFreq(indices = al, window = (om_min,om_max), n_points = n_om) for a,al in self.gf_struct_corr[icrsh]]
            Gloc.append(BlockGf(name_list = b_list, block_list = glist(),make_copies=False))
        for icrsh in xrange(self.n_corr_shells): Gloc[icrsh].zero()                        # initialize to zero

        for ik in xrange(self.n_k):

            G_upfold=self.lattice_gf_realfreq(ik=ik,mu=self.chemical_potential,broadening=broadening,mesh=(om_min,om_max,n_om),with_Sigma=False)
            G_upfold *= self.bz_weights[ik]

            # non-projected DOS
            for iom in range(n_om):
                for sig,gf in G_upfold:
                    asd = gf.data[iom,:,:].imag.trace()/(-3.1415926535)
                    DOS[sig][iom] += asd

            for icrsh in xrange(self.n_corr_shells):
                tmp = Gloc[icrsh].copy()
                for sig,gf in tmp: tmp[sig] << self.downfold(ik,icrsh,sig,G_upfold[sig],gf) # downfolding G
                Gloc[icrsh] += tmp



        if (self.symm_op!=0): Gloc = self.Symm_corr.symmetrize(Gloc)

        if (self.use_rotations):
            for icrsh in xrange(self.n_corr_shells):
                for sig,gf in Gloc[icrsh]: Gloc[icrsh][sig] << self.rotloc(icrsh,gf,direction='toLocal')

        # Gloc can now also be used to look at orbitally resolved quantities
        for ish in range(self.n_inequiv_corr_shells):
            for sig,gf in Gloc[self.invshellmap[ish]]: # loop over spins
                for iom in range(n_om): DOSproj[ish][sig][iom] += gf.data[iom,:,:].imag.trace()/(-3.1415926535)

                DOSproj_orb[ish][sig][:,:,:] += gf.data[:,:,:].imag/(-3.1415926535)

        # output:
        if (mpi.is_master_node()):
            for bn in self.block_names[self.SO]:
                f=open('DOS%s.dat'%bn, 'w')
                for i in range(n_om): f.write("%s    %s\n"%(om_mesh[i],DOS[bn][i]))
                f.close()

                for ish in range(self.n_inequiv_corr_shells):
                    f=open('DOS%s_proj%s.dat'%(bn,ish),'w')
                    for i in range(n_om): f.write("%s    %s\n"%(om_mesh[i],DOSproj[ish][bn][i]))
                    f.close()

                    for i in range(self.corr_shells[self.invshellmap[ish]][3]):
                        for j in range(i,self.corr_shells[self.invshellmap[ish]][3]):
                            Fname = 'DOS'+bn+'_proj'+str(ish)+'_'+str(i)+'_'+str(j)+'.dat'
                            f=open(Fname,'w')
                            for iom in range(n_om): f.write("%s    %s\n"%(om_mesh[iom],DOSproj_orb[ish][bn][iom,i,j]))
                            f.close()




    def read_par_proj_input_from_hdf(self):
        """
        Reads the data for the partial projectors from the HDF file
        """

        thingstoread = ['dens_mat_below','n_parproj','proj_mat_pc','rot_mat_all','rot_mat_all_time_inv']
        read_value = self.read_input_from_hdf(subgrp=self.par_proj_data,things_to_read = thingstoread)
        return read_value



    def dos_partial(self,broadening=0.01):
        """calculates the orbitally-resolved DOS"""

        assert hasattr(self,"Sigma_imp"), "Set Sigma First!!"

        #thingstoread = ['Dens_Mat_below','N_parproj','Proj_Mat_pc','rotmat_all']
        #read_value = self.read_input_from_HDF(SubGrp=self.par_proj_data, thingstoread=thingstoread)
        read_value = self.read_par_proj_input_from_hdf()
        if not read_value: return read_value
        if self.symm_op: self.Symm_par = Symmetry(self.hdf_file,subgroup=self.symm_par_data)

        mu = self.chemical_potential

        gf_struct_proj = [ [ (al, range(self.shells[i][3])) for al in self.block_names[self.SO] ]  for i in xrange(self.n_shells) ]
        Gproj = [BlockGf(name_block_generator = [ (a,GfReFreq(indices = al, mesh = self.Sigma_imp[0].mesh)) for a,al in gf_struct_proj[ish] ], make_copies = False )
                 for ish in xrange(self.n_shells)]
        for ish in range(self.n_shells): Gproj[ish].zero()

        Msh = [x.real for x in self.Sigma_imp[0].mesh]
        n_om = len(Msh)

        DOS = {}
        for bn in self.block_names[self.SO]:
            DOS[bn] = numpy.zeros([n_om],numpy.float_)

        DOSproj     = [ {} for ish in range(self.n_shells) ]
        DOSproj_orb = [ {} for ish in range(self.n_shells) ]
        for ish in range(self.n_shells):
            for bn in self.block_names[self.SO]:
                dl = self.shells[ish][3]
                DOSproj[ish][bn] = numpy.zeros([n_om],numpy.float_)
                DOSproj_orb[ish][bn] = numpy.zeros([n_om,dl,dl],numpy.float_)

        ikarray=numpy.array(range(self.n_k))

        for ik in mpi.slice_array(ikarray):

            S = self.lattice_gf_realfreq(ik=ik,mu=mu,broadening=broadening)
            S *= self.bz_weights[ik]

            # non-projected DOS
            for iom in range(n_om):
                for sig,gf in S: DOS[sig][iom] += gf.data[iom,:,:].imag.trace()/(-3.1415926535)

            #projected DOS:
            for ish in xrange(self.n_shells):
                tmp = Gproj[ish].copy()
                for ir in xrange(self.n_parproj[ish]):
                    for sig,gf in tmp: tmp[sig] << self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                    Gproj[ish] += tmp

        # collect data from mpi:
        for sig in DOS:
            DOS[sig] = mpi.all_reduce(mpi.world,DOS[sig],lambda x,y : x+y)
        for ish in xrange(self.n_shells):
            Gproj[ish] << mpi.all_reduce(mpi.world,Gproj[ish],lambda x,y : x+y)
        mpi.barrier()

        if (self.symm_op!=0): Gproj = self.Symm_par.symmetrize(Gproj)

        # rotation to local coord. system:
        if (self.use_rotations):
            for ish in xrange(self.n_shells):
                for sig,gf in Gproj[ish]: Gproj[ish][sig] << self.rotloc_all(ish,gf,direction='toLocal')

        for ish in range(self.n_shells):
            for sig,gf in Gproj[ish]:
                for iom in range(n_om): DOSproj[ish][sig][iom] += gf.data[iom,:,:].imag.trace()/(-3.1415926535)
                DOSproj_orb[ish][sig][:,:,:] += gf.data[:,:,:].imag / (-3.1415926535)


        if (mpi.is_master_node()):
            # output to files
            for bn in self.block_names[self.SO]:
                f=open('./DOScorr%s.dat'%bn, 'w')
                for i in range(n_om): f.write("%s    %s\n"%(Msh[i],DOS[bn][i]))
                f.close()

                # partial
                for ish in range(self.n_shells):
                    f=open('DOScorr%s_proj%s.dat'%(bn,ish),'w')
                    for i in range(n_om): f.write("%s    %s\n"%(Msh[i],DOSproj[ish][bn][i]))
                    f.close()

                    for i in range(self.shells[ish][3]):
                        for j in range(i,self.shells[ish][3]):
                            Fname = './DOScorr'+bn+'_proj'+str(ish)+'_'+str(i)+'_'+str(j)+'.dat'
                            f=open(Fname,'w')
                            for iom in range(n_om): f.write("%s    %s\n"%(Msh[iom],DOSproj_orb[ish][bn][iom,i,j]))
                            f.close()




    def spaghettis(self,broadening,shift=0.0,plot_range=None, ishell=None, invert_Akw=False, fermi_surface=False):
        """ Calculates the correlated band structure with a real-frequency self energy.
            ATTENTION: Many things from the original input file are are overwritten!!!"""

        assert hasattr(self,"Sigma_imp"), "Set Sigma First!!"
        thingstoread = ['n_k','n_orbitals','proj_mat','hopping','n_parproj','proj_mat_pc']
        read_value = self.read_input_from_hdf(subgrp=self.bands_data,things_to_read=thingstoread)
        if not read_value: return read_value

        if fermi_surface: ishell=None

        # print hamiltonian for checks:
        if ((self.SP==1)and(self.SO==0)):
            f1=open('hamup.dat','w')
            f2=open('hamdn.dat','w')

            for ik in xrange(self.n_k):
                for i in xrange(self.n_orbitals[ik,0]):
                    f1.write('%s    %s\n'%(ik,self.hopping[ik,0,i,i].real))
                for i in xrange(self.n_orbitals[ik,1]):
                    f2.write('%s    %s\n'%(ik,self.hopping[ik,1,i,i].real))
                f1.write('\n')
                f2.write('\n')
            f1.close()
            f2.close()
        else:
            f=open('ham.dat','w')
            for ik in xrange(self.n_k):
                for i in xrange(self.n_orbitals[ik,0]):
                    f.write('%s    %s\n'%(ik,self.hopping[ik,0,i,i].real))
                f.write('\n')
            f.close()


        #=========================================
        # calculate A(k,w):

        mu = self.chemical_potential
        bln = self.block_names[self.SO]

        # init DOS:
        M = [x.real for x in self.Sigma_imp[0].mesh]
        n_om = len(M)

        if plot_range is None:
            om_minplot = M[0]-0.001
            om_maxplot = M[n_om-1] + 0.001
        else:
            om_minplot = plot_range[0]
            om_maxplot = plot_range[1]

        if (ishell is None):
            Akw = {}
            for ibn in bln: Akw[ibn] = numpy.zeros([self.n_k, n_om ],numpy.float_)
        else:
            Akw = {}
            for ibn in bln: Akw[ibn] = numpy.zeros([self.shells[ishell][3],self.n_k, n_om ],numpy.float_)

        if fermi_surface:
            om_minplot = -2.0*broadening
            om_maxplot =  2.0*broadening
            Akw = {}
            for ibn in bln: Akw[ibn] = numpy.zeros([self.n_k,1],numpy.float_)

        if not (ishell is None):
            GFStruct_proj =  [ (al, range(self.shells[ishell][3])) for al in bln ]
            Gproj = BlockGf(name_block_generator = [ (a,GfReFreq(indices = al, mesh = self.Sigma_imp[0].mesh)) for a,al in GFStruct_proj ], make_copies = False)
            Gproj.zero()

        for ik in xrange(self.n_k):

            S = self.lattice_gf_realfreq(ik=ik,mu=mu,broadening=broadening)
            if (ishell is None):
                # non-projected A(k,w)
                for iom in range(n_om):
                    if (M[iom]>om_minplot) and (M[iom]<om_maxplot):
                        if fermi_surface:
                            for sig,gf in S: Akw[sig][ik,0] += gf.data[iom,:,:].imag.trace()/(-3.1415926535) * (M[1]-M[0])
                        else:
                            for sig,gf in S: Akw[sig][ik,iom] += gf.data[iom,:,:].imag.trace()/(-3.1415926535)
                            Akw[sig][ik,iom] += ik*shift                       # shift Akw for plotting in xmgrace


            else:
                # projected A(k,w):
                Gproj.zero()
                tmp = Gproj.copy()
                for ir in xrange(self.n_parproj[ishell]):
                    for sig,gf in tmp: tmp[sig] << self.downfold_pc(ik,ir,ishell,sig,S[sig],gf)
                    Gproj += tmp

                # TO BE FIXED:
                # rotate to local frame
                #if (self.use_rotations):
                #    for sig,gf in Gproj: Gproj[sig] << self.rotloc(0,gf,direction='toLocal')

                for iom in range(n_om):
                    if (M[iom]>om_minplot) and (M[iom]<om_maxplot):
                        for ish in range(self.shells[ishell][3]):
                            for ibn in bln:
                                Akw[ibn][ish,ik,iom] = Gproj[ibn].data[iom,ish,ish].imag/(-3.1415926535)


        # END k-LOOP
        if (mpi.is_master_node()):
            if (ishell is None):

                for ibn in bln:
                    # loop over GF blocs:

                    if (invert_Akw):
                        maxAkw=Akw[ibn].max()
                        minAkw=Akw[ibn].min()


                    # open file for storage:
                    if fermi_surface:
                        f=open('FS_'+ibn+'.dat','w')
                    else:
                        f=open('Akw_'+ibn+'.dat','w')

                    for ik in range(self.n_k):
                        if fermi_surface:
                            if (invert_Akw):
                                Akw[ibn][ik,0] = 1.0/(minAkw-maxAkw)*(Akw[ibn][ik,0] - maxAkw)
                            f.write('%s    %s\n'%(ik,Akw[ibn][ik,0]))
                        else:
                            for iom in range(n_om):
                                if (M[iom]>om_minplot) and (M[iom]<om_maxplot):
                                    if (invert_Akw):
                                        Akw[ibn][ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ibn][ik,iom] - maxAkw)
                                    if (shift>0.0001):
                                        f.write('%s      %s\n'%(M[iom],Akw[ibn][ik,iom]))
                                    else:
                                        f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ibn][ik,iom]))

                            f.write('\n')

                    f.close()

            else:
                for ibn in bln:
                    for ish in range(self.shells[ishell][3]):

                        if (invert_Akw):
                            maxAkw=Akw[ibn][ish,:,:].max()
                            minAkw=Akw[ibn][ish,:,:].min()

                        f=open('Akw_'+ibn+'_proj'+str(ish)+'.dat','w')

                        for ik in range(self.n_k):
                            for iom in range(n_om):
                                if (M[iom]>om_minplot) and (M[iom]<om_maxplot):
                                    if (invert_Akw):
                                        Akw[ibn][ish,ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ibn][ish,ik,iom] - maxAkw)
                                    if (shift>0.0001):
                                        f.write('%s      %s\n'%(M[iom],Akw[ibn][ish,ik,iom]))
                                    else:
                                        f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ibn][ish,ik,iom]))

                            f.write('\n')

                        f.close()


    def constr_Sigma_real_axis(self, filename, hdf=True, hdf_dataset='SigmaReFreq',n_om=0,orb=0, tol_mesh=1e-6):
        """Uses Data from files to construct Sigma (or GF)  on the real axis."""

        if not hdf: # then read sigma from text files

            # first get the mesh out of any one of the files:
            bl = self.gf_struct_solver[orb].items()[0][0] # block name
            ol = self.gf_struct_solver[orb].items()[0][1] # list of orbital indices
            if (len(ol)==1): # if blocks are of size one
                Fname = filename+'_'+bl+'.dat'
            else:
                Fname = filename+'_'+bl+'/'+str(ol[0])+'_'+str(ol[0])+'.dat'

            R = read_fortran_file(Fname)
            mesh = numpy.zeros([n_om],numpy.float_)
            try:
                for i in xrange(n_om):
                    mesh[i] = R.next()
                    sk = R.next()
                    sk = R.next()

            except StopIteration : # a more explicit error if the file is corrupted.
                raise "SumkLDA.read_Sigma_ME : reading mesh failed!"
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
                            raise "SumkLDA.read_Sigma_ME : reading Sigma from file failed!"
                        R.close()


        else: # read sigma from hdf

            omega_min=0.0
            omega_max=0.0
            n_om=0
            if (mpi.is_master_node()):
                ar = HDFArchive(filename)
                SigmaME = ar[hdf_dataset]
                del ar
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



    def partial_charges(self,beta=40):
        """Calculates the orbitally-resolved density matrix for all the orbitals considered in the input.
           The theta-projectors are used, hence case.parproj data is necessary"""


        #thingstoread = ['Dens_Mat_below','N_parproj','Proj_Mat_pc','rotmat_all']
        #read_value = self.read_input_from_HDF(SubGrp=self.par_proj_data,thingstoread=thingstoread)
        read_value = self.read_par_proj_input_from_hdf()
        if not read_value: return read_value
        if self.symm_op: self.Symm_par = Symmetry(self.hdf_file,subgroup=self.symm_par_data)

        # Density matrix in the window
        bln = self.block_names[self.SO]
        ntoi = self.names_to_ind[self.SO]
        self.dens_mat_window = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.n_shells)]
                                 for isp in range(len(bln)) ]    # init the density matrix

        mu = self.chemical_potential
        GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in bln ]  for i in xrange(self.n_shells) ]
        if hasattr(self,"Sigma_imp"):
            Gproj = [BlockGf(name_block_generator = [ (a,GfImFreq(indices = al, mesh = self.Sigma_imp[0].mesh)) for a,al in GFStruct_proj[ish] ], make_copies = False)
                     for ish in xrange(self.n_shells)]
            beta = self.Sigma_imp[0].mesh.beta
        else:
            Gproj = [BlockGf(name_block_generator = [ (a,GfImFreq(indices = al, beta = beta)) for a,al in GFStruct_proj[ish] ], make_copies = False)
                     for ish in xrange(self.n_shells)]

        for ish in xrange(self.n_shells): Gproj[ish].zero()

        ikarray=numpy.array(range(self.n_k))
        #print mpi.rank, mpi.slice_array(ikarray)
        #print "K-Sum starts on node",mpi.rank," at ",datetime.now()

        for ik in mpi.slice_array(ikarray):
            #print mpi.rank, ik, datetime.now()
            S = self.lattice_gf_matsubara(ik=ik,mu=mu,beta=beta)
            S *= self.bz_weights[ik]

            for ish in xrange(self.n_shells):
                tmp = Gproj[ish].copy()
                for ir in xrange(self.n_parproj[ish]):
                    for sig,gf in tmp: tmp[sig] << self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                    Gproj[ish] += tmp

        #print "K-Sum done on node",mpi.rank," at ",datetime.now()
        #collect data from mpi:
        for ish in xrange(self.n_shells):
            Gproj[ish] << mpi.all_reduce(mpi.world,Gproj[ish],lambda x,y : x+y)
        mpi.barrier()

        #print "Data collected on node",mpi.rank," at ",datetime.now()

        # Symmetrisation:
        if (self.symm_op!=0): Gproj = self.Symm_par.symmetrize(Gproj)
        #print "Symmetrisation done on node",mpi.rank," at ",datetime.now()

        for ish in xrange(self.n_shells):

            # Rotation to local:
            if (self.use_rotations):
                for sig,gf in Gproj[ish]: Gproj[ish][sig] << self.rotloc_all(ish,gf,direction='toLocal')

            isp = 0
            for sig,gf in Gproj[ish]: #dmg.append(Gproj[ish].density()[sig])
                self.dens_mat_window[isp][ish] = Gproj[ish].density()[sig]
                isp+=1

        # add Density matrices to get the total:
        dens_mat = [ [ self.dens_mat_below[ntoi[bln[isp]]][ish]+self.dens_mat_window[isp][ish] for ish in range(self.n_shells)]
                     for isp in range(len(bln)) ]

        return dens_mat
