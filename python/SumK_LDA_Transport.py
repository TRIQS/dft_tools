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

#=======================================================================================================================
# #################################################################
# Code for Transport/Optic calculations based on SumK_LDA... class
# by Xiaoyu Deng <xiaoyu.deng@gmail.com>
# #################################################################
#=======================================================================================================================


from types import *
import numpy


# NEW
import pytriqs.utility.dichotomy as Dichotomy
# OLD
#import pytriqs.Base.Utility.Dichotomy as Dichotomy

# NEW
#from pytriqs.gf.local.descriptors import A_Omega_Plus_B
from pytriqs.gf.local import *
# OLD
#from pytrigs.Base.GF_Local.GF import GF
#from pytriqs.Base.GF_local.GFBloc_ImFreq import GFBloc_ImFreq
#from pytriqs.Base.GF_Local.GFBloc_ReFreq import GFBloc_ReFreq
#from pytriqs.Base.GF_Local.GFBloc_ImTime import GFBloc_ImTime
#from pytriqs.Base.GF_Local import GF_Initializers

# NEW
#from pytriqs.operators.operators2 import *
# OLD
#from pytriqs.Solvers.Operators import *

# NEW
# NOT FOUND
# OLD
#from pytriqs.Base.Utility.myUtils import Sum

# NEW
import pytriqs.utility.mpi as myMPI
# OLD
#import pytriqs.Base.Utility.MPI as myMPI

from datetime import datetime

#from pytriqs.Wien2k.Symmetry import *

# NEW
from pytriqs.applications.dft.sumk_lda import *
from pytriqs.applications.dft.sumk_lda_tools import *
#k in xrange(self OLD
#from pytriqs.Wien2k.SumK_LDA import SumK_LDA
#from pytriqs.Wien2k.SumK_LDA_tools import SumK_LDA_tools

import string
import copy

import SumK_LDA_Transport_Wien2k_input as Wien



def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists" % filename
    for line in open(filename, 'r') :
	for x in line.replace('D', 'E').split() : 
	    yield string.atof(x)

def Read_Fortran_File2(filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists" % filename
    for line in open(filename, 'r') :
        for x in line.replace('D', 'E').split() :
            try:
                yield string.atof(x)
            except GeneratorExit:
                raise
            except:
                yield x.replace('E', 'D')


def fermidis(x):
    return 1.0/(numpy.exp(x)+1)

# OLD
#class TransportEtOptic(SumK_LDA_tools): 
# NEW
class TransportEtOptic(SumkLDATools):
    """transport and optic related functions
    calculates distributions: Tr A(k,w) v(k) A(k, w) v(k) or Tr A(k,w) v(k) A(k, w+Omega) v(k)
    .based on thi for ik in xrange(selfs distribution other properties could be obtained.
	!!! worked for cases with spin orbital interaction.
	!!! for non-SOI, need check! Be careful!.
    """

    def initbyWien(self, wiencase):
        """ read in necessary parameters from wien file.
            symmetries:
            volume:
            velocities:
        """
        
        # k_dep_projection is the general case!.
        assert self.k_dep_projection == 1, "Not implemented!"        
        
        self.Nspinblocs = self.SP + 1 - self.SO
        # suffix of wien2k output file  
        self.blocksuffix=[[""],["up","dn"]]
        
        self.velocities=None
        self.bandwin=None
        self.Vol = None
        self.symm = None
        self.nsymm = None
        if myMPI.is_master_node():
            CASE = Wien.WienStruct(wiencase)
            CASE.readSGsymm(wiencase)
            if(self.Nspinblocs == 1): # paramagnetic , without SOI; or spin polarized with SOI
                self.velocities = [Wien.Velocities(wiencase,bln) for bln in self.blocksuffix[0]]
            else:
                self.velocities = [Wien.Velocities(wiencase, bln) for bln in self.blocksuffix[1] ]
# read in band window for each k points. "CASE.oubwin".
            bandwin=self.bandwinfromwiencase(wiencase=wiencase)
        
            self.bandwin=bandwin
            self.Vol = CASE.VolumePC
            self.symm = CASE.symmcartesian
            self.nsymm = CASE.nsymm
        myMPI.barrier() 
        self.velocities=myMPI.bcast(self.velocities)
        self.bandwin=myMPI.bcast(self.bandwin)
        self.Vol=myMPI.bcast(self.Vol)     
        self.symm=myMPI.bcast(self.symm)
        self.nsymm=myMPI.bcast(self.nsymm)                   


    def Transportdistribution_Boltz(self, wiencase, mshape=None, broadening=0.01, energywindow=None, Nomega=1000, loadpw=False):
        """This is for transport calculation of  Boltzmann theory
            calculate \sum_k Tr vv delta(\omega-enk) which is transportdistribution in Boltzmann theory
            Just use the LDA hamiltonian and velocity, DMFT self energy is not needed.
            
            mshape defines the indices of directions. xx,yy,zz,xy,yz,zx. 
            mshape is 3x3 matrix, mshape[0,0]=1 --> calculate xx,  mshape[1,1]=1 --> calculate yy, mshape[1,2]=1 --> calculate xy,
            by default, xx is calculated.
        """
        if mshape==None:
            mshape=numpy.zeros(9).reshape(3,3)
            mshape[0,0]=1            
        assert mshape.shape == (3, 3), "mshape should be 3x3"
        velocities = self.velocities
        mu = 0
        assert self.k_dep_projection == 1, "k_dep_projection = 0 is NOT implemented!"

        if energywindow is None:
            omminplot = -1.0
            ommaxplot = 1.0
        else:
            omminplot = energywindow[0]
            ommaxplot = energywindow[1]
        deltaomega = (ommaxplot - omminplot) / Nomega

        # transport distribution :output P(\omega)_xy should has the same dimension as defined in mshape.
        self.Pw = numpy.zeros((mshape.sum(), Nomega), dtype=numpy.float_)

        mlist = []
        for ir in xrange(3):
            for ic in xrange(3):
                if(mshape[ir][ic] == 1):
                    mlist.append((ir, ic))
        if loadpw:
            with open("TD_Boltz_mp.dat", "r") as pwin:
                for iw in xrange(Nomega):
                    fstr = pwin.readline().split()
                    aomega = iw * deltaomega + omminplot
                    assert abs(float(fstr[0]) - aomega) <= 1e-8, "mesh not match when load transportdistribution"
                    for ipw in range(mshape.sum()): self.Pw[ipw, iw] = float(fstr[ipw + 1])
                print "Blotz Pw loaded"
                return

        # for ik=0
        ik = 0
        # block_names for green function and self energy
        bln = self.block_names[self.SO]
        ntoi = self.names_to_ind[self.SO]
        
        ikarray = numpy.array(range(self.n_k))
        for isp in range(self.Nspinblocs):
            for ik in myMPI.slice_array(ikarray):
                n_orb = self.n_orbitals[ik][isp]
		mupat = numpy.ones(n_orb, numpy.float_) * mu
                Ms = copy.deepcopy(mupat)
                ind = ntoi[bln[isp]]
		Ms = numpy.diag(self.hopping[ik,ind,0:n_orb,0:n_orb].real) - mupat
                if(ik%100==0):
                    print "ik,isp", ik, isp
                kvel = velocities[isp].vks[ik]

                # in general, bandwindows for Annkw and for velocities are not the same.
                # one should make sure the same window are using before go further. otherwise the matrix size are not match.
                Pwtem = numpy.zeros((mshape.sum(), Nomega), dtype=numpy.float_)
                #symmetry loop
                # this symmetrization could be done first to speed up... To be done.
                for Rmat in self.symm:
                    # get new velocity.
                    Rkvel = copy.deepcopy(kvel.vel)
                    for vnb1 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                        Rkvel[vnb1][vnb1][:] = numpy.dot(Rmat, Rkvel[vnb1][vnb1][:])

                    ipw = 0
                    for (ir, ic) in mlist:
                        for iw in xrange(Nomega):
                            omega = deltaomega * iw + omminplot
                            bmin = max(self.bandwin[isp][ik, 0], kvel.bandwin[0])
                            bmax = min(self.bandwin[isp][ik, 1], kvel.bandwin[1])
                            for ib in xrange(bmax - bmin + 1):
                                ibb = ib + bmin - self.bandwin[isp][ik, 0]
                                ibv = ib + bmin - kvel.bandwin[0]
                                enk = Ms[ibb] - omega
                                Ag = 1.0 / numpy.sqrt(2 * numpy.pi) / broadening * numpy.exp(-enk ** 2 / 2.0 / broadening ** 2)
                                vkr = Rkvel[ibv][ibv][ir]
                                vkc = Rkvel[ibv][ibv][ic]

                                #Pwtem[ipw,iw]+=(vkr*vkc*Ag).real
                                Pwtem[ipw, iw] += vkr * vkc * Ag
                        ipw += 1

                self.Pw += Pwtem * self.bz_weights[ik] / self.nsymm
        myMPI.barrier()
        self.Pw = myMPI.all_reduce(myMPI.world, self.Pw, lambda x, y : x + y)
        # for non-magnetic case, the total weight is doubled because of spin degeneracy. 
        self.Pw *= (2 - self.SP)
       
        # scattering is needed here.
        self.Pw *= 1.0/numpy.pi/2.0/broadening


        if myMPI.is_master_node():
            with open("TD_Boltz.dat", "w") as pwout:
                for iw in xrange(Nomega):
                    omega = deltaomega * iw + omminplot
                    pwout.write(str(omega) + "   ")
                    for i in range(self.Pw.shape[0]):
                        pwout.write(str(self.Pw[i, iw]) + "    ")
                    pwout.write("\n")

    def Transportdistribution(self, wiencase, mshape=None, broadening=0.00, energywindow=None, loadpw=False):
        """calculate Tr A(k,w) v(k) A(k, w) v(k). 
        mshape defines the indices of directions. xx,yy,zz,xy,yz,zx. 
        mshape is 3x3 matrix, mshape[0,0]=1 --> calculate xx,  mshape[1,1]=1 --> calculate yy, mshape[1,2]=1 --> calculate xy,
        by default, xx is calculated.
        """
        if mshape==None:
            mshape=numpy.zeros(9).reshape(3,3)
            mshape[0,0]=1
        assert mshape.shape == (3, 3), "mshape should be 3x3"
        assert hasattr(self, "Sigma_imp"), "Set Sigma First!"        
        velocities = self.velocities
        mu = self.chemical_potential
        if myMPI.is_master_node():
            print "Chemical_Potential", mu
        # use k-dependent-projections.
        assert self.k_dep_projection == 1, "Not implemented!"

        # form self energy from impurity self energy and double counting term.
        stmp = self.add_dc()

        #set mesh and energy range for spectrals functions
        # one could construct self energy within only a small energy range when calculating transports.
        M = [x for x in self.Sigma_imp[0].mesh]
        N_om = len(M)
        if energywindow is None:
            omminplot = M[0].real - 0.001
            ommaxplot = M[N_om - 1].real + 0.001
        else:
            omminplot = energywindow[0]
            ommaxplot = energywindow[1]
        # set mesh for Pw, only mesh in focused energyrange is needed. Mpw is just the index of mesh need in M
        Mpw = [i for i in xrange(len(M)) if (M[i].real > omminplot and  M[i].real < ommaxplot)]

        # output P(\omega)_xy should has the same dimension as defined in mshape.
        self.Pw = numpy.zeros((mshape.sum(), N_om), dtype=numpy.float)

        mlist = []
        for ir in xrange(3):
            for ic in xrange(3):
                if(mshape[ir][ic] == 1):
                    mlist.append((ir, ic))

        ik = 0
        # will construct G in the end; don't be mislead by
        # nomenclature; S is sometimes sigma, sometimes sigma^-1, etc.


        bln = self.block_names[self.SO]
        ntoi = self.names_to_ind[self.SO]

        S = BlockGf(name_block_generator=[(bln[isp], GfReFreq(indices=range(self.n_orbitals[ik][isp]), mesh=self.Sigma_imp[0].mesh)) for isp in range(self.Nspinblocs) ], make_copies=False)
        mupat = [numpy.identity(self.n_orbitals[ik][isp], numpy.complex_) * mu for isp in range(self.Nspinblocs)]   # construct mupat
        Annkw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], N_om), dtype=numpy.complex_) for isp in range(self.Nspinblocs)]

        ikarray = numpy.array(range(self.n_k))
        for ik in myMPI.slice_array(ikarray):
           
        #for ik in xrange(self.n_k):
            unchangesize = all([ self.n_orbitals[ik][isp] == mupat[isp].shape[0] for isp in range(self.Nspinblocs)])
            if (not unchangesize):
                # recontruct green functions.
                S = BlockGf(name_block_generator=[(bln[isp], GfReFreq(indices=range(self.n_orbitals[ik][isp]),
                                                                        mesh=self.Sigma_imp[0].mesh))
                                                       for isp in range(self.Nspinblocs) ],
                               make_copies=False)
                                   # change size of mupat
                mupat = [numpy.identity(self.n_orbitals[ik][isp], numpy.complex_) * mu for isp in range(self.Nspinblocs)]   # construct mupat
                #Annkw=numpy.zeros((self.n_orbitals[ik],self.n_orbitals[ik],N_om),dtype=numpy.complex_)
                Annkw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], N_om), dtype=numpy.complex_) for isp in range(self.Nspinblocs)]
            # get lattice green functions.
           # S <<= A_Omega_Plus_B(A=1, B=1j * broadening)
	   
	    S <<= 1*Omega+1j*broadening

            Ms = copy.deepcopy(mupat)
            for ibl in range(self.Nspinblocs):
                n_orb = self.n_orbitals[ik][ibl]
		ind = ntoi[bln[ibl]]
		Ms[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb].real - mupat[ibl]
            S -= Ms

            tmp = S.copy()
            for icrsh in xrange(self.n_corr_shells):
                for sig, gf in tmp: tmp[sig] <<= self.upfold(ik, icrsh, sig, stmp[icrsh][sig], gf)
                S -= tmp

            S.invert()
            
            #hence we have A(k,\omega)_nn' for a special k points.
            for isp in range(self.Nspinblocs):
               # Annkw[isp].real = -copy.deepcopy(S[self.block_names[self.SO][isp]]._data.array).imag / numpy.pi
		Annkw[isp].real = -copy.deepcopy(S[self.block_names[self.SO][isp]].data.swapaxes(0,1).swapaxes(1,2)).imag / numpy.pi
       
#               A=-1/pi*Im G

            # for different spin velocties might be different
            for isp in range(self.Nspinblocs):
                if(ik%100==0):
                    print "ik,isp", ik, isp                 
                kvel = velocities[isp].vks[ik]

                # in general, bandwindows for Annkw and for velocities are not the same.
                # one should make sure the same window are using before go further. otherwise the matrix size are not match.

                Pwtem = numpy.zeros((mshape.sum(), N_om), dtype=numpy.float_)

                #symmetry loop
                # how to symmetrize this part???
                for Rmat in self.symm:
                    # get new velocity.
                    Rkvel = copy.deepcopy(kvel.vel)
                    for vnb1 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                        for vnb2 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                            Rkvel[vnb1][vnb2][:] = numpy.dot(Rmat, Rkvel[vnb1][vnb2][:])
                    ipw = 0
                    bmin = max(self.bandwin[isp][ik, 0], kvel.bandwin[0])
                    bmax = min(self.bandwin[isp][ik, 1], kvel.bandwin[1])
                    Astart = bmin - self.bandwin[isp][ik, 0]
                    Aend = bmax - self.bandwin[isp][ik, 0] + 1
                    vstart = bmin - kvel.bandwin[0]
                    vend = bmax - kvel.bandwin[0] + 1
                    for (ir, ic) in mlist:
                        #for iw in xrange(N_om):
                        for iw in Mpw:
                            #if (M[iw]>omminplot) and (M[iw]<ommaxplot):
                                # here use bandwin to construct match matrix for A and velocity.
                            Annkwt = Annkw[isp][Astart:Aend, Astart:Aend, iw]
                            Rkveltr = Rkvel[vstart:vend, vstart:vend, ir]
                            Rkveltc = Rkvel[vstart:vend, vstart:vend, ic]
                            #print Annkwt.shape,Rkvel[...,ir].shape
                            Pwtem[ipw, iw] += numpy.dot(numpy.dot(numpy.dot(Rkveltr, Annkwt), Rkveltc), Annkwt).trace().real
                        ipw += 1

                # k sum and spin sum.
                self.Pw += Pwtem * self.bz_weights[ik] / self.nsymm



        self.Pw = myMPI.all_reduce(myMPI.world, self.Pw, lambda x, y : x + y)
        # for non-magnetic case, the total weight is doubled because of spin degeneracy. 
        self.Pw *= (2 - self.SP)
 
        if myMPI.is_master_node():
            with open("TD_DMFT.dat", "w") as pwout:
                for iw in xrange(N_om):
                    if (M[iw].real > omminplot) and (M[iw].real < ommaxplot):
                        pwout.write(str(M[iw].real) + "   ")
                        for i in range(self.Pw.shape[0]):
                            pwout.write(str(self.Pw[i, iw]) + "    ")
                        pwout.write("\n")

    def OpticDistribution(self, wiencase, mshape=None, broadening=0.01, energywindow=None, Qmesh=[0.5], Beta=50, loadpw=False):
        """calculate Tr A(k,w) v(k) A(k, w+q) v(k) and optics.
        energywindow is the regime for omega integral
        Qmesh contains the frequencies of the optic conductivitity. I repin the Qmesh to the self-energy mesh, 
        so the exact value might not exactly the same as given in the list.
        
        mshape defines the indices of directions. xx,yy,zz,xy,yz,zx. 
        mshape is 3x3 matrix, mshape[0,0]=1 --> calculate xx,  mshape[1,1]=1 --> calculate yy, mshape[1,2]=1 --> calculate xy,
        by default, xx is calculated.
        """
        assert mshape.shape == (3, 3), "mshape should be 3x3"

        assert hasattr(self, "Sigma_imp"), "Set Sigma First!"
        assert ((self.SP == 0) and (self.SO == 0)), "For SP and SO implementation of spaghettis has to be changed!"
        velocities = self.velocities
        # calculate A(k,w):
        mu = self.chemical_potential

        #we need this somehow for k_dep_projections. So we have to face the problem that the size of A(k,\omega) will
        #change, and more, the band index for the A(k,\omega) matrix is not known yet.

        # use k-dependent-projections.
        assert self.k_dep_projection == 1, "Not implemented!"

        # form self energy from impurity self energy and double counting term.
        stmp = self.add_dc()

        #set mesh and energyrange.
        M = [x for x in self.Sigma_imp[0].mesh]
        deltaM = numpy.abs(M[0] - M[1])
        N_om = len(M)
        if energywindow is None:
            omminplot = M[0].real - 0.001
            ommaxplot = M[N_om - 1].real + 0.001
        else:
            omminplot = energywindow[0]
            ommaxplot = energywindow[1]

        # define exact mesh for optic conductivity
        Qmesh_ex = [int(x / deltaM) for x in Qmesh]
        if myMPI.is_master_node():
            print "Qmesh   ", Qmesh
            print "mesh interval in self energy  ", deltaM
            print "Qmesh / mesh interval  ", Qmesh_ex

        # output P(\omega)_xy should has the same dimension as defined in mshape.
        self.Pw_optic = numpy.zeros((mshape.sum(), len(Qmesh), N_om), dtype=numpy.float_)

        mlist = []
        for ir in xrange(3):
            for ic in xrange(3):
                if(mshape[ir][ic] == 1):
                    mlist.append((ir, ic))
        ik = 0
        
        bln = self.block_names[self.SO]
        ntoi = self.names_to_ind[self.SO]

        S = BlockGf(name_block_generator=[(bln[isp], GfReFreq(indices=range(self.n_orbitals[ik][isp]),
                                                                     mesh=self.Sigma_imp[0].mesh))
                                                    for isp in range(self.Nspinblocs) ],
                            make_copies=False)
        mupat = [numpy.identity(self.n_orbitals[ik][isp], numpy.complex_) * mu for isp in range(self.Nspinblocs)]   # construct mupat
        Annkw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], N_om), dtype=numpy.complex_) for isp in range(self.Nspinblocs)]

        ikarray = numpy.array(range(self.n_k))
        for ik in myMPI.slice_array(ikarray):
            unchangesize = all([ self.n_orbitals[ik][isp] == mupat[isp].shape[0] for isp in range(self.Nspinblocs)])
            if (not unchangesize):
                # recontruct green functions.
                S = BlockGf(name_block_generator=[(bln[isp], GfReFreq(indices=range(self.n_orbitals[ik][isp]),
                                                                        mesh=self.Sigma_imp[0].mesh))
                                                       for isp in range(self.Nspinblocs) ],
                               make_copies=False)

                # S = GF(name_block_generator=[ (s, GFBloc_ReFreq(Indices=BS, Mesh=self.Sigma_imp[0].mesh)) for s in ['up', 'down'] ], Copy=False)
                # mupat = numpy.identity(self.n_orbitals[ik], numpy.complex_)                                     # change size of mupat
                mupat = [numpy.identity(self.n_orbitals[ik][isp], numpy.complex_) * mu for isp in range(self.Nspinblocs)]   # construct mupat
                # mupat *= mu

            #set a temporary array storing spectral functions with band index. Note, usually we should have spin index
                #Annkw=numpy.zeros((self.n_orbitals[ik],self.n_orbitals[ik],N_om),dtype=numpy.complex_)
                Annkw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], N_om), dtype=numpy.complex_) for isp in range(self.Nspinblocs)]

            # get lattice green functions.
            # S <<= A_Omega_Plus_B(A=1, B=1j * broadening)
            
	    S <<= 1*Omega + 1j*broadening

            Ms = copy.deepcopy(mupat)
            for ibl in range(self.Nspinblocs):
                ind = ntoi[bln[ibl]]
		n_orb = self.n_orbitals[ik][ibl]
		Ms[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb].real - mupat[ibl]
            S -= Ms

            tmp = S.copy()    # init temporary storage
            ## substract self energy
            for icrsh in xrange(self.n_corr_shells):
                for sig, gf in tmp: tmp[sig] <<= self.upfold(ik, icrsh, sig, stmp[icrsh][sig], gf)
                S -= tmp

            S.invert()

            for isp in range(self.Nspinblocs):
                Annkw[isp].real = -copy.deepcopy(S[self.block_names[self.SO][isp]].data.swapaxes(0,1).swapaxes(1,2)).imag / numpy.pi

            for isp in range(self.Nspinblocs):
                if(ik%100==0):
                    print "ik,isp", ik, isp
                kvel = velocities[isp].vks[ik]

                Pwtem = numpy.zeros((mshape.sum(), len(Qmesh_ex), N_om), dtype=numpy.float_)

                #symmetry loop
                for Rmat in self.symm:
                    # get new velocity.
                    Rkvel = copy.deepcopy(kvel.vel)
                    for vnb1 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                        for vnb2 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                            Rkvel[vnb1][vnb2][:] = numpy.dot(Rmat, Rkvel[vnb1][vnb2][:])
                    ipw = 0
                    for (ir, ic) in mlist:
                        for iw in xrange(N_om):
				
                            if(M[iw].real > 5.0 / Beta):
                                continue
                            for iq in range(len(Qmesh_ex)):
                                #if(Qmesh_ex[iq]==0 or iw+Qmesh_ex[iq]>=N_om ):
                                # here use fermi distribution to truncate self energy mesh.
				if(Qmesh_ex[iq] == 0 or iw + Qmesh_ex[iq] >= N_om  
                                   or M[iw].real + Qmesh[iq] < -10.0 / Beta or M[iw].real >10.0 / Beta):
                                    continue
				
                                if (M[iw].real > omminplot) and (M[iw].real < ommaxplot):
                                # here use bandwin to construct match matrix for A and velocity.
                                    bmin = max(self.bandwin[isp][ik, 0], kvel.bandwin[0])
                                    bmax = min(self.bandwin[isp][ik, 1], kvel.bandwin[1])
                                    Astart = bmin - self.bandwin[isp][ik, 0]
                                    Aend = bmax - self.bandwin[isp][ik, 0] + 1
                                    vstart = bmin - kvel.bandwin[0]
                                    vend = bmax - kvel.bandwin[0] + 1
                                    Annkwl = Annkw[isp][Astart:Aend, Astart:Aend, iw]
                                    Annkwr = Annkw[isp][Astart:Aend, Astart:Aend, iw + Qmesh_ex[iq]]
                                    Rkveltr = Rkvel[vstart:vend, vstart:vend, ir]
                                    Rkveltc = Rkvel[vstart:vend, vstart:vend, ic]
                                #print Annkwt.shape,Rkvel[...,ir].shape
                                    Pwtem[ipw, iq, iw] += numpy.dot(numpy.dot(numpy.dot(Rkveltr, Annkwl), Rkveltc), Annkwr).trace().real
                        ipw += 1

                # k sum and spin sum.
                self.Pw_optic += Pwtem * self.bz_weights[ik] / self.nsymm



        self.Pw_optic = myMPI.all_reduce(myMPI.world, self.Pw_optic, lambda x, y : x + y)
        self.Pw_optic *= (2 - self.SP)




        # just back up TD_optic data
        if myMPI.is_master_node():
            with open("TD_Optic_DMFT.dat", "w") as pwout:
                #shape
                L1,L2,L3=self.Pw_optic.shape
                pwout.write("%s  %s  %s\n"%(L1,L2,L3))
                #dump Qmesh
                Qmeshr=[i*deltaM for i in Qmesh_ex]
                #dump self energy mesh
                #dump Pw_optic
                for iq in xrange(L2):
                    pwout.write(str(Qmeshr[iq])+"  ")
                pwout.write("\n")
                for iw in xrange(L3):
                    pwout.write(str(M[iw].real)+"  ")
                pwout.write("\n")               
                for i in xrange(L1):
                    for iq in xrange(L2):
                        for iw in xrange(L3):
                            pwout.write(str(self.Pw_optic[i, iq, iw]) + "    ")
                pwout.write("\n")

        # sum over omega to get optic conductivity for ik in xrange(self
        if myMPI.is_master_node():
            OpticConductivity = numpy.zeros((mshape.sum(), len(Qmesh)), dtype=numpy.float_)
            for im in range(mshape.sum()):
                for iq in range(len(Qmesh)):
                    for iw in xrange(N_om):
                        omegaT = M[iw].real * Beta
                        omega_aug = Qmesh_ex[iq] * deltaM
                        OpticConductivity[im, iq] += self.Pw_optic[im, iq, iw] * (fermidis(omegaT) - fermidis(omegaT + omega_aug * Beta)) / omega_aug
            OpticConductivity *= deltaM
            OpticConductivity *= 10700 / self.Vol
            with open("Optic_con.dat", "wOptic_con") as opt:
                for iq in range(len(Qmesh_ex)):
                    opt.write(str(Qmesh_ex[iq] * deltaM) + "   ")
                    for im in range(mshape.sum()):
                        opt.write(str(OpticConductivity[im, iq]) + "   ")
                    opt.write("\n")
    
    def loadOpticTD(self,OpticTDFile="TD_Optic_DMFT.dat",Beta=40):
        """ load optic conductivity distribution and calculate Optical Conductivty
        """
        if myMPI.is_master_node():
            with open(OpticTDFile,"r") as pw:
                L1,L2,L3=(int(i) for i in pw.readline().split())
                #QMeshr=numpy.zeros(L2,dtype=numpy.float)
                #M=numpy.zeros(L3,dtype=numpy.float)
                #Pw_optic=numpy.zeros((L1,L2,L3), dtype=numpy.float)
                
                QMeshr=numpy.array([float(i) for i in pw.readline().split()])
                M=numpy.array([float(i) for i in pw.readline().split()])
                Pw_optic=numpy.array([float(i) for i in pw.readline().split()]).reshape(L1,L2,L3)                            
            
            OpticConductivity = numpy.zeros((L1, L2), dtype=numpy.float)
            deltaM=M[1]-M[0]
            for im in xrange(L1):
                for iq in xrange(L2):
                    for iw in xrange(L3):
                        omegaT = M[iw] * Beta
                        omega_aug = QMeshr[iq]
                        OpticConductivity[im, iq] += Pw_optic[im, iq, iw] * (fermidis(omegaT) - fermidis(omegaT + omega_aug * Beta)) / omega_aug
            OpticConductivity *= deltaM
            ## transform to standard unit as in resistivity
            OpticConductivity *= 10700 / self.Vol
            ##
            with open("Optic_con.dat", "w") as opt:
                for iq in xrange(L2):
                    opt.write(str(QMeshr[iq]) + "   ")
                    for im in xrange(L1):
                        opt.write(str(OpticConductivity[im, iq]) + "   ")
                    opt.write("\n")        
        


    def OpticDistribution_LDA(self, wiencase, mshape=None, broadening=0.01, energywindow=None, Qmesh=[0.5], Beta=50, loadpw=False):
        """calculate Tr A(k,w) v(k) A(k, w+q) v(k) and optics. A constant self-energy is used to mimick noninteracting case.
        It is not the best way to calculate optic for LDA. Just to compare.
        energywindow is the regime for omega integral
        Qmesh contains the frequencies of the optic conductivitity. I repin the Qmesh to the self-energy mesh, 
        so the exact value might not exactly the same as given in the list.
        
        mshape defines the indices of directions. xx,yy,zz,xy,yz,zx. 
        mshape is 3x3 matrix, mshape[0,0]=1 --> calculate xx,  mshape[1,1]=1 --> calculate yy, mshape[1,2]=1 --> calculate xy,
        by default, xx is calculated.
        """
        assert mshape.shape == (3, 3), "mshape should be 3x3"

        assert hasattr(self, "Sigma_imp"), "Set Sigma First!"
        assert ((self.SP == 0) and (self.SO == 0)), "For SP and SO implementation of spaghettis has to be changed!"
        velocities = self.velocities
        # calculate A(k,w):
        mu = self.chemical_potential

        #we need this somehow for k_dep_projections. So we have to face the problem that the size of A(k,\omega) will
        #change, and more, the band index for the A(k,\omega) matrix is not known yet.

        # use k-dependent-projections.
        assert self.k_dep_projection == 1, "Not implemented!"

        # form self energy from impurity self energy and double counting term.
        stmp = self.add_dc()

        #set mesh and energyrange.
        M = [x for x in self.Sigma_imp[0].mesh]
        deltaM = numpy.abs(M[0] - M[1])
        N_om = len(M)
        if energywindow is None:
            omminplot = M[0] - 0.001
            ommaxplot = M[N_om - 1] + 0.001
        else:
            omminplot = energywindow[0]
            ommaxplot = energywindow[1]

        # define exact mesh for optic conductivity
        Qmesh_ex = [int(x / deltaM) for x in Qmesh]
        if myMPI.is_master_node():
            print "Qmesh   ", Qmesh
            print "mesh interval in self energy  ", deltaM
            print "Qmesh / mesh interval  ", Qmesh_ex

        # output P(\omega)_xy should has the same dimension as defined in mshape.
        self.Pw_optic = numpy.zeros((mshape.sum(), len(Qmesh), N_om), dtype=numpy.float_)

        mlist = []
        for ir in xrange(3):
            for ic in xrange(3):
                if(mshape[ir][ic] == 1):
                    mlist.append((ir, ic))
        ik = 0
        
        bln = self.block_names[self.SO]
        ntoi = self.names_to_ind[self.SO]

        S = BlockGf(name_block_generator=[(bln[isp], GfReFreq(indices=range(self.n_orbitals[ik][isp]),
                                                                     mesh=self.Sigma_imp[0].mesh))
                                                    for isp in range(self.Nspinblocs) ],
                            make_copies=False)
        mupat = [numpy.identity(self.n_orbitals[ik][isp], numpy.complex_) * mu for isp in range(self.Nspinblocs)]   # construct mupat
        Annkw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], N_om), dtype=numpy.complex_) for isp in range(self.Nspinblocs)]

        ikarray = numpy.array(range(self.n_k))
        for ik in myMPI.slice_array(ikarray):
            unchangesize = all([ self.n_orbitals[ik][isp] == mupat[isp].shape[0] for isp in range(self.Nspinblocs)])
            if (not unchangesize):
                # recontruct green functions.
                S = BlockGf(name_block_generator=[(bln[isp], GfReFreq(indices=range(self.n_orbitals[ik][isp]),
                                                                        mesh=self.Sigma_imp[0].mesh))
                                                       for isp in range(self.Nspinblocs) ],
                               make_copies=False)

                # S = GF(name_block_generator=[ (s, GFBloc_ReFreq(Indices=BS, Mesh=self.Sigma_imp[0].mesh)) for s in ['up', 'down'] ], Copy=False)
                # mupat = numpy.identity(self.n_orbitals[ik], numpy.complex_)                                     # change size of mupat
                mupat = [numpy.identity(self.n_orbitals[ik][isp], numpy.complex_) * mu for isp in range(self.Nspinblocs)]   # construct mupat
                # mupat *= mu
 
            #set a temporary array storing spectral functions with band index. Note, usually we should have spin index
                #Annkw=numpy.zeros((self.n_orbitals[ik],self.n_orbitals[ik],N_om),dtype=numpy.complex_)
                Annkw = [numpy.zeros((self.n_orbitals[ik][isp], self.n_orbitals[ik][isp], N_om), dtype=numpy.complex_) for isp in range(self.Nspinblocs)]

            # get lattice green functions.
            # S <<= A_Omega_Plus_B(A=1, B=1j * broadening)
             
	    S <<= 1*Omega + 1j*broadening

            Ms = copy.deepcopy(mupat)
            for ibl in range(self.Nspinblocs):
                ind = ntoi[bln[ibl]]
		n_orb = self.n_orbitals[ik][ibl]
		Ms[ibl] = self.hopping[ik,ind,0:n_orb,0:n_orb].real - mupat[ibl]
            S -= Ms
#            tmp = S.copy()    # init temporary storage
#            ## substract self energy
#            for icrsh in xrange(self.n_corr_shells):
#                for sig, gf in tmp: tmp[sig] <<= self.upfold(ik, icrsh, sig, stmp[icrsh][sig], gf)
#                S -= tmp

            S.invert()

            for isp in range(self.Nspinblocs):
                Annkw[isp].real = -copy.deepcopy(S[self.block_names[self.SO][isp]].data.swapaxes(0,1).swapaxes(1,2)).imag / numpy.pi

            for isp in range(self.Nspinblocs):
                if(ik%100==0):
                    print "ik,isp", ik, isp
                kvel = velocities[isp].vks[ik]

                Pwtem = numpy.zeros((mshape.sum(), len(Qmesh_ex), N_om), dtype=numpy.float_)

                #symmetry loop
                for Rmat in self.symm:
                    # get new velocity.
                    Rkvel = copy.deepcopy(kvel.vel)
                    for vnb1 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                        for vnb2 in xrange(kvel.bandwin[1] - kvel.bandwin[0] + 1):
                            Rkvel[vnb1][vnb2][:] = numpy.dot(Rmat, Rkvel[vnb1][vnb2][:])
                    ipw = 0
                    for (ir, ic) in mlist:
                        for iw in xrange(N_om):
                            if(M[iw].real > 5.0 / Beta):
                                continue
                            for iq in range(len(Qmesh_ex)):
                                #if(Qmesh_ex[iq]==0 or iw+Qmesh_ex[iq]>=N_om ):
                                # here use fermi distribution to truncate self energy mesh.
                                if(Qmesh_ex[iq] == 0 or iw + Qmesh_ex[iq] >= N_om  
                                   or M[iw].real + Qmesh[iq] < -10.0 / Beta or M[iw].real >10.0 / Beta):
                                    continue

                                if (M[iw].real > omminplot) and (M[iw].real < ommaxplot):
                                # here use bandwin to construct match matrix for A and velocity.
                                    bmin = max(self.bandwin[isp][ik, 0], kvel.bandwin[0])
                                    bmax = min(self.bandwin[isp][ik, 1], kvel.bandwin[1])
                                    Astart = bmin - self.bandwin[isp][ik, 0]
                                    Aend = bmax - self.bandwin[isp][ik, 0] + 1
                                    vstart = bmin - kvel.bandwin[0]
                                    vend = bmax - kvel.bandwin[0] + 1
                                    Annkwl = Annkw[isp][Astart:Aend, Astart:Aend, iw]
                                    Annkwr = Annkw[isp][Astart:Aend, Astart:Aend, iw + Qmesh_ex[iq]]
                                    Rkveltr = Rkvel[vstart:vend, vstart:vend, ir]
                                    Rkveltc = Rkvel[vstart:vend, vstart:vend, ic]
                                #print Annkwt.shape,Rkvel[...,ir].shape
                                    Pwtem[ipw, iq, iw] += numpy.dot(numpy.dot(numpy.dot(Rkveltr, Annkwl), Rkveltc), Annkwr).trace().real
                        ipw += 1

                # k sum and spin sum.
                self.Pw_optic += Pwtem * self.bz_weights[ik] / self.nsymm



        self.Pw_optic = myMPI.all_reduce(myMPI.world, self.Pw_optic, lambda x, y : x + y)
        self.Pw_optic *= (2 - self.SP)


        # just back up TD_optic data        # just back up TD_optic data
        if myMPI.is_master_node():
            with open("TD_Optic_LDA.dat", "w") as pwout:
                L1,L2,L3=self.Pw_optic.shape
                pwout.write("%s  %s  %s\n"%(L1,L2,L3))
                for i in xrange(L1):
                    for iq in xrange(L2):
                        for iw in xrange(L3):
                            pwout.write(str(i)+"   "+str(Qmesh_ex[iq] * deltaM) + "   " + str(M[iw]) + "   ")
                            pwout.write(str(self.Pw_optic[i, iq, iw]) + "    ")
                            pwout.write("\n")
                    pwout.write("\n")

        # sum over omega to get optic conductivity
        if myMPI.is_master_node():
            OpticConductivity = numpy.zeros((mshape.sum(), len(Qmesh)), dtype=numpy.float_)
            for im in range(mshape.sum()):
                for iq in range(len(Qmesh)):
                    for iw in xrange(N_om):
                        omegaT = M[iw] * Beta
                        omega_aug = Qmesh_ex[iq] * deltaM
                        OpticConductivity[im, iq] += self.Pw_optic[im, iq, iw] * (fermidis(omegaT) - fermidis(omegaT + omega_aug * Beta)) / omega_aug
            OpticConductivity *= deltaM
            OpticConductivity *= 10700 / self.Vol
            with open("Optic_con_LDA.dat", "w") as opt:
                for iq in range(len(Qmesh_ex)):
                    opt.write(str(Qmesh_ex[iq] * deltaM) + "   ")
                    for im in range(mshape.sum()):
                        opt.write(str(OpticConductivity[im, iq]) + "   ")
                    opt.write("\n")


    # load transportdistribution Pw from file.
    def loadTD(self, filename, fermishift=0.0):
        """ load transport distribution from file. Assume energy mesh is uniform
             the first column is energy mesh. The others are TD values.
             fermishift is used to shift the TD mesh to mimick a rigid band shift.
        """
        myMPI.barrier()
        self.TD = numpy.loadtxt(filename)
        self.TD[:, 0] -= fermishift


    # seebeck is just intetral of omega \int Pw f(\omega)f(-\omega)(\beta w).
    def Seebeck(self, Beta, index=0):
        """ get -A1/A0, that is Seebeck in unit k_B/e. index is used to select the column in self.tdintegral.
            Note: for nodiagonal element in Sxy this might not be right. so take care.
        """
        if myMPI.is_master_node():
            print "A0, A1  %.5e  %.5e " % (self.tdintegral(Beta, 0)[index], self.tdintegral(Beta, 1)[index])
            seb = -self.tdintegral(Beta, 1)[index] / self.tdintegral(Beta, 0)[index]
            print "Seebeck%d    %.4f  k_B/e        %.4f  x 10^(-6)V/K" % (index, seb, seb * 86.17)
            return seb

    def Conductivity(self, Beta, index=0):
        """ #return 1/T*A0, that is Conductivity in unit 1/V
            return Conductivity
        """
	if myMPI.is_master_node():
            Cond = Beta * self.tdintegral(Beta, 0)[index]
        #print "Beta*A0     ", Cond
            print "V in bohr^3          ", self.Vol
            Cond *=  10700.0 / self.Vol
            print "Conductivity%d       %.4f  x 10^4 Ohm^-1 cm^-1" % (index, Cond)
            print "Resistivity%d        %.4f  x 10^-4 Ohm cm" % (index, 1.0 / Cond)
            return Cond

    def tdintegral(self, Beta, pn=0):
        """calculate { \pi *\int Pw f(omega)f(-omega)(\beta\omega)^(pn)d\omega }
        """
        M = self.TD[:, 0]
        domega = abs(M[1] - M[0])
        pwint = numpy.zeros(self.TD.shape[1] - 1)
        for ipw in range(self.TD.shape[1] - 1):
            for iw in xrange(len(M)):
                x = M[iw] * Beta
                pwint[ipw] += fermidis(x) * fermidis(-x) * self.TD[iw, ipw + 1] * numpy.float(x) ** pn * domega

        return pwint


    def tdintcore(self, Beta):
        """calculate {  Pw f(omega)f(-omega) for check data)
        """
        M = self.TD[:, 0]
        domega = abs(M[1] - M[0])
        pwint = numpy.zeros((self.TD.shape[1] - 1, M.size))
        for ipw in range(self.TD.shape[1] - 1):
            for iw in xrange(M.size):
                x = M[iw] * Beta
                pwint[ipw, iw] = self.TD[iw, ipw + 1] * fermidis(x) * fermidis(-x)
                #if(pwint[ipw,iw]>=0.3):
                #    self.Pw[ipw,iw]=0.0

        if myMPI.is_master_node():
            with open("tdintcore.dat", "w") as pwout:
                for iw in xrange(M.size):
                    pwout.write(str(M[iw]) + "   ")
                    for i in range(pwint.shape[0]):
                        pwout.write(str(pwint[i, iw]) + "    ")
                    pwout.write("\n")
        return pwint


    def bandwinfromwiencase(self, wiencase):
        """ read in the band window from wiencase.outbwin file.
        """
        bandwin = [numpy.zeros(self.n_k * 2, dtype=int).reshape(self.n_k, 2) for isp in range(self.SP + 1 - self.SO)]
        for isp in range(self.SP + 1 - self.SO):
            if(self.SP == 0 or self.SO == 1):
                winfile = Read_Fortran_File2(wiencase + ".oubwin")
            elif self.SP == 1 and isp == 0:
                winfile = Read_Fortran_File2(wiencase + ".oubwinup")
            elif self.SP == 1 and isp == 1:
                winfile = Read_Fortran_File2(wiencase + ".oubwindn")
            else:
                assert 0, "Reading bandwin error! Check self.SP and self.SO!"
            Nk = int(winfile.next())
            assert Nk == self.n_k, "Number of K points is unconsistent in case.oubwin"
            SO = int(winfile.next())
            assert SO == self.SO, "SO is unconsistent in case.oubwin"

            for i in xrange(self.n_k):
                winfile.next()
                bandwin[isp][i, 0] = winfile.next()
                bandwin[isp][i, 1] = winfile.next()
                winfile.next()
        return bandwin
                 
