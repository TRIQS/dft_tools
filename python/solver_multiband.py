
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


from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb_matrix import Solver
from pytriqs.applications.dft.U_matrix import *

import pytriqs.utility.mpi as mpi
from types import *
import numpy

def sum_list(L):
    """ Can sum any list"""
    return reduce(lambda x, y: x+y, L) if len(L)>0 else []

#########################################
#
#  Solver for the Multi-Band problem
#
#########################################


class SolverMultiBand(Solver):
    """ 
    This is a general solver for a multiband local Hamiltonian. 
    Calling arguments: 
    beta = inverse temperature
    n_orb = Number of local orbitals
    U_interact = Average Coulomb interaction
    J_hund     = Hund coupling
    use_spinflip = true/false
    use_pairhop  = true/false
    use_matrix: Use the interaction matrix calculated from the Slater integrals
    is use_matrix, you need also:
        l: angular momentum of the orbital, l=2 is d
        T: Transformation matrix for U vertex. If not present, use standard complex harmonics
               
    """
    
    def __init__(self, beta, n_orb, gf_struct = False, map = False, omega_max = None):

        self.n_orb = n_orb

        # either get or construct gf_struct
        if (gf_struct):
            assert map, "give also the mapping!"
            self.map = map
        else:
            # standard gf_struct and map
            gf_struct = [ ('%s'%(ud),[n for n in range(n_orb)]) for ud in ['up','down'] ]
            self.map = {'up' : ['up' for v in range(n_orb)], 'down' : ['down' for v in range(n_orb)]}

        # now initialize the solver with the stuff given above:
        if (omega_max is None):
            Solver.__init__(self, beta = beta, gf_struct = gf_struct)
        else:
            n_w = int(omega_max*beta/(2.*numpy.pi))
            Solver.__init__(self, beta = beta, gf_struct = gf_struct, n_w = n_w)
            


    def solve(self, U_interact=None, J_hund=None, use_spinflip=False,
                 use_matrix = True, l=2, T=None, dim_reps=None, irep=None, deg_orbs = [], sl_int = None, **params):

        self.use_spinflip = use_spinflip
        self.U, self.Up, self.U4ind, self.offset = set_U_matrix(U_interact,J_hund,self.n_orb,l,use_matrix,T,sl_int,use_spinflip,dim_reps,irep) 

        # define mapping of indices:
        self.map_ind={}
        for nm in self.map:
            bl_names = self.map[nm]
            block = []
            for a,al in self.gf_struct:
                if a in bl_names: block.append(al)
                        
            self.map_ind[nm] = range(self.n_orb)
            i = 0
            for al in block:
                cnt = 0
                for a in range(len(al)):
                    self.map_ind[nm][i] = cnt
                    i = i+1
                    cnt = cnt+1
        
        # set the Hamiltonian
        if (use_spinflip==False):
            Hamiltonian = self.__set_hamiltonian_density()
        else:
            if (use_matrix):
                #Hamiltonian = self.__set_full_hamiltonian_slater()
                Hamiltonian = self.__set_spinflip_hamiltonian_slater()
            else:
                Hamiltonian = self.__set_full_hamiltonian_kanamori(J_hund = J_hund)

        # set the Quantum numbers
        Quantum_Numbers = self.__set_quantum_numbers(self.gf_struct)
    
        # Determine if there are only blocs of size 1:
        self.blocssizeone = True
        for ib in self.gf_struct:
            if (len(ib[1])>1): self.blocssizeone = False

        nc = params.pop("n_cycles",10000)
        if ((self.blocssizeone) and (self.use_spinflip==False)):
            use_seg = True
        else:
            use_seg = False
        #gm = self.set_global_moves(deg_orbs)

        Solver.solve(self,H_local = Hamiltonian, quantum_numbers = Quantum_Numbers, n_cycles = nc, use_segment_picture = use_seg, **params)


    def set_global_moves(self, deg_orbs, factor=0.05):
        # Sets some global moves given orbital degeneracies:
        
        strbl  = ''
        strind = ''
        inddone = []

        for orbs in deg_orbs:
            ln = len(orbs)
            orbsorted = sorted(orbs)
            for ii in range(ln):
                if (strbl!=''): strbl += ','
                bl1 = orbsorted[ii]
                bl2 = orbsorted[(ii+1)%ln]
                ind1 = [ll for ll in self.Sigma[bl1].indices ]
                ind2 = [ll for ll in self.Sigma[bl2].indices ]

                strbl += "'" + bl1 + "':'" + bl2 + "'"
                for kk, ind in enumerate(ind1):
                    if not (ind in inddone):
                        if (strind!=''): strind += ','
                        strind += '%s:%s'%(ind1[kk],ind2[kk])
                        inddone.append(ind)
                

        if len(deg_orbs)>0:
            str = 'Global_Moves = [ (%s, lambda (a,alpha,dag) : ({ '%factor + strbl + ' }[a], {' + strind + '}[alpha], dag) )]'
            exec str
            return Global_Moves
        else:
            return []
        

    def __set_hamiltonian_density(self):
        # density-density Hamiltonian:
        
        spinblocs = [v for v in self.map]
        #print spinblocs
        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it

        for sp1 in spinblocs:
            for sp2 in spinblocs:
                for i in range(self.n_orb):
                    for j in range(self.n_orb):
                        if (sp1==sp2):
                            Hamiltonian += 0.5 * self.U[self.offset+i,self.offset+j] * N(self.map[sp1][i],self.map_ind[sp1][i]) * N(self.map[sp2][j],self.map_ind[sp2][j]) 
                        else:
                            Hamiltonian += 0.5 * self.Up[self.offset+i,self.offset+j] * N(self.map[sp1][i],self.map_ind[sp1][i]) * N(self.map[sp2][j],self.map_ind[sp2][j]) 

        Hamiltonian -= N(self.map[spinblocs[0]][0],0)      # substract the initializing value

        return Hamiltonian


    def __set_full_hamiltonian_slater(self):
      
        spinblocs = [v for v in self.map]
        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it
        #print "Starting..."
        # use the full 4-index U-matrix:
        #for sp1 in spinblocs:
        #    for sp2 in spinblocs:
        for m1 in range(self.n_orb):
            for m2 in range(self.n_orb):
                for m3 in range(self.n_orb):
                    for m4 in range(self.n_orb):
                        if (abs(self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4])>0.00001):
                            for sp1 in spinblocs:
                                for sp2 in spinblocs:
                                    #print sp1,sp2,m1,m2,m3,m4
                                    Hamiltonian += 0.5 * self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4] * \
                                        Cdag(self.map[sp1][m1],self.map_ind[sp1][m1]) * Cdag(self.map[sp2][m2],self.map_ind[sp2][m2]) * C(self.map[sp2][m4],self.map_ind[sp2][m4]) * C(self.map[sp1][m3],self.map_ind[sp1][m3])
        #print "end..."
        Hamiltonian -= N(self.map[spinblocs[0]][0],0)      # substract the initializing value
                        
        return Hamiltonian


    def __set_spinflip_hamiltonian_slater(self):
        """Takes only spin-flip and pair-hopping terms"""
        
        spinblocs = [v for v in self.map]
        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it
        #print "Starting..."
        # use the full 4-index U-matrix:
        #for sp1 in spinblocs:
        #    for sp2 in spinblocs:
        for m1 in range(self.n_orb):
            for m2 in range(self.n_orb):
                for m3 in range(self.n_orb):
                    for m4 in range(self.n_orb):
                        if ((abs(self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4])>0.00001) and
                            ( ((m1==m2)and(m3==m4)) or ((m1==m3)and(m2==m4)) or ((m1==m4)and(m2==m3)) ) ):
                            for sp1 in spinblocs:
                                for sp2 in spinblocs:
                                    #print sp1,sp2,m1,m2,m3,m4
                                    Hamiltonian += 0.5 * self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4] * \
                                        Cdag(self.map[sp1][m1],self.map_ind[sp1][m1]) * Cdag(self.map[sp2][m2],self.map_ind[sp2][m2]) * C(self.map[sp2][m4],self.map_ind[sp2][m4]) * C(self.map[sp1][m3],self.map_ind[sp1][m3])
        #print "end..."
        Hamiltonian -= N(self.map[spinblocs[0]][0],0)      # substract the initializing value
                        
        return Hamiltonian


            
    def __set_full_hamiltonian_kanamori(self,J_hund):

        spinblocs = [v for v in self.map]
        assert len(spinblocs)==2,"spinflips in Kanamori representation only implemented for up/down structure!"

        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it

        # density terms:
        for sp1 in spinblocs:
            for sp2 in spinblocs:
                for i in range(self.n_orb):
                    for j in range(self.n_orb):
                        if (sp1==sp2):
                            Hamiltonian += 0.5 * self.U[self.offset+i,self.offset+j] * N(self.map[sp1][i],self.map_ind[sp1][i]) * N(self.map[sp2][j],self.map_ind[sp2][j]) 
                        else: 
                            Hamiltonian += 0.5 * self.Up[self.offset+i,self.offset+j] * N(self.map[sp1][i],self.map_ind[sp1][i]) * N(self.map[sp2][j],self.map_ind[sp2][j]) 

        # spinflip term:
        sp1 = spinblocs[0]
        sp2 = spinblocs[1]
        for i in range(self.n_orb-1):
            for j in range(i+1,self.n_orb):
                Hamiltonian -= J_hund * ( Cdag(self.map[sp1][i],self.map_ind[sp1][i]) * C(self.map[sp2][i],self.map_ind[sp2][i]) * Cdag(self.map[sp2][j],self.map_ind[sp2][j]) * C(self.map[sp1][j],self.map_ind[sp1][j]) )     # first term
                Hamiltonian -= J_hund * ( Cdag(self.map[sp2][i],self.map_ind[sp2][i]) * C(self.map[sp1][i],self.map_ind[sp1][i]) * Cdag(self.map[sp1][j],self.map_ind[sp1][j]) * C(self.map[sp2][j],self.map_ind[sp2][j]) )     # second term

        # pairhop terms:
        for i in range(self.n_orb-1):
            for j in range(i+1,self.n_orb):
                Hamiltonian -= J_hund * ( Cdag(self.map[sp1][i],self.map_ind[sp1][i]) * Cdag(self.map[sp2][i],self.map_ind[sp2][i]) * C(self.map[sp1][j],self.map_ind[sp1][j]) * C(self.map[sp2][j],self.map_ind[sp2][j]) )     # first term
                Hamiltonian -= J_hund * ( Cdag(self.map[sp2][j],self.map_ind[sp2][j]) * Cdag(self.map[sp1][j],self.map_ind[sp1][j]) * C(self.map[sp2][i],self.map_ind[sp2][i]) * C(self.map[sp1][i],self.map_ind[sp1][i]) )     # second term  

        Hamiltonian -= N(self.map[spinblocs[0]][0],0)       # substract the initializing value
                        
        return Hamiltonian
   

    def __set_quantum_numbers(self,gf_struct):
    
        QN = {}
        spinblocs = [v for v in self.map]

        # Define the quantum numbers:
        if (self.use_spinflip) :            
            Ntot = sum_list( [ N(self.map[s][i],self.map_ind[s][i]) for s in spinblocs for i in range(self.n_orb) ] )
            QN['NtotQN'] = Ntot
            #QN['Ntot'] = sum_list( [ N(self.map[s][i],i) for s in spinblocs for i in range(self.n_orb) ] )
            if (len(spinblocs)==2):
                # Assuming up/down structure:
                Sz = sum_list( [ N(self.map[spinblocs[0]][i],self.map_ind[spinblocs[0]][i])-N(self.map[spinblocs[1]][i],self.map_ind[spinblocs[1]][i]) for i in range(self.n_orb) ] )
                QN['SzQN'] = Sz
                # new quantum number: works only if there are only spin-flip and pair hopping, not any more complicated things
                for i in range(self.n_orb):
                    QN['Sz2_%s'%i] = (N(self.map[spinblocs[0]][i],self.map_ind[spinblocs[0]][i])-N(self.map[spinblocs[1]][i],self.map_ind[spinblocs[1]][i])) * (N(self.map[spinblocs[0]][i],self.map_ind[spinblocs[0]][i])-N(self.map[spinblocs[1]][i],self.map_ind[spinblocs[1]][i]))

        else :
            for ibl in range(len(gf_struct)):
                QN['N%s'%gf_struct[ibl][0]] = sum_list( [ N(gf_struct[ibl][0],gf_struct[ibl][1][i]) for i in range(len(gf_struct[ibl][1])) ] )

        return QN


    def fit_tails(self): 
	"""Fits the tails using the constant value for the Re Sigma calculated from F=Sigma*G.
           Works only for blocks of size one."""
	
	#if (len(self.gf_struct)==2*self.n_orb):
        if (self.blocssizeone):
            spinblocs = [v for v in self.map]
            mpi.report("Fitting tails manually")
	
            known_coeff = numpy.zeros([1,1,2],numpy.float_)
            msh = [x.imag for x in self.G[self.map[spinblocs[0]][0]].mesh ]
            fit_start = msh[self.fitting_Frequency_Start]
            fit_stop = msh[self.N_Frequencies_Accumulated]	
            
            # Fit the tail of G just to get the density
            for n,g in self.G:
                g.fitTail([[[0,0,1]]],7,fit_start,2*fit_stop) 
            densmat = self.G.density()

            for sig1 in spinblocs:
                for i in range(self.n_orb):

                    coeff = 0.0

                    for sig2 in spinblocs:
                        for j in range(self.n_orb):
                            if (sig1==sig2):
                                coeff += self.U[self.offset+i,self.offset+j] * densmat[self.map[sig1][j]][0,0].real
                            else:
                                coeff += self.Up[self.offset+i,self.offset+j] * densmat[self.map[sig2][j]][0,0].real

                    known_coeff[0,0,1] = coeff
                    self.Sigma[self.map[sig1][i]].fitTail(fixed_coef = known_coeff, order_max = 3, fit_start = fit_start, fit_stop = fit_stop)

        else:

            for n,sig in self.Sigma:

                known_coeff = numpy.zeros([sig.N1,sig.N2,1],numpy.float_)
                msh = [x.imag for x in sig.mesh]
                fit_start = msh[self.fitting_Frequency_Start]
                fit_stop  = msh[self.N_Frequencies_Accumulated]
            
                sig.fitTail(fixed_coef = known_coeff, order_max = 3, fit_start = fit_start, fit_stop = fit_stop)

		


class SolverMultiBandOld(SolverMultiBand):
    """ 
    Old MultiBand Solver construct
    """
    
    def __init__(self, Beta, Norb, U_interact=None, J_Hund=None, GFStruct=False, map=False, use_spinflip=False,
                 useMatrix = True, l=2, T=None, dimreps=None, irep=None, deg_orbs = [], Sl_Int = None):

        SolverMultiBand.__init__(self, beta=Beta, n_orb=Norb, gf_struct=GFStruct, map=map)
        self.U_interact = U_interact
        self.J_Hund = J_Hund
        self.use_spinflip = use_spinflip
        self.useMatrix = useMatrix
        self.l = l
        self.T = T
        self.dimreps = dimreps
        self.irep = irep
        self.deg_orbs = deg_orbs
        self.Sl_Int = Sl_Int
        self.gen_keys = copy.deepcopy(self.__dict__)

        msg = """
**********************************************************************************
 Warning: You are using the old constructor for the solver. Beware that this will
 be deprecated in future versions. Please check the documentation.
**********************************************************************************
"""
        mpi.report(msg)


    def Solve(self):

        params = copy.deepcopy(self.__dict__)
        for i in self.gen_keys: self.params.pop(i)
        self.params.pop("gen_keys")
        self.solve(self, U_interact=self.U_interact, J_hund=self.J_Hund, use_spinflip=self.use_spinflip,
                   use_matrix = self.useMatrix, l=self.l, T=self.T, dim_reps=self.dimreps, irep=self.irep,
                   deg_orbs = self.deg_orbs, sl_int = self.Sl_Int, **params)




	
def set_U_matrix(U_interact,J_hund,n_orb,l,use_matrix=True,T=None,sl_int=None,use_spinflip=False,dim_reps=None,irep=None):
    """ Set up the interaction vertex""" 

    offset = 0
    U4ind = None
    U = None
    Up = None
    if (use_matrix):
        if not (sl_int is None):
            Umat = Umatrix(l=l)
            assert len(sl_int)==(l+1),"sl_int has the wrong length"
            if (type(sl_int)==ListType):
                Rcl = numpy.array(sl_int)
            else:
                Rcl = sl_int
            Umat(T=T,rcl=Rcl)
        else:
            if ((U_interact==None)and(J_hund==None)):
                mpi.report("Give U,J or Slater integrals!!!")
                assert 0
            Umat = Umatrix(U_interact=U_interact, J_hund=J_hund, l=l)
            Umat(T=T)
            
        Umat.reduce_matrix()
        if (Umat.N==Umat.Nmat):
            # Transformation T is of size 2l+1
            U = Umat.U
            Up = Umat.Up
        else:
            # Transformation is of size 2(2l+1)
            U = Umat.U
         # now we have the reduced matrices U and Up, we need it for tail fitting anyways

        if (use_spinflip):
            #Take the 4index Umatrix
            # check for imaginary matrix elements:
            if (abs(Umat.Ufull.imag)>0.0001).any():
                mpi.report("WARNING: complex interaction matrix!! Ignoring imaginary part for the moment!")
                mpi.report("If you want to change this, look into Wien2k/solver_multiband.py")
            U4ind = Umat.Ufull.real
    
        # this will be changed for arbitrary irep:
        # use only one subgroup of orbitals?
        if not (irep is None):
            #print irep, dim_reps
            assert not (dim_reps is None), "Dimensions of the representatives are missing!"
            assert n_orb==dim_reps[irep-1],"Dimensions of dimrep and n_orb do not fit!"
            for ii in range(irep-1):
                offset += dim_reps[ii]
    else:
        if ((U_interact==None)and(J_hund==None)):
            mpi.report("For Kanamori representation, give U and J!!")
            assert 0
        U  = numpy.zeros([n_orb,n_orb],numpy.float_)
        Up = numpy.zeros([n_orb,n_orb],numpy.float_)
        for i in range(n_orb):
            for j in range(n_orb):
	        if (i==j):
	            Up[i,i] = U_interact 
	        else:
	       	    Up[i,j] = U_interact - 2.0*J_hund
		    U[i,j]  = U_interact - 3.0*J_hund

    return U, Up, U4ind, offset
