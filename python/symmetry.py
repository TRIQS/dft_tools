
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


import copy,numpy
import string
from types import *
from pytriqs.gf.local import *
from pytriqs.archive import *
import pytriqs.utility.mpi as mpi


class Symmetry:
    """This class provides the routines for applying symmetry operations for the k sums.
       It contains the permutations of the atoms in the unti cell, and the corresponding
       rotational matrices for each symmetry operation."""

    def __init__(self, hdf_file, subgroup = None):
        """Initialises the class.
           Reads the permutations and rotation matrizes from the file, and constructs the mapping for
           the given orbitals. For each orbit a matrix is read!!!
           SO: Flag for SO coupled calculations.
           SP: Spin polarisation yes/no
           """

        assert type(hdf_file)==StringType,"hdf_file must be a filename"; self.hdf_file = hdf_file
        things_to_read = ['n_s','n_atoms','perm','orbits','SO','SP','time_inv','mat','mat_tinv']
        for it in things_to_read: setattr(self,it,0)

        if (mpi.is_master_node()):
            #Read the stuff on master:
            ar = HDFArchive(hdf_file,'a')
            if (subgroup is None):
                ar2 = ar
            else:
                ar2 = ar[subgroup]

            for it in things_to_read: setattr(self,it,ar2[it])
            del ar2
            del ar

        # Broadcasting
        for it in things_to_read: setattr(self,it,mpi.bcast(getattr(self,it)))

        # now define the mapping of orbitals:
        # self.map[iorb]=jorb gives the permutation of the orbitals as given in the list, when the
        # permutation of the atoms is done:
        self.n_orbits = len(self.orbits)

        self.map = [ [0 for iorb in range(self.n_orbits)] for in_s in range(self.n_s) ]
        for in_s in range(self.n_s):
            for iorb in range(self.n_orbits):

                srch = copy.deepcopy(self.orbits[iorb])
                srch[0] = self.perm[in_s][self.orbits[iorb][0]-1]
                self.map[in_s][iorb] = self.orbits.index(srch)



    def symmetrize(self,obj):

        assert isinstance(obj,list),"obj has to be a list of objects!"
        assert len(obj)==self.n_orbits,"obj has to be a list of the same length as defined in the init"

        if (isinstance(obj[0],BlockGf)):
            symm_obj = [ obj[i].copy() for i in range(len(obj)) ]        # here the result is stored, it is a BlockGf!
            for iorb in range(self.n_orbits): symm_obj[iorb].zero()      # set to zero
        else:
            # if not a BlockGf, we assume it is a matrix (density matrix), has to be complex since self.mat is complex!
            #symm_obj = [ numpy.zeros([self.orbits[iorb][3],self.orbits[iorb][3]],numpy.complex_) for iorb in range(self.n_orbits) ]
            symm_obj = [ copy.deepcopy(obj[i]) for i in range(len(obj)) ]

            for iorb in range(self.n_orbits):
                if (type(symm_obj[iorb])==DictType):
                    for ii in symm_obj[iorb]: symm_obj[iorb][ii] *= 0.0
                else:
                    symm_obj[iorb] *= 0.0


        for in_s in range(self.n_s):

            for iorb in range(self.n_orbits):

                l = self.orbits[iorb][2]         # s, p, d, or f
                dim = self.orbits[iorb][3]
                jorb = self.map[in_s][iorb]


                if (isinstance(obj[0],BlockGf)):

                    tmp = obj[iorb].copy()
                    if (self.time_inv[in_s]): tmp << tmp.transpose()
                    for bname,gf in tmp: tmp[bname].from_L_G_R(self.mat[in_s][iorb],tmp[bname],self.mat[in_s][iorb].conjugate().transpose())
                    tmp *= 1.0/self.n_s
                    symm_obj[jorb] += tmp

                else:

                    if (type(obj[iorb])==DictType):

                        for ii in obj[iorb]:
                            if (self.time_inv[in_s]==0):
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[in_s][iorb],obj[iorb][ii]),
                                                                self.mat[in_s][iorb].conjugate().transpose()) / self.n_s
                            else:
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[in_s][iorb],obj[iorb][ii].conjugate()),
                                                                self.mat[in_s][iorb].conjugate().transpose()) / self.n_s



                    else:
                        if (self.time_inv[in_s]==0):
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[in_s][iorb],obj[iorb]),self.mat[in_s][iorb].conjugate().transpose()) / self.n_s
                        else:
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[in_s][iorb],obj[iorb].conjugate()),
                                                        self.mat[in_s][iorb].conjugate().transpose()) / self.n_s


# Markus: This does not what it is supposed to do, check how this should work (keep for now)
#        if ((self.SO==0) and (self.SP==0)):
#            # add time inv:
            #mpi.report("Add time inversion")
#            for iorb in range(self.n_orbits):
#                if (isinstance(symm_obj[0],BlockGf)):
#                    tmp = symm_obj[iorb].copy()
#                    tmp << tmp.transpose()
#                    for bname,gf in tmp: tmp[bname].from_L_G_R(self.mat_tinv[iorb],tmp[bname],self.mat_tinv[iorb].transpose().conjugate())
#                    symm_obj[iorb] += tmp
#                    symm_obj[iorb] /= 2.0
#
#                else:
#                    if (type(symm_obj[iorb])==DictType):
#                        for ii in symm_obj[iorb]:
#                            symm_obj[iorb][ii] += numpy.dot(numpy.dot(self.mat_tinv[iorb],symm_obj[iorb][ii].conjugate()),
#                                                            self.mat_tinv[iorb].transpose().conjugate())
#                            symm_obj[iorb][ii] /= 2.0
#                    else:
#                        symm_obj[iorb] += numpy.dot(numpy.dot(self.mat_tinv[iorb],symm_obj[iorb].conjugate()),
#                                                    self.mat_tinv[iorb].transpose().conjugate())
#                        symm_obj[iorb] /= 2.0


        return symm_obj




