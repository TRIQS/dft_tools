
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2018 by G. J. Kraberger
# Copyright (C) 2018 by Simons Foundation
# Authors: G. J. Kraberger, O. Parcollet
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


import copy
import numpy as np
from pytriqs.gf import GfImFreq, BlockGf
from ast import literal_eval
import pytriqs.utility.mpi as mpi
from warnings import warn


class BlockStructure(object):
    """ Contains information about the Green function structure.

    This class contains information about the structure of the solver
    and sumk Green functions and the mapping between them.

    Parameters
    ----------
    gf_struct_sumk : list of list of tuple
        gf_struct_sumk[ish][idx] = (block_name,list of indices in block)

        for correlated shell ish; idx is just a counter in the list
    gf_struct_solver : list of dict
        gf_struct_solver[ish][block] = list of indices in that block

        for *inequivalent* correlated shell ish
    solver_to_sumk : list of dict
        solver_to_sumk[ish][(from_block,from_idx)] = (to_block,to_idx)

        maps from the solver block and index to the sumk block and index
        for *inequivalent* correlated shell ish
    sumk_to_solver : list of dict
        sumk_to_solver[ish][(from_block,from_idx)] = (to_block,to_idx)

        maps from the sumk block and index to the solver block and index
        for *inequivalent* correlated shell ish
    solver_to_sumk_block : list of dict
        solver_to_sumk_block[ish][from_block] = to_block

        maps from the solver block to the sumk block
        for *inequivalent* correlated shell ish
    """

    def __init__(self, gf_struct_sumk=None,
                 gf_struct_solver=None,
                 solver_to_sumk=None,
                 sumk_to_solver=None,
                 solver_to_sumk_block=None,
                 deg_shells=None):
        self.gf_struct_sumk = gf_struct_sumk
        self.gf_struct_solver = gf_struct_solver
        self.solver_to_sumk = solver_to_sumk
        self.sumk_to_solver = sumk_to_solver
        self.solver_to_sumk_block = solver_to_sumk_block
        self.deg_shells = deg_shells

    @property
    def gf_struct_solver_list(self):
        """ The structure of the solver Green's function

        This is returned as a
            list (for each shell)
        of lists (for each block)
        of tuples (block_name, block_indices).

        That is,
        ``gf_struct_solver_list[ish][b][0]``
        is the name of the block number ``b`` of shell ``ish``, and
        ``gf_struct_solver_list[ish][b][1]``
        is a list of its indices.

        The list for each shell is sorted alphabetically by block name.
        """
        # we sort by block name in order to get a reproducible result
        return [sorted([(k, v) for k, v in gfs.iteritems()], key=lambda x: x[0])
                for gfs in self.gf_struct_solver]

    @property
    def gf_struct_sumk_list(self):
        """ The structure of the sumk Green's function

        This is returned as a
            list (for each shell)
        of lists (for each block)
        of tuples (block_name, block_indices)

        That is,
        ``gf_struct_sumk_list[ish][b][0]``
        is the name of the block number ``b`` of shell ``ish``, and
        ``gf_struct_sumk_list[ish][b][1]``
        is a list of its indices.
        """
        return self.gf_struct_sumk

    @property
    def gf_struct_solver_dict(self):
        """ The structure of the solver Green's function

        This is returned as a
            list (for each shell)
        of dictionaries.

        That is,
        ``gf_struct_solver_dict[ish][bname]``
        is a list of the indices of block ``bname`` of shell ``ish``.
        """
        return self.gf_struct_solver

    @property
    def gf_struct_sumk_dict(self):
        """ The structure of the sumk Green's function

        This is returned as a
            list (for each shell)
        of dictionaries.

        That is,
        ``gf_struct_sumk_dict[ish][bname]``
        is a list of the indices of block ``bname`` of shell ``ish``.
        """
        return [{block: indices for block, indices in gfs}
                for gfs in self.gf_struct_sumk]

    @classmethod
    def full_structure(cls,gf_struct,corr_to_inequiv):
        """ Construct structure that maps to itself.

        This has the same structure for sumk and solver, and the
        mapping solver_to_sumk and sumk_to_solver is one-to-one.

        Parameters
        ----------
        gf_struct : list of dict
            gf_struct[ish][block] = list of indices in that block

            for (inequivalent) correlated shell ish
        corr_to_inequiv : list
            gives the mapping from correlated shell csh to inequivalent
            correlated shell icsh, so that corr_to_inequiv[csh]=icsh
            e.g. SumkDFT.corr_to_inequiv

            if None, each inequivalent correlated shell is supposed to
            be correspond to just one correlated shell with the same
            index; there is not default, None has to be set explicitly!
        """

        solver_to_sumk = []
        s2sblock = []
        gs_sumk = []
        for ish in range(len(gf_struct)):
            so2su = {}
            so2sublock = {}
            gss = []
            for block in gf_struct[ish]:
                so2sublock[block]=block
                for ind in gf_struct[ish][block]:
                    so2su[(block,ind)]=(block,ind)
                gss.append((block,gf_struct[ish][block]))
            solver_to_sumk.append(so2su)
            s2sblock.append(so2sublock)
            gs_sumk.append(gss)

        # gf_struct_sumk is not given for each inequivalent correlated
        # shell, but for every correlated shell!
        if corr_to_inequiv is not None:
            gs_sumk_all = [None]*len(corr_to_inequiv)
            for i in range(len(corr_to_inequiv)):
                gs_sumk_all[i] = gs_sumk[corr_to_inequiv[i]]
        else:
            gs_sumk_all = gs_sumk

        return cls(gf_struct_solver=copy.deepcopy(gf_struct),
                gf_struct_sumk = gs_sumk_all,
                solver_to_sumk = copy.deepcopy(solver_to_sumk),
                sumk_to_solver = solver_to_sumk,
                solver_to_sumk_block = s2sblock,
                deg_shells = [[] for ish in range(len(gf_struct))])

    def pick_gf_struct_solver(self,new_gf_struct):
        """ Pick selected orbitals within blocks.

        This throws away parts of the Green's function that (for some
        reason - be sure that you know what you're doing) shouldn't be
        included in the calculation.

        To drop an entire block, just don't include it.
        To drop a certain index within a block, just don't include it.

        If it was before:

        [{'up':[0,1],'down':[0,1],'left':[0,1]}]

        to choose the 0th index of the up block and the 1st index of
        the down block and drop the left block, the new_gf_struct would
        have to be

        [{'up':[0],'down':[1]}]

        Note that the indices will be renamed to be a 0-based
        sequence of integers, i.e. the new structure will actually
        be  [{'up':[0],'down':[0]}].

        For dropped indices, sumk_to_solver will map to (None,None).

        Parameters
        ----------
        new_gf_struct : list of dict
            formatted the same as gf_struct_solver:

            new_gf_struct[ish][block]=list of indices in that block.
        """

        for ish in range(len(self.gf_struct_solver)):
            gf_struct = new_gf_struct[ish]

            # create new solver_to_sumk
            so2su={}
            so2su_block = {}
            for blk,idxs in gf_struct.items():
                for i in range(len(idxs)):
                    so2su[(blk,i)]=self.solver_to_sumk[ish][(blk,idxs[i])]
                    so2su_block[blk]=so2su[(blk,i)][0]
            self.solver_to_sumk[ish] = so2su
            self.solver_to_sumk_block[ish] = so2su_block
            # create new sumk_to_solver
            for k,v in self.sumk_to_solver[ish].items():
                blk,ind=v
                if blk in gf_struct and ind in gf_struct[blk]:
                    new_ind = gf_struct[blk].index(ind)
                    self.sumk_to_solver[ish][k]=(blk,new_ind)
                else:
                    self.sumk_to_solver[ish][k]=(None,None)
            # reindexing gf_struct so that it starts with 0
            for k in gf_struct:
                gf_struct[k]=range(len(gf_struct[k]))
            self.gf_struct_solver[ish]=gf_struct

    def pick_gf_struct_sumk(self,new_gf_struct):
        """ Pick selected orbitals within blocks.

        This throws away parts of the Green's function that (for some
        reason - be sure that you know what you're doing) shouldn't be
        included in the calculation.

        To drop an entire block, just don't include it.
        To drop a certain index within a block, just don't include it.

        If it was before:

        [{'up':[0,1],'down':[0,1],'left':[0,1]}]

        to choose the 0th index of the up block and the 1st index of
        the down block and drop the left block, the new_gf_struct would
        have to be

        [{'up':[0],'down':[1]}]

        Note that the indices will be renamed to be a 0-based
        sequence of integers.

        For dropped indices, sumk_to_solver will map to (None,None).

        Parameters
        ----------
        new_gf_struct : list of dict
            formatted the same as gf_struct_solver:

            new_gf_struct[ish][block]=list of indices in that block.

            However, the indices are not according to the solver Gf
            but the sumk Gf.
        """


        gfs = []
        # construct gfs, which is the equivalent of new_gf_struct
        # but according to the solver Gf, by using the sumk_to_solver
        # mapping
        for ish in range(len(new_gf_struct)):
            gfs.append({})
            for block in new_gf_struct[ish].keys():
                for ind in new_gf_struct[ish][block]:
                    ind_sol = self.sumk_to_solver[ish][(block,ind)]
                    if not ind_sol[0] in gfs[ish]:
                        gfs[ish][ind_sol[0]]=[]
                    gfs[ish][ind_sol[0]].append(ind_sol[1])
        self.pick_gf_struct_solver(gfs)


    def map_gf_struct_solver(self,mapping):
        """ Map the Green function structure from one struct to another.

        Parameters
        ----------
        mapping : list of dict
            the dict consists of elements
            (from_block,from_index) : (to_block,to_index)
            that maps from one structure to the other
        """

        for ish in range(len(mapping)):
            gf_struct = {}
            so2su = {}
            su2so = {}
            so2su_block = {}
            for frm,to in mapping[ish].iteritems():
                if not to[0] in gf_struct:
                    gf_struct[to[0]]=[]
                gf_struct[to[0]].append(to[1])

                so2su[to]=self.solver_to_sumk[ish][frm]
                su2so[self.solver_to_sumk[ish][frm]]=to
                if to[0] in so2su_block:
                    if so2su_block[to[0]] != \
                        self.solver_to_sumk_block[ish][frm[0]]:
                            warn("solver block '{}' maps to more than one sumk block: '{}', '{}'".format(
                                    to[0],so2su_block[to[0]],self.solver_to_sumk_block[ish][frm[0]]))
                else:
                    so2su_block[to[0]]=\
                        self.solver_to_sumk_block[ish][frm[0]]
            for k in self.sumk_to_solver[ish].keys():
                if not k in su2so:
                    su2so[k] = (None,None)
            self.gf_struct_solver[ish]=gf_struct
            self.solver_to_sumk[ish]=so2su
            self.sumk_to_solver[ish]=su2so
            self.solver_to_sumk_block[ish]=so2su_block

    def create_gf(self,ish=0,gf_function=GfImFreq,**kwargs):
        """ Create a zero BlockGf having the gf_struct_solver structure.

        When using GfImFreq as gf_function, typically you have to
        supply beta as keyword argument.

        Parameters
        ----------
        ish : int
            shell index
        gf_function : constructor
            function used to construct the Gf objects constituting the
            individual blocks; default: GfImFreq
        **kwargs :
            options passed on to the Gf constructor for the individual
            blocks
        """

        names = self.gf_struct_solver[ish].keys()
        blocks=[]
        for n in names:
            G = gf_function(indices=self.gf_struct_solver[ish][n],**kwargs)
            blocks.append(G)
        G = BlockGf(name_list = names, block_list = blocks)
        return G

    def check_gf(self, G, ish=None):
        """ check whether the Green's function G has the right structure

        This throws an error if the structure of G is not the same
        as ``gf_struct_solver``.

        Parameters
        ----------
        G : BlockGf or list of BlockGf
            Green's function to check
            if it is a list, there should be as many entries as there
            are shells, and the check is performed for all shells (unless
            ish is given).
        ish : int
            shell index
            default: 0 if G is just one Green's function is given,
            check all if list of Green's functions is given
        """

        if isinstance(G, list):
            assert len(G) == len(self.gf_struct_solver),\
                "list of G does not have the correct length"
            if ish is None:
                ishs = range(len(self.gf_struct_solver))
            else:
                ishs = [ish]
            for ish in ishs:
                self.check_gf(G[ish], ish=ish)
            return

        if ish is None:
            ish = 0

        for block in self.gf_struct_solver[ish]:
            assert block in G.indices,\
                "block " + block + " not in G (shell {})".format(ish)
        for block, gf in G:
            assert block in self.gf_struct_solver[ish],\
                "block " + block + " not in struct (shell {})".format(ish)
            assert list(gf.indices) == 2 * [map(str, self.gf_struct_solver[ish][block])],\
                "block " + block + " has wrong indices (shell {})".format(ish)

    def convert_gf(self, G, G_struct, ish=0, show_warnings=True,
                   G_out=None, **kwargs):
        """ Convert BlockGf from its structure to this (solver) structure.

        .. warning::

            Elements that are zero in the new structure due to
            the new block structure will be just ignored, thus
            approximated to zero.

        Parameters
        ----------
        G : BlockGf
            the Gf that should be converted
        G_struct : BlockStructure or str
            the structure of that G or 'sumk' (then, the structure of
            the sumk Green's function of this BlockStructure is used)
        ish : int
            shell index
        show_warnings : bool or float
            whether to show warnings when elements of the Green's
            function get thrown away
            if float, set the threshold for the magnitude of an element
            about to be thrown away to trigger a warning
            (default: 1.e-10)
        G_out : BlockGf
            the output Green's function (if not given, a new one is
            created)
        **kwargs :
            options passed to the constructor for the new Gf
        """

        warning_threshold = 1.e-10
        if isinstance(show_warnings, float):
            warning_threshold = show_warnings
            show_warnings = True

        # we offer the possibility to convert to convert from sumk_dft
        from_sumk = False
        if isinstance(G_struct, str) and G_struct == 'sumk':
            gf_struct_in = self.gf_struct_sumk_dict[ish]
            from_sumk = True
        else:
            gf_struct_in = G_struct.gf_struct_solver[ish]

        # create a Green's function to hold the result
        if G_out is None:
            G_out = self.create_gf(ish=ish, **kwargs)
        else:
            self.check_gf(G_out, ish=ish)

        for block in gf_struct_in.keys():
            for i1 in gf_struct_in[block]:
                for i2 in gf_struct_in[block]:
                    if from_sumk:
                        i1_sumk = (block, i1)
                        i2_sumk = (block, i2)
                    else:
                        i1_sumk = G_struct.solver_to_sumk[ish][(block, i1)]
                        i2_sumk = G_struct.solver_to_sumk[ish][(block, i2)]

                    i1_sol = self.sumk_to_solver[ish][i1_sumk]
                    i2_sol = self.sumk_to_solver[ish][i2_sumk]
                    if i1_sol[0] is None or i2_sol[0] is None:
                        if show_warnings:
                            if mpi.is_master_node():
                                warn(('Element {},{} of block {} of G is not present ' +
                                      'in the new structure').format(i1, i2, block))
                        continue
                    if i1_sol[0] != i2_sol[0]:
                        if (show_warnings and
                                np.max(np.abs(G[block][i1, i2].data)) > warning_threshold):
                            if mpi.is_master_node():
                                warn(('Element {},{} of block {} of G is approximated ' +
                                      'to zero to match the new structure. Max abs value: {}').format(
                                    i1, i2, block, np.max(np.abs(G[block][i1, i2].data))))
                        continue
                    G_out[i1_sol[0]][i1_sol[1], i2_sol[1]] = \
                        G[block][i1, i2]
        return G_out

    def approximate_as_diagonal(self):
        """ Create a structure for a GF with zero off-diagonal elements.

        .. warning::

            In general, this will throw away non-zero elements of the
            Green's function. Be sure to verify whether this approximation
            is justified.
        """

        self.gf_struct_solver=[]
        self.solver_to_sumk=[]
        self.solver_to_sumk_block=[]
        for ish in range(len(self.sumk_to_solver)):
            self.gf_struct_solver.append({})
            self.solver_to_sumk.append({})
            self.solver_to_sumk_block.append({})
            for frm,to in self.sumk_to_solver[ish].iteritems():
                if to[0] is not None:
                    self.gf_struct_solver[ish][frm[0]+'_'+str(frm[1])]=[0]
                    self.sumk_to_solver[ish][frm]=(frm[0]+'_'+str(frm[1]),0)
                    self.solver_to_sumk[ish][(frm[0]+'_'+str(frm[1]),0)]=frm
                    self.solver_to_sumk_block[ish][frm[0]+'_'+str(frm[1])]=frm[0]

    def __eq__(self,other):
        def compare(one,two):
            if type(one)!=type(two):
                if not (isinstance(one, (bool, np.bool_)) and isinstance(two, (bool, np.bool_))):
                    return False
            if one is None and two is None:
                return True
            if isinstance(one,list) or isinstance(one,tuple):
                if len(one) != len(two):
                    return False
                for x,y in zip(one,two):
                    if not compare(x,y):
                        return False
                return True
            elif isinstance(one,(int,bool, str, np.bool_)):
                return one==two
            elif isinstance(one,np.ndarray):
                return np.all(one==two)
            elif isinstance(one,dict):
                if set(one.keys()) != set(two.keys()):
                    return False
                for k in set(one.keys()).intersection(two.keys()):
                    if not compare(one[k],two[k]):
                        return False
                return True
            warn('Cannot compare {}'.format(type(one)))
            return False

        for prop in [ "gf_struct_sumk", "gf_struct_solver",
                "solver_to_sumk", "sumk_to_solver", "solver_to_sumk_block",
                "deg_shells"]:
            if not compare(getattr(self,prop),getattr(other,prop)):
                return False
        return True

    def copy(self):
        return copy.deepcopy(self)

    def __reduce_to_dict__(self):
        """ Reduce to dict for HDF5 export."""

        ret = {}
        for element in [ "gf_struct_sumk", "gf_struct_solver",
                         "solver_to_sumk_block","deg_shells"]:
            ret[element] = getattr(self,element)

        def construct_mapping(mapping):
            d = []
            for ish in range(len(mapping)):
                d.append({})
                for k,v in mapping[ish].iteritems():
                    d[ish][repr(k)] = repr(v)
            return d

        ret['solver_to_sumk']=construct_mapping(self.solver_to_sumk)
        ret['sumk_to_solver']=construct_mapping(self.sumk_to_solver)
        return ret

    @classmethod
    def __factory_from_dict__(cls,name,D) :
        """ Create from dict for HDF5 import."""

        def reconstruct_mapping(mapping):
            d = []
            for ish in range(len(mapping)):
                d.append({})
                for k,v in mapping[ish].iteritems():
                    # literal_eval is a saje alternative to eval
                    d[ish][literal_eval(k)] = literal_eval(v)
            return d

        D['solver_to_sumk']=reconstruct_mapping(D['solver_to_sumk'])
        D['sumk_to_solver']=reconstruct_mapping(D['sumk_to_solver'])
        return cls(**D)

    def __str__(self):
        s=''
        s+= "gf_struct_sumk "+str( self.gf_struct_sumk)+'\n'
        s+= "gf_struct_solver "+str(self.gf_struct_solver)+'\n'
        s+= "solver_to_sumk_block "+str(self.solver_to_sumk_block)+'\n'
        for el in ['solver_to_sumk','sumk_to_solver']:
            s+=el+'\n'
            element=getattr(self,el)
            for ish in range(len(element)):
                s+=' shell '+str(ish)+'\n'
                def keyfun(el):
                    return '{}_{:05d}'.format(el[0],el[1])
                keys = sorted(element[ish].keys(),key=keyfun)
                for k in keys:
                    s+='  '+str(k)+str(element[ish][k])+'\n'
        s += "deg_shells\n"
        for ish in range(len(self.deg_shells)):
            s+=' shell '+str(ish)+'\n'
            for l in range(len(self.deg_shells[ish])):
                s+='  equivalent group '+str(l)+'\n'
                if isinstance(self.deg_shells[ish][l],dict):
                    for key, val in self.deg_shells[ish][l].iteritems():
                        s+='   '+key+('*' if val[1] else '')+':\n'
                        s+='    '+str(val[0]).replace('\n','\n    ')+'\n'
                else:
                    for key in self.deg_shells[ish][l]:
                        s+='   '+key+'\n'
        return s

from pytriqs.archive.hdf_archive_schemes import register_class
register_class(BlockStructure)
