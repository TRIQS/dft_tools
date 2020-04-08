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
    def __init__(self,gf_struct_sumk=None,
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
            for blk,idxs in list(gf_struct.items()):
                for i in range(len(idxs)):
                    so2su[(blk,i)]=self.solver_to_sumk[ish][(blk,idxs[i])]
                    so2su_block[blk]=so2su[(blk,i)][0]
            self.solver_to_sumk[ish] = so2su
            self.solver_to_sumk_block[ish] = so2su_block
            # create new sumk_to_solver
            for k,v in list(self.sumk_to_solver[ish].items()):
                blk,ind=v
                if blk in gf_struct and ind in gf_struct[blk]:
                    new_ind = gf_struct[blk].index(ind)
                    self.sumk_to_solver[ish][k]=(blk,new_ind)
                else:
                    self.sumk_to_solver[ish][k]=(None,None)
            # reindexing gf_struct so that it starts with 0
            for k in gf_struct:
                gf_struct[k]=list(range(len(gf_struct[k])))
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
            for block in list(new_gf_struct[ish].keys()):
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
            for frm,to in mapping[ish].items():
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
            for k in list(self.sumk_to_solver[ish].keys()):
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

        names = list(self.gf_struct_solver[ish].keys())
        blocks=[]
        for n in names:
            G = gf_function(indices=self.gf_struct_solver[ish][n],**kwargs)
            blocks.append(G)
        G = BlockGf(name_list = names, block_list = blocks)
        return G


    def convert_gf(self,G,G_struct,ish=0,show_warnings=True,**kwargs):
        """ Convert BlockGf from its structure to this structure.

        .. warning::

            Elements that are zero in the new structure due to
            the new block structure will be just ignored, thus
            approximated to zero.

        Parameters
        ----------
        G : BlockGf
            the Gf that should be converted
        G_struct : GfStructure
            the structure of that G
        ish : int
            shell index
        show_warnings : bool or float
            whether to show warnings when elements of the Green's
            function get thrown away
            if float, set the threshold for the magnitude of an element
            about to be thrown away to trigger a warning
            (default: 1.e-10)
        **kwargs :
            options passed to the constructor for the new Gf
        """

        warning_threshold = 1.e-10
        if isinstance(show_warnings, float):
            warning_threshold = show_warnings
            show_warnings = True

        G_new = self.create_gf(ish=ish,**kwargs)
        for block in list(G_struct.gf_struct_solver[ish].keys()):
            for i1 in G_struct.gf_struct_solver[ish][block]:
                for i2 in G_struct.gf_struct_solver[ish][block]:
                    i1_sumk = G_struct.solver_to_sumk[ish][(block,i1)]
                    i2_sumk = G_struct.solver_to_sumk[ish][(block,i2)]
                    i1_sol = self.sumk_to_solver[ish][i1_sumk]
                    i2_sol = self.sumk_to_solver[ish][i2_sumk]
                    if i1_sol[0] is None or i2_sol[0] is None:
                        if show_warnings:
                            if mpi.is_master_node():
                                warn(('Element {},{} of block {} of G is not present '+
                                    'in the new structure').format(i1,i2,block))
                        continue
                    if i1_sol[0]!=i2_sol[0]:
                        if show_warnings and np.max(np.abs(G[block][i1,i2].data)) > warning_threshold:
                            if mpi.is_master_node():
                                warn(('Element {},{} of block {} of G is approximated '+
                                    'to zero to match the new structure. Max abs value: {}').format(
                                        i1,i2,block,np.max(np.abs(G[block][i1,i2].data))))
                        continue
                    G_new[i1_sol[0]][i1_sol[1],i2_sol[1]] = \
                            G[block][i1,i2]
        return G_new

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
            for frm,to in self.sumk_to_solver[ish].items():
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
                for k in set(one.keys()).intersection(list(two.keys())):
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
                for k,v in mapping[ish].items():
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
                for k,v in mapping[ish].items():
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
                keys = sorted(list(element[ish].keys()),key=keyfun)
                for k in keys:
                    s+='  '+str(k)+str(element[ish][k])+'\n'
        s += "deg_shells\n"
        for ish in range(len(self.deg_shells)):
            s+=' shell '+str(ish)+'\n'
            for l in range(len(self.deg_shells[ish])):
                s+='  equivalent group '+str(l)+'\n'
                if isinstance(self.deg_shells[ish][l],dict):
                    for key, val in self.deg_shells[ish][l].items():
                        s+='   '+key+('*' if val[1] else '')+':\n'
                        s+='    '+str(val[0]).replace('\n','\n    ')+'\n'
                else:
                    for key in self.deg_shells[ish][l]:
                        s+='   '+key+'\n'
        return s

from h5.formats import register_class
register_class(BlockStructure)
