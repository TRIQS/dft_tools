 
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
#
# DFT tools: Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# PLOVasp: Copyright (C) 2015 by O. E. Peil
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
from pytriqs.archive import *
from converter_tools import *
import os.path
try:
    import simplejson as json
except ImportError:
    import json

class VaspConverter(ConverterTools):
    """
    Conversion from VASP output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, filename, hdf_filename = None,
                       dft_subgrp = 'dft_input', symmcorr_subgrp = 'dft_symmcorr_input',
                       parproj_subgrp='dft_parproj_input', symmpar_subgrp='dft_symmpar_input',
                       bands_subgrp = 'dft_bands_input', misc_subgrp = 'dft_misc_input',
                       transp_subgrp = 'dft_transp_input', repacking = False):
        """
        Init of the class. Variable filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so on. 
        """

        assert type(filename)==StringType, "Please provide the DFT files' base name as a string."
        if hdf_filename is None: hdf_filename = filename+'.h5'
        self.hdf_file = hdf_filename
        self.basename = filename
        self.ctrl_file = filename+'.ctrl'
#        self.pmat_file = filename+'.pmat'
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.parproj_subgrp = parproj_subgrp
        self.symmpar_subgrp = symmpar_subgrp
        self.bands_subgrp = bands_subgrp
        self.misc_subgrp = misc_subgrp
        self.transp_subgrp = transp_subgrp

        # Checks if h5 file is there and repacks it if wanted:
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)


    def read_data(self, fh):
        """
        Generator for reading plain data.
        """
        for line in fh:
            line_ = line.strip()
            if not line or (line_ == '' or line_[0] == '#'):
                continue

            for val in map(float, line.split()):
                yield val

    def read_header_and_data(self, filename):
        """
        Opens a file and returns a JSON-header and the generator for the plain data.
        """
        fh = open(filename, 'rt')
        header = ""
        for line in fh:
            if not "#END" in line:
                header += line
            else:
                break

        f_gen = self.read_data(fh)

        return header, f_gen

    def convert_dft_input(self):
        """
        Reads the input files, and stores the data in the HDFfile
        """
        energy_unit = 1.0 # VASP interface always uses eV
        k_dep_projection = 1
# Symmetries are switched off for the moment
# TODO: implement symmetries
        symm_op = 0                                   # Use symmetry groups for the k-sum

        # Read and write only on the master node
        if not (mpi.is_master_node()): return
        mpi.report("Reading input from %s..."%self.ctrl_file)

        # R is a generator : each R.Next() will return the next number in the file
        jheader, rf = self.read_header_and_data(self.ctrl_file)
        print jheader
        ctrl_head = json.loads(jheader)

        ng = ctrl_head['ngroups']
        n_k = ctrl_head['nk']
# Note the difference in name conventions!
        SP = ctrl_head['ns'] - 1
        SO = ctrl_head['nc_flag']

        kpts = numpy.zeros((n_k, 3))
        bz_weights = numpy.zeros(n_k)
        try:
            for ik in xrange(n_k):
                kx, ky, kz = rf.next(), rf.next(), rf.next()
                kpts[ik, :] = kx, ky, kz
                bz_weights[ik] = rf.next()
        except StopIteration:
            raise "VaspConverter: error reading %s"%self.ctrl_file

#        if nc_flag:
## TODO: check this
#            n_spin_blocs = 1
#        else:
#            n_spin_blocs = ns
        n_spin_blocs = SP + 1 - SO

# Read PLO groups
# First, we read everything into a temporary data structure
# TODO: think about multiple shell groups and how to map them on h5 structures
        assert ng == 1, "Only one group is allowed at the moment"

        try:
            for ig in xrange(ng):
                gr_file = self.basename + '.pg%i'%(ig + 1)
                jheader, rf = self.read_header_and_data(gr_file)
                gr_head = json.loads(jheader)

                e_win = gr_head['ewindow']
                nb_max = gr_head['nb_max']
                p_shells = gr_head['shells']
                density_required = gr_head['nelect']
                charge_below = 0.0 # This is not defined in VASP interface

# Note that in the DftTools convention each site gives a separate correlated shell!
                n_corr_shells = sum([len(sh['ion_list']) for sh in p_shells])

                corr_shells = []
                shion_to_corr_shell = [[] for ish in xrange(len(p_shells))]
                icsh = 0
                for ish, sh in enumerate(p_shells):
                    ion_list = sh['ion_list']
                    for i, ion in enumerate(ion_list):
                        pars = {}
                        pars['atom'] = ion
# We set all sites inequivalent
                        pars['sort'] = sh['ion_sort'][i]
                        pars['l'] = sh['lorb']
                        pars['dim'] = sh['ndim']
                        pars['SO'] = SO
# TODO: check what 'irep' entry does (it seems to be very specific to dmftproj)
                        pars['irep'] = 0
                        corr_shells.append(pars)
                        shion_to_corr_shell[ish].append(i)

# TODO: generalize this to the case of multiple shell groups
                n_shells = n_corr_shells # No non-correlated shells at the moment
                shells = corr_shells

# FIXME: atomic sorts in Wien2K are not the same as in VASP.
#        A symmetry analysis from OUTCAR or symmetry file should be used
#        to define equivalence classes of sites.
            n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = ConverterTools.det_shell_equivalence(self, corr_shells)

            if mpi.is_master_node():
                print "  No. of inequivalent shells:", n_inequiv_shells

# NB!: these rotation matrices are specific to Wien2K! Set to identity in VASP
            use_rotations = 1
            rot_mat = [numpy.identity(corr_shells[icrsh]['dim'],numpy.complex_) for icrsh in range(n_corr_shells)]
            rot_mat_time_inv = [0 for i in range(n_corr_shells)]

# TODO: implement transformation matrices
            n_reps = [1 for i in range(n_inequiv_shells)]
            dim_reps = [0 for i in range(n_inequiv_shells)]
            T = []
            for ish in range(n_inequiv_shells):
                n_reps[ish] = 1   # Always 1 in VASP
                ineq_first = inequiv_to_corr[ish]
                dim_reps[ish] = [corr_shells[ineq_first]['dim']]   # Just the dimension of the shell
            
                # The transformation matrix:
                # is of dimension 2l+1 without SO, and 2*(2l+1) with SO!
                ll = 2 * corr_shells[inequiv_to_corr[ish]]['l']+1
                lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
# TODO: at the moment put T-matrices to identities
                T.append(numpy.identity(lmax, numpy.complex_))
                
#            if nc_flag:
## TODO: implement the noncollinear part
#                raise NotImplementedError("Noncollinear calculations are not implemented")
#            else:
            hopping = numpy.zeros([n_k, n_spin_blocs, nb_max, nb_max], numpy.complex_)
            f_weights = numpy.zeros([n_k, n_spin_blocs, nb_max], numpy.complex_)
            band_window = [numpy.zeros((n_k, 2), dtype=int) for isp in xrange(n_spin_blocs)]
            n_orbitals = numpy.zeros([n_k, n_spin_blocs], numpy.int)

            for isp in xrange(n_spin_blocs):
                for ik in xrange(n_k):
                    ib1, ib2 = int(rf.next()), int(rf.next())
                    band_window[isp][ik, :2] = ib1, ib2
                    nb = ib2 - ib1 + 1
                    n_orbitals[ik, isp] = nb
                    for ib in xrange(nb):
                        hopping[ik, isp, ib, ib] = rf.next()
                        f_weights[ik, isp, ib] = rf.next()

# Projectors
#            print n_orbitals
#            print [crsh['dim'] for crsh in corr_shells]
            proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max([crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], numpy.complex_)

# TODO: implement reading from more than one projector group
# In 'dmftproj' each ion represents a separate correlated shell.
# In my interface a 'projected shell' includes sets of ions.
# How to reconcile this? Two options:
#
# 1. Redefine 'projected shell' in my interface to make it correspond to one site only.
#    In this case the list of ions must be defined at the level of the projector group.
#
# 2. Split my 'projected shell' to several 'correlated shells' here in the converter.
#
# At the moment I choose i.2 for its simplicity. But one should consider possible
# use cases and decide which solution is to be made permanent.
#
            for ish, sh in enumerate(p_shells):
                for isp in xrange(n_spin_blocs):
                    for ik in xrange(n_k):
                        for ion in xrange(len(sh['ion_list'])):
                            icsh = shion_to_corr_shell[ish][ion]
                            for ilm in xrange(sh['ndim']):
                                for ib in xrange(n_orbitals[ik, isp]):
                                    # This is to avoid confusion with the order of arguments
                                    pr = rf.next()
                                    pi = rf.next()
                                    proj_mat[ik, isp, icsh, ilm, ib] = complex(pr, pi)

            things_to_set = ['n_shells','shells','n_corr_shells','corr_shells','n_spin_blocs','n_orbitals','n_k','SO','SP','energy_unit'] 
            for it in things_to_set:
#                print "%s:"%(it), locals()[it]
                setattr(self,it,locals()[it])

        except StopIteration:
           raise "VaspConverter: error reading %s"%self.gr_file

        rf.close()

        
        # Save it to the HDF:
        ar = HDFArchive(self.hdf_file,'a')
        if not (self.dft_subgrp in ar): ar.create_group(self.dft_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created. If it exists, the data is overwritten!
        things_to_save = ['energy_unit','n_k','k_dep_projection','SP','SO','charge_below','density_required',
                          'symm_op','n_shells','shells','n_corr_shells','corr_shells','use_rotations','rot_mat',
                          'rot_mat_time_inv','n_reps','dim_reps','T','n_orbitals','proj_mat','bz_weights','hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
        for it in things_to_save: ar[self.dft_subgrp][it] = locals()[it]

# Store Fermi weights to 'dft_misc_input'
        if not (self.misc_subgrp in ar): ar.create_group(self.misc_subgrp)
        ar[self.misc_subgrp]['dft_fermi_weights'] = f_weights
        ar[self.misc_subgrp]['band_window'] = band_window
        del ar
        # Symmetries are used, so now convert symmetry information for *correlated* orbitals:
        self.convert_symmetry_input(ctrl_head, orbits=self.corr_shells, symm_subgrp=self.symmcorr_subgrp)
# TODO: Implement misc_input
#        self.convert_misc_input(bandwin_file=self.bandwin_file,struct_file=self.struct_file,outputs_file=self.outputs_file,
#                                misc_subgrp=self.misc_subgrp,SO=self.SO,SP=self.SP,n_k=self.n_k)


    def convert_misc_input(self, bandwin_file, struct_file, outputs_file, misc_subgrp, SO, SP, n_k):
        """
        Reads input for the band window from bandwin_file, which is case.oubwin,
                            structure from struct_file, which is case.struct,
                            symmetries from outputs_file, which is case.outputs.
        """

        if not (mpi.is_master_node()): return
        things_to_save = []

        # Read relevant data from .oubwin/up/dn files
        #############################################
        # band_window: Contains the index of the lowest and highest band within the
        #              projected subspace (used by dmftproj) for each k-point.

        if (SP == 0 or SO == 1):
            files = [self.bandwin_file]
        elif SP == 1:
            files = [self.bandwin_file+'up', self.bandwin_file+'dn']
        else: # SO and SP can't both be 1
            assert 0, "convert_transport_input: Reding oubwin error! Check SP and SO!"
        
        band_window = [numpy.zeros((n_k, 2), dtype=int) for isp in range(SP + 1 - SO)]
        for isp, f in enumerate(files):
            if os.path.exists(f):
                mpi.report("Reading input from %s..."%f)
                R = ConverterTools.read_fortran_file(self, f, self.fortran_to_replace)
                assert int(R.next()) == n_k, "convert_misc_input: Number of k-points is inconsistent in oubwin file!"
                assert int(R.next()) == SO, "convert_misc_input: SO is inconsistent in oubwin file!"
                for ik in xrange(n_k):
                    R.next()
                    band_window[isp][ik,0] = R.next() # lowest band
                    band_window[isp][ik,1] = R.next() # highest band
                    R.next()
                things_to_save.append('band_window')

        R.close() # Reading done!

        # Read relevant data from .struct file
        ######################################
        # lattice_type: bravais lattice type as defined by Wien2k
        # lattice_constants: unit cell parameters in a. u.
        # lattice_angles: unit cell angles in rad

        if (os.path.exists(self.struct_file)):
            mpi.report("Reading input from %s..."%self.struct_file)
        
            with open(self.struct_file) as R:
                try:
                    R.readline()
                    lattice_type = R.readline().split()[0]
                    R.readline()
                    temp = R.readline()
#                    print temp
                    lattice_constants = numpy.array([float(temp[0+10*i:10+10*i].strip()) for i in range(3)])
                    lattice_angles = numpy.array([float(temp[30+10*i:40+10*i].strip()) for i in range(3)]) * numpy.pi / 180.0
                    things_to_save.extend(['lattice_type', 'lattice_constants', 'lattice_angles'])
                except IOError:
                    raise "convert_misc_input: reading file %s failed" %self.struct_file

        # Read relevant data from .outputs file
        #######################################
        # rot_symmetries: matrix representation of all (space group) symmetry operations
        
        if (os.path.exists(self.outputs_file)):
            mpi.report("Reading input from %s..."%self.outputs_file)
        
            rot_symmetries = []
            with open(self.outputs_file) as R:
                try:
                    while 1:
                        temp = R.readline().strip(' ').split()
                        if (temp[0] =='PGBSYM:'):
                            n_symmetries = int(temp[-1])
                            break
                    for i in range(n_symmetries):
                        while 1:
                            if (R.readline().strip().split()[0] == 'Symmetry'): break
                        sym_i = numpy.zeros((3, 3), dtype = float)
                        for ir in range(3):
                            temp = R.readline().strip().split()
                            for ic in range(3):
                                sym_i[ir, ic] = float(temp[ic])
                        R.readline()
                        rot_symmetries.append(sym_i)
                    things_to_save.extend(['n_symmetries', 'rot_symmetries'])
                    things_to_save.append('rot_symmetries')
                except IOError:
                    raise "convert_misc_input: reading file %s failed" %self.outputs_file

        # Save it to the HDF:
        ar=HDFArchive(self.hdf_file,'a')
        if not (misc_subgrp in ar): ar.create_group(misc_subgrp)
        for it in things_to_save: ar[misc_subgrp][it] = locals()[it]
        del ar


    def convert_symmetry_input(self, ctrl_head, orbits, symm_subgrp):
        """
        Reads input for the symmetrisations from symm_file, which is case.sympar or case.symqmc.
        """

# In VASP interface the symmetries are read directly from *.ctrl file
# For the moment the symmetry parameters are just stubs
        n_symm = 0
        n_atoms = 1
        perm = [0]
        n_orbits = len(orbits)
        SP = ctrl_head['ns']
        SO = ctrl_head['nc_flag']
        time_inv = [0]
        mat = [numpy.identity(1)]
        mat_tinv = [numpy.identity(1)]

        # Save it to the HDF:
        ar=HDFArchive(self.hdf_file,'a')
        if not (symm_subgrp in ar): ar.create_group(symm_subgrp)
        things_to_save = ['n_symm','n_atoms','perm','orbits','SO','SP','time_inv','mat','mat_tinv']
        for it in things_to_save:
#            print "%s:"%(it), locals()[it]
            ar[symm_subgrp][it] = locals()[it]
        del ar
