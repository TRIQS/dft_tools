
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
from h5 import *
from .converter_tools import *
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
                       transp_subgrp = 'dft_transp_input', repacking = False,
                       proj_or_hk='proj'):
        """
        Init of the class. Variable filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so on.

        Parameters
        ----------
        filename : string
                   Base name of DFT files.
        hdf_filename : string, optional
                       Name of hdf5 archive to be created.
        dft_subgrp : string, optional
                     Name of subgroup storing necessary DFT data.
        symmcorr_subgrp : string, optional
                          Name of subgroup storing correlated-shell symmetry data.
        parproj_subgrp : string, optional
                         Name of subgroup storing partial projector data.
        symmpar_subgrp : string, optional
                         Name of subgroup storing partial-projector symmetry data.
        bands_subgrp : string, optional
                       Name of subgroup storing band data.
        misc_subgrp : string, optional
                      Name of subgroup storing miscellaneous DFT data.
        transp_subgrp : string, optional
                        Name of subgroup storing transport data.
        repacking : boolean, optional
                    Does the hdf5 archive need to be repacked to save space?
        proj_or_hk : string, optional
                    Select scheme to convert between KS bands and localized orbitals.

        """

        assert isinstance(filename, str), "Please provide the DFT files' base name as a string."
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
        assert (proj_or_hk == 'proj') or (proj_or_hk == 'hk'), "proj_or_hk has to be 'proj' of 'hk'"
        self.proj_or_hk = proj_or_hk

        # Checks if h5 file is there and repacks it if wanted:
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)

        # this is to test pull request
    def read_data(self, fh):
        """
        Generator for reading plain data.

        Parameters
        ----------
        fh : file object
              file object which is read in.
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

        Parameters
        ----------
        filename : string
              file name of the file to read.
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
        print(jheader)
        ctrl_head = json.loads(jheader)

        ng = ctrl_head['ngroups']
        n_k = ctrl_head['nk']
# Note the difference in name conventions!
        SP = ctrl_head['ns'] - 1
        SO = ctrl_head['nc_flag']

        # load reciprocal basis
        kpt_basis = numpy.zeros((3,3))
        kpt_basis[:,0] = ctrl_head['kvec1']
        kpt_basis[:,1] = ctrl_head['kvec2']
        kpt_basis[:,2] = ctrl_head['kvec3']

        kpts = numpy.zeros((n_k, 3))
        kpts_cart = numpy.zeros((n_k, 3))
        bz_weights = numpy.zeros(n_k)
        kpt_weights = numpy.zeros(n_k)
        try:
            for ik in range(n_k):
                kx, ky, kz = next(rf), next(rf), next(rf)
                kpts[ik, :] = kx, ky, kz
                bz_weights[ik] = next(rf)
                # bz_weights soon to be removed, and replaced by kpt_weights
                kpt_weights[ik] = bz_weights[ik]
            for ik in range(n_k):
                kx, ky, kz = next(rf), next(rf), next(rf)
                kpts_cart[ik, :] = kx, ky, kz
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
            for ig in range(ng):
                gr_file = self.basename + '.pg%i'%(ig + 1)
                jheader, rf = self.read_header_and_data(gr_file)
                gr_head = json.loads(jheader)


                nb_max = gr_head['nb_max']
                p_shells = gr_head['shells']
                density_required = gr_head['nelect']
                charge_below = 0.0 # This is not defined in VASP interface

# Note that in the DftTools convention each site gives a separate correlated shell!
                n_shells = sum([len(sh['ion_list']) for sh in p_shells])
                n_corr_shells = sum([len(sh['ion_list']) for sh in p_shells])

                shells = []
                corr_shells = []
                shion_to_shell = [[] for ish in range(len(p_shells))]
                cr_shion_to_shell = [[] for ish in range(len(p_shells))]
                shorbs_to_globalorbs = [[] for ish in range(len(p_shells))]
                last_dimension = 0
                crshorbs_to_globalorbs = []
                icsh = 0
                for ish, sh in enumerate(p_shells):
                    ion_list = sh['ion_list']
                    for i, ion in enumerate(ion_list):
                        pars = {}
                        pars['atom'] = ion
# We set all sites inequivalent
                        pars['sort'] = sh['ion_sort'][i]
                        pars['l'] = sh['lorb']
                        #pars['corr'] = sh['corr']
                        pars['dim'] = sh['ndim']
                        #pars['ion_list'] = sh['ion_list']
                        pars['SO'] = SO
# TODO: check what 'irep' entry does (it seems to be very specific to dmftproj)
                        pars['irep'] = 0
                        shells.append(pars)
                        shion_to_shell[ish].append(i)
                        shorbs_to_globalorbs[ish].append([last_dimension,
                                                 last_dimension+sh['ndim']])
                        last_dimension = last_dimension+sh['ndim']
                        if sh['corr']:
                            corr_shells.append(pars)


# TODO: generalize this to the case of multiple shell groups
                n_corr_shells = len(corr_shells)

            n_orbs = sum([sh['dim'] for sh in shells])

# FIXME: atomic sorts in Wien2K are not the same as in VASP.
#        A symmetry analysis from OUTCAR or symmetry file should be used
#        to define equivalence classes of sites.
            n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = ConverterTools.det_shell_equivalence(self, corr_shells)

            if mpi.is_master_node():
                print("  No. of inequivalent shells:", n_inequiv_shells)

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
            band_window = [numpy.zeros((n_k, 2), dtype=int) for isp in range(n_spin_blocs)]
            n_orbitals = numpy.zeros([n_k, n_spin_blocs], numpy.int)


            for isp in range(n_spin_blocs):
                for ik in range(n_k):
                    ib1, ib2 = int(next(rf)), int(next(rf))
                    band_window[isp][ik, :2] = ib1, ib2
                    nb = ib2 - ib1 + 1
                    n_orbitals[ik, isp] = nb
                    for ib in range(nb):
                        hopping[ik, isp, ib, ib] = next(rf)
                        f_weights[ik, isp, ib] = next(rf)

            if self.proj_or_hk == 'hk':
                hopping = numpy.zeros([n_k, n_spin_blocs, n_orbs, n_orbs], numpy.complex_)
                # skip header lines
                hk_file = self.basename + '.hk%i'%(ig + 1)
                f_hk = open(hk_file, 'rt')
                # skip the header (1 line for n_kpoints, n_electrons, n_shells)
                # and one line per shell
                count = 0
                while count <  3 + n_shells:
                    f_hk.readline()
                    count += 1
                rf_hk = self.read_data(f_hk)
                for isp in range(n_spin_blocs):
                    for ik in range(n_k):
                        n_orbitals[ik, isp] = n_orbs
                        for ib in range(n_orbs):
                            for jb in range(n_orbs):
                                hopping[ik, isp, ib, jb] = next(rf_hk)
                        for ib in range(n_orbs):
                            for jb in range(n_orbs):
                                hopping[ik, isp, ib, jb] += 1j*next(rf_hk)
                rf_hk.close()

# Projectors
#            print n_orbitals
#            print [crsh['dim'] for crsh in corr_shells]
            proj_mat_csc = numpy.zeros([n_k, n_spin_blocs, sum([sh['dim'] for sh in shells]), numpy.max(n_orbitals)], numpy.complex_)

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
                for isp in range(n_spin_blocs):
                    for ik in range(n_k):
                        for ion in range(len(sh['ion_list'])):
                            for ilm in range(shorbs_to_globalorbs[ish][ion][0],shorbs_to_globalorbs[ish][ion][1]):
                                for ib in range(n_orbitals[ik, isp]):
                                    # This is to avoid confusion with the order of arguments
                                    pr = next(rf)
                                    pi = next(rf)
                                    proj_mat_csc[ik, isp, ilm, ib] = complex(pr, pi)

# now save only projectors with flag 'corr' to proj_mat
            proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max([crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], numpy.complex_)
            if self.proj_or_hk == 'proj':
                for ish, sh in enumerate(p_shells):
                    if sh['corr']:
                        for isp in range(n_spin_blocs):
                            for ik in range(n_k):
                                for ion in range(len(sh['ion_list'])):
                                    icsh = shion_to_shell[ish][ion]
                                    for iclm,ilm in enumerate(range(shorbs_to_globalorbs[ish][ion][0],shorbs_to_globalorbs[ish][ion][1])):
                                        for ib in range(n_orbitals[ik, isp]):
                                            proj_mat[ik,isp,icsh,iclm,ib] = proj_mat_csc[ik,isp,ilm,ib]
            elif self.proj_or_hk == 'hk':

                for ish, sh in enumerate(p_shells):
                    if sh['corr']:
                        for ion in range(len(sh['ion_list'])):
                            icsh = shion_to_shell[ish][ion]
                            for isp in range(n_spin_blocs):
                                for ik in range(n_k):
                                    for iclm,ilm in enumerate(range(shorbs_to_globalorbs[ish][ion][0],shorbs_to_globalorbs[ish][ion][1])):
                                        proj_mat[ik,isp,icsh,iclm,ilm] = 1.0

            #corr_shell.pop('ion_list')
            things_to_set = ['n_shells','shells','n_corr_shells','corr_shells','n_spin_blocs','n_orbitals','n_k','SO','SP','energy_unit']
            for it in things_to_set:
#                print "%s:"%(it), locals()[it]
                setattr(self,it,locals()[it])

        except StopIteration:
           raise "VaspConverter: error reading %s"%self.gr_file

        rf.close()


        proj_or_hk = self.proj_or_hk

        # Save it to the HDF:
        with HDFArchive(self.hdf_file,'a') as ar:
            if not (self.dft_subgrp in ar): ar.create_group(self.dft_subgrp)
            # The subgroup containing the data. If it does not exist, it is created. If it exists, the data is overwritten!
            things_to_save = ['energy_unit','n_k','k_dep_projection','SP','SO','charge_below','density_required',
                              'symm_op','n_shells','shells','n_corr_shells','corr_shells','use_rotations','rot_mat',
                              'rot_mat_time_inv','n_reps','dim_reps','T','n_orbitals','proj_mat','bz_weights',
                              'hopping','n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr','proj_or_hk',
                              'kpts','kpt_weights', 'kpt_basis']
            if self.proj_or_hk == 'hk' or self.proj_or_hk ==  True:
                things_to_save.append('proj_mat_csc')
            for it in things_to_save: ar[self.dft_subgrp][it] = locals()[it]

            # Store Fermi weights to 'dft_misc_input'
            if not (self.misc_subgrp in ar): ar.create_group(self.misc_subgrp)
            ar[self.misc_subgrp]['dft_fermi_weights'] = f_weights
            ar[self.misc_subgrp]['kpts_cart'] = kpts_cart
            ar[self.misc_subgrp]['band_window'] = band_window

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

        Parameters
        ----------
        bandwin_file : string
                    filename of .oubwin/up/dn file.
        struct_file : string
                    filename of .struct file.
        outputs_file : string
                    filename of .outputs file.
        misc_subgrp : string
                    name of the subgroup in which to save
        SO : boolean
            spin-orbit switch
        SP : int
            spin
        n_k : int
            number of k-points

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
                assert int(next(R)) == n_k, "convert_misc_input: Number of k-points is inconsistent in oubwin file!"
                assert int(next(R)) == SO, "convert_misc_input: SO is inconsistent in oubwin file!"
                for ik in range(n_k):
                    next(R)
                    band_window[isp][ik,0] = next(R) # lowest band
                    band_window[isp][ik,1] = next(R) # highest band
                    next(R)
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
        with HDFArchive(self.hdf_file,'a') as ar:
            if not (misc_subgrp in ar): ar.create_group(misc_subgrp)
            for it in things_to_save: ar[misc_subgrp][it] = locals()[it]


    def convert_symmetry_input(self, ctrl_head, orbits, symm_subgrp):
        """
        Reads input for the symmetrisations from symm_file, which is case.sympar or case.symqmc.

        Parameters
        ----------
        ctrl_head : dict
                dictionary of header of .ctrl file
        orbits : list of shells
                contains all shells
        symm_subgrp : name of symmetry group in h5 archive

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
        with HDFArchive(self.hdf_file,'a') as ar:
            if not (symm_subgrp in ar): ar.create_group(symm_subgrp)
            things_to_save = ['n_symm','n_atoms','perm','orbits','SO','SP','time_inv','mat','mat_tinv']
            for it in things_to_save:
                ar[symm_subgrp][it] = locals()[it]
