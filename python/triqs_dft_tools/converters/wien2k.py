
##########################################################################
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
##########################################################################

from types import *
import numpy
from h5 import *
from .converter_tools import *
import os.path


class Wien2kConverter(ConverterTools):
    """
    Conversion from Wien2k output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, filename, hdf_filename=None,
                 dft_subgrp='dft_input', symmcorr_subgrp='dft_symmcorr_input',
                 parproj_subgrp='dft_parproj_input', symmpar_subgrp='dft_symmpar_input',
                 bands_subgrp='dft_bands_input', misc_subgrp='dft_misc_input',
                 transp_subgrp='dft_transp_input', repacking=False):
        """
        Initialise the class.

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

        """

        assert isinstance(filename, str), "Wien2kConverter: Please provide the DFT files' base name as a string."
        if hdf_filename is None:
            hdf_filename = filename + '.h5'
        self.hdf_file = hdf_filename
        self.dft_file = filename + '.ctqmcout'
        self.symmcorr_file = filename + '.symqmc'
        self.parproj_file = filename + '.parproj'
        self.symmpar_file = filename + '.sympar'
        self.band_file = filename + '.outband'
        self.bandwin_file = filename + '.oubwin'
        self.struct_file = filename + '.struct'
        self.outputs_file = filename + '.outputs'
        self.pmat_file = filename + '.pmat'
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.parproj_subgrp = parproj_subgrp
        self.symmpar_subgrp = symmpar_subgrp
        self.bands_subgrp = bands_subgrp
        self.misc_subgrp = misc_subgrp
        self.transp_subgrp = transp_subgrp
        self.fortran_to_replace = {'D': 'E'}

        # Checks if h5 file is there and repacks it if wanted:
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)

    def convert_dft_input(self):
        """
        Reads the appropriate files and stores the data for the

        - dft_subgrp
        - symmcorr_subgrp
        - misc_subgrp

        in the hdf5 archive. 

        """

        # Read and write only on the master node
        if not (mpi.is_master_node()):
            return
        mpi.report("Reading input from %s..." % self.dft_file)

        # R is a generator : each R.Next() will return the next number in the
        # file
        R = ConverterTools.read_fortran_file(
            self, self.dft_file, self.fortran_to_replace)
        try:
            energy_unit = next(R)                        # read the energy convertion factor
            # read the number of k points
            n_k = int(next(R))
            k_dep_projection = 1
            # flag for spin-polarised calculation
            SP = int(next(R))
            # flag for spin-orbit calculation
            SO = int(next(R))
            charge_below = next(R)                       # total charge below energy window
            # total density required, for setting the chemical potential
            density_required = next(R)
            symm_op = 1                                   # Use symmetry groups for the k-sum

            # the information on the non-correlated shells is not important
            # here, maybe skip:
            # number of shells (e.g. Fe d, As p, O p) in the unit cell,
            n_shells = int(next(R))
            # corresponds to index R in formulas
            # now read the information about the shells (atom, sort, l, dim):
            shell_entries = ['atom', 'sort', 'l', 'dim']
            shells = [{name: int(val) for name, val in zip(
                shell_entries, R)} for ish in range(n_shells)]

            # number of corr. shells (e.g. Fe d, Ce f) in the unit cell,
            n_corr_shells = int(next(R))
            # corresponds to index R in formulas
            # now read the information about the shells (atom, sort, l, dim, SO
            # flag, irep):
            corr_shell_entries = ['atom', 'sort', 'l', 'dim', 'SO', 'irep']
            corr_shells = [{name: int(val) for name, val in zip(
                corr_shell_entries, R)} for icrsh in range(n_corr_shells)]

            # determine the number of inequivalent correlated shells and maps,
            # needed for further reading
            n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = ConverterTools.det_shell_equivalence(
                self, corr_shells)

            use_rotations = 1
            rot_mat = [numpy.identity(
                corr_shells[icrsh]['dim'], numpy.complex_) for icrsh in range(n_corr_shells)]

            # read the matrices
            rot_mat_time_inv = [0 for i in range(n_corr_shells)]

            for icrsh in range(n_corr_shells):
                for i in range(corr_shells[icrsh]['dim']):    # read real part:
                    for j in range(corr_shells[icrsh]['dim']):
                        rot_mat[icrsh][i, j] = next(R)
                # read imaginary part:
                for i in range(corr_shells[icrsh]['dim']):
                    for j in range(corr_shells[icrsh]['dim']):
                        rot_mat[icrsh][i, j] += 1j * next(R)

                if (SP == 1):             # read time inversion flag:
                    rot_mat_time_inv[icrsh] = int(next(R))

            # Read here the info for the transformation of the basis:
            n_reps = [1 for i in range(n_inequiv_shells)]
            dim_reps = [0 for i in range(n_inequiv_shells)]
            T = []
            for ish in range(n_inequiv_shells):
                # number of representatives ("subsets"), e.g. t2g and eg
                n_reps[ish] = int(next(R))
                dim_reps[ish] = [int(next(R)) for i in range(
                    n_reps[ish])]   # dimensions of the subsets

                # The transformation matrix:
                # is of dimension 2l+1 without SO, and 2*(2l+1) with SO!
                ll = 2 * corr_shells[inequiv_to_corr[ish]]['l'] + 1
                lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
                T.append(numpy.zeros([lmax, lmax], numpy.complex_))

                # now read it from file:
                for i in range(lmax):
                    for j in range(lmax):
                        T[ish][i, j] = next(R)
                for i in range(lmax):
                    for j in range(lmax):
                        T[ish][i, j] += 1j * next(R)

            # Spin blocks to be read:
            n_spin_blocs = SP + 1 - SO

            # read the list of n_orbitals for all k points
            n_orbitals = numpy.zeros([n_k, n_spin_blocs], numpy.int)
            for isp in range(n_spin_blocs):
                for ik in range(n_k):
                    n_orbitals[ik, isp] = int(next(R))

            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k, n_spin_blocs, n_corr_shells, max(
                [crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], numpy.complex_)

            # Read the projectors from the file:
            for ik in range(n_k):
                for icrsh in range(n_corr_shells):
                    n_orb = corr_shells[icrsh]['dim']
                    # first Real part for BOTH spins, due to conventions in
                    # dmftproj:
                    for isp in range(n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik][isp]):
                                proj_mat[ik, isp, icrsh, i, j] = next(R)
                    # now Imag part:
                    for isp in range(n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik][isp]):
                                proj_mat[ik, isp, icrsh, i, j] += 1j * next(R)

            # now define the arrays for weights and hopping ...
            # w(k_index),  default normalisation
            bz_weights = numpy.ones([n_k], numpy.float_) / float(n_k)
            hopping = numpy.zeros([n_k, n_spin_blocs, numpy.max(
                n_orbitals), numpy.max(n_orbitals)], numpy.complex_)

            # weights in the file
            for ik in range(n_k):
                bz_weights[ik] = next(R)

            # if the sum over spins is in the weights, take it out again!!
            sm = sum(bz_weights)
            bz_weights[:] /= sm

            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian -- convention
            # for Wien2K.
            for isp in range(n_spin_blocs):
                for ik in range(n_k):
                    n_orb = n_orbitals[ik, isp]
                    for i in range(n_orb):
                        hopping[ik, isp, i, i] = next(R) * energy_unit

            # keep some things that we need for reading parproj:
            things_to_set = ['n_shells', 'shells', 'n_corr_shells', 'corr_shells',
                             'n_spin_blocs', 'n_orbitals', 'n_k', 'SO', 'SP', 'energy_unit']
            for it in things_to_set:
                setattr(self, it, locals()[it])
        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("wien2k : reading file %s failed!" % self.dft_file)

        R.close()
        # Reading done!

        # Save it to the HDF:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.dft_subgrp in ar):
                ar.create_group(self.dft_subgrp)
            # The subgroup containing the data. If it does not exist, it is
            # created. If it exists, the data is overwritten!
            things_to_save = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                          'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                          'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                          'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr']
            for it in things_to_save:
                ar[self.dft_subgrp][it] = locals()[it]

        # Symmetries are used, so now convert symmetry information for
        # *correlated* orbitals:
        self.convert_symmetry_input(orbits=self.corr_shells, symm_file=self.symmcorr_file,
                                    symm_subgrp=self.symmcorr_subgrp, SO=self.SO, SP=self.SP)
        self.convert_misc_input()

    def convert_parproj_input(self):
        """
        Reads the appropriate files and stores the data for the 

        - parproj_subgrp
        - symmpar_subgrp

        in the hdf5 archive. 

        """

        if not (mpi.is_master_node()):
            return

        # get needed data from hdf file
        with HDFArchive(self.hdf_file, 'a') as ar:
            things_to_read = ['SP', 'SO', 'n_shells',
                          'n_k', 'n_orbitals', 'shells']

            for it in things_to_read:
                if not hasattr(self, it):
                    setattr(self, it, ar[self.dft_subgrp][it])
            self.n_spin_blocs = self.SP + 1 - self.SO

        mpi.report("Reading input from %s..." % self.parproj_file)

        dens_mat_below = [[numpy.zeros([self.shells[ish]['dim'], self.shells[ish]['dim']], numpy.complex_) for ish in range(self.n_shells)]
                          for isp in range(self.n_spin_blocs)]

        R = ConverterTools.read_fortran_file(
            self, self.parproj_file, self.fortran_to_replace)

        n_parproj = [int(next(R)) for i in range(self.n_shells)]
        n_parproj = numpy.array(n_parproj)

        # Initialise P, here a double list of matrices:
        proj_mat_all = numpy.zeros([self.n_k, self.n_spin_blocs, self.n_shells, max(
            n_parproj), max([sh['dim'] for sh in self.shells]), numpy.max(self.n_orbitals)], numpy.complex_)

        rot_mat_all = [numpy.identity(
            self.shells[ish]['dim'], numpy.complex_) for ish in range(self.n_shells)]
        rot_mat_all_time_inv = [0 for i in range(self.n_shells)]

        for ish in range(self.n_shells):
            # read first the projectors for this orbital:
            for ik in range(self.n_k):
                for ir in range(n_parproj[ish]):

                    for isp in range(self.n_spin_blocs):
                        # read real part:
                        for i in range(self.shells[ish]['dim']):
                            for j in range(self.n_orbitals[ik][isp]):
                                proj_mat_all[ik, isp, ish, ir, i, j] = next(R)

                    for isp in range(self.n_spin_blocs):
                        # read imaginary part:
                        for i in range(self.shells[ish]['dim']):
                            for j in range(self.n_orbitals[ik][isp]):
                                proj_mat_all[ik, isp, ish,
                                             ir, i, j] += 1j * next(R)

            # now read the Density Matrix for this orbital below the energy
            # window:
            for isp in range(self.n_spin_blocs):
                for i in range(self.shells[ish]['dim']):    # read real part:
                    for j in range(self.shells[ish]['dim']):
                        dens_mat_below[isp][ish][i, j] = next(R)
            for isp in range(self.n_spin_blocs):
                # read imaginary part:
                for i in range(self.shells[ish]['dim']):
                    for j in range(self.shells[ish]['dim']):
                        dens_mat_below[isp][ish][i, j] += 1j * next(R)
                if (self.SP == 0):
                    dens_mat_below[isp][ish] /= 2.0

            # Global -> local rotation matrix for this shell:
            for i in range(self.shells[ish]['dim']):    # read real part:
                for j in range(self.shells[ish]['dim']):
                    rot_mat_all[ish][i, j] = next(R)
            for i in range(self.shells[ish]['dim']):    # read imaginary part:
                for j in range(self.shells[ish]['dim']):
                    rot_mat_all[ish][i, j] += 1j * next(R)

            if (self.SP):
                rot_mat_all_time_inv[ish] = int(next(R))

        R.close()
        # Reading done!

        # Save it to the HDF:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.parproj_subgrp in ar):
                ar.create_group(self.parproj_subgrp)
            # The subgroup containing the data. If it does not exist, it is
            # created. If it exists, the data is overwritten!
            things_to_save = ['dens_mat_below', 'n_parproj',
                          'proj_mat_all', 'rot_mat_all', 'rot_mat_all_time_inv']
            for it in things_to_save:
                ar[self.parproj_subgrp][it] = locals()[it]

        # Symmetries are used, so now convert symmetry information for *all*
        # orbitals:
        self.convert_symmetry_input(orbits=self.shells, symm_file=self.symmpar_file,
                                    symm_subgrp=self.symmpar_subgrp, SO=self.SO, SP=self.SP)

    def convert_bands_input(self):
        """
        Reads the appropriate files and stores the data for the bands_subgrp in the hdf5 archive. 

        """

        if not (mpi.is_master_node()):
            return

        try:
            # get needed data from hdf file
            with HDFArchive(self.hdf_file, 'a') as ar:
                things_to_read = ['SP', 'SO', 'n_corr_shells',
                              'n_shells', 'corr_shells', 'shells', 'energy_unit']

                for it in things_to_read:
                    if not hasattr(self, it):
                        setattr(self, it, ar[self.dft_subgrp][it])
                self.n_spin_blocs = self.SP + 1 - self.SO

            mpi.report("Reading input from %s..." % self.band_file)
            R = ConverterTools.read_fortran_file(
                self, self.band_file, self.fortran_to_replace)
            n_k = int(next(R))

            # read the list of n_orbitals for all k points
            n_orbitals = numpy.zeros([n_k, self.n_spin_blocs], numpy.int)
            for isp in range(self.n_spin_blocs):
                for ik in range(n_k):
                    n_orbitals[ik, isp] = int(next(R))

            # Initialise the projectors:
            proj_mat = numpy.zeros([n_k, self.n_spin_blocs, self.n_corr_shells, max(
                [crsh['dim'] for crsh in self.corr_shells]), numpy.max(n_orbitals)], numpy.complex_)

            # Read the projectors from the file:
            for ik in range(n_k):
                for icrsh in range(self.n_corr_shells):
                    n_orb = self.corr_shells[icrsh]['dim']
                    # first Real part for BOTH spins, due to conventions in
                    # dmftproj:
                    for isp in range(self.n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik, isp]):
                                proj_mat[ik, isp, icrsh, i, j] = next(R)
                    # now Imag part:
                    for isp in range(self.n_spin_blocs):
                        for i in range(n_orb):
                            for j in range(n_orbitals[ik, isp]):
                                proj_mat[ik, isp, icrsh, i, j] += 1j * next(R)

            hopping = numpy.zeros([n_k, self.n_spin_blocs, numpy.max(
                n_orbitals), numpy.max(n_orbitals)], numpy.complex_)

            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian!!!!
            for isp in range(self.n_spin_blocs):
                for ik in range(n_k):
                    n_orb = n_orbitals[ik, isp]
                    for i in range(n_orb):
                        hopping[ik, isp, i, i] = next(R) * self.energy_unit

            # now read the partial projectors:
            n_parproj = [int(next(R)) for i in range(self.n_shells)]
            n_parproj = numpy.array(n_parproj)

            # Initialise P, here a double list of matrices:
            proj_mat_all = numpy.zeros([n_k, self.n_spin_blocs, self.n_shells, max(n_parproj), max(
                [sh['dim'] for sh in self.shells]), numpy.max(n_orbitals)], numpy.complex_)

            for ish in range(self.n_shells):
                for ik in range(n_k):
                    for ir in range(n_parproj[ish]):
                        for isp in range(self.n_spin_blocs):

                            # read real part:
                            for i in range(self.shells[ish]['dim']):
                                for j in range(n_orbitals[ik, isp]):
                                    proj_mat_all[ik, isp, ish,
                                                 ir, i, j] = next(R)

                            # read imaginary part:
                            for i in range(self.shells[ish]['dim']):
                                for j in range(n_orbitals[ik, isp]):
                                    proj_mat_all[ik, isp, ish,
                                                 ir, i, j] += 1j * next(R)

            R.close()

        except KeyError:
            raise IOError("convert_bands_input : Needed data not found in hdf file. Consider calling convert_dft_input first!")
        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("wien2k : reading file %s failed!" % self.band_file)

        # Reading done!

        # Save it to the HDF:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.bands_subgrp in ar):
                ar.create_group(self.bands_subgrp)
            # The subgroup containing the data. If it does not exist, it is
            # created. If it exists, the data is overwritten!
            things_to_save = ['n_k', 'n_orbitals', 'proj_mat',
                          'hopping', 'n_parproj', 'proj_mat_all']
            for it in things_to_save:
                ar[self.bands_subgrp][it] = locals()[it]

    def convert_misc_input(self):
        """
        Reads additional information on:

        - the band window from :file:`case.oubwin`,
        - lattice parameters from :file:`case.struct`,
        - symmetries from :file:`case.outputs`,

        if those Wien2k files are present and stores the data in the hdf5 archive.
        This function is automatically called by :meth:`convert_dft_input <triqs_dft_tools.converters.wien2k.Wien2kConverter.convert_dft_input>`. 

        """

        if not (mpi.is_master_node()):
            return

        # Check if SP, SO and n_k are already in h5
        with HDFArchive(self.hdf_file, 'r') as ar:
            if not (self.dft_subgrp in ar):
                raise IOError("convert_misc_input: No %s subgroup in hdf file found! Call convert_dft_input first." % self.dft_subgrp)
            SP = ar[self.dft_subgrp]['SP']
            SO = ar[self.dft_subgrp]['SO']
            n_k = ar[self.dft_subgrp]['n_k']

        things_to_save = []

        # Read relevant data from .oubwin/up/dn files
        #############################################
        # band_window: Contains the index of the lowest and highest band within the
        #              projected subspace (used by dmftproj) for each k-point.

        if (SP == 0 and SO == 0): # read .oubwin file
            files = [self.bandwin_file]
        elif (SP == 1 and SO == 0): # read .oubwinup and .oubwindn
            files = [self.bandwin_file + 'up', self.bandwin_file + 'dn']
        elif (SP == 1 and SO == 1): # read either .oubwinup or .oubwindn
            if os.path.exists(self.bandwin_file + 'up'):
                files = [self.bandwin_file + 'up']
            elif os.path.exists(self.bandwin_file + 'dn'):
                files = [self.bandwin_file + 'dn']
            else:
                assert 0, "convert_misc_input: If SO and SP are 1 provide either .oubwinup or .oubwindn file"
        else:
            assert 0, "convert_misc_input: Reading oubwin error! Check SP and SO, if SO=1 SP must be 1."

        band_window = [None for isp in range(SP + 1 - SO)]
        for isp, f in enumerate(files):
            if os.path.exists(f):
                mpi.report("Reading input from %s..." % f)
                R = ConverterTools.read_fortran_file(
                    self, f, self.fortran_to_replace)
                n_k_oubwin = int(next(R))
                if (n_k_oubwin != n_k):
                    mpi.report(
                        "convert_misc_input : WARNING : n_k in case.oubwin is different from n_k in case.klist")
                assert int(
                    next(R)) == SO, "convert_misc_input: SO is inconsistent in oubwin file!"

                band_window[isp] = numpy.zeros((n_k_oubwin, 2), dtype=int)
                for ik in range(n_k_oubwin):
                    next(R)
                    band_window[isp][ik, 0] = next(R)  # lowest band
                    band_window[isp][ik, 1] = next(R)  # highest band
                    next(R)
                things_to_save.append('band_window')

                R.close()  # Reading done!

        # Read relevant data from .struct file
        ######################################
        # lattice_type: bravais lattice type as defined by Wien2k
        # lattice_constants: unit cell parameters in a. u.
        # lattice_angles: unit cell angles in rad

        if (os.path.exists(self.struct_file)):
            mpi.report("Reading input from %s..." % self.struct_file)

            with open(self.struct_file) as R:
                try:
                    R.readline()
                    lattice_type = R.readline().split()[0]
                    R.readline()
                    temp = R.readline()
                    lattice_constants = numpy.array(
                        [float(temp[0 + 10 * i:10 + 10 * i].strip()) for i in range(3)])
                    lattice_angles = numpy.array(
                        [float(temp[30 + 10 * i:40 + 10 * i].strip()) for i in range(3)]) * numpy.pi / 180.0
                    things_to_save.extend(
                        ['lattice_type', 'lattice_constants', 'lattice_angles'])
                except IOError:
                    raise IOError("convert_misc_input: reading file %s failed" % self.struct_file)

        # Read relevant data from .outputs file
        #######################################
        # rot_symmetries: matrix representation of all (space group) symmetry
        # operations

        if (os.path.exists(self.outputs_file)):
            mpi.report("Reading input from %s..." % self.outputs_file)

            rot_symmetries = []
            with open(self.outputs_file) as R:
                try:
                    while 1:
                        temp = R.readline().strip(' ').split()
                        if (temp[0] == 'PGBSYM:'):
                            n_symmetries = int(temp[-1])
                            break
                    for i in range(n_symmetries):
                        while 1:
                            if (R.readline().strip().split()[0] == 'Symmetry'):
                                break
                        sym_i = numpy.zeros((3, 3), dtype=float)
                        for ir in range(3):
                            temp = R.readline().strip().split()
                            for ic in range(3):
                                sym_i[ir, ic] = float(temp[ic])
                        R.readline()
                        rot_symmetries.append(sym_i)
                    things_to_save.extend(['n_symmetries', 'rot_symmetries'])
                    things_to_save.append('rot_symmetries')
                except IOError:
                    raise IOError("convert_misc_input: reading file %s failed" % self.outputs_file)

        # Save it to the HDF:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.misc_subgrp in ar):
                ar.create_group(self.misc_subgrp)
            for it in things_to_save:
                ar[self.misc_subgrp][it] = locals()[it]

    def convert_transport_input(self):
        """ 
        Reads the necessary information for transport calculations on:

        - the optical band window and the velocity matrix elements from :file:`case.pmat`

        and stores the data in the hdf5 archive.

        """

        if not (mpi.is_master_node()):
            return

        # Check if SP, SO and n_k are already in h5
        with HDFArchive(self.hdf_file, 'r') as ar:
            if not (self.dft_subgrp in ar):
                raise IOError("convert_transport_input: No %s subgroup in hdf file found! Call convert_dft_input first." % self.dft_subgrp)
            SP = ar[self.dft_subgrp]['SP']
            SO = ar[self.dft_subgrp]['SO']
            n_k = ar[self.dft_subgrp]['n_k']

        # Read relevant data from .pmat/up/dn files
        ###########################################
        # band_window_optics: Contains the index of the lowest and highest band within the
        #                     band window (used by optics) for each k-point.
        # velocities_k: velocity (momentum) matrix elements between all bands in band_window_optics
        #               and each k-point.

        if (SP == 0 and SO == 0): # read .pmat file
            files = [self.pmat_file]
        elif (SP == 1 and SO == 0): # read .pmatup and pmatdn
            files = [self.pmat_file + 'up', self.pmat_file + 'dn']
        elif (SP == 1 and SO == 1): # read either .pmatup or .pmatdn
            if os.path.exists(self.pmat_file + 'up'):
                files = [self.pmat_file + 'up']
            elif os.path.exists(self.pmat_file + 'dn'):
                files = [self.pmat_file + 'dn']
            else:
                assert 0, "convert_transport_input: If SO and SP are 1 provide either .pmatup or .pmatdn file"
        else:
            assert 0, "convert_transport_input: Reading velocity file error! Check SP and SO, if SO=1 SP must be 1."

        velocities_k = [[] for f in files]
        band_window_optics = []
        for isp, f in enumerate(files):
            if not os.path.exists(f):
                raise IOError("convert_transport_input: File %s does not exist" % f)
            mpi.report("Reading input from %s..." % f)

            R = ConverterTools.read_fortran_file(
                self, f, {'D': 'E', '(': '', ')': '', ',': ' '})
            band_window_optics_isp = []
            for ik in range(n_k):
                next(R)
                nu1 = int(next(R))
                nu2 = int(next(R))
                band_window_optics_isp.append((nu1, nu2))
                n_bands = nu2 - nu1 + 1
                for _ in range(4):
                    next(R)
                if n_bands <= 0:
                    velocity_xyz = numpy.zeros((1, 1, 3), dtype=complex)
                else:
                    velocity_xyz = numpy.zeros(
                        (n_bands, n_bands, 3), dtype=complex)
                    for nu_i in range(n_bands):
                        for nu_j in range(nu_i, n_bands):
                            for i in range(3):
                                velocity_xyz[nu_i][nu_j][
                                    i] = next(R) + next(R) * 1j
                                if (nu_i != nu_j):
                                    velocity_xyz[nu_j][nu_i][i] = velocity_xyz[
                                        nu_i][nu_j][i].conjugate()
                velocities_k[isp].append(velocity_xyz)
            band_window_optics.append(numpy.array(band_window_optics_isp))
            R.close()  # Reading done!

        # Put data to HDF5 file
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (self.transp_subgrp in ar):
                ar.create_group(self.transp_subgrp)
            # The subgroup containing the data. If it does not exist, it is
            # created. If it exists, the data is overwritten!!!
            things_to_save = ['band_window_optics', 'velocities_k']
            for it in things_to_save:
                ar[self.transp_subgrp][it] = locals()[it]

    def convert_symmetry_input(self, orbits, symm_file, symm_subgrp, SO, SP):
        """
        Reads and stores symmetrisation data from symm_file, which can be is case.sympar or case.symqmc.

        Parameters
        ----------
        orbits : list of dicts
                 This is either shells or corr_shells depending on whether the symmetry 
                 information is for correlated shells or partial projectors.
        symm_file : string
                    Name of the file containing symmetry data. 
                    This is case.symqmc for correlated shells and case.sympar for partial projectors.
        symm_subgrp : string, optional
                      Name of subgroup storing symmetry data.
        SO : integer
             Is spin-orbit coupling considered?
        SP : integer
             Is the system spin-polarised?

        """

        if not (mpi.is_master_node()):
            return
        mpi.report("Reading input from %s..." % symm_file)

        n_orbits = len(orbits)

        R = ConverterTools.read_fortran_file(
            self, symm_file, self.fortran_to_replace)

        try:
            n_symm = int(next(R))           # Number of symmetry operations
            n_atoms = int(next(R))       # number of atoms involved
            perm = [[int(next(R)) for i in range(n_atoms)]
                    for j in range(n_symm)]    # list of permutations of the atoms
            if SP:
                # time inversion for SO coupling
                time_inv = [int(next(R)) for j in range(n_symm)]
            else:
                time_inv = [0 for j in range(n_symm)]

            # Now read matrices:
            mat = []
            for i_symm in range(n_symm):

                mat.append([numpy.zeros([orbits[orb]['dim'], orbits[orb][
                           'dim']], numpy.complex_) for orb in range(n_orbits)])
                for orb in range(n_orbits):
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            # real part
                            mat[i_symm][orb][i, j] = next(R)
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            mat[i_symm][orb][i, j] += 1j * \
                                next(R)      # imaginary part

            mat_tinv = [numpy.identity(orbits[orb]['dim'], numpy.complex_)
                        for orb in range(n_orbits)]

            if ((SO == 0) and (SP == 0)):
                # here we need an additional time inversion operation, so read
                # it:
                for orb in range(n_orbits):
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            # real part
                            mat_tinv[orb][i, j] = next(R)
                    for i in range(orbits[orb]['dim']):
                        for j in range(orbits[orb]['dim']):
                            mat_tinv[orb][i, j] += 1j * \
                                next(R)      # imaginary part

        except StopIteration:  # a more explicit error if the file is corrupted.
            raise IOError("wien2k : reading file %s failed!" %symm_file)

        R.close()
        # Reading done!

        # Save it to the HDF:
        with HDFArchive(self.hdf_file, 'a') as ar:
            if not (symm_subgrp in ar):
                ar.create_group(symm_subgrp)
            things_to_save = ['n_symm', 'n_atoms', 'perm',
                          'orbits', 'SO', 'SP', 'time_inv', 'mat', 'mat_tinv']
            for it in things_to_save:
                ar[symm_subgrp][it] = locals()[it]
