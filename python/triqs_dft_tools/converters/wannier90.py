
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

###
#  Wannier90 to HDF5 converter for the SumkDFT class of dfttools/TRIQS;
#
#   written by Gabriele Sclauzero (Materials Theory, ETH Zurich), Dec 2015 -- Jan 2016,
#   updated by Maximilian Merkel (Materials Theory, ETH Zurich), Aug 2020 -- Jan 2021,
#   and by Sophie Beck (Materials Theory, ETH Zurich), Sep 2020 -- Apr 2021,
#   under the supervision of Claude Ederer (Materials Theory).
#   Partially based on previous work by K. Dymkovski and the DFT_tools/TRIQS team.
#
#  Limitations of the current implementation:
# - the T rotation matrices are not used in this implementation
# - in bloch_basis mode, the number of wannier functions per shell has to be equal
#   to the dim of the shell
#
#  Things to be improved/checked:
# - the case with SP=1 might work, but was never tested (do we need to define
#   rot_mat_time_inv also if symm_op = 0?)
# - the calculation of rot_mat in find_rot_mat() relies on the eigenvalues of H(0);
#   this might fail in presence of degenerate eigenvalues (now just prints warning)
# - make the code more MPI safe (error handling): if we run with more than one process
#   and an error occurs on the masternode, the calculation does not abort
# - in case of disentanglement, the outer window being close to Kohn-Sham energies
#   can cause a problem in creating the udis_mat in read_wannier90data
# - add_lambda does not work for multiple impurities
###
"""
Wannier90 converter
"""

import numpy
import os.path
from itertools import product

from h5 import HDFArchive
from .converter_tools import ConverterTools
import triqs.utility.mpi as mpi

class Wannier90Converter(ConverterTools):
    """
    Conversion from Wannier90 output to an hdf5 file that can be used as input for the SumkDFT class.
    """

    def __init__(self, seedname, hdf_filename=None, dft_subgrp='dft_input',
                 symmcorr_subgrp='dft_symmcorr_input', misc_subgrp='dft_misc_input',
                 repacking=False, rot_mat_type='hloc_diag', bloch_basis=False, add_lambda=None,
                 w90zero=2.e-6):
        """
        Initialise the class.

        Parameters
        ----------
        seedname : string
            Base name of Wannier90 files
        hdf_filename : string, optional
            Name of hdf5 archive to be created
        dft_subgrp : string, optional
            Name of subgroup storing necessary DFT data
        symmcorr_subgrp : string, optional
            Name of subgroup storing correlated-shell symmetry data
        misc_subgrp : string, optional
            Name of subgroup storing miscellaneous DFT data.
        repacking : boolean, optional
            Does the hdf5 archive need to be repacked to save space?
        rot_mat_type : string, optional
            Type of rot_mat used
            Can be 'hloc_diag', 'wannier', 'none'
        bloch_basis : boolean, optional
            Should the Hamiltonian be written in Bloch rather than Wannier basis?
        add_lambda : list of floats, optional
            add local spin-orbit term
        w90zero : float, optional
            threshold on symmetry checks of Hamiltonian and rot_mat
        """

        self._name = "Wannier90Converter"
        assert isinstance(seedname, str), self._name + \
            ": Please provide the DFT files' base name as a string."
        if hdf_filename is None:
            hdf_filename = seedname + '.h5'
        self.hdf_file = hdf_filename
        # if the w90 output is seedname_hr.dat, the input file for the
        # converter must be called seedname.inp
        self.inp_file = seedname + '.inp'
        self.w90_seed = seedname
        self.dft_subgrp = dft_subgrp
        self.symmcorr_subgrp = symmcorr_subgrp
        self.misc_subgrp = misc_subgrp
        self.fortran_to_replace = {'D': 'E'}
        # threshold below which matrix elements from wannier90 should be
        # considered equal
        self._w90zero = w90zero
        self.rot_mat_type = rot_mat_type
        self.bloch_basis = bloch_basis
        self.add_lambda = add_lambda
        if self.rot_mat_type not in ('hloc_diag', 'wannier', 'none'):
            raise ValueError('Parameter rot_mat_type invalid, should be one of'
                             + '"hloc_diag", "wannier", "none"')

        # Checks if h5 file is there and repacks it if wanted:
        if (os.path.exists(self.hdf_file) and repacking):
            ConverterTools.repack(self)

    def convert_dft_input(self):
        """
        Reads the appropriate files and stores the data for the

        - dft_subgrp
        - symmcorr_subgrp

        in the hdf5 archive.

        """

        mpi.report("\nReading input from %s..." % self.inp_file)

        # R is a generator : each R.Next() will return the next number in the
        # file
        R = ConverterTools.read_fortran_file(
            self, self.inp_file, self.fortran_to_replace)
        shell_entries = ['atom', 'sort', 'l', 'dim']
        corr_shell_entries = ['atom', 'sort', 'l', 'dim', 'SO', 'irep']
        # First, let's read the input file with the parameters needed for the
        # conversion
        try:
            # read k - point mesh generation option
            kmesh_mode = int(next(R))
            if kmesh_mode >= 0:
                # read k-point mesh size from input
                nki = [int(next(R)) for idir in range(3)]
            else:
                # some default grid, if everything else fails...
                nki = [8, 8, 8]
            # read the total number of electrons per cell if not in bloch basis
            # in bloch basis, this is later calculated from the partial occupations
            density_required = float(next(R))

            # we do not read shells, because we have no additional shells beyond correlated ones,
            # and the data will be copied from corr_shells into shells (see below)
            # number of corr. shells (e.g. Fe d, Ce f) in the unit cell,
            n_corr_shells = int(next(R))
            # now read the information about the correlated shells (atom, sort,
            # l, dim, SO flag, irep):
            corr_shells = [{name: int(val) for name, val in zip(
                corr_shell_entries, R)} for icrsh in range(n_corr_shells)]

            try:
                self.fermi_energy = float(next(R))
            except:
                self.fermi_energy = 0.
        except StopIteration:  # a more explicit error if the file is corrupted.
            mpi.report(self._name + ": reading input file %s failed!" %
                       self.inp_file)
        # close the input file
        R.close()

        # Set or derive some quantities
        # Wannier90 does not use symmetries to reduce the k-points
        # the following might change in future versions
        symm_op = 0
        # copy corr_shells into shells (see above)
        n_shells = n_corr_shells
        shells = []
        for ish in range(n_shells):
            shells.append({key: corr_shells[ish].get(
                key, None) for key in shell_entries})

        # Determine if any shell requires SO
        if any([corr_shell['SO']==1 for corr_shell in corr_shells]):
            SO = 1
            SP = 1
            mpi.report('Spin-orbit interaction turned on')
        else:
            SO = 0
            SP = 0
        # Only one block supported - either non-spin-polarized or spin-orbit coupled
        assert SP == SO, 'Spin-polarized calculations not implemented'
        if self.add_lambda:
            assert n_shells == 1, 'add_lambda not implemented for more than one t2g shell'
            assert [sh['dim'] for sh in corr_shells] == [3 for sh in corr_shells], 'add_lambda only implemented for t2g shell'
            assert SO == SP == 0, 'add_lambda not implemented for SO = SP = 1'
            assert self.bloch_basis == False, 'add_lambda not implemented for bloch_basis = True'
            # now setting SO and SP to 1
            SO = SP = 1

        charge_below = 0                # total charge below energy window NOT used for now
        energy_unit = 1.0               # should be understood as eV units

        # this is more general
        n_spin_blocs = SP + 1 - SO
        assert n_spin_blocs > 0, 'Input error, if SO=1, SP must be 1.'

        if SO == 1:
            for shell_list in [shells, corr_shells]:
                for entry in shell_list:
                    entry['dim'] *= 2
                    if 'SO' in entry.keys() and self.add_lambda: entry['SO'] = 1

        dim_corr_shells = sum([sh['dim'] for sh in corr_shells])
        mpi.report('Total number of WFs expected in the correlated shells: {0:d}'.format(dim_corr_shells))

        # build the k-point mesh, if its size was given on input (kmesh_mode >= 0),
        # otherwise it is built according to the data in the hr file (see below)
        # If output is in bloch_basis, we use k mesh from seedname_u.mat for consistency
        if kmesh_mode >= 0 and not self.bloch_basis:
            n_k, kpts, kpt_weights = self.kmesh_build(nki, kmesh_mode)
            self.n_k = n_k
            self.kpts = kpts

        # determine the number of inequivalent correlated shells and maps,
        # needed for further processing
        n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = ConverterTools.det_shell_equivalence(
            self, corr_shells)
        mpi.report("Number of inequivalent shells: %d" % n_inequiv_shells)
        mpi.report("Shell representatives: " + format(inequiv_to_corr))
        shells_map = [inequiv_to_corr[corr_to_inequiv[ish]]
                      for ish in range(n_corr_shells)]
        mpi.report("Mapping: " + format(shells_map))

        # not used in this version: reset to dummy values?
        n_reps = [1 for i in range(n_inequiv_shells)]
        dim_reps = [0 for i in range(n_inequiv_shells)]
        T = []
        for ish in range(n_inequiv_shells):
            ll = 2 * corr_shells[inequiv_to_corr[ish]]['l'] + 1
            lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
            T.append(numpy.zeros([lmax, lmax], dtype=complex))

        spin_w90name = ['_up', '_down']
        hamr_full = []
        umat_full = []
        udismat_full = []
        bandmat_full = []

        # TODO: generalise to SP=1 (only partially done)
        rot_mat_time_inv = [0 for i in range(n_corr_shells)]

        # Second, let's read the file containing the Hamiltonian in WF basis
        # produced by Wannier90
        for isp in range(n_spin_blocs):
            # begin loop on isp

            # build filename according to wannier90 conventions
            if SP == 1 and SO == 0:
                mpi.report(
                    "Reading information for spin component n. %d" % isp)
                file_seed = self.w90_seed + spin_w90name[isp]
            else:
                file_seed = self.w90_seed
            # now grab the data from the H(R) file
            mpi.report(
                "\nThe Hamiltonian in MLWF basis is extracted from %s files..." % file_seed)
            nr = rvec = rdeg = nw = hamr = u_mat = udis_mat = band_mat = k_mesh_from_umat = None
            if (mpi.is_master_node()):
                (nr, rvec, rdeg, nw, hamr, u_mat, udis_mat,
                 band_mat, k_mesh_from_umat) = self.read_wannier90data(file_seed, kmesh_mode)
            mpi.barrier()
            nr = mpi.bcast(nr)
            rvec = mpi.bcast(rvec)
            rdeg = mpi.bcast(rdeg)
            nw = mpi.bcast(nw)
            hamr = mpi.bcast(hamr)
            u_mat = mpi.bcast(u_mat)
            udis_mat = mpi.bcast(udis_mat)
            band_mat = mpi.bcast(band_mat)
            k_mesh_from_umat = mpi.bcast(k_mesh_from_umat)
            # number of R vectors, their indices, their degeneracy, number of WFs, H(R),
            # U matrices, U(dis) matrices, band energies, k_mesh of U matrices
            mpi.report('\n... done: {} R vectors, {} WFs found'.format(nr, nw))

            if self.add_lambda:
                mpi.report('Adding local spin-orbit term to Hamiltonian (assuming dxz, dyz, dxy as orbital order)')
                # upscaling quantities
                nw *= 2
                # scale Hamiltonian by 2 to account for spin DOF
                hamr = [numpy.kron(numpy.eye(2), hamr[ir]) for ir in range(nr)]
                # scale lambda matrix by number of correlated shells to account for shells
                # FIXME: does not give the correct order for multiple impurities!
                hamr[nr//2] += numpy.kron(numpy.eye(n_corr_shells), self.lambda_matrix_w90_t2g())
                with numpy.printoptions(linewidth=100, formatter={'complexfloat': '{:+.3f}'.format}):
                    mpi.report('Local Hamiltonian including spin-orbit coupling:')
                    mpi.report(hamr[nr//2])

            if isp == 0:
                # set or check some quantities that must be the same for both
                # spins
                self.nrpt = nr

                # k-point grid: (if not defined before)
                if self.bloch_basis:
                    mpi.report('Reading k mesh from seedname_u.mat file')
                    kpts = k_mesh_from_umat
                    n_k = len(kpts)
                    kpt_weights = numpy.full(n_k, 1/n_k)
                    self.n_k = n_k
                    self.kpts = kpts
                elif kmesh_mode == -1:
                    # the size of the k-point mesh is determined from the
                    # largest R vector
                    nki = [2 * rvec[:, idir].max() + 1 for idir in range(3)]
                    # it will be the same as in the win only when nki is odd, because of the
                    # wannier90 convention: if we have nki k-points along the i-th direction,
                    # then we should get 2*(nki/2)+nki%2 R points along that
                    # direction
                    n_k, kpts, kpt_weights = self.kmesh_build(nki)
                    self.n_k = n_k
                    self.kpts = kpts


                # set the R vectors and their degeneracy
                self.rvec = rvec
                self.rdeg = rdeg

                self.nwfs = nw
                # check that the total number of WFs makes sense
                if self.nwfs < dim_corr_shells:
                    mpi.report(
                        "ERROR: number of WFs in the file smaller than number of correlated orbitals!")
                elif self.nwfs > dim_corr_shells:
                    # NOTE: correlated shells must appear before uncorrelated
                    # ones inside the file
                    mpi.report("Number of WFs larger than correlated orbitals:\n" +
                               "WFs from %d to %d treated as uncorrelated" % (dim_corr_shells + 1, self.nwfs))
                else:
                    mpi.report(
                        "Number of WFs equal to number of correlated orbitals")

                # we assume spin up and spin down always have same total number
                # of WFs
                # get second dimension of udis_mat which corresponds to number of bands in window
                # n_bands_max corresponds to numpy.max(n_orbitals)
                n_bands_max = udis_mat.shape[1] if not self.add_lambda else 2*udis_mat.shape[1]
                n_orbitals = numpy.full([self.n_k, n_spin_blocs], n_bands_max)
            else:
                # consistency check between the _up and _down file contents
                if nr != self.nrpt:
                    mpi.report(
                        "Different number of R vectors for spin-up/spin-down!")
                if nw != self.nwfs:
                    mpi.report(
                        "Different number of WFs for spin-up/spin-down!")

            hamr_full.append(hamr)
            umat_full.append(u_mat)
            udismat_full.append(udis_mat)
            bandmat_full.append(band_mat)

            for ir in range(nr):
                # checks if the Hamiltonian is real (it should, if
                # wannierisation worked fine)
                if numpy.abs((hamr[ir].imag.max()).max()) > self._w90zero:
                    mpi.report(
                        "H(R) has large complex components at R %d" % ir)
                # copy the R=0 block corresponding to the correlated shells
                # into another variable (needed later for finding rot_mat)
                if rvec[ir, 0] == 0 and rvec[ir, 1] == 0 and rvec[ir, 2] == 0:
                    ham_corr0 = hamr[ir][0:dim_corr_shells, 0:dim_corr_shells]

            # checks if ham0 is Hermitian
            if not numpy.allclose(ham_corr0.transpose().conjugate(), ham_corr0, atol=self._w90zero, rtol=0):
                raise ValueError("H(R=0) matrix is not Hermitian!")

            # find rot_mat symmetries by diagonalising the on-site Hamiltonian
            # of the first spin
            if isp == 0:
                use_rotations, rot_mat = self.find_rot_mat(
                    n_corr_shells, corr_shells, shells_map, ham_corr0)
            else:
                # consistency check
                use_rotations_, rot_mat_ = self.find_rot_mat(
                    n_corr_shells, corr_shells, shells_map, ham_corr0)
                if (use_rotations and not use_rotations_):
                    mpi.report(
                        "Rotations cannot be used for spin component n. %d" % isp)
                for icrsh in range(n_corr_shells):
                    if not numpy.allclose(rot_mat_[icrsh], rot_mat[icrsh], atol=self._w90zero, rtol=0):
                        mpi.report(
                            "Rotations for spin component n. %d do not match!" % isp)
        # end loop on isp

        # Reads misc input needed for CSC calculations
        if self.bloch_basis:
            if os.path.isfile(self.w90_seed + '.nscf.out'):
                fermi_weight_file = self.w90_seed + '.nscf.out'
                mpi.report('Reading DFT band occupations from Quantum Espresso output {}'.format(fermi_weight_file))
            elif os.path.isfile('LOCPROJ'):
                assert os.path.isfile('OUTCAR')
                fermi_weight_file = 'LOCPROJ'
                mpi.report('Reading DFT band occupations from Vasp output {}'.format(fermi_weight_file))
            else:
                raise IOError('seedname.nscf.out or LOCPROJ required in bloch_basis mode')

            f_weights = band_window = self.fermi_energy = kpt_basis = None
            if (mpi.is_master_node()):
                f_weights, band_window, self.fermi_energy, kpt_basis = self.convert_misc_input(fermi_weight_file,
                                                                                    self.w90_seed + '.nnkp', n_spin_blocs)
            mpi.barrier()
            f_weights = mpi.bcast(f_weights)
            band_window = mpi.bcast(band_window)
            self.fermi_energy = mpi.bcast(self.fermi_energy)
            kpt_basis = mpi.bcast(kpt_basis)
            # Get density from k-point weighted average and sum over all spins and bands
            density_required = numpy.sum(f_weights.T * kpt_weights) * (2 - SP)
            mpi.report('Overwriting required density with DFT result {:.5f}'.format(density_required))
            mpi.report('and using the DFT Fermi energy {:.5f} eV\n'.format(self.fermi_energy))

        if not self.bloch_basis:
            mpi.report("The k-point grid has dimensions: %d, %d, %d" % tuple(nki))

        # if calculations are spin-polarized, then renormalize k-point weights
        if SP == 1 and SO == 0:
            kpt_weights *= 0.5

        # Third, initialise the projectors
        k_dep_projection = 0   # at the moment not really used, but might get important
        proj_mat = numpy.zeros([self.n_k, n_spin_blocs, n_corr_shells, max(
            [crsh['dim'] for crsh in corr_shells]), numpy.max(n_orbitals)], dtype=complex)
        iorb = 0
        # Projectors are either identity matrix blocks to use with Wannier basis
        # OR correspond to the overlap between Kohn-Sham and Wannier orbitals as
        # P_{nu,alpha](k) = <w_{alpha,k}|psi_{nu,k}>
        # NOTE: we assume that the correlated orbitals appear at the beginning of the H(R)
        # file and that the ordering of MLWFs matches the corr_shell info from
        # the input.
        for isp in range(n_spin_blocs):
            # now combine udismat and umat
            u_total = numpy.einsum('abc,acd->abd',udismat_full[isp],umat_full[isp])
            # transpose and write into proj_mat
            u_temp = numpy.transpose(u_total.conj(),(0,2,1))
            # scale unitary U by 2 to account for spin DOF
            if self.add_lambda: u_temp = numpy.kron(numpy.eye(2), u_temp)
            for icrsh in range(n_corr_shells):
                dim = corr_shells[icrsh]['dim']
                proj_mat[:, isp, icrsh, 0:dim, :] = u_temp[:,iorb:iorb+dim,:]
                iorb += dim

        # Then, compute the hoppings in reciprocal space
        hopping = numpy.zeros([self.n_k, n_spin_blocs, numpy.max(n_orbitals), numpy.max(n_orbitals)], dtype=complex)
        for isp in range(n_spin_blocs):
            # if bloch_basis is True, use Kohn-Sham eigenvalues as hamk
            # this ensures that the calculation of the band-correlation energy
            # is consistent with SumkDFT's calc_density_correction
            if self.bloch_basis:
                # diagonal Kohn-Sham bands
                # TODO: test for system with self.nwfs > dim_corr_shells
                hamk = [numpy.diag(bandmat_full[isp][ik]) for ik in range(self.n_k)]

                # Sanity check if the local Hamiltonian with the projector method
                # corresponds to W90 result
                wannier_ham = self.fourier_ham(hamr_full[isp])
                for ik in range(self.n_k):
                    proj_mat_flattened = numpy.zeros((numpy.max(n_orbitals), numpy.max(n_orbitals)), dtype=complex)
                    iorb = 0
                    for icrsh in range(n_corr_shells):
                        dim = corr_shells[icrsh]['dim']
                        proj_mat_flattened[iorb:iorb+dim,:] = proj_mat[ik, isp][icrsh,0:dim,:].reshape(dim, numpy.max(n_orbitals))
                        iorb += dim
                    downfolded_ham = proj_mat_flattened.dot(hamk[ik].dot(proj_mat_flattened.conj().T))
                    if dim_corr_shells < numpy.max(n_orbitals):
                        downfolded_ham = downfolded_ham[:dim_corr_shells,:dim_corr_shells]
                        wannier_ham[ik] = wannier_ham[ik][:dim_corr_shells,:dim_corr_shells]

                    if not numpy.allclose(downfolded_ham, wannier_ham[ik], atol=1e-4, rtol=0):
                        mpi.report('WARNING: mismatch between downfolded Hamiltonian and '
                                   + f'Fourier transformed H(R). First occurred at kpt {ik}:')

                        with numpy.printoptions(formatter={'complexfloat': '{:+.4f}'.format}):
                            mpi.report('Downfolded Hamiltonian, P H_eig P')
                            mpi.report(downfolded_ham)
                            mpi.report('\nWannier Hamiltonian, Fourier(H(r))')
                            mpi.report(wannier_ham[ik])
                        break
            # else for an isolated set of bands use fourier transform of H(R)
            else:
                # make Fourier transform H(R) -> H(k) : it can be done one spin at a time
                hamk = self.fourier_ham(hamr_full[isp])

                # Sanity check if imaginary diagonal elements are zero, otherwise instabilties in lattice Gf!
                diag_iterator = range(0,dim_corr_shells)
                for ik in range(self.n_k):
                    if not numpy.allclose(hamk[ik][diag_iterator, diag_iterator].imag, 0, atol=1e-10):
                        mpi.report('ERROR: Wannier Hamiltonian has complex diagonal entries. '
                                   + f'First occurred at kpt {ik}:')
                        with numpy.printoptions(formatter={'float': '{:+.10f}'.format}):
                            mpi.report('\nWannier Hamiltonian diagonal, Fourier(H(r)), imaginary')
                            mpi.report(hamk[ik][diag_iterator, diag_iterator].imag)
                        mpi.MPI.COMM_WORLD.Abort(1)
                    # set imaginary part to zero
                    hamk[ik][diag_iterator, diag_iterator] = hamk[ik][diag_iterator, diag_iterator].real + 0*1j
            # finally write hamk into hoppings
            for ik in range(self.n_k):
                hopping[ik, isp] = hamk[ik] - numpy.identity(numpy.max(n_orbitals)) * self.fermi_energy
            hopping *= energy_unit
        mpi.report("Subtracting {:.5f} eV from the Fermi level.".format(self.fermi_energy))

        # bz_weights required by triqs h5 standard but soon to be replaced by kpt_weights
        bz_weights = kpt_weights

        # Finally, save all required data into the HDF archive:
        # use_rotations is supposed to be an int = 0, 1, no bool
        use_rotations = int(use_rotations)
        if mpi.is_master_node():
            with HDFArchive(self.hdf_file, 'a') as ar:
                if not (self.dft_subgrp in ar):
                    ar.create_group(self.dft_subgrp)
                # The subgroup containing the data. If it does not exist, it is
                # created. If it exists, the data is overwritten!
                things_to_save = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                              'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                              'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                              'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr', 'kpt_weights', 'kpts']
                if self.bloch_basis: numpy.append(things_to_save, 'kpt_basis')
                for it in things_to_save:
                    ar[self.dft_subgrp][it] = locals()[it]

                # Store Fermi weights to 'dft_misc_input'
                if not (self.misc_subgrp in ar):
                    ar.create_group(self.misc_subgrp)
                ar[self.misc_subgrp]['dft_fermi_energy'] = self.fermi_energy
                if self.bloch_basis:
                    ar[self.misc_subgrp]['dft_fermi_weights'] = f_weights
                    ar[self.misc_subgrp]['band_window'] = band_window+1 # Change to 1-based index
                    ar[self.misc_subgrp]['kpts_cart'] = numpy.dot(kpts, kpt_basis.T)
        mpi.barrier()

    def read_wannier90data(self, wannier_seed="wannier", kmesh_mode=0):
        """
        Method for reading the seedname_hr.dat file produced by Wannier90 (http://wannier.org)

        Parameters
        ----------
        wannier_seed : string
            seedname to read H(R) file produced by Wannier90 (usually seedname_hr.dat)

        Returns
        -------
        nrpt : integer
            number of R vectors found in the file
        rvec_idx : numpy.array of integers
            Miller indices of the R vectors
        rvec_deg : numpy.array of floats
            weight of the R vectors
        num_wf : integer
            number of Wannier functions found
        h_of_r : list of numpy.array
            <w_i|H(R)|w_j> = Hamilonian matrix elements in the Wannier basis
        u_mat : numpy.array
            U_mn^k = unitary matrix elements which mix the Kohn-Sham states
        udis_mat : numpy.array
            U^dis(k) = rectangular matrix for entangled bands
        band_mat : numpy.array
            \epsilon_nk = Kohn-Sham eigenvalues (in eV) needed for entangled bands
            If not self.bloch_basis unused and therefore None
        k_mesh : numpy.array
            The k mesh read from the seedname_u.mat file to ensure consistency
            If not self.bloch_basis unused and therefore None
        """

        hr_filename = wannier_seed + '_hr.dat'
        try:
            with open(hr_filename, 'r') as hr_filedesc:
                hr_data = hr_filedesc.readlines()
        except IOError:
            mpi.report("The file %s could not be read!" % hr_filename)

        mpi.report('reading {:20}...{}'.format(hr_filename,hr_data[0].strip('\n')))

        try:
            # reads number of Wannier functions per spin
            num_wf = int(hr_data[1])
            nrpt = int(hr_data[2])
        except ValueError:
            mpi.report("Could not read number of WFs or R vectors")

        k_mesh = None

        if not self.bloch_basis:
            # For kmesh_mode == -1, size of automatic k mesh known when rvec_idx
            # have been read. For kmesh_mode >= 0, kpts have been determined already
            if kmesh_mode >= 0:
                n_k = self.n_k
        else:
            # first, read u matrices from 'seedname_u.mat'
            u_filename = wannier_seed + '_u.mat'
            with open(u_filename,'r') as u_file:
                u_data = u_file.readlines()
            # reads number of kpoints and number of wannier functions
            n_k, num_wf_u, _ = map(int, u_data[1].split())
            assert num_wf_u == num_wf, '#WFs must be identical for *_u.mat and *_hr.dat'
            mpi.report('reading {:20}...{}'.format(u_filename,u_data[0].strip('\n')))

            # Reads k mesh from all lines with 3 floats
            k_mesh = numpy.loadtxt((line for line in u_data if line.count('.') == 3))
            assert k_mesh.shape == (n_k, 3)
            # Reads u matrices from all lines with 2 floats
            u_mat = numpy.loadtxt((line for line in u_data if line.count('.') == 2))
            assert u_mat.shape == (n_k*num_wf*num_wf, 2)

            mpi.report('Writing h5 archive in projector formalism: H(k) defined in KS Bloch basis')

            try:
                # read 'seedname_u_dis.mat'
                udis_filename = wannier_seed + '_u_dis.mat'
                # if it exists the Kohn-Sham eigenvalues and the window are needed
                band_filename = wannier_seed + '.eig'
                wout_filename = wannier_seed + '.wout'

                with open(udis_filename,'r') as udis_file:
                    udis_data = udis_file.readlines()
                disentangle = True
            except IOError:
                disentangle = False
                mpi.report('WARNING: File {} missing.'.format(udis_filename))
                mpi.report('Assuming an isolated set of bands. Check if this is what you want!')

            # read Kohn-Sham eigenvalues from 'seedname.eig'
            mpi.report('Reading {}'.format(band_filename))
            band_data = numpy.loadtxt(band_filename, usecols=2)

            if disentangle:
                # reads number of kpoints, number of wannier functions and bands
                num_k_udis, num_wf_udis, num_ks_bands = map(int, udis_data[1].split())

                assert num_k_udis == n_k, '#k points must be identical for *.inp and *_u_dis.mat'
                assert num_wf_udis == num_wf, '#WFs must be identical for *_u.mat and *_hr.dat'

                mpi.report('Found {:22}...{}, '.format(udis_filename,udis_data[0].strip('\n')))
                udis_data = numpy.loadtxt(udis_data, usecols=(0, 1), skiprows=2)

                # read disentanglement window from 'seedname.wout'
                with open(wout_filename) as wout_file:
                    for line in wout_file:
                        if 'Outer:' in line:
                            content = line.split()
                            index = content.index('Outer:') + 1
                            dis_window_min = float(content[index])
                            dis_window_max = float(content[index+2])
                            break
                mpi.report('and {} for disentanglement energy window.'.format(wout_filename))
            else:
                num_ks_bands = num_wf

        # allocate arrays to save the R vector indexes and degeneracies and the
        # Hamiltonian
        rvec_idx = numpy.zeros((nrpt, 3), dtype=int)
        rvec_deg = numpy.zeros(nrpt, dtype=int)
        h_of_r = [numpy.zeros((num_wf, num_wf), dtype=complex)
                  for n in range(nrpt)]

        # variable currpos points to the current line in the file
        currpos = 2
        try:
            ir = 0
            # read the degeneracy of the R vectors (needed for the Fourier
            # transform)
            while ir < nrpt:
                currpos += 1
                for x in hr_data[currpos].split():
                    if ir >= nrpt:
                        raise IndexError("wrong number of R vectors??")
                    rvec_deg[ir] = int(x)
                    ir += 1
            # for each direct lattice vector R read the block of the
            # Hamiltonian H(R)
            for ir, jj, ii in product(range(nrpt), range(num_wf), range(num_wf)):
                # advance one line, split the line into tokens
                currpos += 1
                cline = hr_data[currpos].split()
                # check if the orbital indexes in the file make sense
                if int(cline[3]) != ii + 1 or int(cline[4]) != jj + 1:
                    mpi.report(
                        "Inconsistent indices at %s%s of R n. %s" % (ii, jj, ir))
                rcurr = numpy.array([int(cline[0]), int(cline[1]), int(cline[2])])
                if ii == 0 and jj == 0:
                    rvec_idx[ir] = rcurr
                    rprec = rcurr
                else:
                    # check if the vector indices are consistent
                    if not numpy.array_equal(rcurr, rprec):
                        mpi.report(
                            "Inconsistent indices for R vector n. %s" % ir)

                # fill h_of_r with the matrix elements of the Hamiltonian
                h_of_r[ir][ii, jj] = complex(float(cline[5]), float(cline[6]))

        except ValueError:
            mpi.report("Wrong data or structure in file %s" % hr_filename)

        # first, get the input for u_mat
        if self.bloch_basis:
            u_mat = u_mat[:, 0] + 1j * u_mat[:, 1]
            u_mat = u_mat.reshape((n_k, num_wf, num_wf)).transpose((0, 2, 1))
        else:
            if kmesh_mode == -1:
                n_k = numpy.prod([2 * rvec_idx[:, idir].max() + 1 for idir in range(3)])
            # Wannier basis; fill u_mat with identity
            u_mat = numpy.zeros([n_k, num_wf, num_wf], dtype=complex)
            for ik in range(n_k):
                u_mat[ik,:,:] = numpy.identity(num_wf,dtype=complex)

        # now, check what is needed in the case of disentanglement:
        # The file seedname_u_dis.mat contains only the bands inside the window
        # and then fills the rest up with zeros. Therefore, we need to put the
        # entries from udis_data in the correct position in udis_mat, i.e.
        # shifting by the number of bands below dis_window_min
        if self.bloch_basis:
            # reshape band_data
            band_mat = band_data.reshape(n_k, num_ks_bands)
        else:
            # Not used in wannier basis
            band_mat = None

        if self.bloch_basis and disentangle:
            # Determine which bands are inside the band window
            inside_window = numpy.logical_and(band_mat >= dis_window_min,
                                              band_mat <= dis_window_max)
            n_inside_per_k = numpy.sum(inside_window, axis=1)

            # Reformats udis_data as complex, without header
            udis_data = udis_data[:, 0] + 1j * udis_data[:, 1]
            udis_data = udis_data.reshape((n_k, num_wf*num_ks_bands+1))[:, 1:]
            udis_data = udis_data.reshape((n_k, num_wf, num_ks_bands))

            #initiate U disentanglement matrices and fill from file "seedname_u_dis.mat"
            udis_mat = numpy.zeros([n_k, num_ks_bands, num_wf], dtype=complex)
            for ik in range(n_k):
                udis_mat[ik, inside_window[ik]] = udis_data[ik, :, :n_inside_per_k[ik]].T
                if not numpy.allclose(udis_data[ik, :, n_inside_per_k[ik]:], 0):
                    raise ValueError('This error could come from rounding of the band window in the seedname.wout. '
                                     + 'Never use default outer window but something wider and '
                                     + 'check that your outer window is not close to any band energy.')
        else:
            # no disentanglement; fill udis_mat with identity
            udis_mat = numpy.array([numpy.identity(num_wf,dtype=complex)] * n_k)

        # return the data into variables
        return nrpt, rvec_idx, rvec_deg, num_wf, h_of_r, u_mat, udis_mat, band_mat, k_mesh

    def find_rot_mat(self, n_sh, sh_lst, sh_map, ham0):
        """
        Method for finding the matrices that bring from local to global coordinate systems
        (and viceversa), based on the eigenvalues of H(R=0)

        Parameters
        ----------
        n_sh : integer
            number of shells
        sh_lst : list of shells-type dictionaries
            contains the shells (could be correlated or not)
        sh_map : list of integers
            mapping between shells
        ham0 : numpy.array of floats
            local Hamiltonian matrix elements

        Returns
        -------
        succeeded : integer
            if 0, something failed in the construction of the matrices
        rot_mat : list of numpy.array
            rotation matrix for each of the shell

        """

        # initialize the rotation matrices to identities
        rot_mat = [numpy.identity(sh_lst[ish]['dim'], dtype=complex)
                   for ish in range(n_sh)]
        succeeded = True

        hs = ham0.shape
        if hs[0] != hs[1] or hs[0] != sum([sh['dim'] for sh in sh_lst]):
            mpi.report(
                "find_rot_mat: wrong block structure of input Hamiltonian!")
            # this error will lead into troubles later... early return
            succeeded = False
            return succeeded, rot_mat

        # Method none as physically unsound option for testing
        # Returns identity matrices as rotation matrices
        if self.rot_mat_type == 'none':
            mpi.report('WARNING: using the method "none" leads to physically wrong results. '
                       + 'Only use for testing if other methods fail.')
            succeeded = True
            return succeeded, rot_mat

        # TODO: better handling of degenerate eigenvalue case
        eigval_lst = [None] * n_sh
        eigvec_lst = [None] * n_sh
        ham0_lst = [None] * n_sh
        iwf = 0
        # loop over shells
        for ish in range(n_sh):
            # nw = number of orbitals in this shell
            nw = sh_lst[ish]["dim"]
            # save the sub-block of H(0) corresponding to this shell
            ham0_lst[ish] = ham0[iwf:iwf+nw, iwf:iwf+nw]
            # diagonalize the sub-block for this shell
            eigval, eigvec = numpy.linalg.eigh(ham0_lst[ish])
            eigval_lst[ish] = eigval
            eigvec_lst[ish] = eigvec
            iwf += nw
            # TODO: better handling of degenerate eigenvalue case
            if sh_map[ish] != ish:  # issue warning only when there are equivalent shells
                for i in range(nw):
                    for j in range(i + 1, nw):
                        if abs(eigval[j] - eigval[i]) < self._w90zero:
                            mpi.report("WARNING: degenerate eigenvalue of H(0) detected for shell %d: " % (ish) +
                                       "global-to-local transformation might not work!")

        for ish in range(n_sh):
            try:
                # build rotation matrices either...
                if self.rot_mat_type == 'hloc_diag':
                    # using the unitary transformations that diagonalize H(0)
                    rot_mat[ish] = eigvec_lst[ish]
                elif self.rot_mat_type == 'wannier':
                    # or by combining those transformations (i.e. for each group,
                    # the representative site is chosen as the global frame of reference)
                    rot_mat[ish] = numpy.dot(eigvec_lst[ish],
                                             eigvec_lst[sh_map[ish]].conjugate().transpose())
            except ValueError:
                mpi.report(
                    "Global-to-local rotation matrices cannot be constructed!")

            # check that eigenvalues are the same (within accuracy) for
            # equivalent shells
            if not numpy.allclose(eigval_lst[ish], eigval_lst[sh_map[ish]],
                    atol=self._w90zero, rtol=0):
                mpi.report(f'ERROR: eigenvalue mismatch between equivalent shells! {ish:d}')
                eigval_diff = eigval_lst[ish] - eigval_lst[sh_map[ish]]
                mpi.report(f'Eigenvalue difference {eigval_diff}, but threshold set to {self._w90zero:.1e}.')
                mpi.report('Consider lowering threshold if you are certain the mapping is correct.')
                succeeded = False

            # check that rotation matrices are unitary
            # nw = number of orbitals in this shell
            nw = sh_lst[ish]["dim"]
            tmp_mat = numpy.dot(rot_mat[ish],rot_mat[ish].conjugate().transpose())
            if not numpy.allclose(tmp_mat, numpy.identity(nw),
                                  atol=self._w90zero, rtol=0):
                mpi.report(f'ERROR: rot_mat for shell {ish:d} is not unitary!')
                succeeded = False

            # check that rotation matrices map equivalent H(0) blocks as they should
            # (assuming representative shell as global frame of reference)
            if self.rot_mat_type == 'hloc_diag':
                tmp_mat = numpy.dot( rot_mat[ish],
                        rot_mat[sh_map[ish]].conjugate().transpose() )
            elif self.rot_mat_type == 'wannier':
                tmp_mat = rot_mat[ish]
            tmp_mat = numpy.dot(tmp_mat.conjugate().transpose(),
                    numpy.dot(ham0_lst[ish],tmp_mat))
            if not numpy.allclose(tmp_mat, ham0_lst[sh_map[ish]],
                                  atol=self._w90zero, rtol=0):
                mpi.report(f'ERROR: rot_mat does not map H(0) correctly! {ish:d}')
                succeeded = False

        # abort in case the rot_mat was not found correctly to avoid the user to ignore ERRORS
        if not succeeded: mpi.MPI.COMM_WORLD.Abort(1)

        return succeeded, rot_mat

    def kmesh_build(self, msize, mmode=0):
        """
        Method for the generation of the k-point mesh.
        Right now it only supports the option for generating a full grid containing k=0,0,0.

        Parameters
        ----------
        msize : list of 3 integers
            the dimensions of the mesh
        mmode : integer
            mesh generation mode (right now, only full grid available)

        Returns
        -------
        nkpt : integer
            total number of k-points in the mesh
        kpts : numpy.array[nkpt,3] of floats
            the coordinates of all k-points
        wk : numpy.array[nkpt] of floats
            the weight of each k-point

        """

        if mmode != 0:
            raise ValueError("Mesh generation mode not supported: %s" % mmode)

        assert len(msize) == 3

        # a regular mesh including Gamma point
        # total number of k-points
        nkpt = numpy.prod(msize)
        kpts = numpy.array(list(product(range(msize[0]), range(msize[1]), range(msize[2]))))
        kpts = kpts / numpy.array(msize)
        # weight is equal for all k-points because wannier90 uses uniform grid on whole BZ
        # (normalization is always 1 and takes into account spin degeneracy)
        wk = numpy.full(nkpt, 1/nkpt)

        return nkpt, kpts, wk

    def fourier_ham(self, h_of_r):
        """
        Method for obtaining H(k) from H(R) via Fourier transform
        The R vectors and k-point mesh are read from global module variables

        Parameters
        ----------
        h_of_r : list of numpy.array[norb,norb]
            Hamiltonian H(R) in Wannier basis

        Returns
        -------
        h_of_k : list of numpy.array[norb,norb]
            transformed Hamiltonian H(k) in Wannier basis

        """

        h_of_k = [numpy.zeros((self.nwfs, self.nwfs), dtype=complex)
                  for ik in range(self.n_k)]
        h_of_k_array = numpy.array(h_of_k, dtype=complex)
        ridx = numpy.array(list(range(self.nrpt)))

        for ir in mpi.slice_array(ridx):
            for ik in list(range(self.n_k)):
                rdotk = numpy.dot(self.kpts[ik], self.rvec[ir])
                factor = numpy.exp(2j * numpy.pi * rdotk) / self.rdeg[ir]
                h_of_k_array[ik, :, :] += factor * h_of_r[ir][:, :]
        h_of_k_array = mpi.all_reduce(mpi.world, h_of_k_array, lambda x, y: x + y)
        mpi.barrier()

        h_of_k = list(h_of_k_array[:])

        return h_of_k

    def convert_misc_input(self, out_filename, nnkp_filename, n_spin_blocs):
        """
        Reads input from DFT code calculations to get occupations

        Parameters
        ----------
        out_filename : string
            filename of DFT output file containing occupation data
        nnkp_file : string
            filename of Wannier postproc_setup run
        n_spin_blocs : int
            SP + 1 - SO

        Returns
        -------
        fermi_weights : numpy.array[self.n_k, n_spin_blocs ,n_orbitals]
            occupations from DFT calculation
        band_window : numpy.array[self.n_k, n_spin_blocs ,n_orbitals]
            band indices of correlated subspace
        fermi_energy : float
            the DFT Kohn-Sham Fermi energy
        kpt_basis: numpy.array[3, 3]
            the basis vectors in reciprocal space
        """

        assert n_spin_blocs == 1, 'spin-polarized not implemented'
        assert 'nscf.out' in out_filename or out_filename == 'LOCPROJ'

        occupations = []
        reading_kpt_basis = False
        lines_read_kpt_basis = 0
        kpt_basis = numpy.zeros((3, 3))

        if 'nscf.out' in out_filename:
            occupations = []
            with open(out_filename,'r') as out_file:
                out_data = out_file.readlines()
            # Reads number of Kohn-Sham states and Fermi energy
            for line in out_data:
                if 'number of Kohn-Sham states' in line:
                    n_ks = int(line.split()[-1])
                elif 'Fermi energy' in line:
                    fermi_energy = float(line.split()[-2])
                elif 'reciprocal axes' in line:
                    reading_kpt_basis = True
                    continue
                elif reading_kpt_basis and lines_read_kpt_basis < 3:
                    kpt_basis[lines_read_kpt_basis, :] = line.split()[3:6]
                    lines_read_kpt_basis +=1

                # get occupations
            for ct, line in enumerate(out_data):
                if line.strip() == 'End of band structure calculation':
                    break

            assert 'k =' in out_data[ct + 2], 'Cannot read occupations. Set verbosity = "high" in {}'.format(out_filename)
            out_data = out_data[ct+2:]

            # block size of eigenvalues + occupations per k-point
            n_block = int(2*numpy.ceil(n_ks/8)+5)

            for ik in range(self.n_k):
                # get data
                k_block = [line.split() for line in out_data[ik*n_block+2:ik*n_block+n_block-1]]
                # second half corresponds to occupations
                occs = k_block[int(len(k_block)/2)+1:]
                flattened_occs = [float(item) for sublist in occs for item in sublist]
                occupations.append(flattened_occs)
        else:
            # Reads LOCPROJ
            with open(out_filename, 'r') as file:
                header = file.readline()
                n_ks = int(header.split()[2])
                fermi_energy = float(header.split()[4])

                occupations = numpy.loadtxt((line for line in file if 'orbital' in line), usecols=5)
            occupations = occupations.reshape((self.n_k, n_ks))

            # Read reciprocal vectors from OUTCAR
            with open('OUTCAR', 'r') as file:
                for line in file:
                    if 'reciprocal lattice vectors' in line:
                        reading_kpt_basis = True
                    elif reading_kpt_basis:
                        kpt_basis[lines_read_kpt_basis, :] = line.split()[3:6]
                        lines_read_kpt_basis += 1
                        if lines_read_kpt_basis == 3:
                            break

        # assume that all bands contribute, then remove from exclude_bands; python indexing
        corr_bands = list(range(n_ks))
        # read exclude_bands from "seedname.nnkp" file
        with open(nnkp_filename, 'r') as nnkp_file:
            read = False
            skip = False
            for line in nnkp_file:
                if line.strip() == 'begin exclude_bands':
                    read = True
                    # skip one more line that contains total number of excluded bands
                    skip = True
                    continue
                elif line.strip() == 'end exclude_bands':
                    read = False
                    continue
                elif skip:
                    skip = False
                    continue
                elif read:
                    # wannier index -1
                    corr_bands.remove(int(line)-1)

        # For now, it is only supported to exclude the lowest and highest bands
        # If bands in the middle are supposed to be excluded, this doesn't work with the band_window
        #     We'd need to manually add rows of zeros to the projectors
        if numpy.any(numpy.diff(corr_bands) != 1):
            raise NotImplementedError('Can only exclude the lowest or highest bands')
        band_window = numpy.array([[(min(corr_bands), max(corr_bands))]*self.n_k]*n_spin_blocs)

        included_occupations = numpy.array(occupations)[:, corr_bands]
        # Adds spin dimension without changing the array
        f_weights = included_occupations.reshape(included_occupations.shape[0], 1,
                                                 included_occupations.shape[1])

        return f_weights, band_window, fermi_energy, kpt_basis

    def lambda_matrix_w90_t2g(self):
        """
        Adds local spin-orbit interaction term to the t2g subspace. Orbital order is assumed to be
        xz_up, yz_up, xy_up, xz_dn, yz_dn, xy_dn as given by Wannier90 per default.
        Parameters are defined as self.add_lambda = [lambda_x, lambda_y, lambda_z], representative of
        the orbital coupling terms perpendicular to [x, y, z] i.e. [d_yz, d_xz, d_xy], respectively.

        Returns
        -------
        lambda_matrix : numpy.array[6, 6]
            local spin-orbit term to be added to H(0)

        """

        lambda_x, lambda_y, lambda_z = self.add_lambda

        lambda_matrix = numpy.zeros((6,6), dtype=complex)
        lambda_matrix[0,1] = -1j*lambda_z/2.0
        lambda_matrix[0,5] =  1j*lambda_x/2.0
        lambda_matrix[1,5] =    -lambda_y/2.0
        lambda_matrix[2,3] = -1j*lambda_x/2.0
        lambda_matrix[2,4] =     lambda_y/2.0
        lambda_matrix[3,4] =  1j*lambda_z/2.0
        lambda_matrix += numpy.transpose(numpy.conjugate(lambda_matrix))

        return lambda_matrix
