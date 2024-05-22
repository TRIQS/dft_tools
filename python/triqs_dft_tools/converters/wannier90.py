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

"""
Wannier90 to HDF5 converter for TRIQS/dft_tools

Written by Gabriele Sclauzero (Materials Theory, ETH Zurich), Dec 2015 -- Jan 2016,
updated by Maximilian Merkel (Materials Theory, ETH Zurich), Aug 2020 -- Feb 2022,
and by Sophie Beck (Materials Theory, ETH Zurich), Sep 2020 -- Apr 2021,
under the supervision of Claude Ederer (Materials Theory).
Partially based on previous work by K. Dymkovski and the TRIQS/dft_tools team.

Limitations of the current implementation:
- the T rotation matrices are not used in this implementation

Things to be improved/checked:
- the case with SP=1 is only half implemented and never tested (do we need to
  define rot_mat_time_inv also if symm_op = 0?)
- the calculation of rot_mat in find_rot_mat() relies on the eigenvalues of H(0);
  this might fail in presence of degenerate eigenvalues (now just prints warning)
- make the code more MPI safe (error handling): if we run with more than one process
  and an error occurs on the masternode, the calculation does not abort
- in case of disentanglement, the outer window being close to Kohn-Sham energies
  can cause a problem in creating the udis_mat_spin in read_wannier90data
- bloch_basis on SO coupled calculations has never been tested but might work
- would be helpful to read the order of orbitals from the nnkp or wout file
  and save it to, e.g., misc_subgrp for codes working on the generated h5
"""

import os.path
from itertools import product
import numpy as np

from h5 import HDFArchive
from triqs.utility import mpi
from .converter_tools import ConverterTools

class Wannier90Converter(ConverterTools):
    """
    Conversion from Wannier90 output to an hdf5 file that can be used as input
    for the SumkDFT class.
    """

    def __init__(self, seedname, hdf_filename=None, dft_subgrp='dft_input',
                 symmcorr_subgrp='dft_symmcorr_input', misc_subgrp='dft_misc_input',
                 repacking=False, rot_mat_type='hloc_diag', bloch_basis=False, add_lambda=None,
                 w90zero=2e-6, reorder_orbital_and_spin_vasp5=False):
        r"""
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
            Add local spin-orbit term
        w90zero : float, optional
            Threshold on symmetry checks of Hamiltonian and rot_mat
        reorder_orbital_and_spin_vasp5 : bool, optional
            Is false for output from VASP 6 and Quantum Espresso
            Reorder orbitals and spins from the VASP 5 convention of first all
            orbitals with up and then all orbitals with down to the "usual"
            convention of every up orbital immediately being followed by its
            corresponding down orbital
        """
        self._name = 'Wannier90Converter'
        assert isinstance(seedname, str), self._name + \
            ': Please provide the DFT files\' base name as a string.'
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
        self.reorder_orbital_and_spin_vasp5 = reorder_orbital_and_spin_vasp5
        if self.add_lambda is not None and len(self.add_lambda) != 3:
            raise ValueError('If specifying add_lambda, give three values.')
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
        - misc_subgrp

        in the hdf5 archive.
        """
        # Reads in inp file
        input_params = None
        if mpi.is_master_node():
            input_params = read_input_file(self.inp_file, self.fortran_to_replace)
        (kmesh_mode, kmesh_size, density_required, n_corr_shells,
         corr_shells, n_shells, shells, fermi_energy) = mpi.bcast(input_params)
        if density_required is None and not self.bloch_basis:
            raise ValueError('Required density necessary if not in bloch basis')

        shells, corr_shells, SO, SP, n_spin_blocks = check_and_adapt_for_soc(shells, corr_shells,
                                                                             self.bloch_basis, self.add_lambda)

        # Determines sum of dimensions of all impurities
        dim_corr_shells = sum([sh['dim'] for sh in corr_shells])
        mpi.report('Total number of WFs expected in the correlated shells: {}'.format(dim_corr_shells))

        # determine the number of inequivalent correlated shells and maps,
        # needed for further processing
        n_inequiv_shells, corr_to_inequiv, inequiv_to_corr = self.det_shell_equivalence(corr_shells)
        mpi.report('Number of inequivalent shells: {}'.format(n_inequiv_shells))
        mpi.report('Shell representatives: {}'.format(inequiv_to_corr))
        shells_map = [inequiv_to_corr[corr_to_inequiv[icrsh]]
                      for icrsh in range(n_corr_shells)]
        mpi.report('Mapping: {}'.format(shells_map))

        # Second, let's read the file containing the Hamiltonian in WF basis
        # produced by Wannier90
        (wannier_hr, u_total, ks_eigenvals, r_vector, r_degeneracy, n_wannier, n_bands,
         k_mesh_from_umat, wan_centres) = read_all_wannier90_data(n_spin_blocks, dim_corr_shells,
                                                                  self.w90_seed, self.add_lambda,
                                                                  self.bloch_basis)

        # Read high-symmetry k-path from _band.kpt
        w90_kpath_results = None
        if mpi.is_master_node():
            w90_kpath_results = read_wannier90_symm_kpath(self.w90_seed)
        symm_kpath_info = mpi.bcast(w90_kpath_results)

        # Builds kmesh or uses kmesh from _u.mat
        if self.bloch_basis:
            # If output is in bloch_basis, we use k mesh from seedname_u.mat for consistency
            kpts = k_mesh_from_umat
            n_k = len(kpts)
            kpt_weights = np.full(n_k, 1/n_k)
            mpi.report('Using k mesh from seedname_u.mat with {} k points.'.format(n_k))
        else:
            if kmesh_mode == -1:
                # The size of the k-point mesh is determined from the largest R vector.
                # It will only be the same as in the win when kmesh_size is odd, because of the
                # wannier90 convention: if we have kmesh_size k-points along the i-th direction,
                # then we should get 2*(kmesh_size/2)+kmesh_size%2 R points along that direction
                kmesh_size = [2 * r_vector[:, idir].max() + 1 for idir in range(3)]
            mpi.report('Building k-point grid with dimension {} x {} x {}.'.format(*kmesh_size))
            n_k, kpts, kpt_weights = build_kmesh(kmesh_size, kmesh_mode if kmesh_mode >=0 else 0)

        # Reads misc input needed for CSC calculations
        if self.bloch_basis:
            misc_results = None
            if mpi.is_master_node():
                misc_results = read_misc_input(self.w90_seed, n_spin_blocks, n_k)
            f_weights, band_window, fermi_energy, kpt_basis = mpi.bcast(misc_results)

            # Get density from k-point weighted average and sum over all spins and bands
            density_required = np.sum(f_weights.T * kpt_weights) * (2 - SP)
            mpi.report('Overwriting required density with DFT result {:.5f}'.format(density_required))
            mpi.report('and using the DFT Fermi energy {:.5f} eV\n'.format(fermi_energy))

        # Switches from Vasp 5 + wannier90 order to triqs order
        if SO and not self.add_lambda and self.reorder_orbital_and_spin_vasp5:
            wannier_hr, u_total = reorder_orbital_and_spin(n_wannier, wannier_hr, u_total)

        # Finds R=0 index
        r_zero_index = np.nonzero(np.all(r_vector == 0, axis=1))[0][0]

        # Increases Hamiltonian size and adds lambda term
        if self.add_lambda:
            # scale Hamiltonian by 2 to account for spin DOF
            wannier_hr = np.array([[np.kron(hr, np.eye(2)) for hr in wannier_hr[0]]])
            wannier_hr[0, r_zero_index] += generate_local_so_matrix_t2g(self.add_lambda, n_corr_shells, n_wannier)

            with np.printoptions(linewidth=100, formatter={'complexfloat': '{:+.3f}'.format}):
                mpi.report('Local Hamiltonian including spin-orbit coupling:')
                mpi.report(wannier_hr[0, r_zero_index])

        # Runs tests in H(R)
        check_hr(wannier_hr, self._w90zero, r_zero_index)

        # Determines rot_mat
        wannier_hr0 = wannier_hr[:, r_zero_index]
        rot_mat = find_rot_mat(n_corr_shells, corr_shells, shells_map, wannier_hr0,
                               self.rot_mat_type, self._w90zero, n_spin_blocks)
        if rot_mat is None:
            raise ValueError('Something went wrong in the creation of rotation matrices.')

        # Renormalizes k-point weights if calculations are spin-polarized
        if SP == 1 and SO == 0:
            kpt_weights *= 0.5

        # Sets the projectors
        # Projectors are either identity matrix blocks to use with Wannier basis
        # OR correspond to the overlap between Kohn-Sham and Wannier orbitals as
        # P_{nu,alpha](k) = <w_{alpha,k}|psi_{nu,k}>
        # NOTE: we assume that the correlated orbitals appear at the beginning of the H(R)
        # file and that the ordering of MLWFs matches the corr_shell info from
        # the input.
        proj_mat = np.zeros([n_k, n_spin_blocks, n_corr_shells,
                             max(crsh['dim'] for crsh in corr_shells), n_bands], dtype=complex)
        if not self.bloch_basis:
            u_total = np.array([[np.identity(n_wannier)] * n_k] * n_spin_blocks)
        for isp in range(n_spin_blocks):
            iorb = 0
            for icrsh in range(n_corr_shells):
                dim = corr_shells[icrsh]['dim']
                proj_mat[:, isp, icrsh, :dim, :] = u_total[isp, :, iorb:iorb+dim, :]
                iorb += dim

        # Then, compute the hoppings in reciprocal space
        wannier_hk = fourier_transform_hamiltonian(wannier_hr, r_vector, r_degeneracy, kpts)

        diag_iterator = range(n_bands)
        if self.bloch_basis:
            # if bloch_basis is True, use Kohn-Sham eigenvalues as hamk
            # this ensures that the calculation of the band-correlation energy
            # is consistent with SumkDFT's calc_density_correction
            hopping = np.zeros([n_k, n_spin_blocks, n_bands, n_bands], dtype=complex)
            hopping[:, :, diag_iterator, diag_iterator] = ks_eigenvals.transpose((1, 0, 2))
            check_bloch_basis_hk(n_corr_shells, corr_shells, n_k, n_spin_blocks, n_bands,
                                 proj_mat, dim_corr_shells, wannier_hk, hopping)
        else:
            hopping = wannier_hk.transpose((1, 0, 2, 3))
            check_wannier_basis_hk(hopping, dim_corr_shells)
            hopping[:, :, diag_iterator, diag_iterator] = hopping[:, :, diag_iterator, diag_iterator].real + 0j

        hopping[:, :, diag_iterator, diag_iterator] -= fermi_energy
        mpi.report('Subtracting {:.5f} eV from the Fermi level.'.format(fermi_energy))

        # Prints Hamiltonian at first k point
        for icrsh in range(n_corr_shells):
            ik = 0
            isp = 0
            dim = corr_shells[icrsh]['dim']
            hamiltonian = np.einsum('ji,jk,kl,ml,mn->in', rot_mat[icrsh].conj(), proj_mat[ik][isp][icrsh, :dim],
                                    hopping[ik, isp], proj_mat[ik][isp][icrsh, :dim].conj(), rot_mat[icrsh])
            with np.printoptions(floatmode='fixed', precision=3, suppress=True, linewidth=100):
                mpi.report(f'Hamiltonian at first k point for corr. shell {icrsh}', hamiltonian)

        # Sets variables that are the same for every input
        symm_op = 0 # Wannier90 does not use symmetries to reduce the k-points
        charge_below = 0 # total charge below energy window NOT used for now
        energy_unit = 1.0 # should be understood as eV units
        hopping *= energy_unit
        # not used in this version: reset to dummy values?
        n_reps = [1 for i in range(n_inequiv_shells)]
        dim_reps = [0 for i in range(n_inequiv_shells)]
        T = []
        for ish in range(n_inequiv_shells):
            ll = 2 * corr_shells[inequiv_to_corr[ish]]['l'] + 1
            lmax = ll * (corr_shells[inequiv_to_corr[ish]]['SO'] + 1)
            T.append(np.zeros([lmax, lmax], dtype=complex))
        # TODO: generalise to SP=1 (only partially done)
        rot_mat_time_inv = [0 for i in range(n_corr_shells)]
        k_dep_projection = 0   # at the moment not really used, but might get important
        use_rotations = 1

        # bz_weights required by triqs h5 standard but soon to be replaced by kpt_weights
        bz_weights = kpt_weights
        # n_orbitals required by triqs h5 standard, which actually contains the number of bands
        n_orbitals = np.full((n_k, n_spin_blocks), n_bands)

        #new variable: dft_code - this determines which DFT code the inputs come from.
        #used for certain routines within dft_tools if treating the inputs differently is required.
        dft_code = 'w90'

        # Finally, save all required data into the HDF archive:
        if mpi.is_master_node():
            with HDFArchive(self.hdf_file, 'a') as archive:
                if self.dft_subgrp not in archive:
                    archive.create_group(self.dft_subgrp)
                # The subgroup containing the data. If it does not exist, it is
                # created. If it exists, the data is overwritten!
                things_to_save = ['energy_unit', 'n_k', 'k_dep_projection', 'SP', 'SO', 'charge_below', 'density_required',
                              'symm_op', 'n_shells', 'shells', 'n_corr_shells', 'corr_shells', 'use_rotations', 'rot_mat',
                              'rot_mat_time_inv', 'n_reps', 'dim_reps', 'T', 'n_orbitals', 'proj_mat', 'bz_weights', 'hopping',
                              'n_inequiv_shells', 'corr_to_inequiv', 'inequiv_to_corr', 'kpt_weights', 'kpts', 'dft_code']
                if wan_centres is not None:
                    things_to_save.append('wan_centres')
                if self.bloch_basis:
                    things_to_save.append('kpt_basis')
                for it in things_to_save:
                    archive[self.dft_subgrp][it] = locals()[it]

                # Store Fermi weights to 'dft_misc_input'
                if self.misc_subgrp not in archive:
                    archive.create_group(self.misc_subgrp)
                archive[self.misc_subgrp]['dft_fermi_energy'] = fermi_energy
                if symm_kpath_info is not None:
                    archive[self.misc_subgrp].create_group('symm_kpath')
                    kpath_grp = archive[self.misc_subgrp]['symm_kpath']
                    kpath_grp['kpts'] = symm_kpath_info[0]
                    kpath_grp['labels'] = symm_kpath_info[1]
                    kpath_grp['label_idx'] = symm_kpath_info[2]
                if self.bloch_basis:
                    archive[self.misc_subgrp]['dft_fermi_weights'] = f_weights
                    archive[self.misc_subgrp]['band_window'] = band_window+1 # Change to 1-based index
                    archive[self.misc_subgrp]['kpts_cart'] = np.dot(kpts, kpt_basis.T)
        mpi.barrier()

        # Makes Fermi energy a class variable for testing
        self.fermi_energy = fermi_energy


def read_input_file(inp_file, fortran_to_replace):
    """
    Reads the input file.

    We do not read shells, because we have no additional shells beyond correlated
    ones, and the data will be copied from corr_shells into shells (see below).

    Raises NotImplementedError if an unknown kmesh mode is chosen.
    """
    mpi.report('\nReading input from {}'.format(inp_file))

    def filter_and_replace(file, to_replace):
        """
        Filters for comments (#) and empty lines and replaces the to_replace entries.
        Similar to triqs_dft_tools.converters.converter_tools.read_fortran_file.
        """
        for line in file:
            line = line.split('#')[0].strip()
            if not line:
                continue
            for old, new in to_replace.items():
                line = line.replace(old, new)
            yield line.split()

    # Reads in filtered input file
    with open(inp_file, 'r') as file:
        file_content = list(filter_and_replace(file, fortran_to_replace))

    # Defines the properties for the shells and corr_shells lists
    shell_entries = ['atom', 'sort', 'l', 'dim']
    corr_shell_entries = ['atom', 'sort', 'l', 'dim', 'SO', 'irep']

    # Reads k-point mesh generation option
    kmesh_mode = int(file_content[0][0])
    if kmesh_mode not in (-1, 0):
        raise NotImplementedError(f'kmesh_mode {kmesh_mode} not supported')

    if kmesh_mode >= 0:
        # Read k-point mesh size from input
        assert len(file_content[0]) == 4, 'specify k grid dimensions'
        kmesh_size = [int(s) for s in file_content[0][1:]]
    else:
        assert len(file_content[0]) == 1, 'automatic k grid, no specifications needed'
        kmesh_size = None

    # Reads the total number of electrons per cell if not in bloch basis
    # in bloch basis, this is later calculated from the partial occupations
    assert len(file_content[1]) == 1
    if file_content[1][0].lower() == 'none':
        density_required = None
    else:
        density_required = float(file_content[1][0])

    # Reads in number of corr. shells (e.g. Fe d, Ce f) in the unit cell
    assert len(file_content[2]) == 1
    n_corr_shells = int(file_content[2][0])

    # Now reads the information about the correlated shells
    corr_shells = []
    for line in file_content[3:3+n_corr_shells]:
        assert len(line) == len(corr_shell_entries)
        corr_shells.append({name: int(val) for name, val in zip(corr_shell_entries, line)})

    # Reads in the Fermi energy if there is an additional line after the correlated shells
    if len(file_content) > 3+n_corr_shells:
        assert len(file_content[3+n_corr_shells]) == 1
        fermi_energy = float(file_content[3+n_corr_shells][0])
    else:
        fermi_energy = 0.

    # Checks length of file
    assert len(file_content) in (3+n_corr_shells, 4+n_corr_shells), 'input file contains additional lines'

    # Copies corr_shells into shells
    n_shells = n_corr_shells
    shells = [{key: corr_shells[ish][key] for key in shell_entries}
              for ish in range(n_shells)]

    return (kmesh_mode, kmesh_size, density_required, n_corr_shells,
            corr_shells, n_shells, shells, fermi_energy)


def check_and_adapt_for_soc(shells, corr_shells, bloch_basis, add_lambda):
    """
    Checks compatibilities, modifies shells and corr_shells and sets variables
    needed for spin-orbit coupled systems.
    """
    # Determines if any shell requires SO
    if any(sh['SO'] == 1 for sh in corr_shells):
        SO = 1
        SP = 1
        mpi.report('Spin-orbit interaction turned on')
    else:
        SO = 0
        SP = 0

    # Only one block supported - either non-spin-polarized or spin-orbit coupled
    assert SP == SO, 'Spin-polarized calculations not implemented'

    # If adding a local SOC term, turn on SOC
    if add_lambda:
        assert all(sh['dim'] == 3 for sh in corr_shells), 'add_lambda only implemented for t2g shell'
        assert SO == SP == 0, 'add_lambda not implemented for SO = SP = 1'
        assert not bloch_basis, 'add_lambda not implemented for bloch_basis = True'
        # now setting SO and SP to 1
        SO = SP = 1

    # this is more general
    n_spin_blocks = SP + 1 - SO
    assert n_spin_blocks > 0, 'Input error, if SO=1, SP must be 1.'

    # Doubles dimension of all shells if SOC is on
    if SO == 1:
        for shell_list in [shells, corr_shells]:
            for entry in shell_list:
                entry['dim'] *= 2
                if 'SO' in entry.keys() and add_lambda:
                    entry['SO'] = 1

    return shells, corr_shells, SO, SP, n_spin_blocks


def read_wannier90_hr_data(wannier_seed):
    """
    Method for reading the seedname_hr.dat file produced by Wannier90 (http://wannier.org).
    If spin polarized, reads in the properties for a given spin

    Parameters
    ----------
    wannier_seed : string
        seedname to read H(R) file produced by Wannier90, seedname_hr.dat

    Returns
    -------
    n_r_spin : integer
        number of R vectors found in the file
    r_vector_spin : np.ndarray of integers
        Miller indices of the R vectors
    r_degeneracy_spin : np.ndarray of floats
        weight of the R vectors
    n_wannier_spin : integer
        number of Wannier functions found
    wannier_hr_spin : np.ndarray, # R points x # wannier funcs x # wannier funcs
        <w_i|H(R)|w_j> = Hamilonian matrix elements in the Wannier basis
    """
    hr_filename = wannier_seed + '_hr.dat'
    with open(hr_filename, 'r') as file:
        hr_data = file.readlines()

    mpi.report('Reading {}: {}'.format(hr_filename, hr_data[0].strip()))

    n_wannier_spin = int(hr_data[1])
    n_r_spin = int(hr_data[2])

    r_vector_spin = np.zeros((n_r_spin, 3), dtype=int)
    r_degeneracy_spin = np.zeros(n_r_spin, dtype=int)
    wannier_hr_spin = np.zeros((n_r_spin, n_wannier_spin, n_wannier_spin), dtype=complex)

    currpos = 2
    ir = 0
    # read the degeneracy of the R vectors (needed for the Fourier
    # transform)
    while ir < n_r_spin:
        currpos += 1
        for x in hr_data[currpos].split():
            if ir >= n_r_spin:
                raise IndexError('wrong number of R vectors??')
            r_degeneracy_spin[ir] = int(x)
            ir += 1
    # for each direct lattice vector R read the block of the
    # Hamiltonian H(R)
    for ir, jj, ii in product(range(n_r_spin), range(n_wannier_spin), range(n_wannier_spin)):
        # advance one line, split the line into tokens
        currpos += 1
        cline = hr_data[currpos].split()
        # check if the orbital indexes in the file make sense
        assert int(cline[3]) == ii + 1 and int(cline[4]) == jj + 1, 'Inconsistent indices for R vector n. {}'.format(ir)
        rcurr = np.array([int(cline[0]), int(cline[1]), int(cline[2])])
        if ii == 0 and jj == 0:
            r_vector_spin[ir] = rcurr
            rprev = rcurr
        else:
            # check if the vector indices are consistent
            assert np.all(rcurr == rprev), 'Inconsistent indices for R vector n. {}'.format(ir)

        # fill wannier_hr_spin with the matrix elements of the Hamiltonian
        wannier_hr_spin[ir, ii, jj] = complex(float(cline[5]), float(cline[6]))

    return n_r_spin, r_vector_spin, r_degeneracy_spin, n_wannier_spin, wannier_hr_spin


def read_wannier90_blochbasis_data(wannier_seed, n_wannier_spin):
    """
    Method for reading the files needed in the bloch_basis: seedname_u.mat,
    seedname.eig and potentially seedname_u_dis.mat.

    Parameters
    ----------
    wannier_seed : string
        seedname to Wannier90 output

    Returns
    -------
    u_mat_spin : np.ndarray
        U_mn^k = unitary matrix elements which mix the Kohn-Sham states
    udis_mat_spin : np.ndarray
        U^dis(k) = rectangular matrix for entangled bands
    ks_eigenvals_spin : np.ndarray
        \epsilon_nk = Kohn-Sham eigenvalues (in eV) needed for entangled bands
    k_mesh : np.ndarray
        The k mesh read from the seedname_u.mat file to ensure consistency
    """
    mpi.report('Writing h5 archive in projector formalism: H(k) defined in KS Bloch basis')

    # ------------- Reads seedname_u.mat
    # first, read u matrices from 'seedname_u.mat'
    u_filename = wannier_seed + '_u.mat'
    with open(u_filename,'r') as u_file:
        u_data = u_file.readlines()
    # reads number of kpoints and number of wannier functions
    n_k, n_wannier_spin_u, _ = map(int, u_data[1].split())
    assert n_wannier_spin_u == n_wannier_spin, '#WFs must be identical for *_u.mat and *_hr.dat'
    mpi.report('Reading {}: {}'.format(u_filename, u_data[0].strip()))

    # Reads k mesh from all lines with 3 floats
    k_mesh = np.loadtxt((line for line in u_data if line.count('.') == 3), ndmin=2)
    assert k_mesh.shape == (n_k, 3)
    # Reads u matrices from all lines with 2 floats
    u_mat_spin = np.loadtxt((line for line in u_data if line.count('.') == 2))
    assert u_mat_spin.shape == (n_k*n_wannier_spin*n_wannier_spin, 2)

    # first, get the input for u_mat_spin
    u_mat_spin = u_mat_spin[:, 0] + 1j * u_mat_spin[:, 1]
    u_mat_spin = u_mat_spin.reshape((n_k, n_wannier_spin, n_wannier_spin)).transpose((0, 2, 1))

    # ------------- Reading seedname.eig
    band_filename = wannier_seed + '.eig'
    # read Kohn-Sham eigenvalues from 'seedname.eig'
    mpi.report('Reading {}'.format(band_filename))
    band_data = np.loadtxt(band_filename, usecols=2)
    ks_eigenvals_spin = band_data.reshape(n_k, -1)

    # ------------- Reading seedname_u_dis.mat
    udis_filename = wannier_seed + '_u_dis.mat'
    disentangle = os.path.isfile(udis_filename)
    if disentangle:
        with open(udis_filename,'r') as udis_file:
            udis_data = udis_file.readlines()
    else:
        mpi.report('WARNING: File {} missing.'.format(udis_filename))
        mpi.report('Assuming an isolated set of bands. Check if this is what you want!')

    if disentangle:
        # if disentangle the window is needed
        wout_filename = wannier_seed + '.wout'
        # reads number of kpoints, number of wannier functions and bands
        num_k_udis, n_wannier_spin_udis, num_ks_bands = map(int, udis_data[1].split())

        assert num_k_udis == n_k, '#k points must be identical for *.inp and *_u_dis.mat'
        assert n_wannier_spin_udis == n_wannier_spin, '#WFs must be identical for *_u.mat and *_hr.dat'

        mpi.report('Reading {}: {}, '.format(udis_filename, udis_data[0].strip()))
        udis_data = np.loadtxt(udis_data, usecols=(0, 1), skiprows=2)

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
        num_ks_bands = n_wannier_spin

    # now, check what is needed in the case of disentanglement:
    # The file seedname_u_dis.mat contains only the bands inside the window
    # and then fills the rest up with zeros. Therefore, we need to put the
    # entries from udis_data in the correct position in udis_mat_spin, i.e.
    # shifting by the number of bands below dis_window_min
    # reshape band_data
    assert ks_eigenvals_spin.shape[1] == num_ks_bands, '.eig and u_dis.mat data inconsistent'

    if disentangle:
        # In case the disentanglement window is not set by the user, change manually both limits to 
        # larger window to avoid possible counting error in next line
        dis_tol = 1e-5
        shift_dis_down = np.any(np.isclose(np.min(ks_eigenvals_spin, axis=1), dis_window_min, atol=dis_tol, rtol=0.) == True)
        shift_dis_up = np.any(np.isclose(np.max(ks_eigenvals_spin, axis=1), dis_window_max, atol=dis_tol, rtol=0.) == True)
        if shift_dis_down or shift_dis_up:
            if shift_dis_down:
                mpi.report('WARNING: dis_win_min too close to value in .eig, manually shifting it down by '
                           + '{} to avoid counting error'.format(2*dis_tol))
                dis_window_min -= 2*dis_tol
            if shift_dis_up:
                mpi.report('WARNING: dis_win_max too close to value in .eig, manually shifting it up by '
                           + '{} to avoid counting error'.format(2*dis_tol))
                dis_window_max += 2*dis_tol
            mpi.report('This can happen when the user did not specify disentanglement windows in .win. '
                       + 'Check if this is the case here.')

        # Determine which bands are inside the band window
        inside_window = np.logical_and(ks_eigenvals_spin >= dis_window_min,
                                       ks_eigenvals_spin <= dis_window_max)
        n_inside_per_k = np.sum(inside_window, axis=1)

        # Reformats udis_data as complex, without header
        udis_data = udis_data[:, 0] + 1j * udis_data[:, 1]
        udis_data = udis_data.reshape((n_k, n_wannier_spin*num_ks_bands+1))[:, 1:]
        udis_data = udis_data.reshape((n_k, n_wannier_spin, num_ks_bands))

        #initiate U disentanglement matrices and fill from file 'seedname_u_dis.mat'
        udis_mat_spin = np.zeros([n_k, num_ks_bands, n_wannier_spin], dtype=complex)
        for ik in range(n_k):
            # this assumes that u_dis_mat first index is the first state in the dis window
            # which is different from the eig file which starts with the first KS ev given to w90
            udis_mat_spin[ik, inside_window[ik]] = udis_data[ik, :, :n_inside_per_k[ik]].T
            if not np.allclose(udis_data[ik, :, n_inside_per_k[ik]:], 0):
                raise ValueError('This error could come from rounding of the band window in the seedname.wout. '
                                 + 'Never use default outer window but something wider and '
                                 + 'check that your outer window is not close to any band energy.')
    else:
        # no disentanglement; fill udis_mat_spin with identity
        udis_mat_spin = np.array([np.identity(n_wannier_spin, dtype=complex)] * n_k)

    # return the data into variables
    return u_mat_spin, udis_mat_spin, ks_eigenvals_spin, k_mesh


def read_wannier90_centres(wannier_seed):
    centres_filename = wannier_seed + '_centres.xyz'
    if not os.path.isfile(centres_filename):
        mpi.report('Wannier centres file, {}, does not exist. Please set: '
                   'write_xyz = true, translate_home_cell = false'.format(centres_filename))
        return None

    centres = []
    with open(centres_filename, 'r') as c_file:
        c_data = c_file.readlines()
        mpi.report('Reading {}: {}'.format(centres_filename, c_data[1].strip()))
        for l in c_data:
            if l[0] == 'X':
                x = np.asarray(l.split())
                R = x[1:].astype(float)
                centres.append(R)
    centres = np.asarray(centres)

    return centres


def read_wannier90_symm_kpath(wannier_seed):
    kpath_filename = wannier_seed + '_band.kpt'
    label_filename = wannier_seed + '_band.labelinfo.dat'
    if not os.path.isfile(kpath_filename) or not os.path.isfile(label_filename):
        return None

    labels = ''
    label_idx = []
    with open(label_filename, 'r') as label_file:
        label_data = label_file.readlines()
        for line in label_data:
            l, idx = line.split()[:2]
            labels += l
            label_idx.append(int(idx))
    label_idx = np.asarray(label_idx)

    kpts_interpolate = []
    with open(kpath_filename, 'r') as path_file:
        linelist = path_file.readlines()
        mpi.report('Reading {}: high-symmetry path \'{}\' for '
                   'band structure'.format(kpath_filename, labels))
        for line in linelist[1:]:
            x = np.asarray(line.split()[:-1])
            x = x.astype(float)
            kpts_interpolate.append(x)
    kpts_interpolate = np.asarray(kpts_interpolate)

    return kpts_interpolate, labels, label_idx


def read_all_wannier90_data(n_spin_blocks, dim_corr_shells, w90_seed, add_lambda, bloch_basis):
    """
    Reads in all the wannier90 data using the functions read_wannier90_hr_data
    and read_wannier90_blochbasis_data. Reads in everything for each spin
    channel (into the variables marked with _spin) and runs consistency checks
    on them or combines them.

    Returns
    -------
    wannier_hr : np.ndarray[n_spin_blocks, n_r, n_wannier, n_wannier] of complex
        Hamilonian matrix elements in the Wannier basis
    u_total : np.ndarray[n_spin_blocks, n_k, n_wannier, n_bands] of complex
        The projection matrix as read from wannier. None if not bloch_basis
    ks_eigenvals : np.ndarray[n_spin_blocks, n_k, n_bands] of float
        The KS eigenvalues of the bands per k point. None if not bloch_basis
    r_vector : np.ndarray[n_r, 3] of int
        The r vectors of the wannier Hamiltonian
    r_degeneracy : np.ndarray[n_r] of int
        The degeneracy per r point
    n_wannier: int
        Number of wannier functions
    n_bands : int
        Number of bands
    k_mesh_from_umat : np.ndarray[n_k, 3] of float
        The k points as used in wannier for consistency. None if not bloch_basis
    centres: np.ndarray[n_spin_blocks, 3, 3] of float or None
        Centres of wannier functions
    """
    spin_w90name = ['_up', '_down']
    wannier_hr = []
    centres = []
    if bloch_basis:
        u_mat = []
        udis_mat = []
        ks_eigenvals = []

    for isp in range(n_spin_blocks):
        # build filename according to wannier90 conventions
        if n_spin_blocks == 2:
            mpi.report('Reading information for spin component n. {}'.format(isp))
            file_seed = w90_seed + spin_w90name[isp]
        else:
            file_seed = w90_seed
        # now grab the data from the wannier90 output files
        mpi.report('\nThe Hamiltonian in MLWF basis is extracted from {} files'.format(file_seed))

        w90_hr_results = None
        if mpi.is_master_node():
            w90_hr_results = read_wannier90_hr_data(file_seed)
        n_r_spin, r_vector_spin, r_degeneracy_spin, n_wannier_spin, wannier_hr_spin = mpi.bcast(w90_hr_results)

        if bloch_basis:
            w90_results = None
            if mpi.is_master_node():
                w90_results = read_wannier90_blochbasis_data(file_seed, n_wannier_spin)
            # number of R vectors, their indices, their degeneracy, number of WFs, H(R),
            # U matrices, U(dis) matrices, band energies, k_mesh of U matrices
            u_mat_spin, udis_mat_spin, ks_eigenvals_spin, k_mesh_from_umat = mpi.bcast(w90_results)

        w90_centres_results = None
        if mpi.is_master_node():
            w90_centres_results = read_wannier90_centres(file_seed)
        centres_spin = mpi.bcast(w90_centres_results)

        mpi.report('\n... done: {} R vectors, {} WFs found'.format(n_r_spin, n_wannier_spin))

        if add_lambda:
            mpi.report('Adding local spin-orbit term to Hamiltonian (assuming dxz, dyz, dxy as orbital order)')
            # doubling number of wannier functions
            n_wannier_spin *= 2

        if isp == 0:
            # set the R vectors and their degeneracy
            n_r = n_r_spin
            r_vector = r_vector_spin
            r_degeneracy = r_degeneracy_spin
            n_wannier = n_wannier_spin

            # check that the total number of WFs makes sense
            if n_wannier < dim_corr_shells:
                mpi.report('ERROR: number of WFs in the file smaller than number of correlated orbitals!')
            elif n_wannier > dim_corr_shells:
                # NOTE: correlated shells must appear before uncorrelated
                # ones inside the file
                mpi.report('Number of WFs larger than correlated orbitals:',
                           'WFs from {} to {} treated as uncorrelated'.format(dim_corr_shells + 1, n_wannier))
            else:
                mpi.report('Number of WFs equal to number of correlated orbitals')

            # we assume spin up and spin down always have same total number
            # of WFs
            # get second dimension of udis_mat_spin which corresponds to number of bands in window
            if bloch_basis:
                n_bands = udis_mat_spin.shape[1]
                if add_lambda:
                    n_bands *= 2
            else:
                n_bands = n_wannier
        else:
            # consistency check between the _up and _down file contents
            assert n_r_spin == n_r, 'Different number of R vectors for spin-up/spin-down!'
            assert n_wannier_spin == n_wannier, 'Different number of WFs for spin-up/spin-down!'
            assert np.all(r_vector_spin == r_vector), 'R vectors different between spin components'
            assert np.all(r_degeneracy_spin == r_degeneracy), 'R vec. degeneracy different between spin components'

        wannier_hr.append(wannier_hr_spin)
        centres.append(centres_spin)
        if bloch_basis:
            u_mat.append(u_mat_spin)
            udis_mat.append(udis_mat_spin)
            ks_eigenvals.append(ks_eigenvals_spin)

    if bloch_basis:
        # Definition of projectors in Wannier and Triqs different by Hermitian conjugate
        u_total = np.einsum('skab,skbc->skca', udis_mat, u_mat).conj()
        ks_eigenvals = np.array(ks_eigenvals)
    else:
        u_total = None
        ks_eigenvals = None
        k_mesh_from_umat = None
    wannier_hr = np.array(wannier_hr)
    centres = np.array(centres) if centres[0] is not None else None

    return (wannier_hr, u_total, ks_eigenvals, r_vector, r_degeneracy,
            n_wannier, n_bands, k_mesh_from_umat, centres)


def build_kmesh(kmesh_size, kmesh_mode=0):
    """
    Method for the generation of the k-point mesh. Right now it only supports
    the option for generating a full grid containing k=0,0,0.

    Parameters
    ----------
    kmesh_size : list of 3 integers
        the dimensions of the mesh
    kmesh_mode : integer
        mesh generation mode (right now, only full grid available)

    Returns
    -------
    n_k : integer
        total number of k-points in the mesh
    kpts : np.ndarray[n_k, 3] of floats
        the coordinates of all k-points
    kpt_weights : np.ndarray[n_k] of floats
        the weight of each k-point
    """
    if kmesh_mode != 0:
        raise ValueError('Mesh generation mode not supported: {}'.format(kmesh_mode))

    assert len(kmesh_size) == 3

    # a regular mesh including Gamma point
    # total number of k-points
    n_k = np.prod(kmesh_size)
    kpts = np.array(list(product(range(kmesh_size[0]), range(kmesh_size[1]), range(kmesh_size[2]))))
    kpts = kpts / np.array(kmesh_size)
    # weight is equal for all k-points because wannier90 uses uniform grid on whole BZ
    # (normalization is always 1 and takes into account spin degeneracy)
    kpt_weights = np.full(n_k, 1/n_k)

    return n_k, kpts, kpt_weights


def read_misc_input(w90_seed, n_spin_blocks, n_k):
    """
    Reads input from DFT code calculations to get occupations, the band window,
    the Fermi energy and the basis for the k points.

    Parameters
    ----------
    w90_seed : string
        seed of wannier90 calculations
    n_spin_blocks : int
        SP + 1 - SO
    n_k : int
        Number of k points

    Returns
    -------
    fermi_weights : np.ndarray[n_k, n_spin_blocks, n_bands]
        occupations from DFT calculation
    band_window : np.ndarray[n_k, n_spin_blocks, n_bands]
        band indices of correlated subspace
    fermi_energy : float
        the DFT Kohn-Sham Fermi energy
    kpt_basis: np.ndarray[3, 3]
        the basis vectors in reciprocal space
    """
    w90_seed_dir = os.path.dirname(w90_seed)
    nscf_filename = w90_seed + '.nscf.out'
    nnkp_filename = w90_seed + '.nnkp'
    locproj_filename = os.path.join(w90_seed_dir, 'LOCPROJ')
    outcar_filename = os.path.join(w90_seed_dir, 'OUTCAR')

    if os.path.isfile(nscf_filename):
        read_from = 'qe'
        mpi.report('Reading DFT band occupations from Quantum Espresso output {}'.format(nscf_filename))
    elif os.path.isfile(locproj_filename) and os.path.isfile(outcar_filename):
        read_from = 'vasp'
        mpi.report('Reading DFT band occupations from Vasp output {}'.format(locproj_filename))
    else:
        raise IOError('seedname.nscf.out or LOCPROJ and OUTCAR required in bloch_basis mode')

    assert n_spin_blocks == 1, 'spin-polarized not implemented'
    assert read_from in ('qe', 'vasp')

    occupations = []
    reading_kpt_basis = False
    lines_read_kpt_basis = 0
    kpt_basis = np.zeros((3, 3))

    if read_from == 'qe':
        occupations = []
        with open(nscf_filename,'r') as out_file:
            out_data = out_file.readlines()
        # Reads number of Kohn-Sham states and Fermi energy
        for line in out_data:
            if 'number of Kohn-Sham states' in line:
                n_ks = int(line.split()[-1])
            elif 'the Fermi energy is' in line:
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

        assert 'k =' in out_data[ct + 2], 'Cannot read occupations. Set verbosity = "high" in {}'.format(nscf_filename)
        out_data = out_data[ct+2:]

        # block size of eigenvalues + occupations per k-point
        n_block = int(2*np.ceil(n_ks/8)+5)

        for ik in range(n_k):
            # get data
            k_block = [line.split() for line in out_data[ik*n_block+2:ik*n_block+n_block-1]]
            # second half corresponds to occupations
            occs = k_block[int(len(k_block)/2)+1:]
            flattened_occs = [float(item) for sublist in occs for item in sublist]
            occupations.append(flattened_occs)
    else:
        # Reads LOCPROJ
        with open(locproj_filename, 'r') as file:
            header = file.readline()
            n_ks = int(header.split()[2])
            fermi_energy = float(header.split()[4])

            occupations = np.loadtxt((line for line in file if 'orbital' in line), usecols=5)
        occupations = occupations.reshape((n_k, n_ks))

        # Read reciprocal vectors from OUTCAR
        with open(outcar_filename, 'r') as file:
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
    # read exclude_bands from 'seedname.nnkp' file
    with open(nnkp_filename, 'r') as nnkp_file:
        read = False
        skip = False
        for line in nnkp_file:
            if line.strip() == 'begin exclude_bands':
                read = True
                # skip one more line that contains total number of excluded bands
                skip = True
                continue
            if line.strip() == 'end exclude_bands':
                read = False
                continue
            if skip:
                skip = False
                continue
            if read:
                # wannier index -1
                corr_bands.remove(int(line)-1)

    # For now, it is only supported to exclude the lowest and highest bands
    # If bands in the middle are supposed to be excluded, this doesn't work with the band_window
    #     We'd need to manually add rows of zeros to the projectors
    if np.any(np.diff(corr_bands) != 1):
        raise NotImplementedError('Can only exclude the lowest or highest bands')
    band_window = np.array([[(min(corr_bands), max(corr_bands))]*n_k]*n_spin_blocks)

    included_occupations = np.array(occupations)[:, corr_bands]
    # Adds spin dimension without changing the array
    f_weights = included_occupations.reshape(included_occupations.shape[0], 1,
                                             included_occupations.shape[1])

    return f_weights, band_window, fermi_energy, kpt_basis


def reorder_orbital_and_spin(nwfs, wannier_hr, u_total):
    """
    Changes order from VASP5 + wannier90 (first all up spin wannier orbitals,
    then all down spin) to usual order (every up orbital is followed directly
    by the corresponding down orbital). The order between orbital degrees of
    freedom is not changed.

    Returns
    -------
    wannier_hr : np.ndarray
        The re-ordered wannier real-space Hamiltonian.
    u_total : list of np.ndarray
        The re-ordered projection matrix.
    """
    reorder = [i+shift for i in range(nwfs//2) for shift in (0, nwfs//2)]
    mpi.report('Reordering the orbitals because of SOC, new order:', reorder)

    wannier_hr = wannier_hr[:, :, reorder][:, :, :, reorder]
    if u_total is not None:
        u_total = u_total[:, :, reorder, :]

    return wannier_hr, u_total


def generate_local_so_matrix_t2g(add_lambda, n_corr_shells, nwfs):
    """
    Adds local spin-orbit interaction term to the t2g subspace. Orbital order
    is assumed to be the wannier90/triqs order xz_up, xz_dn, yz_up, yz_dn,
    xy_up, xy_dn.
    Parameters are defined as add_lambda = [lambda_x, lambda_y, lambda_z],
    representative of the orbital coupling terms perpendicular to [x, y, z],
    i.e., [d_yz, d_xz, d_xy], respectively.

    Returns
    -------
    lambda_matrix : np.ndarray[6, 6] of complex
        local spin-orbit term to be added to H(0)
    """
    lambda_x, lambda_y, lambda_z = add_lambda

    lambda_matrix = np.zeros((nwfs, nwfs), dtype=complex)

    for icrsh in range(n_corr_shells):
        ind = 6*icrsh
        lambda_matrix[ind  , ind+2] = -1j*lambda_z/2.0
        lambda_matrix[ind  , ind+5] =  1j*lambda_x/2.0
        lambda_matrix[ind+2, ind+5] =    -lambda_y/2.0
        lambda_matrix[ind+4, ind+1] = -1j*lambda_x/2.0
        lambda_matrix[ind+4, ind+3] =     lambda_y/2.0
        lambda_matrix[ind+1, ind+3] =  1j*lambda_z/2.0
    lambda_matrix += lambda_matrix.T.conj()

    return lambda_matrix


def check_hr(wannier_hr, w90_zero, r_zero_index):
    """
    Checks of the real-space Hamiltonian. Prints a warning if Hamiltonian is
    not real and raises a ValueError if Hamiltonian at R=0 is not Hermitian.
    """
    # Tests if Hamiltonian is real
    imag_components_hr = np.isclose(np.imag(wannier_hr), 0, atol=w90_zero, rtol=0)
    index_imag_components = np.nonzero(np.logical_not(np.all(imag_components_hr, axis=(2, 3))))
    for ind in np.transpose(index_imag_components):
        mpi.report('H(R) has large complex components at R {}'.format(ind[1])
                    + ' and spin component {}'.format(ind[0]))

    # Checks if Hamiltonian at R=0 is hermitian
    wannier_hr0 = wannier_hr[:, r_zero_index]
    if not np.allclose(wannier_hr0.transpose((0, 2, 1)).conj(), wannier_hr0, atol=w90_zero, rtol=0):
        raise ValueError('H(R=0) matrix is not Hermitian')


def find_rot_mat(n_corr_shells, corr_shells, shells_map, wannier_hr0,
                 rot_mat_type, w90_zero, n_spin_blocks):
    """
    Method for finding the matrices that bring from local to global coordinate
    systems, based on the eigenvalues of H(R=0).

    Returns
    -------
    rot_mat : list of list of np.ndarray[dim, dim] of complex
        Rotation matrix for each of the shell. None if construction failed
    """
    rot_mat_full = [None] * n_spin_blocks

    for isp in range(n_spin_blocks):
        # initialize the rotation matrices to identities
        rot_mat = [np.identity(corr_shells[icrsh]['dim'], dtype=complex)
                   for icrsh in range(n_corr_shells)]

        # Method none only physical if no equivalent impurities, otherwise
        # potentially useful for debugging
        # Returns identity matrices as rotation matrices
        if rot_mat_type == 'none':
            mpi.report('WARNING: using the method "none" leads to physically wrong results '
                       + 'if there is a mapping of multiple correlated shells on one impurity.')
            rot_mat_full[isp] = rot_mat
            continue

        # TODO: better handling of degenerate eigenvalue case
        hr0_list = [None] * n_corr_shells
        eigval_lst = [None] * n_corr_shells
        eigvec_lst = [None] * n_corr_shells
        ind = 0
        for icrsh in range(n_corr_shells):
            dim = corr_shells[icrsh]['dim']
            # Saves the sub-block of H(0) corresponding to this shell
            hr0_list[icrsh] = np.zeros((dim, dim), dtype=complex)
            hr0_list[icrsh] = wannier_hr0[isp, ind:ind+dim, ind:ind+dim]
            ind += dim
            # Diagonalizes the sub-block for this shell
            eigval_lst[icrsh], eigvec_lst[icrsh] = np.linalg.eigh(hr0_list[icrsh])

        # Checks for degenerate eigenvalues if there are equivalent shells
        # TODO: better handling of degenerate eigenvalue case
        for iineq in set(shells_map):
            if shells_map.count(iineq) > 1:
                icrsh = shells_map.index(iineq)
                dim = corr_shells[icrsh]['dim']
                if any(abs(eigval_lst[icrsh][j] - eigval_lst[icrsh][i]) < w90_zero
                       for i in range(dim) for j in range(i+1, dim)):
                    mpi.report('WARNING: degenerate eigenvalue of H(0) detected for shell {}: '.format(icrsh) +
                               'global-to-local transformation might not work!')

        for icrsh in range(n_corr_shells):
            # build rotation matrices either...
            if rot_mat_type == 'hloc_diag':
                # using the unitary transformations that diagonalize H(0)
                rot_mat[icrsh] = eigvec_lst[icrsh]
            elif rot_mat_type == 'wannier':
                # or by combining those transformations (i.e. for each group,
                # the representative site is chosen as the global frame of reference)
                rot_mat[icrsh] = np.dot(eigvec_lst[icrsh], eigvec_lst[shells_map[icrsh]].T.conj())

            # check that eigenvalues are the same (within accuracy) for
            # equivalent shells
            if not np.allclose(eigval_lst[icrsh], eigval_lst[shells_map[icrsh]], atol=w90_zero, rtol=0):
                mpi.report(f'ERROR: eigenvalue mismatch between equivalent shells! {icrsh}, {shells_map[icrsh]}')
                eigval_diff = eigval_lst[icrsh] - eigval_lst[shells_map[icrsh]]
                mpi.report(f'Eigenvalue difference {eigval_diff}, but threshold set to {w90_zero:.1e}.')
                mpi.report('Consider lowering threshold if you are certain the mapping is correct.')
                return None

            # check that rotation matrices are unitary
            # dim = number of orbitals in this shell
            dim = corr_shells[icrsh]['dim']
            tmp_mat = np.dot(rot_mat[icrsh], rot_mat[icrsh].conj().T)
            if not np.allclose(tmp_mat, np.identity(dim), atol=w90_zero, rtol=0):
                mpi.report(f'ERROR: rot_mat for shell {icrsh:d} is not unitary!')
                return None

            # check that rotation matrices map equivalent H(0) blocks as they should
            # (assuming representative shell as global frame of reference)
            if rot_mat_type == 'hloc_diag':
                tmp_mat = np.dot(rot_mat[icrsh], rot_mat[shells_map[icrsh]].conj().T)
            elif rot_mat_type == 'wannier':
                tmp_mat = rot_mat[icrsh]
            tmp_mat = np.dot(tmp_mat.conj().T, np.dot(hr0_list[icrsh], tmp_mat))
            if not np.allclose(tmp_mat, hr0_list[shells_map[icrsh]], atol=w90_zero, rtol=0):
                mpi.report(f'ERROR: rot_mat does not map H(0) correctly! {icrsh:d}')
                return None

        rot_mat_full[isp] = rot_mat

    # Equality check between rot_mats for different spins
    if n_spin_blocks == 2 and not all(np.allclose(r0, r1, atol=w90_zero, rtol=0)
                                     for r0, r1 in zip(rot_mat_full[0], rot_mat_full[1])):
        mpi.report('Rotations between spin components do not match!')
        return None

    return rot_mat_full[0]


def fourier_transform_hamiltonian(wannier_hr, r_vector, r_degeneracy, kpts):
    """
    Method for obtaining H(k) from H(R) via Fourier transform.

    Parameters
    ----------
    wannier_hr : np.ndarray[n_spin_blocks, n_r, n_wannier, n_wannier] of complex
        Hamiltonian H(R) in Wannier basis
    r_vector : np.ndarray[n_r, 3] of float
        R vectors on which wannier real-space Hamiltonian is defined
    r_degeneracy : np.ndarray[n_r] of int
        Degeneracy of R vector
    kpts : np.ndarray[n_k, 3] of float
        k points where the Fourier transform is executed on

    Returns
    -------
    wannier_hk : np.ndarray[n_spin_blocks, n_k, n_wannier, n_wannier]
        Transformed Hamiltonian H(k) in Wannier basis
    """
    factors = np.exp(2j * np.pi * np.matmul(r_vector, kpts.T)) / r_degeneracy.reshape(-1, 1)
    wannier_hk = np.einsum('rk,srab->skab', factors, wannier_hr)

    return wannier_hk


def check_bloch_basis_hk(n_corr_shells, corr_shells, n_k, n_spin_blocks, n_bands,
                         proj_mat, dim_corr_shells, wannier_hk, hopping):
    """
    Check of the reciprocal-space Hamiltonian in bloch basis. Prints warning if
    the local downfolded Hamiltonian with the projector method does not
    correspond to the W90 result.
    """
    proj_mat_flattened = np.zeros((n_k, n_spin_blocks, dim_corr_shells, n_bands), dtype=complex)
    iorb = 0
    for icrsh in range(n_corr_shells):
        dim = corr_shells[icrsh]['dim']
        proj_mat_flattened[:, :, iorb:iorb+dim, :] = proj_mat[:, :, icrsh, :dim, :]
        iorb += dim

    downfolded_ham = np.einsum('ksab,ksbc,ksdc->skad', proj_mat_flattened, hopping,
                               proj_mat_flattened.conj())
    if dim_corr_shells < n_bands:
        wannier_ham_corr = wannier_hk[:, :, :dim_corr_shells, :dim_corr_shells]
    else:
        wannier_ham_corr = wannier_hk

    hks_are_equal = np.isclose(downfolded_ham, wannier_ham_corr, atol=1e-4, rtol=0)
    if not np.all(hks_are_equal):
        index_difference = np.nonzero(np.logical_not(np.all(hks_are_equal, axis=2)))
        isp = index_difference[0][0]
        ik = index_difference[1][0]
        mpi.report('WARNING: mismatch between downfolded Hamiltonian and Fourier transformed '
                   + 'H(R). First occurred at kpt {} and spin {}:'.format(ik, isp))

        with np.printoptions(formatter={'complexfloat': '{:+.4f}'.format}):
            mpi.report('Downfolded Hamiltonian, P H_eig P')
            mpi.report(downfolded_ham[isp, ik])
            mpi.report('\nWannier Hamiltonian, Fourier(H(r))')
            mpi.report(wannier_ham_corr[isp, ik])


def check_wannier_basis_hk(hopping, dim_corr_shells):
    """
    Check of the reciprocal-space Hamiltonian in wannier basis. Raises error
    if imaginary diagonal elements are not zero, because otherwise there can be
    instabilties in lattice Gf.
    """
    #TODO: do we want to apply this on bloch_basis downfolded Hamiltonian as well?
    diag_iterator = range(dim_corr_shells)
    hk_has_imag_diag = np.isclose(hopping[:, :, diag_iterator, diag_iterator].imag, 0, atol=1e-10, rtol=0)

    if not np.all(hk_has_imag_diag):
        index_imag = np.nonzero(np.logical_not(np.all(hk_has_imag_diag, axis=2)))
        ik, isp = np.transpose(index_imag)[0]
        mpi.report('ERROR: Wannier Hamiltonian has complex diagonal entries. '
                   + 'First occurred at kpt {} and spin {}:'.format(ik, isp))
        with np.printoptions(formatter={'float': '{:+.10f}'.format}):
            mpi.report('\nWannier Hamiltonian diagonal, Fourier(H(r)), imaginary')
            mpi.report(hopping[ik, isp, diag_iterator, diag_iterator].imag)
        raise ValueError
