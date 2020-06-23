
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
r"""
    plovasp.proj_shell
    ==================

    Storage and manipulation on projector shells.
"""
def issue_warning(message):
    """
    Issues a warning.
    """
    print()
    print("  !!! WARNING !!!: " + message)
    print()

import itertools as it
import numpy as np
from . import atm

np.set_printoptions(suppress=True)

################################################################################
################################################################################
#
# class ProjectorShell
#
################################################################################
################################################################################
class ProjectorShell:
    """
    Container of projectors related to a specific shell.

    The constructor pre-selects a subset of projectors according to
    the shell parameters passed from the config-file.

    Parameters:

    - sh_pars (dict) : shell parameters from the config-file
    - proj_raw (numpy.array) : array of raw projectors

    """
    def __init__(self, sh_pars, proj_raw, proj_params, kmesh, structure, nc_flag):
        self.lorb = sh_pars['lshell']
        self.ions = sh_pars['ions']
        self.user_index = sh_pars['user_index']
        self.corr = sh_pars['corr']
        self.ion_sort = [sh_pars['ion_sort']]
        self.nc_flag = nc_flag
#        try:
#            self.tmatrix = sh_pars['tmatrix']
#        except KeyError:
#            self.tmatrix = None

        self.lm1 = self.lorb**2
        self.lm2 = (self.lorb+1)**2

        self.nion = self.ions['nion']
# Extract ion list and equivalence classes (ion sorts)
        self.ion_list = sorted(it.chain(*self.ions['ion_list']))

        if self.ion_sort[0] is None:
            self.ion_sort = []
            for ion in self.ion_list:
                for icl, eq_cl in enumerate(self.ions['ion_list']):
                    if ion in eq_cl:
                        self.ion_sort.append(icl + 1) # Enumerate classes starting from 1
                        break

        self.ndim = self.extract_tmatrices(sh_pars)

        self.extract_projectors(proj_raw, proj_params, kmesh, structure)

################################################################################
#
# extract_tmatrices
#
################################################################################
    def extract_tmatrices(self, sh_pars):
        """
        Extracts and interprets transformation matrices provided by the
        config-parser.
        There are two relevant options in 'sh_pars':

          'tmatrix'  : a transformation matrix applied to all ions in the shell
          'tmatrices': interpreted as a set of transformation matrices for each ion.

        If both of the options are present a warning is issued and 'tmatrices'
        supersedes 'tmatrix'.

        Flag 'self.do_transform' is introduced for the optimization purposes
        to avoid superfluous matrix multiplications.
        """
        nion = self.nion
        nm = self.lm2 - self.lm1

        if 'tmatrices' in sh_pars:
            self.do_transform = True

            if 'tmatrix' in sh_pars:
                mess = "Both TRANSFORM and TRANSFILE are specified, TRANSFORM will be ignored."
                issue_warning(mess)

            raw_matrices = sh_pars['tmatrices']
            nrow, ncol = raw_matrices.shape

            assert nrow%nion == 0, "Number of rows in TRANSFILE must be divisible by the number of ions"
            assert ncol%nm == 0, "Number of columns in TRANSFILE must be divisible by the number of orbitals 2*l + 1"

            nr = nrow // nion
            nsize = ncol // nm
            assert nsize in (1, 2, 4), "Number of columns in TRANSFILE must be divisible by either 1, 2, or 4"
#
# Determine the spin-dimension and whether the matrices are real or complex
#
#            if nsize == 1 or nsize == 2:
# Matrices (either real or complex) are spin-independent
#                nls_dim = nm
#                if msize == 2:
#                    is_complex = True
#                else:
#                    is_complex = False
#            elif nsize = 4:
# Matrices are complex and spin-dependent
#                nls_dim = 2 * nm
#                is_complex = True
#
            is_complex = nsize > 1
            ns_dim = max(1, nsize // 2)

# Dimension of the orbital subspace
            assert nr%ns_dim == 0, "Number of rows in TRANSFILE is not compatible with the spin dimension"
            ndim = nr // ns_dim

            self.tmatrices = np.zeros((nion, nr, nm * ns_dim), dtype=np.complex128)

            if is_complex:
                raw_matrices = raw_matrices[:, ::2] + raw_matrices[:, 1::2] * 1j

            for io in range(nion):
                i1 = io * nr
                i2 = (io + 1) * nr
                self.tmatrices[io, :, :] = raw_matrices[i1:i2, :]

            return ndim

        if 'tmatrix' in sh_pars:
            self.do_transform = True

            raw_matrix = sh_pars['tmatrix']
            nrow, ncol = raw_matrix.shape

            assert ncol%nm == 0, "Number of columns in TRANSFORM must be divisible by the number of orbitals 2*l + 1"

# Only spin-independent matrices are expected here
            nsize = ncol // nm
            assert nsize in (1, 2), "Number of columns in TRANSFORM must be divisible by either 1 or 2"

            is_complex = nsize > 1
            if is_complex:
                matrix = raw_matrix[:, ::2] + raw_matrix[:, 1::2] * 1j
            else:
                matrix = raw_matrix

            ndim = nrow

            self.tmatrices = np.zeros((nion, nrow, nm), dtype=np.complex128)
            for io in range(nion):
                self.tmatrices[io, :, :] = raw_matrix

            return ndim

# If no transformation matrices are provided define a default one
        self.do_transform = False

        ns_dim = 2 if self.nc_flag else 1
        ndim = nm * ns_dim

# We still need the matrices for the output
        self.tmatrices = np.zeros((nion, ndim, ndim), dtype=np.complex128)
        for io in range(nion):
            self.tmatrices[io, :, :] = np.identity(ndim, dtype=np.complex128)

        return ndim

################################################################################
#
# extract_projectors
#
################################################################################
    def extract_projectors(self, proj_raw, proj_params, kmesh, structure):
        """
        Extracts projectors for the given shell.

        Projectors are selected from the raw-projector array 'proj_raw'
        according to the shell parameters.
        If necessary the projectors are transformed usin 'self.tmatrices'.
        """
        assert self.nc_flag == False, "Non-collinear case is not implemented"

#        nion = len(self.ion_list)
        nion = self.nion
        nlm = self.lm2 - self.lm1
        _, ns, nk, nb = proj_raw.shape

        if self.do_transform:
# TODO: implement a non-collinear case
#       for a non-collinear case 'ndim' is 'ns * nm'
            ndim = self.tmatrices.shape[1]
            self.proj_arr = np.zeros((nion, ns, nk, ndim, nb), dtype=np.complex128)
            for ik in range(nk):
                kp = kmesh['kpoints'][ik]
                for io, ion in enumerate(self.ion_list):
                    proj_k = np.zeros((ns, nlm, nb), dtype=np.complex128)
                    qcoord = structure['qcoords'][ion]
#                    kphase = np.exp(-2.0j * np.pi * np.dot(kp, qcoord))
#                    kphase = 1.0
                    for m in range(nlm):
# Here we search for the index of the projector with the given isite/l/m indices
                        for ip, par in enumerate(proj_params):
                            if par['isite'] - 1 == ion and par['l'] == self.lorb and par['m'] == m:
                                proj_k[:, m, :] = proj_raw[ip, :, ik, :]  #* kphase
                                break
                    for isp in range(ns):
                        self.proj_arr[io, isp, ik, :, :] = np.dot(self.tmatrices[io, :, :], proj_k[isp, :, :])

        else:
# No transformation: just copy the projectors as they are
            self.proj_arr = np.zeros((nion, ns, nk, nlm, nb), dtype=np.complex128)
            for io, ion in enumerate(self.ion_list):
                qcoord = structure['qcoords'][ion]
                for m in range(nlm):
# Here we search for the index of the projector with the given isite/l/m indices
                    for ip, par in enumerate(proj_params):
                        if par['isite'] - 1 == ion and par['l'] == self.lorb and par['m'] == m:
                            self.proj_arr[io, :, :, m, :] = proj_raw[ip, :, :, :]
#                            for ik in range(nk):
#                                kp = kmesh['kpoints'][ik]
##                                kphase = np.exp(-2.0j * np.pi * np.dot(kp, qcoord))
#                                kphase = 1.0
#                                self.proj_arr[io, :, :, m, :] = proj_raw[ip, :, :, :] # * kphase
                            break


################################################################################
#
# select_projectors
#
################################################################################
    def select_projectors(self, ib_win, ib_min, ib_max):
        """
        Selects a subset of projectors corresponding to a given energy window.
        """
        self.ib_win = ib_win
        self.ib_min = ib_min
        self.ib_max = ib_max
        nb_max = ib_max - ib_min + 1

# Set the dimensions of the array
        nion, ns, nk, nlm, nbtot = self.proj_arr.shape
# !!! Note that the order of the two last indices is different !!!
        self.proj_win = np.zeros((nion, ns, nk, nlm, nb_max), dtype=np.complex128)

# Select projectors for a given energy window
        ns_band = self.ib_win.shape[1]
        for isp in range(ns):
            for ik in range(nk):
# TODO: for non-collinear case something else should be done here
                is_b = min(isp, ns_band)
                ib1 = self.ib_win[ik, is_b, 0]
                ib2 = self.ib_win[ik, is_b, 1] + 1
                ib_win = ib2 - ib1
                self.proj_win[:, isp, ik, :, :ib_win] = self.proj_arr[:, isp, ik, :, ib1:ib2]

################################################################################
#
# density_matrix
#
################################################################################
    def density_matrix(self, el_struct, site_diag=True, spin_diag=True):
        """
        Returns occupation matrix/matrices for the shell.
        """
        nion, ns, nk, nlm, nbtot = self.proj_win.shape

#        assert site_diag, "site_diag = False is not implemented"
        assert spin_diag, "spin_diag = False is not implemented"

        if site_diag:
            occ_mats = np.zeros((ns, nion, nlm, nlm), dtype=np.float64)
            overlaps = np.zeros((ns, nion, nlm, nlm), dtype=np.float64)
        else:
            ndim = nion * nlm
            occ_mats = np.zeros((ns, 1, ndim, ndim), dtype=np.float64)
            overlaps = np.zeros((ns, 1, ndim, ndim), dtype=np.float64)

#        self.proj_win = np.zeros((nion, ns, nk, nlm, nb_max), dtype=np.complex128)
        kweights = el_struct.kmesh['kweights']
        occnums = el_struct.ferw
        ib1 = self.ib_min
        ib2 = self.ib_max + 1
        if site_diag:
            for isp in range(ns):
                for ik, weight, occ in zip(it.count(), kweights, occnums[isp, :, :]):
                    for io in range(nion):
                        proj_k = self.proj_win[io, isp, ik, ...]
                        occ_mats[isp, io, :, :] += np.dot(proj_k * occ[ib1:ib2],
                                                     proj_k.conj().T).real * weight
                        overlaps[isp, io, :, :] += np.dot(proj_k,
                                                     proj_k.conj().T).real * weight
        else:
            proj_k = np.zeros((ndim, nbtot), dtype=np.complex128)
            for isp in range(ns):
                for ik, weight, occ in zip(it.count(), kweights, occnums[isp, :, :]):
                    for io in range(nion):
                        i1 = io * nlm
                        i2 = (io + 1) * nlm
                        proj_k[i1:i2, :] = self.proj_win[io, isp, ik, ...]
                    occ_mats[isp, 0, :, :] += np.dot(proj_k * occ[ib1:ib2],
                                                 proj_k.conj().T).real * weight
                    overlaps[isp, 0, :, :] += np.dot(proj_k,
                                                 proj_k.conj().T).real * weight

#        if not symops is None:
#            occ_mats = symmetrize_matrix_set(occ_mats, symops, ions, perm_map)

        return occ_mats, overlaps

################################################################################
#
# local_hamiltonian
#
################################################################################
    def local_hamiltonian(self, el_struct, site_diag=True, spin_diag=True):
        """
        Returns occupation matrix/matrices for the shell.
        """
        nion, ns, nk, nlm, nbtot = self.proj_win.shape

        assert site_diag, "site_diag = False is not implemented"
        assert spin_diag, "spin_diag = False is not implemented"

        loc_ham = np.zeros((ns, nion, nlm, nlm), dtype=np.complex128)

#        self.proj_win = np.zeros((nion, ns, nk, nlm, nb_max), dtype=np.complex128)
        kweights = el_struct.kmesh['kweights']
        occnums = el_struct.ferw
        ib1 = self.ib_min
        ib2 = self.ib_max + 1
        for isp in range(ns):
            for ik, weight, occ, eigk in zip(it.count(), kweights, occnums[isp, :, :],
                                          el_struct.eigvals[:, ib1:ib2, isp]):
                for io in range(nion):
                    proj_k = self.proj_win[io, isp, ik, ...]
                    loc_ham[isp, io, :, :] += np.dot(proj_k * (eigk - el_struct.efermi),
                                                 proj_k.conj().T) * weight

#        if not symops is None:
#            occ_mats = symmetrize_matrix_set(occ_mats, symops, ions, perm_map)

        return loc_ham

################################################################################
#
# density_of_states
#
################################################################################
    def density_of_states(self, el_struct, emesh):
        """
        Returns projected DOS for the shell.
        """
        nion, ns, nk, nlm, nbtot = self.proj_win.shape

# There is a problem with data storage structure of projectors that will
# make life more complicated. The problem is that band-indices of projectors
# for different k-points do not match because we store 'nb_max' values starting
# from 0.
        nb_max = self.ib_max - self.ib_min + 1
        ns_band = self.ib_win.shape[1]

        ne = len(emesh)
        dos = np.zeros((ne, ns, nion, nlm))
        w_k = np.zeros((nk, nb_max, ns, nion, nlm), dtype=np.complex128)
        for isp in range(ns):
            for ik in range(nk):
                is_b = min(isp, ns_band)
                ib1 = self.ib_win[ik, is_b, 0]
                ib2 = self.ib_win[ik, is_b, 1] + 1
                for ib_g in range(ib1, ib2):
                    for io in range(nion):
# Note the difference between 'ib' and 'ibn':
#  'ib'  counts from 0 to 'nb_k - 1'
#  'ibn' counts from 'ib1 - ib_min' to 'ib2 - ib_min'
                        ib = ib_g - ib1
                        ibn = ib_g - self.ib_min
                        proj_k = self.proj_win[io, isp, ik, :, ib]
                        w_k[ik, ib, isp, io, :] = proj_k * proj_k.conj()

#        eigv_ef = el_struct.eigvals[ik, ib, isp] - el_struct.efermi
        itt = el_struct.kmesh['itet'].T
# k-indices are starting from 0 in Python
        itt[1:, :] -= 1
        for isp in range(ns):
            for ib, eigk in enumerate(el_struct.eigvals[:, self.ib_min:self.ib_max+1, isp].T):
                for ie, e in enumerate(emesh):
                    eigk_ef = eigk - el_struct.efermi
                    cti = atm.dos_tetra_weights_3d(eigk_ef, e, itt)
                    for im in range(nlm):
                        for io in range(nion):
                            dos[ie, isp, io, im] += np.sum((cti * w_k[itt[1:, :], ib, isp, io, im].real).sum(0) * itt[0, :])

        dos *= 2 * el_struct.kmesh['volt']
#        for isp in range(ns):
#            for ik, weight, occ in zip(it.count(), kweights, occnums[isp, :, :]):
#                for io in range(nion):
#                    proj_k = self.proj_win[isp, io, ik, ...]
#                    occ_mats[isp, io, :, :] += np.dot(proj_k * occ[ib1:ib2],
#                                                 proj_k.conj().T).real * weight
#                    overlaps[isp, io, :, :] += np.dot(proj_k,
#                                                 proj_k.conj().T).real * weight

#        if not symops is None:
#            occ_mats = symmetrize_matrix_set(occ_mats, symops, ions, perm_map)

        return dos

################################################################################
################################################################################
#
# class ProjectorShell
#
################################################################################
################################################################################
class ComplementShell(ProjectorShell):
    """
    Container of projectors related to a complement shell.


    Parameters:

    - sh_pars (dict) : shell parameters from the config-file
    - proj_compl (numpy.array) : array of complement projectors

    """
    def __init__(self, sh_pars, proj_compl, nc_flag):
        self.lorb = sh_pars['lshell']
        self.ions = sh_pars['ions']
        self.user_index = sh_pars['user_index']
        self.corr = sh_pars['corr']
        self.nc_flag = nc_flag

        self.ib_min = sh_pars['ib_min']
        self.ib_max = sh_pars['ib_max']
        self.ib_win = sh_pars['ib_win']


        #self.lm1 = self.lorb**2
        #self.lm2 = (self.lorb+1)**2

        self.nion = self.ions['nion']
# Extract ion list and equivalence classes (ion sorts)
        self.ion_list = sorted(it.chain(*self.ions['ion_list']))
        self.ion_sort = []
        for ion in self.ion_list:
            for icl, eq_cl in enumerate(self.ions['ion_list']):
                if ion in eq_cl:
                    self.ion_sort.append(icl + 1) # Enumerate classes starting from 1
                    break

        self.ndim = proj_compl.shape[3]
        self.proj_win = proj_compl

    def extract_tmatrices(self, sh_pars):
        raise Exception('not implemented')

    def local_hamiltonian(self, el_struct, site_diag=True, spin_diag=True):
        raise Exception('not implemented')

    def density_matrix(self, el_struct, site_diag=True, spin_diag=True):
        raise Exception('not implemented')

    #def density_of_states(self, el_struct, emesh):
    #    raise Exception('not implemented')
