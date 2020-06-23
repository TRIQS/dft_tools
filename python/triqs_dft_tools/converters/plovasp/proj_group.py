
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
    plovasp.proj_group
    ==================

    Storage and manipulation of projector groups.
"""
import numpy as np
from .proj_shell import ComplementShell
np.set_printoptions(suppress=True)

################################################################################
################################################################################
#
# class ProjectorGroup
#
################################################################################
################################################################################
class ProjectorGroup:
    """
    Container of projectors defined within a certain energy window.

    The constructor selects a subset of projectors according to
    the parameters from the config-file (passed in `pars`).

    Parameters:
        - gr_pars (dict) : group parameters from the config-file
        - shells ([ProjectorShell]) : array of ProjectorShell objects
        - eigvals (numpy.array) : array of KS eigenvalues

    """
    def __init__(self, gr_pars, shells, eigvals):
        """
        Constructor
        """
        self.emin, self.emax = gr_pars['ewindow']
        self.ishells = gr_pars['shells']
        self.ortho = gr_pars['normalize']
        self.normion = gr_pars['normion']
        self.complement = gr_pars['complement']

        self.shells = shells

# Determine the minimum and maximum band numbers
        if 'bands' in gr_pars:
            nk, nband, ns_band = eigvals.shape
            ib_win = np.zeros((nk, ns_band, 2), dtype=np.int32)
            ib_win[:,:,0] = gr_pars['bands'][0]-1
            ib_win[:,:,1] = gr_pars['bands'][1]-1
            ib_min = gr_pars['bands'][0] - 1
            ib_max = gr_pars['bands'][1] - 1

        else:
            ib_win, ib_min, ib_max = self.select_bands(eigvals)
        self.ib_win = ib_win
        self.ib_min = ib_min
        self.ib_max = ib_max
        self.nb_max = ib_max - ib_min + 1



        if self.complement:
            n_bands = self.ib_win[:,:,1] - self.ib_win[:,:,0]+1
            n_orbs = sum([x.ndim for x in self.shells])
            assert np.all( n_bands == n_bands[0,0] ), "At each band the same number of bands has to be selected for calculating the complement (to end up with an equal number of orbitals at each k-point)."
            if n_orbs == n_bands[0,0]:
                self.complement = False
                print("\nWARNING: The total number of orbitals in this group is  ")
                print("equal to the number of bands. Setting COMPLEMENT to FALSE!\n")


# Select projectors within the energy window
        for ish in self.ishells:
            shell = self.shells[ish]
            shell.select_projectors(ib_win, ib_min, ib_max)



################################################################################
#
# nelect_window
#
################################################################################
    def nelect_window(self, el_struct):
        """
        Determines the total number of electrons within the window.
        """
        self.nelect = 0
        nk, ns_band, _ = self.ib_win.shape
        rspin = 2.0 if ns_band == 1 else 1.0
        for isp in range(ns_band):
            for ik in range(nk):
                ib1 = self.ib_win[ik, isp, 0]
                ib2 = self.ib_win[ik, isp, 1]+1
                occ = el_struct.ferw[isp, ik, ib1:ib2]
                kwght = el_struct.kmesh['kweights'][ik]
                self.nelect += occ.sum() * kwght * rspin

        return self.nelect

################################################################################
#
# orthogonalize
#
################################################################################
    def orthogonalize(self):
        """
        Orthogonalize a group of projectors.

        There are two options for orthogonalizing projectors:
          1. one ensures orthogonality on each site (NORMION = True);
          2. one ensures orthogonality for subsets of sites (NORMION = False),
             as, e.g., in cluster calculations.

        In order to handle various cases the strategy is first to build a
        mapping that selects appropriate blocks of raw projectors, forms a
        matrix consisting of these blocks, orthogonalize the matrix, and use
        the mapping again to write the orthogonalized projectors back to the
        projector arrays. Note that the blocks can comprise several projector arrays
        contained in different projector shells.

        The construction of block maps is performed in 'self.get_block_matrix_map()'.
        """
# Quick exit if no normalization is requested
        if not self.ortho:
            return

        block_maps, ndim = self.get_block_matrix_map()

        _, ns, nk, _, _ = self.shells[0].proj_win.shape
        p_mat = np.zeros((ndim, self.nb_max), dtype=np.complex128)
# Note that 'ns' and 'nk' are the same for all shells
        for isp in range(ns):
            for ik in range(nk):
                nb = self.ib_win[ik, isp, 1] - self.ib_win[ik, isp, 0] + 1
# Combine all projectors of the group to one block projector
                for bl_map in block_maps:
                    p_mat[:, :] = 0.0j  # !!! Clean-up from the last k-point and block!
                    for ibl, block in enumerate(bl_map):
                        i1, i2 = block['bmat_range']
                        ish, ion = block['shell_ion']
                        nlm = i2 - i1 + 1
                        shell = self.shells[ish]
                        p_mat[i1:i2, :nb] = shell.proj_win[ion, isp, ik, :nlm, :nb]
# Now orthogonalize the obtained block projector
                    ibl_max = i2
                    p_orth, overl, eig = self.orthogonalize_projector_matrix(p_mat[:ibl_max, :nb])
# Distribute projectors back using the same mapping
                    for ibl, block in enumerate(bl_map):
                        i1, i2 = block['bmat_range']
                        ish, ion = block['shell_ion']
                        nlm = i2 - i1 + 1
                        shell = self.shells[ish]
                        shell.proj_win[ion, isp, ik, :nlm, :nb] = p_orth[i1:i2, :nb]


################################################################################
#
# calc_hk
#
################################################################################
    def calc_hk(self, eigvals):
        """
        Calculate H(k) for a group by applying the projectors P
        to the eigenvalues eps.

        H_ij(k) = sum_l P*_il eps_l P_lj

        """

# here we abuse the get_block_matrix_map(), however, it only works
# if self.normion is false
        temp = self.normion
        self.normion = False
        block_maps, ndim = self.get_block_matrix_map()
        self.normion = temp

        _, ns, nk, _, _ = self.shells[0].proj_win.shape

        self.hk = np.zeros((ns,nk,ndim,ndim), dtype=np.complex128)
# Note that 'ns' and 'nk' are the same for all shells
        for isp in range(ns):
            for ik in range(nk):
                bmin = self.ib_win[ik, isp, 0]
                bmax = self.ib_win[ik, isp, 1]+1

                nb = bmax - bmin
                p_mat = np.zeros((ndim, nb), dtype=np.complex128)
                #print(bmin,bmax,nb)
# Combine all projectors of the group to one block projector
                for bl_map in block_maps:
                    p_mat[:, :] = 0.0j  # !!! Clean-up from the last k-point and block!
                    for ibl, block in enumerate(bl_map):
                        i1, i2 = block['bmat_range']
                        ish, ion = block['shell_ion']
                        nlm = i2 - i1 + 1
                        shell = self.shells[ish]
                        p_mat[i1:i2, :nb] = shell.proj_win[ion, isp, ik, :nlm, :nb]

                self.hk[isp,ik,:,:] = np.dot(p_mat*eigvals[ik,bmin:bmax,isp],
                                        p_mat.transpose().conjugate())


################################################################################
#
# complement
#
################################################################################
    def calc_complement(self,eigvals):
        """
        Calculate the complement for a group of projectors.

        This leads to quadtratic projectors P = <l|n> by using a Gram-Schmidt.

        The projector on the orthogonal complement of the existing projectors
        |l> is P^u = 1 - sum_l |l><l|
        We get candidates for complement projectors by applying P^u to a Bloch
        state |n>: |l*> = P^u |n>. For numerical stability we select that Bloch
        state which leads to the |l*> with the largest norm (that corresponds to
        that Bloch state with the smallest overlap with the space spanned by |l>)
        We normalize |l*> and add it to |l>. We do so untill we have as many
        |l> states as we have |n> states.

        """

        print('\nCalculating complement\n')

        block_maps, ndim = self.get_block_matrix_map()
        _, ns, nk, _, _ = self.shells[0].proj_win.shape
        p_mat = np.zeros((ndim, self.nb_max), dtype=np.complex128)
        p_full = np.zeros((1,ns,nk,self.nb_max, self.nb_max), dtype=np.complex128)

# Note that 'ns' and 'nk' are the same for all shells


        for isp in range(ns):
            for ik in range(nk):
                bmin = self.ib_win[ik, isp, 0]
                bmax = self.ib_win[ik, isp, 1]+1

                nb = bmax - bmin
# Combine all projectors of the group to one block projector
                for bl_map in block_maps:
                    p_mat[:, :] = 0.0j  # !!! Clean-up from the last k-point and block!
                    for ibl, block in enumerate(bl_map):
                        i1, i2 = block['bmat_range']
                        ish, ion = block['shell_ion']
                        nlm = i2 - i1 + 1
                        shell = self.shells[ish]
                        p_mat[i1:i2, :nb] = shell.proj_win[ion, isp, ik, :nlm, :nb]
                orbs_done = 1*ndim
                p_full[0,isp,ik,:ndim,:] = p_mat
                while orbs_done < self.nb_max:
#We calculate the overlap of all bloch states: sum_l <n|l><l|m>
                    overlap = np.dot(p_full[0,isp,ik,:orbs_done,:].transpose().conjugate(),p_full[0,isp,ik,:orbs_done,:])
# work is the projector onto the orthogonal complment <n| ( 1 - sum_l |l><l| ) |m>
                    work = np.eye(self.nb_max) - overlap
# calculate the norm of the projected bloch function
                    norm = np.sqrt(np.sum(work*work.transpose(),axis=1))
# select the bloch function leading to the largest norm
                    max_ind = np.argmax(norm)
# normalize and put it to the projectors
                    p_full[0,isp,ik,orbs_done,:] = work[:,max_ind].conjugate()/norm[max_ind]

                    orbs_done += 1

        sh_pars = {}
        sh_pars['lshell'] = -1
        sh_pars['ions'] = {'nion':1,'ion_list':[[1]]}
        sh_pars['user_index'] = 'complement'
        sh_pars['corr']  = False
        sh_pars['ib_min']  = bmin
        sh_pars['ib_max']  = bmax
        sh_pars['ib_win']  = self.ib_win

        self.shells.append(ComplementShell(sh_pars,p_full[:,:,:,ndim:,:],False))
        self.ishells.append(self.ishells[-1]+1)


################################################################################
#
# gen_block_matrix_map
#
################################################################################
    def get_block_matrix_map(self):
        """
        Generates a map from a set of projectors belonging to different shells
        and ions onto a set of block projector matrices, each of which is
        orthonormalized.

        Returns the map and the maximum orbital dimension of the block projector
        matrix.


        Mapping is defined as a list of 'block_maps' corresponding to subsets
        of projectors to be orthogonalized. Each subset corresponds to a subset of sites
        and spans all orbital indices. defined by 'bl_map' as

           bl_map = [((i1_start, i1_end), (i1_shell, i1_ion)),
                     ((i2_start, i2_end), (i2_shell, i2_ion)),
                     ...],

        where `iX_start`, `iX_end` is the range of indices of the block matrix
        (in Python convention `iX_end = iX_last + 1`, with `iX_last` being the last index
        of the range),
        `iX_shell` and `iX_ion` the shell and site indices. The length of the range
        should be consistent with 'nlm' dimensions of a corresponding shell, i.e.,
        `iX_end - iX_start = nlm[iX_shell]`.

        Consider particular cases:
          1. Orthogonality is ensured on each site (NORMION = True).
             For each site 'ion' we have the following mapping:

                 block_maps = [bl_map[ion] for ion in range(shell.nion)
                                           for shell in shells]

                 bl_map = [((i1_start, i1_end), (i1_shell, ion)),
                           ((i2_start, i2_end), (i2_shell, ion)),
                           ...],

          2. Orthogonality is ensured on all sites within the group (NORMION = False).
             The mapping:

                 block_maps = [bl_map]

                 bl_map = [((i1_start, i1_end), (i1_shell, i1_shell.ion1)),
                           ((i1_start, i1_end), (i1_shell, i1_shell.ion2)),
                           ...
                           ((i2_start, i2_end), (i2_shell, i2_shell.ion1)),
                           ((i2_start, i2_end), (i2_shell, i2_shell.ion2)),
                           ...],
        """
        if self.normion:
# Projectors for each site are mapped onto a separate block matrix
            block_maps = []
            ndim = 0
            for ish in self.ishells:
                _shell = self.shells[ish]
                nion, ns, nk, nlm, nb_max = _shell.proj_win.shape
                ndim = max(ndim, nlm)
                for ion in range(nion):
                    i1_bl = 0
                    i2_bl = nlm
                    block = {'bmat_range': (i1_bl, i2_bl)}
                    block['shell_ion'] = (ish, ion)
                    bl_map = [block]
                    block_maps.append(bl_map)

        else:
# All projectors within a group are combined into one big block matrix
            block_maps = []
            bl_map = []
            i1_bl = 0
            for ish in self.ishells:
                _shell = self.shells[ish]
                nion, ns, nk, nlm, nb_max = _shell.proj_win.shape
                for ion in range(nion):
                    i2_bl = i1_bl + nlm
                    block = {'bmat_range': (i1_bl, i2_bl)}
                    block['shell_ion'] = (ish, ion)
                    bl_map.append(block)
                    i1_bl = i2_bl

            ndim = i2_bl
            block_maps.append(bl_map)

        return block_maps, ndim

################################################################################
#
# orthogonalize_projector_matrix()
#
################################################################################
    def orthogonalize_projector_matrix(self, p_matrix):
        """
        Orthogonalizes a projector defined by a rectangular matrix `p_matrix`.

        Parameters
        ----------

        p_matrix (numpy.array[complex]) : matrix `Nm x Nb`, where `Nm` is
          the number of orbitals, `Nb` number of bands

        Returns
        -------

        Orthogonalized projector matrix, initial overlap matrix and its eigenvalues.
        """
# TODO: check the precision of the calculations below,
#       it seems to be inferior to that of Fortran implementation
# Overlap matrix O_{m m'} = \sum_{v} P_{m v} P^{*}_{v m'}
        overlap = np.dot(p_matrix, p_matrix.conj().T)
# Calculate [O^{-1/2}]_{m m'}
        eig, eigv = np.linalg.eigh(overlap)
        assert np.all(eig > 0.0), ("Negative eigenvalues of the overlap matrix:"
           "projectors are ill-defined")
        sqrt_eig = 1.0 / np.sqrt(eig)
        shalf = np.dot(eigv * sqrt_eig, eigv.conj().T)
# Apply \tilde{P}_{m v} = \sum_{m'} [O^{-1/2}]_{m m'} P_{m' v}
        p_ortho = np.dot(shalf, p_matrix)

        return (p_ortho, overlap, eig)

################################################################################
#
# select_bands()
#
################################################################################
    def select_bands(self, eigvals):
        """
        Select a subset of bands lying within a given energy window.
        The band energies are assumed to be sorted in an ascending order.

        Parameters
        ----------

        eigvals (numpy.array) : all eigenvalues
        emin, emax (float) : energy window

        Returns
        -------

        ib_win, nb_min, nb_max : lowest and highest indices of the selected bands

        """
# Sanity check
        if self.emin > eigvals.max() or self.emax < eigvals.min():
            raise Exception("Energy window does not overlap with the band structure")

        nk, nband, ns_band = eigvals.shape
        ib_win = np.zeros((nk, ns_band, 2), dtype=np.int32)

        ib_min = 10000000
        ib_max = 0
        for isp in range(ns_band):
            for ik in range(nk):
                for ib in range(nband):
                    en = eigvals[ik, ib, isp]
                    if en >= self.emin:
                        break
                ib1 = ib
                for ib in range(ib1, nband):
                    en = eigvals[ik, ib, isp]
                    if en > self.emax:
                        break
                else:
# If we reached the last band add 1 to get the correct bound
                    ib += 1
                ib2 = ib - 1

                assert ib1 <= ib2, "No bands inside the window for ik = %s"%(ik)

                ib_win[ik, isp, 0] = ib1
                ib_win[ik, isp, 1] = ib2

                ib_min = min(ib_min, ib1)
                ib_max = max(ib_max, ib2)

        return ib_win, ib_min, ib_max
