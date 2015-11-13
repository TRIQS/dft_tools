
import numpy as np

np.set_printoptions(suppress=True)

################################################################################
#
# orthogonalize_projector_matrix()
#
################################################################################
def orthogonalize_projector_matrix(p_matrix):
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
# Overlap matrix O_{m m'} = \sum_{v} P_{m v} P^{*}_{v m'}
    overlap = np.dot(p_matrix, p_matrix.conj().T)
# Calculate [O^{-1/2}]_{m m'}
    eig, eigv = np.linalg.eigh(overlap)
    assert np.all(eig > 0.0), ("Negative eigenvalues of the overlap matrix:"
       "projectors are ill-defined")
    sqrt_eig = np.diag(1.0 / np.sqrt(eig))
    shalf = np.dot(eigv, np.dot(sqrt_eig, eigv.conj().T))
# Apply \tilde{P}_{m v} = \sum_{m'} [O^{-1/2}]_{m m'} P_{m' v}
    p_ortho = np.dot(shalf, p_matrix)

    return (p_ortho, overlap, eig)

################################################################################
#
# select_bands()
#
################################################################################
def select_bands(eigvals, emin, emax):
    """
    Select a subset of bands lying within a given energy window.
    The band energies are assumed to be sorted in an ascending order.

    Parameters
    ----------
    
    eigvals (numpy.array) : all eigenvalues
    emin, emax (float) : energy window

    Returns
    -------

    ib_win, nb_min, nb_max : 
    """
# Sanity check
    if emin > eigvals.max() or emax < eigvals.min():
        raise Exception("Energy window does not overlap with the band structure")

    nk, nband, ns_band = eigvals.shape
    ib_win = np.zeros((nk, ns_band, 2), dtype=np.int32)

    ib_min = 10000000
    ib_max = 0
    for isp in xrange(ns_band):
        for ik in xrange(nk):
            for ib in xrange(nband):
                en = eigvals[ik, ib, isp]
                if en >= emin:
                    break
            ib1 = ib
            for ib in xrange(ib1, nband):
                en = eigvals[ik, ib, isp]
                if en > emax:
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
    def __init__(self, gr_pars, shells, eigvals, ferw):
        """
        Constructor
        """
        self.emin, self.emax = gr_pars['ewindow']
        self.ishells = gr_pars['shells']
        self.ortho = gr_pars['normalize']
        self.normion = gr_pars['normion']

        self.shells = shells

# Determine the minimum and maximum band numbers
        ib_win, ib_min, ib_max = select_bands(eigvals, self.emin, self.emax)
        self.ib_win = ib_win
        self.ib_min = ib_min
        self.ib_max = ib_max
        self.nb_max = ib_max - ib_min + 1

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
        for isp in xrange(ns_band):
            for ik in xrange(nk):
                ib1 = self.ib_win[ik, isp, 0]
                ib2 = self.ib_win[ik, isp, 1]
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

                 block_maps = [bl_map[ion] for ion in xrange(shell.nion) 
                                           for shell in shells]

                 bl_map = [((i1_start, i1_end), (i1_shell, ion)),
                           ((i2_start, i2_end), (i2_shell, ion)),
                           ...],

          2. Orthogonality is ensured on all sites within the group (NORMION = True).
             The mapping:

                 block_maps = [bl_map]

                 bl_map = [((i1_start, i1_end), (i1_shell, i1_shell.ion1)),
                           ((i1_start, i1_end), (i1_shell, i1_shell.ion2)),
                           ...
                           ((i2_start, i2_end), (i2_shell, i2_shell.ion1)),
                           ((i2_start, i2_end), (i2_shell, i2_shell.ion2)),
                           ...],
             
        """
# Quick exit if no normalization is requested
        if not self.ortho:
            return

# TODO: add the case of 'normion = True'
        assert not self.normion, "'NORMION = True' is not yet implemented"

# Determine the dimension of the projector matrix
# and map the blocks to the big matrix
        i1_bl = 0
        bl_map = [{} for ish in self.ishells]
        for ish in self.ishells:
            _shell = self.shells[ish]
            nion, ns, nk, nlm, nb_max = _shell.proj_win.shape
            bmat_bl = []  # indices corresponding to a big block matrix
            for ion in xrange(nion):
                i2_bl = i1_bl + nlm
                bmat_bl.append((i1_bl, i2_bl))
                i1_bl = i2_bl
            bl_map[ish]['bmat_blocks'] = bmat_bl

        ndim = i2_bl
        p_mat = np.zeros((ndim, nb_max), dtype=np.complex128)
        for isp in xrange(ns):
            for ik in xrange(nk):
                nb = self.ib_win[ik, isp, 1] - self.ib_win[ik, isp, 0] + 1
# Combine all projectors of the group to one block projector
                for ish in self.ishells:
                    shell = self.shells[ish]
                    blocks = bl_map[ish]['bmat_blocks']
                    for ion in xrange(nion):
                        i1, i2 = blocks[ion]
                        p_mat[i1:i2, :nb] = shell.proj_win[ion, isp, ik, :nlm, :nb]
# Now orthogonalize the obtained block projector
                p_orth, overl, eig = orthogonalize_projector_matrix(p_mat)
#                print "ik = ", ik
#                print overl.real
# Distribute back projectors in the same order
                for ish in self.ishells:
                    shell = self.shells[ish]
                    blocks = bl_map[ish]['bmat_blocks']
                    for ion in xrange(nion):
                        i1, i2 = blocks[ion]
                        shell.proj_win[ion, isp, ik, :nlm, :nb] = p_orth[i1:i2, :nb]



