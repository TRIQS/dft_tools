
import numpy as np

class Projector:
    """
    Class describing a local-orbital projector.
    """
  
    def __init__(self, matrix, ib1=1, ib2=None, nion=1):
        self.p_matrix = matrix.astype(np.complex128)
        self.norb, self.nb = matrix.shape
        self.nion = nion
        self.ib1 = ib1 - 1
        if not ib2 is None:
            self.ib2 = ib2 - 1
        else:
            self.ib2 = self.nb - 1
  
    def project_up(self, mat):
        return np.dot(self.p_matrix.conj().T, np.dot(mat, self.p_matrix))
  
    def project_down(self, mat):
        assert mat.shape == (self.nb, self.nb), "  Matrix must match projector in size"
        return np.dot(self.p_matrix, np.dot(mat, self.p_matrix.conj().T))
  
    def orthogonalize(self):
        """
        Orthogonalizes a projector.
        Returns an overlap matrix and its eigenvalues for initial projectors.
        """
        self.p_matrix, overlap, over_eig = orthogonalize_projector(self.p_matrix)

        return (overlap, over_eig)
        
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
# check_vasp_data_consistency()
################################################################################
def check_vasp_data_consistency(vasp_data):
    """
    Check the consistency of the VASP data.
    """
    pass

################################################################################
# select_bands()
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
    print nk, nband, ns_band
    print emin, emax
    ib_win = np.zeros((nk, ns_band, 2), dtype=np.int32)

    nb_min = 10000000
    nb_max = 0
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

            ib_win[ik, isp, 0] = ib1
            ib_win[ik, isp, 1] = ib2

            nb_min = min(nb_min, ib1)
            nb_max = max(nb_max, ib2)

    return ib_win, nb_min, nb_max

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
        self.emin = gr_pars['emin']
        self.emax = gr_pars['emax']
        self.ishells = gr_pars['shells']
        self.ortho = gr_pars['normalize']
        self.normion = gr_pars['normion']

        self.shells = shells

# Determine the minimum and maximum band numbers
        ib_win, nb_min, nb_max = select_bands(eigvals, self.emin, self.emax)
        self.ib_win = ib_win
        self.nb_min = nb_min
        self.nb_max = nb_max

# Select projectors within the energy window
        for ish in self.ishells:
            shell = self.shells[ish]
            shell.select_projectors(ib_win, nb_min, nb_max)

################################################################################
#
# orthogonalize
#
################################################################################
    def orthogonalize(self):
        """
        Orthogonalize a group of projectors.
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

        ndim = i2
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
# Distribute back projectors in the same order
                for ish in self.ishells:
                    shell = self.shells[ish]
                    blocks = bl_map[ish]['bmat_blocks']
                    for ion in xrange(nion):
                        i1, i2 = blocks[ion]
                        shell.proj_win[ion, isp, ik, :nlm, :nb] = p_mat[i1:i2, :nb]


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
    def __init__(self, sh_pars, proj_raw):
        self.lorb = sh_pars['lshell']
        self.ion_list = sh_pars['ion_list']
        try:
            self.tmatrix = sh_pars['tmatrix']
        except KeyError:
            self.tmatrix = None

        self.lm1 = self.lorb**2
        self.lm2 = (self.lorb+1)**2

# Pre-select a subset of projectors (this should be an array view => no memory is wasted)
        self.proj_arr = proj_raw[self.ion_list, :, :, :, self.lm1:self.lm2]

################################################################################
#
# select_projectors
#
################################################################################
    def select_projectors(self, ib_win, nb_min, nb_max):
        """
        Selects a subset of projectors corresponding to a given energy window.
        """
        self.ib_win = ib_win
        self.nb_min = nb_min
        self.nb_max = nb_max

# Set the dimensions of the array
        nb_win = self.nb_max - self.nb_min + 1
        nion, ns, nk, nbtot, nlm = self.proj_arr.shape
# !!! Note that the order is changed below !!!
        self.proj_win = np.zeros((nion, ns, nk, nb_win, nlm), dtype=np.complex128)

# Select projectors for a given energy window
        ns_band = self.ib_win.shape[1]
        for isp in xrange(ns):
            for ik in xrange(nk):
# TODO: for non-collinear case something else should be done here
                is_b = min(isp, ns_band)
                ib1 = self.ib_win[ik, is_b, 0]
                ib2 = self.ib_win[ik, is_b, 1] + 1
                ib1_win = ib1 - self.nb_min
                ib2_win = ib2 - self.nb_min
                self.proj_win[:, isp, ik, ib1_win:ib2_win, :] = self.proj_arr[:, isp, ik, ib1:ib2, :]

# !!! This sucks but I have to change the order of 'ib' and 'ilm' indices here
# This should perhaps be done right after the projector array is read from PLOCAR
        self.proj_win.transpose((0, 1, 2, 4, 3))


def generate_ortho_plos(conf_pars, vasp_data):
    """
    Parameters
    ----------

      conf_pars (dict) : dictionary of input parameters (from conf-file)
      vasp_data (dict) : dictionary of object representing various VASP files
    """

    check_vasp_data_consistency(vasp_data)

    proj_raw = vaps_data['plocar'].plo
    try:
        efermi = conf_pars.general['efermi']
    except KeyError:
        efermi = vasp_data['doscar'].efermi

# eigvals(nktot, nband, ispin) are defined with respect to the Fermi level
    eigvals = vasp_data['eigenval'].eigs - efermi 

    shells = []
    for sh_par in conf_pars.shells:
        shells.append(ProjectorShell(sh_par, proj_raw))

    groups = []
    for gr_par in conf_pars.groups:
        group = ProjectorGroup(gr_par, shells, eigvals)


