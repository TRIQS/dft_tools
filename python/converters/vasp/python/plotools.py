
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
# orthogonalize_projector()
################################################################################
def orthogonalize_projector(p_matrix):
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

    overlap = np.dot(p_matrix, p_matrix.conj().T)
    eig, eigv = np.linalg.eigh(overlap)
    assert np.all(eig > 0.0), ("  Negative eigenvalues of the overlap matrix:"
       "projectors are ill-defined")
    sqrt_eig = np.diag(1.0 / np.sqrt(eig))
    shalf = np.dot(eigv, np.dot(sqrt_eig, eigv.conj().T))
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
    nk, nband, ns_band = eigvals.shape
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
            for ib in xrange(ib1, nb_max):
                en = eigvals[ik, ib, isp]
                if en <= emax:
                    break
            ib2 = ib

            ib_win[ik, isp, 0] = ib1
            ib_win[ik, isp, 1] = ib2

            nb_min = min(nb_min, ib1)
            nb_max = max(nb_max, ib2)

    return ib_win, nb_min, nb_max

################################################################################
#
# class ProjectorGroup
#
################################################################################
class ProjectorGroup:
    """
    Container of projectors defined within a certain energy window.

    The constructor selects a subset of projectors according to
    the parameters from the config-file (passed in `pars`).

    Parameters:

    - pars (dict) : dictionary of parameters from the config-file for a given PLO group
    - proj_raw (numpy.array) : array of raw projectors
    - eigvals (numpy.array) : array of KS eigenvalues

    """
#    def __init__(self, proj_set, nb_min, nb_max, ib_win):
#        """
#        Constructor.
#
#        Parameters
#        ----------
#
#        proj_set (numpy.array) : projector array
#        nb_min (int) : the lowest absolute band index 
#        nb_max (int) : the lowest absolute band index
#        ib_win (numpy.array((nk, ns, 2), dtype=int)) : the lowest and highest band indices
#          for a given `k`-point
#        """
#        self.proj_set = proj_set
#        self.nb_min = nb_min
#        self.nb_max = nb_max
#        self.ib_win = ib_win

#################################################################################
# __init__()
#################################################################################
    def __init__(self, pars, proj_raw, eigvals):
        """
        Constructor
        """
        ns = proj_raw.shape[1]
        nk, nband, ns_band = eigvals.shape

        self.lorb = pars['lshell']
        self.lm_l = range(lorb**2, (lorb+1)**2)
        nlm = len(self.lm_l)

        self.emin = pars['emin']
        self.emax = pars['emax']

# Determine the minimum and maximum band numbers
        ib_win, nb_min, nb_max = select_bands(eigvals, self.emin, self.emax)
        self.ib_win = ib_win
        self.nb_min = nb_min
        self.nb_max = nb_max

# Set the dimensions of the array
        nb_win = self.nb_max - self.nb_min + 1
        nion_sel = pars['ion_list'].shape[0]

        self.proj_set = np.zeros((nion_sel, ns, nk, nb_win, nlm), dtype=np.complex128)
# Select projectors for a given energy window
        for isp in xrange(ns):
            for ik in xrange(nk):
# TODO: for non-collinear case something else should be done here
                is_b = min(isp, ns_band)
                ib1 = self.ib_win[ik, is_b, 0]
                ib2 = self.ib_win[ik, is_b, 1] + 1
                ib1_win = ib1 - self.nb_min
                ib2_win = ib2 - self.nb_min
                for ion, ion_sel in enumerate(pars['ion_list']):
                    self.proj_set[ion, isp, ik, ib1_win:ib2_win, :] = proj_raw[ion_sel, isp, ik, ib1:ib2, self.lm_l]

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

    - pars (dict) : dictionary of parameters from the config-file for a given PLO group
    - proj_raw (numpy.array) : array of raw projectors
    - eigvals (numpy.array) : array of KS eigenvalues

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

# Pre-select a subset of projectors
        self.proj_arr = proj_raw[self.ion_list, :, :, :, lm1:lm2]

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


