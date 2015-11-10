
import itertools as it
import numpy as np
import vasp.atm.c_atm_dos as c_atm_dos

np.set_printoptions(suppress=True)

# 'simplejson' is supposed to be faster than 'json' in stdlib.
try:
    import simplejson as json
except ImportError:
    import json

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
# check_data_consistency()
################################################################################
def check_data_consistency(pars, el_struct):
    """
    Check the consistency of the VASP data.
    """
# Check that ions inside each shell are of the same sort
    for sh in pars.shells:
        assert max(sh['ion_list']) <= el_struct.natom, "Site index in the projected shell exceeds the number of ions in the structure"
        sorts = set([el_struct.type_of_ion[io] for io in sh['ion_list']])
        assert len(sorts) == 1, "Each projected shell must contain only ions of the same sort"

# Check that ion and orbital lists in shells match those of projectors
        ion_list = sh['ion_list']
        lshell = sh['lshell']
        for ion in ion_list:
            for par in el_struct.proj_params:
                if par['isite'] - 1 == ion and par['l'] == lshell:
                    break
            else:
                errmsg = "Projector for isite = %s, l = %s does not match PROJCAR"%(ion + 1, lshell)
                raise Exception(errmsg)
        

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
        self.emin = gr_pars['emin']
        self.emax = gr_pars['emax']
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
    def __init__(self, sh_pars, proj_raw, proj_params):
        self.lorb = sh_pars['lshell']
        self.ion_list = sh_pars['ion_list']
        self.user_index = sh_pars['user_index']
        try:
            self.tmatrix = sh_pars['tmatrix']
        except KeyError:
            self.tmatrix = None

        self.lm1 = self.lorb**2
        self.lm2 = (self.lorb+1)**2

        if self.tmatrix is None:
            self.ndim = self.lm2 - self.lm1
        else:
# TODO: generalize this to a tmatrix for every ion
            self.ndim = self.tmatrix.shape[0]

# Pre-select a subset of projectors (this should be an array view => no memory is wasted)
# !!! This sucks but I have to change the order of 'ib' and 'ilm' indices here
# This should perhaps be done right after the projector array is read from PLOCAR
#        self.proj_arr = proj_raw[self.ion_list, :, :, :, self.lm1:self.lm2].transpose((0, 1, 2, 4, 3))
# We want to select projectors from 'proj_raw' and form an array
#   self.proj_arr[nion, ns, nk, nlm, nb]
# TODO: think of a smart way of copying the selected projectors
#       perhaps, by redesigning slightly the structure of 'proj_arr' and 'proj_win'
#       or by storing only a mapping between site/orbitals and indices of 'proj_raw'
#        iproj_l = []
        nion = len(self.ion_list)
        nlm = self.lm2 - self.lm1
        _, ns, nk, nb = proj_raw.shape
        self.proj_arr = np.zeros((nion, ns, nk, nlm, nb), dtype=np.complex128)
        for io, ion in enumerate(self.ion_list):
            for m in xrange(nlm):
# Here we search for the index of the projector with the given isite/l/m indices
                for ip, par in enumerate(proj_params):
                    if par['isite'] - 1 == ion and par['l'] == self.lorb and par['m'] == m:
#                        iproj_l.append(ip)
                        self.proj_arr[io, :, :, m, :] = proj_raw[ip, :, :, :]
                        break

#        self.proj_arr = proj_raw[iproj_l, :, :, :].transpose((1, 2, 0, 3))

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
        for isp in xrange(ns):
            for ik in xrange(nk):
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

        assert site_diag, "site_diag = False is not implemented"
        assert spin_diag, "spin_diag = False is not implemented"

        occ_mats = np.zeros((ns, nion, nlm, nlm), dtype=np.float64)
        overlaps = np.zeros((ns, nion, nlm, nlm), dtype=np.float64)

#        self.proj_win = np.zeros((nion, ns, nk, nlm, nb_max), dtype=np.complex128)
        kweights = el_struct.kmesh['kweights']
        occnums = el_struct.ferw
        ib1 = self.ib_min
        ib2 = self.ib_max + 1
        for isp in xrange(ns):
            for ik, weight, occ in it.izip(it.count(), kweights, occnums[isp, :, :]):
                for io in xrange(nion):
                    proj_k = self.proj_win[io, isp, ik, ...]
                    occ_mats[isp, io, :, :] += np.dot(proj_k * occ[ib1:ib2],
                                                 proj_k.conj().T).real * weight
                    overlaps[isp, io, :, :] += np.dot(proj_k,
                                                 proj_k.conj().T).real * weight

#        if not symops is None:
#            occ_mats = symmetrize_matrix_set(occ_mats, symops, ions, perm_map)

        return occ_mats, overlaps

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
        for isp in xrange(ns):
            for ik in xrange(nk):
                is_b = min(isp, ns_band)
                ib1 = self.ib_win[ik, is_b, 0]
                ib2 = self.ib_win[ik, is_b, 1] + 1
                for ib_g in xrange(ib1, ib2):
                    for io in xrange(nion):
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
        for isp in xrange(ns):
            for ib, eigk in enumerate(el_struct.eigvals[:, self.ib_min:self.ib_max+1, isp].T):
                for ie, e in enumerate(emesh):
                    eigk_ef = eigk - el_struct.efermi
                    cti = c_atm_dos.dos_weights_3d(eigk_ef, e, itt)
                    for im in xrange(nlm):
                        for io in xrange(nion):
                            dos[ie, isp, io, im] += np.sum((cti * w_k[itt[1:, :], ib, isp, io, im].real).sum(0) * itt[0, :])

        dos *= 2 * el_struct.kmesh['volt']
#        for isp in xrange(ns):
#            for ik, weight, occ in it.izip(it.count(), kweights, occnums[isp, :, :]):
#                for io in xrange(nion):
#                    proj_k = self.proj_win[isp, io, ik, ...]
#                    occ_mats[isp, io, :, :] += np.dot(proj_k * occ[ib1:ib2],
#                                                 proj_k.conj().T).real * weight
#                    overlaps[isp, io, :, :] += np.dot(proj_k,
#                                                 proj_k.conj().T).real * weight

#        if not symops is None:
#            occ_mats = symmetrize_matrix_set(occ_mats, symops, ions, perm_map)

        return dos



################################################################################
#
# generate_plo()
#
################################################################################
def generate_plo(conf_pars, el_struct):
    """
    Parameters
    ----------

      conf_pars (dict) : dictionary of input parameters (from conf-file)
      el_struct : ElectronicStructure object
    """

    check_data_consistency(conf_pars, el_struct)

    proj_raw = el_struct.proj_raw
    try:
        efermi = conf_pars.general['efermi']
    except (KeyError, AttributeError):
        efermi = el_struct.efermi

# eigvals(nktot, nband, ispin) are defined with respect to the Fermi level
    eigvals = el_struct.eigvals - efermi

    nshell = len(conf_pars.shells)
    print
    print "  Generating %i shell%s..."%(nshell, '' if nshell == 1 else 's')
    pshells = []
    for sh_par in conf_pars.shells:
        pshell = ProjectorShell(sh_par, proj_raw, el_struct.proj_params)
        print
        print "    Shell         : %s"%(pshell.user_index)
        print "    Orbital l     : %i"%(pshell.lorb)
        print "    Number of ions: %i"%(len(pshell.ion_list))
        pshells.append(pshell)

    pgroups = []
    for gr_par in conf_pars.groups:
        pgroup = ProjectorGroup(gr_par, pshells, eigvals, el_struct.ferw)
        pgroup.orthogonalize()
        print "Density matrix:"
        dm, ov = pshells[pgroup.ishells[0]].density_matrix(el_struct)
        print dm
        print
        print "Overlap:"
        print ov
        if 'dosmesh' in conf_pars.general:
            print
            print "Evaluating DOS..."
            mesh_pars = conf_pars.general['dosmesh']
            if np.isnan(mesh_pars['emin']):
                dos_emin = pgroup.emin
                dos_emax = pgroup.emax
            else:
                dos_emin = mesh_pars['emin']
                dos_emax = mesh_pars['emax']
            n_points = mesh_pars['n_points']

            emesh = np.linspace(dos_emin, dos_emax, n_points)
            dos = pshells[pgroup.ishells[0]].density_of_states(el_struct, emesh)
            de = emesh[1] - emesh[0]
            ntot = (dos[1:,...] + dos[:-1,...]).sum(0) / 2 * de
            print "  Total number of states:", ntot
            for io in xrange(dos.shape[2]):
                np.savetxt('pdos_%i.dat'%(io), np.vstack((emesh.T, dos[:, 0, io, :].T)).T)

        pgroups.append(pgroup)

    return pshells, pgroups


# TODO: k-points with weights should be stored once and for all
################################################################################
#
# kpoints_output
#
################################################################################
def kpoints_output(basename, el_struct):
    """
    Outputs k-point data into a text file.
    """

    kmesh = el_struct.kmesh
    fname = basename + '.kpoints'
    with open(fname, 'wt') as f:
        f.write("# Number of k-points: nktot\n")
        nktot = kmesh['nktot']
        f.write("%i\n"%(nktot))
# TODO: add the output of reciprocal lattice vectors
        f.write("# List of k-points with weights\n")
        for ik in xrange(nktot):
            kx, ky, kz = kmesh['kpoints'][ik, :]
            kwght = kmesh['kweights'][ik]
            f.write("%15.10f%15.10f%15.10f%20.10f\n"%(kx, ky, kz, kwght))

# Check if there are tetrahedra defined and if they are, output them
        try:
            ntet = kmesh['ntet']
            volt = kmesh['volt']
            f.write("\n# Number of tetrahedra and volume: ntet, volt\n")
            f.write("%i %s\n"%(ntet, volt))
            f.write("# List of tetrahedra: imult, ik1, ..., ik4\n")
            for it in xrange(ntet):
                f.write('  '.join(map("{0:d}".format, *kmesh['itet'][it, :])) + '\n')
        except KeyError:
            pass


################################################################################
#
# ctrl_output
#
################################################################################
def ctrl_output(conf_pars, el_struct, ng):
    """
    Outputs a ctrl-file.
    """
    ctrl_fname = conf_pars.general['basename'] + '.ctrl'
    head_dict = {}

# TODO: Add output of tetrahedra
# Construct the header dictionary
    head_dict['ngroups'] = ng
    head_dict['nk'] = el_struct.kmesh['nktot']
    head_dict['ns'] = el_struct.nspin
    head_dict['nc_flag'] = 1 if el_struct.nc_flag else 0
#    head_dict['efermi'] = conf_pars.general['efermi']  # We probably don't need Efermi

    header = json.dumps(head_dict, indent=4, separators=(',', ': '))

    print "  Storing ctrl-file..."
    with open(ctrl_fname, 'wt') as f:
        f.write(header + "\n")
        f.write("#END OF HEADER\n")

        f.write("# k-points and weights\n")
        labels = ['kx', 'ky', 'kz', 'kweight']
        out = "".join(map(lambda s: s.center(15), labels))
        f.write("#" + out + "\n")
        for ik, kp in enumerate(el_struct.kmesh['kpoints']):
            tmp1 = "".join(map("{0:15.10f}".format, kp))
            out = tmp1 + "{0:16.10f}".format(el_struct.kmesh['kweights'][ik])
            f.write(out + "\n")


################################################################################
#
# plo_output
#
################################################################################
def plo_output(conf_pars, el_struct, pshells, pgroups):
    """
    Outputs PLO groups into text files.

    Filenames are defined by <basename> that is passed from config-file.
    All necessary general parameters are stored in a file '<basename>.ctrl'.

    Each group is stored in a '<basename>.plog<Ng>' file. The format is the
    following:

    # Energy window: emin, emax
    ib_min, ib_max
    nelect
    # Eigenvalues
    isp, ik1, kx, ky, kz, kweight
    ib1, ib2
    eig1
    eig2
    ...
    eigN
    ik2, kx, ky, kz, kweight
    ...

    # Projected shells
    Nshells
    # Shells: <shell indices>
    # Shell <1>
    Shell 1
    ndim
    # complex arrays: plo(ns, nion, ndim, nb)
    ...
    # Shells: <shell indices>
    # Shell <2>
    Shell 2
    ...

    """
    for ig, pgroup in enumerate(pgroups):
        plo_fname = conf_pars.general['basename'] + '.pg%i'%(ig + 1)
        print "  Storing PLO-group file '%s'..."%(plo_fname)
        head_dict = {}

        head_dict['ewindow'] = (pgroup.emin, pgroup.emax)
        head_dict['nb_max'] = pgroup.nb_max

# Number of electrons within the window
        head_dict['nelect'] = pgroup.nelect_window(el_struct)
        print "  Density within window:", head_dict['nelect']

        head_shells = []
        for ish in pgroup.ishells:
            shell = pgroup.shells[ish]
            sh_dict = {}
            sh_dict['shell_index'] = ish
            sh_dict['lorb'] = shell.lorb
            sh_dict['ndim'] = shell.ndim
# Convert ion indices from the internal representation (starting from 0)
# to conventional VASP representation (starting from 1)
            ion_output = [io + 1 for io in shell.ion_list]
            sh_dict['ion_list'] = ion_output
            sh_dict['ion_sort'] = el_struct.type_of_ion[shell.ion_list[0]]

# TODO: add the output of transformation matrices

            head_shells.append(sh_dict)

        head_dict['shells'] = head_shells

        header = json.dumps(head_dict, indent=4, separators=(',', ': '))

        with open(plo_fname, 'wt') as f:
            f.write(header + "\n")
            f.write("#END OF HEADER\n")
            
# Eigenvalues within the window
            f.write("# Eigenvalues within the energy window: %s, %s\n"%(pgroup.emin, pgroup.emax))
            nk, nband, ns_band = el_struct.eigvals.shape
            for isp in xrange(ns_band):
                f.write("# is = %i\n"%(isp + 1))
                for ik in xrange(nk):
                    ib1, ib2 = pgroup.ib_win[ik, isp, 0], pgroup.ib_win[ik, isp, 1]
                    f.write(" %i  %i\n"%(ib1, ib2))
                    for ib in xrange(ib1, ib2 + 1):
                        eigv_ef = el_struct.eigvals[ik, ib, isp] - el_struct.efermi
                        f.write("%15.7f\n"%(eigv_ef))

# Projected shells
            f.write("# Projected shells\n")
            f.write("# Shells: %s\n"%(pgroup.ishells))
            for ish in pgroup.ishells:
                shell = pgroup.shells[ish]
                f.write("# Shell %i\n"%(ish))

                nion, ns, nk, nlm, nb = shell.proj_win.shape
                for isp in xrange(ns):
                    f.write("# is = %i\n"%(isp + 1))
                    for ik in xrange(nk):
                        f.write("# ik = %i\n"%(ik + 1))
                        for ion in xrange(nion):
                            for ilm in xrange(nlm):
                                ib1, ib2 = pgroup.ib_win[ik, isp, 0], pgroup.ib_win[ik, isp, 1]
                                ib_win = ib2 - ib1 + 1
                                for ib in xrange(ib_win):
                                    p = shell.proj_win[ion, isp, ik, ilm, ib]
                                    f.write("{0:16.10f}{1:16.10f}\n".format(p.real, p.imag))
                                f.write("\n")


################################################################################
#
# output_as_text
#
################################################################################
def output_as_text(pars, el_struct, pshells, pgroups):
    """
    Output all information necessary for the converter as text files.
    """
    ctrl_output(pars, el_struct, len(pgroups))
    plo_output(pars, el_struct, pshells, pgroups)

