
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
    plovasp.plotools
    ================

    Set of routines for processing and outputting PLOs.

    This is the main module containing routines responsible for checking
    the consistency of the input data, generation of projected localized
    orbitals (PLOs) out of raw VASP projectors, and outputting data
    required by DFTTools.

    The first step of PLO processing is to select subsets of projectors
    corresponding to PLO groups. Each group contains a set of shells. Each
    projector shell is represented by an object 'ProjectorShell' that contains
    an array of projectors and information on the shell itself (orbital number,
    ions, etc.). 'ProjectorShell's are contained in both a list of shells
    (according to the original list as read from config-file) and in a
    'ProjectorGroup' object, the latter also providing information about the
    energy window.

    Order of operations:
     - transform projectors (all bands) in each shell
     - select transformed shell projectors for a given group within the window
     - orthogonalize if necessary projectors within a group by performing
       the following operations for each k-point:
       * combine all projector shells into a single array
       * orthogonalize the array
       * distribute back the arrays assuming that the order is preserved

"""
import itertools as it
import numpy as np
from .proj_group import ProjectorGroup
from .proj_shell import ProjectorShell
from .proj_shell import ComplementShell

np.set_printoptions(suppress=True)

# 'simplejson' is supposed to be faster than 'json' in stdlib.
try:
    import simplejson as json
except ImportError:
    import json

def issue_warning(message):
    """
    Issues a warning.
    """
    print()
    print("  !!! WARNING !!!: " + message)
    print()

################################################################################
# check_data_consistency()
################################################################################
def check_data_consistency(pars, el_struct):
    """
    Check the consistency of the VASP data.
    """
# Check that ions inside each shell are of the same sort
    for sh in pars.shells:
        max_ion_index = max([max(gr) for gr in sh['ions']['ion_list']])
        assert max_ion_index < el_struct.natom, "Site index in the projected shell exceeds the number of ions in the structure"
        ion_list = list(it.chain(*sh['ions']['ion_list']))

        sorts = set([el_struct.type_of_ion[io] for io in ion_list])
        assert len(sorts) == 1, "Each projected shell must contain only ions of the same sort"

# Check that ion and orbital lists in shells match those of projectors
        lshell = sh['lshell']
        for ion in ion_list:
            for par in el_struct.proj_params:
                if par['isite'] - 1 == ion and par['l'] == lshell:
                    break
            else:
                errmsg = "Projector for isite = %s, l = %s does not match PROJCAR"%(ion + 1, lshell)
                raise Exception(errmsg)


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
# check if at least one shell is correlated
    assert np.any([shell['corr'] for shell in conf_pars.shells]), 'at least one shell has be CORR = True'
    nshell = len(conf_pars.shells)
    print()
    print("  Generating %i shell%s..."%(nshell, '' if nshell == 1 else 's'))
    pshells = []
    for sh_par in conf_pars.shells:
        pshell = ProjectorShell(sh_par, proj_raw, el_struct.proj_params, el_struct.kmesh, el_struct.structure, el_struct.nc_flag)
        print()
        print("    Shell         : %s"%(pshell.user_index))
        print("    Orbital l     : %i"%(pshell.lorb))
        print("    Number of ions: %i"%(pshell.nion))
        print("    Dimension     : %i"%(pshell.ndim))
        print("    Correlated    : %r"%(pshell.corr))
        print("    Ion sort      : %r"%(pshell.ion_sort))
        pshells.append(pshell)


    pgroups = []
    for gr_par in conf_pars.groups:
        pgroup = ProjectorGroup(gr_par, pshells, eigvals)
        pgroup.orthogonalize()
        if pgroup.complement:
            pgroup.calc_complement(eigvals)
        if conf_pars.general['hk']:
            pgroup.calc_hk(eigvals)
            #testout = 'hk.out.h5'
            #from h5 import HDFArchive
            #with HDFArchive(testout, 'w') as h5test:
            #    h5test['hk'] = pgroup.hk
# DEBUG output
        print("Density matrix:")
        nimp = 0.0
        ov_all = []
        for ish in pgroup.ishells:
            if  not isinstance(pshells[pgroup.ishells[ish]],ComplementShell):
                print("  Shell %i"%(ish + 1))
                dm_all, ov_all_ = pshells[ish].density_matrix(el_struct)
                ov_all.append(ov_all_[0])
                spin_fac = 2 if dm_all.shape[0] == 1 else 1
                for io in range(dm_all.shape[1]):
                    print("    Site %i"%(io + 1))
                    dm = spin_fac * dm_all[:, io, : ,:].sum(0)
                    for row in dm:
                        print(''.join(map("{0:14.7f}".format, row)))
                    ndm = dm.trace()
                    if pshells[ish].corr:
                        nimp += ndm
                    print("      trace: ", ndm)
        print()
        print("  Impurity density:", nimp)
        print()
        print("Overlap:")
        for io, ov in enumerate(ov_all):
            print("  Site %i"%(io + 1))
            print(ov[0,...])
        print()
        print("Local Hamiltonian:")
        for ish in pgroup.ishells:
            if  not isinstance(pshells[pgroup.ishells[ish]],ComplementShell):
                print("  Shell %i"%(ish + 1))
                loc_ham = pshells[pgroup.ishells[ish]].local_hamiltonian(el_struct)
                for io in range(loc_ham.shape[1]):
                    print("    Site %i (real | complex part)"%(io + 1))
                    for row in loc_ham[:, io, :, :].sum(0):
                        print(''.join(map("{0:14.7f}".format, row.real))+' |'+''.join(map("{0:14.7f}".format, row.imag)))
# END DEBUG output
        if 'dosmesh' in conf_pars.general:
            print()
            print("Evaluating DOS...")
            mesh_pars = conf_pars.general['dosmesh']
            if np.isnan(mesh_pars['emin']):
                dos_emin = pgroup.emin
                dos_emax = pgroup.emax
            else:
                dos_emin = mesh_pars['emin']
                dos_emax = mesh_pars['emax']
            n_points = mesh_pars['n_points']

            emesh = np.linspace(dos_emin, dos_emax, n_points)
            for ish in pgroup.ishells:
                if  not isinstance(pshells[pgroup.ishells[ish]],ComplementShell) or True:
                    print("  Shell %i"%(ish + 1))
                    dos = pshells[pgroup.ishells[ish]].density_of_states(el_struct, emesh)
                    de = emesh[1] - emesh[0]
                    ntot = (dos[1:,...] + dos[:-1,...]).sum(0) / 2 * de
                    print("    Total number of states:", ntot)
                    for io in range(dos.shape[2]):
                        np.savetxt('pdos_%i_%i.dat'%(ish,io), np.vstack((emesh.T, dos[:, 0, io, :].T)).T)

        pgroups.append(pgroup)

    return pshells, pgroups

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
    if pars.general['hk']:
        hk_output(pars, el_struct, pgroups)


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
        for ik in range(nktot):
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
            for it in range(ntet):
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

    Control file format
    """"""""""""""""""""""""""""""

    Filename '<namebase>.ctrl'. Contains the data shared between all shells.
    The JSON-header consists of the following elements:

     * *nk*: number of `k`-points

     * *ns*: number of spin channels

     * *nc_flag*: collinear/noncollinear case (False/True)

     * *ng*: number of projector groups

     * Symmetry information (list of symmetry operations)

     * *efermi*: Fermi level (optional)

     * Lattice information

    """
    ctrl_fname = conf_pars.general['basename'] + '.ctrl'
    head_dict = {}

# TODO: Add output of tetrahedra
# Construct the header dictionary
    head_dict['ngroups'] = ng
    head_dict['nk'] = el_struct.kmesh['nktot']
    head_dict['ns'] = el_struct.nspin
    head_dict['kvec1'] = list(el_struct.structure['kpt_basis'][:,0])
    head_dict['kvec2'] = list(el_struct.structure['kpt_basis'][:,1])
    head_dict['kvec3'] = list(el_struct.structure['kpt_basis'][:,2])
    head_dict['nc_flag'] = 1 if el_struct.nc_flag else 0
#    head_dict['efermi'] = conf_pars.general['efermi']  # We probably don't need Efermi

    header = json.dumps(head_dict, indent=4, separators=(',', ': '))

    print("  Storing ctrl-file...")
    with open(ctrl_fname, 'wt') as f:
        f.write(header + "\n")
        f.write("#END OF HEADER\n")

        f.write("# k-points and weights\n")
        labels = ['kx', 'ky', 'kz', 'kweight']
        out = "".join([s.center(15) for s in labels])
        f.write("#" + out + "\n")
        for ik, kp in enumerate(el_struct.kmesh['kpoints']):
            tmp1 = "".join(map("{0:15.10f}".format, kp))
            out = tmp1 + "{0:16.10f}".format(el_struct.kmesh['kweights'][ik])
            f.write(out + "\n")
        f.write("# k-points and weights cartesian\n")
        labels = ['kx', 'ky', 'kz']
        out = "".join([s.center(15) for s in labels])
        f.write("#" + out + "\n")
        for ik, kp in enumerate(el_struct.kmesh['kpoints_cart']):
            out = "".join(map("{0:15.10f}".format, kp))
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

    | # Energy window: emin, emax
    | ib_min, ib_max
    | nelect
    | # Eigenvalues
    | isp, ik1, kx, ky, kz, kweight
    | ib1, ib2
    | eig1
    | eig2
    | ...
    | eigN
    | ik2, kx, ky, kz, kweight
    | ...

    | # Projected shells
    | Nshells
    | # Shells: <shell indices>
    | # Shell <1>
    | Shell 1
    | ndim
    | # complex arrays: plo(ns, nion, ndim, nb)
    | ...
    | # Shells: <shell indices>
    | # Shell <2>
    | Shell 2
    | ...

    """
    for ig, pgroup in enumerate(pgroups):
        plo_fname = conf_pars.general['basename'] + '.pg%i'%(ig + 1)
        print("  Storing PLO-group file '%s'..."%(plo_fname))
        head_dict = {}


        head_dict['nb_max'] = pgroup.nb_max

        if 'bands' in conf_pars.groups[ig]:
            head_dict['bandwindow'] = (pgroup.ib_min, pgroup.ib_max)
        else:
            head_dict['ewindow'] = (pgroup.emin, pgroup.emax)

# Number of electrons within the window
        head_dict['nelect'] = pgroup.nelect_window(el_struct)
        print("  Density within window:", head_dict['nelect'])

        head_shells = []
        for ish in pgroup.ishells:
            shell = pgroup.shells[ish]
            sh_dict = {}
            sh_dict['shell_index'] = ish
            sh_dict['lorb'] = shell.lorb
            sh_dict['ndim'] = shell.ndim
            sh_dict['corr'] = shell.corr
# Convert ion indices from the internal representation (starting from 0)
# to conventional VASP representation (starting from 1)
            ion_output = [io + 1 for io in shell.ion_list]
# Derive sorts from equivalence classes
            sh_dict['ion_list'] = ion_output
            sh_dict['ion_sort'] = shell.ion_sort

# TODO: add the output of transformation matrices

            head_shells.append(sh_dict)

        head_dict['shells'] = head_shells

        header = json.dumps(head_dict, indent=4, separators=(',', ': '))

        with open(plo_fname, 'wt') as f:
            f.write(header + "\n")
            f.write("#END OF HEADER\n")

# Eigenvalues within the window
            if 'bands' in conf_pars.groups[ig]:
                f.write("# Eigenvalues within the band window: %s, %s\n"%(pgroup.ib_min+1, pgroup.ib_max+1))
            else:
                f.write("# Eigenvalues within the energy window: %s, %s\n"%(pgroup.emin, pgroup.emax))

            nk, nband, ns_band = el_struct.eigvals.shape
            for isp in range(ns_band):
                f.write("# is = %i\n"%(isp + 1))
                for ik in range(nk):
                    ib1, ib2 = pgroup.ib_win[ik, isp, 0], pgroup.ib_win[ik, isp, 1]
# Output band indices in Fortran convention!
                    f.write(" %i  %i\n"%(ib1 + 1, ib2 + 1))
                    for ib in range(ib1, ib2 + 1):
                        eigv_ef = el_struct.eigvals[ik, ib, isp] - el_struct.efermi
                        f_weight = el_struct.ferw[isp, ik, ib]
                        f.write("%13.8f %12.7f\n"%(eigv_ef, f_weight))

# Projected shells
            f.write("# Projected shells\n")
            f.write("# Shells: %s\n"%(pgroup.ishells))
            for ish in pgroup.ishells:
                shell = pgroup.shells[ish]
                f.write("# Shell %i\n"%(ish))

                nion, ns, nk, nlm, nb = shell.proj_win.shape
                for isp in range(ns):
                    f.write("# is = %i\n"%(isp + 1))
                    for ik in range(nk):
                        f.write("# ik = %i\n"%(ik + 1))
                        for ion in range(nion):
                            for ilm in range(nlm):
                                ib1, ib2 = pgroup.ib_win[ik, isp, 0], pgroup.ib_win[ik, isp, 1]
                                ib_win = ib2 - ib1 + 1
                                for ib in range(ib_win):
                                    p = shell.proj_win[ion, isp, ik, ilm, ib]
                                    f.write("{0:16.10f}{1:16.10f}\n".format(p.real, p.imag))
                                f.write("\n")

################################################################################
#
# plo_output
#
################################################################################
def hk_output(conf_pars, el_struct, pgroups):
    """
    Outputs HK into text file.

    Filename is defined by <basename> that is passed from config-file.

    The Hk for each groups is stored in a '<basename>.hk<Ng>' file. The format is
    similar as defined in the Hk dft_tools format, but does not store info
    about correlated shells and irreps

    nk                  # number of k-points
    n_el                # electron density
    n_sh                # number of total atomic shells
    at sort l dim       # atom, sort, l, dim
    at sort l dim       # atom, sort, l, dim

    After these header lines, the file has to contain the Hamiltonian matrix
    in orbital space. The standard convention is that you give for each k-point
    first the matrix of the real part, then the matrix of the imaginary part,
    and then move on to the next k-point.

    """


    for ig, pgroup in enumerate(pgroups):

        hk_fname = conf_pars.general['basename'] + '.hk%i'%(ig + 1)
        print("  Storing HK-group file '%s'..."%(hk_fname))

        head_shells = []
        for ish in pgroup.ishells:

            shell = pgroup.shells[ish]

            ion_output = [io + 1 for io in shell.ion_list]

            for iion in ion_output:
                sh_dict = {}
                sh_dict['shell_index'] = ish
                sh_dict['lorb'] = shell.lorb
                sh_dict['ndim'] = shell.ndim
    # Convert ion indices from the internal representation (starting from 0)
    # to conventional VASP representation (starting from 1)

    # Derive sorts from equivalence classes
                sh_dict['ion_list'] = ion_output
                sh_dict['ion_sort'] = shell.ion_sort


                head_shells.append(sh_dict)

        with open(hk_fname, 'wt') as f:
            # Eigenvalues within the window
            nk, nband, ns_band = el_struct.eigvals.shape
            f.write('%i          # number of kpoints\n'%nk)
            f.write('{0:0.4f}      # electron density\n'.format(pgroup.nelect_window(el_struct)))
            f.write('%i            # number of shells\n'%len(head_shells))
            for head in head_shells:
                f.write('%i %i %i %i      # atom sort l dim\n'%(head['ion_list'][0],head['ion_sort'][0],head['lorb'],head['ndim']))

            norbs = pgroup.hk.shape[2]
            for isp in range(ns_band):
                for ik in range(nk):
                    for io in range(norbs):
                        for iop in range(norbs):
                            f.write(" {0:14.10f}".format(pgroup.hk[isp,ik,io,iop].real))
                        f.write("\n")
                    for io in range(norbs):
                        for iop in range(norbs):
                            f.write(" {0:14.10f}".format(pgroup.hk[isp,ik,io,iop].imag))
                        f.write("\n")
