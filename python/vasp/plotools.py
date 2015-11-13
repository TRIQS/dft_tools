
import itertools as it
import numpy as np
from proj_group import ProjectorGroup
from proj_shell import ProjectorShell

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
    print
    print "  !!! WARNING !!!: " + message
    print

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
        pshell = ProjectorShell(sh_par, proj_raw, el_struct.proj_params, el_struct.nc_flag)
        print
        print "    Shell         : %s"%(pshell.user_index)
        print "    Orbital l     : %i"%(pshell.lorb)
        print "    Number of ions: %i"%(len(pshell.ion_list))
        print "    Dimension     : %i"%(pshell.ndim)
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

