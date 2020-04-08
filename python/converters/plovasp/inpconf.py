
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
  plovasp.inpconfig
  =================

  Module for parsing and checking an input config-file.
"""
import configparser
import numpy as np
import re
import sys
import itertools as it
from . import vaspio

def issue_warning(message):
    """
    Issues a warning.
    """
    print()
    print("  !!! WARNING !!!: " + message)
    print()

################################################################################
################################################################################
#
# class ConfigParameters
#
################################################################################
################################################################################
class ConfigParameters:
    r"""
    Class responsible for parsing of the input config-file.

    Parameters:

    - *sh_required*, *sh_optional* : required and optional parameters of shells
    - *gr_required*, *gr_optional* : required and optional parameters of groups

    The dictionary contains a mapping of conf-file keywords to
    a pair of objects:

      1. internal name of a parameter
      2. function used to convert an input string into data for a given parameter
    """
################################################################################
#
# __init__()
#
################################################################################
    def __init__(self, input_filename, verbosity=1):
        self.verbosity = verbosity
        self.cp = configparser.SafeConfigParser()
        self.cp.readfp(open(input_filename, 'r'))

        self.parameters = {}

        self.sh_required = {
            'ions': ('ions', self.parse_string_ion_list),
            'lshell': ('lshell', int)}

        self.sh_optional = {
            'transform': ('tmatrix', lambda s: self.parse_string_tmatrix(s, real=True)),
            'transfile': ('tmatrices', self.parse_file_tmatrix),
            'sort': ('ion_sort', self.parse_string_int,None),
            'corr': ('corr', self.parse_string_logical, True)}

        self.gr_required = {
            'shells': ('shells', lambda s: list(map(int, s.split()))),
            'ewindow': ('ewindow', self.parse_energy_window)}

        self.gr_optional = {
            'normalize' : ('normalize', self.parse_string_logical, True),
            'normion' : ('normion', self.parse_string_logical, True),
            'complement' : ('complement', self.parse_string_logical, False),
            'bands': ('bands', self.parse_band_window)}


        self.gen_optional = {
            'basename' : ('basename', str, 'vasp'),
            'efermi' : ('efermi', float),
            'dosmesh': ('dosmesh', self.parse_string_dosmesh),
            'hk': ('hk', self.parse_string_logical, False)}

#
# Special parsers
#
################################################################################
#
# parse_string_ion_list()
#
################################################################################
    def parse_string_ion_list(self, par_str):
        """
        The ion list accepts the following formats:
          1). A list of ion indices according to POSCAR.
              The list can be defined as a range '9..20'.

          2). A list of ion groups (e.g. '[1  4]  [2  3]') in which
              case each group defines a set of equivalent sites.

          3). An element name, in which case all ions with
              this name are included. NOT YET IMPLEMENTED.

        The second option requires an input from POSCAR file.
        """
        ion_info = {}

# First check if a range is given
        patt = '([0-9]+)\.\.([0-9]+)'
        match = re.match(patt, par_str)
        if match:
            i1, i2 = tuple(map(int, match.groups()[:2]))
            mess = "First index of the range must be smaller or equal to the second"
            assert i1 <= i2, mess
# Note that we need to subtract 1 from VASP indices
            ion_info['ion_list'] = [[ion - 1] for ion in range(i1, i2 + 1)]
            ion_info['nion'] = len(ion_info['ion_list'])
        else:
# Check if a set of indices is given
            try:
                l_tmp = list(map(int, par_str.split()))
                l_tmp.sort()
# Subtract 1 so that VASP indices (starting with 1) are converted
# to Python indices (starting with 0)
                ion_info['ion_list'] = [[ion - 1] for ion in l_tmp]
                ion_info['nion'] = len(ion_info['ion_list'])
            except ValueError:
                pass
# Check if equivalence classes are given

        if not ion_info:
            try:
                patt = '[0-9][0-9,\s]*'
                patt2 = '[0-9]+'
                classes = re.findall(patt, par_str)
                ion_list = []
                nion = 0
                for cl in classes:
                    ions = list(map(int, re.findall(patt2, cl)))
                    ion_list.append([ion - 1 for ion in ions])
                    nion += len(ions)

                if not ion_list:
                    raise ValueError

                ion_info['ion_list'] = ion_list
                ion_info['nion'] = nion
            except ValueError:
                err_msg = "Error parsing list of ions"
                raise NotImplementedError(err_msg)

        if 'ion_list' in ion_info:
            ion_list = ion_info['ion_list']

            assert all([all([ion >= 0 for ion in gr]) for gr in ion_list]), (
               "Lowest ion index is smaller than 1 in '%s'"%(par_str))

        return ion_info

################################################################################
#
# parse_string_logical()
#
################################################################################
    def parse_string_logical(self, par_str):
        """
        Logical parameters are given by string 'True' or 'False'
        (case does not matter). In fact, only the first symbol matters so that
        one can write 'T' or 'F'.
        """
        first_char = par_str[0].lower()
        assert first_char in 'tf', "Logical parameters should be given by either 'True' or 'False'"
        return first_char == 't'


################################################################################
#
# parse_string_int()
#
################################################################################
    def parse_string_int(self, par_str):
        """
        int parameters
        """
        return int(par_str)

################################################################################
#
# parse_energy_window()
#
################################################################################
    def parse_energy_window(self, par_str):
        """
        Energy window is given by two floats, with the first one being smaller
        than the second one.
        """
        ftmp = list(map(float, par_str.split()))
        assert len(ftmp) == 2, "EWINDOW must be specified by exactly two floats"
        assert ftmp[0] < ftmp[1], "The first float in EWINDOW must be smaller than the second one"
        return tuple(ftmp)

################################################################################
#
# parse_band_window()
#
################################################################################
    def parse_band_window(self, par_str):
        """
        Band window is given by two ints, with the first one being smaller
        than the second one.
        """
        ftmp = list(map(int, par_str.split()))
        assert len(ftmp) == 2, "BANDS must be specified by exactly two ints"
        assert ftmp[0] < ftmp[1], "The first int in BANDS must be smaller than the second one"
        return tuple(ftmp)

################################################################################
#
# parse_string_tmatrix()
#
################################################################################
    def parse_string_tmatrix(self, par_str, real):
        """
        Transformation matrix is defined as a set of rows separated
        by a new line symbol.
        """
        str_rows = par_str.split('\n')
        try:
            rows = [list(map(float, s.split())) for s in str_rows]
        except ValueError:
            err_mess = "Cannot parse a matrix string:\n%s"%(par_str)
            raise ValueError(err_mess)

        nr = len(rows)
        nm = len(rows[0])

        err_mess = "Number of columns must be the same:\n%s"%(par_str)
        for row in rows:
            assert len(row) == nm, err_mess

        if real:
            mat = np.array(rows)
        else:
            err_mess = "Complex matrix must contain 2*M values:\n%s"%(par_str)
            assert 2 * (nm // 2) == nm, err_mess

            tmp = np.array(rows, dtype=np.complex128)
            mat = tmp[:, 0::2] + 1.0j * tmp[:, 1::2]

        return mat

################################################################################
#
# parse_file_tmatrix()
#
################################################################################
    def parse_file_tmatrix(self, filename):
        """
        Parses a file 'filename' containing transformation matrices
        for each ion. The parser returns a raw matrix that will be
        interpreted elsewhere because the interpretation depends on
        shell parameters.
        """
        tmatrices = np.loadtxt(filename)
        return tmatrices

################################################################################
#
# parse_string_dosmesh()
#
################################################################################
    def parse_string_dosmesh(self, par_str):
        """
        Two formats are accepted:

          1. Two floats (energy range) and an integer (number of energy points).

          2. One integer (number of energy points). In this case the energy
             range is taken to be equal to EMIN, EMAX of a shell.

        The parser returns a dictionary:
          {'n_points': int,
           'emin': float,
           'emax': float}

        If the second option is used, 'emin' and 'emax' are undefined
        and set to 'nan'.
        """
        stmp = par_str.split()
        if len(stmp) == 3:
            emin, emax = float(stmp[0]), float(stmp[1])
            n_points = int(stmp[2])
        elif len(stmp) == 1:
            n_points = int(stmp[0])
            emin = emax = float('nan')
        else:
            err_mess = "DOSMESH must be either 'EMIN EMAX NPOINTS' or 'NPOINTS'"
            raise ValueError(err_mess)

        dos_pars = {
            'n_points': n_points,
            'emin': emin,
            'emax': emax}

        return dos_pars

################################################################################
#
# parse_parameter_set()
#
################################################################################
    def parse_parameter_set(self, section, param_set, exception=False, defaults=True):
        """
        Parses required or optional parameter set from a section.
        For required parameters `exception=True` must be set.
        """
        parsed = {}
        for par in list(param_set.keys()):
            key = param_set[par][0]
            try:
                par_str = self.cp.get(section, par)
            except (configparser.NoOptionError, configparser.NoSectionError):
                if exception:
                    message = "Required parameter '%s' not found in section [%s]"%(par, section)
                    raise Exception(message)
                else:
# Use the default value if there is one
                    if defaults and len(param_set[par]) > 2:
                        parsed[key] = param_set[par][2]
                    continue

            if self.verbosity > 0:
                print("  %s = %s"%(par, par_str))

            parse_fun = param_set[par][1]
            parsed[key] = parse_fun(par_str)

        return parsed


################################################################################
#
# parse_shells()
#
################################################################################
    def parse_shells(self):
        """
        Parses all [Shell] sections.
        """
# Find all [Shell] sections
# (note that ConfigParser transforms all names to lower case)
        sections = self.cp.sections()

        sh_patt1 = re.compile('shell +.*', re.IGNORECASE)
        sec_shells = list(filter(sh_patt1.match, sections))

        self.nshells = len(sec_shells)
        assert self.nshells > 0, "No projected shells found in the input file"

        if self.verbosity > 0:
            print()
            if self.nshells > 1:
                print("  Found %i projected shells"%(self.nshells))
            else:
                print("  Found 1 projected shell")

# Get shell indices
        sh_patt2 = re.compile('shell +([0-9]*)$', re.IGNORECASE)
        try:
            get_ind = lambda s: int(sh_patt2.match(s).groups()[0])
            sh_inds = list(map(get_ind, sec_shells))
        except (ValueError, AttributeError):
            raise ValueError("Failed to extract shell indices from a list: %s"%(sec_shells))

        self.sh_sections = {ind: sec for ind, sec in zip(sh_inds, sec_shells)}

# Check that all indices are unique
# In principle redundant because the list of sections will contain only unique names
        assert len(sh_inds) == len(set(sh_inds)), "There must be no shell with the same index!"

# Ideally, indices should run from 1 to <nshells>
# If it's not the case, issue a warning
        sh_inds.sort()
        if sh_inds != list(range(1, len(sh_inds) + 1)):
            issue_warning("Shell indices are not uniform or not starting from 1. "
               "This might be an indication of a incorrect setup.")

# Parse shell parameters and put them into a list sorted according to the original indices
        self.shells = []
        for ind in sh_inds:
            shell = {}
# Store the original user-defined index
            shell['user_index'] = ind
            section = self.sh_sections[ind]

            if self.verbosity > 0:
                print()
                print("  Shell parameters:")
# Shell required parameters
            parsed = self.parse_parameter_set(section, self.sh_required, exception=True)
            shell.update(parsed)

# Shell optional parameters
            parsed = self.parse_parameter_set(section, self.sh_optional, exception=False)
            shell.update(parsed)

# Group required parameters
# Must be given if no group is explicitly specified
# If in conflict with the [Group] section, the latter has a priority
            parsed = self.parse_parameter_set(section, self.gr_required, exception=False)
            shell.update(parsed)

# Group optional parameters
            parsed = self.parse_parameter_set(section, self.gr_optional, exception=False, defaults=False)
            shell.update(parsed)

            self.shells.append(shell)

################################################################################
#
# parse_groups()
#
################################################################################
    def parse_groups(self):
        """
        Parses [Group] sections.
        """
# Find group sections
        sections = self.cp.sections()

        gr_patt = re.compile('group +(.*)', re.IGNORECASE)
        sec_groups = list(filter(gr_patt.match, sections))

        self.ngroups = len(sec_groups)

        self.groups = []
# Parse group parameters
        for section in sec_groups:
            group = {}

# Extract group index (FIXME: do we really need it?)
            gr_patt2 = re.compile('group +([0-9]*)$', re.IGNORECASE)
            try:
                gr_ind = int(gr_patt2.match(section).groups()[0])
            except (ValueError, AttributeError):
                raise ValueError("Failed to extract group index from a group name: %s"%(section))
            group['index'] = gr_ind

            if self.verbosity > 0:
                print()
                print("  Group parameters:")
# Group required parameters
            parsed = self.parse_parameter_set(section, self.gr_required, exception=True)
            group.update(parsed)

# Group optional parameters
            parsed = self.parse_parameter_set(section, self.gr_optional, exception=False)
            group.update(parsed)

            self.groups.append(group)

# Sort groups according to indices defined in the config-file
        if self.ngroups > 0:
            self.groups.sort(key=lambda g: g['index'])

################################################################################
#
# groups_shells_consistency()
#
################################################################################
    def groups_shells_consistency(self):
        """
        Ensures consistency between groups and shells. In particular:
            - if no groups are explicitly defined and only shell is defined create a group automatically
            - check the existance of all shells referenced in the groups
            - check that all shells are referenced in the groups
            
        """
# Special case: no groups is defined
        if self.ngroups == 0:
# Check that 'nshells = 1'
            assert self.nshells == 1, "At least one group must be defined if there are more than one shells."

# Otherwise create a single group taking group information from [Shell] section
            self.groups.append({})
            self.groups[0]['index'] = '1'
# Check that the single '[Shell]' section contains enough information
# (required group parameters except 'shells')
# and move it to the `groups` dictionary
            sh_gr_required = dict(self.gr_required)
            sh_gr_required.pop('shells')
            try:
                for par in list(sh_gr_required.keys()):
                    key = sh_gr_required[par][0]
                    value = self.shells[0].pop(key)
                    self.groups[0][key] = value
            except KeyError:
                message = "One [Shell] section is specified but no explicit [Group] section is provided."
                message += " In this case the [Shell] section must contain all required group information.\n"
                message += "  Required parameters are: %s"%(list(sh_gr_required.keys()))
                raise KeyError(message)

# Do the same for optional group parameters, but do not raise an exception
            for par in list(self.gr_optional.keys()):
                try:
                    key = self.gr_optional[par][0]
                    value = self.shells[0].pop(key)
                    self.groups[0][key] = value
                except KeyError:
                    if len(self.gr_optional[par]) > 2:
                        self.groups[0][key] = self.gr_optional[par][2]
                    continue
# Add the index of the single shell into the group
            self.groups[0].update({'shells': [1]})

#
# Consistency checks
#
# Check the existence of shells referenced in the groups
        def find_shell_by_user_index(uindex):
            for ind, shell in enumerate(self.shells):
                if shell['user_index'] == uindex:
                    return ind, shell
            raise KeyError

        sh_all_inds = []
        for group in self.groups:
            gr_shells = group['shells']
            sh_inds = []
            for user_ind in gr_shells:
                try:
                    ind, shell = find_shell_by_user_index(user_ind)
                except KeyError:
                    raise Exception("Shell %i referenced in group '%s' does not exist"%(user_ind, group['index']))
                sh_inds.append(ind)

# If [Shell] section contains (potentially conflicting) group parameters
# remove them and issue a warning.
#
# First, required group parameters
                for par in list(self.gr_required.keys()):
                    try:
                        key = self.gr_required[par][0]
                        value = shell.pop(key)
                        mess = ("  Redundant group parameter '%s' in [Shell] section"
                                " %i is discarded"%(par, user_ind))
                        issue_warning(mess)
                    except KeyError:
                        continue

# Second, optional group parameters
                for par in list(self.gr_optional.keys()):
                    try:
                        key = self.gr_optional[par][0]
                        value = shell.pop(key)
                        mess = ("  Redundant group parameter '%s' in [Shell] section"
                                " %i is discarded"%(par, user_ind))
                        issue_warning(mess)
                    except KeyError:
                        continue

            sh_all_inds += sh_inds
# Replace user shell indices with internal ones
            group['shells'] = sh_inds

        sh_refs_used = list(set(sh_all_inds))
        sh_refs_used.sort()

# Check that all shells are referenced in the groups
        assert sh_refs_used == list(range(self.nshells)), "Some shells are not inside any of the groups"


################################################################################
#
# parse_general()
#
################################################################################
    def parse_general(self):
        """
        Parses [General] section.
        """
        self.general = {}
        sections = self.cp.sections()
        gen_section = [s for s in sections if s.lower() == 'general']
# If no [General] section is found parse a dummy section name to the parser
# to reset parameters to their default values
        if len(gen_section) > 1:
            raise Exception("More than one section [General] is found")
        if len(gen_section) == 0:
            gen_section = 'general'
        gen_section = gen_section[0]
        parsed = self.parse_parameter_set(gen_section, self.gen_optional, exception=False)
        self.general.update(parsed)

################################################################################
#
# Main parser function
#
################################################################################
    def parse_input(self):
        """
        Parses input conf-file.
        """
        self.parse_general()
        self.parse_shells()
        self.parse_groups()

        self.groups_shells_consistency()

#
# Obsolete part
#
if __name__ == '__main__':
    narg = len(sys.argv)
    if narg < 2:
        raise SystemExit("  Usage: python pyconf.py <conf-file> [<path-to-vasp-calcultaion>]")
    else:
        filename = sys.argv[1]
        if narg > 2:
            vasp_dir = sys.argv[2]
            if vasp_dir[-1] != '/':
                vasp_dir += '/'
        else:
            vasp_dir = './'


#    plocar = vaspio.Plocar()
#    plocar.from_file(vasp_dir)
#    poscar = vaspio.Poscar()
#    poscar.from_file(vasp_dir)
#    kpoints = vaspio.Kpoints()
#    kpoints.from_file(vasp_dir)
    eigenval = vaspio.Eigenval()
    eigenval.from_file(vasp_dir)
    doscar = vaspio.Doscar()
    doscar.from_file(vasp_dir)
#    pars = parse_input(filename)
