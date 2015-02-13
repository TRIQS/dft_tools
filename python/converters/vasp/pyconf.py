
import ConfigParser
import numpy as np
import re
import sys
import vaspio

def issue_warning(message):
    """
    Issues a warning.
    """
    print
    print "  !!! WARNING !!!: " + message
    print
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
    self.sh_required = {
        'ions': ('ion_list', self.parse_ion_list),
        'lshell': ('lshell', int)}

    self.sh_optional = {
        'rtransform': ('tmatrix', lambda s: self.parse_tmatrix(s, real=True)),
        'ctransform': ('tmatrix', lambda s: self.parse_tmatrix(s, real=False))}

    self.gr_required = {
        'emin': ('emin', float),
        'emax': ('emax', float)}

    self.gr_optional = {
        'normalize' : ('normalize', self.parse_logical),
        'normion' : ('normion', self.parse_logical)}


################################################################################
#
# __init__()
#
################################################################################
    def __init__(self, input_filename, verbosity=1):
        self.verbosity = verbosity
        self.cp = ConfigParser.ConfigParser()
        self.cp.readfp(open(input_filename, 'r'))

        self.conf_pars = {}

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
        The ion list accepts two formats:
          1). A list of ion indices according to POSCAR.
          2). An element name, in which case all ions with
              this name are included.

        The second option requires an input from POSCAR file.
        """
        try:
            l_tmp = map(int, par_str.split())
# Subtract 1 so that VASP indices (starting with 1) are converted
# to Python indices (starting with 0)
            ion_list = np.array(l_tmp) - 1
        except ValueError:
            err_msg = "Only an option with a list of ion indices is implemented"
            raise NotImplementedError(err_msg)

        return ion_list

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
        assert par_str[0] in 'tf', "Logical parameters should be given by either 'True' or 'False'"
        return par_str[0] == 't'

################################################################################
#
# parse_parameter_set()
#
################################################################################
    def parse_parameter_set(self, section, param_set, exception=False):
        """
        Parses required or optional parameter set from a section.
        For required parameters `exception=True` must be set.
        """
        parsed = {}
        for par in param_set.keys():
            try:
                par_str = self.cp.get(section, par)
            except ConfigParser.NoOptionError:
                if exception:
                    message = "Required parameter '%s' not found in section [%s]"%(par, section)
                    raise ConfigParser.NoOptionError(message)
                else:
                    continue

            if self.verbosity > 0:
                print "  %s = %s"%(par, par_str)

            key = param_set[par][0]
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

        sh_patt = 'shell *([0-9]*)'
        ismatch = lambda s: not re.match(sh_patt, s) is None
        sec_shells = filter(ismatch, sections)

        self.nshells = len(sec_shells)
        assert self.nshells > 0, "No projected shells found in the input file"

        if self.verbosity > 0:
            print
            if self.nshells > 1:
                print "  Found %i projected shells"%(self.nshells)
            else:
                print "  Found 1 projected shell"

# Get shell indices
        get_ind = lambda s: int(re.match(sh_patt, s).groups()[0])
        try:
            sh_inds = map(get_ind, sec_shells)
        except ValueError:
            raise ValueError("Failed to extract shell indices from a list: %s"%(sec_shells))

        self.sh_sections = {ind: sec for ind, sec in sh_inds, sec_shells}

# Check that all indices are unique
        assert len(sh_inds) == len(set(sh_inds)), "There must be no shell with the same index!"

# Ideally, indices should run from 1 to <nshells>
# If it's not the case, issue a warning
        sh_inds.sort()
        if sh_inds != range(1, len(sh_inds) + 1):
            issue_warning("Shell indices are not uniform or not starting from 1. "
               "This might be an indication of a incorrect setup."

# Parse shell parameters
        self.shells = {}
        for ind in sh_inds:
            self.shells[ind] = {}
            section = self.sh_sections[ind]

# Shell required parameters
            if self.verbosity > 0:
                print
                print "  Required shell parameters:"
            parsed = self.parse_parameter_set(section, self.sh_required, exception=True)
            self.shells[ind].update(parsed)

# Shell optional parameters
            if self.verbosity > 0:
                print
                print "  Optional shell parameters:"
            parsed = self.parse_parameter_set(section, self.sh_optional, exception=False)
            self.shells[ind].update(parsed)

# Group required parameters
# Must be given if no group is explicitly specified
# If in conflict with the [Group] section, the latter has a priority
            if self.verbosity > 0:
                print
                print "  Required group parameters:"
            parsed = self.parse_parameter_set(section, self.gr_required, exception=False)
            self.shells[ind].update(parsed)

# Group optional parameters
            if self.verbosity > 0:
                print
                print "  Optional group parameters:"
            parsed = self.parse_parameter_set(section, self.gr_optional, exception=False)
            self.shells[ind].update(parsed)

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

        gr_patt = 'group *([0-9]*)'
        ismatch = lambda s: not re.match(gr_patt, s) is None
        sec_groups = filter(ismatch, sections)

        self.ngroups = len(sec_groups)

# Special case: no groups is defined
        if self.ngroups == 0:
# Check that 'nshells = 1'
            assert self.nshells == 1, "At least one group must be defined if there are more than one shells."

# Otherwise create a single group taking group information from [Shell] section
            self.groups = [{}]
# Check that the single '[Shell]' section contains enough information
# and move it to the `groups` dictionary
            ind = self.sh_sections.keys()[0]
            try:
                for par in self.gr_required.keys():
                    key = self.gr_required[par][0]
                    value = self.shells[ind].pop(key)
                    self.groups[0][key] = value
            except KeyError:
                message = "One [Shell] section is specified but no explicit [Group] section is provided."
                message += " In this case the [Shell] section must contain all required group information.\n"
                message += "  Required parameters are: %s"%(self.gr_required.keys())
                raise KeyError(message)

# Do the same for optional group parameters, but do not raise an exception
            for par in self.gr_required.keys():
                try:
                    key = self.gr_required[par][0]
                    value = self.shells[ind].pop(key)
                    self.groups[0][key] = value
                except KeyError:
                    continue
                
            self.groups.update({'shells': [self.shells[ind]]})
            

################################################################################
#
# Main parser
# parse_logical()
#
################################################################################
    def parse_input(self):
        """
        Parses input conf-file.
        """
        self.parse_shells()
        self.parse_groups()


# Output list of dictionaries
        output_pars = [{} for isec in xrange(nsections)]
        for isec, section in enumerate(sections):
            print "Section: %s"%(section)
            for par in required.keys():
                try:
                    par_str = cp.get(section, par)
                except ConfigParser.NoOptionError:
                    raise SystemExit("*** Error: Required entry '%s' not found in the input file"%(par))

                print "  %s: %s"%(par, par_str)
                key = required[par][0]
                parse_fun = required[par][1]
                output_pars[isec][key] = parse_fun(par_str)

        print output_pars
        print cp.get(section, 'rtransform').strip().split('\n')

        return output_pars

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

