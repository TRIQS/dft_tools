.. _convVASP:

===================
Interface with VASP
===================

The VASP interface relies on new options introduced since version 5.4.x In
particular, a new INCAR-option `LOCPROJ
<https://cms.mpi.univie.ac.at/wiki/index.php/LOCPROJ>`_, the new `LORBIT` modes
13 and 14 have been added, and the new `ICHARG` mode 5 for charge
self-consistent DFT+DMFT calculations have been added.

The VASP interface methodologically builds on the so called projection on
localized orbitals (PLO) scheme, where the resulting KS states from DFT are
projected on localized orbitals, which defines a basis for setting up a
Hubbard-like model Hamiltonian. Resulting in lattice object stored in `SumkDFT`.
The implementation is presented in `M. Schüler et al. 2018 J. Phys.: Condens.
Matter 30 475901 <https://doi.org/10.1088/1361-648X/aae80a>`_.

The interface consists of two parts, :ref:`PLOVASP<refPLOVASP>`, a collection of
python classes and functions converting the raw VASP output to proper projector
functions, and the python based :ref:`VaspConverter<refVASPconverter>`, which
creates a h5 archive from the :ref:`PLOVASP<refPLOVASP>` output readable by
`SumkDFT`. Therefore, the conversion consist always of two steps.

Here, we will present a guide how the interface `can` be used to create input for a DMFT calculation, using SrVO3 as an example. Full examples can be found in the :ref:`tutorial section of DFTTools<tutorials>`.

Limitations of the interface
============================

* The interface works correctly only if the k-point symmetries
  are turned off during the VASP run (ISYM=-1).
* Generation of projectors for k-point lines (option `Lines` in KPOINTS)
  needed for Bloch spectral function calculations is not possible at the moment.
* The interface currently supports only collinear-magnetism calculation
  (this implies no spin-orbit coupling) and spin-polarized projectors have not
  been tested.
* The converter needs the correct Fermi energy from VASP, which is read from
  the LOCPROJ file. However, VASP by default does not output this information.
  Please see `Remarks on the VASP version`_.

VASP: generating raw projectors
===============================

The VASP **INCAR** option `LOCPROJ` selects a set of localized projectors that
will be written to the file **LOCPROJ** after a successful VASP run. A projector set is specified by site indices, labels of the target local states, and projector type:

  | `LOCPROJ = <sites> : <shells> : <projector type>`

where `<sites>` represents a list of site indices separated by spaces, with the
indices corresponding to the site position in the **POSCAR** file; `<shells>`
specifies local states (see below); `<projector type>` chooses a particular type
of the local basis function. The recommended projector type is `Pr 2`. This will
perform a projection of the Kohn-Sham states onto the VASP PAW projector
functions. The number specified behind `Pr` is selecting a specific PAW channel,
see the `VASP wiki page <https://cms.mpi.univie.ac.at/wiki/index.php/LOCPROJ>`_
for more information. The formalism for this type of projectors is presented in
`M. Schüler et al. 2018 J. Phys.: Condens. Matter 30 475901
<https://doi.org/10.1088/1361-648X/aae80a>`_. For further details on the
`LOCPROJ` flag also have a look in the `VASP wiki
<https://cms.mpi.univie.ac.at/wiki/index.php/LOCPROJ>`_.

The allowed labels of the local states defined in terms of cubic
harmonics are (mind the order):

 * Entire shells: `s`, `p`, `d`, `f`

 * `p`-states: `py`, `pz`, `px`

 * `d`-states: `dxy`, `dyz`, `dz2`, `dxz`, `dx2-y2`

 * `f`-states: `fy(3x2-y2)`, `fxyz`, `fyz2`, `fz3`,
   `fxz2`, `fz(x2-y2)`, `fx(x2-3y2)`.

For projector type `Pr`, one should ideally also set `LORBIT = 14` in the
INCAR file and provide parameters `EMIN`, `EMAX`, defining, in this case, an
energy range (energy window) corresponding to the valence states. Note that,
as in the case of a DOS calculation, the position of the valence states
depends on the Fermi level, which can usually be found at the end of the
OUTCAR file. Setting `LORBIT=14` will perform an automatic optimization of
the PAW projector channel as described in `M. Schüler et al. 2018 J. Phys.:
Condens. Matter 30 475901 <https://doi.org/10.1088/1361-648X/aae80a>`_, by
using a linear combination of the PAW channels, to maximize the overlap in
the chosen energy window between the projector and the Kohn-Sham state.
Therefore, setting `LORBIT=14` will let VASP ignore the channel specified
after `Pr`. This optimization is only performed for the projector type `Pr`,
not for `Ps` and obviously not for `Hy`. We recommend to specify the PAW
channel anyway, in case one forgets to set `LORBIT=14`.

In case of SrVO3 one may first want to perform a self-consistent
calculation to know the Fermi level and the rough position of the target states.
In the next step one sets `ICHARG = 1` and adds the following additional lines
into INCAR (provided that V is the second ion in POSCAR):

  | `EMIN = 3.0`
  | `EMAX = 8.0`
  | `LORBIT = 14`
  | `LOCPROJ = 2 : d : Pr 2`

The energy range does not have to be precise. Important is that it has a large
overlap with valence bands and no overlap with semi-core or high unoccupied
states. This **INCAR** will calculate and write-out projections for all five d-orbitals.


VASP input-output
-----------------

The calculated projections :math:`\langle \chi_L | \Psi_\mu \rangle` are written
into files **PROJCAR** and **LOCPROJ**. The difference between these two files
is that **LOCPROJ** contains raw matrices without any reference to
sites/orbitals, while **PROJCAR** is more detailed. In particular, the
information that can be obtained for each projector from **PROJCAR** is the
following:

  * site (and species) index
  * for each `k`-point and band: a set of complex numbers for labeled orbitals

At the same time, **LOCPROJ** contains the total number of projectors (as well
as the number of `k`-points, bands, and spin channels) in the first line, which
can be used to allocate the arrays before parsing.

Conversion for the DMFT self-consistency cycle
==============================================

The projectors generated by VASP require certain post-processing before they can
be used for DMFT calculations. The most important step is to (ortho-)normalize
them within an energy window that selects band states relevant for the impurity
problem. This will create proper Wannier functions of the initial projections
produced by VASP. Note that this energy window is different from the one
described above and it must be chosen independently of the energy range given by
`EMIN, EMAX` in the **INCAR** VASP input file. This part is done in `PLOVASP`.


PLOVASP: converting VASP output
--------------------------------

:ref:`PLOVASP<refPLOVASP>` is a collection of python functions and classes, post-processing the raw VASP `LOCPROJ` output creating proper projector functions.

The following VASP files are used by PLOVASP:
  * PROJCAR, LOCPROJ: raw projectors generated by VASP-PLO interface
  * EIGENVAL: Kohn-Sham eigenvalues as well as `k`-points with weights and Fermi weights
  * IBZKPT: `k`-point data (:math:`\Gamma`)
  * POSCAR: crystal structure data

To run `PLOVASP`, one first prepares an input file `<name>.cfg` (default name `plo.cfg`) that describes the definition of the correlated subspace. For SrVO3 this input file would look like this:

.. literalinclude:: ../tutorials/svo_vasp/plo.cfg

In the [section] general, the `DOSMESH` defines an energy window and number of
data points, which lets the converter calculate the density of states of the
created projector functions in a given energy window. Each projector shell is
defined by a section `[Shell 1]` where the number can be arbitrary and used only
for user convenience. Several parameters are required

- **IONS**: list of site indices which must be a subset of indices given earlier
  in the VASP INCAR `LOCPROJ` flag. Note: If projections are performed for
  multiple sites one can specify symmetry equivalent sites with brackets: `[2
  3]`. Here the projector are generated for ions 2 and 3, but they will be
  marked as symmetry equivalent later in 'SumkDFT'.
- **LSHELL**: :math:`l`-quantum number of the projector shell; the corresponding
  orbitals must be present in `LOCPROJ`.
- **EWINDOW**: energy window in which the projectors are normalized;
  note that the energies are defined with respect to the Fermi level.

The Option **TRANSFORM** is optional here, and it is specified to extract
only the three :math:`t_{2g}` orbitals out of the five `d` orbitals given by
:math:`l = 2`. A detailed explanation of all input parameters can be found
further below `PLOVASP detailed guide`_.

Next, the converter is executed. This can be done by calling :program:`PLOVASP` directly in the command line with the input file as an argument, e.g.:
     | `plovasp plo.cfg`

or embedded in a python script as::

    import triqs_dft_tools.converters.plovasp.converter as plo_converter
    # Generate and store PLOs
    plo_converter.generate_and_output_as_text('plo.cfg', vasp_dir='./')

This will create the xml files `vasp.ctrl` and `vasp.pg1` containing the orthonormalized projector functions readable by the :ref:`VaspConverter<refVASPconverter>`. Moreover, `PLOVASP` will output important information of the orthonormalization process, such as the density matrix of the correlated shell and the local Hamiltonian.

Running the VASP converter
-------------------------------------

The actual conversion to a h5-file is performed with the orthonormalized projector functions readable by the :ref:`VaspConverter<refVASPconverter>` in the same fashion as with the other `DFTTools` converters::

    from triqs_dft_tools.converters.vasp import *
    Converter = VaspConverter(filename = 'vasp')
    Converter.convert_dft_input()

As usual, the resulting h5-file can then be used with the SumkDFT class::
    sk = SumkDFTTools(hdf_file='vasp.h5')

Note that the automatic detection of the correct block structure might fail for
VASP inputs. This can be circumvented by setting a bigger value of the threshold
in :class:`SumkDFT <dft.sumk_dft.SumkDFT>`, e.g.::

    SK.analyse_block_structure(threshold = 1e-4)

However, this should only be done after a careful study of the density matrix and the projected DOS in the localized basis. For the complete process for SrVO3 see the tutorial for the VASP interface `here <../tutorials/svo_vasp/svo_notebook.html>`_.

PLOVASP detailed guide
======================

The general purpose of the PLOVASP tool is to transform raw, non-normalized
projectors generated by VASP into normalized projectors corresponding to
user-defined projected localized orbitals (PLO). To enhance the performance
parsing the raw VASP output files, the parser is implemented in plain C. The
idea is that the python part of the parser first reads the first line of
**LOCPROJ** and then calls the C-routine with necessary parameters to parse
**PROJCAR**. The resulting PLOs can then be used for DFT+DMFT calculations with
or without charge self-consistency. PLOVASP also provides some utilities for
basic analysis of the generated projectors, such as outputting density matrices,
local Hamiltonians, and projected density of states.

PLOs are determined by the energy window in which the raw projectors are
normalized. This allows to define either atomic-like strongly localized Wannier
functions (large energy window) or extended Wannier functions focusing on
selected low-energy states (small energy window).

In PLOVASP, all projectors sharing the same energy window are combined into a
`projector group`. This allows one in principal to define several groups with
different energy windows for the same set of raw projectors. Note: multiple groups are not yet implemented.

A set of projectors defined on sites related to each other either by symmetry
or by an atomic sort, along with a set of :math:`l`, :math:`m` quantum numbers,
forms a `projector shell`. There could be several projectors shells in a
projector group, implying that they will be normalized within the same energy
window.

Projector shells and groups are specified by a user-defined input file whose
format is described below. Additionally, each shell can be marked correlated or non-correlated, to tell `SumkDFT` whether or not these should be treated in the DMFT impurity problem.

Input file format
-----------------

The input file is written in the standard config-file format.
Parameters (or 'options') are grouped into sections specified as
`[Section name]`. All parameters must be defined inside some section.

A PLOVASP input file can contain three types of sections:

#. **[General]**: includes parameters that are independent
   of a particular projector set, such as the Fermi level, additional
   output (e.g. the density of states), etc.
#. **[Group <Ng>]**: describes projector groups, i.e. a set of
   projectors sharing the same energy window and normalization type.
   At the moment, DFTtools support only one projector group, therefore
   there should be no more than one projector group.
#. **[Shell <Ns>]**: contains parameters of a projector shell labelled
   with `<Ns>`. If there is only one group section and one shell section,
   the group section can be omitted but in this case, the group required
   parameters must be provided inside the shell section.

Section [General]
"""""""""""""""""

The entire section is optional and it contains four parameters:

*  **BASENAME** (string): provides a base name for output files.
   Default filenames are :file:`vasp.*`.
*  **DOSMESH** ([float float] integer): if this parameter is given,
   the projected density of states for each projected orbital will be
   evaluated and stored to files :file:`pdos_<s>_<n>.dat`, where `s` is the
   shell index and `n` the ion index. The energy mesh is defined by three
   numbers: `EMIN`  `EMAX`  `NPOINTS`. The first two
   can be omitted in which case they are taken to be equal to the projector
   energy window. **Important note**: at the moment this option works
   only if the tetrahedron integration method (`ISMEAR = -4` or `-5`)
   is used in VASP to produce `LOCPROJ`.
*  **EFERMI** (float): provides the Fermi level. This value overrides
   the one extracted from VASP output files.
*  **HK** (True/False): If True, the projectors are applied the the Kohn-Sham
   eigenvalues which results in a Hamitlonian H(k) in orbital basis. The H(k)
   is written for each group to a file :file:`Basename.hk<Ng>`. It is recommended
   to also set `COMPLEMENT = True` (see below). Default is False.

There are no required parameters in this section.

Section [Shell]
"""""""""""""""

This section specifies a projector shell. Each `[Shell]` section must be
labeled by an index, e.g. `[Shell 1]`. These indices can then be referenced
in a `[Group]` section.

In each `[Shell]` section two parameters are required:

*  **IONS** (list of integer): indices of sites included in the shell.
   The sites can be given either by a list of integers `IONS = 5 6 7 8`
   or by a range `IONS = 5..8`. The site indices must be compatible with
   the POSCAR file. Morever, sites can be marked to be identical by 
   grouping them with brackets, i.e. `IONS = [5 6] [7 8]` will mark the
   sites 5 and 6 in the POSCAR (and of course also 7 and 8) to be idential.
   This will mark these correlated site as equivalent, and only one 
   impurity problem per bracket group is generated.
*  **LSHELL** (integer): :math:`l` quantum number of the desired local states.

It is important that a given combination of site indices and local states
given by `LSHELL` must be present in the LOCPROJ file.

There are additional optional parameters that allow one to transform
the local states:

*  **CORR** (True/False): Determines if shell is correlated or not. At least one
   shell has to be correlated. Default is True.
*  **SORT** (integer): Overrides the default detection of ion sorts by supplying
   an integer. Default is `None`, for which the default behavior is retained.
*  **TRANSFORM** (matrix): local transformation matrix applied to all states
   in the projector shell. The matrix is defined by a (multiline) block
   of floats, with each line corresponding to a row. The number of columns
   must be equal to :math:`2 l + 1`, with :math:`l` given by `LSHELL`. Only real matrices
   are allowed. This parameter can be useful to select certain subset of
   orbitals or perform a simple global rotation.
*  **TRANSFILE** (string): name of the file containing transformation
   matrices for each site. This option allows for a full-fledged functionality
   when it comes to local state transformations. The format of this file
   is described :ref:`below <transformation_file>`.

Section [Group]
"""""""""""""""

Each defined projector shell must be part of a projector group. In the current
implementation of DFTtools only a single group (labelled by any integer, e.g. `[Group 1]`)
is supported. This implies that all projector shells
must be included in this group.

Required parameters for any group are the following:

*  **SHELLS** (list of integers): indices of projector shells included in the group.
   All defined shells must be grouped.
*  **EWINDOW** (float float): the energy window specified by two floats: bottom
   and top. All projectors in the current group are going to be normalized within
   this window. *Note*: This option must be specified inside the `[Shell]` section
   if only one shell is defined and the `[Group]` section is omitted.

Optional group parameters:

*  **NORMALIZE** (True/False): specifies whether projectors in the group are
   to be normalized. The default value is **True**.
*  **NORMION** (True/False): specifies whether projectors are normalized on
   a per-site (per-ion) basis. That is, if `NORMION = True`, the orthogonality
   condition will be enforced on each site separately but the Wannier functions
   on different sites will not be orthogonal. If `NORMION = False`, the Wannier functions
   on different sites included in the group will be orthogonal to each other. The default value is **True**
*  **BANDS** (int int): the energy window specified by two ints: band index of
   lowest band and band index of highest band. Using this overrides the selection
   in `EWINDOW`.
*  **COMPLEMENT** (True/False). If True, the orthogonal complement is calculated
   resulting in unitary (quadratic) projectors, i.e., the same number of orbitals
   as bands. It is required to have an equal number of bands in the energy window
   at each k-point. Default is False.


.. _transformation_file:

File of transformation matrices
"""""""""""""""""""""""""""""""

.. warning::
  The description below applies only to collinear cases (i.e., without spin-orbit
  coupling). In this case, the matrices are spin-independent.

The file specified by option `TRANSFILE` contains transformation matrices
for each ion.  Each line must contain a series of floats whose number is either equal to
the number of orbitals :math:`N_{orb}` (in this case the transformation matrices
are assumed to be real) or to :math:`2 N_{orb}` (for the complex transformation matrices).
The total number of lines :math:`N` must be a multiple of the number of ions :math:`N_{ion}`
and the ratio :math:`N / N_{ion}`, then, gives the dimension of the transformed
orbital space. The lines with floats can be separated by any number of empty or
comment lines (starting from `#`), which are ignored.

A very simple example is a transformation matrix that selects the :math:`t_{2g}` manifold.
For two correlated sites, one can define the file as follows:
::

  # Site 1
    1.0   0.0   0.0   0.0   0.0
    0.0   1.0   0.0   0.0   0.0
    0.0   0.0   0.0   1.0   0.0

  # Site 2
    1.0   0.0   0.0   0.0   0.0
    0.0   1.0   0.0   0.0   0.0
    0.0   0.0   0.0   1.0   0.0

Remarks on the VASP version
===============================

In the current version of the interface the Fermi energy is extracted from the
`LOCPROJ` file. The file should contain the Fermi energy in the header. One can
either copy the Fermi energy manually there after a successful VASP run, or
modify the VASP source code slightly, by replacing the following line in
`locproj.F` (around line 695):
::
  WRITE(99,'(4I6,"  # of spin, # of k-points, # of bands, # of proj" )') NS,NK,NB,NF

with:
::
  WRITE(99,'(4I6,F12.7,"  # of spin, # of k-points, # of bands, # of proj, Efermi" )') W%WDES%NCDIJ,NK,NB,NF,EFERMI

Another critical point for CSC calculations is the function call of
`LPRJ_LDApU` in VASP. This function is not needed, and was left there for debug
purposes, but is called every iteration. Removing the call to this function in `electron.F` in line 644 speeds up the calculation significantly in the `ICHARG=5` mode. Moreover, this prevents VASP from generating the `GAMMA` file, which should ideally only be done by the DMFT code after a successful DMFT step, and then be read by VASP.


Furthermore, there is a bug in `fileio.F` around line 1710 where VASP tries to
print "reading the density matrix from Gamma". This should be done only by the
master node, and VASP gets stuck sometimes. Adding a
::
  IF (IO%IU0>=0) THEN
  ...
  ENDIF
statement resolves this issue. A similar problem occurs, when VASP writes the
`OSZICAR` file and a buffer is stuck. Adding a `flush` to the buffer in
`electron.F` around line 580 after
::
  CALL STOP_TIMING("G",IO%IU6,"DOS")
  flush(17)
  print *, ' '

resolves this issue. Otherwise the OSZICAR file is not written properly.
