Input Config-file
=================

A config-file describes subsets of PLOs that are to be generated.
The PLOs are defined in terms of `shells` determined uniquely by an orbital
number `l` and a set of ions (nomrmally, of the same sort).
The shells are further combined into `groups` such that PLO in each group
are orthogonalized together. This, for instance, allows to create a group of several
mutually orthogonal subsets of PLOs.
A group is characterized by a single projection energy window.

A config-file contains three types of sections:

 - **[General]** : providing information applicable to all projected shells
   (e.g. Fermi level)
 - **[Shell <Ns>]** : each section like this describes a PLO shell, with the index
   `Ns` used for referencing
 - **[Group <tag>]** : describes shell groups

The format requirements are relatively flexible. A config-file must contain
at least one `[Shell]` section. If there is only one shell defined, it is possible
to specify the energy window by providing parameters `EMIN`, `EMAX` (see below)
right in this section, in which case a group
will be created automatically and the `[Group]` section can be omitted.
If, nevertheless, a group referencing the single shell is explicitly given
the energy window parameters provided in the `[Group]` have higher priority
and in case of conflict with `[Shell]` section a warning is issued.

An example of a config-file:

.. literalinclude:: adv_example.cfg

Here two shells, one corresponding to `d` states on ions 5-8, another to `p`
states of ions 9-20, are created. They form a single group that, by default, will be
orthogonalized within a window `[-7.6, 2.7]` eV. Also Fermi level is explicitly
specified, which might be necessary sometimes, e.g., for non-self-consistent calcuation
of the band structure.

Below, the sections and their parameters are described.
All parameter names are case-insensitive.

Section [General]
-----------------

**Required parameters:**

In principle, there are unconditionally required parameters in this section.
However, if VASP data file do not contain a meaningful value of the Fermi level
it must be given here using parameter *EFERMI*. Note that if this parameter
is given it will override a value that might have been read from VASP files.

**Optional parameters:**
 - *BASENAME* (string): a basename for output files. The filenames will be
   of the form '<basename>.*'.
 - *DOSMESH* (float, float, int): if this parameter is provided the projected
   DOS for each projected shell will be generated, with the energy mesh parameters
   given by the energy range (two floats) and the number of points (int). It is also
   possible to omit the energy range, in which case it will be set to the energy window
   of the corresponding projector group.
 
Section [Shell <Ns>]
--------------------

Defines a projected shell with an integer index `<Ns>`. Ideally, the indices should
run uniformly starting from 1. However, since they are used only to reference
shells in group sections, their values are not important. One should only
make sure that there are no sections with the same name, in which case one
of them will be ignored by config parser.

**Required parameters:**

 - *IONS* ([int]): provides a list of ions. The list can be either given
   by a explicit enumeration of ion indices or by a range `N1..N2` (ions `N1` and `N2`
   will be included). 
 - *LSHELL* (int): orbital number `l` of the shell (`l = 0,1,2,3` for `s,p,d,f`, respectively).


**Optional parameters:**

 - *TRANSFORM* (matrix of floats): transformation matrices
   (real or complex) that are applied to projectors before orthogonalization.
   The number of columns `Nc` must be consistent with the number of orbitals
   (`Nc = 2l+1` for real matrices and `Nc = 2(2l+1)` for complex ones).
   The dimension of the resulting orbital subspace is determined by the number of rows.
 - *TRANSFILE* (str): file containing transformation matrices for individual ions.
   The file must contain rows of floats. Empty lines and lines starting with '#' are ignored.
   The data is interpreted as follows:

   * The number of rows is divided by the number of ions to give the number
     of rows per ion, `Nr`.
   * The number of columns is divided by `Nm = 2 * l + 1` to give `nsize`.
     There are, then, three options:

     #. if `nc_flag = 0`, i.e. a calculation is collinear (no spin-orbit coupling),
        and `nsize = 1` the matrices are considered to be real;
     #. if `nc_flag = 0` and `nsize = 2` the matrices are considered to be complex;
     #. if `nc_flag = 1`, i.e. a calculation is non-collinear, it is expected
        that `nsize = 4`, and the matrices are considered to be complex and 
        spin-dependent.

   * The subspace dimension is determined simply as `Ndim = Nr / nsize`.

   In all cases when a division is performed the result must be integer. Otherwise,
   the matrices are considered to be inconsistent.
 - *EWINDOW* (float, float): energy window. Should be given only if no excplicit groups
   is specified. Otherwise, the values are overriden by group parameters.


Section [Group <tag>]
---------------------

Defines a group of shells. Note that the group tag can be any string without whitespaces.
It will be used to tag intermediate output files.

**Required parameters:**

 - *SHELLS* ([int]): indices refering to shells forming the group.
 - *EWINDOW* (float, float): the bottom and top of the energy window with respect to the Fermi level.

**Optional parameters:**
 - NORMALIZE (True/False): if True, orthogonalizetion is performed (default behavior).
 - *NORMION* (True/False): if True, orthogonalization is performed on each site
   separately; if False, all projectors of the group are orthogonalized together.

