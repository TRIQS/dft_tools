Input Config-file
=================

A config-file describes subsets of PLOs that are to be generated.
The PLOs are defined in terms of `shells` determined uniquely by an orbital
number `l` and a set of ions (nomrmally, of the same sort).
The shells are further combined into `groups` such that PLO in each group
are orthogonalized together. This allows to create several mutually orthogonal
subsets of PLOs. A group is characterized by a single projection energy window.

A config-file contains three types of sections:

 - **[General]** : providing information applicable to all projected shells
   (e.g. Fermi level)
 - **[Shell <Ns>]** : each section like this describes a PLO shell, with the index
   `Ns` used for referencing
 - **[Group <Ng>]** : describes shell groups

..
  It must contain at least one group defined by a section `[PLO Group 1]`.
  There is also an optional section `[General]` containing options that concern
  all PLO groups (e.g. `k`-mesh properties).

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

A config file must contain at least on group of PLOs. Each group is described
by a set of ions to which the projectors are applied, an atomic shell number
(:math:`l = 0,1,2,3` for `s,p,d,f`, respectively), and an energy window defining
the subset of bands from which the projectors are constructed.

In addition, one can define a real or complex transformation, which allows one
to produce projectors in a basis set different from the original one.

Below, the format of config-file is described.
All option names are case-insensitive.

Required parameters
-------------------

- **IONS**: ion indices as defined in POSCAR files
- **LSHELL**: atomic shell (values 0, 1, 2, 3 for `s,p,d,f` orbitals, respectively)
- **EMIN**, **EMAX**: the bottom and top of the energy window with respect to the Fermi level

Optional parameters
-------------------

- **RTRANSFORM**, **CTRANSFORM**: real or complex transformation matrix used to produce projectors
  in a different basis; the number of columns is determined by the size of the atomic shell

 
