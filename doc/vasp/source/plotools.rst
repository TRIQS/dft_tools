.. highlight:: python

#########
PLO tools
#########

Introduction
************

This set of tools is intended for processing of raw projectors read
from VASP. One of the main tasks is to generate an orthonormalized subset
of PLOs constructed according to the :doc:`config-file </config>`.

To produce the final output the following steps are undertaken:

 * Parse input config-file

 * Input raw VASP data 

 * Convert the raw VASP data into an internal representaion and check it for consistency.
 
 * Generate a set of projector shells according to the parameters of the config-file

 * Create a set of projector groups

 * Perform necessary group operations (such as :ref:`orthogonalization<ortho>`) on the constituing shells 

 * Calculate and output some useful quantities (bare density matrix, DOS, etc.)


Initial Processing
******************

The raw data from VASP files is initially read in simple objects (containers).
Then these objects are combined in an another object containing all necessary
electronic structure data. At this stage simple consistency checks are performed:

 * the basic dimensions, such as the number of bands, number of `k`-points, etc.,
   are consistent for all data

 * the `k`-point mesh are read both the IBZKPT and EIGENVAL and it is worth checking
   that both sets are coinciding

 * in case tetrahedron data is read from IBZKPT, the tetrahedron volume must be related
   to the total volume of the unit cell as derived from POSCAR

All electronic structure from VASP is stored in a class ElectronicStructure:

.. autoclass:: elstruct.ElectronicStructure
   :members:


Consistency with parameters

 * parameters in the config-file should pass trivial checks such as that the ion
   list does not contain non-existing ions (boundary check for ion indices)

.. function:: check_vasp_data_consistency(conf_pars, vasp_data)

   **Parameters**:

   - *conf_pars* (dict) : dictionary of input parameters from conf-file
   - *vasp_data* (dict) : dictionary containing all VASP data

   **Returns**:
   
   *None*

   **Raises**:
   
   A meaningful exception indicating an inconsistency in the input data


Selecting projector subsets
===========================

The first step of PLO processing is to select subsets of projectors
corresponding to PLO groups. Each group contains a set of shells.
Each projector shell is represented by an object 'ProjectorShell'
that contains an array of projectors and information on the shell itself
(orbital number, ions, etc.). 'ProjectorShell's are contained in
both a list of shells (according to the original list as read
from config-file) and in a 'ProjectorGroup' object, the latter
also providing information about the energy window.
`[In fact, shell container can be a simple dictionary.]`

Order of operations:

  - transform projectors (all bands) in each shell
  - select transformed shell projectors for a given group within the window
  - orthogonalize if necessary projectors within a group by performing
    the following operations for each k-point:
    * combine all projector shells into a single array
    * orthogonalize the array
    * distribute back the arrays assuming that the order is preserved
    

.. autoclass:: proj_shell.ProjectorShell
   :members:


The class is using a helper function `select_bands()` for selecting a subset of bands.

.. function:: select_bands(eigvals, emin, emax)

   **Parameters**:

   - *eigvals* (numpy.array) : array of eigenvalues
   - *emin*, *emax* (float) : energy window

   **Returns**:
   
   *ib_win*, *nb_min*, *nb_max* (numpy.array[int], int, int) :
     lowest and highest indices of the selected bands
 

.. _ortho:

Orthogonalization
-----------------

At the second stage the selected projectors are orthogonalized (orthonormalized). 
Orthogonalization can be performed in different ways if projection is made
on several ions or if several correlated shells per ion are considered.
In the case of several correlated ions per unit cell (and one correlated shell per ion)
at least two options can be considered:

 #. Projectors are normalized for each ion separetely. In this case, corresponding
    Wannier functions for different ions are generally not orthogonal.

 #. Projectors are normalized for all ions in the unit cell simultaneously. This
    ensures that the Wannier functions for different ions are mutually orthogonal.

The way the normalization of a PLO group is done is controlled by two group parameters:

  - *NORMALIZE* (True/False) : indicates whether the PLO group is normalized (True by default)
  - *NORMION* (True/False) : indicates whether the PLO group is normalized on a per-ion basis
    (False by default)


Storing generated projectors
****************************

After the PLOs are generated they are stored to text files which can be subsequently
converted to TRIQS h5-files (using the converter). The general format of the file
is a JSON-header containing all necessary parameters followed by a set of arrays.
There is always one (control) file containing general information (`k`-kpoints, lattice vectors etc.)
and `at least` one file containing correlated groups (one file for each group).

Control file format
===================

Filename '<namebase>.ctrl'. Contains the data shared between all shells.
The JSON-header consists of the following elements:

  * *nk*: number of `k`-points

  * *ns*: number of spin channels

  * *nc_flag*: collinear/noncollinear case (False/True)

  * *ng*: number of projector groups

  * Symmetry information (list of symmetry operations)

  * *efermi*: Fermi level (optional)

  * Lattice information

Projector-group file format
===========================

Projector-group files have names '<namebase>.plog<Ng>'.
They essentially contain serialized objects of class 'ProjectorGroup'.
The JSON-header has, thus, the following elements:

  * *shells*: list of shells

  * each shell is a dictionary:

    - *lshell*: orbital number `l`

    - *nion*: number of ions

    - *ndim*: number of orbitals/ion

    - *nbmax*: maxmimum number of bands (needed for array allocations)

  * *emin*, *emax*: energy window


