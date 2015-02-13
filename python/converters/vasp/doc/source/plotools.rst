.. highlight:: python

PLO tools
#########

Introduction
============

This set of tools is intended for processing of raw projectors read
from VASP. One of the main tasks is to generate an orthonormalized subset
of PLOs constructed according to the :doc:`config-file </config>`.

As an input this sub-library accepts a dictionary with parameters from
the config-file and a dictionary containing VASP-data objects.
Generally, the following workflow is adopted:

 * The VASP data is checked for consistency. 
 
 * For each group of PLOs a corresponding subset of projectors and eigenvalues
   is selected according to the config-file.

 * The selected subsets of projectors are orthogonalized. In general, there are different ways
   how it can be done (see :ref:`Orthogonalization<ortho>`).


Initial Processing
==================

Consistency check
-----------------

The purpose of a consistency check is to make sure that the VASP data passed
to PLOtools are complete and originating from the same calculation.
In particular, the following things are supposed to be checked:

 * the basic dimensions, such as the number of bands, number of `k`-points, etc.,
   are consistent for all data

 * the `k`-point mesh are read both the IBZKPT and EIGENVAL and it is worth checking
   that both sets are coinciding

 * in case tetrahedron data is read from IBZKPT, the tetrahedron volume must be related
   to the total volume of the unit cell as derived from POSCAR

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
---------------------------

.. autoclass:: plotools.ProjectorSet
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
=================

At the second stage the selected projectors are orthogonalized (orthonormalized). 
Orthogonalization can be performed in different ways if projection is made
on several ions or if several correlated shells per ion are considered.
In the case of several correlated ions per unit cell (and one correlated shell per ion)
at least two options can be considered:

 #. Projectors are normalized for each ion separetely. In this case, corresponding
    Wannier functions for different ions are generally not orthogonal.

 #. Projectors are normalized for all ions in the unit cell simultaneously. This
    ensures that the Wannier functions for different ions are mutually orthogonal.

If more than one shells is considered (say, `p` and `d` orbitals), the
normalization can be imposed either for a combined set of shells or for each shell
separately.

The way the normalization of a PLO group is done is controlled by two flags:

  - **NORMALIZE** (True/False) : indicates whether the PLO group is normalized (True by default)
  - **NORMION** (True/False) : indicates whether the PLO group is normalized on a per-ion basis

If there are several PLO groups defined, the convention is the following: All PLO groups
marked with `NORMALIZE = True` are orthogonalized with respect to each other.

