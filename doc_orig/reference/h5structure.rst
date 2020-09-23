.. _hdfstructure:

hdf5 structure
==============

All the data is stored using the hdf5 standard, as described also in the
documentation of the TRIQS package itself. In order to do a DMFT calculation,
using input from DFT applications, a converter is needed on order to provide
the necessary data in the hdf5 format.

groups and their formats
------------------------

In order to be used with the DMFT routines, the following data needs to be
provided in the hdf5 file. It contains a lot of information in order to perform
DMFT calculations for all kinds of situations, e.g. d-p Hamiltonians, more than
one correlated atomic shell, or using symmetry operations for the k-summation.
We store all data in subgroups of the hdf5 archive:

Main data
^^^^^^^^^
There needs to be one subgroup for the main data of the
calculation. The default name of this group is `dft_input`. Its contents are

=================  ======================================================================  =====================================================================================
Name               Type                                                                    Meaning
=================  ======================================================================  =====================================================================================
energy_unit        numpy.float                                                             Unit of energy used for the calculation.
n_k                numpy.int                                                               Number of k-points used for the BZ integration.
k_dep_projection   numpy.int                                                               1 if the dimension of the projection operators depend on the k-point,
                                                                                           0 otherwise.
SP                 numpy.int                                                               1 for spin-polarised Hamiltonian, 0 for paramagnetic Hamiltonian.
SO                 numpy.int                                                               1 if spin-orbit interaction is included, 0 otherwise.
charge_below       numpy.float                                                             Number of electrons in the crystal below the correlated orbitals.
                                                                                           Note that this is for compatibility with dmftproj, otherwise set to 0
density_required   numpy.float                                                             Required total electron density. Needed to determine the chemical potential.
                                                                                           The density in the projection window is then `density_required`-`charge_below`.
symm_op            numpy.int                                                               1 if symmetry operations are used for the BZ sums,
                                                                                           0 if all k-points are directly included in the input.
n_shells           numpy.int                                                               Number of atomic shells for which post-processing is possible.
                                                                                           Note: this is `not` the number of correlated orbitals!
                                                                                           If there are two equivalent atoms in the unit cell, `n_shells` is 2.
shells             list of dict {string:int}, dim n_shells x 4                             Atomic shell information.
                                                                                           For each shell, have a dict with keys ['atom', 'sort', 'l', 'dim'].
                                                                                           'atom' is the atom index, 'sort' defines the equivalency of the atoms,
                                                                                           'l' is the angular quantum number, 'dim' is the dimension of the atomic shell.
                                                                                           e.g. for two equivalent atoms in the unit cell, `atom` runs from 0 to 1,
                                                                                           but `sort` can take only one value 0.
n_corr_shells      numpy.int                                                               Number of correlated atomic shells.
                                                                                           If there are two correlated equivalent atoms in the unit cell, `n_corr_shells` is 2.
n_inequiv_shells   numpy.int                                                               Number of inequivalent atomic shells. Needs to be smaller than `n_corr_shells`.
                                                                                           The up / downfolding routines mediate between all correlated shells and the
                                                                                           actual inequivalent shells, by using the self-energy etc. for all equal shells
                                                                                           belonging to the same class of inequivalent shells. The mapping is performed with
                                                                                           information stored in `corr_to_inequiv` and `inequiv_to_corr`.
corr_to_inequiv    list of numpy.int, dim `n_corr_shells`                                  mapping from correlated shells to inequivalent correlated shells.
                                                                                           A list of length `n_corr_shells` containing integers, where same numbers mark
                                                                                           equivalent sites.
inequiv_to_corr    list of numpy.int, dim `n_inequiv_shells`                               A list of length `n_inequiv_shells` containing list indices as integers pointing
                                                                                           to the corresponding sites in `corr_to_inequiv`.
corr_shells        list of dict {string:int}, dim n_corr_shells x 6                        Correlated orbital information.
                                                                                           For each correlated shell, have a dict with keys
                                                                                           ['atom', 'sort', 'l', 'dim', 'SO', 'irrep'].
                                                                                           'atom' is the atom index, 'sort' defines the equivalency of the atoms,
                                                                                           'l' is the angular quantum number, 'dim' is the dimension of the atomic shell.
                                                                                           'SO' is one if spin-orbit is included, 0 otherwise, 'irep' is a dummy integer 0.
use_rotations      numpy.int                                                               1 if local and global coordinate systems are used, 0 otherwise.
rot_mat            list of numpy.array.complex,                                            Rotation matrices for correlated shells, if `use_rotations`.
                   dim n_corr_shells x [corr_shells['dim'],corr_shells['dim']]             These rotations are automatically applied for up / downfolding.
                                                                                           Set to the unity matrix if no rotations are used.
rot_mat_time_inv   list of numpy.int, dim n_corr_shells                                    If `SP` is 1, 1 if the coordinate transformation contains inversion, 0 otherwise.
                                                                                           If `use_rotations` or `SP` is 0, give a list of zeros.
n_reps             numpy.int                                                               Number of irreducible representations of the correlated shell.
                                                                                           e.g. 2 if eg/t2g splitting is used.
dim_reps           list of numpy.int, dim n_reps                                           Dimension of the representations.
                                                                                           e.g. [2,3] for eg/t2g subsets.
T                  list of numpy.array.complex,                                            Transformation matrix from the spherical harmonics to impurity problem basis
                   dim n_inequiv_corr_shell x                                              normally the real cubic harmonics).
                   [max(corr_shell['dim']),max(corr_shell['dim'])]                         This matrix can be used to calculate the 4-index U matrix, not automatically done.
n_orbitals         numpy.array.int, dim [n_k,SP+1-SO]                                      Number of Bloch bands included in the projection window for each k-point.
                                                                                           If SP+1-SO=2, the number of included bands may depend on the spin projection up/down.
proj_mat           numpy.array.complex,                                                    Projection matrices from Bloch bands to Wannier orbitals.
                   dim [n_k,SP+1-SO,n_corr_shells,max(corr_shell['dim']),max(n_orbitals)]  For efficient storage reasons, all matrices must be of the same size
                                                                                           (given by last two indices).
                                                                                           For k-points with fewer bands, only the first entries are used, the rest are zero.
                                                                                           e.g. if number of Bloch bands ranges from 4-6, all matrices are of size 6.
bz_weights         numpy.array.float, dim n_k                                              Weights of the k-points for the k summation. Soon be replaced by `kpt_weights`
hopping            numpy.array.complex,                                                    Non-interacting Hamiltonian matrix for each k point.
                   dim [n_k,SP+1-SO,max(n_orbitals),max(n_orbitals)]                       As for `proj_mat`, all matrices have to be of the same size.
=================  ======================================================================  =====================================================================================

Converter specific data
^^^^^^^^^^^^^^^^^^^^^^^

This data is specific to the different converters and stored in the `dft_input`
group as well.

For the Vasp converter:

=================  ======================================================================  =====================================================================================
Name               Type                                                                    Meaning
=================  ======================================================================  =====================================================================================
kpt_basis          numpy.array.float, dim 3x3                                              Basis for the k-point mesh, reciprocal lattice vectors.
kpts               numpy.array.float, dim n_k x 3                                          k-points given in reciprocal coordinates.
kpt_weights        numpy.array.float, dim n_k                                              Weights of the k-points for the k summation.
proj_or_hk         string                                                                  Switch determining whether the Vasp converter is running in projection mode `proj`, or
                                                                                           in Hamiltonian mode `hk`. In Hamiltonian mode, the hopping matrix is written in
                                                                                           orbital basis, whereas in projection mode hopping is written in band basis.
proj_mat_csc       numpy.array.complex,                                                    Projection matrices from Bloch bands to Wannier orbitals for Hamiltonian based `hk`
                   dim                                                                     approach. No site index is given, since hk is written in orbital basis. The last to
                   [n_k,SP+1-SO, n_corr_shells x max(corr_shell['dim']), max(n_orbitals)]  indices are a square matrix rotating from orbital to band space.
dft_fermi_weights  numpy.array.float, dim n_k x 1 x max(n_orbitals)                        DFT fermi weights (occupations) of KS eigenstates for each k-point for calculation
                   (stored in dft_misc_input)                                              of density matrix correction.
band_window        list of numpy.array.int , dim(SP+1-SO)x n_k x 2                         Band windows as KS band indices in Vasp for each spin channel, and k-point. Needed for
                   (stored in dft_misc_input)                                              writing out the GAMMA file.
=================  ======================================================================  =====================================================================================



Symmetry operations
^^^^^^^^^^^^^^^^^^^
In this subgroup we store all the data for applying the symmetry operations in
the DMFT loop (in case you want to use symmetry operations). The default name
of this subgroup is `dft_symmcorr_input`. This information is needed only if symmetry
operations are used to do the k summation. To be continued...

.. warning::
   TO BE COMPLETED!

General and simple H(k) Converter
---------------------------------

The above described converter of the Wien2k input is quite involved, since
Wien2k provides a lot of information, e.g. about symmetry operations, that can
be used in the calculation. However, sometimes we want to use a light
implementation where the input consists basically only of the Hamiltonian
matrix in Wannier basis, given at a grid of k points in the first Brillouin
zone. For this purpose, a simple converter is included in the package, called
:class:`HkConverter`, which is implemented for the simplest case of
paramagnetic DFT calculations without spin-orbit coupling. It reads a simple,
easy to construct text file, and produces an archive that can be used for the
DMFT calculations. An example input file for a structure with one correlated
site with 3 t2g orbitals in the unit cell contains the following:

  10               <- n_k

  1.0              <- density_required

  1                <- n_shells

  1 1 2 3          <- shells, as above: atom, sort, l, dim

  1                <- n_corr_shells

  1 1 2 3 0 0      <- corr_shells, as above: atom, sort, l, dim, SO, dummy

  2 2 3            <- n_reps, dim_reps (length 2, because eg/t2g splitting) for each inequivalent correlated shell

After this header, we give the Hamiltonian matrices for al the k-points. for
each k-point we give first the matrix of the real part, then the matrix of the
imaginary part. The projection matrices are set automatically to unity
matrices, no rotations, no symmetry operations are used. That means that the
symmetry sub group in the hdf5 archive needs not be set, since it is not used.
It is furthermore assumed that all k-points have equal weight in the k-sum.
Note that the input file should contain only the numbers, not the comments
given in above example.

The Hamiltonian matrices can be taken, e.g., from Wannier90, which constructs
the Hamiltonian in a maximally localized Wannier basis.

Note that with this simplified converter, no full charge self consistent
calculations are possible!
